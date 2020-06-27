/*	
* noise_remover.cpp
*
* This program removes noise from an image based on Speckle Reducing Anisotropic Diffusion
* Y. Yu, S. Acton, Speckle reducing anisotropic diffusion, 
* IEEE Transactions on Image Processing 11(11)(2002) 1260-1270 <http://people.virginia.edu/~sc5nf/01097762.pdf>
* Original implementation is Modified by Burak BASTEM
*/

// srun: srun -A users -p short -N 1 -n 1 --gres=gpu:1 --qos=users --pty $SHELL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define DIM_THREAD_BLOCK_X 32
#define DIM_THREAD_BLOCK_Y 32

#define CUDA_CHECK(call) \
	if((call) != cudaSuccess) { \
		cudaError_t err = cudaGetLastError(); \
		printf("CUDA error calling method \""#call"\" - err: %s\n", cudaGetErrorString(err)); \
	}

#define CUDA_GET_LAST_ERR(verbose_num) \
	GLOBAL_ERR = cudaGetLastError(); \
	if ( cudaSuccess != GLOBAL_ERR ){ \
		fprintf(stderr, "%d-cudaCheckError() failed: %s\n", verbose_num, cudaGetErrorString(GLOBAL_ERR)); \
	}

#define MATCH(s) (!strcmp(argv[ac], (s)))

cudaError_t GLOBAL_ERR; // used by macro

// returns the current time
static const double kMicro = 1.0e-6;
double get_time() {
	struct timeval TV;
	struct timezone TZ;
	const int RC = gettimeofday(&TV, &TZ);
	if(RC == -1) {
		printf("ERROR: Bad call to gettimeofday\n");
		return(-1);
	}
	return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );
}
 
// KERNELS
// Reduction: Reduce 2D image into to arrays reduarr_sum and reduarr_sum2, both 1D.
__global__ void reduction_kernel_1(unsigned char* image, int width, int height, int n_pixels, float* reduarr_sum1, float* reduarr_sum2) {
	unsigned int j = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int i = blockIdx.y * blockDim.y + threadIdx.y;
	__shared__ unsigned char sdata_1[DIM_THREAD_BLOCK_X * DIM_THREAD_BLOCK_Y];
	__shared__ unsigned char sdata_2[DIM_THREAD_BLOCK_X * DIM_THREAD_BLOCK_Y];
	if (j < width && i < height) {
		// Reduce (level 1) image (size n) -> reduction_array (size blockDim.x * blockDim.y)
		int k = i * width + j; // convert to 1D index on image
		int tid = threadIdx.y * blockDim.x + threadIdx.x; // convert to 1D index within block
		float image_k = image[k]; // each thread loads one pixel of image
		sdata_1[tid] = image_k; // for sum1
		sdata_2[tid] = image_k * image_k; // for sum2
		__syncthreads(); // wait for the load to complete
		for (unsigned int s = (DIM_THREAD_BLOCK_X * DIM_THREAD_BLOCK_Y)/2; s>0; s>>=1) {
			if (tid < s) {
				sdata_1[tid] += sdata_1[tid + s]; // sum
				sdata_2[tid] += sdata_2[tid + s] * sdata_2[tid + s]; // sum2
			}
			__syncthreads();
		}
		if (tid == 0) {
			reduarr_sum1[blockIdx.y * blockDim.x + blockIdx.x] = sdata_1[0];
			reduarr_sum2[blockIdx.y * blockDim.x + blockIdx.x] = sdata_2[0];
		}		
	}		
}

 // Reduce reduarr_sum and reduarr_sum2
__global__ void reduction_kernel_2(float* reduarr_sum1, float* reduarr_sum2, int n, float* sums) {
	unsigned int i = blockIdx.x * blockDim.x * 2 + threadIdx.x; // notice x2
	__shared__ unsigned char sdata_1[DIM_THREAD_BLOCK_X * DIM_THREAD_BLOCK_Y];
	__shared__ unsigned char sdata_2[DIM_THREAD_BLOCK_X * DIM_THREAD_BLOCK_Y];
	if (i <= n) {
		int tid = threadIdx.x; 
		sdata_1[tid] = reduarr_sum1[i] + reduarr_sum1[i + blockDim.x];
		sdata_2[tid] = reduarr_sum2[i] + reduarr_sum2[i + blockDim.x]; 
		__syncthreads();

		for (unsigned int s = blockDim.x/2; s>0; s>>=1) {
			if (tid < s) { 
				sdata_1[tid] += sdata_1[tid + s]; 
				sdata_2[tid] += sdata_2[tid + s];
			}
			__syncthreads();
		}
		if (tid == 0) {
			// At this point we only have few values left to sum, so we should just use atomic add rather than launching another reduction kernel
			atomicAdd(&sums[0], sdata_1[0]); // 1 global access
			atomicAdd(&sums[1], sdata_2[0]); // 1 global access
		}
	}		
}


// Statistics: After reduction is finished, sum parameter will be holding the sum basically. The result will be written to std_dev.
__global__ void statistics_kernel(int n_pixels, float* sums, float* std_dev) { 
	float mean = sums[0] / n_pixels;  // 1 global access
	float variance = (sums[1] / n_pixels) - mean * mean;  // 1 global access
	std_dev[0] = variance / (mean * mean);  // 1 global access
	// reset sums for next iteration (is this better to do it here or issue a memsetasync?)
	//sums[0] = 0;
	//sums[1] = 0;
}

// Compute 1: 
__global__ void compute1_kernel_v1(unsigned char* image, int width, int height, int n_pixels, float* std_dev, float* north_deriv, float* south_deriv, float* west_deriv, float* east_deriv, float* diff_coef) {
	unsigned int col = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int row = blockIdx.y * blockDim.y + threadIdx.y;
	__shared__ unsigned char sdata[DIM_THREAD_BLOCK_Y][DIM_THREAD_BLOCK_X+1]; // +1 to avoid bank conflicts
	float stddev = std_dev[0];
	long k;
	float image_k;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	// Load image to shared memory
	if (col < width && row < height) {	
		k = row * width + col;	// position of current element
		image_k = image[k];
		sdata[ty][tx] = image_k;
	}
	__syncthreads();
	if (col != 0 && col < width - 1 && row != 0 && row < height - 1) {		 
		float gradient_square, num, den, std_dev2, laplacian;
		float north, south, west, east, diff_coef_k;
		if (tx == 0) {
			west = image[k - 1] - image_k;	// west derivative on register --- 1 floating point arithmetic operations, 1 global access
		} else {
			west = sdata[ty][tx-1] - image_k;	// west derivative on register --- 1 floating point arithmetic operations
		} 
		if (ty == 0) {
			north = image[k - width] - image_k;	// north derivative on register --- 1 floating point arithmetic operations, 1 global access
		} else {
			north = sdata[ty-1][tx] - image_k;	// north derivative on register --- 1 floating point arithmetic operations
		}
		if (tx == blockDim.x - 1) {
			east = image[k + 1] - image_k;	// east derivative on register --- 1 floating point arithmetic operations, 1 global access
		} else {
			east = sdata[ty][tx+1] - image_k;	// east derivative on  register --- 1 floating point arithmetic operations
		}
		if (ty == blockDim.y - 1) {
			south = image[k + width] - image_k;	// south derivative on register --- 1 floating point arithmetic operations, 1 global access
		} else {
			south = sdata[ty+1][tx] - image_k;	// south derivative on register --- 1 floating point arithmetic operations
		}
		gradient_square = (north * north + south * south + west * west + east * east) / (image_k * image_k); // 9 floating point arithmetic operations
		laplacian = (north + south + west + east) / image_k; // 4 floating point arithmetic operations
		num = (0.5 * gradient_square) - ((1.0 / 16.0) * (laplacian * laplacian)); // 5 floating point arithmetic operations
		den = 1 + (.25 * laplacian); // 2 floating point arithmetic operations
		std_dev2 = num / (den * den); // 2 floating point arithmetic operations
		den = (std_dev2 - stddev) / (stddev * (1 + stddev)); // 4 floating point arithmetic operations
		diff_coef_k = 1.0 / (1.0 + den); // 2 floating point arithmetic operations
		if (diff_coef_k < 0) {
			diff_coef_k = 0;
		} else if (diff_coef_k > 1) {
			diff_coef_k = 1;
		}
		north_deriv[k] = north; // north derivative on register back to global memory, 1 global access
		south_deriv[k] = south; // south derivative on register back to global memory, 1 global access
		west_deriv[k] = west; // west derivative on register back to global memory, 1 global access
		east_deriv[k] = east; // east derivative on register back to global memory, 1 global access
		diff_coef[k] = diff_coef_k; // diff coef on register back to global memory, 1 global access
	}	
}

// Compute 1: 
__global__ void compute1_kernel_v2(unsigned char* image, int width, int height, int n_pixels, float* std_dev, float* north_deriv, float* south_deriv, float* west_deriv, float* east_deriv, float* diff_coef) {
	unsigned int col = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int row = blockIdx.y * blockDim.y + threadIdx.y;
	__shared__ unsigned char sdata[DIM_THREAD_BLOCK_Y][DIM_THREAD_BLOCK_X+1]; // +1 to avoid bank conflicts
	float stddev = std_dev[0];
	long k;
	float image_k;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	// Load image to shared memory
	if (col < width && row < height) {
		k = row * width + col;	// position of current element
		image_k = image[k];
		sdata[ty][tx] = image_k;
	}
	__syncthreads();
	if (col != 0 && col < width - 1 && row != 0 && row < height - 1) {		 
		float gradient_square, num, den, std_dev2, laplacian;
		float north, south, west, east, diff_coef_k;
		if (tx == 0 || ty == 0 || tx == blockDim.x - 1 || ty == blockDim.y - 1) {
			// doing halo checks here to whether use global memory or shared data actually slows things down
			north = image[k - width] - image_k;	// north derivative on register --- 1 floating point arithmetic operations, 1 global access
			south = image[k + width] - image_k;	// south derivative on register --- 1 floating point arithmetic operations, 1 global access
			west = image[k - 1] - image_k;	// west derivative on register --- 1 floating point arithmetic operations, 1 global access
			east = image[k + 1] - image_k;	// east derivative on register --- 1 floating point arithmetic operations, 1 global access
		} else {
			north = sdata[ty-1][tx] - image_k;	// north derivative on register --- 1 floating point arithmetic operations, 1 global access
			south = sdata[ty+1][tx] - image_k;	// south derivative on register --- 1 floating point arithmetic operations, 1 global access
			west = sdata[ty][tx-1] - image_k;	// west derivative on register --- 1 floating point arithmetic operations, 1 global access
			east = sdata[ty][tx+1] - image_k;	// east derivative on register --- 1 floating point arithmetic operations, 1 global access
		} 
		gradient_square = (north * north + south * south + west * west + east * east) / (image_k * image_k); // 9 floating point arithmetic operations
		laplacian = (north + south + west + east) / image_k; // 4 floating point arithmetic operations
		num = (0.5 * gradient_square) - ((1.0 / 16.0) * (laplacian * laplacian)); // 5 floating point arithmetic operations
		den = 1 + (.25 * laplacian); // 2 floating point arithmetic operations
		std_dev2 = num / (den * den); // 2 floating point arithmetic operations
		den = (std_dev2 - stddev) / (stddev * (1 + stddev)); // 4 floating point arithmetic operations
		diff_coef_k = 1.0 / (1.0 + den); // 2 floating point arithmetic operations
		if (diff_coef_k < 0) {
			diff_coef_k = 0;
		} else if (diff_coef_k > 1) {
			diff_coef_k = 1;
		}
		north_deriv[k] = north; // north derivative on register back to global memory, 1 global access
		south_deriv[k] = south; // south derivative on register back to global memory, 1 global access
		west_deriv[k] = west; // west derivative on register back to global memory, 1 global access
		east_deriv[k] = east; // east derivative on register back to global memory, 1 global access
		diff_coef[k] = diff_coef_k; // diff coef on register back to global memory, 1 global access
	}	
}

__global__ void compute1_kernel_v3(unsigned char* image, int width, int height, int n_pixels, float* std_dev, float* north_deriv, float* south_deriv, float* west_deriv, float* east_deriv, float* diff_coef) {
	unsigned int col = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int row = blockIdx.y * blockDim.y + threadIdx.y;
	__shared__ unsigned char sdata[DIM_THREAD_BLOCK_Y+2][DIM_THREAD_BLOCK_X+2+1]; // +1 to avoid bank conflicts
	float stddev = std_dev[0];
	long k;
	float image_k;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	// Load image to shared memory
	if (col <= width - 1 && row <= height - 1) {	
	   k = row * width + col;	// position of current element
	   image_k = image[k];
	   sdata[ty+1][tx+1] = image_k;
	}
	if (row != 0 && ty == 0) {
	   sdata[0][tx+1] = image[k - width]; // fill upper halo if needed
	}
	if (col != 0 && tx == 0) {
	   sdata[ty+1][0] = image[k - 1];
	}
	if (row != height - 1 && ty == blockDim.y - 1) {
	   sdata[blockDim.y+1][tx+1] = image[k + width];
	}
	if (col != width - 1 && tx == blockDim.x - 1) {
	   sdata[ty+1][blockDim.x+1] = image[k + 1];
	}
	__syncthreads();
	if (col != 0 && col < width - 1 && row != 0 && row < height - 1) {		 
		float gradient_square, num, den, std_dev2, laplacian;
		float north, south, west, east, diff_coef_k;
		ty = ty + 1;
		tx = tx + 1;
	   north = sdata[ty-1][tx] - image_k;	// north derivative on register --- 1 floating point arithmetic operations
	   south = sdata[ty+1][tx] - image_k;	// south derivative on register --- 1 floating point arithmetic operations
	   west = sdata[ty][tx-1] - image_k;	// west derivative on register --- 1 floating point arithmetic operations
	   east = sdata[ty][tx+1] - image_k;	// east derivative on register --- 1 floating point arithmetic operations
	   gradient_square = (north * north + south * south + west * west + east * east) / (image_k * image_k); // 9 floating point arithmetic operations
	   laplacian = (north + south + west + east) / image_k; // 4 floating point arithmetic operations
	   num = (0.5 * gradient_square) - ((1.0 / 16.0) * (laplacian * laplacian)); // 5 floating point arithmetic operations
	   den = 1 + (.25 * laplacian); // 2 floating point arithmetic operations
	   std_dev2 = num / (den * den); // 2 floating point arithmetic operations
	   den = (std_dev2 - stddev) / (stddev * (1 + stddev)); // 4 floating point arithmetic operations
	   diff_coef_k = 1.0 / (1.0 + den); // 2 floating point arithmetic operations
	   if (diff_coef_k < 0) {
		   diff_coef_k = 0;
	   } else if (diff_coef_k > 1) {
		   diff_coef_k = 1;
	   }
	   north_deriv[k] = north; // north derivative on register back to global memory, 1 global access
	   south_deriv[k] = south; // south derivative on register back to global memory, 1 global access
	   west_deriv[k] = west; // west derivative on register back to global memory, 1 global access
	   east_deriv[k] = east; // east derivative on register back to global memory, 1 global access
	   diff_coef[k] = diff_coef_k; // diff coef on register back to global memory, 1 global access
	}	
}

__global__ void compute1_kernel_v4(unsigned char* image, int width, int height, int n_pixels, float* std_dev, float* north_deriv, float* south_deriv, float* west_deriv, float* east_deriv, float* diff_coef) {
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	if (j > 0 && j < width - 1 && i > 0 && i < height - 1) {		 
		float gradient_square, num, den, std_dev2, laplacian;
		long k = i * width + j;	// position of current element
		// from global memory to registers 
		float std_dev_0 = std_dev[0]; // 1 global access
		float image_k = image[k]; // 1 global access
		float north, south, west, east;
		float diff_coef_k;
		// actual process
		north = image[k - width] - image_k;	// north derivative on register --- 1 floating point arithmetic operations, 1 global access
		south = image[k + width] - image_k;	// south derivative on register --- 1 floating point arithmetic operations, 1 global access
		west = image[k - 1] - image_k;	// west derivative on register --- 1 floating point arithmetic operations, 1 global access
		east = image[k + 1] - image_k;	// east derivative on register --- 1 floating point arithmetic operations, 1 global access
		north_deriv[k] = north;	// north derivative on register back to global memory, 1 global access
		south_deriv[k] = south;	// south derivative on register back to global memory, 1 global access
		west_deriv[k] = west;	// west derivative on register back to global memory, 1 global access
		east_deriv[k] = east;	// east derivative on register back to global memory, 1 global access
		gradient_square = (north * north + south * south + west * west + east * east) / (image_k * image_k); // 9 floating point arithmetic operations
		laplacian = (north + south + west + east) / image_k; // 4 floating point arithmetic operations
		num = (0.5 * gradient_square) - ((1.0 / 16.0) * (laplacian * laplacian)); // 5 floating point arithmetic operations
		den = 1 + (.25 * laplacian); // 2 floating point arithmetic operations
		std_dev2 = num / (den * den); // 2 floating point arithmetic operations
		den = (std_dev2 - std_dev_0) / (std_dev_0 * (1 + std_dev_0)); // 4 floating point arithmetic operations
		diff_coef_k = 1.0 / (1.0 + den); // 2 floating point arithmetic operations
		if (diff_coef_k < 0) {
		   diff_coef_k = 0;
		} else if (diff_coef_k > 1) {
		   diff_coef_k = 1;
		}
		diff_coef[k] = diff_coef_k; // 1 global access
	}	
}

// Compute 2:
__global__ void compute2_kernel_v1(float* diff_coef, int width, int height, float lambda, float* north_deriv, float* south_deriv, float* west_deriv, float* east_deriv, unsigned char* image) {
	unsigned int col = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int row = blockIdx.y * blockDim.y + threadIdx.y;
	__shared__ unsigned char sdata[DIM_THREAD_BLOCK_Y][DIM_THREAD_BLOCK_X+1]; // +1 to avoid bank conflicts
	float diff_coef_k;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	long k;
	// Load image to shared memory
	if (col < width && row < height) {	
		k = row * width + col;	// position of current element
		diff_coef_k = diff_coef[k];
		sdata[ty][tx] = diff_coef_k;
	}
	__syncthreads();	
	if (col > 0 && col < width - 1 && row > 0 && row < height - 1) {
		float divergence;
		if (tx == blockDim.x - 1 && ty == blockDim.y - 1) {
			divergence = diff_coef_k * north_deriv[k] + diff_coef[k + width] * south_deriv[k] + diff_coef_k * west_deriv[k] + diff_coef[k + 1] * east_deriv[k];
		} else if (tx == blockDim.x - 1) {
			divergence = diff_coef_k * north_deriv[k] + sdata[ty+1][tx] * south_deriv[k] + diff_coef_k * west_deriv[k] + diff_coef[k + 1] * east_deriv[k];			
		} else if (ty == blockDim.y - 1) {			
			divergence = diff_coef_k * north_deriv[k] + diff_coef[k + width] * south_deriv[k] + diff_coef_k * west_deriv[k] + sdata[ty][tx+1] * east_deriv[k];
		} else {
			divergence = diff_coef_k * north_deriv[k] + sdata[ty+1][tx] * south_deriv[k] + diff_coef_k * west_deriv[k] + sdata[ty][tx+1] * east_deriv[k];
		}
		image[k] = image[k] + 0.25 * lambda * divergence; // --- 3 floating point arithmetic operations, 2 global access (read and write)
	}
}
 
__global__ void compute2_kernel_v2(float* diff_coef, int width, int height, float lambda, float* north_deriv, float* south_deriv, float* west_deriv, float* east_deriv, unsigned char* image) {
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	if (j > 0 && j < width - 1 && i > 0 && i < height - 1) {
		long k = i * width + j; // get position of current element
		float diff_coef_k = diff_coef[k]; // 1 global access
		float divergence = diff_coef_k * north_deriv[k] + diff_coef[k + width] * south_deriv[k] + diff_coef_k * west_deriv[k] + diff_coef[k + 1] * east_deriv[k]; // --- 7 floating point arithmetic operations, 6 global accesses
		image[k] = image[k] + 0.25 * lambda * divergence; // --- 3 floating point arithmetic operations, 2 global access (read and write)
	}
}
int main(int argc, char *argv[]) {
	// (HOST) Part I: Allocate and initialize variables
	double time_0, time_1, time_2, time_3, time_4, time_5, time_6, time_7, time_8;	// time variables
	time_0 = get_time();
	const char *filename = "input.pgm";
	const char *outputname = "output.png";
	int width, height, pixelWidth, n_pixels;
	int n_iter = 50;
	float lambda = 0.5;
	float *d_north_deriv, *d_south_deriv, *d_west_deriv, *d_east_deriv; // directional derivatives (TO DEVICE)
	float *d_diff_coef;	// diffusion coefficient (TO DEVICE)
	short compute1_kernel_type = 4; // 1 = partially shared directionwise, 2 = partially shared halowise, 3 = fully shared, 4 = no shared
	short compute2_kernel_type = 2; // 1 = partially shared directionwise, 2 = no shared
	bool use_pinned = true;// not 1 = dont use it, 1 = use pinned memory
	time_1 = get_time();
	cudaStream_t stream;
	
	// (HOST) Part II: Parse command line arguments
	if(argc<2) {
		printf("Usage: %s [-i < filename>] [-iter <n_iter>] [-l <lambda>] [-o <outputfilename>] [-conek 1|2|3|4] [-ctwok 1|2]\n",argv[0]);
		return(-1);
	}
	for(int ac=1;ac<argc;ac++) {
		if(MATCH("-i")) {
			filename = argv[++ac];
		} else if(MATCH("-iter")) {
			n_iter = atoi(argv[++ac]);
		} else if(MATCH("-l")) {
			lambda = atof(argv[++ac]);
		} else if(MATCH("-o")) {
			outputname = argv[++ac];			 
		} else if(MATCH("-conek")) { // Compute 1 Kernel type
	 		compute1_kernel_type = atoi(argv[++ac]);
		} else if(MATCH("-ctwok")) { // Compute 2 Kernel type
			compute2_kernel_type = atoi(argv[++ac]);
		} else if(MATCH("-usepinned")) { // Use Pinned Memory or not
			use_pinned = atoi(argv[++ac]) == 1;
		} else {
			printf("Usage: %s [-i < filename>] [-iter <n_iter>] [-l <lambda>] [-o <outputfilename>] [-ck 1|2|3|4]\n",argv[0]);
			return(-1);
		}
	}
	time_2 = get_time();
 
	// (HOST) Part III: Read image	
	printf("Reading image...\n");
	unsigned char *image = stbi_load(filename, &width, &height, &pixelWidth, 0);
	if (!image) {
		fprintf(stderr, "Couldn't load image.\n");
		return (-1);
	}
	printf("Image Read. Width : %d, Height : %d, nComp: %d\n",width,height,pixelWidth);
	n_pixels = height * width; // W * H is the number of pixels
	time_3 = get_time();

	// (HOST) Extra: Make a pinned image at host
	unsigned char *pinned_host_image;
	if (use_pinned == 1) {
		CUDA_CHECK(cudaMallocHost((void**)&pinned_host_image, n_pixels * sizeof(unsigned char)));
		memcpy(pinned_host_image, image, n_pixels * sizeof(unsigned char));
	}
	
	// (HOST & DEVICE) Part IV: Allocate variables
	unsigned char *d_image;
	float* d_sums; // 0=sum, 1=sum2
	float* d_std_dev;
	CUDA_CHECK(cudaMalloc(&d_image, n_pixels * sizeof(unsigned char))); // image is also allocated there
	CUDA_CHECK(cudaMalloc(&d_north_deriv, n_pixels * sizeof(float))); // north derivative allocate on device
	CUDA_CHECK(cudaMalloc(&d_south_deriv, n_pixels * sizeof(float))); // south derivative allocate on device
	CUDA_CHECK(cudaMalloc(&d_west_deriv, n_pixels * sizeof(float))); // west derivative allocate on device
	CUDA_CHECK(cudaMalloc(&d_east_deriv, n_pixels * sizeof(float))); // east derivative allocate on device
	CUDA_CHECK(cudaMalloc(&d_diff_coef, n_pixels * sizeof(float))); // diffusion coefficient derivative allocate on device
	CUDA_CHECK(cudaMalloc(&d_sums, 2 * sizeof(float))); // scalar
	CUDA_CHECK(cudaMalloc(&d_std_dev, 1 * sizeof(float))); // scalar

	// (HOST) Extra: Calculate block and grid parameters.	
	dim3 block(DIM_THREAD_BLOCK_X, DIM_THREAD_BLOCK_Y, 1);
	dim3 grid((unsigned int) ceil(((float) width) / ((float) block.x)), (unsigned int) ceil(((float) height) / ((float) block.y)),1); // +1 blocks for fringing

	// (HOST) Extra: Reduction allocations
	/// 2D image to first level reduction done by reduction_1_grid which is same as grid
	// First level reduction to second level reduction
	dim3 block_1D(DIM_THREAD_BLOCK_X * DIM_THREAD_BLOCK_Y, 1, 1);
	dim3 reduction_2_grid((unsigned int) ceil(((float) (grid.x * grid.y)) / (float) block_1D.x));
	printf("\n\tBlocks per grid: (%d, %d, 1)\n\tThreads per block: (%d, %d, 1)\n\n",grid.x, grid.y, block.x, block.y);
	printf("\tReduction Info: \n\t\tLvl. 1: (%d,%d) image reduced to (%d * %d) array.\n\t\tLvl. 2: That (%d) reduced to a (%d) array.\n\n", width, height, grid.x, grid.y, grid.x * grid.y, reduction_2_grid.x);
	float* reduarr_sum1;
	float* reduarr_sum2;
	CUDA_CHECK(cudaMalloc(&reduarr_sum1, grid.x * grid.y * sizeof(float)));
	CUDA_CHECK(cudaMalloc(&reduarr_sum2, grid.x * grid.y * sizeof(float)));

	// (DEVICE) Extra: Copy everything to device.
	CUDA_CHECK(cudaStreamCreate(&stream)); // create stream
	if (use_pinned) {
		CUDA_CHECK(cudaMemcpy(d_image, pinned_host_image, n_pixels * sizeof(unsigned char), cudaMemcpyHostToDevice)); // transfer image
	} else {
		CUDA_CHECK(cudaMemcpy(d_image, image, n_pixels * sizeof(unsigned char), cudaMemcpyHostToDevice)); // transfer image
	}	
	time_4 = get_time();	

	// (DEVICE) Part V: Compute --- n_iter * (3 * height * width + 42 * (height-1) * (width-1) + 6) floating point arithmetic operations in totaL
	for (int iter = 0; iter < n_iter; iter++) {
		CUDA_CHECK(cudaMemsetAsync(d_sums, 0, 2 * sizeof(float), stream)); // set sums to 0 at first
		reduction_kernel_1<<<grid, block, 0, stream>>>(d_image, width, height, n_pixels, reduarr_sum1, reduarr_sum2);
		CUDA_GET_LAST_ERR(11);
		reduction_kernel_2<<<reduction_2_grid, block_1D, 0, stream>>>(reduarr_sum1, reduarr_sum2, grid.x * grid.y, d_sums);
		CUDA_GET_LAST_ERR(12);
		statistics_kernel<<<1, 1, 0, stream>>>(n_pixels, d_sums, d_std_dev);
		CUDA_GET_LAST_ERR(2);
		if (compute1_kernel_type == 1) {
			// Direction wise if-else
			compute1_kernel_v1<<<grid, block, 0, stream>>>(d_image, width, height, n_pixels, d_std_dev, d_north_deriv, d_south_deriv, d_west_deriv, d_east_deriv, d_diff_coef);
		} else if (compute1_kernel_type == 2) {
			// Halo-wise if-else
			compute1_kernel_v2<<<grid, block, 0, stream>>>(d_image, width, height, n_pixels, d_std_dev, d_north_deriv, d_south_deriv, d_west_deriv, d_east_deriv, d_diff_coef);
		} else if (compute1_kernel_type == 3) {
			// Fully shared memory
			compute1_kernel_v3<<<grid, block, 0, stream>>>(d_image, width, height, n_pixels, d_std_dev, d_north_deriv, d_south_deriv, d_west_deriv, d_east_deriv, d_diff_coef);
		} else if (compute1_kernel_type == 4) {
			// Same as noiseRemover v2, no shared memory
			compute1_kernel_v4<<<grid, block, 0, stream>>>(d_image, width, height, n_pixels, d_std_dev, d_north_deriv, d_south_deriv, d_west_deriv, d_east_deriv, d_diff_coef);
		} else {
			printf("Invalid compute 1 kernel type parameter.");
			exit(-1);
		}
		CUDA_GET_LAST_ERR(3);
		if (compute2_kernel_type == 1) {
			// Shared memory, same as v1 of compute 1. Direction-wise
			compute2_kernel_v1<<<grid, block, 0, stream>>>(d_diff_coef, width, height, lambda, d_north_deriv, d_south_deriv, d_west_deriv, d_east_deriv, d_image);
		} else if (compute2_kernel_type == 2) {
			// Same as noiseRemover v2, no shared memory
			compute2_kernel_v2<<<grid, block, 0, stream>>>(d_diff_coef, width, height, lambda, d_north_deriv, d_south_deriv, d_west_deriv, d_east_deriv, d_image);
		} else {
			printf("Invalid compute 2 kernel type parameter.");
			exit(-1);
		}
		CUDA_GET_LAST_ERR(4);
	}
	// (DEVICE) Extra: Retrieve image from device
	CUDA_CHECK(cudaDeviceSynchronize());
	CUDA_CHECK(cudaStreamDestroy(stream)); // destroy stream
	time_5 = get_time();

	// (DEVICE & HOST) Part VI: Retrieve and write image to file
	if (use_pinned) {
		CUDA_CHECK(cudaMemcpy(pinned_host_image, d_image, n_pixels * sizeof(unsigned char), cudaMemcpyDeviceToHost)); // transfer image
	} else {
		CUDA_CHECK(cudaMemcpy(image, d_image, n_pixels * sizeof(unsigned char), cudaMemcpyDeviceToHost)); // transfer image
	}
	
	stbi_write_png(outputname, width, height, pixelWidth, image, 0);
	time_6 = get_time();

	// (HOST) Part VII: Get the average of sum of pixels for testing and calculate GFLOPS
	// FOR VALIDATION - DO NOT PARALLELIZE
	float test = 0;
	if (use_pinned) { // we do this out here to avoid too many if checks in the loop
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				test += pinned_host_image[i * width + j];
			}
		}
	} else {
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				test += image[i * width + j];
			}
		}
	}
	test /= n_pixels;
 
	float gflops = (float) (n_iter * 1E-9 * (3 * height * width + 42 * (height-1) * (width-1) + 6)) / (time_5 - time_4);
	time_7 = get_time();

	// (HOST & DEVICE) Part VII: Deallocate variables
	if (use_pinned) {
		CUDA_CHECK(cudaFreeHost(pinned_host_image));
	}
	stbi_image_free(image);
	CUDA_CHECK(cudaFree(d_image)); // device's image
	CUDA_CHECK(cudaFree(d_north_deriv)); // north derivative allocate on device
	CUDA_CHECK(cudaFree(d_south_deriv)); // south derivative allocate on device
	CUDA_CHECK(cudaFree(d_west_deriv)); // west derivative allocate on device
	CUDA_CHECK(cudaFree(d_east_deriv)); // east derivative allocate on device
	CUDA_CHECK(cudaFree(d_diff_coef)); // diffusion coefficient derivative allocate on device
	CUDA_CHECK(cudaFree(d_sums)); // scalar
	CUDA_CHECK(cudaFree(d_std_dev)); // scalar
	CUDA_CHECK(cudaFree(reduarr_sum1)); // reduction auxiliary
	CUDA_CHECK(cudaFree(reduarr_sum2)); // reduction auxiliary
	time_8 = get_time();

	// (HOST) Final: Print
	printf("Time spent in different stages of the application:\n");
	printf("%9.6f s => Part I: allocate and initialize variables\n", (time_1 - time_0));
	printf("%9.6f s => Part II: parse command line arguments\n", (time_2 - time_1));
	printf("%9.6f s => Part III: read image\n", (time_3 - time_2));
	printf("%9.6f s => Part IV: allocate variables\n", (time_4 - time_3));
	printf("%9.6f s => Part V: compute\n", (time_5 - time_4));
	printf("%9.6f s => Part VI: write image to file\n", (time_6 - time_5));
	printf("%9.6f s => Part VII: get average of sum of pixels for testing and calculate GFLOPS\n", (time_7 - time_6));
	printf("%9.6f s => Part VIII: deallocate variables\n", (time_7 - time_6));
	printf("Total time: %9.6f s\n", (time_8 - time_0));
	printf("Average of sum of pixels: %9.6f\n", test);
	printf("GFLOPS: %f\n", gflops);
	return 0;
}
 
 