/*	
 * noise_remover.cpp
 *
 * This program removes noise from an image based on Speckle Reducing Anisotropic Diffusion
 * Y. Yu, S. Acton, Speckle reducing anisotropic diffusion, 
 * IEEE Transactions on Image Processing 11(11)(2002) 1260-1270 <http://people.virginia.edu/~sc5nf/01097762.pdf>
 * Original implementation is Modified by Burak BASTEM
 */

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
 // Reduction: Most naive version possible. Everything is in global memory.
 __global__ void reduction_kernel(unsigned char* image, int width, int height, int n_pixels, float* sum, float* sum2) {
	 int j = blockIdx.x * blockDim.x + threadIdx.x;
	 int i = blockIdx.y * blockDim.y + threadIdx.y;
	 if (j >= 0 && j < width && i >= 0 && i < height) {
		 long p = (i * width + j); // pixel index
		 float tmp = image[p]; // 1 global access
		 atomicAdd(&sum[0], tmp); // 1 global access
		 atomicAdd(&sum2[0], tmp * tmp); // 1 global access
	 }	
 }
 // Statistics: After reduction is finished, sum parameter will be holding the sum basically. The result will be written to std_dev.
 __global__ void statistics_kernel(float* sum, int n_pixels, float* sum2, float* std_dev) { 
	 //printf("Sum: %f\tSum2: %f\n", sum[0], sum2[0]);
	 float mean = sum[0] / n_pixels;  // 1 global access
	 float variance = (sum2[0] / n_pixels) - mean * mean;  // 1 global access
	 std_dev[0] = variance / (mean * mean);  // 1 global access
 }
 // Compute 1: 
 __global__ void compute1_kernel(unsigned char* image, int width, int height, int n_pixels, float* std_dev, float* north_deriv, float* south_deriv, float* west_deriv, float* east_deriv, float* diff_coef) {
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
 __global__ void compute2_kernel(float* diff_coef, int width, int height, float lambda, float* north_deriv, float* south_deriv, float* west_deriv, float* east_deriv, unsigned char* image) {
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
	 time_1 = get_time();
	 cudaStream_t stream;
	 
	 // (HOST) Part II: Parse command line arguments
	 if(argc<2) {
	   printf("Usage: %s [-i < filename>] [-iter <n_iter>] [-l <lambda>] [-o <outputfilename>]\n",argv[0]);
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
		 } else {
		 printf("Usage: %s [-i < filename>] [-iter <n_iter>] [-l <lambda>] [-o <outputfilename>]\n",argv[0]);
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
 
	 // (HOST & DEVICE) Part IV: Allocate variables
	 unsigned char *d_image;
	 float* d_sum;
	 float* d_sum2;
	 float* d_std_dev;
	 CUDA_CHECK(cudaMalloc(&d_image, n_pixels * sizeof(unsigned char))); // image is also allocated there
	 CUDA_CHECK(cudaMalloc(&d_north_deriv, n_pixels * sizeof(float))); // north derivative allocate on device
	 CUDA_CHECK(cudaMalloc(&d_south_deriv, n_pixels * sizeof(float))); // south derivative allocate on device
	 CUDA_CHECK(cudaMalloc(&d_west_deriv, n_pixels * sizeof(float))); // west derivative allocate on device
	 CUDA_CHECK(cudaMalloc(&d_east_deriv, n_pixels * sizeof(float))); // east derivative allocate on device
	 CUDA_CHECK(cudaMalloc(&d_diff_coef, n_pixels * sizeof(float))); // diffusion coefficient derivative allocate on device
	 CUDA_CHECK(cudaMalloc(&d_sum, 1 * sizeof(float))); // scalar
	 CUDA_CHECK(cudaMalloc(&d_sum2, 1 * sizeof(float))); // scalar
	 CUDA_CHECK(cudaMalloc(&d_std_dev, 1 * sizeof(float))); // scalar
 
	 // (HOST) Extra: Calculate block and grid parameters.	
	 dim3 block(DIM_THREAD_BLOCK_X, DIM_THREAD_BLOCK_Y, 1);
	 dim3 grid((unsigned int) ceil(((float) width) / ((float) block.x)), (unsigned int) ceil(((float) height) / ((float) block.y)),1); // +1 blocks for fringing
 
	 printf("\n\tBlocks per grid: (%d, %d, 1)\n\tThreads per block: (%d, %d, 1)\n\n",grid.x, grid.y, block.x, block.y);
 
	 // (DEVICE) Extra: Copy everything to device.
	 CUDA_CHECK(cudaStreamCreate(&stream)); // create stream
	 CUDA_CHECK(cudaMemcpy(d_image, image, n_pixels * sizeof(unsigned char), cudaMemcpyHostToDevice)); // transfer image
	 time_4 = get_time();	
 
	 // (DEVICE) Part V: Compute --- n_iter * (3 * height * width + 42 * (height-1) * (width-1) + 6) floating point arithmetic operations in totaL
	 for (int iter = 0; iter < n_iter; iter++) {
		 CUDA_CHECK(cudaMemsetAsync(d_sum, 0, 1 * sizeof(float), stream)); // set sum 0
		 CUDA_CHECK(cudaMemsetAsync(d_sum2, 0, 1 * sizeof(float), stream)); // set sum2 0
		 reduction_kernel<<<grid, block, 0, stream>>>(d_image, width, height, n_pixels, d_sum, d_sum2);
		 CUDA_GET_LAST_ERR(1);
		 statistics_kernel<<<1, 1, 0, stream>>>(d_sum, n_pixels, d_sum2, d_std_dev);
		 CUDA_GET_LAST_ERR(2);
		 compute1_kernel<<<grid, block, 0, stream>>>(d_image, width, height, n_pixels, d_std_dev, d_north_deriv, d_south_deriv, d_west_deriv, d_east_deriv, d_diff_coef);
		 CUDA_GET_LAST_ERR(3);
		 compute2_kernel<<<grid, block, 0, stream>>>(d_diff_coef, width, height, lambda, d_north_deriv, d_south_deriv, d_west_deriv, d_east_deriv, d_image);
		 CUDA_GET_LAST_ERR(4);
	 }
	 // (DEVICE) Extra: Retrieve image from device
	 CUDA_CHECK(cudaDeviceSynchronize());
	 CUDA_CHECK(cudaStreamDestroy(stream)); // destroy streams
	 CUDA_CHECK(cudaMemcpy(image, d_image, n_pixels * sizeof(unsigned char), cudaMemcpyDeviceToHost)); // transfer image
	 time_5 = get_time();
 
	 // (DEVICE & HOST) Part VI: Write image to file
	 stbi_write_png(outputname, width, height, pixelWidth, image, 0);
	 time_6 = get_time();
 
	 // (HOST) Part VII: Get the average of sum of pixels for testing and calculate GFLOPS
	 // FOR VALIDATION - DO NOT PARALLELIZE
	 float test = 0;
	 for (int i = 0; i < height; i++) {
		 for (int j = 0; j < width; j++) {
			 test += image[i * width + j];
		 }
	 }
	 test /= n_pixels;
 
	 float gflops = (float) (n_iter * 1E-9 * (3 * height * width + 42 * (height-1) * (width-1) + 6)) / (time_5 - time_4);
	 time_7 = get_time();
 
	 // (HOST & DEVICE) Part VII: Deallocate variables
	 stbi_image_free(image);
	 CUDA_CHECK(cudaFree(d_image)); // device's image
	 CUDA_CHECK(cudaFree(d_north_deriv)); // north derivative allocate on device
	 CUDA_CHECK(cudaFree(d_south_deriv)); // south derivative allocate on device
	 CUDA_CHECK(cudaFree(d_west_deriv)); // west derivative allocate on device
	 CUDA_CHECK(cudaFree(d_east_deriv)); // east derivative allocate on device
	 CUDA_CHECK(cudaFree(d_diff_coef)); // diffusion coefficient derivative allocate on device
	 CUDA_CHECK(cudaFree(d_sum)); // scalar
	 CUDA_CHECK(cudaFree(d_sum2)); // scalar
	 CUDA_CHECK(cudaFree(d_std_dev)); // scalar
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
 
 