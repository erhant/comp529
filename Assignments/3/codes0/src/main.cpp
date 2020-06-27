// SpMV Code
// Created: 03-12-2019
// Author: Najeeb Ahmad
// Updated: 13-05-2020
// Author: Muhammad Aditya Sasongko
// Updated: 20-05-2020
// Author: Erhan Tezcan

#include <bits/stdc++.h>
#include "common.h"
#include "matrix.h"
#include "mmio.h"

using namespace std;

#define VERBOSE 0 // 1 for specifics, 2 for sum verification

int main(int argc, char **argv)
{
	// Initial argument check
    if(argc < 3)
    {
    	std::cout << "Error: Missing arguments\n";
    	std::cout << "Usage: " << argv[0] << " <matrix.mtx> <time_step_count>";
    	return EXIT_FAILURE;
    }
    csr_matrix matrix;
    string matrix_name;
    printf("Reading .mtx file...\n");
    int retCode = 0;
    int time_steps = atoi(argv[2]);
    matrix_name = argv[1];
    cout << matrix_name << endl;
    double *rhs;
    double *result;
    retCode = mm_read_unsymmetric_sparse(argv[1], &matrix.m, &matrix.n, &matrix.nnz, &matrix.csrVal, &matrix.csrRowPtr, &matrix.csrColIdx); // loads the CSR format
    rhs = (double *)malloc(sizeof(double) * matrix.n); // Allocation for vector rhs    
    result = (double *)malloc(sizeof(double) * matrix.n); // Allocation for vector result
    if(retCode == -1) {
		cout << "Error reading input .mtx file\n";
		return EXIT_FAILURE;
	}
    
    printf("Matrix Rows: %d\n", matrix.m);
    printf("Matrix Cols: %d\n", matrix.n);
    printf("Matrix nnz: %d\n", matrix.nnz);
    coo2csr_in(matrix.m, matrix.nnz, matrix.csrVal, matrix.csrRowPtr, matrix.csrColIdx); 
    printf("Done reading file!\n");

    // Initialize right-hand-side
    for(int i = 0; i < matrix.n; i++)
      rhs[i] = (double) 1.0/matrix.n;

    clock_t start, end;

    start = clock();

    for(int k = 0; k < time_steps; k++) {
		// For each timestep

    	for(int i = 0; i < matrix.m; i++) {
			// For each row
			
        result[i]=0.0;
        for(int j = matrix.csrRowPtr[i]; j < matrix.csrRowPtr[i+1]; j++)
        {
            result[i] += matrix.csrVal[j] * rhs[matrix.csrColIdx[j]]; // from lec23, p. 24 |	y[i] += A[i][j] * x[j]
        }
      }

    	for(int i = 0; i < matrix.m; i++) {
		  	rhs[i] = result[i];
      }
      
    }

    end = clock(); 

    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    cout << "Time taken by program is : " << fixed  
         << time_taken << setprecision(5); 
    cout << " sec " << endl;  

    #if VERBOSE == 2
    // Verify sum
		long sum = 0;
		for (int i = 0; i<matrix.n && i<100; i++) {
			sum += result[i];
		}
		cout << "Sum is: " << sum << endl;
    #endif

    return EXIT_SUCCESS;
}
