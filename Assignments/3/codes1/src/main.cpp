// SpMV Code
// Created: 03-12-2019
// Author: Najeeb Ahmad
// Updated: 13-05-2020
// Author: Muhammad Aditya Sasongko
// Updated: 20-05-2020
// Author: Erhan Tezcan

/*
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)

int* csrRowPtr = new int[numtasks];
(unsigned int) ceil(((float) width) / ((float) block.x)
*/

#include <bits/stdc++.h>
#include <math.h>
#include <mpi.h> // MPI!
#include "common.h"
#include "matrix.h"
#include "mmio.h"

using namespace std;

#define DO_NOT_PROCEED -1
#define DO_PROCEED 1
#define RO_MSG_TAG 2 // row offsets
#define CI_MSG_TAG 3 // column indices
#define CV_MSG_TAG 4 // column values
#define VERBOSE 0

void printArray(int* arr, int n, string name) {
	cout << name << ":\n[";
	for (int i = 0; i<n-1; i++) {
		printf("%d ", arr[i]);
	}
	printf("%d]\n",arr[n-1]);
}

void printArray(double* arr, int n, string name) {
	cout << name << ":\n[";
	for (int i = 0; i<n-1; i++) {
		printf("%f ", arr[i]);
	}
	printf("%f]\n",arr[n-1]);
}

int main(int argc, char **argv) {
	// Initial argument check
	if (argc < 3) {
		cout << "Error: Missing arguments\n";
		cout << "Usage: " << argv[0] << " <matrix.mtx> <time_step_count>";
		return EXIT_FAILURE;
	}
	
	// MPI Initializations
	int numtasks, taskid;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	/*
	if(numtasks <= 1) {
		printf("Error: Number of tasks should be greater than 1. You have: %d\n",numtasks);
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	*/
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	// printf("MPI task %d has started...\n", taskid); // we don't really need this ;)
	int master = numtasks - 1;

	// Initializations for all processes
	int* recieveCounts = new int[numtasks]; // for allgatherv
	int* offsets = new int[numtasks]; // for all gatherv
	double *rhs;
	double *result;
	int time_steps = atoi(argv[2]);
	int i, j, k, t, myRow, rpw;
	short proceedCode = DO_PROCEED;

	// Above runs for everyone
	if (taskid == master) {
		/*** MASTER ***/

		// Master reads the matrix		
		int retCode = 0;
		csr_matrix matrix;
    	string matrix_name;
		matrix_name = argv[1];
		printf("Reading .mtx file...\n");
		cout << matrix_name << endl;
		retCode = mm_read_unsymmetric_sparse(argv[1], &matrix.m, &matrix.n, &matrix.nnz, &matrix.csrVal, &matrix.csrRowPtr, &matrix.csrColIdx); // loads the matrix
		// If we couldn't read the input, abort.
		if (retCode == -1) {
			printf("Error: Couldn't read input .mtx file!\n");
			proceedCode = DO_NOT_PROCEED;
			MPI_Bcast(&proceedCode, 1, MPI_SHORT, master, MPI_COMM_WORLD); // notify workers about abort
			MPI_Finalize();
			return EXIT_FAILURE;
		}
		// If this is not a square matrix, abort.
		if (matrix.m != matrix.n) {
			printf("Error: Expected square matrix, but input is not a square matrix!\n");
			proceedCode = DO_NOT_PROCEED;
			MPI_Bcast(&proceedCode, 1, MPI_SHORT, master, MPI_COMM_WORLD); // notify workers about abort
			MPI_Finalize();
			return EXIT_FAILURE;
		}
		// If matrix dimension is less than number of processors, abort.
		if (matrix.m < numtasks) {
			printf("Error: The matrix dimensions (%d,%d) are too small for %d processors... \n", matrix.m, matrix.n, numtasks);
			proceedCode = DO_NOT_PROCEED;
			MPI_Bcast(&proceedCode, 1, MPI_SHORT, master, MPI_COMM_WORLD); // notify workers about abort
			MPI_Finalize();
			return EXIT_FAILURE;
		}
		// If we didn't return by now, everything is ok.
		MPI_Bcast(&proceedCode, 1, MPI_SHORT, master, MPI_COMM_WORLD); // notify workers that they can continue

		// Broadcast matrix dimension
		MPI_Bcast(&(matrix.n), 1, MPI_INT, master, MPI_COMM_WORLD);

		// Read CSR
		printf("Matrix Rows: %d\nMatrix Cols: %d\nMatrix nnz: %d\n", matrix.m, matrix.n, matrix.nnz); // M is row count, N is col count, NNZ is number of non-zero elements
		coo2csr_in(matrix.m, matrix.nnz, matrix.csrVal, matrix.csrRowPtr, matrix.csrColIdx); 
		printf("Done reading file!\n");

		// Allocate and initialize rhs
		rhs = new double[matrix.n]; 
		for(i = 0; i < matrix.n; ++i) {
			rhs[i] = (double) 1.0/matrix.n;
		}
		
		// Calculate rpw
		rpw = int(floor(float(matrix.n) / float(numtasks))); // rows per worker

		// Construct recieve count array for the allgatherv
		for (i = 0; i<numtasks - 1; ++i) {
			recieveCounts[i] = rpw;
		}
		recieveCounts[numtasks - 1] = matrix.n - rpw * master; // master gets the fringes too

		// Construct offsets for the allgatherv
		offsets[0] = 0;
		for (i = 1; i<numtasks; ++i) {
			offsets[i] = offsets[i-1] + recieveCounts[i-1];
		}

		#if VERBOSE == 1
		// Print recieve counts
		printArray(recieveCounts, numtasks, "Recieve Counts");
		printArray(offsets, numtasks, "Offsets");
		#endif

		// Distribute RO (row offsets as csrRowPtr), CI (column indices as csrColIdx), CV (column values as csrVal)
		MPI_Request* sendRequests = new MPI_Request[3 * (numtasks-1)]; // 3 messages sent per processor except the master
		MPI_Status* sendStatuses = new MPI_Status[3 * (numtasks-1)];
		for (i = 0; i < numtasks - 1; ++i) {
			#if VERBOSE == 1
			printf("\nFor task %d\n", i);
			printArray(&(matrix.csrRowPtr[rpw * i]), rpw+1, "Row Pointers");
			printArray(&(matrix.csrColIdx[matrix.csrRowPtr[rpw * i]]), matrix.csrRowPtr[rpw * i + rpw] - matrix.csrRowPtr[rpw * i], "Column Indices");
			printArray(&(matrix.csrVal[matrix.csrRowPtr[rpw * i]]), matrix.csrRowPtr[rpw * i + rpw] - matrix.csrRowPtr[rpw * i], "Column Values");
			#endif
			// Send row pointers to slave
			MPI_Isend(&(matrix.csrRowPtr[rpw * i]), rpw+1, MPI_INT, i, RO_MSG_TAG, MPI_COMM_WORLD, &sendRequests[3*i]);
			// Send column indices to slave
			MPI_Isend(&(matrix.csrColIdx[matrix.csrRowPtr[rpw * i]]), matrix.csrRowPtr[rpw * i + rpw] - matrix.csrRowPtr[rpw * i], MPI_INT, i, CI_MSG_TAG, MPI_COMM_WORLD, &sendRequests[3*i+1]);			
			// Send column values to slave
			MPI_Isend(&(matrix.csrVal[matrix.csrRowPtr[rpw * i]]), matrix.csrRowPtr[rpw * i + rpw] - matrix.csrRowPtr[rpw * i], MPI_DOUBLE, i, CV_MSG_TAG, MPI_COMM_WORLD, &sendRequests[3*i+2]);
		}
		MPI_Waitall(3 * (numtasks-1), sendRequests, sendStatuses);		

		// Allocate result array
		result = new double[matrix.n]; // this will be the overall result

		// Calculations at master		
		myRow = rpw * taskid;
		#if VERBOSE == 1
		printf("\nFor MASTER task\n");
		printArray(&(matrix.csrRowPtr[myRow]), matrix.m - myRow + 1, "Row Pointers");
		printArray(&(matrix.csrColIdx[matrix.csrRowPtr[myRow]]), matrix.csrRowPtr[matrix.m] - matrix.csrRowPtr[myRow], "Column Indices");
		printArray(&(matrix.csrVal[matrix.csrRowPtr[myRow]]), matrix.csrRowPtr[matrix.m] - matrix.csrRowPtr[myRow], "Column Values");
		#endif

		clock_t start, end;
		start = clock();
		for (t = 0; t < time_steps; ++t) {
			k = matrix.csrRowPtr[myRow];
			// Compute locals for the MASTER (it also gets the fringe part)
			for (i = 0; i < matrix.m - myRow; ++i) { 
				// For each row I have:			
				result[myRow + i] = 0.0;
				for (j = 0; j < matrix.csrRowPtr[myRow + i + 1] - matrix.csrRowPtr[myRow + i]; ++j) { // rowPtr is length m+1
					result[myRow + i] += matrix.csrVal[k] * rhs[matrix.csrColIdx[k]]; // from lec23, p. 24 |	y[i] += A[i][j] * x[j]
					k++;
				}

			}

			#if VERBOSE == 1
			printArray(result, matrix.n, "Result at master BEFORE gather");
			#endif

			// Allgather the results 
			MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, result, recieveCounts, offsets, MPI_DOUBLE, MPI_COMM_WORLD); // 2nd and 3rd params are ignored like this

			#if VERBOSE == 1
			printArray(result, matrix.n, "Result at master AFTER gather");
			#endif

			// Swap result and rhs
			double* tmp = rhs; rhs = result; result = tmp;
			
		}

		// Master reports the calculation time
		end = clock(); 
		double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
		cout << "Time taken by program is : " << fixed << time_taken << setprecision(5) << " sec " << endl; 

		#if VERBOSE == 2
		// Verify sum
		long sum = 0;
		for (i = 0; i<matrix.n && i<100; ++i) {
			sum += rhs[i]; // we use rhs instead of result here because we swap the pointers at the end of the loop, so result is stored at rhs in the end
		}
		cout << "Sum is: " << sum << endl;
		#endif

		// Memory free for master
		delete[] sendRequests;
		delete[] sendStatuses;
	} else {
		/*** SLAVE ***/

		// Check if there were any errors that require termination
		MPI_Bcast(&proceedCode, 1, MPI_SHORT, master, MPI_COMM_WORLD);
		if (proceedCode == DO_NOT_PROCEED) {
			MPI_Finalize();
			return EXIT_FAILURE;
		}

		// Master is bcasting matrix.n, receive it
		int n;
		MPI_Bcast(&n, 1, MPI_INT, master, MPI_COMM_WORLD);

		// Allocate and initialize rhs
		rhs = new double[n]; 
		for(i = 0; i < n; ++i) {
			rhs[i] = (double) 1.0/n;
		}

		// Calculate rpw
		rpw = int(floor(float(n) / float(numtasks))); // rows per worker

		// Construct recieve count array for the allgatherv
		for (i = 0; i<numtasks - 1; ++i) {
			recieveCounts[i] = rpw;
		}
		recieveCounts[numtasks - 1] = n - rpw * master; // master gets the fringes too
		// Construct offsets for the allgatherv
		offsets[0] = 0;
		for (i = 1; i<numtasks; ++i) {
			offsets[i] = offsets[i-1] + recieveCounts[i-1];
		}

		// Allocate row pointer array (RO)
		int* csrRowPtr = new int[rpw+1]; // with 1 extra

		// Master is sending row pointers, receive it
		MPI_Recv(csrRowPtr, rpw+1, MPI_INT, master, RO_MSG_TAG, MPI_COMM_WORLD, &status);

		// Allocate column value and column index arrays
		const int valCount = csrRowPtr[rpw] - csrRowPtr[0]; // this will give the amount of values for my portion of rows
		int* csrColIdx = new int[valCount];
		double* csrVal = new double[valCount];

		// Master is sending column indices, receive it
		MPI_Recv(csrColIdx, valCount, MPI_INT, master, CI_MSG_TAG, MPI_COMM_WORLD, &status);

		// Master is sending column values, receive it
		MPI_Recv(csrVal, valCount, MPI_DOUBLE, master, CV_MSG_TAG, MPI_COMM_WORLD, &status);

		// Allocate result array
		result = new double[n]; // this will be the overall result

		// Calculations at slave
		myRow = rpw * taskid;
		for (t = 0; t < time_steps; ++t) {
			k = 0;
			// Calculations
			for (i = 0; i < rpw; ++i) {
				result[myRow + i] = 0.0;
				for (j = 0; j < csrRowPtr[i+1] - csrRowPtr[i]; ++j) {
					result[myRow + i] += csrVal[k] * rhs[csrColIdx[k]];
					k++;
				}
			}

			#if VERBOSE == 1
			printArray(result, n, "Result BEFORE gather at task" + to_string(taskid));
			#endif

			// Allgather the results 
			MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, result, recieveCounts, offsets, MPI_DOUBLE, MPI_COMM_WORLD); // 2nd and 3rd params are ignored like this

			#if VERBOSE == 1
			printArray(result, n, "Result AFTER gather at task" + to_string(taskid));
			#endif

			// Swap result and rhs
			double* tmp = rhs; rhs = result; result = tmp;
		}	

		// Memory free for workers
		delete[] csrRowPtr;
		delete[] csrColIdx;
		delete[] csrVal;
	}

    // Memory free for everyone
	delete[] rhs;
	delete[] result;   
	delete[] recieveCounts;
	delete[] offsets;

	// Finalize
	MPI_Finalize();
    return EXIT_SUCCESS;
}