#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h> 
#include <string.h>
#include <omp.h>

#define VERBOSE 0
#define N 14 // board size

bool can_be_placed(int board[N][N], int row, int col);
void print_solution(int board[N][N]);
void solve_NQueens(int board[N][N], int col);

int CUTOFF_DEPTH = N; // if >= N there will be no CUTOFF
int NUM_SOLUTIONS = 0;

#define MATCH(s) (!strcmp(argv[ac], (s)))

bool can_be_placed(int board[N][N], int row, int col) 
{ 
    int i, j;   
    /* Check this row on left side */
    for (i = 0; i < col; i++) 
        if (board[row][i]) 
            return false; 
  
    /* Check upper diagonal on left side */
    for (i=row, j=col; i>=0 && j>=0; i--, j--) 
        if (board[i][j]) 
            return false; 
  
    /* Check lower diagonal on left side */
    for (i=row, j=col; j>=0 && i<N; i++, j--) 
        if (board[i][j]) 
            return false; 
  
    return true; 
} 

void print_solution(int board[N][N]) 
{ 
    
    static int k = 1; 
    NUM_SOLUTIONS++;

    #if VERBOSE == 1
    //printf("Solution (from thread %d) #%d-\n",omp_get_thread_num(),k++); // debug purposes
    printf("Solution #%d-\n",k++);
    int i,j;
    for (i = 0; i < N; i++) 
    { 
        for (j = 0; j < N; j++) 
            printf(" %d ", board[i][j]); 
        printf("\n"); 
    } 
    printf("\n"); 
    #endif
} 

void solve_NQueens(int board[N][N], int col) 
{ 

    if (col == N) 
    { 
        // only one thread should print the solution
        #pragma omp critical 
        print_solution(board); 
        //SOLUTION_EXISTS = true;
        //return true; 
    } 
    else 
    {
        int i;
        for (i = 0; i < N; i++) 
        {             
            if (can_be_placed(board, i, col)) 
            { 
                // If i can place it, I place it and launch the next check to another task with the placed board
                board[i][col] = 1; 
                // Concurrent tasks cannot operate on the same board, so we have to sacrifice some memory for the task.
                int boardCopy[N][N];
                memcpy(boardCopy, board, sizeof(int) * N * N);
                // I just launch the task here, let another thread deal with the solution. 
                // as condition we can have if(false) or final(true) and this will execute on the point.
                // https://stackoverflow.com/questions/53641799/openmp-task-if-or-final
                // http://www.red-ricap.org/documents/1071192/1486573/OpenMP_02_TaskingInDepth.pdf/3b11710c-486a-47b8-907e-f1f3e0f9b605
                // we can use column number as tree depth
                #pragma omp task firstprivate(col, boardCopy) final(col > CUTOFF_DEPTH)
                solve_NQueens(boardCopy, col + 1);

                // Then I remove that place on my board.
                board[i][col] = 0;  // TODO: perhaps a small optimization to decrease this 2 writes + memcpy to 1 write + memcpy
            } 
        } 
        //return false; // i dont care about what I return (or do I?)
    }
  
    
}


int main(int argc,char **argv)
{ 
    int numthreads = 1;
    int ac;
    for(ac = 1; ac < argc; ac++)
    {
        if(MATCH("-c")) {CUTOFF_DEPTH = atoi(argv[++ac]);}
        else if(MATCH("-t")) {numthreads = atoi(argv[++ac]);}
        else {
            printf("\nUsage: %s [-c <cutoff depth>]\n",argv[0]);
            return(-1);
        }
    }
    printf("Using cutoff depth of %d with %d threads.\n",CUTOFF_DEPTH,numthreads);
    int board[N][N];

    memset(board, 0, sizeof(board)); 
    double time1 = omp_get_wtime();
    
    // create the team here, use export OMP_NUM_THREADS for the thread count
    #pragma omp parallel num_threads(numthreads)
    #pragma omp single
    {
        // how to wait for the completion of this?
        // https://docs.oracle.com/cd/E37069_01/html/E37081/gozsi.html
        #pragma omp taskgroup // we need taskgroup because we require child-level sychronization, as opposed to scope-limited sync of taskwait.
        solve_NQueens(board, 0);
    }
    
    printf("%d solutions found.\n",NUM_SOLUTIONS);
    printf("Elapsed time: %0.6lf\n", omp_get_wtime() - time1); 
    return 0;
} 