/* **********************************************************
 * Sample Serial Code : Conways' game of life
 *
 * 
 *  Author : Urvashi R.V. [04/06/2004]
 *     Modified by Scott Baden [10/8/06]
 *     Modified by Pietro Cicotti [10/8/08]
 *     Modified by Didem Unat [03/06/15]
 *************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

#define MATCH(s) (!strcmp(argv[ac], (s)))

int MeshPlot(int t, int m, int n, char **mesh);

double real_rand();
int seed_rand(long sd);

static char **currWorld=NULL, **nextWorld=NULL, **tmesh=NULL;
static int maxiter = 200; /* number of iteration timesteps */
static int population[2] = {0,0}; /* number of live cells */

int nx = 100;      /* number of mesh points in the x dimension (default 100)*/
int ny = 100;      /* number of mesh points in the y dimension (default 100)*/

static int w_update = 0;
static int w_plot = 1;

double getTime();
extern FILE *gnu;

FILE* fp; // File pointer for which we write the results to (then compare using diff)

int main(int argc,char **argv)
{
    int i,j,ac;

    /* Set default input parameters */
    float prob = 0.5;   /* Probability of placing a cell */
    long seedVal = 0;
    int game = 0;
    int s_step = 0;
    int numthreads = 1;
    int disable_display= 0;
    int disable_log = 0; 

    /* Over-ride with command-line input parameters (if any) */
    // ./life -i MAXITER -t NUMTHREAD -p PROB -s SEEDVAL -step SINGLESTEP -g GAMENO -d
    // -d takes no parameters but it will disable display, I set it to be -d by default. 
    for(ac = 1; ac < argc; ac++)
    {
        if(MATCH("-n")) {nx = atoi(argv[++ac]);}
        else if(MATCH("-i")) {maxiter = atoi(argv[++ac]);}
          else if(MATCH("-t"))  {numthreads = atof(argv[++ac]);}
        else if(MATCH("-p"))  {prob = atof(argv[++ac]);}
        else if(MATCH("-s"))  {seedVal = atof(argv[++ac]);}
        else if(MATCH("-step"))  {s_step = 1;}
        else if(MATCH("-d"))  {disable_display = 1;}
        else if(MATCH("-l"))  {disable_log = 1;}
        else if(MATCH("-g"))  {game = atoi(argv[++ac]);}
        else {
            printf("Usage: %s [-n < meshpoints>] [-i <iterations>] [-s seed] [-p prob] [-t numthreads] [-step] [-g <game #>] [-d]\n",argv[0]);
            return(-1);
        }
    }

    int rs = seed_rand(seedVal);

    /* Increment sizes to account for boundary ghost cells */
    nx = nx+2;
    ny = nx; 
    
    /* Allocate contiguous memory for two 2D arrays of size nx*ny.
     * Two arrays are required because in-place updates are not
     * possible with the simple iterative scheme listed below */
    currWorld = (char**)malloc(sizeof(char*)*nx + sizeof(char)*nx*ny); // update whole 2D at once
    for(i=0;i<nx;i++) 
      currWorld[i] = (char*)(currWorld+nx) + i*ny; // assign pointers here (instead of mallocing here)
    
    nextWorld = (char**)malloc(sizeof(char*)*nx + sizeof(char)*nx*ny);
    for(i=0;i<nx;i++) 
      nextWorld[i] = (char*)(nextWorld+nx) + i*ny;
    
    /* Set the boundary ghost cells to hold 'zero' */
    for(i=0;i<nx;i++)
    {
        currWorld[i][0]=0;
        currWorld[i][ny-1]=0;
        nextWorld[i][0]=0;
        nextWorld[i][ny-1]=0;
    }
    for(i=0;i<ny;i++)
    {
        currWorld[0][i]=0;
        currWorld[nx-1][i]=0;
        nextWorld[0][i]=0;
        nextWorld[nx-1][i]=0;
    }

    // Generate a world
    if (game == 0){ // Use Random input
        // Due to first-touch policy it would be sensible to parallelize this but the assignment PDF tells us not to do so.
        for(i=1;i<nx-1;i++) {
            for(j=1;j<ny-1;j++) {   
                currWorld[i][j] = (real_rand() < prob);
                population[w_plot] += currWorld[i][j];
            }
        }
            
    }
    else if (game == 1) { //  Block, still life
        printf("2x2 Block, still life\n");
        int nx2 = nx/2;
        int ny2 = ny/2;
        currWorld[nx2+1][ny2+1] = currWorld[nx2][ny2+1] = currWorld[nx2+1][ny2] = currWorld[nx2][ny2] = 1;
        population[w_plot] = 4;
    }
    else if (game == 2){ //  Glider (spaceship)
        printf("Glider (spaceship)\n");
        // Your code codes here
    }
    else{
        printf("Unknown game %d\n",game);
        exit(-1);
    }
    
    /* Plot the initial data */
    if(!disable_display) {
        MeshPlot(0,nx,ny,currWorld);
    }
      
    
    // Open the file for logging
    if (!disable_log) {
        char *filename;
        filename = (char *) malloc(30 * sizeof(char)); // 30 is arbitrary, but it should be enough to hold "log_SEEDNUMBER.txt"
        if (numthreads == 1) {
            sprintf(filename, "log_s_%d.txt", rs); // modify filename to hold the seed itself
        } else {
            sprintf(filename, "log_p_%d.txt", rs); // modify filename to hold the seed itself
        }
        
        fp = fopen(filename,"w"); // write mode
        free(filename); // free filename here just in case!
        if (fp == NULL) {
            // Error from fopen
            printf("Error occured during file opening.");
            exit(-2);
        }
    }
    

    // Auxillaries
    int t; 
    int sum = 0;
    short proceedToSwap; // flag variable
    double t0;
    
    // RUN
    if (numthreads == 1) {
        // SERIAL IMPLEMENTATION
        printf("Serial\n\tProbability: %f\n\tSeed: %d\n\tThreads: %d\n\tIterations: %d\n\tProblem Size: %d\n",prob,rs,numthreads,maxiter,nx);
        t0 = getTime();
        
        for(t=0;t<maxiter && population[w_plot];t++)
        {
            /* Use currWorld to compute the updates and store it in nextWorld */
            population[w_update] = 0;
            for(i=1;i<nx-1;i++) {
                for(j=1;j<ny-1;j++) {
                    // Calculate neighbor count
                    int nn = currWorld[i+1][j] + currWorld[i-1][j] + 
                            currWorld[i][j+1] + currWorld[i][j-1] + 
                            currWorld[i+1][j+1] + currWorld[i-1][j-1] + 
                            currWorld[i-1][j+1] + currWorld[i+1][j-1];
                    // If alive check if you die, if dead check if you can produce.
                    nextWorld[i][j] = currWorld[i][j] ? (nn == 2 || nn == 3) : (nn == 3);
                    // Update population
                    population[w_update] += nextWorld[i][j];
                }
            }  
            
            /* Pointer Swap : nextWorld <-> currWorld */
            tmesh = nextWorld;
            nextWorld = currWorld;
            currWorld = tmesh;

            if(!disable_display) {
                // Actually plot the data
                MeshPlot(t,nx,ny,currWorld);
            } else if (!disable_log) {
                // Write the values to the file, we actually dont care about them in terms of visualization, we just want them to be the same! Therefore this is ran by a single thread.
                for(i=1;i<nx-1;i++) {
                    // Current world only
                    for(j=1;j<ny-1;j++) {            
                        fprintf(fp, "%d", currWorld[i][j] ? 1 : 0);
                    }
                    fprintf(fp, "\n");
                }
            } 

            // Is the singlestep option set?
            if (s_step) {
                printf("Finished with step %d\n",t);
                printf("Press enter to continue.\n");
                getchar();
            }
        }

    } else {
        // PARALLEL IMPLEMENTATION
        printf("Parallel\n\tProbability: %f\n\tSeed: %d\n\tThreads: %d\n\tIterations: %d\n\tProblem Size: %d\n",prob,rs,numthreads,maxiter,nx);
        t0 = getTime();

        // Do one iteration outside
        population[w_update] = 0;
        #pragma omp parallel num_threads(numthreads)
        {
            #pragma omp for reduction(+: sum) private(j) // we could use collapse(2) but then it would be 2D decomposition, so we go for only one loop parallelization.
            for(i=1;i<nx-1;i++) {
                for(j=1;j<ny-1;j++) {
                    // Calculate neighbor count
                    int nn = currWorld[i+1][j] + currWorld[i-1][j] + 
                            currWorld[i][j+1] + currWorld[i][j-1] + 
                            currWorld[i+1][j+1] + currWorld[i-1][j-1] + 
                            currWorld[i-1][j+1] + currWorld[i+1][j-1];
                    // If alive: check if you die, if dead: check if you can produce.
                    nextWorld[i][j] = currWorld[i][j] ? (nn == 2 || nn == 3) : (nn == 3);
                    // Update population (CRITICAL)   
                    sum  += nextWorld[i][j];
                }
            } 
            #pragma omp single nowait
            population[w_update] += sum;
        }

        /* Pointer Swap : nextWorld <-> currWorld */
        tmesh = nextWorld;
        nextWorld = currWorld;
        currWorld = tmesh;    

        // now inside the loop we will print the world, while calculating the next one
        // is it possible to perhaps combine these two parallel regions into one, while also including this iterative loop inside?
        for(t=1; t<maxiter && population[w_plot]; t++)
        {
            /* Use currWorld to compute the updates and store it in nextWorld */
            population[w_update] = 0; // do this outside and at the end (instead of start)
            sum = 0;
            proceedToSwap = 0; // flag variable

            #pragma omp parallel num_threads(2) if(numthreads> 1)
            {
                #pragma omp single 
                {
                    // Launch plot task
                    #pragma omp task private(i,j)
                    {
                        // We will plot the current one while the other task is calculating the next.
                        //int i_p, j_p; // index variable for (p)rinting.
                        /* Start the new plot */
                        if(!disable_display) {
                            // Actually plot the data
                            MeshPlot(t,nx,ny,currWorld);
                        } else  if (!disable_log) {
                            // Write the values to the file, we actually dont care about them in terms of visualization, we just want them to be the same! Therefore this is ran by a single thread.
                            for(i=1;i<nx-1;i++) {
                                // Current world only
                                for(j=1;j<ny-1;j++) {            
                                    fprintf(fp, "%d", currWorld[i][j] ? 1 : 0);
                                }
                                fprintf(fp, "\n");
                            }
                        }
                    } // end of plot task
                    
                    // Launch compute task
                    #pragma omp task 
                    {
                        //printf("3 Thread: %d\n",omp_get_thread_num());
                        #pragma omp parallel num_threads(numthreads-1) if(numthreads> 2) // nested parallel enable
                        {
                            //printf("4 Thread: %d\n",omp_get_thread_num()); // why is this always 0?
                            #pragma omp for reduction(+: sum) private(j) // we could use collapse(2) but then it would be 2D decomposition, so we go for only one loop parallelization.
                            for(i=1;i<nx-1;i++) {
                                for(j=1;j<ny-1;j++) {
                                    //printf("5 Thread: %d\n",omp_get_thread_num());
                                    // Calculate neighbor count
                                    int nn = currWorld[i+1][j] + currWorld[i-1][j] + 
                                            currWorld[i][j+1] + currWorld[i][j-1] + 
                                            currWorld[i+1][j+1] + currWorld[i-1][j-1] + 
                                            currWorld[i-1][j+1] + currWorld[i+1][j-1];
                                    // If alive check if you die, if dead check if you can produce.
                                    nextWorld[i][j] = currWorld[i][j] ? (nn == 2 || nn == 3) : (nn == 3);
                                    // Update population (CRITICAL)
                                    sum += nextWorld[i][j];
                                }
                            } 

                            #pragma omp single nowait
                            {
                                //printf("Thread %d reporting in to swap.\n",omp_get_thread_num()); // debug purposes
                                population[w_update] += sum;
                            }
                        }
                    } // end of compute task
                    
                } // end of single
            } // end of 2 threads

            // Swap
            /* Pointer Swap : nextWorld <-> currWorld */
            tmesh = nextWorld;
            nextWorld = currWorld;
            currWorld = tmesh;    
            
            // Is the singlestep option set?
            if (s_step) {
                printf("Finished with step %d\n",t);
                printf("Press enter to continue.\n");
                getchar();
            }
        }

        // We have print one more, because this was calculated at the last iteration
        if(!disable_display) {
            // Actually plot the data
            MeshPlot(t,nx,ny,currWorld);
        } else  if (!disable_log) {
            // Write the values to the file, we actually dont care about them in terms of visualization, we just want them to be the same! Therefore this is ran by a single thread.
            for(i=1;i<nx-1;i++) {
                // Current world only
                for(j=1;j<ny-1;j++) {            
                    fprintf(fp, "%d", currWorld[i][j] ? 1 : 0);
                }
                fprintf(fp, "\n");
            }
        }
    }

    double t1 = getTime(); 
    printf("Running time for the iterations: %f sec.\n",t1-t0);
    //printf("Press enter to end.\n");
    //getchar();
    
    if(gnu != NULL)
      pclose(gnu);

    // Close the file
    if (!disable_log) {
        fprintf(fp, "Population: {%d, %d}\n",population[0], population[1]);
        fclose(fp);
    }
    

    /* Free resources */
    free(nextWorld);
    free(currWorld);

    return(0);
}
