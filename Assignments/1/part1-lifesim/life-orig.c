/* **********************************************************
 * Sample Serial Code : Conways' game of life
 *
 * Usage : life [-n <meshpoints>]
 *              [-p <prob>]  [-t <timesteps>]
 *
 * 
 *  Author : Urvashi R.V. [04/06/2004]
 *     Modified by Scott Baden [10/8/06]
 *     Modified by Pietro Cicotti [10/8/08]
 *     Modified by Scott Baden [10/13/06]
 *     Modified by Didem Unat [03/06/15]
 *************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#define MATCH(s) (!strcmp(argv[ac], (s)))

int MeshPlot(int t, int m, int n, char **mesh);

double real_rand();
int seed_rand(long sd);

static char **currWorld=NULL, **nextWorld=NULL, **tmesh=NULL;
static int maxiter = 200; /* number of iteration timesteps */
static int population[2] = {0,0}; /* number of live cells */

int nx = 16;      /* number of mesh points in the x dimension */
int ny = 16;      /* number of mesh points in the y dimension */

static int w_update = 0;
static int w_plot = 1;

double getTime();
extern FILE *gnu;

int main(int argc,char **argv)
{
    int i,j,ac;

    /* Set default input parameters */
    
    float prob = 0.5;   /* Probability of placing a cell */
    long seedVal = 0;
    int game = 0;
    int s_step = 0;
    int numthreads = 1;
    int disable_display = 0;

    /* Over-ride with command-line input parameters (if any) */
    
    for(ac=1;ac<argc;ac++)
    {
        if(MATCH("-n")) {nx = atoi(argv[++ac]);}
        else if(MATCH("-i")) {maxiter = atoi(argv[++ac]);}
        else if(MATCH("-t"))  {numthreads = atof(argv[++ac]);}
        else if(MATCH("-p"))  {prob = atof(argv[++ac]);}
        else if(MATCH("-s"))  {seedVal = atof(argv[++ac]);}
        else if(MATCH("-d"))  {disable_display = 1;}
        else if(MATCH("-step"))  {s_step = 1;}
        else if(MATCH("-g"))  {game = atoi(argv[++ac]);}
        else {
            printf("Usage: %s [-n < meshpoints>] [-i <iterations>] [-p prob] [-step] [-game <game #>]\n",argv[0]);
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
    
    currWorld = (char**)malloc(sizeof(char*)*nx + sizeof(char)*nx*ny);
    for(i=0;i<nx;i++) 
      currWorld[i] = (char*)(currWorld+nx) + i*ny;
    
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
        for(i=1;i<nx-1;i++)
            for(j=1;j<ny-1;j++) {
                currWorld[i][j] = (real_rand() < prob);
                population[w_plot] += currWorld[i][j];
            }
    }
    else if (game == 1){ //  Block, still life
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
    
    printf("probability: %f\n",prob);
    printf("Random # generator seed: %d\n", rs);

    /* Plot the initial data */
    
    /* Perform updates for maxiter iterations */
    double t0 = getTime();
    int t;
    for(t=0;t<maxiter && population[w_plot];t++)
    {
        /* Use currWorld to compute the updates and store it in nextWorld */
        population[w_update] = 0;
        for(i=1;i<nx-1;i++)
            for(j=1;j<ny-1;j++) {
                int nn = currWorld[i+1][j] + currWorld[i-1][j] + 
		         currWorld[i][j+1] + currWorld[i][j-1] + 
                         currWorld[i+1][j+1] + currWorld[i-1][j-1] + 
		         currWorld[i-1][j+1] + currWorld[i+1][j-1];
                
		nextWorld[i][j] = currWorld[i][j] ? (nn == 2 || nn == 3) : (nn == 3);
                population[w_update] += nextWorld[i][j];
            }

        /* Wait for previous plot to end */
            
        /* Pointer Swap : nextWorld <-> currWorld */
        tmesh = nextWorld;
        nextWorld = currWorld;
        currWorld = tmesh;
        /* Update the ghost cells */

        /* Start the new plot */
        if(!disable_display) {
            // Actually plot the data
            MeshPlot(t,nx,ny,currWorld);
        } 

    if (s_step){
       printf("Finished with step %d\n",t);
       printf("Press enter to continue.\n");
       getchar();
    }

    }
    double t1 = getTime(); 
    printf("Running time for the iterations: %f sec.\n",t1-t0);
    printf("Press enter to end.\n");
    getchar();
    if(gnu != NULL)
      pclose(gnu);

    /* Wait for plot of final data */

    /* Free resources */
    free(nextWorld);
    free(currWorld);

    return(0);
}
