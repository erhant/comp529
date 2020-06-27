/* Random number generation routines based on system 48-bit LC RNG:
 * Written by Thomas Kammeyer, UCSD
 *
 * double real_rand() automatically seeds itself if need be
 *                    and returns random numbers from [0,1)
 *
 * 05/22/96 jcm Correction to use of gettimeofday
 *
 * Requirements: Be sure and load with -lm
 */

/* 
 * This scheme uses gettmeofday() to seed the random number generator,
 * and has obvious flaws.
 *
 */
 
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#if 0
extern "C" {
   extern double drand48();
   extern void srand48(long);
}

#endif


int seed_rand(long sd)
{
static int seed_me = 0;  /* "...Seymour, seed me all night long..." */
   if (sd){
      seed_me = sd;
      }
   else{
      struct timeval tp;
      gettimeofday(&tp,NULL);
      long hashed = ((tp.tv_sec & ~tp.tv_usec) | (~tp.tv_sec & tp.tv_usec));
      seed_me = hashed;
   }
   srand48( seed_me );
   return(seed_me);
}

double real_rand()
{

   return drand48();
}
