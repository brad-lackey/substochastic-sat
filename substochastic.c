/** @file  substochastic.c
 * @brief Main file for Substochastic Monte Carlo.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 4/2/16.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "macros.h"
#include "bitstring.h"
#include "sat.h"
#include "population.h"

//On machines with very old versions of glibc (e.g. the Raritan cluster)
//we need to define gnu_source in order to avoid warnings about implicit
//declaration of getline() and round().
#define _GNU_SOURCE

#define TEST 1

extern int nbts;
extern int arraysize;

static int popsize;

void update(double a, double b, double mean, Population P, int parity);
int parseCommand(int argc, char **argv, Population *Pptr, double *weight, double *runtime, int *trials, double *starttime);


/*Added by SPJ 3/17
 *Returns the average number of variables per clause.
 */
double avgLength(SAT instance) {
  int c;
  int total;
  total = 0;
  for(c = 0; c < instance->num_clauses; c++) total += instance->clause_length[c];
  return (double)total/(double)instance->num_clauses;
}

/*Added by SPJ 3/17
 *Automatically chooses parameters, such as runtime and number of trials.
 */
void autoparam(double varsperclause, int vars, double *weight, double *runtime, int *trials, double *starttime) {
  int k;
  k = (int)round(varsperclause);
  //*weight = 1.0;
  *weight = 100.0; // Changed since weight now percent -- Michael 3/30/16
  *trials = 5; //default
  *starttime = 0.0; // default -- Michael 3/30/16
  if(k == 2) {
    //These values are determined using the experimental results in
    //summary_2sat120v.txt and summary_2sat200v.txt.
    *trials = 10;
    *runtime = 3.2E4;
  }
  if(k == 3) {
    //These values are determined using the experimental results in
    //summary_3sat70v.txt and summary3sat110vlong.txt.
    *trials = 5;
    *runtime = 223.0*exp(0.07*vars);
  }
  if(k == 4 || vars > 200) {
    //This overrides the above.
    //Just go for broke and use up all the time.
    *trials = 5;
    *runtime = 3E6;
  }
}

int main(int argc, char **argv){
  int parity, trials, try, err;
  double weight, runtime, mean, starttime;
  double a, b, t, dt;
  Population pop;
  double min = -1.0;    //the best minimum from different trials
  Bitstring solution;   //the corresponding bitstring
  clock_t beg, end;     //for code timing
  double time_spent;    //for code timing
  
  
  
  beg = clock();
  
  if ( (err = parseCommand(argc, argv, &pop, &weight, &runtime, &trials, &starttime)) ){
    return err;
  }
  
  if ( (err = initBitstring(&solution)) ){
    fprintf(stderr, "Could not initialize answerspace.\n");
    return err;
  }
  
  
  for(try = 0; try < trials; try++) {
    
    t = starttime;
    parity = 0;
    randomPopulation(pop,popsize);
    
    while (t < runtime) {
      
      // The annealing schedule
      a = weight*(1.0 - t/runtime)/100.0; // Turned weight into percent -- Michael 3/30/16
      b = (t/runtime);
      
      mean = pop->avg_v + (pop->max_v - pop->min_v)*(popsize - pop->psize)/(2*popsize);
      
      if ( (pop->max_v - mean) > (mean - pop->min_v) )
        dt = 0.9/(a + b*(pop->max_v - mean));
      else
        dt = 0.9/(a + b*(mean - pop->min_v));
      if (t + dt > runtime)
        dt = runtime - t;
      
      update(a*dt, b*dt, mean, pop, parity);
      
      
      t += dt;
      parity ^= 1;
    }
    
#if TEST
    printBits(stdout, pop->winner);
    if ((min<0) || (pop->winner->potential < min)) {
      min = pop->winner->potential;
      copyBitstring(solution, pop->winner);
    }
#else
    end = clock();
    time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
    if ((min<0) || (pop->winner->potential < min)) {
      min = pop->winner->potential;
      copyBitstring(solution, pop->winner);
      printBits(stdout, pop->winner);
      printf("c Walltime: %f seconds\n", time_spent);
      fflush(stdout);
    }
    if ( time_spent > 240 )
      break;
#endif
    
  }
  
  freeBitstring(&solution);
  
#if TEST
  end = clock();
  time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
  printf("c Walltime: %f seconds\n", time_spent);
#endif
  
  return 0;
}


int parseCommand(int argc, char **argv, Population *Pptr, double *weight, double *runtime, int *trials, double *starttime){
  SAT sat;
  int seed;
  double varsperclause; //average number of variables per clause
  FILE *fp;
  Population pop;
  
  if (argc != 7 && argc != 2) {
    fprintf(stderr, "Usage: %s instance.cnf [<step weight> <runtime> <population size> <trials> <start time>]\n",argv[0]);
    return 2;
  }
  
  if ( (fp = fopen(argv[1], "r")) == NULL ){
    fprintf(stderr,"Could not open file %s\n",argv[1]);
    return IO_ERROR;
  }
  
  if ( loadDIMACSFile(fp,&sat) ){
    fprintf(stderr,"Error reading in DIMACS SAT file %s\n",argv[1]);
    return IO_ERROR;
  }
  
  fclose(fp);
  
  setBitLength(sat->num_vars);
  varsperclause = avgLength(sat);
  
  if ( argc == 2 )
    autoparam(varsperclause, nbts, weight, runtime, trials, starttime);
  
  if( argc == 7 ) {
    *weight = (double) atoi(argv[2]);
    sscanf(argv[3], "%lf", runtime); //allows scientific notation unlike atoi
    popsize = atoi(argv[4]);
    *trials = atoi(argv[5]);
    *starttime = atoi(argv[6]);
  } else {
    popsize = 64;
  }
  
  arraysize = 10*popsize;
  
  if ( initPopulation(&pop, sat) ) {
    fprintf(stderr,"Could not initialize potential.\n");
    return MEMORY_ERROR;
  }
  
  seed = time(0);
  //seed = 0; // for testing
  srand48(seed);
  printf("c ------------------------------------------------------\n");
  printf("c Substochastic Monte Carlo, version 1.0                \n");
  printf("c Brad Lackey, Stephen Jordan, and Michael Jarret, 2016.\n");
  printf("c ------------------------------------------------------\n");
  printf("c Input: %s\n", argv[1]);
  printf("c Bits: %d\n", nbts);
  printf("c Clauses (after tautology removal): %d\n", pop->sat->num_clauses);
  printf("c Step weight: %f\n", *weight);
  printf("c Population size: %d\n", popsize);
  printf("c Runtime: %e\n", *runtime);
  printf("c Start time: %e\n", *starttime);
  printf("c Trials: %i\n", *trials);
  printf("c Variables per clause: %f\n", varsperclause);
  printf("c Seed: %i\n", seed);
  
  *Pptr = pop;
  
  return 0;
}





void update(double a, double b, double mean, Population P, int parity){
  int i,j,k;
  double p,e;
  double min = P->sat->num_clauses;
  double avg = 0.0;
  double max = -(P->sat->num_clauses);
  int old = parity*arraysize;
  int new = (1-parity)*arraysize;
  
  for (i=j=0; i<P->psize; ++i) {   // Loop over each walker (i) and set target position (j) to zero.
    p = drand48();
    
    // First potential event: walker steps.
    if ( p < a ) {
      k = randomBitFlip(P->walker[new+j], P->walker[old+i]);
      e = P->walker[old+i]->potential + ((k>0)-(k<0))*getPotential(P->walker[new+j],P->ds->der[abs(k)-1]);
      P->walker[new+j]->potential = e;
      if ( e < min ){
        min = e;
        if ( e < P->winner->potential )
          copyBitstring(P->winner, P->walker[new+j]);
      }
      if ( e > max ) max = e;
      avg += e;
      ++j;
      continue;
    }
    p -= a;
    
    // Second potential event: walker spawns/dies.
    e = b*(mean - P->walker[old+i]->potential);
    if ( p < e ) {  // Spawn.
      copyBitstring(P->walker[new+j], P->walker[old+i]);
      copyBitstring(P->walker[new+j+1], P->walker[old+i]);
      e = P->walker[old+i]->potential;
      if ( e < min ) min = e;
      if ( e > max ) max = e;
      avg += 2*e;
      j += 2;
      continue;
    }
    
    if ( p < -e ) { // Die.
      continue;
    }
    
    // Third potential event: walker stays.
    copyBitstring(P->walker[new+j], P->walker[old+i]);
    e = P->walker[old+i]->potential;
    if ( e < min ) min = e;
    if ( e > max ) max = e;
    avg += e;
    ++j;
  }
  
  P->psize = j;
  P->avg_v = avg/j;
  P->min_v = min;
  P->max_v = max;
}
