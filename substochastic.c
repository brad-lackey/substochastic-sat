/** @file  substochastic.c
 * @brief Main file for Substochastic Monte Carlo.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 4/7/16.
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

extern int blen;
extern int nbts;
extern int arraysize;
extern int problem_type;

static int popsize,runmode;
static double weight, runtime, runstep;
static int optimal;




void update(double a, double b, double mean, Population P, int parity);
int parseCommand(int argc, char **argv, Population *Pptr);


int main(int argc, char **argv){
  int parity, try, err;
  double mean;
  double a, b, t, dt;
  Population pop;
  double min = -1;      //the best minimum from different trials
  Bitstring solution;   //the corresponding bitstring
  clock_t beg, end;     //for code timing
  double time_spent;    //for code timing
  
#if TRACK_GLOBAL_BIASES
  int i;
  word_t u;
#endif
  
  
  beg = clock();
  
  if ( (err = parseCommand(argc, argv, &pop)) ){
    return err;
  }
  
  if ( (err = initBitstring(&solution)) ){
    fprintf(stderr, "Could not initialize answerspace.\n");
    return err;
  }
  
  
  try = 1;
  while (1) {
    
    t = 0.0;
    parity = 0;
    randomPopulation(pop,popsize);
    
    while (t < runtime) {
      
      // The annealing schedule
      a = weight*(1.0 - t/runtime); // Turned weight into percent -- Michael 3/30/16
      b = (t/runtime);
      
      mean = pop->avg_v + (pop->max_v - pop->min_v)*(popsize - pop->psize)/(2.0*popsize);
      
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
    
    end = clock();
    time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
    
    if ((min<0) || (pop->winner->potential < min)) {
      printBits(stdout, pop->winner);
      printf("c Walltime: %f seconds, %d loops\n", time_spent, try);
      fflush(stdout);
      min = pop->winner->potential;
      copyBitstring(solution, pop->winner);
      if (min == optimal) {
        break;
      }
    }
    if ( time_spent > 60 )
      break;
    
#if TRACK_GLOBAL_BIASES
    
    for (i=0; i<pop->sat->num_vars; ++i) {
      pop->sat->global_bias[i] *= REMIX_PERCENTAGE;
      u = pop->winner->node[i/BITS_PER_WORD] >> (i % BITS_PER_WORD);
      pop->sat->global_bias[i] += (1-REMIX_PERCENTAGE)*(UPDATE_RELAXATION - (2*UPDATE_RELAXATION -1)*(u%2));
    }
    
#endif
    
    runtime += runstep;
    try += 1;
  }
  
  freeBitstring(&solution);

  return 0;
}


int parseCommand(int argc, char **argv, Population *Pptr){
  SAT sat;
  int seed;
  FILE *fp;
  Population pop;
  
  if ( argc != 2 && argc != 3 && argc != 4 && argc != 6 ) {
    fprintf(stderr, "Usage: %s <instance.cnf> \n",argv[0]);
    fprintf(stderr, "Usage: %s <instance.cnf> <optimum>\n",argv[0]);
    fprintf(stderr, "Usage: %s <instance.cnf> <optimum> <seed>\n",argv[0]);
    fprintf(stderr, "Usage: %s <instance.cnf> <optimum> <seed> <step weight> <runtime>\n",argv[0]);
    return 2;
  }
  
  printf("c ------------------------------------------------------\n");
  printf("c Substochastic Monte Carlo, version 1.0                \n");
  printf("c Brad Lackey, Stephen Jordan, and Michael Jarret, 2016.\n");
  printf("c ------------------------------------------------------\n");
  printf("c Input: %s\n", argv[1]);
  
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
  
  printf("c Bits: %d\n", nbts);
  printf("c Clauses (after tautology removal): %d\n", sat->num_clauses);
  printf("c Problem type: %d\n", problem_type);

  if ( argc <= 4 ) {
    
    
    if ( (problem_type == UNKNOWN) || (problem_type == UNWEIGHTED_4_SAT) ) {
      weight = 0.50;
      runtime = 10000.0;
      runstep = 0.0;
      popsize = 16;
      runmode = 0;
    }
    
    if ( (problem_type == UNWEIGHTED_2_SAT) || (problem_type == WEIGHTED_2_SAT) ){
      weight = 0.18;
      runtime = exp(0.041*sat->num_vars + 3.5)/weight;
      runstep = exp(0.025*sat->num_vars + 2.5)/weight;
      popsize = 16;
      runmode = 1;
    }
    
    if ( (problem_type == UNWEIGHTED_3_SAT) || (problem_type == WEIGHTED_3_SAT) ) {
      weight = 0.16;
      runtime = exp(0.010*sat->num_vars + 8.9)/weight;
      runstep = exp(0.010*sat->num_vars + 5.8)/weight;
      popsize = 16;
      runmode = 1;
    }
    
    if ( problem_type == UNWEIGHTED_4_SAT ) {
      weight = 0.50;
      runtime = 1000000.0;
      runstep = 0.0;
      popsize = 16;
      runmode = 1;
    }
    
    
    if ( argc >= 3 ) {
      
      optimal = atoi(argv[2]);
      
      if ( argc == 4 ) {
        seed = atoi(argv[3]);
      } else {
        seed = time(0);
      }
      
      
    } else {
      
      optimal = 0;
      seed = time(0);
      
    }
    
  } else { // This is the timing trials case.
    
    optimal = atoi(argv[2]);
    seed = atoi(argv[3]);
    sscanf(argv[4], "%lf", &weight);
    sscanf(argv[5], "%lf", &runtime);
    weight /= 100.0;
    runstep = 0.0;
    popsize = 16;
    runmode = 1;
    
  }
  
  printf("c Population size: %d\n", popsize);
  printf("c Starting runtime: %.0f\n", runtime);
  printf("c Runtime step per loop: %.0f\n", runstep);
  printf("c Step weight: %.3f\n", weight);
  printf("c Target potential: %d\n", optimal);
  
  arraysize = 10*popsize;
  
  srand48(seed);
  printf("c Seed: %i\n", seed);

  if ( initPopulation(&pop, sat, runmode) ) {
    fprintf(stderr,"Could not initialize potential.\n");
    return MEMORY_ERROR;
  }
  
  *Pptr = pop;
  
  return 0;
}





void update(double a, double b, double mean, Population P, int parity){
  int i,j,k;
  double p,e;
  int min = P->sat->num_clauses;
  double avg = 0.0;
  int max = -(P->sat->num_clauses);
  int old = parity*arraysize;
  int new = (1-parity)*arraysize;
  
  for (i=j=0; i<P->psize; ++i) {   // Loop over each walker (i) and set target position (j) to zero.
    p = drand48();
    
    // First potential event: walker steps.
    if ( p < a ) {
      k = randomBitFlip(P->walker[new+j], P->walker[old+i]);
      if ( P->ds != NULL )
        e = P->walker[old+i]->potential + ((k>0)-(k<0))*getPotential(P->walker[new+j],P->ds->der[abs(k)-1]);
      else {
        if ( P->tbl != NULL )
          e = getPotential2(P->walker[new+j],P->tbl);
        else
          e = getPotential(P->walker[new+j],P->sat);
      }

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
