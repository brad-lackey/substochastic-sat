/** @file  greedy.c
 * @brief Main file for a greedy descent.
 *
 * Created by Brad Lackey on 5/24/17. Last modified 5/24/17.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
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
extern potential_t topweight;

static int popsize,runmode;
static potential_t optimal;


int parseCommand(int argc, char **argv, Population *Pptr);
int descend(Population P);


int main(int argc, char **argv){
  int err;
  Population pop;
  potential_t min = -1;      //the best minimum from different trials
  Bitstring solution;        //the corresponding bitstring
  clock_t beg, end;          //for code timing
  double time_spent;         //for code timing
  
  beg = clock();
  
  if ( (err = parseCommand(argc, argv, &pop)) ){
    return err;
  }
  
  end = clock();
  time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
  printf("c Problem loaded: %f seconds\n", time_spent);
  
  if ( (err = initBitstring(&solution)) ){
    fprintf(stderr, "Could not initialize answerspace.\n");
    return err;
  }
  
  
  randomPopulation(pop,popsize);
  end = clock();
  min = pop->winner->potential;
  if ( min < topweight ) {
    printf("o %ld\n", min);
    fflush(stdout);
    printBits(stdout, pop->winner);
    fflush(stdout);
    time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
    printf("c Walltime: %f seconds, 0 loops\n", time_spent);
  }
  fflush(stdout);
  if (min <= optimal) exit(0);
  
  if ( (err = descend(pop)) ){
    return err;
  }
  
  end = clock();
  if ( pop->winner->potential < min ) {
    min = pop->winner->potential;
    if ( min < topweight ) {
      printf("o %ld\n", min);
      fflush(stdout);
      printBits(stdout, pop->winner);
      fflush(stdout);
      time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
      printf("c Walltime: %f seconds, 0 loops\n", time_spent);
    }
    fflush(stdout);
    if (min <= optimal) exit(0);
  }
  
  freeBitstring(&solution);
  
  return 0;
}


int parseCommand(int argc, char **argv, Population *Pptr) {
  SAT sat;
  int seed;
  FILE *fp;
  Population pop;
  
  if ( argc < 2 || argc > 4 ) {
    fprintf(stderr, "Usage: %s <instance.cnf> \n",argv[0]);
    fprintf(stderr, "Usage: %s <instance.cnf> [<target optimum> [<seed>]]\n",argv[0]);
    return 2;
  }
  
  printf("c ------------------------------------------------------\n");
  printf("c ----------------------------------------------------------\n");
  printf("c Substochastic Monte Carlo, version 1.1, 2017\n");
  printf("c Michael Jarret[1], Stephen Jordan[2], Brad Lackey[2], and Bryan Luu[3]\n");
  printf("c [1] Perimeter Institute for Theoretical Physics\n");
  printf("c [2] Joint Center for Quantum Information and Computer Science\n");
  printf("c       University of Maryland, College Park.\n");
  printf("c [3] University of British Columbia\n");
  printf("c ----------------------------------------------------------\n");
  printf("c Input: %s\n", argv[1]);
  
  if ( (fp = fopen(argv[1], "r")) == NULL ){
    fprintf(stderr,"Could not open file %s\n",argv[2]);
    return IO_ERROR;
  }
  
  if ( loadDIMACSFile(fp,&sat) ){
    fprintf(stderr,"Error reading in DIMACS SAT file %s\n",argv[2]);
    return IO_ERROR;
  }
  
  fclose(fp);
  
  setBitLength(sat->num_vars);
  
  printf("c Bits: %d\n", nbts);
  printf("c Clauses (after tautology removal): %d\n", sat->num_clauses);
  printf("c Problem type: %d\n", problem_type);
  
  if ( problem_type == UNKNOWN ) {
    popsize = 1024;
    if ( NUM_VARIABLE_WORDS*vlen*tlen*sizeof(word_t) + NUM_CLAUSE_WORDS*clen*sizeof(int) < (1<<30) )
      runmode = 2;
    else
      runmode = 1;
  }
  
  if ( problem_type == UNWEIGHTED_2_SAT ){
    popsize = 16;
    runmode = 1;
  }
  
  if (  problem_type == PARTIAL_2_SAT ){
    popsize = 16;
    runmode = 1;
  }
  
  if ( problem_type == WEIGHTED_2_SAT ){
    popsize = 16;
    runmode = 1;
  }
  
  if ( problem_type == UNWEIGHTED_3_SAT ) {
    popsize = 16;
    runmode = 1;
  }
  
  if ( problem_type == PARTIAL_3_SAT ) {
    popsize = 16;
    runmode = 1;
  }
  
  if ( problem_type == WEIGHTED_3_SAT ){
    popsize = 16;
    runmode = 1;
  }
  
  if ( (problem_type == UNWEIGHTED_4_SAT) || (problem_type == WEIGHTED_4_SAT) ) {
    popsize = 16;
    runmode = 2;
  }
  
  
  if ( argc >= 3 ) {
    
    optimal = atoi(argv[2]);
    
    if ( argc == 4 )
      seed = atoi(argv[3]);
    else
      seed = time(0);
    
  } else {
    
    optimal = 0;
    seed = time(0);
    
  }
  
  
  printf("c Target potential: %ld\n", optimal);
  printf("c Top potential: %ld\n", topweight);
  arraysize = 2000;
  
  srand48(seed);
  printf("c Seed: %i\n", seed);
  
  if ( initPopulation(&pop, sat, runmode) ) {
    fprintf(stderr,"Could not initialize potential.\n");
    return MEMORY_ERROR;
  }
  
  *Pptr = pop;
  
  return 0;
}


int descend(Population P){
  int i,j,k,argmin;
  int err;
  potential_t e,min;
  static Bitstring * stack;
  
  if ( (stack = (Bitstring *) malloc(P->sat->num_vars*sizeof(Bitstring))) == NULL ) return MEMORY_ERROR;
  for (i=0; i<P->sat->num_vars; ++i) {
    if ( (err = initBitstring(stack + i)) ) {
      return err;
    }
  }
  
  for (i=0; i<P->psize; ++i) {
    
    while(1){
      
      argmin = -1;
      min = P->walker[i]->potential;
      
      for (j=0; j<P->sat->num_vars; ++j) {
        k = flipBit(stack[j], P->walker[i], j);
        if ( P->ds != NULL )
          e = P->walker[i]->potential + ((k>0)-(k<0))*getPotential(stack[j],P->ds->der[j]);
        else {
          if ( P->tbl != NULL )
            e = getPotential2(stack[j],P->tbl);
          else
            e = getPotential(stack[j],P->sat);
        }
        
        if ( e < min ) {
          argmin = j;
          min = e;
          stack[j]->potential = e;
        }
      }
      
      if ( argmin >= 0 ) {
        copyBitstring(P->walker[i], stack[argmin]);
        if ( stack[argmin]->potential < P->winner->potential ) {
          copyBitstring(P->winner, stack[argmin]);
        }
      } else {
        break;
      }
    }
  }
  
  for (i=0; i<P->sat->num_vars; ++i)
    freeBitstring(stack + i);
  
  return 0;
}

