/** @file  population.c
 * @brief Source file for a population.
 *
 * Created by Brad Lackey on 3/30/16. Last modified 6/7/17.
 */

#include <stdio.h>
#include <stdlib.h>
#include "population.h"



/**
 * This creates a population from the passed SAT instance.
 * If the low bit of parameter \a mode is set, it also computes and stores the derivative of that SAT instance.
 * If the next low bit of \a mode is set, it also computes a look table for fast evaluation, and stores it.
 * The amount of memory assigned to the population is an global integer \a arraysize.
 * Note the population size is dynamics, and so enough memory must be assigned.
 * WARNING: there are currently no checks to ensure the population does not die off entirely, or overflow the array!
 * However in the current implementation of the \a update routine, it is highly unlikely that the population size grows/shrinks by more that a few percent of its initial size.
 * @param Pptr points to the population to be created.
 * @param sat is the underlying SAT instance.
 * @param mode indicated what additional structures to create.
 * @return Zero if successful, error code(s) if failed.
 */
int initPopulation(Population *Pptr, SAT sat, int mode){
  int i,err;
  Population P;
  
  if ( (P = (Population) malloc(sizeof(struct population_st))) == NULL) {
    *Pptr = NULL;
    return MEMORY_ERROR;
  }
  
  P->sat = sat;
  
  if ( (mode & 1) )
    createSATDerivative(&(P->ds),sat);
  else
    P->ds = NULL;
  
  if ( (mode & 2) ){
    createIncidenceTable(&(P->tbl),sat);
  }
  else
    P->tbl = NULL;

  P->psize = 0;

  if ( (P->walker = (Bitstring *) malloc((2*arraysize)*sizeof(Bitstring))) == NULL ){
    freePopulation(&P);
    *Pptr = NULL;
    return MEMORY_ERROR;
  }
  
  if ( (err = initBitstring(&(P->winner))) ) {
    freePopulation(&P);
    *Pptr = NULL;
    return err;
  }
  
  for (i=0; i<(2*arraysize); ++i) {
    if ( (err = initBitstring(P->walker + i)) ) {
      freePopulation(&P);
      *Pptr = NULL;
      return err;
    }
  }
  
  P->avg_v = 0.0;
  P->max_v = 0;
  P->min_v = 0;
  
  *Pptr = P;
  return 0;
}


/**
 * @param Pptr points to the population instance to be destroyed.
 * @return None.
 */
void freePopulation(Population *Pptr){
  int j;
  if ( *Pptr != NULL ) {
    if ( (*Pptr)->sat != NULL )
      freeSAT(&((*Pptr)->sat));
    if ( (*Pptr)->ds != NULL )
      freeSATDerivative(&((*Pptr)->ds));
    if ( (*Pptr)->winner != NULL )
      freeBitstring(&((*Pptr)->winner));
    if ( (*Pptr)->walker != NULL ) {
      for (j=0; j<(2*arraysize); ++j)
        freeBitstring((*Pptr)->walker + j);
      free((*Pptr)->walker);
    }
    *Pptr = NULL;
  }
}


/**
 * This routine initializes an already created population of specified size with random walkers.
 * The potential of each walker is computed and stored, so updates can be performed using derivatives.
 * The walker with the best potential is stored off as the winner.
 * Note: this initial distribution is uniform, which is the ground state for the hypercube Laplacian.
 * If one changes the driving Hamiltonian, then this routine will likely need to be modified.
 * @param P is the population to be initialized.
 * @param size is the initial size of the population.
 * @return None.
 */
void randomPopulation(Population P, int size){
  int i, argmin;
  potential_t e, min, max;
  double avg;
  
#if TRACK_GLOBAL_BIASES
  int j;
  word_t u;
  
  for (j=0; j<blen; ++j) {
    P->walker[0]->node[j] = (word_t) 0;
  }
  for (j=0; j<nbts; ++j) {
    u = (word_t) (drand48() > P->sat->global_bias[j]);
    P->walker[0]->node[j/VARIABLE_NUMB_BITS] ^= u << (j % VARIABLE_NUMB_BITS);
  }
#else
  randomBitstring(P->walker[0]);
#endif

  if ( P->tbl != NULL )
    e = getPotential2(P->walker[0],P->tbl);
  else
    e = getPotential(P->walker[0], P->sat);

  avg = e;
  min = e;
  argmin = 0;
  max = e;
  P->walker[0]->potential = e;
  
  for (i=1; i<size; ++i){
    
#if TRACK_GLOBAL_BIASES
    for (j=0; j<blen; ++j) {
      P->walker[i]->node[j] = (word_t) 0;
    }
    for (j=0; j<nbts; ++j) {
      u = (word_t) (drand48() > P->sat->global_bias[j]);
      P->walker[i]->node[j/VARIABLE_NUMB_BITS] ^= u << (j % VARIABLE_NUMB_BITS);
    }
#else
    randomBitstring(P->walker[i]);
#endif
    
    if ( P->tbl != NULL )
      e = getPotential2(P->walker[i],P->tbl);
    else
      e = getPotential(P->walker[i], P->sat);
    avg += e;
    if ( e < min ){
      min = e;
      argmin = i;
    }
    if ( e > max )
      max = e;
    P->walker[i]->potential = e;
  }
  
  P->psize = size;
  P->avg_v = avg/size;
  P->max_v = max;
  P->min_v = min;
  
  copyBitstring(P->winner, P->walker[argmin]);
}


int reallocatePopulation(Population P, int new_size, int parity){
  int i,j,err,new_arraysize;
  
  // First check if we need more memory (and move walkers into correct position if needed).
  if ( arraysize < 3*new_size ) {
    
    // Allocate memory.
    new_arraysize = 3*new_size;
    if ( (P->walker = (Bitstring *) realloc(P->walker, (2*new_arraysize)*sizeof(Bitstring))) == NULL ){
      freePopulation(&P);
      return MEMORY_ERROR;
    }
    
    // Initialize new walker locations.
    for (i=2*arraysize; i<(2*new_arraysize); ++i) {
      if ( (err = initBitstring(P->walker + i)) ) {
        freePopulation(&P);
        return err;
      }
    }
    
    // Move walkers to new position in array if necessary.
    if ( parity )
      for (i=0; i<P->psize; ++i)
        copyBitstring(P->walker[new_arraysize+i], P->walker[arraysize+i]);
    
    // Adjust arraysize
    arraysize = new_arraysize;
  }
  
  // Now we can grow/shrink the population to correct size.
  err = parity*arraysize;                              // Reuse the err variable as the parity offset.
  for (i=P->psize; i<new_size; ++i) {                  // Add walkers if needed.
    j = lrand48() % P->psize;                          // Choose a random walker.
    copyBitstring(P->walker[err+i], P->walker[err+j]); // Clone it into the next position.
  }
  P->psize = new_size;
  
  return 0;
}

