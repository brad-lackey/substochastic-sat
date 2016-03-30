/** @file  population.c
 * @brief Source file for a population.
 *
 * Created by Brad Lackey on 3/30/16. Last modified 3/30/16.
 */

#include <stdio.h>
#include <stdlib.h>
#include "population.h"



/**
 * This creates a population, reads in its associated SAT problem, and computes the derivative.
 * The problem instance must be a standard DIMACS formatted file, pointed to by \a fp.
 * The amount of memory assigned to the population is an global integer \a arraysize.
 * Note the population size is dynamics, and so enough memory must be assigned.
 * WARNING: there are currently no checks to ensure the population does not die off entirely, or overflow the array!
 * However in the current implementation of the \a update routine, it is highly unlikely that the population size grows/shrinks by more that a few percent of its initial size.
 * @param Pptr points to the population to be created.
 * @param fp points to the file to be read.
 * @return Zero if successful, error code if failed.
 */
int initPopulation(Population *Pptr, FILE *fp){
  int i,err;
  Population P;
  
  if ( (P = (Population) malloc(sizeof(struct population_st))) == NULL) {
    *Pptr = NULL;
    return MEMORY_ERROR;
  }
  
  if ( (i = loadDIMACSFile(fp,&(P->sat))) == 0 ){
    freePopulation(&P);
    *Pptr = NULL;
    return IO_ERROR;
  }
  createSATDerivative(&(P->ds),i,P->sat);
  setBitLength(i);
  
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
  P->max_v = 0.0;
  P->min_v = 0.0;
  
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
 * @return Zero.
 */
int randomPopulation(Population P, int size){
  int i, argmin;
  double e, avg, min, max;
  
  randomBitstring(P->walker[0]);
  e = getPotential(P->walker[0], P->sat);
  avg = e;
  min = e;
  argmin = 0;
  max = e;
  P->walker[0]->potential = e;
  
  for (i=1; i<size; ++i){
    randomBitstring(P->walker[i]);
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
  
  return 0;
}

