/** @file  population.h
 * @brief Header file for a population.
 *
 * Created by Brad Lackey on 3/30/16. Last modified 3/30/16.
 */

#ifndef population_h
#define population_h

#include <stdio.h>
#include <stdlib.h>

#include "macros.h"
#include "bitstring.h"
#include "sat.h"

/// The integer \a arraysize is the amount of memory assigned to the population.
int arraysize;

/// The underlying type for a population.
/**
 * This data type also holds a reference to the SAT instance.
 */
struct population_st;
typedef struct population_st * Population;

struct population_st {
  SAT sat;             ///< Reference to the underlying SAT problem.
  DSAT ds;             ///< Reference to the SAT derivative for updates.
  int psize;           ///< The current population size.
  Bitstring *walker;   ///< Array of bitstrings that form the population.
  double avg_v;        ///< Average potential for the population.
  double max_v;        ///< Maximum potential for the population.
  double min_v;        ///< Minimum potential for the population.
};

// Memory management routines.
int initPopulation(Population *Pptr, FILE *fp);  ///< Create a population, and load its SAT instance from a file and compute the derivative.
void freePopulation(Population *Pptr);           ///< Deallocation routine for a population.
int randomPopulation(Population P, int size);    ///< Initialize a population of size \a size with random walkers.

#endif /* population_h */
