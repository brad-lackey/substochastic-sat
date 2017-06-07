/** @file  population.h
 * @brief Header file for a population.
 *
 * Created by Brad Lackey on 3/30/16. Last modified 6/7/17.
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
  Table tbl;
  int psize;           ///< The current population size.
  Bitstring *walker;   ///< Array of bitstrings that form the population.
  Bitstring winner;    ///< A copy of the best walker yet found.
  double avg_v;        ///< Average potential for the population.
  potential_t max_v;           ///< Maximum potential for the population.
  potential_t min_v;           ///< Minimum potential for the population.
};

// Memory management routines.
int initPopulation(Population *Pptr, SAT sat, int mode);  ///< Create a population from its SAT instance.
void freePopulation(Population *Pptr);                    ///< Deallocation routine for a population.
void randomPopulation(Population P, int size);            ///< Initialize a population of size \a size with random walkers.

// Adjust a population to a new size.
int reallocatePopulation(Population P, int new_size, int parity);

#endif /* population_h */
