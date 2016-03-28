/** @file  sat.h
 * @brief Header file for a weighted SAT potential.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 3/14/16.
 */

#ifndef sat_h
#define sat_h

#include <stdio.h>
#include <stdlib.h>

#include "macros.h"
#include "bitstring.h"



/// The underlying type for a SAT instance
/**
 * This data type just holds the SAT clauses in DIMACS format.
 */
struct sat_st;
typedef struct sat_st * SAT;

struct sat_st {
  int num_clauses;        ///< Number of problems in the instance.
  int *clause_weight;     ///< Array holding the clause weights.
  int *clause_length;     ///< Array holding the length of each clause.
  int **clause;           ///< Array of clauses.
};

// I/O routines, since potentials are expected to be loaded from files.
int initSAT(SAT *sat_ptr, int ncls);
int loadDIMACSFile(FILE *fp, SAT *sat_ptr);   ///< Create a SAT instance from a file.
void freeSAT(SAT *sat_ptr);                   ///< Deallocation routine for a SAT instance.
void printSAT(FILE *fp, int nvars, SAT sat);  ///< Print in DIMACS format.

struct diff_sat_st;
typedef struct diff_sat_st * DSAT;

struct diff_sat_st {
  int num_vars;
  SAT *der;
};

int createSATDerivative(DSAT *dsat_ptr, int nvars, SAT sat);
void freeSATDerivative(DSAT *dsat_ptr);

#endif /* sat_h */
