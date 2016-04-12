/** @file  sat.h
 * @brief Header file for a weighted SAT potential.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 4/7/16.
 */

#ifndef sat_h
#define sat_h

#include <stdio.h>
#include <stdlib.h>

#include "macros.h"
#include "bitstring.h"

#if GMP
#include "gmp.h"
#endif

int problem_type;
int clen;           /// The length of clause strings (in words).
int tlen;           /// The length of clause strings (in numbs).
int vlen;           /// The length of variable strings (in words).

/// The underlying type for a SAT instance.
/**
 * This data type just holds the SAT clauses in DIMACS format.
 */
struct sat_st;
typedef struct sat_st * SAT;

struct sat_st {
  int num_vars;           ///< Number of variables in the instance.
  int num_clauses;        ///< Number of problems in the instance.
  int *clause_weight;     ///< Array holding the clause weights.
  int *clause_length;     ///< Array holding the length of each clause.
  int **clause;           ///< Array of clauses.
#if TRACK_GLOBAL_BIASES
  double *global_bias;    ///< Array of biases on each bit.
#endif
};

// Constuctors and I/O routines. Potentials are expected to be loaded from files.
int initSAT(SAT *sat_ptr, int nvars, int ncls); ///< Initialize memory for a sat instance.
int loadDIMACSFile(FILE *fp, SAT *sat_ptr);     ///< Create a SAT instance from a file.
void freeSAT(SAT *sat_ptr);                     ///< Deallocation routine for a SAT instance.
void printSAT(FILE *fp, SAT sat);               ///< Print in DIMACS format.

int getPotential(Bitstring bts, SAT sat);       ///< Evaluates the SAT instance on the passed bitstring.


/// This is the underlying type for the derivative of a SAT instance.
/**
 * A crucial speedup in the algorithm is owed to not evaluating the SAT instance, but rather simply computing the difference based on the bit flipped.
 * This is really on efficient for very sparse SAT instances.
 */
struct diff_sat_st;
typedef struct diff_sat_st * DSAT;

struct diff_sat_st {
  int num_vars;
  SAT *der;
};

// Constuctors and I/O routines.
int createSATDerivative(DSAT *dsat_ptr, SAT sat);  ///< Creates and computes the derivative of the passed SAT instance.
void freeSATDerivative(DSAT *dsat_ptr);            ///< Deallocation routine for a SAT derivative.

struct incidence_table_st;
typedef struct incidence_table_st * Table;

struct incidence_table_st {
#if GMP
  int num_words;
  int num_bits;
  mpz_t *temp[1<<(8*VARIABLE_WORD_SIZE)];
  mpz_t buffer2;
#else
  word_t **incident[NUM_VARIABLE_WORDS];
  word_t *buffer;
#endif
  int *weight[NUM_CLAUSE_WORDS];
};

int createIncidenceTable(Table *t_ptr, SAT sat);
void freeIncidenceTable(Table *t_ptr);

int getPotential2(Bitstring bts, Table tbl);

#endif /* sat_h */
