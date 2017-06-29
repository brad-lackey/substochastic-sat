/** @file  sat.h
 * @brief Header file for a weighted SAT potential.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 5/18/16.
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
potential_t topweight;

/// The underlying type for a SAT instance.
/**
 * This data type just holds the SAT clauses in DIMACS format.
 */
struct sat_st;
typedef struct sat_st * SAT;

struct sat_st {
  int num_vars;           ///< Number of variables in the instance.
  int num_clauses;        ///< Number of problems in the instance.
  int total_weight;
  potential_t *clause_weight;     ///< Array holding the clause weights.
  int *clause_length;     ///< Array holding the length of each clause.
  int **clause;           ///< Array of clauses.
#if TRACK_GLOBAL_BIASES
  double *global_bias;    ///< Array of biases on each bit.
#endif
};

struct map;
typedef struct map * SATMAP;

struct map {
  int num_vars;
  int * clause_map;
  int indices;
};

// Constuctors and I/O routines. Potentials are expected to be loaded from files.
int loadSATMAP(FILE * fp, SATMAP *map_ptr);     ///< Creates a SAT map from file.
int initSAT(SAT *sat_ptr, int nvars, int ncls); ///< Initialize memory for a sat instance.
int loadDIMACSFile(FILE *fp, SAT *sat_ptr);     ///< Create a SAT instance from a file.
int removeSoftClauses(SAT *sat_ptr, SAT *removed);     ///< Removes Soft Clauses from a SAT instance, placing them in removed.
void freeSAT(SAT *sat_ptr);                     ///< Deallocation routine for a SAT instance.
void printSAT(FILE *fp, SAT sat);               ///< Print in DIMACS format.
void printToCNF(FILE *fp, SAT * sat_ptr);               ///< Print in DIMACS format.
int recombineSAT(SAT * sat_one, SAT * sat_two, SAT * sum_ptr); ///< Combine SAT instances into a single instance.
void mapSAT(SAT *sat_ptr, SATMAP *map_ptr);


potential_t getPotential(Bitstring bts, SAT sat);       ///< Evaluates the SAT instance on the passed bitstring.


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
  potential_t *weight[NUM_CLAUSE_WORDS];
};

int createIncidenceTable(Table *t_ptr, SAT sat);
void freeIncidenceTable(Table *t_ptr);

potential_t getPotential2(Bitstring bts, Table tbl);

#endif /* sat_h */
