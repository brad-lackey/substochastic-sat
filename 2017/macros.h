/** @file macros.h
 * @brief Header file defining all the macros.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 5/18/16.
 */

#ifndef _macros_h
#define _macros_h

#define GMP 0

#define MEMORY_ERROR 1         ///< Error code for failure to allocate memory.
#define IO_ERROR 2             ///< Error code for read/write failure (including formatting).
#define OUT_OF_BOUNDS_ERROR 4  ///< Error code for passing bad bounds to an array.

typedef unsigned long word_t;
typedef long potential_t;

#define VARIABLE_WORD_SIZE 1   ///< Number of bytes per word in generating bitstring lookup table.

#define CLAUSE_LIMB_BITS 64
#define CLAUSE_NUMB_BITS 64
#define CLAUSE_WORD_BITS 16
#define NUM_CLAUSE_WORDS (1<<CLAUSE_WORD_BITS)

#define VARIABLE_LIMB_BITS 64
#define VARIABLE_NUMB_BITS 64
#define VARIABLE_WORD_BITS 8
#define NUM_VARIABLE_WORDS (1<<VARIABLE_WORD_BITS)


//#define BYTES_PER_WORD sizeof(word_t)
//#define BITS_PER_WORD (8*sizeof(word_t))

#define TRACK_GLOBAL_BIASES 0   ///< Turn on/off biasing of individual bits based on problem structure/previous solutions.
#define GREEDY_DESCENT 0

#if TRACK_GLOBAL_BIASES
#define INITIAL_BUILD_RELAXATION 0.95  ///< Relaxation on the initial distribution based on examining the SAT instance.
#define UPDATE_RELAXATION        0.90  ///< Relaxation on the update distribution given by the previous winner(s).
#define REMIX_PERCENTAGE         0.60  ///< What percentage of the previous distribution is kept in the next step.
#endif

// Problem types.
#define UNKNOWN 0

#define UNWEIGHTED_2_SAT 10
#define UNWEIGHTED_3_SAT 11
#define UNWEIGHTED_4_SAT 12

#define PARTIAL_2_SAT 20
#define PARTIAL_3_SAT 21

#define WEIGHTED_2_SAT 30
#define WEIGHTED_3_SAT 31
#define WEIGHTED_4_SAT 32


#endif
