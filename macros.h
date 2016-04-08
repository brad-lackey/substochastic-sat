/** @file macros.h
 * @brief Header file defining all the macros.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 4/7/16.
 */

#ifndef _macros_h
#define _macros_h

#define GMP 0

#define MEMORY_ERROR 1         ///< Error code for failure to allocate memory.
#define IO_ERROR 2             ///< Error code for read/write failure (including formatting).
#define OUT_OF_BOUNDS_ERROR 4  ///< Error code for passing bad bounds to an array.

typedef unsigned long word_t;

#define CHUNK_SIZE 1           ///< Number of bytes used in lookup table.

#define BYTES_PER_WORD sizeof(word_t)
#define BITS_PER_WORD (8*sizeof(word_t))

#define TRACK_GLOBAL_BIASES 1   ///< Turn on/off biasing of individual bits based on problem structure/previous solutions.

#if TRACK_GLOBAL_BIASES
#define INITIAL_BUILD_RELAXATION 0.9   ///< Relaxation on the initial distribution based on examining the SAT instance.
#define UPDATE_RELAXATION        0.9   ///< Relaxation on the update distribution given by the previous winner(s).
#define REMIX_PERCENTAGE         0.95  ///< What percentage of the previous distribution is kept in the next step.
#endif

// Problem types.
#define UNKNOWN 0

#define UNWEIGHTED_2_SAT 10
#define UNWEIGHTED_3_SAT 11
#define UNWEIGHTED_4_SAT 12

#define WEIGHTED_2_SAT 30
#define WEIGHTED_3_SAT 31


#endif
