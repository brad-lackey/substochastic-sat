/** @file macros.h
 * @brief Header file defining all the macros.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 3/16/16.
 */

#ifndef _macros_h
#define _macros_h

#define MEMORY_ERROR 1         ///< Error code for failure to allocate memory.
#define IO_ERROR 2             ///< Error code for read/write failure (including formatting).
#define OUT_OF_BOUNDS_ERROR 4  ///< Error code for passing bad bounds to an array.

typedef unsigned long word_t;

#define BYTES_PER_WORD sizeof(word_t)
#define BITS_PER_WORD (8*sizeof(word_t))
#endif
