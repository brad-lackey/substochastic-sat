/** @file  LUT.h
 * @brief Header file for a LUT.
 *
 * Created by Bryan Luu on 5/3/17.
 */

#ifndef LUT_h
#define LUT_h

#include <stdio.h>
#include <stdlib.h>

/// The underlying type for a LUT.
/**
 * This data type also holds a reference to the LUT instance.
 */
struct lut;
typedef struct lut * LUT;

struct lut {
    int nrows;  ///< Number of rows.
    double total_time;
    double *times;   ///< Array of times.
    double *vals;    ///< Array of values.
};

// Memory management routines.
int initLUT(FILE *fp, LUT *lut);  ///< Create a population from its SAT instance.
void freeLUT(LUT *lut);                    ///< Deallocation routine for a LUT.

#endif //LUT_h