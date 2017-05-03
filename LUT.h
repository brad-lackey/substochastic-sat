/** @file  LUT.h
 * @brief Header file for a LUT.
 *
 * Created by Bryan Luu on 5/3/17.
 */

#ifndef LUT_h
#define LUT_h

#include <stdio.h>

/// The underlying type for a LUT.
/**
 * This data type also holds a reference to the LUT instance.
 */
struct lut;
typedef struct lut * LUT;

struct lut {
    int *cols;     ///< Array for column values.
    int **rows;   ///< Array for row values.
};

// Memory management routines.
int initLUT(FILE *fp, LUT *tbl);  ///< Create a population from its SAT instance.
void freeLUT(LUT *tbl);                    ///< Deallocation routine for a LUT.

#endif //LUT_h