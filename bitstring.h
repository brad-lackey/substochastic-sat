/** @file  bitstring.h
 * @brief Header file for a weighted bitstring.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 3/16/16.
 */

#ifndef bitstring_h
#define bitstring_h

#include <stdlib.h>
#include "macros.h"

/// The underlying type for a weighted bitstring.
/**
 * The basic container for a bitstring is an array of unsigned integers.
 * As the length of all bitstrings are the same, this is stored externally.
 * Additionally, a double (the weight) is stored.
 */
struct bitstring_st {
  word_t *node;                ///< The array that holds the bits.
  double potential;            ///< The weight associated to this bitstring.
};

int nbts;      ///< The length of bitstrings (in number of bits).
int blen;      ///< The length of bitstrings (in number of words).

typedef struct bitstring_st * Bitstring;

int setBitLength(int num_bits);                         ///< Sets the global variables appropriately.

int initBitstring(Bitstring *bst_ptr);                  ///< Routine for initializing a bitstring.
void freeBitstring(Bitstring *bst_ptr);                 ///< Routine for freeing a bitstring.
int randomBitstring(Bitstring bst);                     ///< Routine for filling with random bits.
int copyBitstring(Bitstring bst_out, Bitstring bst_in); ///< Routine for copying a bitstring.

int randomBitFlip(Bitstring bst_out, Bitstring bst_in); ///< Routine for stepping on the hypercube.

void printBits(Bitstring bst); //routine for printing a bitstring in MaxSAT 2015 format

#endif /* bitstring_h */
