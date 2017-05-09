//
// Created by Bryan Luu on 5/3/2017.
//

#include "LUT.h"
#include "macros.h"

/**
 * @param fp file pointer to LUT
 * @param lut points to the LUT instance to be initialized.
 * @return Zero if successful, error code(s) if failed.
 */
int initLUT(FILE *fp, LUT *lut){
    char *line = NULL;
    size_t linecap = 0;
    ssize_t linelen;
    int nrows = 0;
    LUT _lut;

    // Read first line
    linelen = getline(&line, &linecap, fp);

    // Scan the line and if it does not have the right form: reject.
    if ( (sscanf(line, "%d", &nrows)) < 1 ){
        *lut = NULL;
        return IO_ERROR;
    }

    // Allocate memory for the LUT
    if ( (_lut = (LUT) malloc(sizeof(struct lut))) == NULL) {
        *lut = NULL;
        return MEMORY_ERROR;
    }

    _lut->nrows = nrows;

    if ( (_lut->times = (double*) malloc(nrows * sizeof(double))) == NULL) {
        *lut = NULL;
        return MEMORY_ERROR;
    }

    if ( (_lut->vals = (double*) malloc(nrows * sizeof(double))) == NULL) {
        freeLUT(&_lut);
        *lut = NULL;
        return MEMORY_ERROR;
    }

    double time;
    double val;
    int i=0;

    while( (linelen = getline(&line, &linecap, fp)) > 0)
    {
        if( sscanf(line, "%lf\t%lf", &time, &val) < 2 ){
            freeLUT(&_lut);
            *lut = NULL;
            return IO_ERROR;
        }

        _lut->times[i] = time;
        _lut->vals[i] = val;
        i++;
    }

    if(i != nrows){
        freeLUT(&_lut);
        return IO_ERROR;
    }

    *lut = _lut;

    return 0;
}

/**
 * @param lut points to the LUT instance to be destroyed.
 * @return None.
 */
void freeLUT(LUT *lut){
    if ( (*lut) != NULL ) {
        if ( (*lut)->times != NULL ){
            free((*lut)->times);
        }
        if ( (*lut)->vals != NULL ){
            free((*lut)->vals);
        }
        free(*lut);
        (*lut) = NULL;
    }
}