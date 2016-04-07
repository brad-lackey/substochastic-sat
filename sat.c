/** @file  sat.c
 * @brief Source file for a SAT potential type in the Substochastic library.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 4/7/16.
 */

#include <string.h>
#include <math.h>
#include "sat.h"


/**
 * A SAT instance is initialized.
 * Memory is allocated based on the number of clauses passed by \a ncls.
 * The number of variables is stored in the structure but not used here.
 * @param sat_ptr points to the SAT structure.
 * @param nvars is the number of variables in the instance.
 * @param ncls is the number of clauses in the instance.
 * @return Zero if successful, error code(s) if failed.
 */
int initSAT(SAT *sat_ptr, int nvars, int ncls){
  SAT sat;
  
  if ( (sat = (SAT) malloc(sizeof(struct sat_st))) == NULL ) {
    (*sat_ptr) = NULL;
    return MEMORY_ERROR;
  }
  
  sat->num_vars = nvars;
  sat->num_clauses = ncls;
  if ( (sat->clause_weight = (int *) calloc(ncls,sizeof(int))) == NULL ) {
    freeSAT(&sat);
    *sat_ptr = NULL;
    return MEMORY_ERROR;
  }
  if ( (sat->clause_length = (int *) calloc(ncls,sizeof(int))) == NULL ) {
    freeSAT(&sat);
    *sat_ptr = NULL;
    return MEMORY_ERROR;
  }
  if ( (sat->clause = (int **) calloc(ncls,sizeof(int *))) == NULL ) {
    freeSAT(&sat);
    *sat_ptr = NULL;
    return MEMORY_ERROR;
  }
  
#if TRACK_GLOBAL_BIASES
  if ( (sat->global_bias = (double *) calloc(nvars,sizeof(double))) == NULL ) {
    freeSAT(&sat);
    *sat_ptr = NULL;
    return MEMORY_ERROR;
  }
#endif

  (*sat_ptr) = sat;
  return 0;
}


// Comparison routine for quicksort.
int abscompare(const int *a, const int *b){
  int i = *a;
  int j = *b;
  return (abs(i) < abs(j)) ? -1 : 1;
}

// Hack job for deduping repeated variables in clauses using quicksort.
int dedupe(int *buffer, int size){
  int i,j;
  qsort(buffer,size,sizeof(int),abscompare);
  for (i=0,j=1; j<size; ++j) {
    if (buffer[i] == -buffer[j]) {
      return 0;
    }
    if (buffer[i] != buffer[j]) {
      buffer[i+1] = buffer[j];
      ++i;
    }
  }
  return i+1;
}

// Parse out the type of problem and parameters.
int parseHeader(char *line, int *nv, int *nc){
  int argc;
  char prob[50];
  int nvars, ncls, maxw;
  
  // Scan the line and if it does not have the right form: reject.
  if ( (argc = sscanf(line, "p %s %d %d %d", prob, &nvars, &ncls, &maxw)) < 3 ){
    return -1;
  }
  
  *nv = nvars;
  *nc = ncls;

  if ( strcmp(prob,"cnf") == 0 ){ // Then this is a unweighted max-sat instance.
    return 0;
  }
  
  if ( strcmp(prob,"wcnf") == 0 ){
    if ( argc == 3 ) { // Then this is a weighted max-sat instance.
      return 2;
    } else {           // Then this is a partial max-sat instance.
      return 1;
    }
  }
  
  return -1;
}


/**
 * This ((un)weighted/partial) SAT instance must be given in DIMACS-CNF format.
 * @param fp points to the file to be read.
 * @param sat_ptr points to the SAT instance to be created.
 * @return Number of variables if successful, zero if failed.
 */
int loadDIMACSFile(FILE *fp, SAT *sat_ptr){
  int i,j,k,w,off,type;
  int *buf;
  SAT sat = NULL;
  char *line = NULL;
  size_t linecap = 0;
  ssize_t linelen;
  int nvars,ncls;
  int avg_cls_length = 0;
#if TRACK_GLOBAL_BIASES
  int *total_weight;
#endif
  
  // First skip down until the parameter line is found.
  while ( (linelen = getline(&line, &linecap, fp)) > 0 ){
    if (line[0] != 'p') continue;
    
    // Read in the parameter line.
    if ( (type = parseHeader(line,&nvars,&ncls)) < 0 ) {
      *sat_ptr = NULL;
      return IO_ERROR;
    }
    
    // Now create the SAT instance with the given number of variables and clauses.
    if ( initSAT(&sat, nvars, ncls) ) {
      *sat_ptr = NULL;
      return MEMORY_ERROR;
    }
    
    // Create a large enough buffer to read in clause lines for the instance.
    if ( (buf = (int *) malloc(nvars*sizeof(int))) == NULL ) {
      *sat_ptr = NULL;
      return MEMORY_ERROR;
    }
    
#if TRACK_GLOBAL_BIASES
    if ( (total_weight = (int *) calloc(nvars,sizeof(int))) == NULL ) {
      free(buf);
      *sat_ptr = NULL;
      return MEMORY_ERROR;
    }
#endif
    
    
    break;
  }
  
  // Now start reading in clause lines.
  for (i=0; i<sat->num_clauses; ++i) {
    if ( (linelen = getline(&line, &linecap, fp)) <= 0 ) {
      freeSAT(&sat);
      free(buf);
#if TRACK_GLOBAL_BIASES
      free(total_weight);
#endif
      *sat_ptr = NULL;
      return IO_ERROR;
    }
    
    // Skip over a comment line. (Why have these in the body of the file anyway?)
    if (line[0] == 'c') {
      --i;
      continue;
    }
    
    
    if ( type == 0 ) {                // We can read from the beginning of the line.
      w = 1;
      off = 0;
    } else {                          // We need to read off the variable weight first.
      sscanf(line,"%d%n",&w,&off);
    }

    // Now read in the terms.
    for (j=0; j<nvars; ++j) {
      sscanf(line+off,"%d%n",buf+j,&k);
      off += k;
      if (buf[j] == 0) break;
    }
    
    // Some files were not preprocessed and so need tautology removal.
    j = dedupe(buf,j);
    
    
    if (j == 0) {                 // Then this clause is tautology and can be removed.
      --i;
      --(sat->num_clauses);
    } else {                      // Otherwise copy the buffer into the instance.
      sat->clause_weight[i] = w;
      sat->clause_length[i] = j;
      avg_cls_length += j;
      sat->clause[i] = (int *) malloc(j*sizeof(int));
      for (j=0; j<sat->clause_length[i]; ++j) {
        sat->clause[i][j] = buf[j];
#if TRACK_GLOBAL_BIASES
        if( buf[j] > 0 ){ // then this variable prefers to be true in this clause...
          sat->global_bias[buf[j]-1] += -w; // so it gets a penalty if it is set to false.
          total_weight[buf[j]-1] += abs(w);
        } else {          // then it would rather be false...
          sat->global_bias[(-buf[j])-1] += w; // so it gets a penalty if it is set to true.
          total_weight[(-buf[j])-1] += abs(w);
        }
#endif
      }
    }
  }
  
#if TRACK_GLOBAL_BIASES
  for (j=0; j<sat->num_vars; ++j) {
    if ( total_weight[j] > 0 ) {
      sat->global_bias[j] = (total_weight[j] + INITIAL_BUILD_RELAXATION*sat->global_bias[j])/(2*total_weight[j]);
    } else{
      sat->global_bias[j] = 0.5;
    }
  }
  free(total_weight);
#endif
  
  // Point to our newly minted instance and free the buffer memory.
  *sat_ptr = sat;
  free(buf);
  
  // Finalize our problem type.
  
  avg_cls_length = (int) round(((double) avg_cls_length)/sat->num_clauses);
  problem_type = UNKNOWN;
  
  if ( type == 0 ) { // This is an unweighted max-sat problem.
    
    if ( avg_cls_length <= 4 ) {
      problem_type = UNWEIGHTED_4_SAT;
    }

    if ( avg_cls_length <= 3 ) {
      problem_type = UNWEIGHTED_3_SAT;
    }

    if ( avg_cls_length <= 2) {
      problem_type = UNWEIGHTED_2_SAT;
    }
    
  }
  
  if ( type == 2 ){ // This is a weighted max-sat problem.
    
    if ( avg_cls_length <= 3) {
      problem_type = WEIGHTED_3_SAT;
    }

    if ( avg_cls_length <= 2) {
      problem_type = WEIGHTED_2_SAT;
    }
    
  }

  return 0;
}

/**
 * @param sat_ptr points to the SAT instance to be destroyed.
 * @return None.
 */
void freeSAT(SAT *sat_ptr){
  int j;
  if ( *sat_ptr != NULL ) {
    if ( (*sat_ptr)->clause != NULL ) {
      for (j=0; j<(*sat_ptr)->num_clauses; ++j)
        free((*sat_ptr)->clause[j]);
      free((*sat_ptr)->clause);
    }
    if ( (*sat_ptr)->clause_length != NULL ) {
      free((*sat_ptr)->clause_length);
    }
    if ( (*sat_ptr)->clause_weight != NULL ) {
      free((*sat_ptr)->clause_weight);
    }
#if TRACK_GLOBAL_BIASES
    if ( (*sat_ptr)->global_bias != NULL ) {
      free((*sat_ptr)->global_bias);
    }
#endif
    free(*sat_ptr);
    *sat_ptr = NULL;
  }
}

// Don't look at this.
void printSAT(FILE *fp, SAT sat){
  int i,j;
  
  fprintf(fp, "p wcnf %d %d\n", sat->num_vars, sat->num_clauses);
  for (i=0; i<sat->num_clauses; ++i) {
    fprintf(fp, "%d ",sat->clause_weight[i]);
    for (j=0; j<sat->clause_length[i]; ++j)
      fprintf(fp, "%d ", sat->clause[i][j]);
    fprintf(fp,"0 \n");
  }
}

/**
 * We need to return the potential as the sum of the weights of the failed clauses.
 * This is in the innermost loop of the algorithm so needs to be heavily optimized.
 * @param bst is the bitstring on which we compute the potential.
 * @param sat is the SAT instance that holds the potential.
 * @return The potential value.
 */
int getPotential(Bitstring bts, SAT sat){
  int i,j;
  int k,l;
  int p,q,v;
  word_t r;
  
  for (i=v=0; i<sat->num_clauses; ++i) {   // Loop through the clauses; set the output to zero.
    l = sat->clause_length[i];             // Store off the length of the i-th clause.
    for (j=0, k=1; j<l; ++j) {             // Loop though the terms in the i-th clause; set truth value to false. (Note: k is really ~k)
      p = sat->clause[i][j];               // Store off the j-th term of the i-th clause.
      q = abs(p)-1;                        // This the index of the variable in this term.
      r = bts->node[q/BITS_PER_WORD];      // Get the correct word from the bitstring.
      k &= (p > 0)^(r>>(q%BITS_PER_WORD)); // If the term is positive the value is kept, otherwise it is negated. (Note: k is really ~k)
    }
    k &= 1;
    v += k*sat->clause_weight[i];          // If the clause is _false_ then add in the weight as penalty. (Note: k is really ~k)
  }
  
  return v;
}

/**
 * A SAT derivative instance is initialized, and computed from the passed SAT instance.
 * The underlying SAT instance may or may not be weighted, but the derivative always is weighted.
 * @param dsat_ptr points to the SAT derivative structure to be created.
 * @param sat is the SAT instance to be differentiated.
 * @return Zero if successful, error code(s) if failed.
 */
int createSATDerivative(DSAT *dsat_ptr, SAT sat){
  int i,j,k,l;
  int *clen;
  SAT temp;
  DSAT dsat;
  
  if ( (clen = (int *) calloc(sat->num_vars,sizeof(int))) == NULL ) {
    (*dsat_ptr) = NULL;
    return MEMORY_ERROR;
  }
  if ( (dsat = (DSAT) malloc(sizeof(struct diff_sat_st))) == NULL) {
    free(clen);
    (*dsat_ptr) = NULL;
    return MEMORY_ERROR;
  }
  
  dsat->num_vars = sat->num_vars;
  if ( (dsat->der = (SAT *) malloc(sat->num_vars*sizeof(SAT))) == NULL) {
    free(clen);
    free(dsat);
    (*dsat_ptr) = NULL;
    return MEMORY_ERROR;
  }
  
  
  for (i=0; i<sat->num_clauses; ++i) {
    for (j=0; j<sat->clause_length[i]; ++j) {
      k = abs(sat->clause[i][j])-1;
      ++(clen[k]);
    }
  }
  
  for (i=0; i<sat->num_vars; ++i){
    initSAT(dsat->der + i, 0, clen[i]);
    clen[i] = 0;
  }
  
  for (i=0; i<sat->num_clauses; ++i) {
    for (j=0; j<sat->clause_length[i]; ++j) {
      if (sat->clause[i][j] < 0) {
        l = clen[-sat->clause[i][j]-1]++;
        temp = dsat->der[-sat->clause[i][j]-1];
        temp->clause_weight[l] = sat->clause_weight[i];
      } else {
        l = clen[sat->clause[i][j]-1]++;
        temp = dsat->der[sat->clause[i][j]-1];
        temp->clause_weight[l] = -(sat->clause_weight[i]);
      }
      
      temp->clause_length[l] = sat->clause_length[i]-1;
      temp->clause[l] = (int *) malloc(temp->clause_length[l]*sizeof(int));
      for (k=0; k<j; ++k)
        temp->clause[l][k] = sat->clause[i][k];
      for (; k<temp->clause_length[l]; ++k)
        temp->clause[l][k] = sat->clause[i][k+1];
    }
  }
  
  (*dsat_ptr) = dsat;
  free(clen);
  
  return 0;
}

/**
 * @param dsat_ptr points to the SAT derivative instance to be destroyed.
 * @return None.
 */
void freeSATDerivative(DSAT *dsat_ptr){
  int i;
  
  if ( (*dsat_ptr) != NULL ) {
    if ( (*dsat_ptr)->der != NULL ) {
      for (i=0; i<(*dsat_ptr)->num_vars; ++i)
        freeSAT((*dsat_ptr)->der + i);
      free((*dsat_ptr)->der);
    }
    free(*dsat_ptr);
    (*dsat_ptr) = NULL;
  }
}

