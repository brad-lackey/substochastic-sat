/** @file  sat.c
 * @brief Source file for a SAT potential type in the Substochastic library.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 5/18/16.
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
  if ( (sat->clause_weight = (potential_t *) calloc(ncls,sizeof(potential_t))) == NULL ) {
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
int parseHeader(char *line, int *nv, int *nc, potential_t *mw){
  int argc;
  char prob[50];
  int nvars, ncls;
  potential_t maxw;
  
  // Scan the line and if it does not have the right form: reject.
  if ( (argc = sscanf(line, "p %s %d %d %ld", prob, &nvars, &ncls, &maxw)) < 3 ){
    return -1;
  }
  
  *nv = nvars;
  *nc = ncls;
  
  if ( strcmp(prob,"cnf") == 0 ){ // Then this is a unweighted max-sat instance.
    *mw = 1;
    return 0;
  }
  
  if ( strcmp(prob,"wcnf") == 0 ){
    if ( argc == 3 ) { // Then this is a weighted max-sat instance.
      *mw = 10;
      return 2;
    } else {           // Then this is a partial max-sat instance.
      *mw = maxw;
      return 1;
    }
  }
  
  return -1;
}


int removeSoftClauses(SAT * sat_ptr, SAT * removed){
  SAT sat = *sat_ptr;
  SAT leftovers = NULL;
  int i, j;
  int numc = sat->num_clauses;

  for(i=0; i<sat->num_clauses; i++){
    // if soft clause
    if(sat->clause_weight[i] < numc){
      potential_t tmp_cls_weight = sat->clause_weight[i];
      int * tmp_cls = sat->clause[i];
      int tmp_cls_len = sat->clause_length[i];

      for(j=i; j<numc-1; j++){
        // remove j-th clause
        sat->clause_weight[j] = sat->clause_weight[j+1];
        sat->clause[j] = sat->clause[j+1];
        sat->clause_length[j] = sat->clause_length[j+1];
      }

      sat->clause_weight[sat->num_clauses - 1] = tmp_cls_weight;
      sat->clause[sat->num_clauses - 1] = tmp_cls;
      sat->clause_length[sat->num_clauses - 1] = tmp_cls_len;

      sat->num_clauses--;
      i--;
    }
  }

  // Now create the SAT instance with the given number of variables and clauses.
  if ( initSAT(&leftovers, sat->num_vars, numc-sat->num_clauses) ) {
    *removed = NULL;
    return MEMORY_ERROR;
  }

  for(i=0; i<leftovers->num_clauses; i++){
    leftovers->clause_length[i] = sat->clause_length[numc - 1 - i];
    leftovers->clause_weight[i] = sat->clause_weight[numc - 1 - i];
    leftovers->clause[i] = sat->clause[numc - 1 - i];
  }

  *removed = leftovers;

  return 0;
}


int recombineSAT(SAT * sat_one, SAT * sat_two, SAT * sum_ptr){
  SAT sat1 = *sat_one;
  SAT sat2 = *sat_two;
  SAT result = NULL;

  // Now create the SAT instance with the given number of variables and clauses.
  if ( initSAT(&result, sat1->num_vars, sat1->num_clauses + sat2->num_clauses) ) {
    *sum_ptr = NULL;
    return MEMORY_ERROR;
  }

  int i;

  for(i=0; i<sat1->num_clauses; i++){
    result->clause_weight[i] = sat1->clause_weight[i];
    result->clause_length[i] = sat1->clause_length[i];
    result->clause[i] = sat1->clause[i];
  }

  for(i=sat1->num_clauses; i<result->num_clauses; i++){
    result->clause_weight[i] = sat2->clause_weight[i-sat1->num_clauses];
    result->clause_length[i] = sat2->clause_length[i-sat1->num_clauses];
    result->clause[i] = sat2->clause[i-sat1->num_clauses];
  }


  *sum_ptr = result;

  return 0;
}


/**
 * This ((un)weighted/partial) SAT instance must be given in DIMACS-CNF format.
 * @param fp points to the file to be read.
 * @param sat_ptr points to the SAT instance to be created.
 * @return Number of variables if successful, zero if failed.
 */
int loadDIMACSFile(FILE *fp, SAT *sat_ptr){
  int i,j,k,off,type;
  int *buf;
  SAT sat = NULL;
  char *line = NULL;
  size_t linecap = 0;
  ssize_t linelen;
  int nvars,ncls;
  int max_cls_length = 0;
  int min_cls_length;
  potential_t weight_sum = 0;
  double avg_cls_length = 0.0;
  potential_t w,max_weight;
#if TRACK_GLOBAL_BIASES
  potential_t *total_weight;
#endif
  
  // First skip down until the parameter line is found.
  while ( (linelen = getline(&line, &linecap, fp)) > 0 ){
    if (line[0] != 'p') continue;
    
    // Read in the parameter line.
    if ( (type = parseHeader(line,&nvars,&ncls,&max_weight)) < 0 ) {
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
    if ( (total_weight = (potential_t *) calloc(nvars,sizeof(potential_t))) == NULL ) {
      free(buf);
      *sat_ptr = NULL;
      return MEMORY_ERROR;
    }
#endif
    
    min_cls_length = nvars;
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
      sscanf(line,"%ld%n",&w,&off);
    }
    
    // Now read in the terms.
    for (j=0; j<nvars; ++j) {
      sscanf(line+off,"%d%n",buf+j,&k);
      off += k;
      if (buf[j] == 0) break;
    }
    
    // Some files were not preprocessed and so need tautology removal.
    if ( j < min_cls_length )
      min_cls_length = j;
    j = dedupe(buf,j);
    
    
    if (j == 0) {                 // Then this clause is tautology and can be removed.
      --i;
      --(sat->num_clauses);
    } else {                      // Otherwise copy the buffer into the instance.
      sat->clause_weight[i] = w;
      weight_sum += w;
      sat->clause_length[i] = j;
      avg_cls_length += (double) j;
      if ( j > max_cls_length )
        max_cls_length = j;
      sat->clause[i] = (int *) malloc(j*sizeof(int));
      for (j=0; j<sat->clause_length[i]; ++j) {
        sat->clause[i][j] = buf[j];
        
#if TRACK_GLOBAL_BIASES
        if( buf[j] > 0 ){ // then this variable prefers to be true in this clause...
          sat->global_bias[buf[j]-1] += -w; // so it gets a penalty if it is set to false.
          total_weight[buf[j]-1] += labs(w);
        } else {          // then it would rather be false...
          sat->global_bias[(-buf[j])-1] += w; // so it gets a penalty if it is set to true.
          total_weight[(-buf[j])-1] += labs(w);
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
  sat->total_weight = weight_sum;
  *sat_ptr = sat;
  free(buf);
  
  // Finalize our problem type.
  
  problem_type = UNKNOWN;
  avg_cls_length /= ncls;
  
  if ( type == 0 ) { // This is an unweighted max-sat problem.
    
    topweight = weight_sum;
    
    if ( avg_cls_length < 4.01)
      problem_type = UNWEIGHTED_4_SAT;
    
    if ( avg_cls_length < 3.01)
      problem_type = UNWEIGHTED_3_SAT;
    
    if ( avg_cls_length < 2.01)
      problem_type = UNWEIGHTED_2_SAT;
    
  }
  
  if ( type == 1 ) { // This is a partial max-sat problem.
    
    topweight = max_weight;
    
    // Min-sat problems create some strange issues with average clause density.
    if ( avg_cls_length <= 3.1 ){
      problem_type = PARTIAL_3_SAT;
    }
    
    if ( avg_cls_length <= 2.1 ){
      problem_type = PARTIAL_2_SAT;
    }
    
  }
  
  if ( type == 2 ){ // This is a weighted max-sat problem.
    
    topweight = weight_sum;
    
    if ( max_cls_length <= 4) {
      problem_type = WEIGHTED_4_SAT;
    }
    
    if ( max_cls_length <= 3) {
      problem_type = WEIGHTED_3_SAT;
    }
    
    if ( max_cls_length <= 2) {
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
    fprintf(fp, "%ld ",sat->clause_weight[i]);
    for (j=0; j<sat->clause_length[i]; ++j)
      fprintf(fp, "%d ", sat->clause[i][j]);
    fprintf(fp,"0 \n");
  }
}


void printToCNF(FILE *fp, SAT * sat_ptr){
  int i,j;
  SAT sat = *sat_ptr;

  fprintf(fp, "p cnf %d %d\n", sat->num_vars, sat->num_clauses);
  for (i=0; i<sat->num_clauses; ++i) {
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
potential_t getPotential(Bitstring bts, SAT sat){
  int i,j;
  int k,l;
  int p,q;
  word_t r;
  potential_t v;
  
  for (i=v=0; i<sat->num_clauses; ++i) {        // Loop through the clauses; set the output to zero.
    l = sat->clause_length[i];                  // Store off the length of the i-th clause.
    for (j=0, k=1; j<l; ++j) {                  // Loop though the terms in the i-th clause; set truth value to false. (Note: k is really ~k)
      p = sat->clause[i][j];                    // Store off the j-th term of the i-th clause.
      q = abs(p)-1;                             // This the index of the variable in this term.
      r = bts->node[q/VARIABLE_NUMB_BITS];      // Get the correct word from the bitstring.
      k &= (p > 0)^(r>>(q%VARIABLE_NUMB_BITS)); // If the term is positive the value is kept, otherwise it is negated. (Note: k is really ~k)
    }
    k &= 1;
    v += k*sat->clause_weight[i];               // If the clause is _false_ then add in the weight as penalty. (Note: k is really ~k)
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
/// Need to remove clen if we want to do both derivatives and incidence tables.
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


int createIncidenceTable(Table *t_ptr, SAT sat){
  int i,j,k,l;
  Table tbl;
  
  if ( (tbl = (Table) malloc(sizeof(struct incidence_table_st))) == NULL ) {
    *t_ptr = NULL;
    return MEMORY_ERROR;
  }
  
  clen = (sat->num_clauses-1)/CLAUSE_WORD_BITS + 1;
  tlen =(sat->num_clauses-1)/CLAUSE_NUMB_BITS + 1;
  vlen = (sat->num_vars-1)/VARIABLE_WORD_BITS + 1;
//  printf("c Number of clauses per word: %d\n",NUM_CLAUSE_WORDS);
//  printf("c Number of clause words: %d\n",clen);
//  printf("c Number of variables per word: %d\n",NUM_VARIABLE_WORDS);
//  printf("c Number of variable words: %d\n",vlen);
  
  
#if GMP
  tbl->num_bits = sat->num_clauses;
  for (i=0; i<(1<<(8*VARIABLE_WORD_SIZE)); ++i) {
    if ( (tbl->temp[i] = (mpz_t *) malloc(tbl->num_words*sizeof(mpz_t))) == NULL ) {
      freeIncidenceTable(&tbl);
      *t_ptr = NULL;
      return MEMORY_ERROR;
    }
    for (j=0; j<tbl->num_words; ++j)
      mpz_init2(tbl->temp[i][j],sat->num_clauses);
  }
#else
  for (i=0; i<NUM_VARIABLE_WORDS; ++i) {
    if ( (tbl->incident[i] = (word_t **) malloc(vlen*sizeof(word_t *))) == NULL ) {
      freeIncidenceTable(&tbl);
      *t_ptr = NULL;
      return MEMORY_ERROR;
    }
    for (j=0; j<vlen; ++j){
      if ( (tbl->incident[i][j] = (word_t *) calloc(tlen, sizeof(word_t))) == NULL ) {
        freeIncidenceTable(&tbl);
        *t_ptr = NULL;
        return MEMORY_ERROR;
      }
    }
  }
//  printf("c Incidence table size: %lu bytes\n",NUM_VARIABLE_WORDS*vlen*tlen*sizeof(word_t));
#endif
  
  for (i=0; i<NUM_CLAUSE_WORDS; ++i) {
    if ( (tbl->weight[i] = (potential_t *) calloc(clen, sizeof(potential_t))) == NULL ) {
      freeIncidenceTable(&tbl);
      *t_ptr = NULL;
      return MEMORY_ERROR;
    }
  }
//  printf("c Weight table size: %lu bytes\n",NUM_CLAUSE_WORDS*clen*sizeof(int));
  
#if GMP
  tbl->num_bits = sat->num_clauses;
  for (i=0; i<(1<<(8*VARIABLE_WORD_SIZE)); ++i) {
    if ( (tbl->temp[i] = (mpz_t *) malloc(tbl->num_words*sizeof(mpz_t))) == NULL ) {
      freeIncidenceTable(&tbl);
      *t_ptr = NULL;
      return MEMORY_ERROR;
    }
    for (j=0; j<tbl->num_words; ++j)
      mpz_init2(tbl->temp[i][j],sat->num_clauses);
  }
#else
  if ( (tbl->buffer = (word_t *) malloc(tlen*sizeof(word_t))) == NULL ) {
    freeIncidenceTable(&tbl);
    *t_ptr = NULL;
    return MEMORY_ERROR;
  }
#endif
  
  for (i=0; i<sat->num_clauses; ++i) {                         // Loop over clauses.
    
    // ------------- Build incidence table -----------
    for (j=0; j<sat->clause_length[i]; ++j) {                  // Loop over variables in each clause.
      k = sat->clause[i][j];                                   // Let k be the current variable (1-up).
      if ( k > 0 ) {                                           // Then this variable is not negated.
        for (l=0; l<NUM_VARIABLE_WORDS; ++l)                   // Loop over variable set patterns.
          if (  ((l>>((k-1)%VARIABLE_WORD_BITS))&1) == 1 ) {   // If in this pattern, this variable is set to true...
#if GMP
            mpz_setbit(tbl->temp[l][(k-1)/(8*VARIABLE_WORD_SIZE)],i);
#else                                                          // then mark the i-th bit in the incidence table as true.
            tbl->incident[l][(k-1)/VARIABLE_WORD_BITS][i/CLAUSE_NUMB_BITS] |= ((word_t) 1) << (i%CLAUSE_NUMB_BITS);
#endif
          }
      } else {                                                 // Then this variable is negated.
        for (l=0; l<NUM_VARIABLE_WORDS; ++l)                   // Loop over variable set patterns.
          if (  ((l>>((-k-1)%VARIABLE_WORD_BITS))&1) == 0 ) {  // If in this pattern, this variable is set to false...
#if GMP
            mpz_setbit(tbl->temp[l][(-k-1)/(8*VARIABLE_WORD_SIZE)],i);
#else                                                          // then mark the i-th bit in the incidence table as true.
            tbl->incident[l][(-k-1)/VARIABLE_WORD_BITS][i/CLAUSE_NUMB_BITS] |= ((word_t) 1) << (i%CLAUSE_NUMB_BITS);
#endif
          }
      }
    }
    
    // ------------- Build weight table -----------
    for (j=0; j<NUM_CLAUSE_WORDS; ++j) {                             // Loop over clause patterns.
      if ( ((j>>(i%CLAUSE_WORD_BITS)) & 1) == 0 ) {                  // If this clause is false...
        tbl->weight[j][i/CLAUSE_WORD_BITS] += sat->clause_weight[i]; // then add its weight to the table.
      }
    }
    
  }
  
  *t_ptr = tbl;
  return 0;
}

void freeIncidenceTable(Table *t_ptr){
  int j,k;
  if ( (*t_ptr) != NULL ) {
    
    for (j=0; j<VARIABLE_WORD_SIZE; ++j) {
      
#if GMP
      if ( (*t_ptr)->temp[j] != NULL ){
        for (k=0; k<(*t_ptr)->num_words; ++k)
          mpz_clear((*t_ptr)->temp[j][k]);
        free((*t_ptr)->temp[j]);
      }
#else
      if ( (*t_ptr)->incident[j] != NULL ){
        for (k=0; k<vlen; ++k){
          if ( (*t_ptr)->incident[j][k] != NULL)
            free((*t_ptr)->incident[j][k]);
        }
        free((*t_ptr)->incident[j]);
      }
#endif
    }
    
#if GMP
    mpz_clear((*t_ptr)->buffer2);
#else
    if ( (*t_ptr)->buffer != NULL )
      free((*t_ptr)->buffer);
#endif
    free(*t_ptr);
    
    for (j=0; j<NUM_CLAUSE_WORDS; ++j)
      if ( (*t_ptr)->weight[j] != NULL )
        free((*t_ptr)->weight[j]);
  }
  (*t_ptr) = NULL;
}


potential_t getPotential2(Bitstring bts, Table tbl){
  int i,j,k;
  potential_t w;
  
  
#if GMP
  mpz_set_ui(tbl->buffer2,0);
#else
  for (j=0; j<tlen; ++j) tbl->buffer[j] = 0;
#endif
  
  for (i=0; i<vlen; ++i) {
    j = VARIABLE_WORD_BITS*i;
    k = (bts->node[j/VARIABLE_NUMB_BITS] >> (j%VARIABLE_NUMB_BITS)) % NUM_VARIABLE_WORDS;
    
#if GMP
    mpz_ior(tbl->buffer2,tbl->temp[k][i],tbl->buffer2);
#else
    
    for (j=0; j<tlen; ++j)
      tbl->buffer[j] |= tbl->incident[k][i][j];
#endif
  }
  
#if GMP
  for (i=j=0; i<tbl->num_bits; ++i){
    if ( mpz_tstbit(tbl->buffer2,i) == 0 )
      j += tbl->weight[i];
  }
#else
  for (i=w=0; i<clen; ++i) {
    j = CLAUSE_WORD_BITS*i;
    k = (tbl->buffer[j/CLAUSE_NUMB_BITS] >> (j%CLAUSE_NUMB_BITS)) % NUM_CLAUSE_WORDS;
    w += tbl->weight[k][i];
  }
#endif
  
  return w;
}

