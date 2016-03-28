#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "macros.h"
#include "bitstring.h"
#include "sat.h"

//On machines with very old versions of glibc (e.g. the Raritan cluster)
//we need to define gnu_source in order to avoid warnings about implicit
//declaration of getline() and round().
#define _GNU_SOURCE

extern int nbts;
static int arraysize;

struct population_st {
  SAT sat;
  DSAT ds;
  int psize;
  Bitstring *walker;
  double avg_v;
  double max_v;
  double min_v;
};

typedef struct population_st * Population;

int initPopulation(Population *Pptr, FILE *fp);
void freePopulation(Population *Pptr);
int randomPopulation(Population P, int size);

double getPotential(Bitstring bts, SAT sat);
void update(double a, double b, double mean, Population P, int parity);

/*Added by SPJ 3/17
 *Print optimal bitstring out of the current population. If there are multiple
 *optima with the same potential, just print one of them. Returns the index of
 *the optimum printed.
 */
int printOpt(Population P) {
  int i;
  double potential;
  i = -1;
  do {
    i++;
    potential = getPotential(P->walker[i], P->sat);
    if(potential == P->min_v) printBits(P->walker[i]);
  }while(potential != P->min_v);
  return i;
}

/*Added by SPJ 3/17
 *Returns the average number of variables per clause.
 */
double avgLength(SAT instance) {
  int c;
  int total;
  total = 0;
  for(c = 0; c < instance->num_clauses; c++) total += instance->clause_length[c];
  return (double)total/(double)instance->num_clauses;
}

/*Added by SPJ 3/17
 *Automatically chooses parameters, such as runtime and number of trials.
 */
void autoparam(double varsperclause, int vars, double *weight, double *runtime, int *trials) {
  int k;
  k = (int)round(varsperclause);
  *weight = 1.0;
  *trials = 5; //default
  if(k == 2) {
    //These values are determined using the experimental results in
    //summary_2sat120v.txt and summary_2sat200v.txt.
    *trials = 10;
    *runtime = 3.2E4;
  }
  if(k == 3) {
    //These values are determined using the experimental results in
    //summary_3sat70v.txt and summary3sat110vlong.txt.
    *trials = 5;
    *runtime = 223.0*exp(0.07*vars);
  }
  if(k == 4 || vars > 200) {
    //This overrides the above.
    //Just go for broke and use up all the time.
    *trials = 5;
    *runtime = 3E6;
  }
}

int main(int argc, char **argv){
  int i, parity, trials, try, popsize;
  double weight, runtime, mean, mark;
  double a, b, t, dt;
  FILE *fp;
  Population pop;
  int optindex;         //index of optimum found
  int *mins;            //the minima from different trials
  Bitstring *solutions; //the corresponding bitstrings
  int best_try;         //the trial in which the best minimum was found
  int initfail;         //a flag for memory allocation failure
  clock_t beg, end;     //for code timing
  double time_spent;    //for code timing
  double varsperclause; //average number of variables per clause
  long seed;            //the seed for the RNG
  beg = clock();
  if (argc != 6 && argc != 2) {
    printf("Usage: %s instance.cnf <step weight> <runtime> <population size> <trials>\n",argv[0]);
    return 2;
  }
  if ( (fp = fopen(argv[1], "r")) == NULL ){
    fprintf(stderr,"Could not open file %s\n",argv[1]);
    return IO_ERROR;
  }
  //If we want to autoparam popsize we will have to separate loading the DIMACS
  //file and allocating the population memory into separate steps. Right now they
  //are done together in initPopulation.
  if(argc == 6) {
    weight = (double) atoi(argv[2]);
    sscanf(argv[3], "%lf", &runtime); //allows scientific notation unlike atoi
    popsize = atoi(argv[4]);
    trials = atoi(argv[5]);
  }
  else popsize = 128;
  arraysize = 2*popsize;
  if ( initPopulation(&pop, fp) ) {
    fprintf(stderr,"Could not initialize potential.\n");
    return MEMORY_ERROR;
  }
  fclose(fp);
  varsperclause = avgLength(pop->sat);
  seed = time(0);
  //seed = 0 // for testing
  srand48(seed);
  if(argc == 2) autoparam(varsperclause, nbts, &weight, &runtime, &trials);
  printf("c ------------------------------------------------------\n");
  printf("c Substochastic Monte Carlo, version 1.0                \n");
  printf("c Brad Lackey, Stephen Jordan, and Michael Jarret, 2016.\n");
  printf("c ------------------------------------------------------\n");
  printf("c Input: %s\n", argv[1]); 
  printf("c Bits: %d\n", nbts);
  printf("c Clauses (after tautology removal): %d\n", pop->sat->num_clauses);
  printf("c Step weight: %f\n", weight);
  printf("c Population size: %d\n", popsize);
  printf("c Runtime: %e\n", runtime);
  printf("c Trials: %i\n", trials);
  printf("c Variables per clause: %f\n", varsperclause);
  printf("c Seed: %li\n", seed);
  mark = runtime/100.0;
  mins = (int *)malloc(trials*sizeof(int));
  solutions = (Bitstring *)malloc(trials*sizeof(Bitstring));
  if(mins == NULL || solutions == NULL) {
    fprintf(stderr, "Could not initialize trials.\n");
    return MEMORY_ERROR;
  }
  initfail = 0;
  for(try = 0; try < trials && !initfail; try++) initfail = initBitstring(&solutions[try]);
  if(initfail) {
    fprintf(stderr, "Could not initialize answerspace.\n");
    return MEMORY_ERROR;
  }
  for(try = 0; try < trials; try++) {
    t = 0.0;
    parity = 0;
    randomPopulation(pop,arraysize/2);
    while (t < runtime) {
      a = weight*(1.0 - t/runtime);
      b = (t/runtime);
      mean = pop->avg_v + ((pop->max_v - pop->min_v)/arraysize)*((arraysize/2) - pop->psize);
      if ( (pop->max_v - mean) > (mean - pop->min_v) )
        dt = 0.9/(a + b*(pop->max_v - mean));
      else
        dt = 0.9/(a + b*(mean - pop->min_v));
      if (t + dt > runtime)
        dt = runtime - t;      
      update(a*dt, b*dt, mean, pop, parity);      
      t += dt;
      parity ^= 1;
      if ( (trials==0) && (t >= mark) ){
        fprintf(stderr,"time=%f: size=%3d    min=%f   max=%f\n", t, pop->psize, pop->min_v, pop->max_v);
        mark += runtime/100.0;
      }
    }
    printf("o %i\n",(int)pop->min_v);
    printf("v ");
    optindex = printOpt(pop);
    mins[try] = (int)pop->min_v;
    copyBitstring(solutions[try], pop->walker[optindex]);
  }
  best_try = 0;
  for(try = 1; try < trials; try++) if(mins[try] < mins[best_try]) best_try = try;
  printf("c Final answer: \n");
  printf("o %i\n", mins[best_try]);
  printf("v ");
  printBits(solutions[best_try]);
  for(try = 0; try < trials; try++) freeBitstring(&solutions[try]);
  free(solutions);
  free(mins);
  end = clock();
  time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
  printf("c Walltime: %f seconds\n", time_spent);
  return 0;
}


int initPopulation(Population *Pptr, FILE *fp){
  int i,err;
  Population P;

  if ( (P = (Population) malloc(sizeof(struct population_st))) == NULL) {
    *Pptr = NULL;
    return MEMORY_ERROR;
  }

  if ( (i = loadDIMACSFile(fp,&(P->sat))) == 0 ){
    freePopulation(&P);
    *Pptr = NULL;
    return IO_ERROR;
  }
  createSATDerivative(&(P->ds),i,P->sat);
  setBitLength(i);

  P->psize = 0;
  
  if ( (P->walker = (Bitstring *) malloc((2*arraysize)*sizeof(Bitstring))) == NULL ){
    freePopulation(&P);
    *Pptr = NULL;
    return MEMORY_ERROR;
  }

  for (i=0; i<(2*arraysize); ++i) {
    if ( (err = initBitstring(P->walker + i)) ) {
      freePopulation(&P);
      *Pptr = NULL;
      return err;
    }
  }

  P->avg_v = 0.0;
  P->max_v = 0.0;
  P->min_v = 0.0;
  
  *Pptr = P;
  return 0;
}

void freePopulation(Population *Pptr){
  int j;
  if ( *Pptr != NULL ) {
    if ( (*Pptr)->sat != NULL )
      freeSAT(&((*Pptr)->sat));
    if ( (*Pptr)->ds != NULL )
      freeSATDerivative(&((*Pptr)->ds));
    if ( (*Pptr)->walker != NULL ) {
      for (j=0; j<(2*arraysize); ++j)
        freeBitstring((*Pptr)->walker + j);
      free((*Pptr)->walker);
    }
    *Pptr = NULL;
  }
}


int randomPopulation(Population P, int size){
  int i;
  double e, avg, min, max;
  
  randomBitstring(P->walker[0]);
  e = getPotential(P->walker[0], P->sat);
  avg = e;
  min = e;
  max = e;
  P->walker[0]->potential = e;
  
  for (i=1; i<size; ++i){
    randomBitstring(P->walker[i]);
    e = getPotential(P->walker[i], P->sat);
    avg += e;
    if ( e < min )
      min = e;
    if ( e > max )
      max = e;
    P->walker[i]->potential = e;
  }
  
  P->psize = size;
  P->avg_v = avg/size;
  P->max_v = max;
  P->min_v = min;

  return 0;
}


double getPotential(Bitstring bts, SAT sat){
  int i,j;
  int k,l;
  int p,q,v;
  word_t r;
  
  for (i=v=0; i<sat->num_clauses; ++i) { // Loop through the clauses; set the output to zero.
    l = sat->clause_length[i];           // Store off the length of the i-th clause.
    for (j=k=0; j<l; ++j) {              // Loop though the terms in the i-th clause; set truth value to false.
      p = sat->clause[i][j];             // Store off the j-th term of the i-th clause.
      q = abs(p);                        // This the index (1-up) of the variable in this term.
      --q;                               // This is the index (0-up) of the varible in this term.
      r = bts->node[q/BITS_PER_WORD];    // Get the correct word from the bitstring.
      r >>= q%BITS_PER_WORD;             // Shift the correct variable value to the lowest bit.
      r &= 1;                            // Zeroize the other bits.
      k |= (p > 0) ? r : 1-r;            // If the term is positive the value is kept, otherwise it is negated.
    }
    v += (1-k)*sat->clause_weight[i];    // If the clause is _false_ then add in the weight as penalty.
  }
  
  return (double) v;
}


void update(double a, double b, double mean, Population P, int parity){
  int i,j,k;
  double p,e;
  double min = P->sat->num_clauses;
  double avg = 0.0;
  double max = -(P->sat->num_clauses);
  int old = parity*arraysize;
  int new = (1-parity)*arraysize;
  
  for (i=j=0; i<P->psize; ++i) {   // Loop over each walker (i) and set target position (j) to zero.
    p = drand48();
    
    // First potential event: walker steps.
    if ( p < a ) {
      k = randomBitFlip(P->walker[new+j], P->walker[old+i]);
      e = P->walker[old+i]->potential + ((k>0)-(k<0))*getPotential(P->walker[new+j],P->ds->der[abs(k)-1]);
      if ( e < min ) min = e;
      if ( e > max ) max = e;
      avg += e;
      P->walker[new+j]->potential = e;
      ++j;
      continue;
    }
    p -= a;
   
    // Second potential event: walker spawns/dies.
    if ( mean > P->walker[old+i]->potential ) { // Case of spawning.
      if ( p < b*(mean - P->walker[old+i]->potential) ) {
        copyBitstring(P->walker[new+j], P->walker[old+i]);
        copyBitstring(P->walker[new+j+1], P->walker[old+i]);
        e = P->walker[old+i]->potential;
        if ( e < min ) min = e;
        if ( e > max ) max = e;
        avg += 2*e;
        j += 2;
        continue;
      }
    } else {
      if ( p < b*(P->walker[old+i]->potential - mean) ) { // Case of dying.
        continue;
      }
    }
    
    // Third potential event: walker stays.
    copyBitstring(P->walker[new+j], P->walker[old+i]);
    e = P->walker[old+i]->potential;
    if ( e < min ) min = e;
    if ( e > max ) max = e;
    avg += e;
    ++j;
  }
  
  P->psize = j;
  P->avg_v = avg/j;
  P->min_v = min;
  P->max_v = max;
}
