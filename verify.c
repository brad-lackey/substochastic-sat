/*-----------------------------------------------------------------
  This software verifies kSAT solutions. The kSAT instances are
  read in DIMACS cnf format. The solution is given in the format
  specified for the MaxSAT 2015 contest. That is, the bits are
  numbered 1,2,3,4... If the xth bit is 1 then output x, and if
  the xth bit is 0 then output -x. This software was written in
  2016 as part of a collaboration between Stephen Jordan, Brad
  Lackey, and Michael Jarret.
  -----------------------------------------------------------------*/

//On machines with very old versions of glibc (e.g. the Raritan cluster)
//we need to define gnu_source in order to avoid warnings about implicit
//declaration of getline().
#define _GNU_SOURCE

#include <stdio.h>
//#include <malloc.h> //omit on mac
#include <string.h>
#include <stdlib.h> //for abs
#include <ctype.h>  //for isdigit and isspace

//variables are numbered 1,2,3,...
//negative means logical negation in the clause
typedef struct {
  int vars[4];
  int nots[4];
  int numvars;
}clause;

typedef struct {
  clause *clauses;   //the clauses
  int numclauses;    //how many clauses there are
  int B;             //how many bits there are
}instance;

void printbits(int *bits, int B) {
  int i;
  for(i = 0; i < B; i++) printf("%i ", bits[i]);
  printf("\n");
}

int violated(clause c, int *bits) {
  int i;
  for(i = 0; i < c.numvars; i++) if(c.nots[i]^bits[c.vars[i]]) return 0;
  return 1;
}

int numviolated(int *bits, instance *sat) {
  int i;
  int returnval;
  returnval = 0;
  for(i = 0; i < sat->numclauses; i++) returnval += violated(sat->clauses[i], bits);
  return returnval;
}

int loadbits(char *filename, int **bits) {
  char *current;
  int vars;
  int index;
  int scanned;
  int *localbits;
  int i;
  int unassigned;
  char *string;
  FILE *fp;
  size_t nbytes;
  int bytes_read;
  fp = fopen(filename, "r");
  if(fp == NULL) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }
  string = (char *)malloc(16384*sizeof(char));
  nbytes = 16384;
  bytes_read = getline(&string, &nbytes, fp);
  fclose(fp);
  if(bytes_read < 1 || bytes_read > 16383) {
    printf("%s is not a valid bitstring file.\n", filename);
    return 0;
  }
  vars = 0;
  current = string;
  do {
    scanned = sscanf(current, "%i", &index);
    if(scanned == 1) {
      if(abs(index) > vars) vars = abs(index);
      while(isspace(current[0])) current++;
      while(isdigit(current[0])||current[0] == '-') current++;
    }
  }while(scanned == 1);
  *bits = (int *)malloc(vars*sizeof(int));
  localbits = *bits;
  for(i = 0; i < vars; i++) localbits[i] = 3;
  current = string;
  do {
    scanned = sscanf(current, "%i", &index);
    if(scanned == 1) {
      if(index < 0) localbits[abs(index)-1] = 0;
      if(index > 0) localbits[abs(index)-1] = 1;
      while(isspace(current[0])) current++;
      while(isdigit(current[0])||current[0] == '-') current++;
    }
  }while(scanned == 1);
  i = 0;
  unassigned = 0;
  for(i = 0; i < vars && !unassigned; i++) {
    if(localbits[i] == 3) {
      fprintf(stderr, "Unassigned bits are present!\n");
      unassigned = 1;
    }
  }
  free(string);
  return vars;
}

//Here we load an instance of SAT in the DIMACS file format.
//This returns 0 on failure, 1 on success. We can handle
//clauses with up to 4 variables.
int loadsat(char *filename, instance *sat) {
  FILE *fp;
  size_t nbytes;
  int bytes_read;
  char *line;
  char junk1;
  char junk2[64];
  int vars, clauses;
  int i, j;
  int success;
  int numread;
  int x[5];
  nbytes = 255;
  line = (char *)malloc(256*sizeof(char));
  fp = fopen(filename, "r");
  if(fp == NULL) {
    fprintf(stderr, "Error: unable to open %s\n", filename);
    return 0;
  }
  success = 0;
  do {
    bytes_read = getline(&line, &nbytes, fp);
    //lines starting with c are comments
    //lines starting with p specify the parameters
    //lines with two or fewer bytes are either carriage return or %
    if(line[0] == 'p') {
      sscanf(line, "%c %s %i %i", &junk1, junk2, &vars, &clauses);
      success = 1;
    }
  }while(!success && bytes_read > 0);
  if(!success) {
    fprintf(stderr, "Finished scanning file without finding parameters.\n");
    return 0;
  }
  sat->clauses = (clause *)malloc(clauses*sizeof(clause));
  if(sat->clauses == NULL) {
    fprintf(stderr, "Memory allocation error in loadsat.\n");
    return 0;
  }
  sat->B = vars;
  sat->numclauses = clauses;
  fseek(fp, 0, SEEK_SET); //return to beginning
  i = 0;
  do {
    bytes_read = getline(&line, &nbytes, fp);    
    if(line[0] != 'c' && line[0] != 'p' && bytes_read > 2) {
      numread = sscanf(line, "%i %i %i %i %i", &x[0], &x[1], &x[2], &x[3], &x[4]);
      if(x[numread-1] != 0) printf("Warning: line %s not terminated with 0.\n", line);
      sat->clauses[i].numvars = numread-1;
      for(j = 0; j < numread-1; j++) {
        sat->clauses[i].vars[j] = abs(x[j])-1;
        sat->clauses[i].nots[j] = 0;
        if(x[j] < 0) sat->clauses[i].nots[j] = 1;
      }
      i++;
    }
  }while(bytes_read > 0);
  fclose(fp);
  free(line);
  return 1;
}

//print the SAT instance to stdout
void printsat(instance *sat) {
  int i,j;
  printf("c CNF formula in DIMACS format.\n");
  printf("p cnf %i %i\n", sat->B, sat->numclauses);
  for(i = 0; i < sat->numclauses; i++) {
    for(j = 0; j < sat->clauses[i].numvars; j++) {
      if(sat->clauses[i].nots[j]) printf("%i ", -(sat->clauses[i].vars[j]+1));
      else printf("%i ", sat->clauses[i].vars[j]+1);
    }
    printf("0\n");
  }
}

int main(int argc, char *argv[]) {
  instance sat;
  int *bits;
  int satsuccess;
  int vars;
  if(argc != 3) {
    printf("Usage: verify bitstring.txt instance.cnf\n");
    return 0;
  }
  satsuccess = loadsat("example.cnf", &sat);
  if(!satsuccess) return 0;
  vars = loadbits("example.txt", &bits);
  if(!vars) return 0;
  //printsat(&sat);
  printf("%i clauses violated\n", numviolated(bits, &sat));
  free(sat.clauses);
  free(bits);
  return 0;
}
