#include <stdio.h>
//#include <malloc.h>
#include <string.h>

#define MAXRUNS 5000

typedef struct {
  double walltime;
  char filename[256];
  int trials;
  int bits;
  int answer;
  int achieved[200];
  double runtime;
  int popsize;
}run;

void printrun(run *r) {
  int i;
  printf("Filename: %s\n", r->filename);
  printf("trials: %i\n", r->trials-1);
  printf("bits: %i\n", r->bits);
  printf("answer: %i\n", r->answer);
  printf("runtime: %e\n", r->runtime);
  printf("walltime: %f\n", r->walltime);
  printf("popsize: %i\n", r->popsize);
  for(i = 0; i < r->trials-1; i++)
    printf("optimum from trial %i: %i\n", i, r->achieved[i]);
  printf("final optimum: %i\n", r->achieved[r->trials-1]);
}

int loadruns(run *runs, FILE *optima, FILE *input) {
  int i, thisrun, scanned;
  int answers;
  char *line;
  int bytes_read;
  size_t nbytes;
  int numruns;  
  char junk1[256];
  char junk2[256];
  char opt_filename[1000][256];
  int opt[1000];
  char junkchar1, junkchar2;
  char *tmp;
  int answer_found;
  //load the optima that were obtained from the maxsat website
  i = 0;
  do {
    scanned = fscanf(optima, "%s %c %c %i", opt_filename[i], &junkchar1, &junkchar2, &opt[i]);
    i++;
  }while(scanned == 4);
  answers = i-1;
  //load the data about the runs
  line = (char *)malloc(512*sizeof(char));
  nbytes = 512;
  numruns = 0;
  do {
    bytes_read = getline(&line, &nbytes, input);
    if(bytes_read > 0) {
      if(strstr(line, "c Input")) { //filename
        tmp = strrchr(line, '/');
        for(i = 0; i < strlen(tmp)-2; i++) runs[numruns].filename[i] = tmp[i+1];
        runs[numruns].filename[i] = '\0';
        runs[numruns].trials = 0;
      }
      if(line[0] == 'o') {//output of optimum
	sscanf(line+2, "%i", &runs[numruns].achieved[runs[numruns].trials]);
        runs[numruns].trials++;
      }
      if(strstr(line, "c Walltime")) { //total walltime
        sscanf (line, "%c %s %lf", &junkchar1, junk1, &runs[numruns].walltime);
        numruns++;
      }
      if(strstr(line, "c Bits")) //number of bits
        sscanf(line, "%c %s %i", &junkchar1, junk1, &runs[numruns].bits);
      if(strstr(line, "c Runtime")) //physical runtime (hbar = 1)
        sscanf(line, "%c %s %lf", &junkchar1, junk1, &runs[numruns].runtime);
      if(strstr(line, "c Population")) //size of population of walkers
        sscanf(line, "%c %s %s %i", &junkchar1, junk1, junk2, &runs[numruns].popsize);
    }
  }while(bytes_read > 0);
  //find out which runs found the optima
  for(thisrun = 0; thisrun < numruns; thisrun++) {
    answer_found = 0;
    for(i = 0; i < answers && !answer_found; i++) {
      if(strcmp(runs[thisrun].filename, opt_filename[i]) == 0) {
        runs[thisrun].answer = opt[i];
        answer_found = 1;
      }
    }
    if(!answer_found) printf("Error: answer not found for %s!\n", runs[thisrun].filename);
  }
  free(line);
  return numruns;
}

void binnit(int numruns, run *r) {
  double runtimelist[100];
  int poplist[100];
  int runtimes;
  int pops;
  int found;
  int i,j;
  int popbin[MAXRUNS];
  int timebin[MAXRUNS];
  int rawsuccess[100][100];
  int rawruns[100][100];
  int finalruns[100][100];
  int finalsuccess[100][100];
  double finalwalltime[100][100];
  runtimes = 0;
  pops = 0;
  for(i = 0; i < numruns; i++) {
    found = 0;
    for(j = 0; j < runtimes; j++) if(runtimelist[j] == r[i].runtime) found = 1;
    if(!found) runtimelist[runtimes++] = r[i].runtime;
    found = 0;
    for(j = 0; j < pops; j++) if(poplist[j] == r[i].popsize) found = 1;
    if(!found) poplist[pops++] = r[i].popsize;
  }
  printf("%i populations used:\n", pops);
  for(i = 0; i < pops; i++) printf("%i\n", poplist[i]);
  printf("%i runtimes used:\n", runtimes);
  for(i = 0; i < runtimes; i++) printf("%e\n", runtimelist[i]);
  for(i = 0; i < numruns; i++) {
    for(j = 0; j < pops; j++) if(r[i].popsize == poplist[j]) popbin[i] = j;
    for(j = 0; j < runtimes; j++) if(r[i].runtime == runtimelist[j]) timebin[i] = j;
  }
  //for(i = 0; i < numruns; i++)
  //  printf("run %i: pop %i: runtime %e: popbin: %i->%i: timebin: %i->%e\n", 
  //  i, r[i].popsize, r[i].runtime, popbin[i], poplist[popbin[i]], timebin[i], runtimelist[timebin[i]]);
  for(i = 0; i < pops; i++) {
    for(j = 0; j < runtimes; j++) {
      rawsuccess[i][j] = 0;
      finalsuccess[i][j] = 0;
      rawruns[i][j] = 0;
      finalruns[i][j] = 0;
      finalwalltime[i][j] = 0;
    }
  }
  for(i = 0; i < numruns; i++) {
    for(j = 0; j < r[i].trials-1; j++) {
      rawruns[popbin[i]][timebin[i]] += 1;
      if(r[i].achieved[j] == r[i].answer) rawsuccess[popbin[i]][timebin[i]] += 1;
    }
    finalwalltime[popbin[i]][timebin[i]] += r[i].walltime;
    finalruns[popbin[i]][timebin[i]] += 1;
    if(r[i].achieved[r[i].trials-1] == r[i].answer) finalsuccess[popbin[i]][timebin[i]] += 1;
  }
  for(i = 0; i < pops; i++) {
    printf("population: %i------------------------\n", poplist[i]);
    for(j = 0; j < runtimes; j++) {
      printf("runtime: %e\n", runtimelist[j]);
      printf("Raw success = %i/%i = %f\n", rawsuccess[i][j], rawruns[i][j], (double)rawsuccess[i][j]/(double)rawruns[i][j]);
      printf("Final success = %i/%i = %f\n", finalsuccess[i][j], finalruns[i][j], (double)finalsuccess[i][j]/(double)finalruns[i][j]);
      printf("Average walltime = %e\n", finalwalltime[i][j]/(double)rawruns[i][j]);
    }
  }
}

int main(int argc, char *argv[]) {
  FILE *input;
  FILE *optima;
  run runs[MAXRUNS];
  int numruns;
  if(argc != 2) {
    printf("Usage: process data.txt\n");
    return 0;
  }
  input = fopen(argv[1], "r");
  if(input == NULL) {
    printf("Unable to open %s\n", argv[1]);
    return 0;
  }
  optima = fopen("optima.dat", "r");
  if(optima == NULL) {
    printf("Unable to open optima.dat\n");
    return 0;
  }
  numruns = loadruns(runs, optima, input);
  //for(i = 0; i < numruns; i++) printrun(&runs[i]);
  binnit(numruns, runs);
  fclose(optima);
  fclose(input);
  return 0;
}
