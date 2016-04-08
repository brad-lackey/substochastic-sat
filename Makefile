CC=gcc
#CFLAGS=-O3 -fvectorize -funroll-loops
CFLAGS=-O3 -funroll-loops -Wall
INC=-I/opt/local/include
LIB=-L/opt/local/lib
LINK=-lm -lgmp

all : substochastic verify process

substochastic : substochastic.c bitstring.o sat.o population.o
	$(CC) $(CFLAGS) $(INC) substochastic.c bitstring.o sat.o population.o $(LIB) $(LINK) -o substochastic
	strip substochastic

verify : verify.c
	$(CC) -Wall verify.c -o verify

process : process.c
	$(CC) -Wall process.c -o process

bitstring.o : bitstring.c bitstring.h macros.h
	$(CC) $(CFLAGS) $(INC) -c bitstring.c 

sat.o : sat.c sat.h macros.h bitstring.h
	$(CC) $(CFLAGS) $(INC) -c sat.c

population.o: population.c population.h macros.h bitstring.h sat.h
	$(CC) $(CFLAGS) $(INC) -c population.c

clean :
	rm -f *~ verify substochastic process *.o
