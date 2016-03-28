CC=gcc
#CFLAGS=-O3 -fvectorize -funroll-loops
CFLAGS=-O3 -funroll-loops -Wall

all : substochastic verify process

substochastic : substochastic.c bitstring.o sat.o
	$(CC) $(CFLAGS) substochastic.c bitstring.o sat.o -lm -o substochastic

verify : verify.c
	$(CC) -Wall verify.c -o verify

process : process.c
	$(CC) -Wall process.c -o process

bitstring.o : bitstring.c bitstring.h macros.h
	$(CC) $(CFLAGS) -c bitstring.c 

sat.o : sat.c sat.h macros.h
	$(CC) $(CFLAGS) -c sat.c

clean :
	rm -f *~ verify substochastic process *.o
