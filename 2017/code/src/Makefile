CC=gcc
CFLAGS=-g -O3 -funroll-loops -Wall
LINK=-lm
INC=
LIB=

BINPATH = ../../bin

all : ssmc

ssmc : substochastic.c bitstring.o sat.o population.o LUT.o
	$(CC) $(CFLAGS) $(INC) substochastic.c bitstring.o sat.o population.o LUT.o $(LIB) $(LINK) -o $(BINPATH)/ssmc	

greedy : greedy.c bitstring.o sat.o population.o
	$(CC) $(CFLAGS) $(INC) greedy.c bitstring.o sat.o population.o $(LIB) $(LINK) -o greedy

bitstring.o : bitstring.c bitstring.h macros.h
	$(CC) $(CFLAGS) $(INC) -c bitstring.c 

sat.o : sat.c sat.h macros.h bitstring.h
	$(CC) $(CFLAGS) $(INC) -c sat.c

population.o: population.c population.h macros.h bitstring.h sat.h
	$(CC) $(CFLAGS) $(INC) -c population.c

LUT.o : LUT.c LUT.h macros.h
	$(CC) $(CFLAGS) $(INC) -c LUT.c

clean :
	rm -f *~ ssmc *.o
