CC=gcc
CFLAGS=-O3 -funroll-loops -Wall
LINK=-lm
INC=
LIB=

all : ssmc

ssmc : substochastic.c bitstring.o sat.o population.o LUT.o
	$(CC) $(CFLAGS) $(INC) substochastic.c bitstring.o sat.o population.o lut.o $(LIB) $(LINK) -o ssmc
	strip ssmc

bitstring.o : bitstring.c bitstring.h macros.h
	$(CC) $(CFLAGS) $(INC) -c bitstring.c 

sat.o : sat.c sat.h macros.h bitstring.h
	$(CC) $(CFLAGS) $(INC) -c sat.c

population.o: population.c population.h macros.h bitstring.h sat.h
	$(CC) $(CFLAGS) $(INC) -c population.c

LUT.o : lut.c lut.h macros.h
	$(CC) $(CFLAGS) $(INC) -c lut.c

clean :
	rm -f *~ ssmc *.o
