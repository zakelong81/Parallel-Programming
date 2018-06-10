CFLAGS= -Wall -W -Werror -g -lm

all: tmin ptrap

tmin: tmin.c
	gcc -o tmin tmin.c $(CFLAGS)
ptrap: ptrap.c
	mpicc -o ptrap ptrap.c $(CFLAGS)
clean:
					rm tmin mpicc
