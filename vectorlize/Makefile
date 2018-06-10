CC=icc
CFLAGS=-Wall -W -Werror
VEC= -march=core-avx2 -O2 -qopt-report=2 -qopt-report-phase=vec
NOVEC= -march=core-avx2 -O2 -no-vec


all: mmvec mmnovec
mmvec: mm.c
	      $(CC) $(CFLAGS) $(VEC) -o mmvec mm.c
mmnovec: mm.c
			  $(CC) $(CFLAGS) $(NOVEC) -o mmnovec mm.c
clean:
				rm mmvec mmnovec
