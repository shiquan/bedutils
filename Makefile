CC=		gcc
CFLAGS=		-g -Wall -O2
LOBJS=		bedutil.o commons.o
PROG=		bedutils
INCLUDES=	-I.
SUBDIRS=	. 
LIBPATH=        -L.

.SUFFIXES:.c .o
.PHONY: all

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:clean $(PROG) 

.PHONY:all  clean

bedutils: bedutil.o commons.o bedutil.h
		$(CC) $(CFLAGS) -o $@ bedutil.o commons.o $(LDFLAGS) bedtk.c $(LIBPATH) $(INCLUDES) -lz

commons.o:commons.c commons.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) commons.c -o $@	

bedutil.o:bedutil.c bedutil.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) bedutil.c -o $@	

clean:
		rm -fr gmon.out *.o a.out *.exe *.dSYM  $(PROG) *~ *.a 
