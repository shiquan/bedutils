CC=		gcc
CFLAGS=		-g -Wall -O0
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

bedutils:
		$(CC) $(CFLAGS) $(LDFLAGS) $(LIBPATH) $(INCLUDES) -lz -o $@ bedutil.c commons.c kstring.c bedtk.c

clean:
		rm -fr gmon.out *.o a.out *.exe *.dSYM  $(PROG) *~ *.a 
