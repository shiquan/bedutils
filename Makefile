CC=		gcc
CFLAGS= -Wall
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
		$(CC) -O3 $(CFLAGS) $(LDFLAGS) $(LIBPATH) $(INCLUDES) -lz -o $@ bedutil.c commons.c kstring.c bedtk.c number.c

debug:
		$(CC) -g -O0 $(CFLAGS) $(LDFLAGS) $(LIBPATH) $(INCLUDES) -lz -D_DEBUG_MODE -o bedutils bedutil.c commons.c kstring.c bedtk.c number.c

clean:
		rm -fr gmon.out *.o a.out *.exe *.dSYM  $(PROG) *~ *.a 
