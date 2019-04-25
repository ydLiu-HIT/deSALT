CC =	gcc
CFLAGS =	-g -Wall -O2
#SRCS =	main.c hello.c
SRCS =	$(wildcard *.c)
#OBJS =	$(SRCS:.c=.o)
OBJS =	index_build.o seed_ali_core.o load_input.o bit_operation.o LandauVishkin.o ksw.o main.o
LIBS =	-lz -lpthread -lgomp -lm
INCLUDES =
SUBDIRS=	.

BIN_DIR = ..
PROG =	${BIN_DIR}/deBGA

.SUFFIXES:.c .o .cc

.c.o:
	$(CC) -c $(INCLUDES) $< -o $@
#$(CFLAGS)
all:$(PROG)

#dep:$(DEPS)

$(PROG):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	@rm -f *.o $(PROG)
	
depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) -- *.c )
	
index_build.o: index_build.h bit_operation.h load_input.h
seed_ali_core.o: binarys_qsort.h LandauVishkin.h load_input.h seed_ali_p.h bit_operation.h
load_input.o: load_input.h
bit_operation.o: bit_operation.h
LandauVishkin.o: LandauVishkin.h
ksw.o: ksw.h
