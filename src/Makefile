CFLAGS=		-g -Wall -O2 -Wc++-compat
CPPFLAGS=	-DHAVE_KALLOC
#CPPFLAGS=
INCLUDES=
OBJS=		main.o desalt_index.o read_seeding.o bit_operation.o format.o binarys_qsort.o hash_index.o aln_2pass.o splic_junction.o bseq.o load_unipath_size.o graph.o kalloc.o ksw2_ll_sse.o
PROG=		deSALT
# PROG_EXTRA=	sdust minimap2-lite
LIBS=		-lm -lz -lpthread -lgomp

ifeq ($(sse2only),)
	OBJS+=ksw2_extz2_sse41.o ksw2_extd2_sse41.o ksw2_exts2_sse41.o ksw2_extz2_sse2.o ksw2_extd2_sse2.o ksw2_exts2_sse2.o ksw2_dispatch.o
else
	OBJS+=ksw2_extz2_sse.o ksw2_extd2_sse.o ksw2_exts2_sse.o
endif

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@  #$(CFLAGS)

all:$(PROG)

# extra:all $(PROG_EXTRA)
deSALT:$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

ksw2_extz2_sse41.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c -msse4 $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extz2_sse2.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_extd2_sse41.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c -msse4 $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extd2_sse2.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_exts2_sse41.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c -msse4 $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_exts2_sse2.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_dispatch.o:ksw2_dispatch.c ksw2.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

clean:
		@rm -f *.o $(PROG)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# # DO NOT DELETE

# align.o: minimap.h mmpriv.h bseq.h ksw2.h kalloc.h
# bseq.o: bseq.h kvec.h kalloc.h kseq.h
# chain.o: minimap.h mmpriv.h bseq.h kalloc.h
# example.o: minimap.h kseq.h
# format.o: kalloc.h mmpriv.h minimap.h bseq.h
# getopt.o: getopt.h
# hit.o: mmpriv.h minimap.h bseq.h kalloc.h khash.h
# index.o: kthread.h bseq.h minimap.h mmpriv.h kvec.h kalloc.h khash.h
# kalloc.o: kalloc.h
# ksw2_extd2_sse.o: ksw2.h kalloc.h
# ksw2_exts2_sse.o: ksw2.h kalloc.h
# ksw2_extz2_sse.o: ksw2.h kalloc.h
# ksw2_ll_sse.o: ksw2.h kalloc.h
# main.o: bseq.h minimap.h mmpriv.h getopt.h
# map.o: kthread.h kvec.h kalloc.h sdust.h mmpriv.h minimap.h bseq.h khash.h
# misc.o: minimap.h ksort.h
# pe.o: mmpriv.h minimap.h bseq.h kvec.h kalloc.h ksort.h
# sdust.o: kalloc.h kdq.h kvec.h sdust.h
# sketch.o: kvec.h kalloc.h minimap.h
