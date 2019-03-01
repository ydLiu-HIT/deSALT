/*************************************************************************
	> File Name: bseq.c
	> Author: 
	> Mail: 
 ************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include "bseq.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

struct bseq_file_s
{
    gzFile fp;
    kseq_t *ks;
};

bseq_file_t *bseq_open(const char *fn)
{
    bseq_file_t *fp;
    gzFile f;
    f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (f == 0) return 0;
    fp = (bseq_file_t* )calloc(1, sizeof(bseq_file_t));
    fp->fp = f;
    fp->ks = kseq_init(fp->fp);
    return fp;
}

void bseq_close(bseq_file_t *fp)
{
    kseq_destroy(fp->ks);
    gzclose(fp->fp);
    free(fp);
}

uint32_t bseq_read(bseq_file_t *fp, uint32_t batch, READ_t* query_info)
{
    uint32_t seqii = 0;
    int64_t kr1 = 1;
    uint32_t i;
    kseq_t *seq1 = fp->ks;

    for(seqii = 0; (seqii < batch) && ((kr1 = kseq_read(seq1)) > 0); seqii++)
    {
        //strdup, later process like before and compare the memory
        query_info[seqii].name = strdup(seq1->name.s);
        query_info[seqii].read_seq = strdup(seq1->seq.s);
        query_info[seqii].read_length = seq1->seq.l;
        for (i = 0; i < seq1->seq.l; ++i) //convert U to T 
            if (query_info[seqii].read_seq[i] == 'u' || query_info[seqii].read_seq[i] == 'U')
                --query_info[seqii].read_seq[i];
    }

    return seqii;
}

uint32_t bseq_read_2pass(bseq_file_t *fp, uint32_t batch, seq_io* seqio)
{
    uint32_t seqii = 0;
    int64_t kr1 = 1;
    uint32_t i;
    kseq_t *seq1 = fp->ks;

    for(seqii = 0; (seqii < batch) && ((kr1 = kseq_read(seq1)) > 0); seqii++)
    {
        seqio[seqii].name = strdup(seq1->name.s);
        seqio[seqii].read_seq = strdup(seq1->seq.s);
        seqio[seqii].qual = seq1->qual.l? strdup(seq1->qual.s):0;
        seqio[seqii].read_length = seq1->seq.l;
        for (i = 0; i < seq1->seq.l; ++i) //convert U to T
            if (seqio[seqii].read_seq[i] == 'u' || seqio[seqii].read_seq[i] == 'U')
                --seqio[seqii].read_seq[i];
    }

    return seqii;
}
