/*************************************************************************
	> File Name: bseq.h
	> Author: 
	> Mail: 
 ************************************************************************/

#ifndef _BSEQ_H
#define _BSEQ_H

#include <stdio.h>
#include "read_seeding.h"

// KSEQ_INIT(gzFile, gzread)

typedef struct bseq_file_s bseq_file_t;

bseq_file_t *bseq_open(const char *fn);

void bseq_close(bseq_file_t *fp);

uint32_t bseq_read(bseq_file_t *fp, uint32_t batch, READ_t* query_info);

uint32_t bseq_read_2pass(bseq_file_t *fp, uint32_t batch, seq_io* seqio);


#endif
