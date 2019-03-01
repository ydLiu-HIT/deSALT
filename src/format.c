/*************************************************************************
	> File Name: format.c
	> Author: 
	> Mail: 
 ************************************************************************/

#include<stdio.h>
#include<stdint.h>

#include "read_seeding.h"
#include "load_unipath_size.h"

void ff_print_line(_aln_t *aln, char *query_name, char *read_seq, char *qual, uint16_t mapq, uint32_t id, param_map *opt)
{
    fprintf(fp_sam, "%s\t%u\t%s\t%u\t%u\t",
                    query_name, aln->flag, chr_names[aln->chr_n], aln->_1_based_pos,
                    mapq);
    uint32_t c;

    for (c = 0; c < aln->n_cigar; ++c) 
    {
        fprintf(fp_sam, "%d%c", (aln->cigar[c])>>4, "MIDNS"[(aln->cigar[c])&0xf]);// print cigar
    }
    fprintf(fp_sam, "\t*\t0\t0\t%s", read_seq);
    if(qual && opt->with_qual) fprintf(fp_sam, "\t%s", qual);
    else	fprintf(fp_sam, "\t*" );
    

    //fprintf(fp_sam, "\tblen:%d\tmlen:%d\tn_ambi:%d\tNM:i:%d\tms:i:%d\tAS:i:%d\tnn:i:%d", aln->blen, aln->mlen, aln->n_ambi, aln->blen - aln->mlen + aln->n_ambi, aln->dp_max, aln->dp_score, aln->n_ambi);
    fprintf(fp_sam, "\tNM:i:%d\tms:i:%d\tAS:i:%d\tnn:i:%d", aln->blen - aln->mlen + aln->n_ambi, aln->dp_max, aln->dp_score, aln->n_ambi);
    fprintf(fp_sam, "\n");
}

int ff_print_sam (seq_io *seqio, uint32_t read_cnt, param_map *opt)
{
    int mapped_reads = 0;
	uint32_t r_i, r_ii;
    //output sam
    for ( r_i = 0; r_i < read_cnt; ++r_i)
    {
    	// fprintf(stderr, "flag = %u\n", aln->flag);
        if (seqio[r_i].mapable)
        {
            //print the best alignment
            ff_print_line(&seqio[r_i].aln[0], seqio[r_i].name, seqio[r_i].read_seq, seqio[r_i].qual, seqio[r_i].mapq, r_i, opt);
            for (r_ii = 1; r_ii < seqio[r_i].multi_n; ++r_ii)
            {
                ff_print_line(&seqio[r_i].aln[r_ii], seqio[r_i].name, "*", "*", 0, r_i, opt);
            }
            mapped_reads ++;
        }
        else
        {
           	//unmapped read
           	fprintf(fp_sam, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s", seqio[r_i].name, seqio[r_i].read_seq);
        	if(seqio[r_i].qual && opt->with_qual) fprintf(fp_sam, "\t%s", seqio[r_i].qual);
            else	fprintf(fp_sam, "\t*" );
            fprintf(fp_sam, "\tXO:Z:NM\n");
        }
    }
    return mapped_reads;
}