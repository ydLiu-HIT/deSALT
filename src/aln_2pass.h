/*************************************************************************
	> File Name: aln_2pass.h
	> Author: 
	> Mail: 
 ************************************************************************/

#ifndef _ALN_2PASS_H
#define _ALN_2PASS_H

#include "read_seeding.h"
#include "ksw2.h"

typedef struct
{
	uint32_t ts;
	uint32_t te;
	int key; //order in anchor_map2ref
} REF_t;

typedef struct
{
	uint32_t ts_f;
	uint32_t te_f;
	uint32_t ts_r;
	uint32_t te_r;
} EXON_t;

typedef struct
{
    //char chr[1204];
	uint32_t start;
	uint32_t end;
	uint32_t strand;
} Anno_t;

typedef struct
{
    uint32_t start;
    uint32_t end;
    uint32_t upper;
    uint32_t down;
}Anno_range_t;

typedef struct
{
    uint32_t Is;
    uint32_t Ie;
    int32_t strand;
}Ival_anno_t;

typedef struct
{
    int32_t Intron_n;
    Ival_anno_t *Ival;
}Ival_Anno_t;


//pthread
typedef struct THREAD_2PASS
{
	uint8_t tid;
	uint32_t seqn;
	ksw_extz_t ez;
	ksw_extz_t ez2;
	param_map *map;
	dpSkeleton_t *dp_skeleton;
	TARGET_t *anchor_map2ref;
	void *km;
} thread_2pass_t;

REF_t** REF_pos;
QUERY_t** QUERY_pos;
uint8_t **strand_arr;
EXON_t* EXON_T;
seq_io* seqio;


uint32_t merge_anchor_cnt;
uint8_t splice_offset;
uint8_t hash_kmer;
uint8_t e_shift;

void load_fasta_2pass(uint32_t map2ref_cnt, param_map *opt, char *read_fastq, int *load_fasta_2pass);
void test_main(param_map *opt);
void ksw_init();
void ksw_free_mem();

//Ival_Anno_t *read_Annotation(const char *annoDir);
Ival_Anno_t *read_Annotation(param_map *opt, TARGET_t *anchor_map2ref, uint32_t map2ref_cnt, uint32_t *merger_cnt);
void get_refseq(uint8_t *ref, uint32_t len, uint32_t start);
void get_ref_jump(uint8_t *ref, uint32_t start1, uint32_t len1, uint32_t start2, uint32_t len2);
void get_ref_onebyone(uint8_t *ref, uint32_t start, uint32_t len, uint32_t pre_pos);
int find_chr_n_by_name(char *chr_name);
int chromosome_judge(uint32_t position, uint32_t *chromosome_begin);
void mm_append_cigar(_aln_t *aln, uint32_t n_cigar, uint32_t *cigar);
void copy_aln_value(_aln_t *aln1, _aln_t *aln2, uint32_t read_line);
void mm_update_extra(_aln_t *aln, uint8_t *qseq, uint8_t *tseq, uint32_t qlen, uint32_t tlen, uint32_t *qs, uint32_t *ts, uint8_t q, uint8_t e);
uint32_t append_intron_to_cigar(void *km, ksw_extz_t *ez, uint32_t pre_pos, uint32_t intron_pos, uint32_t intron_len);
void check_more_part(uint32_t *cigar, uint32_t *n_cigar, int *score, int q2);
void check_cigar(uint8_t *qseq, uint8_t *tseq, uint32_t *cigar, uint32_t *n_cigar, uint32_t *qs, uint32_t *ts, int *score, uint32_t qlen, uint32_t tlen, uint32_t boundary, uint8_t type, param_map *opt);
int check_realign(void *km, param_map *opt, int bandwith, uint8_t *qseq, uint32_t qlen, uint32_t tlen, ksw_extz_t *ez, int pre_score, uint32_t start_pos, int intron_cnt, int splice_flag);
void align_splic_FOR_REV(void *km, uint8_t *qseq, uint8_t *tseq, uint32_t qlen, uint32_t tlen, int splice_flag, param_map *opt, ksw_extz_t *ez, uint8_t splice_type, uint8_t *junc);
void align_non_splice(void *km, uint8_t *qseq, uint8_t *tseq, uint32_t qlen, uint32_t tlen, param_map *opt, ksw_extz_t *ez, int bandwith, int flag, uint8_t type);


#endif
