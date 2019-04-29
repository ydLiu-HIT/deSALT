/*************************************************************************
	> File Name: load_unipath_size.h
	> Author: 
	> Mail: 
 ************************************************************************/

#ifndef _LOAD_UNIPATH_SIZE_H
#define _LOAD_UNIPATH_SIZE_H

#define MAX_CHR_NAME_LENGTH 200
#define MAX_CHR_NUM 6000000
#define START_POS_REF 0
#define ROUTE_LENGTH_MAX 1024

extern uint64_t* buffer_ref_seq;
extern uint64_t* buffer_seq;

extern uint64_t* buffer_seqf;
extern uint64_t* buffer_off_g;
extern uint64_t* buffer_p;
extern uint64_t* buffer_pp;
extern uint64_t* buffer_hash_g;

extern uint32_t* buffer_kmer_g;//
	
extern uint64_t result_ref_seq;
extern uint64_t result_seq;
extern uint64_t result_seqf;
extern uint64_t result_edge;
extern uint64_t result_p;
extern uint64_t result_pp;
extern uint64_t result_pu;
extern uint64_t result_hash_g;
extern uint64_t result_kmer_g;
extern uint64_t result_off_g;
extern uint64_t result_ref_g;

extern uint32_t chr_end_n[MAX_CHR_NUM];
extern char chr_names[MAX_CHR_NUM][MAX_CHR_NAME_LENGTH];
extern char chr_line_content[MAX_CHR_NAME_LENGTH];
extern int chr_file_n;
extern uint64_t reference_len;

int load_index_file(char *index_dir);

#endif
