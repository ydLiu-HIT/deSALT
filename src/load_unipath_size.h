/*************************************************************************
	> File Name: load_unipath_size.h
	> Author: 
	> Mail: 
	> Created Time: 2017年04月05日 星期三 10时17分17秒
 ************************************************************************/

#ifndef _LOAD_UNIPATH_SIZE_H
#define _LOAD_UNIPATH_SIZE_H


#define UNPIPATH_OFF_K20



#ifdef UNPIPATH_OFF_K20

#define SIZE64
#define MAX_CHR_NAME_LENGTH 200
#define MAX_CHR_NUM 6000000
#define START_POS_REF 0
#define UNI_SEQ64
#define ROUTE_LENGTH_MAX 1024

extern uint64_t* buffer_ref_seq;
#ifdef UNI_SEQ64
extern uint64_t* buffer_seq;
#else
extern uint8_t* buffer_seq;
#endif

extern uint64_t* buffer_seqf;
extern uint64_t* buffer_off_g;
extern uint64_t* buffer_p;
extern uint64_t* buffer_pp;
extern uint64_t* buffer_hash_g;
#else
extern uint32_t* buffer_seqf;
extern uint32_t* buffer_off_g;
extern uint32_t* buffer_p;
extern uint32_t* buffer_pp;
extern uint32_t* buffer_hash_g;
#endif

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

#ifdef	UNPIPATH_OFF_K20
extern uint64_t chr_end_n[MAX_CHR_NUM];
extern char chr_names[MAX_CHR_NUM][MAX_CHR_NAME_LENGTH];
extern char chr_line_content[MAX_CHR_NAME_LENGTH];
extern uint32_t chr_file_n;
extern uint64_t reference_len;
#else
extern uint32_t chr_end_n[MAX_CHR_NUM];
extern uint32_t first_cnt_cl;
#endif

int load_index_file(char *index_dir);

#endif
