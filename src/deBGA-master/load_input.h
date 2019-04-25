
#ifndef LOAD_INPUT_H_
#define LOAD_INPUT_H_

#include <stdio.h>
#include <stdint.h>
#include <getopt.h> 
#include <zlib.h>
#include "kseq.h"

#define	CHR_NAME_SPLIT
#define	HANDLE_DIR

#define UNPIPATH_OFF_K20

#define	UNI_SEQ64
#define PACKAGE_VERSION "0.1"

#define ROUTE_LENGTH_MAX   500

#define	MAX_CHR_NUM 600000
#define	MAX_CHR_NAME_LENGTH	200

#define	START_POS_REF	0

KSEQ_INIT(gzFile, gzread)
extern kseq_t *seq1;
extern kseq_t *seq2;
gzFile fp1;
gzFile fp2;

extern uint64_t* buffer_ref_seq;
#ifdef UNI_SEQ64
extern uint64_t* buffer_seq;
#else
extern uint8_t* buffer_seq;
#endif

#ifdef UNPIPATH_OFF_K20
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

extern uint8_t* buffer_edge;
extern uint32_t* buffer_pupos;
extern uint32_t* buffer_puid;

extern uint32_t* buffer_kmer_g;
	
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

extern const char* refseqs;
extern const char* uniseqs;
extern const char* uniseq_bs;
extern const char* uniseqf_bs;
extern const char* uniedges;
extern const char* unipus;
extern const char* uniposs;
extern const char* uniposps;
extern const char* unistas;
extern const char* unihash_gs;
extern const char* unikmer_gs;
extern const char* unioff_gs;
extern const char* unichrs;
extern const char* f_n;
extern const char* graph;
extern const char* divs;
extern const char* f_size;
extern const char* sys_c_mkdir;
extern const char* sys_c_rm;
extern const char* sys_rm;

extern char sam_result[ROUTE_LENGTH_MAX];
extern char read_fastq1[ROUTE_LENGTH_MAX];
extern char read_fastq2[ROUTE_LENGTH_MAX];
extern char index_route[ROUTE_LENGTH_MAX];
extern char filename_ref[ROUTE_LENGTH_MAX];
extern char filename_div[ROUTE_LENGTH_MAX];
extern char filename_sta[ROUTE_LENGTH_MAX];
extern char ref_seq[ROUTE_LENGTH_MAX];
extern char uniseq[ROUTE_LENGTH_MAX];
extern char uniseq_b[ROUTE_LENGTH_MAX];
extern char uniseqf_b[ROUTE_LENGTH_MAX];
extern char uniedge[ROUTE_LENGTH_MAX];
extern char unipu[ROUTE_LENGTH_MAX];
extern char unipos[ROUTE_LENGTH_MAX];
extern char uniposp[ROUTE_LENGTH_MAX];
//extern char unista[ROUTE_LENGTH_MAX];
extern char unihash_g[ROUTE_LENGTH_MAX];
extern char unikmer_g[ROUTE_LENGTH_MAX];
extern char unioff_g[ROUTE_LENGTH_MAX];
extern char unichr[ROUTE_LENGTH_MAX];
extern char N_route[ROUTE_LENGTH_MAX];
extern char unisize[ROUTE_LENGTH_MAX];
extern char chr_names[MAX_CHR_NUM][MAX_CHR_NAME_LENGTH];
extern char chr_line_content[MAX_CHR_NAME_LENGTH];

//files
extern FILE* fp_ref_seq;
extern FILE* fp_us_b;
extern FILE* fp_usf_b;
extern FILE* fp_ub;
extern FILE* fp_ue;
extern FILE* fp_up;
extern FILE* fp_upp;
extern FILE* fp_pu;
extern FILE* fp_sta;
extern FILE* fp_chr;
extern FILE* fp_uh;
extern FILE* fp_uf;
extern FILE* fp_hash;
extern FILE* fp_kmer;
extern FILE* fp_off;
extern FILE* fp_n;
extern FILE* fp_us;
extern FILE* fp_num;
extern FILE* unipath_debug;

extern const uint8_t f;
extern const uint8_t k;
extern uint8_t k_t;	

#ifdef	UNPIPATH_OFF_K20
extern uint64_t chr_end_n[MAX_CHR_NUM];
extern uint64_t first_cnt_cl;
#else
extern uint32_t chr_end_n[MAX_CHR_NUM];
extern uint32_t first_cnt_cl;
#endif

extern uint32_t chr_file_n;
extern uint8_t thread_n;
extern uint32_t upper_ins;
extern uint32_t floor_ins;
extern uint16_t readlen_max;
extern uint16_t seed_l_max;
extern uint8_t seed_l_l;
extern uint16_t seed_step;
extern uint8_t cir_fix_n;
extern uint16_t pos_n_max;
extern uint16_t length_reduce;
extern float mis_match_r;
extern float score_filt_r;
extern int32_t match;
//extern mismatch;
extern float va_ra;

extern uint16_t cus_ali_n;
extern uint16_t cus_max_output_ali;
extern float max_pair_score_r;
extern float last_circle_rate;
extern uint8_t local_ksw;
extern uint8_t mgn_flag;
extern uint8_t flag_std;

extern float max_single_score_r;
extern uint16_t pr_single_outputn;
extern uint16_t seed_filter_pos_numn;
extern uint16_t seed_filter_pos_num_singlen;
extern float lv_rate;

extern float mis_match_r_single;
extern float lv_rate_anchor;

extern int8_t mat_score;
extern int8_t mis_score;
extern int8_t gapo_score;
extern int8_t gape_score;

int load_input_index(int , char * []);
int load_input_map(int , char * []);
int load_input_map_single_end(int , char * []);
void load_index_file();
int seed_ali();
int seed_ali_single_end();

#endif /* LOAD_INPUT_H_ */
