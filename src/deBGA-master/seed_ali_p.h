
#ifndef SEED_ALI_H_
#define SEED_ALI_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "load_input.h"

//#define RD_FILED

#define	FIX_SV
#define FIX_SA

#define	PAR_OUPUT
//define SAMTOOLS_BUG
//#define PICARD_BUG
//#define POST_DEBUG
#define	REDUCE_ANCHOR
#define ANCHOR_LV_S

#ifdef	REDUCE_ANCHOR
uint8_t anchor_seed_length_thr = 10;//11: 1993234; 12 1992946; 15 1990***; 20 1989***; 30 1988***
#else
uint8_t anchor_seed_length_thr = 0;	
#endif

#define ANCHOR_LV_S_FLOAT 0.5
//0.3 
//0.2

#define	CIGAR_LEN_ERR
#define	READN_RANDOM_SEED

#define	QUAL_20

#define	R_W_LOCK

#define	DEBUG_ALN_64BIT

#define	MERGE_MODIFY
#define	POUND_MODIFY
//#define	CIGAR_S_MODIFY

//#define ALL_ALL_SINGLE
//#define	ALL_ALL
#define	DM_COPY_PAIR

#define	DM_COPY_SINGLE

#define	CIR_CONTINUE_SINGLE_END
#define	CHAR_CP_SINGLE_END
#define	QUAL_FILT_SINGLE_END
#define	QUAL_FILT_LV_MIS_SINGLE_END
#define	QUAL_FILT_LV_OUT_SINGLE_END
#define	MAPPING_QUALITY_SINGLE_END
#define	ALTER_DEBUG_SINGLE_END

#define	ALTER_DEBUG
#define	ALTER_DEBUG_ANCHOR
//#define	MP_DEBUG
//#define	ALTER_DEBUG_PRINT
//#define	A_DEBUG
//#define	S_DEBUG
//#define	NO_S_OFF	

#define	CHAR_CP

//UNMATCH_SINGLE_END: original anchor
//UNMATCH_SINGLE_END_MIS: new non-indel anchor
//CIR_JUMP: If there is matched seeds which don't get a reasonable LV result in current circle, use --cl to set new circle filter and lv rate for next circle
//notice: UNMATCH_SINGLE_END and UNMATCH_SINGLE_END_MIS can execute only one of them 
//notice: when change the condition compilation, retype 'make clean' and 'make' to work

#define UNMATCH_SINGLE_END
#define	SINGLE_END_NOEXECUTE

//#define UNMATCH_SINGLE_END_MIS
#define	CIR_JUMP

#ifdef	R_W_LOCK
int read_seq = 0;
pthread_rwlock_t rwlock;
#endif

#define	MAX_REDUCE_ANCHOR_NUM	1000

uint32_t cir_cnt = 0;
uint32_t extension_cnt = 0;
//uint32_t un_match_cnt = 0;
//uint32_t single_end_num = 0;

#define	QUAL_FILT
#define QUAL_FILT_SINGLE

#define QUAL_FILT_REP
#define	QUAL_FILT_LV
#define	QUAL_FILT_LV_MIS

#define	QUAL_FILT_LV_SINGLE

#define	QUAL_FILT_LV_OUT
#define	QUAL_FILT_SINGLE_OUT

#define	MAPPING_QUALITY

#ifdef	MAPPING_QUALITY
float*** mp_subs1 = NULL;
float*** mp_subs2 = NULL;
float* sub_t = NULL;
#endif

//uint32_t last_cir_n = 0;
//uint32_t mapped_n = 0;
//#define	LAST_CIRCLE_NOPOUND															

#define	SPLIT_LV_SINGLE
#define	KSW_ALN
#define	KSW_ALN_PAIR

//#define	KSW_ALN_PAIR_DEBUG

//#define	SINGLE_PAIR

#ifdef	KSW_ALN_PAIR
int** op_dm_kl2 = NULL;
int** op_dm_kr2 = NULL;
int** ops_dm_kl2 = NULL;
int** ops_dm_kr2 = NULL;
#endif

#ifdef	KSW_ALN
int8_t* mat = NULL;
char** ali_ref_seq2 = NULL;
char** read_char2 = NULL;
int** op_dm_kl1 = NULL;
int** op_dm_kr1 = NULL;
int** ops_dm_kl1 = NULL;
int** ops_dm_kr1 = NULL;
#endif

//#define	SINGLE_DEBUG
//#define	POUND_DEBUG
//#define	TER_J
//#define	QUAL_ANCHOR
//#define	QUAL_PAIR

#define	LAST_CIR
#define	LAST_CIR_SINGLE
#define	POUND_MIS
#define	QUAL_PAIR_SINGLE
#define	OUTPUT_DEBUG

float var_avg = 0.05;
float var_max = 0.2;

#define ANCHOR_HASH_ALI

#ifdef	ANCHOR_HASH_ALI
uint8_t k_anchor = 10;
uint8_t anchor_seed_d = 5;
uint8_t k_anchor_back = 4;
uint8_t anchor_back_mask = 0Xf;
uint8_t k_anchor_re = 22;
uint8_t k_anchor_b = 20;

typedef struct read_hs{
	uint32_t des;
	uint16_t off_set;
}read_h;

typedef struct anchor_seeds{
	int read_left_off;
	int read_right_off;
	uint32_t ref_left_off;
	uint32_t ref_right_off;
	uint16_t seed_length;
}anchor_seed;

uint16_t*** anchor_hash = NULL;
uint8_t*** anchor_array = NULL;
uint16_t*** anchor_point = NULL;
uint16_t*** anchor_pos = NULL;

read_h** rh = NULL;
anchor_seed** anchor_seed_buffer = NULL;

#endif


#ifdef	ALTER_DEBUG

typedef struct seed_length_arrays{
	uint16_t seed_length;
	uint16_t index;
}seed_length_array;

seed_length_array** seed_length_arr = NULL;

uint8_t* rep_go = NULL;

#endif

#define	SEED_FILTER_LENGTH
#define	SEED_FILTER_POS

#ifdef	SEED_FILTER_POS
uint16_t seed_filter_pos_num = 0;
uint16_t seed_filter_pos_num_single = 0;
#endif


#define	LV_CCIGAR
#define	CIGAR_OP

uint32_t lv_cnt = 0;
uint32_t lv_cnt_s = 0;
uint32_t single_lvcigar_t = 0;
uint32_t single_lvcigar_t1 = 0;
uint32_t single_lvcigar_t2 = 0;
uint32_t ssw_n = 0;

#define	PTHREAD_USE

#define	PAIR_RANDOM
#define MERGE_B_R_VS
#define PAIR_B_R_V
#define HASH_KMER_READ_J
#define ALI_B_V_R
#define SPLIT_LV
#define ALI_OUT
#define	PAIR_TRAVERSE
#define ALT_ALL
#define	SINGLE_COV_FILT
#define	SSW_SCORE_FILT
#define	OUPUT_REPEAT
#define	OUPUT_SINGLE_REPEAT
#define	OUPUT_SINGLE_REPEAT2
#define	SINGLE_SCORE_FILT

//#define	DM_COPY_SINGLE
//#define	SINGLE_END_RANDOM

//#define	DM_COPY_REP
//#define LV_CIGAR
//#define	SEED_FILTER
//#define	PAIR_SEED_LENGTH_FILT

#define	MAX_READLENTH	2048
#define	MAX_READLEN	2049
#define SEED_LENGTH	28
#define SEED_HASH_LENGTH	16
#define SEED_KMER_LENGTH	16
#define UNI_SEQ_WRI_ARR	1024
#define MAX_POS_N	300

//#define	MAX_POS_R_N	300
#define	MAX_PTHREAD	32
#define MAX_DIFF_POS	10
#define MAX_COV_A	(((MAX_READLEN - 1) >> 6) + 1)
#define MAX_REF_SEQ_C	(((MAX_READLEN - 1) >> 5) + 3)
#define MAX_Q_NUM	32
#define MAX_LV_REFLEN	(MAX_READLEN + 64)
#define MAX_LV_K 32
#define MAX_LV_CIGAR	(MAX_READLEN + 64)
#define MAX_LV_CIGARCOM	(MAX_READLEN + 64)
#define MAX_EDIT_SCORE	(MAX_READLEN << 1)
#define MAX_OP_SCORE	MAX_READLEN
#define	MAX_OP_SCORE_P1	1024
#define	MAX_OP_SCORE_P2	1025
#define	CIGARMN	6

#define	CUS_SEED_SET	2000
//#define	CUS_MAX_OUTPUT_ALI	150
//#define CUS_MAX_OUTPUT_ALI2 (CUS_MAX_OUTPUT_ALI << 1)
#define CUS_MAX_OUTPUT_ALI2 2000
//1000
//500
//const uint8_t cus_ali_n = 20;


#ifdef	PAIR_RANDOM

#define	PR_COV_FILTER
#define	PAIR_RANDOM_SEED
#define	PR_SINGLE

#ifdef	PR_COV_FILTER
#define	COV_FIT_SI	5
#define	RAN_CIR	4096
#endif

#ifdef	PR_SINGLE
#define	MAX_SPM	0Xfe
#endif

#define	RANDOM_RANGE	2000
#define	RANDOM_RANGE_MAX (RANDOM_RANGE << 1)
#define MINKEY 0
#define MAXKEY 0Xffffffff
#define MAX_COL (RANDOM_RANGE + 1)
#define	PAIR_RAN_INTV	5

#endif

typedef struct seed_MEM
{
	uint64_t cov[MAX_COV_A];//33
	
#ifdef UNPIPATH_OFF_K20
	uint64_t uni_id;
#else
	uint32_t uni_id;
#endif

	uint32_t ref_pos_off;
	uint32_t ref_pos_off_r;
	uint32_t ref_pos_n;
	int s_r_o_l;
	int s_r_o_r;
	uint16_t length;
	uint16_t read_pos;
	uint8_t ui;
	uint8_t tid;
}seed_m;

typedef struct seed_X
{
	uint32_t seed_id_r;
	uint32_t ref_pos_n;
	uint32_t kmer_pos_uni;
	uint8_t read_off;
}seed_x;

typedef struct seed_PARRAY
{
	uint64_t cov[MAX_COV_A];

	uint32_t pos_start;
	uint32_t pos_n;
	uint32_t ref_pos_off;
	uint32_t ref_pos_off_r;

	int s_r_o_l;
	int s_r_o_r;

	uint16_t length;
	uint8_t ui;
}seed_pa;

typedef struct seed_PARRAY_SINGLE
{
	uint64_t cov[MAX_COV_A];

	uint32_t pos_start;
	uint32_t pos_n;
	uint32_t ref_pos_off;
	uint32_t ref_pos_off_r;

	int s_r_o_l;
	int s_r_o_r;

	uint16_t length;
	uint8_t ui;
	uint8_t rc;
	uint8_t pair_flag;
}seed_pa_single;

typedef struct seed_SETS
{
	uint64_t cov[MAX_COV_A];
	
#ifdef UNPIPATH_OFF_K20
	uint64_t seed_set;
#else
	uint32_t seed_set;
#endif

	uint32_t ref_pos_off;
	uint32_t ref_pos_off_r;

	int s_r_o_l;
	int s_r_o_r;

	uint8_t ui;
	uint8_t tid;

}seed_sets;


//for pthread seq input
#define	OUTPUT_ARR

typedef struct SEQIOIN
{
	char* read_seq1;
	char* read_seq2;

	char* qual1;
	char* qual2;

	char* name;
	char* name_other;
	
	uint16_t read_length1;
	uint16_t read_length2;

	uint16_t length_h1;
	uint16_t length_h2;
	uint16_t nm1;
	uint16_t nm2;

	uint16_t v_cnt;
	uint8_t flag1;
	uint8_t flag2;
	int chr_re;
	int chr_re1;
	int chr_re2;

	int64_t pos1;
	int64_t pos2;
	int qualc1;
	int qualc2;
	int64_t cross;
	char* seq1;
	char* seq2;
	char* cigar1;
	char* cigar2;
	uint16_t xa_n;
	uint16_t xa_n_p1;
	uint16_t xa_n_p2;	
	uint16_t xa_n1;
	uint16_t xa_n2;
	
#ifdef	FIX_SV	
	uint16_t xa_n_x1;
	uint16_t xa_n_x2;
	int* chr_res_s1;
	int* chr_res_s2;
#endif
	int* chr_res;
	char* xa_d1s;
	uint32_t* sam_pos1s;
	char* cigar_p1s[CUS_MAX_OUTPUT_ALI2];
	char* cigar_p2s[CUS_MAX_OUTPUT_ALI2];

	int* lv_re1s;
	char* xa_d2s;
	uint32_t* sam_pos2s;

	int* lv_re2s;

	uint8_t* chr_res1;
	char* xa_ds1;
	uint32_t* sam_poss1;
	uint8_t* lv_res1;

	uint8_t* chr_res2;
	char* xa_ds2;
	uint32_t* sam_poss2;
	uint8_t* lv_res2;

}seq_io;

#ifdef	FIX_SA
char** read_rev_buffer_1 = NULL;
#endif

seq_io* seqio = NULL;
char** qual1_buffer = NULL;
char** qual2_buffer = NULL;
char** read_rev_buffer = NULL;
char** pr_cigar1_buffer = NULL;
char** pr_cigar2_buffer = NULL;
int** chr_res_buffer = NULL;

#ifdef	FIX_SV
int** chr_res_buffer1 = NULL;
int** chr_res_buffer2 = NULL;
#endif

char** xa_d1s_buffer = NULL;
char** xa_d2s_buffer = NULL;
uint32_t** sam_pos1s_buffer = NULL;
uint32_t** sam_pos2s_buffer = NULL;
char*** cigar_p1s_buffer = NULL;
char*** cigar_p2s_buffer = NULL;
int** lv_re1s_buffer = NULL;
int** lv_re2s_buffer = NULL;
uint8_t** chr_res1_buffer = NULL;
char** xa_ds1_buffer = NULL;
uint32_t** sam_poss1_buffer = NULL;
uint8_t** lv_res1_buffer = NULL;
uint8_t** chr_res2_buffer = NULL;
char** xa_ds2_buffer = NULL;
uint32_t** sam_poss2_buffer = NULL;
uint8_t** lv_res2_buffer = NULL;
//end for pthread seq input

const uint8_t ali_exl = 32;
const uint8_t max_extension_length = 32;
const char pr_dir[2] = {'+','-'};
const char ssw_cigar[3] = {'M', 'I', 'D'};

//for one circle
int nei_array[7];
uint8_t nei_n = 0;

uint16_t pair_ran_intvp = 0;
uint32_t insert_dis = 500;
int64_t devi = 50 * 3;
uint32_t seed_num = 0;
uint64_t new_seed_cnt = 0;
uint8_t uni_d = 22;//30
uint8_t k_r = 0;
uint8_t re_b = 0;
uint8_t re_bt = 0;
uint8_t re_2bt = 0;
FILE* fp_sam = NULL;
uint32_t* random_buffer = NULL;
uint32_t* seed_r_dup = NULL;

#ifdef	READN_RANDOM_SEED	
#define	RANDOM_RANGE_READN	1024
uint32_t* random_buffer_readn = NULL;
uint64_t readn_cnt = 0;
#endif
	
#ifdef	PTHREAD_USE
typedef struct{
	uint32_t v_cnt;
	uint32_t vs_cnt;
}cnt_re;

typedef struct THREAD_ALI{
	uint8_t tid;
	uint32_t seqn;
} thread_ali_t;

int64_t* g_low = NULL;
int64_t* r_low = NULL;
seed_m** seedm = NULL;
seed_m** seedu = NULL;
seed_sets** seedsets = NULL;
uint32_t** seed_set_off = NULL;

#ifdef UNPIPATH_OFF_K20
uint64_t** seed_set_pos[2][2];
uint64_t** mat_pos1 = NULL;
uint64_t** mat_pos2 = NULL;
uint64_t** op_vector_pos1 = NULL;
uint64_t** op_vector_pos2 = NULL;
uint64_t** ops_vector_pos1 = NULL;
uint64_t** ops_vector_pos2 = NULL;

#ifdef	REDUCE_ANCHOR
uint64_t** poses1 = NULL;
uint64_t** poses2 = NULL;
#endif
					
#else

uint32_t** seed_set_pos[2][2];
uint32_t** mat_pos1 = NULL;
uint32_t** mat_pos2 = NULL;
uint32_t** op_vector_pos1 = NULL;
uint32_t** op_vector_pos2 = NULL;
uint32_t** ops_vector_pos1 = NULL;
uint32_t** ops_vector_pos2 = NULL;

#ifdef	REDUCE_ANCHOR
uint32_t** poses1 = NULL;
uint32_t** poses2 = NULL;
#endif

#endif

seed_pa** seedpa1[2];
seed_pa** seedpa2[2];

#ifdef	REDUCE_ANCHOR
int** ls1 = NULL;
int** ls2 = NULL;
int** rs1 = NULL;
int** rs2 = NULL;
uint8_t** rcs1 = NULL;
uint8_t** rcs2 = NULL;
#endif

int** op_dm_l1 = NULL;
int** op_dm_r1 = NULL;
int** op_dm_l2 = NULL;
int** op_dm_r2 = NULL;
int** ops_dm_l1 = NULL;
int** ops_dm_r1 = NULL;
int** ops_dm_l2 = NULL;
int** ops_dm_r2 = NULL;
int** op_dm_ex1 = NULL;
int** op_dm_ex2 = NULL;
int** ops_dm_ex1 = NULL;
int** ops_dm_ex2 = NULL;
uint64_t*** op_vector_seq1 = NULL;
uint64_t*** op_vector_seq2 = NULL;
uint64_t*** ops_vector_seq1 = NULL;
uint64_t*** ops_vector_seq2 = NULL;

#ifdef ALT_ALL
int** chr_res = NULL;
uint32_t** sam_pos1s = NULL;
uint32_t** sam_pos2s = NULL;
char*** cigar_p1s = NULL;
char*** cigar_p2s = NULL;
char** xa_d1s = NULL;
char** xa_d2s = NULL;
int** lv_re1s = NULL;
int** lv_re2s = NULL;
#endif

uint8_t** op_rc = NULL;
uint8_t** ops_rc = NULL;
uint64_t** ref_seq_tmp1 = NULL;
uint64_t** ref_seq_tmp2 = NULL;

uint8_t** mat_rc = NULL;
uint32_t** seed_no1 = NULL;
uint32_t** seed_no2 = NULL;

#ifdef	PAIR_RANDOM

#ifdef UNPIPATH_OFF_K20
uint64_t*** seed_k_pos = NULL;
uint64_t** b = NULL;
uint64_t** seedposk = NULL;	
uint64_t** seedpos[2][2];
uint64_t** seed_single_pos[2];
#else
uint32_t*** seed_k_pos = NULL;
uint32_t** b = NULL;
uint32_t** seedposk = NULL;	
uint32_t** seedpos[2][2];
uint32_t** seed_single_pos[2];
#endif

uint8_t** seed_posf[2];

int** seed_single_ld[2];
int** seed_single_rd[2];
int** seed_single_dm[2];
uint64_t*** seed_single_refs[2];

int16_t** left_dis = NULL;
int16_t** right_dis = NULL;

#ifdef PR_SINGLE
uint8_t** seedpos_mis[2][2];
uint8_t** pr_chr_res1 = NULL;
uint32_t** pr_sam_pos1 = NULL;
char** pr_xa_d1 = NULL;
uint8_t** pr_lv_re1 = NULL;
uint8_t** pr_chr_res2 = NULL;
uint32_t** pr_sam_pos2 = NULL;
char** pr_xa_d2 = NULL;
uint8_t** pr_lv_re2 = NULL;
#endif


uint32_t** ls = NULL;

uint64_t*** read_bit_pr = NULL;
uint8_t*** read_val_pr = NULL;

uint32_t** kcol = NULL;
uint32_t* seed_posn[2];
uint16_t* seedn = NULL;

void seed_repetitive(uint8_t, uint16_t, uint16_t, uint16_t, uint16_t, uint16_t, uint32_t, uint32_t, cnt_re* , uint32_t, uint16_t, uint16_t, int16_t, int16_t, int16_t, int16_t, int16_t, int16_t, int16_t, int16_t);

#ifdef UNPIPATH_OFF_K20
void seed_repetitive_single64(uint64_t* , uint32_t* , uint64_t**, uint8_t, uint16_t );
cnt_re* pair_op_add64(uint64_t , uint64_t , int , int , int , int , int , uint8_t, int, int, uint8_t, uint32_t, uint32_t, cnt_re* );
uint64_t** seed_set_pos_single = NULL;	
#else
void seed_repetitive_single(uint64_t* , uint32_t* , uint32_t**, uint8_t, uint16_t );
cnt_re* pair_op_add(uint32_t , uint32_t , int , int , int , int , int , uint8_t, int, int, uint8_t, uint32_t, uint32_t, cnt_re* );
uint32_t** seed_set_pos_single = NULL;	
#endif

#endif

//for single mapping changing

uint32_t* set_pos_n_single = NULL;
uint32_t* spa_i_single = NULL;
seed_pa_single** seedpa1_single = NULL;

#ifdef	SINGLE_PAIR
uint16_t* cov_num_front_single = NULL;
uint16_t* cov_num_re_single = NULL;
#endif

int8_t** num = NULL;
int8_t** ref_num = NULL;

//thread global values
uint8_t* rc_cov_f = NULL;
uint8_t** cov_a_n = NULL;
uint8_t* s_uid_f = NULL;
uint8_t* cov_a_n_s = NULL;

uint8_t* min_mis = NULL;
uint8_t* de_m_p = NULL;
uint8_t* de_m_p_o = NULL;

uint32_t* mat_posi = NULL;
uint16_t* end_dis1 = NULL;
uint16_t* end_dis2 = NULL;
int16_t* seed_l = NULL;
uint64_t* low_mask = NULL;
uint64_t* low_mask1 = NULL;
uint64_t* low_mask2 = NULL;
uint32_t* max_lengtht = NULL;

#ifdef	REDUCE_ANCHOR
uint16_t** orders1 = NULL;
uint16_t** orders2 = NULL;
int** dms1 = NULL;
int** dms2 = NULL;
#endif

int* dm_op = NULL;
int* dm_ops = NULL;

uint8_t* max_mismatch = NULL;
uint8_t* max_mismatch_p = NULL;
uint16_t* max_mismatch1 = NULL;
uint16_t* max_mismatch2 = NULL;

uint16_t* max_mismatch1_single = NULL;
uint16_t* max_mismatch2_single = NULL;

#ifdef	PR_SINGLE
uint8_t* pr_o_f = NULL;
uint8_t* seed_re_r = NULL;
#endif

uint8_t** pos_add = NULL;
char** cigar_m = NULL;
char** cigar_m1 = NULL;
char** cigar_m2 = NULL;
uint16_t* pos_ren[2][2];
uint32_t* seedpos_misn[2][2];
uint16_t** sub_mask = NULL;
uint64_t** ex_d_mask = NULL;
uint16_t** sub_mask1 = NULL;
uint64_t** ex_d_mask1 = NULL;
uint16_t** sub_mask2 = NULL;
uint64_t** ex_d_mask2 = NULL;

char** ali_ref_seq = NULL;
char** read_char = NULL;

uint64_t** read_bit_1 = NULL;
uint64_t** read_bit_2 = NULL;
uint8_t** read_val_1 = NULL;
uint8_t** read_val_2 = NULL;

uint8_t*** read_val1 = NULL;
uint8_t*** read_val2 = NULL;
uint64_t read_bit1[MAX_PTHREAD][2][((MAX_READLEN - 1) >> 5) + 1];
uint64_t read_bit2[MAX_PTHREAD][2][((MAX_READLEN - 1) >> 5) + 1];

#ifdef	QUAL_FILT_LV
uint8_t*** qual_filt_lv1 = NULL;
uint8_t*** qual_filt_lv2 = NULL;
#endif


#ifdef	QUAL_FILT
uint64_t qual_filt1[MAX_PTHREAD][2][((MAX_READLEN - 1) >> 5) + 1];
uint64_t qual_filt2[MAX_PTHREAD][2][((MAX_READLEN - 1) >> 5) + 1];
#endif

float** qual_arr_single = NULL;

#ifdef	PR_COV_FILTER
uint8_t* cov_filt_f = NULL;
#endif

float single_random_filter_r = 0.04;

//end thread global values

#ifdef	R_W_LOCK
int seed_ali_core(int, uint8_t);
int seed_ali_core_single_end(int, uint8_t );
#else
int seed_ali_core(uint32_t , uint8_t );
int seed_ali_core_single_end(uint32_t , uint8_t );
#endif

void pair_sam_output(uint8_t, uint16_t, uint16_t, uint16_t, cnt_re*, uint32_t, uint16_t, uint16_t, int16_t, int16_t, int16_t, int16_t, int16_t, int16_t, int16_t, int16_t, uint8_t);


seed_pa_single* single_seed_reduction_core_single(seed_pa_single*, uint64_t*, uint8_t*, uint32_t**, uint8_t, uint16_t, uint16_t*, uint8_t );


#ifdef UNPIPATH_OFF_K20
seed_pa* single_seed_reduction_core_filter64(seed_pa* , uint64_t* , uint8_t* , uint32_t* , uint64_t** , uint16_t*, uint8_t, uint16_t, uint16_t* , uint16_t* );
void merge_pair_end64(uint32_t, uint32_t, seed_pa*, seed_pa*, uint64_t*, uint64_t*, uint8_t, uint8_t );

#else
seed_pa* single_seed_reduction_core_filter(seed_pa* , uint64_t* , uint8_t* , uint32_t* , uint32_t** , uint16_t*, uint8_t, uint16_t, uint16_t* , uint16_t* );
void merge_pair_end(uint32_t, uint32_t, seed_pa*, seed_pa*, uint32_t*, uint32_t*, uint8_t, uint8_t );

#endif

#ifdef	SINGLE_END_RANDOM
void seed_repetitive_single_end(uint64_t* , uint8_t , uint16_t );
#endif

#endif

#endif /* SEED_ALI_H_ */
