
#include <time.h>
#include <ctype.h>
#include <inttypes.h>

#include <omp.h>
#include <math.h>

//#include "global.h"
#include "binarys_qsort.h"
#include "LandauVishkin.h"
#include "seed_ali_p.h"
#include "bit_operation.h"
#include "ksw.h"

#ifdef PTHREAD_USE

#include <pthread.h>

#ifdef	R_W_LOCK
static void *worker(void *data)
{
	thread_ali_t *d = (thread_ali_t*)data;
	int _read_seq;

	while (1)
	{
		pthread_rwlock_wrlock(&rwlock);
		_read_seq = read_seq++;
		pthread_rwlock_unlock(&rwlock);

		if (_read_seq < d->seqn)
		{
			seed_ali_core(_read_seq, d->tid);
		}
		else break;
	}

	return 0;
}

static void *worker_single_end(void *data)
{
	thread_ali_t *d = (thread_ali_t*)data;

	int _read_seq;

	while (1)
	{
		pthread_rwlock_wrlock(&rwlock);
		_read_seq = read_seq++;
		pthread_rwlock_unlock(&rwlock);

		if (_read_seq < d->seqn)
		{
			seed_ali_core_single_end(_read_seq, d->tid);
		}
		else break;
	}

	return 0;
}

#else

static void *worker(void *data)
{
	thread_ali_t *d = (thread_ali_t*)data;
	seed_ali_core(d->seqn, d->tid);
	return 0;
}

static void *worker_single_end(void *data)
{
	thread_ali_t *d = (thread_ali_t*)data;
	seed_ali_core_single_end(d->seqn, d->tid);
	return 0;
}

#endif



#endif

int seed_ali()
{
	fprintf(stderr, "pair end reads mapping\n\n");

	double t = 0.00;
	clock_t start = 0, end = 0;

	uint8_t readlen_name = 200;
	uint16_t v_cnt_i = 0;
	//uint16_t read_bit_char = 0;
	uint32_t r_i = 0;
	uint32_t seqi = 0;
	uint32_t seqii = 0;
	uint32_t read_in = 0X10000;
	int32_t l, m, k1;//maybe changed here
	int64_t kr1 = 1;
	int64_t kr2 = 1;
	char h_chars[6];

#ifdef	PAIR_RANDOM_SEED
	uint32_t r_dup_i = 0;
	uint32_t pos_r_ir = 0;
	int pos_r_nr = 0;
#endif

	start = clock();

	load_index_file();

	end = clock();

	t = (double)(end - start) / CLOCKS_PER_SEC;

	fp1 = gzopen(read_fastq1, "r");
	seq1 = kseq_init(fp1);

	fp2 = gzopen(read_fastq2, "r");
	seq2 = kseq_init(fp2);

	if((fp1 == NULL) || (fp2 == NULL))
	{

		fprintf(stderr, "wrong input file route or name: %s %s\n", read_fastq1, read_fastq2);

		exit(1);
	}

	if(flag_std == 0)
	{
		fp_sam = fopen(sam_result, "w");
		if(fp_sam == NULL)
		{
			fprintf(stderr, "wrong output file route or name: %s\n", sam_result);
			exit(1);
		}
	}
	

	seed_l_max = seed_step * cir_fix_n;

	insert_dis = (upper_ins + floor_ins) >> 1;
	devi = (upper_ins - floor_ins) >> 1;
	pair_ran_intvp = seed_l_max;

#ifdef	SEED_FILTER_POS
	seed_filter_pos_num = seed_filter_pos_numn >> 2;
#endif


	//fprintf(stderr, "number of threads: %u\nkmer size: %u\ninsert size: %u\nstandard varaiance: %"PRId64"\nseed interval: %u\nthe max read length: %u\nmax lv error rate: %.2f\nmax pair error rate: %.2f\nthe minimum seed interval: %u\nthe number of seed circulation: %u\nlv rate on anchor: %.2f\nthe max rate of mismatch alignment on anchor: %.2f\n", thread_n, k_t, insert_dis, (devi / 3), pair_ran_intvp, readlen_max, lv_rate, max_pair_score_r, seed_step, cir_fix_n, lv_rate_anchor, mis_match_r_single);
	fprintf(stderr, "number of threads: %u\nkmer size: %u\nseed interval: %u\nthe max read length: %u\nmax lv error rate: %.2f\nmax pair error rate: %.2f\nthe minimum seed interval: %u\nthe number of seed circulation: %u\nlv rate on anchor: %.2f\nthe max rate of mismatch alignment on anchor: %.2f\n", thread_n, k_t, pair_ran_intvp, readlen_max, lv_rate, max_pair_score_r, seed_step, cir_fix_n, lv_rate_anchor, mis_match_r_single);

	fprintf(stderr, "local alignment score match: %d mismatch: %d  gap open: %d  gap extension: %d\n", mat_score, mis_score, gapo_score, gape_score);

	fprintf(stderr, "%lf seconds is used for loading index\n", t);

	if(local_ksw)	fprintf(stderr, "use local sw alignment\n");
	if(last_circle_rate)	fprintf(stderr, "use last circle end-to-end alignment and the max rate of alignment score is %.2f\n", last_circle_rate);
	if(!mgn_flag)	fprintf(stderr, "use the mode of multi-genomes\n");

#ifdef	LV_MP

	fprintf(stderr, "cal mapping quality on LV\n");

#endif

#ifdef	QUAL_ANCHOR

	fprintf(stderr, "use quality in anchor\n");

#endif
	if(flag_std)
		fprintf(stderr, "output alignments by stdout\n");
	
	k_r = k_t - k;
	re_b = 32 - k_t;
	re_bt = (re_b << 1);
	re_2bt = 64 - k_t;

	//should allocate at beginning
	seed_num = ((readlen_max - k_t) / seed_l_l) + 1;
	new_seed_cnt = seed_num * pos_n_max;

	uint16_t cigar_max_n = MAX_LV_CIGARCOM;

#ifdef	STA_SEED
	fp_read = fopen(file_read, "w");
	fp_seed = fopen(file_seed, "w");

#endif

	//allocation
	start = clock();

#ifdef	PTHREAD_USE

#ifdef	ANCHOR_HASH_ALI

	anchor_hash = (uint16_t*** )calloc(thread_n, sizeof(uint16_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		anchor_hash[r_i] = (uint16_t** )calloc(4, sizeof(uint16_t* ));
		for(m = 0; m < 4; m++)	anchor_hash[r_i][m] = (uint16_t* )calloc(1 << 16, 2);
	}

	anchor_array = (uint8_t*** )calloc(thread_n, sizeof(uint8_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		anchor_array[r_i] = (uint8_t** )calloc(4, sizeof(uint8_t* ));
		for(m = 0; m < 4; m++)	anchor_array[r_i][m] = (uint8_t* )calloc(readlen_max - k_anchor + 1, 1);
	}

	anchor_point = (uint16_t*** )calloc(thread_n, sizeof(uint16_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		anchor_point[r_i] = (uint16_t** )calloc(4, sizeof(uint16_t* ));
		for(m = 0; m < 4; m++)	anchor_point[r_i][m] = (uint16_t* )calloc(readlen_max - k_anchor + 1, 2);
	}
	anchor_pos = (uint16_t*** )calloc(thread_n, sizeof(uint16_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		anchor_pos[r_i] = (uint16_t** )calloc(4, sizeof(uint16_t* ));
		for(m = 0; m < 4; m++)	anchor_pos[r_i][m] = (uint16_t* )calloc(readlen_max - k_anchor + 1, 2);
	}

	rh = (read_h** )calloc(thread_n, sizeof(read_h* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		rh[r_i] = (read_h* )calloc(readlen_max - k_anchor + 1, sizeof(read_h));

	anchor_seed_buffer = (anchor_seed** )calloc(thread_n, sizeof(anchor_seed* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		anchor_seed_buffer[r_i] = (anchor_seed* )calloc((readlen_max - k_anchor + 1) * insert_dis, sizeof(anchor_seed));

#endif

	g_low = (int64_t* )calloc(thread_n, 8);
	r_low = (int64_t* )calloc(thread_n, 8);

	seedm = (seed_m** )calloc(thread_n, sizeof(seed_m* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedm[r_i] = (seed_m* )calloc(seed_num, sizeof(seed_m));

	seedu = (seed_m** )calloc(thread_n, sizeof(seed_m* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedu[r_i] = (seed_m* )calloc(seed_num, sizeof(seed_m));


	seedsets = (seed_sets** )calloc(thread_n, sizeof(seed_sets* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedsets[r_i] = (seed_sets* )calloc(new_seed_cnt, sizeof(seed_sets));

	seed_set_off = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_off[r_i] = (uint32_t* )calloc(new_seed_cnt, sizeof(uint32_t ));

#ifdef UNPIPATH_OFF_K20
	seed_set_pos[0][0] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos[0][0][r_i] = (uint64_t* )calloc(new_seed_cnt, 8);

	seed_set_pos[0][1] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos[0][1][r_i] = (uint64_t* )calloc(new_seed_cnt, 8);

	seed_set_pos[1][0] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos[1][0][r_i] = (uint64_t* )calloc(new_seed_cnt, 8);

	seed_set_pos[1][1] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos[1][1][r_i] = (uint64_t* )calloc(new_seed_cnt, 8);
#else
	seed_set_pos[0][0] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos[0][0][r_i] = (uint32_t* )calloc(new_seed_cnt, 4);

	seed_set_pos[0][1] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos[0][1][r_i] = (uint32_t* )calloc(new_seed_cnt, 4);

	seed_set_pos[1][0] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos[1][0][r_i] = (uint32_t* )calloc(new_seed_cnt, 4);

	seed_set_pos[1][1] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos[1][1][r_i] = (uint32_t* )calloc(new_seed_cnt, 4);
#endif

	seedpa1[0] = (seed_pa** )calloc(thread_n, sizeof(seed_pa* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpa1[0][r_i] = (seed_pa* )calloc(new_seed_cnt, sizeof(seed_pa ));

	seedpa1[1] = (seed_pa** )calloc(thread_n, sizeof(seed_pa* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpa1[1][r_i] = (seed_pa* )calloc(new_seed_cnt, sizeof(seed_pa ));

	seedpa2[0] = (seed_pa** )calloc(thread_n, sizeof(seed_pa* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpa2[0][r_i] = (seed_pa* )calloc(new_seed_cnt, sizeof(seed_pa ));

	seedpa2[1] = (seed_pa** )calloc(thread_n, sizeof(seed_pa* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpa2[1][r_i] = (seed_pa* )calloc(new_seed_cnt, sizeof(seed_pa ));


	fprintf(stderr, "Load seed reduction allocation\n");


#ifdef	ALTER_DEBUG
	seed_length_arr = (seed_length_array** )calloc(thread_n, sizeof(seed_length_array* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_length_arr[r_i] = (seed_length_array* )calloc(cus_max_output_ali, sizeof(seed_length_array));

	rep_go = (uint8_t* )calloc(thread_n, 1);

#endif

#ifdef UNPIPATH_OFF_K20

#ifdef	REDUCE_ANCHOR
	poses1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		poses1[r_i] = (uint64_t* )calloc(cus_max_output_ali << 1, 8);

	poses2 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		poses2[r_i] = (uint64_t* )calloc(cus_max_output_ali << 1, 8);

#endif

	op_vector_pos1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_vector_pos1[r_i] = (uint64_t* )calloc(cus_max_output_ali, 8);

	op_vector_pos2 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_vector_pos2[r_i] = (uint64_t* )calloc(cus_max_output_ali, 8);

	ops_vector_pos1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_vector_pos1[r_i] = (uint64_t* )calloc(cus_max_output_ali, 8);

	ops_vector_pos2 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_vector_pos2[r_i] = (uint64_t* )calloc(cus_max_output_ali, 8);

#else

#ifdef	REDUCE_ANCHOR
	poses1 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		poses1[r_i] = (uint32_t* )calloc(cus_max_output_ali << 1, 4);

	poses2 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		poses2[r_i] = (uint32_t* )calloc(cus_max_output_ali << 1, 4);

#endif

	op_vector_pos1 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_vector_pos1[r_i] = (uint32_t* )calloc(cus_max_output_ali, 4);

	op_vector_pos2 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_vector_pos2[r_i] = (uint32_t* )calloc(cus_max_output_ali, 4);

	ops_vector_pos1 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_vector_pos1[r_i] = (uint32_t* )calloc(cus_max_output_ali, 4);

	ops_vector_pos2 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_vector_pos2[r_i] = (uint32_t* )calloc(cus_max_output_ali, 4);

#endif


#ifdef	REDUCE_ANCHOR
	ls1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ls1[r_i] = (int* )calloc(cus_max_output_ali << 1, 4);

	rs1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		rs1[r_i] = (int* )calloc(cus_max_output_ali << 1, 4);

	ls2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ls2[r_i] = (int* )calloc(cus_max_output_ali << 1, 4);

	rs2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		rs2[r_i] = (int* )calloc(cus_max_output_ali << 1, 4);

#endif

	op_dm_l1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_l1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_r1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_r1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_l2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_l2[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_r2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_r2[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_l1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_l1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_r1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_r1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_l2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_l2[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_r2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_r2[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_ex1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_ex1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_ex2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_ex2[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_ex1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_ex1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_ex2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_ex2[r_i] = (int* )calloc(cus_max_output_ali, 4);


	//for S
#ifdef QUALS_CHECK
	op_err = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_err[r_i] = (uint8_t* )calloc(cus_max_output_ali, 1);

	ops_err = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_err[r_i] = (uint8_t* )calloc(cus_max_output_ali, 1);
#endif

	op_vector_seq1 = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		op_vector_seq1[r_i] = (uint64_t** )calloc(cus_max_output_ali, sizeof(uint64_t* ));
		for(m = 0; m < cus_max_output_ali; m++)
			if((op_vector_seq1[r_i][m] = (uint64_t* )calloc(MAX_REF_SEQ_C, 8)) == NULL)
				exit(1);
	}

	op_vector_seq2 = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		op_vector_seq2[r_i] = (uint64_t** )calloc(cus_max_output_ali, sizeof(uint64_t* ));
		for(m = 0; m < cus_max_output_ali; m++)
			if((op_vector_seq2[r_i][m] = (uint64_t* )calloc(MAX_REF_SEQ_C, 8)) == NULL)
				exit(1);
	}

	ops_vector_seq1 = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		ops_vector_seq1[r_i] = (uint64_t** )calloc(cus_max_output_ali, sizeof(uint64_t* ));
		for(m = 0; m < cus_max_output_ali; m++)
			if((ops_vector_seq1[r_i][m] = (uint64_t* )calloc(MAX_REF_SEQ_C, 8)) == NULL)
				exit(1);
	}

	ops_vector_seq2 = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		ops_vector_seq2[r_i] = (uint64_t** )calloc(cus_max_output_ali, sizeof(uint64_t* ));
		for(m = 0; m < cus_max_output_ali; m++)
			if((ops_vector_seq2[r_i][m] = (uint64_t* )calloc(MAX_REF_SEQ_C, 8)) == NULL)
				exit(1);
	}


	fprintf(stderr, "Load alignment allocation\n");


#ifdef ALT_ALL

	chr_res = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		chr_res[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

	sam_pos1s = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		if((sam_pos1s[r_i] = (uint32_t* )calloc(CUS_MAX_OUTPUT_ALI2, sizeof(uint32_t ))) == NULL)
			exit(1);

	sam_pos2s = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		if((sam_pos2s[r_i] = (uint32_t* )calloc(CUS_MAX_OUTPUT_ALI2, sizeof(uint32_t ))) == NULL)
			exit(1);

	cigar_p1s = (char*** )calloc(thread_n, sizeof(char** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		cigar_p1s[r_i] = (char** )calloc(CUS_MAX_OUTPUT_ALI2, sizeof(char* ));
		for(m = 0; m < CUS_MAX_OUTPUT_ALI2; m++)
			if((cigar_p1s[r_i][m] = (char* )calloc(cigar_max_n, 1)) == NULL)
				exit(1);
	}

	cigar_p2s = (char*** )calloc(thread_n, sizeof(char** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		cigar_p2s[r_i] = (char** )calloc(CUS_MAX_OUTPUT_ALI2, sizeof(char* ));
		for(m = 0; m < CUS_MAX_OUTPUT_ALI2; m++)
			if((cigar_p2s[r_i][m] = (char* )calloc(cigar_max_n, 1)) == NULL)
				exit(1);
	}

	xa_d1s = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		xa_d1s[r_i] = (char* )calloc(CUS_MAX_OUTPUT_ALI2, 1);

	xa_d2s = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		xa_d2s[r_i] = (char* )calloc(CUS_MAX_OUTPUT_ALI2, 1);

	lv_re1s = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		lv_re1s[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, sizeof(int));

	lv_re2s = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		lv_re2s[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, sizeof(int));

#endif


	fprintf(stderr, "Load output allocation\n");


#ifdef	REDUCE_ANCHOR
	rcs1 = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		rcs1[r_i] = (uint8_t* )calloc(cus_max_output_ali << 1, 1);

	rcs2 = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		rcs2[r_i] = (uint8_t* )calloc(cus_max_output_ali << 1, 1);

#endif

	op_rc = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_rc[r_i] = (uint8_t* )calloc(cus_max_output_ali, 1);

	ops_rc = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_rc[r_i] = (uint8_t* )calloc(cus_max_output_ali, 1);

	ref_seq_tmp1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ref_seq_tmp1[r_i] = (uint64_t* )calloc(MAX_REF_SEQ_C, 8);

	ref_seq_tmp2 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ref_seq_tmp2[r_i] = (uint64_t* )calloc(MAX_REF_SEQ_C, 8);

#ifdef UNPIPATH_OFF_K20
	mat_pos1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		mat_pos1[r_i] = (uint64_t* )calloc(CUS_SEED_SET, 8);

	mat_pos2 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		mat_pos2[r_i] = (uint64_t* )calloc(CUS_SEED_SET, 8);
#else
	mat_pos1 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		mat_pos1[r_i] = (uint32_t* )calloc(CUS_SEED_SET, 4);

	mat_pos2 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		mat_pos2[r_i] = (uint32_t* )calloc(CUS_SEED_SET, 4);
#endif

	mat_rc = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		mat_rc[r_i] = (uint8_t* )calloc(CUS_SEED_SET, 1);

	seed_no1 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_no1[r_i] = (uint32_t* )calloc(CUS_SEED_SET, 4);

	seed_no2 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_no2[r_i] = (uint32_t* )calloc(CUS_SEED_SET, 4);


	fprintf(stderr, "Load alignment more allocation\n");


#ifdef	PAIR_RANDOM

#ifdef UNPIPATH_OFF_K20
	seedpos[0][0] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos[0][0][r_i] = (uint64_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 8);

	seedpos[0][1] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos[0][1][r_i] = (uint64_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 8);

	seedpos[1][0] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos[1][0][r_i] = (uint64_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 8);

	seedpos[1][1] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos[1][1][r_i] = (uint64_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 8);

	seed_single_pos[0] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_pos[0][r_i] = (uint64_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 8);

	seed_single_pos[1] = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_pos[1][r_i] = (uint64_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 8);

	seed_k_pos = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ));

	for(r_i = 0; r_i < thread_n; r_i++)   //2000
	{
		seed_k_pos[r_i] = (uint64_t** )calloc(readlen_max / pair_ran_intvp, sizeof(uint64_t* ));
		for(m = 0; m < readlen_max / pair_ran_intvp; m++)
			if((seed_k_pos[r_i][m] = (uint64_t* )calloc(RANDOM_RANGE, 8)) == NULL)
				exit(1);
	}

	seedposk = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedposk[r_i] = (uint64_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 8);

#else
	seedpos[0][0] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos[0][0][r_i] = (uint32_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seedpos[0][1] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos[0][1][r_i] = (uint32_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seedpos[1][0] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos[1][0][r_i] = (uint32_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seedpos[1][1] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos[1][1][r_i] = (uint32_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seed_single_pos[0] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_pos[0][r_i] = (uint32_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seed_single_pos[1] = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_pos[1][r_i] = (uint32_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seed_k_pos = (uint32_t*** )calloc(thread_n, sizeof(uint32_t** ));

	for(r_i = 0; r_i < thread_n; r_i++)   //2000
	{
		seed_k_pos[r_i] = (uint32_t** )calloc(readlen_max / pair_ran_intvp, sizeof(uint32_t* ));
		for(m = 0; m < readlen_max / pair_ran_intvp; m++)
			if((seed_k_pos[r_i][m] = (uint32_t* )calloc(RANDOM_RANGE, 4)) == NULL)
				exit(1);
	}

	seedposk = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedposk[r_i] = (uint32_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

#endif

	seed_posf[0] = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_posf[0][r_i] = (uint8_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 1);

	seed_posf[1] = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_posf[1][r_i] = (uint8_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 1);


	seed_single_ld[0] = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_ld[0][r_i] = (int* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seed_single_ld[1] = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_ld[1][r_i] = (int* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seed_single_rd[0] = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_rd[0][r_i] = (int* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seed_single_rd[1] = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_rd[1][r_i] = (int* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seed_single_dm[0] = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_dm[0][r_i] = (int* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	seed_single_dm[1] = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_single_dm[1][r_i] = (int* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 4);

	if((seed_single_refs[0] = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ))) == NULL)
		exit(1);
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if((seed_single_refs[0][r_i] = (uint64_t** )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, sizeof(uint64_t* ))) == NULL)
			exit(1);
		for(m = 0; m < (readlen_max / pair_ran_intvp) * RANDOM_RANGE; m++)
		{
			if((seed_single_refs[0][r_i][m] = (uint64_t* )calloc((((readlen_max - 1) >> 5) + 3), 8)) == NULL)
				exit(1);
		}
	}

	if((seed_single_refs[1] = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ))) == NULL)
		exit(1);
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if((seed_single_refs[1][r_i] = (uint64_t** )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, sizeof(uint64_t* ))) == NULL)
			exit(1);
		for(m = 0; m < (readlen_max / pair_ran_intvp) * RANDOM_RANGE; m++)
		{
			if((seed_single_refs[1][r_i][m] = (uint64_t* )calloc((((readlen_max - 1) >> 5) + 3), 8)) == NULL)
				exit(1);
		}

	}

#ifdef	PR_SINGLE
	seedpos_mis[0][0] = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos_mis[0][0][r_i] = (uint8_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 1);

	seedpos_mis[0][1] = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos_mis[0][1][r_i] = (uint8_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 1);

	seedpos_mis[1][0] = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos_mis[1][0][r_i] = (uint8_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 1);

	seedpos_mis[1][1] = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpos_mis[1][1][r_i] = (uint8_t* )calloc((readlen_max / pair_ran_intvp) * RANDOM_RANGE, 1);

	pr_chr_res1 = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pr_chr_res1[r_i] = (uint8_t* )calloc(pr_single_outputn, 1);

	pr_sam_pos1 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pr_sam_pos1[r_i] = (uint32_t* )calloc(pr_single_outputn, 4);

	pr_xa_d1 = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pr_xa_d1[r_i] = (char* )calloc(pr_single_outputn, 1);

	pr_lv_re1 = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pr_lv_re1[r_i] = (uint8_t* )calloc(pr_single_outputn, 1);

	pr_chr_res2 = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pr_chr_res2[r_i] = (uint8_t* )calloc(pr_single_outputn, 1);

	pr_sam_pos2 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pr_sam_pos2[r_i] = (uint32_t* )calloc(pr_single_outputn, 4);

	pr_xa_d2 = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pr_xa_d2[r_i] = (char* )calloc(pr_single_outputn, 1);

	pr_lv_re2 = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pr_lv_re2[r_i] = (uint8_t* )calloc(pr_single_outputn, 1);
#endif


	fprintf(stderr, "Load pair random allocation\n");


#ifdef	PAIR_RANDOM_SEED
	seed_r_dup = (uint32_t* )calloc(RANDOM_RANGE, 4);
	srand((unsigned)time(0));
	pos_r_nr = RANDOM_RANGE;

	r_dup_i = 0;
	for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE_MAX; pos_r_ir++)
		if((rand() % (RANDOM_RANGE_MAX - pos_r_ir)) < pos_r_nr)
		{
			seed_r_dup[r_dup_i++] = pos_r_ir;
			pos_r_nr--;
		}

	random_buffer = (uint32_t* )malloc((RAN_CIR) << 2);
	for(pos_r_ir = 0; pos_r_ir < RAN_CIR; pos_r_ir++)
		random_buffer[pos_r_ir] = (uint32_t )rand();
#endif


#endif


#ifdef	READN_RANDOM_SEED
	random_buffer_readn = (uint32_t* )calloc(RANDOM_RANGE_READN, 4);
	for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE_READN; pos_r_ir++)
		random_buffer_readn[pos_r_ir] = (uint32_t )rand();
#endif

	//thread gloabl value
	rc_cov_f = (uint8_t* )calloc(thread_n, 1);
	cov_a_n = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		cov_a_n[r_i] = (uint8_t* )calloc(2,1);

	cov_a_n_s = (uint8_t* )calloc(thread_n, 1);
	s_uid_f = (uint8_t* )calloc(thread_n, 1);

	min_mis = (uint8_t* )calloc(thread_n, 1);
	de_m_p = (uint8_t* )calloc(thread_n, 1);
	de_m_p_o = (uint8_t* )calloc(thread_n, 1);
	seed_l = (int16_t* )calloc(thread_n, 2);//cannot be 0
	max_mismatch1 = (uint16_t* )calloc(thread_n, 2);
	max_mismatch2 = (uint16_t* )calloc(thread_n, 2);

	max_mismatch1_single = (uint16_t* )calloc(thread_n, 2);
	max_mismatch2_single = (uint16_t* )calloc(thread_n, 2);

#ifdef	PR_SINGLE
	pr_o_f = (uint8_t* )calloc(thread_n, 1);
	seed_re_r = (uint8_t* )calloc(thread_n, 1);
#endif

	end_dis1 = (uint16_t* )calloc(thread_n, 2);
	end_dis2 = (uint16_t* )calloc(thread_n, 2);
	mat_posi = (uint32_t* )calloc(thread_n, 4);
#ifdef	PAIR_SEED_LENGTH_FILT
	max_lengtht = (uint32_t* )calloc(thread_n, 4);
#endif

#ifdef	REDUCE_ANCHOR
	orders1 = (uint16_t** )calloc(thread_n, sizeof(uint16_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		orders1[r_i] = (uint16_t* )calloc(cus_max_output_ali << 1, 2);

	orders2 = (uint16_t** )calloc(thread_n, sizeof(uint16_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		orders2[r_i] = (uint16_t* )calloc(cus_max_output_ali << 1, 2);

	dms1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		dms1[r_i] = (int* )calloc(cus_max_output_ali << 1, 4);

	dms2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		dms2[r_i] = (int* )calloc(cus_max_output_ali << 1, 4);
#endif

	dm_op = (int* )calloc(thread_n, 4);
	dm_ops = (int* )calloc(thread_n, 4);
	low_mask1 = (uint64_t* )calloc(thread_n, 8);
	low_mask2 = (uint64_t* )calloc(thread_n, 8);

	pos_add = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pos_add[r_i] = (uint8_t* )calloc(pos_n_max, 1);

	cigar_m1 = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		cigar_m1[r_i] = (char* )calloc(CIGARMN, 1);

	cigar_m2 = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		cigar_m2[r_i] = (char* )calloc(CIGARMN, 1);

	ali_ref_seq = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ali_ref_seq[r_i] = (char* )calloc(readlen_max + 64, 1);

	read_char = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		read_char[r_i] = (char* )calloc(readlen_max + 32, 1);

	pos_ren[0][0] = (uint16_t* )calloc(thread_n, sizeof(uint16_t ));
	pos_ren[0][1] = (uint16_t* )calloc(thread_n, sizeof(uint16_t ));
	pos_ren[1][0] = (uint16_t* )calloc(thread_n, sizeof(uint16_t ));
	pos_ren[1][1] = (uint16_t* )calloc(thread_n, sizeof(uint16_t ));

	sub_mask1 = (uint16_t** )calloc(thread_n, sizeof(uint16_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		sub_mask1[r_i] = (uint16_t* )calloc(33, 2);

	sub_mask2 = (uint16_t** )calloc(thread_n, sizeof(uint16_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		sub_mask2[r_i] = (uint16_t* )calloc(33, 2);

	seedpos_misn[0][0] = (uint32_t* )calloc(thread_n, sizeof(uint32_t ));
	seedpos_misn[0][1] = (uint32_t* )calloc(thread_n, sizeof(uint32_t ));
	seedpos_misn[1][0] = (uint32_t* )calloc(thread_n, sizeof(uint32_t ));
	seedpos_misn[1][1] = (uint32_t* )calloc(thread_n, sizeof(uint32_t ));

	ex_d_mask1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ex_d_mask1[r_i] = (uint64_t* )calloc(33, 8);

	ex_d_mask2 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ex_d_mask2[r_i] = (uint64_t* )calloc(33, 8);

#ifdef	PAIR_RANDOM

#ifdef	UNPIPATH_OFF_K20
	b = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		b[r_i] = (uint64_t* )calloc((readlen_max / pair_ran_intvp) + 1, 8);
#else
	b = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		b[r_i] = (uint32_t* )calloc((readlen_max / pair_ran_intvp) + 1, 4);
#endif
	ls = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ls[r_i] = (uint32_t* )calloc(readlen_max / pair_ran_intvp, 4);

	read_bit_pr = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
		read_bit_pr[r_i] = (uint64_t** )calloc(4, sizeof(uint64_t* ));

	read_val_pr = (uint8_t*** )calloc(thread_n, sizeof(uint8_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
		read_val_pr[r_i] = (uint8_t** )calloc(4, sizeof(uint8_t* ));

	kcol = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		kcol[r_i] = (uint32_t* )calloc(readlen_max / pair_ran_intvp, 4);

	seed_posn[0] = (uint32_t* )calloc(thread_n, 4);
	seed_posn[1] = (uint32_t* )calloc(thread_n, 4);

	seedn = (uint16_t* )calloc(thread_n, 2);

#endif

	read_bit_1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	read_bit_2 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	read_val_1 = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	read_val_2 = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));

	read_val1 = (uint8_t*** )calloc(thread_n, sizeof(uint8_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		read_val1[r_i] = (uint8_t** )calloc(2, sizeof(uint8_t* ));
		read_val1[r_i][0] = (uint8_t* )calloc(readlen_max, 1);
		read_val1[r_i][1] = (uint8_t* )calloc(readlen_max, 1);
	}

	read_val2 = (uint8_t*** )calloc(thread_n, sizeof(uint8_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		read_val2[r_i] = (uint8_t** )calloc(2, sizeof(uint8_t* ));
		read_val2[r_i][0] = (uint8_t* )calloc(readlen_max, 1);
		read_val2[r_i][1] = (uint8_t* )calloc(readlen_max, 1);
	}

#ifdef	QUAL_FILT_LV
	qual_filt_lv1 = (uint8_t*** )calloc(thread_n, sizeof(uint8_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		qual_filt_lv1[r_i] = (uint8_t** )calloc(2, sizeof(uint8_t* ));
		qual_filt_lv1[r_i][0] = (uint8_t* )calloc(readlen_max, 1);
		qual_filt_lv1[r_i][1] = (uint8_t* )calloc(readlen_max, 1);
	}
	qual_filt_lv2 = (uint8_t*** )calloc(thread_n, sizeof(uint8_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		qual_filt_lv2[r_i] = (uint8_t** )calloc(2, sizeof(uint8_t* ));
		qual_filt_lv2[r_i][0] = (uint8_t* )calloc(readlen_max, 1);
		qual_filt_lv2[r_i][1] = (uint8_t* )calloc(readlen_max, 1);
	}

#endif

#ifdef	MAPPING_QUALITY
	mp_subs1 = (float*** )calloc(thread_n, sizeof(float** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		mp_subs1[r_i] = (float** )calloc(2, sizeof(float* ));
		mp_subs1[r_i][0] = (float* )calloc(readlen_max, sizeof(float));
		mp_subs1[r_i][1] = (float* )calloc(readlen_max, sizeof(float));
	}

	mp_subs2 = (float*** )calloc(thread_n, sizeof(float** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		mp_subs2[r_i] = (float** )calloc(2, sizeof(float* ));
		mp_subs2[r_i][0] = (float* )calloc(readlen_max, sizeof(float));
		mp_subs2[r_i][1] = (float* )calloc(readlen_max, sizeof(float));
	}

	sub_t = (float* )calloc(thread_n, sizeof(float ));
#endif

#ifdef	PR_COV_FILTER
	cov_filt_f = (uint8_t* )calloc(thread_n, 1);
#endif

	//end thread gloabl value
#endif


	fprintf(stderr, "successfully allocate memory\n");

	InitiateLVCompute(L, thread_n);

	//for pthread read input

	fprintf(stderr, "Load seq input memory\n");

	seqio = (seq_io* )calloc(read_in, sizeof(seq_io));

	char** read_seq1_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		read_seq1_buffer[r_i] = (char* )calloc(readlen_max, 1);

	char** read_seq2_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		read_seq2_buffer[r_i] = (char* )calloc(readlen_max, 1);

	char** name_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		name_buffer[r_i] = (char* )calloc(readlen_name, 1);

	char** name_buffer_other= (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		name_buffer_other[r_i] = (char* )calloc(readlen_name, 1);

	qual1_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		qual1_buffer[r_i] = (char* )calloc(readlen_max, 1);

	qual2_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		qual2_buffer[r_i] = (char* )calloc(readlen_max, 1);

#ifdef	FIX_SA
	read_rev_buffer_1 = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		read_rev_buffer_1[r_i] = (char* )calloc(readlen_max, 1);
#endif


#ifdef	OUTPUT_ARR
	read_rev_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		read_rev_buffer[r_i] = (char* )calloc(readlen_max, 1);

	pr_cigar1_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		pr_cigar1_buffer[r_i] = (char* )calloc(cigar_max_n, 1);

	pr_cigar2_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		pr_cigar2_buffer[r_i] = (char* )calloc(cigar_max_n, 1);

	chr_res_buffer = (int** )calloc(read_in, sizeof(int* ));
	for(r_i = 0; r_i < read_in; r_i++)
		chr_res_buffer[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

#ifdef	FIX_SV
	chr_res_buffer1 = (int** )calloc(read_in, sizeof(int* ));
	for(r_i = 0; r_i < read_in; r_i++)
		chr_res_buffer1[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

	chr_res_buffer2 = (int** )calloc(read_in, sizeof(int* ));
	for(r_i = 0; r_i < read_in; r_i++)
		chr_res_buffer2[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, 4);
#endif

	xa_d1s_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		xa_d1s_buffer[r_i] = (char* )calloc(CUS_MAX_OUTPUT_ALI2, 1);

	xa_d2s_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		xa_d2s_buffer[r_i] = (char* )calloc(CUS_MAX_OUTPUT_ALI2, 1);

	sam_pos1s_buffer = (uint32_t** )calloc(read_in, sizeof(uint32_t* ));
	for(r_i = 0; r_i < read_in; r_i++)
		sam_pos1s_buffer[r_i] = (uint32_t* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

	sam_pos2s_buffer = (uint32_t** )calloc(read_in, sizeof(uint32_t* ));
	for(r_i = 0; r_i < read_in; r_i++)
		sam_pos2s_buffer[r_i] = (uint32_t* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

	lv_re1s_buffer = (int** )calloc(read_in, sizeof(int* ));
	for(r_i = 0; r_i < read_in; r_i++)
		lv_re1s_buffer[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

	lv_re2s_buffer = (int** )calloc(read_in, sizeof(int* ));
	for(r_i = 0; r_i < read_in; r_i++)
		lv_re2s_buffer[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

#ifdef	PR_SINGLE
	chr_res1_buffer = (uint8_t** )calloc(read_in, sizeof(uint8_t* ));
	for(r_i = 0; r_i < read_in; r_i++)
		chr_res1_buffer[r_i] = (uint8_t* )calloc(pr_single_outputn, 1);

	chr_res2_buffer = (uint8_t** )calloc(read_in, sizeof(uint8_t* ));
	for(r_i = 0; r_i < read_in; r_i++)
		chr_res2_buffer[r_i] = (uint8_t* )calloc(pr_single_outputn, 1);

	xa_ds1_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		xa_ds1_buffer[r_i] = (char* )calloc(pr_single_outputn, 1);

	xa_ds2_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		xa_ds2_buffer[r_i] = (char* )calloc(pr_single_outputn, 1);

	lv_res1_buffer = (uint8_t** )calloc(read_in, sizeof(uint8_t* ));
	for(r_i = 0; r_i < read_in; r_i++)
		lv_res1_buffer[r_i] = (uint8_t* )calloc(pr_single_outputn, 1);

	lv_res2_buffer = (uint8_t** )calloc(read_in, sizeof(uint8_t* ));
	for(r_i = 0; r_i < read_in; r_i++)
		lv_res2_buffer[r_i] = (uint8_t* )calloc(pr_single_outputn, 1);

	sam_poss1_buffer = (uint32_t** )calloc(read_in, sizeof(uint32_t* ));
	for(r_i = 0; r_i < read_in; r_i++)
		sam_poss1_buffer[r_i] = (uint32_t* )calloc(pr_single_outputn, 4);

	sam_poss2_buffer = (uint32_t** )calloc(read_in, sizeof(uint32_t* ));
	for(r_i = 0; r_i < read_in; r_i++)
		sam_poss2_buffer[r_i] = (uint32_t* )calloc(pr_single_outputn, 4);
#endif

#endif

#ifdef KSW_ALN_PAIR
	ali_ref_seq2 = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ali_ref_seq2[r_i] = (char* )calloc(readlen_max + 64, 1);

	read_char2 = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		read_char2[r_i] = (char* )calloc(readlen_max + 32, 1);

	op_dm_kl1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_kl1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_kr1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_kr1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_kl2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_kl2[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_kr2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_kr2[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_kl1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_kl1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_kr1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_kr1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_kl2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_kl2[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_kr2 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_kr2[r_i] = (int* )calloc(cus_max_output_ali, 4);

	mat = (int8_t*)calloc(25, sizeof(int8_t));
	for (l = k1 = 0; l < 4; ++l)
	{
		for (m = 0; m < 4; ++m) mat[k1++] = l == m ? mat_score : -mis_score;	/* weight_match : -weight_mismatch */
		mat[k1++] = -1; // ambiguous base
	}
	for (m = 0; m < 5; ++m) mat[k1++] = 0;
#endif

#ifdef	PAIR_RANDOM

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		read_bit_pr[r_i][0] = read_bit1[r_i][0];
		read_bit_pr[r_i][1] = read_bit2[r_i][1];
		read_bit_pr[r_i][2] = read_bit2[r_i][0];
		read_bit_pr[r_i][3] = read_bit1[r_i][1];

		read_val_pr[r_i][0] = read_val1[r_i][0];
		read_val_pr[r_i][1] = read_val2[r_i][1];
		read_val_pr[r_i][2] = read_val2[r_i][0];
		read_val_pr[r_i][3] = read_val1[r_i][1];
	}

#endif

	end = clock();

	t = (double)(end - start) / CLOCKS_PER_SEC;

	fprintf(stderr, "%lf seconds is used for allocating memory\n", t);

	//write .sam head into result file

	if(flag_std)
		for(r_i = 1; r_i < chr_file_n; r_i++)
			printf("@SQ\tSN:%s\tLN:%u\n", chr_names[r_i], (uint32_t )(chr_end_n[r_i] - chr_end_n[r_i - 1]));
	else
		for(r_i = 1; r_i < chr_file_n; r_i++)
			fprintf(fp_sam, "@SQ\tSN:%s\tLN:%u\n", chr_names[r_i], (uint32_t )(chr_end_n[r_i] - chr_end_n[r_i - 1]));


#ifdef	RD_FILED

	if(flag_std)
	{
		printf("@RG\tID:L1\tPU:SC_1_10\tLB:SC_1\tSM:NA12878\n");
		printf("@RG\tID:L2\tPU:SC_2_12\tLB:SC_2\tSM:NA12878\n");
	}else{
		fprintf(fp_sam, "@RG\tID:L1\tPU:SC_1_10\tLB:SC_1\tSM:NA12878\n");
		fprintf(fp_sam, "@RG\tID:L2\tPU:SC_2_12\tLB:SC_2\tSM:NA12878\n");
	}

#endif

	fprintf(stderr, "begin reading fastq pair-end reads and doing alignment\n");
	//fflush(stdout);

#ifdef	R_W_LOCK
	pthread_rwlock_init(&rwlock, NULL);
	int seqii_i = 0;
#endif

	double dtime = omp_get_wtime(); //value in seconds

	int read_length_tmp1 = 0;
	int read_length_tmp2 = 0;

	char tmp_pos_ch;

#ifdef	FIX_SA
	uint16_t xa_n_p1_tmp = 0;
	uint16_t sam_seq_i = 0;
	char sam_seq1[MAX_READLEN] = {};
	uint16_t flag_tmp = 0;
	char* seq_tmp = NULL;
	char* read_tmp = NULL;
	char* qual_tmp = NULL;
#endif

	while ((kr1 > 0) && (kr2 > 0))
	{
		for(seqii = 0; (seqii < read_in) && (((kr1 = kseq_read(seq1)) > 0) && ((kr2 = kseq_read(seq2)) > 0)); seqii++)   //((kr1 > 0) && (kr2 > 0))
		{
			strncpy(name_buffer[seqii], seq1->name.s, seq1->name.l);
			strncpy(name_buffer_other[seqii], seq2->name.s, seq2->name.l);

			if((seq1->name.s)[seq1->name.l - 2] == '/')
				name_buffer[seqii][seq1->name.l - 2] = '\0';
			else	name_buffer[seqii][seq1->name.l] = '\0';

			if((seq2->name.s)[seq2->name.l - 2] == '/')
				name_buffer_other[seqii][seq2->name.l - 2] = '\0';
			else	name_buffer_other[seqii][seq2->name.l] = '\0';

			if(seq1->seq.l > readlen_max)
			{
				seqio[seqii].read_length1 = readlen_max;
				seqio[seqii].length_h1 = seq1->seq.l - readlen_max;

			}
			else
			{
				seqio[seqii].read_length1 = seq1->seq.l;
				seqio[seqii].length_h1 = 0;
			}

			if(seq2->seq.l > readlen_max)
			{
				seqio[seqii].read_length2 = readlen_max;
				seqio[seqii].length_h2 = seq2->seq.l - readlen_max;
			}
			else
			{
				seqio[seqii].read_length2 = seq2->seq.l;
				seqio[seqii].length_h2 = 0;
			}

			read_length_tmp1 = seqio[seqii].read_length1;
			read_length_tmp2 = seqio[seqii].read_length2;

			strncpy(read_seq1_buffer[seqii], seq1->seq.s, read_length_tmp1);
			strncpy(read_seq2_buffer[seqii], seq2->seq.s, read_length_tmp2);

			read_seq1_buffer[seqii][read_length_tmp1] = '\0';
			read_seq2_buffer[seqii][read_length_tmp2] = '\0';

			strncpy(qual1_buffer[seqii], seq1->qual.s, read_length_tmp1);
			strncpy(qual2_buffer[seqii], seq2->qual.s, read_length_tmp2);

			qual1_buffer[seqii][read_length_tmp1] = '\0';
			qual2_buffer[seqii][read_length_tmp2] = '\0';

			seqio[seqii].read_seq1 = read_seq1_buffer[seqii];
			seqio[seqii].read_seq2 = read_seq2_buffer[seqii];

			seqio[seqii].name = name_buffer[seqii];
			seqio[seqii].name_other = name_buffer_other[seqii];

			seqio[seqii].qual1 = qual1_buffer[seqii];
			seqio[seqii].qual2 = qual2_buffer[seqii];
			
		}

		if(thread_n <= 1)
		{
#ifdef	R_W_LOCK
			for(seqii_i = 0; seqii_i < seqii; seqii_i++)
				seed_ali_core(seqii_i, 0);
#else
			seed_ali_core(seqii, 0);
#endif
		}
		else
		{
			pthread_t* tid;
			thread_ali_t* data;
			pthread_attr_t attr;

			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_ali_t* )calloc(thread_n, sizeof(thread_ali_t ));
			tid = (pthread_t* )calloc(thread_n, sizeof(pthread_t ));
#ifdef	R_W_LOCK
			read_seq = 0;
#endif
			for(r_i = 0; r_i < thread_n; ++r_i)
			{
				data[r_i].tid = r_i;
				data[r_i].seqn = seqii;
				pthread_create(&tid[r_i], &attr, worker, data + r_i);
			}
			for(r_i = 0; r_i < thread_n; ++r_i) pthread_join(tid[r_i], 0);

			free(data);
			free(tid);
		}


		if(flag_std)
		{
			//output
#ifdef	OUTPUT_ARR

			for(seqi = 0; seqi < seqii; seqi++)
			{
				if((seqio[seqi].flag1 == 77) && (seqio[seqi].flag2 == 141))	tmp_pos_ch = '*';
				else	tmp_pos_ch = '=';

				if(seqio[seqi].length_h1 == 0)
				{
					if(seqio[seqi].chr_re1 == seqio[seqi].chr_re2)
					{
						printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t%"PRId64"\t%s\t%s",
						       seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re1], seqio[seqi].pos1,
						       seqio[seqi].qualc1, seqio[seqi].cigar1, tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].cross, seqio[seqi].seq1,
						       seqio[seqi].qual1
						      );
					}
					else
					{
						printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t0\t%s\t%s",
						       seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re1], seqio[seqi].pos1,
						       seqio[seqi].qualc1, seqio[seqi].cigar1, tmp_pos_ch, seqio[seqi].pos1, seqio[seqi].seq1,
						       seqio[seqi].qual1
						      );
					}


					if(seqio[seqi].xa_n_p1 > 0)
					{
						printf("\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p1; v_cnt_i++)
							printf("%s,%c%u,%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d1s[v_cnt_i], seqio[seqi].sam_pos1s[v_cnt_i], seqio[seqi].cigar_p1s[v_cnt_i], seqio[seqi].lv_re1s[v_cnt_i]);			}

#ifdef	PR_SINGLE
					if(seqio[seqi].xa_n1 > 0)
					{
						if((seqio[seqi].xa_n_p1 == 0) && (seqio[seqi].xa_n1 > 0))	printf("\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n1; v_cnt_i++)
							printf("%s,%c%u,%uM,%u;", chr_names[seqio[seqi].chr_res1[v_cnt_i]], seqio[seqi].xa_ds1[v_cnt_i], seqio[seqi].sam_poss1[v_cnt_i], seqio[seqi].read_length1, seqio[seqi].lv_res1[v_cnt_i]);

					}
#endif

#ifdef	FIX_SA
					if(seqio[seqi].xa_n_x1 > 0)
					{

						xa_n_p1_tmp = seqio[seqi].xa_n_p1;
						printf("\tSA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x1; v_cnt_i++)
						{
							printf("%s,%u,%c,%s,%u,%d;", chr_names[seqio[seqi].chr_res_s1[v_cnt_i]], seqio[seqi].sam_pos1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].xa_d1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].cigar_p1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].qualc1, seqio[seqi].lv_re1s[v_cnt_i + xa_n_p1_tmp]);

						}

					}
#endif

#ifdef	RD_FILED

					printf("\tNM:i:%u\tRG:Z:L1\n", seqio[seqi].nm1);

#else

					printf("\tNM:i:%u\n", seqio[seqi].nm1);


#endif
				}
				else
				{
					sprintf(h_chars, "%uH", seqio[seqi].length_h1);
					if(seqio[seqi].chr_re1 == seqio[seqi].chr_re2)
					{
						printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s%s\t%c\t%"PRId64"\t%"PRId64"\t%s\t%s",
						       seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re1], seqio[seqi].pos1,
						       seqio[seqi].qualc1, seqio[seqi].cigar1, h_chars, tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].cross, seqio[seqi].seq1,
						       seqio[seqi].qual1
						      );

					}
					else
					{
						printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s%s\t%c\t%"PRId64"\t0\t%s\t%s",
						       seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re1], seqio[seqi].pos1,
						       seqio[seqi].qualc1, seqio[seqi].cigar1, h_chars, tmp_pos_ch, seqio[seqi].pos1, seqio[seqi].seq1,
						       seqio[seqi].qual1
						      );

					}


					if(seqio[seqi].xa_n_p1 > 0)
					{
						printf("\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p1; v_cnt_i++)
						{
							printf("%s,%c%u,%s%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d1s[v_cnt_i], seqio[seqi].sam_pos1s[v_cnt_i], seqio[seqi].cigar_p1s[v_cnt_i], h_chars, seqio[seqi].lv_re1s[v_cnt_i]);
						}

					}
					printf("NM:i:%u\t", seqio[seqi].nm1);

#ifdef	PR_SINGLE
					if(seqio[seqi].xa_n1 > 0)
					{
						if((seqio[seqi].xa_n_p1 == 0) && (seqio[seqi].xa_n1 > 0))	printf("\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n1; v_cnt_i++)
							printf("%s,%c%u,%uM,%u;", chr_names[seqio[seqi].chr_res1[v_cnt_i]], seqio[seqi].xa_ds1[v_cnt_i], seqio[seqi].sam_poss1[v_cnt_i], seqio[seqi].read_length1, seqio[seqi].lv_res1[v_cnt_i]);

					}
#endif

#ifdef	FIX_SA
					if(seqio[seqi].xa_n_x1 > 0)
					{
						xa_n_p1_tmp = seqio[seqi].xa_n_p1;

						printf("\tSA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x1; v_cnt_i++)
						{
							printf("%s,%u,%c,%s,%u,%d;", chr_names[seqio[seqi].chr_res_s1[v_cnt_i]], seqio[seqi].sam_pos1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].xa_d1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].cigar_p1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].qualc1, seqio[seqi].lv_re1s[v_cnt_i + xa_n_p1_tmp]);
						}
					}
#endif

#ifdef	RD_FILED

					printf("\tNM:i:%u\tRG:Z:L1\n", seqio[seqi].nm1);

#else

					printf("\tNM:i:%u\n", seqio[seqi].nm1);

#endif
				}

#ifdef FIX_SV

#ifndef	FIX_SA

				for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x1; v_cnt_i++)
				{
					if(seqio[seqi].xa_d1s[v_cnt_i + seqio[seqi].xa_n_p1] == '+')
					{
						flag_tmp = 2145;
						seq_tmp = seqio[seqi].read_seq1;
					}
					else
					{
						flag_tmp = 2129;
						read_tmp = seqio[seqi].read_seq1;
						for(sam_seq_i = 0; sam_seq_i < seqio[seqi].read_length1; sam_seq_i++)
							sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[read_tmp[sam_seq_i]] ^ 0X3];
						sam_seq1[sam_seq_i] = '\0';
						strrev1(sam_seq1);
						seq_tmp = sam_seq1;
					}
					qual_tmp = seqio[seqi].qual1;

					printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t%d\t%s\t%s\tNM:i:%u\tRG:Z:L1\n",
					       seqio[seqi].name, flag_tmp, chr_names[seqio[seqi].chr_res_s1[v_cnt_i]], seqio[seqi].sam_pos1s[v_cnt_i + seqio[seqi].xa_n_p1],
					       seqio[seqi].qualc1, seqio[seqi].cigar_p1s[v_cnt_i + seqio[seqi].xa_n_p1], tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].pos2 - seqio[seqi].sam_pos1s[v_cnt_i + seqio[seqi].xa_n_p1], seq_tmp,
					       qual_tmp, seqio[seqi].lv_re1s[v_cnt_i + seqio[seqi].xa_n_p1]
					      );

				}
#endif
				
#endif

				if(seqio[seqi].length_h2 == 0)
				{
					if(seqio[seqi].chr_re1 == seqio[seqi].chr_re2)
					{
						printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t%"PRId64"\t%s\t%s",
						       seqio[seqi].name_other, seqio[seqi].flag2, chr_names[seqio[seqi].chr_re2], seqio[seqi].pos2,
						       seqio[seqi].qualc2, seqio[seqi].cigar2, tmp_pos_ch, seqio[seqi].pos1, -seqio[seqi].cross, seqio[seqi].seq2,
						       seqio[seqi].qual2
						      );
					}
					else
					{
						printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t0\t%s\t%s",
						       seqio[seqi].name_other, seqio[seqi].flag2, chr_names[seqio[seqi].chr_re2], seqio[seqi].pos2,
						       seqio[seqi].qualc2, seqio[seqi].cigar2, tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].seq2,
						       seqio[seqi].qual2
						      );
					}


					if(seqio[seqi].xa_n_p2 > 0)
					{
						printf("\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p2; v_cnt_i++)
							printf("%s,%c%u,%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d2s[v_cnt_i], seqio[seqi].sam_pos2s[v_cnt_i], seqio[seqi].cigar_p2s[v_cnt_i], seqio[seqi].lv_re2s[v_cnt_i]);

					}
#ifdef	PR_SINGLE
					if(seqio[seqi].xa_n2 > 0)
					{
						if((seqio[seqi].xa_n_p2 == 0) && (seqio[seqi].xa_n2 > 0))	printf("\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n2; v_cnt_i++)
							printf("%s,%c%u,%uM,%u;", chr_names[seqio[seqi].chr_res2[v_cnt_i]], seqio[seqi].xa_ds2[v_cnt_i], seqio[seqi].sam_poss2[v_cnt_i], seqio[seqi].read_length2, seqio[seqi].lv_res2[v_cnt_i]);

					}
#endif

#ifdef	FIX_SA
					if(seqio[seqi].xa_n_x2 > 0)
					{
						xa_n_p1_tmp = seqio[seqi].xa_n_p2;
						printf("\tSA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x2; v_cnt_i++)
						{
							printf("%s,%u,%c,%s,%u,%d;", chr_names[seqio[seqi].chr_res_s2[v_cnt_i]], seqio[seqi].sam_pos2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].xa_d2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].cigar_p2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].qualc2, seqio[seqi].lv_re2s[v_cnt_i + xa_n_p1_tmp]);

						}

					}
#endif

#ifdef	RD_FILED

					printf("\tNM:i:%u\tRG:Z:L2\n", seqio[seqi].nm2);
#else

					printf("\tNM:i:%u\n", seqio[seqi].nm2);

#endif

				}
				else
				{
					sprintf(h_chars, "%uH", seqio[seqi].length_h2);
					if(seqio[seqi].chr_re1 == seqio[seqi].chr_re2)
					{
						printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s%s\t%c\t%"PRId64"\t%"PRId64"\t%s\t%s",
						       seqio[seqi].name_other, seqio[seqi].flag2, chr_names[seqio[seqi].chr_re2], seqio[seqi].pos2,
						       seqio[seqi].qualc2, seqio[seqi].cigar2, h_chars, tmp_pos_ch, seqio[seqi].pos1, -seqio[seqi].cross, seqio[seqi].seq2,
						       seqio[seqi].qual2
						      );
					}
					else
					{
						printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s%s\t%c\t%"PRId64"\t0\t%s\t%s",
						       seqio[seqi].name_other, seqio[seqi].flag2, chr_names[seqio[seqi].chr_re2], seqio[seqi].pos2,
						       seqio[seqi].qualc2, seqio[seqi].cigar2, h_chars, tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].seq2,
						       seqio[seqi].qual2
						      );
					}


					if(seqio[seqi].xa_n_p2 > 0)
					{
						printf( "\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p2; v_cnt_i++)
						{
							printf("%s,%c%u,%s%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d2s[v_cnt_i], seqio[seqi].sam_pos2s[v_cnt_i], seqio[seqi].cigar_p2s[v_cnt_i], h_chars, seqio[seqi].lv_re2s[v_cnt_i]);
						}
					}
#ifdef	PR_SINGLE
					if(seqio[seqi].xa_n2 > 0)
					{
						if((seqio[seqi].xa_n_p2 == 0) && (seqio[seqi].xa_n2 > 0))	printf("\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n2; v_cnt_i++)
							printf("%s,%c%u,%uM,%u;", chr_names[seqio[seqi].chr_res2[v_cnt_i]], seqio[seqi].xa_ds2[v_cnt_i], seqio[seqi].sam_poss2[v_cnt_i], seqio[seqi].read_length2, seqio[seqi].lv_res2[v_cnt_i]);

					}
#endif

#ifdef	FIX_SA

					if(seqio[seqi].xa_n_x2 > 0)
					{
						xa_n_p1_tmp = seqio[seqi].xa_n_p2;
						printf("\tSA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x2; v_cnt_i++)
						{
							printf("%s,%u,%c,%s,%u,%d;", chr_names[seqio[seqi].chr_res_s2[v_cnt_i]], seqio[seqi].sam_pos2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].xa_d2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].cigar_p2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].qualc2, seqio[seqi].lv_re2s[v_cnt_i + xa_n_p1_tmp]);

						}
					}
#endif

#ifdef	RD_FILED

					printf("\tNM:i:%u\tRG:Z:L2\n", seqio[seqi].nm2);
#else

					printf("\tNM:i:%u\n", seqio[seqi].nm2);

#endif
				}


#ifdef FIX_SV

#ifndef	FIX_SA
				
				for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x2; v_cnt_i++)
				{
					if(seqio[seqi].xa_d2s[v_cnt_i + seqio[seqi].xa_n_p2] == '+')
					{
						flag_tmp = 2209;
						seq_tmp = seqio[seqi].read_seq2;
					}
					else
					{
						flag_tmp = 2193;
						read_tmp = seqio[seqi].read_seq2;
						for(sam_seq_i = 0; sam_seq_i < seqio[seqi].read_length2; sam_seq_i++)
							sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[read_tmp[sam_seq_i]] ^ 0X3];
						sam_seq1[sam_seq_i] = '\0';
						strrev1(sam_seq1);
						seq_tmp = sam_seq1;
					}

					qual_tmp = seqio[seqi].qual2;

					//seqio[seqi].read_seq1
					//read_rev_buffer[seqi]
					printf("%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t%d\t%s\t%s\tNM:i:%u\tRG:Z:L2\n",
					       seqio[seqi].name_other, flag_tmp, chr_names[seqio[seqi].chr_res_s2[v_cnt_i]], seqio[seqi].sam_pos2s[v_cnt_i + seqio[seqi].xa_n_p2],
					       seqio[seqi].qualc2, seqio[seqi].cigar_p2s[v_cnt_i + seqio[seqi].xa_n_p2], tmp_pos_ch, seqio[seqi].sam_pos2s[v_cnt_i + seqio[seqi].xa_n_p2], seqio[seqi].sam_pos2s[v_cnt_i + seqio[seqi].xa_n_p2] - seqio[seqi].pos1, seq_tmp,
					       qual_tmp, seqio[seqi].lv_re2s[v_cnt_i + seqio[seqi].xa_n_p2]
					      );
				}
				
#endif

#endif
			}
#endif
		}
		else
		{
			//output
#ifdef	OUTPUT_ARR

			for(seqi = 0; seqi < seqii; seqi++)
			{
				if((seqio[seqi].flag1 == 77) && (seqio[seqi].flag2 == 141))	tmp_pos_ch = '*';
				else	tmp_pos_ch = '=';

				if(seqio[seqi].length_h1 == 0)
				{
					if(seqio[seqi].chr_re1 == seqio[seqi].chr_re2)
					{
						fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t%"PRId64"\t%s\t%s",
						        seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re1], seqio[seqi].pos1,
						        seqio[seqi].qualc1, seqio[seqi].cigar1, tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].cross, seqio[seqi].seq1,
						        seqio[seqi].qual1
						       );

					}
					else
					{

						fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t0\t%s\t%s",
						        seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re1], seqio[seqi].pos1,
						        seqio[seqi].qualc1, seqio[seqi].cigar1, tmp_pos_ch, seqio[seqi].pos1, seqio[seqi].seq1,
						        seqio[seqi].qual1
						       );

					}


					if(seqio[seqi].xa_n_p1 > 0)
					{
						fprintf(fp_sam, "\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p1; v_cnt_i++)
							fprintf(fp_sam, "%s,%c%u,%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d1s[v_cnt_i], seqio[seqi].sam_pos1s[v_cnt_i], seqio[seqi].cigar_p1s[v_cnt_i], seqio[seqi].lv_re1s[v_cnt_i]);

					}

#ifdef	PR_SINGLE
					if(seqio[seqi].xa_n1 > 0)
					{

						if((seqio[seqi].xa_n_p1 == 0) && (seqio[seqi].xa_n1 > 0))	fprintf(fp_sam, "\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n1; v_cnt_i++)
							fprintf(fp_sam, "%s,%c%u,%uM,%u;", chr_names[seqio[seqi].chr_res1[v_cnt_i]], seqio[seqi].xa_ds1[v_cnt_i], seqio[seqi].sam_poss1[v_cnt_i], seqio[seqi].read_length1, seqio[seqi].lv_res1[v_cnt_i]);

					}
#endif

#ifdef	FIX_SA
					if(seqio[seqi].xa_n_x1 > 0)
					{

						xa_n_p1_tmp = seqio[seqi].xa_n_p1;

						fprintf(fp_sam, "\tSA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x1; v_cnt_i++)
						{
							fprintf(fp_sam, "%s,%u,%c,%s,%u,%d;", chr_names[seqio[seqi].chr_res_s1[v_cnt_i]], seqio[seqi].sam_pos1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].xa_d1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].cigar_p1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].qualc1, seqio[seqi].lv_re1s[v_cnt_i + xa_n_p1_tmp]);

						}
					}
#endif

#ifdef	RD_FILED

					fprintf(fp_sam, "\tNM:i:%u\tRG:Z:L1\n", seqio[seqi].nm1);

#else

					fprintf(fp_sam, "\tNM:i:%u\n", seqio[seqi].nm1);

#endif
				}
				else
				{
					sprintf(h_chars, "%uH", seqio[seqi].length_h1);
					if(seqio[seqi].chr_re1 == seqio[seqi].chr_re2)
					{
						fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s%s\t%c\t%"PRId64"\t%"PRId64"\t%s\t%s",
						        seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re1], seqio[seqi].pos1,
						        seqio[seqi].qualc1, seqio[seqi].cigar1, h_chars, tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].cross, seqio[seqi].seq1,
						        seqio[seqi].qual1
						       );
					}
					else
					{
						fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s%s\t%c\t%"PRId64"\t0\t%s\t%s",
						        seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re1], seqio[seqi].pos1,
						        seqio[seqi].qualc1, seqio[seqi].cigar1, h_chars, tmp_pos_ch, seqio[seqi].pos1, seqio[seqi].seq1,
						        seqio[seqi].qual1
						       );
					}

					if(seqio[seqi].xa_n_p1 > 0)
					{
						fprintf(fp_sam, "\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p1; v_cnt_i++)
						{
							fprintf(fp_sam, "%s,%c%u,%s%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d1s[v_cnt_i], seqio[seqi].sam_pos1s[v_cnt_i], seqio[seqi].cigar_p1s[v_cnt_i], h_chars, seqio[seqi].lv_re1s[v_cnt_i]);
						}

					}

					fprintf(fp_sam, "NM:i:%u\t", seqio[seqi].nm1);

#ifdef	PR_SINGLE
					if(seqio[seqi].xa_n1 > 0)
					{

						if((seqio[seqi].xa_n_p1 == 0) && (seqio[seqi].xa_n1 > 0))	fprintf(fp_sam, "\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n1; v_cnt_i++)
							fprintf(fp_sam, "%s,%c%u,%uM,%u;", chr_names[seqio[seqi].chr_res1[v_cnt_i]], seqio[seqi].xa_ds1[v_cnt_i], seqio[seqi].sam_poss1[v_cnt_i], seqio[seqi].read_length1, seqio[seqi].lv_res1[v_cnt_i]);

					}
#endif

#ifdef	FIX_SA
					if(seqio[seqi].xa_n_x1 > 0)
					{
						xa_n_p1_tmp = seqio[seqi].xa_n_p1;

						fprintf(fp_sam, "\tSA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x1; v_cnt_i++)
						{
							fprintf(fp_sam, "%s,%u,%c,%s,%u,%d;", chr_names[seqio[seqi].chr_res_s1[v_cnt_i]], seqio[seqi].sam_pos1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].xa_d1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].cigar_p1s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].qualc1, seqio[seqi].lv_re1s[v_cnt_i + xa_n_p1_tmp]);

						}
					}
#endif

#ifdef	RD_FILED

					fprintf(fp_sam, "\tNM:i:%u\tRG:Z:L1\n", seqio[seqi].nm1);

#else

					fprintf(fp_sam, "\tNM:i:%u\n", seqio[seqi].nm1);

#endif
				}

#ifdef FIX_SV

#ifndef	FIX_SA
				
				for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x1; v_cnt_i++)
				{
					if(seqio[seqi].xa_d1s[v_cnt_i + seqio[seqi].xa_n_p1] == '+')
					{
						flag_tmp = 2145;
						seq_tmp = seqio[seqi].read_seq1;
					}
					else
					{
						flag_tmp = 2129;
						read_tmp = seqio[seqi].read_seq1;
						for(sam_seq_i = 0; sam_seq_i < seqio[seqi].read_length1; sam_seq_i++)
							sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[read_tmp[sam_seq_i]] ^ 0X3];
						sam_seq1[sam_seq_i] = '\0';
						strrev1(sam_seq1);
						seq_tmp = sam_seq1;
					}
					qual_tmp = seqio[seqi].qual1;

					fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t%d\t%s\t%s\tNM:i:%u\tRG:Z:L1\n",
					        seqio[seqi].name, flag_tmp, chr_names[seqio[seqi].chr_res_s1[v_cnt_i]], seqio[seqi].sam_pos1s[v_cnt_i + seqio[seqi].xa_n_p1],
					        seqio[seqi].qualc1, seqio[seqi].cigar_p1s[v_cnt_i + seqio[seqi].xa_n_p1], tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].pos2 - seqio[seqi].sam_pos1s[v_cnt_i + seqio[seqi].xa_n_p1], seq_tmp,
					        qual_tmp, seqio[seqi].lv_re1s[v_cnt_i + seqio[seqi].xa_n_p1]
					       );

				}
#endif
				
#endif

				if(seqio[seqi].length_h2 == 0)
				{
					if(seqio[seqi].chr_re1 == seqio[seqi].chr_re2)
					{
						fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t%"PRId64"\t%s\t%s",
						        seqio[seqi].name_other, seqio[seqi].flag2, chr_names[seqio[seqi].chr_re2], seqio[seqi].pos2,
						        seqio[seqi].qualc2, seqio[seqi].cigar2, tmp_pos_ch, seqio[seqi].pos1, -seqio[seqi].cross, seqio[seqi].seq2,
						        seqio[seqi].qual2
						       );
					}
					else
					{
						fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t0\t%s\t%s",
						        seqio[seqi].name_other, seqio[seqi].flag2, chr_names[seqio[seqi].chr_re2], seqio[seqi].pos2,
						        seqio[seqi].qualc2, seqio[seqi].cigar2, tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].seq2,
						        seqio[seqi].qual2
						       );
					}


					if(seqio[seqi].xa_n_p2 > 0)
					{
						fprintf(fp_sam, "\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p2; v_cnt_i++)
							fprintf(fp_sam, "%s,%c%u,%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d2s[v_cnt_i], seqio[seqi].sam_pos2s[v_cnt_i], seqio[seqi].cigar_p2s[v_cnt_i], seqio[seqi].lv_re2s[v_cnt_i]);
					}
#ifdef	PR_SINGLE
					if(seqio[seqi].xa_n2 > 0)
					{
						if((seqio[seqi].xa_n_p2 == 0) && (seqio[seqi].xa_n2 > 0))	fprintf(fp_sam, "\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n2; v_cnt_i++)
							fprintf(fp_sam, "%s,%c%u,%uM,%u;", chr_names[seqio[seqi].chr_res2[v_cnt_i]], seqio[seqi].xa_ds2[v_cnt_i], seqio[seqi].sam_poss2[v_cnt_i], seqio[seqi].read_length2, seqio[seqi].lv_res2[v_cnt_i]);

					}
#endif

#ifdef	FIX_SA
					if(seqio[seqi].xa_n_x2 > 0)
					{
						xa_n_p1_tmp = seqio[seqi].xa_n_p2;

						fprintf(fp_sam, "\tSA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x2; v_cnt_i++)
						{
							fprintf(fp_sam, "%s,%u,%c,%s,%u,%d;", chr_names[seqio[seqi].chr_res_s2[v_cnt_i]], seqio[seqi].sam_pos2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].xa_d2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].cigar_p2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].qualc2, seqio[seqi].lv_re2s[v_cnt_i + xa_n_p1_tmp]);

						}

					}
#endif

#ifdef	RD_FILED

					fprintf(fp_sam, "\tNM:i:%u\tRG:Z:L2\n", seqio[seqi].nm2);

#else

					fprintf(fp_sam, "\tNM:i:%u\n", seqio[seqi].nm2);

#endif

				}
				else
				{
					sprintf(h_chars, "%uH", seqio[seqi].length_h2);
					if(seqio[seqi].chr_re1 == seqio[seqi].chr_re2)
					{
						fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s%s\t%c\t%"PRId64"\t%"PRId64"\t%s\t%s",
						        seqio[seqi].name_other, seqio[seqi].flag2, chr_names[seqio[seqi].chr_re2], seqio[seqi].pos2,
						        seqio[seqi].qualc2, seqio[seqi].cigar2, h_chars, tmp_pos_ch, seqio[seqi].pos1, -seqio[seqi].cross, seqio[seqi].seq2,
						        seqio[seqi].qual2
						       );
					}
					else
					{
						fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s%s\t%c\t%"PRId64"\t0\t%s\t%s",
						        seqio[seqi].name_other, seqio[seqi].flag2, chr_names[seqio[seqi].chr_re2], seqio[seqi].pos2,
						        seqio[seqi].qualc2, seqio[seqi].cigar2, h_chars, tmp_pos_ch, seqio[seqi].pos2, seqio[seqi].seq2,
						        seqio[seqi].qual2
						       );
					}


					if(seqio[seqi].xa_n_p2 > 0)
					{
						fprintf(fp_sam, "\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p2; v_cnt_i++)
						{
							fprintf(fp_sam, "%s,%c%u,%s%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d2s[v_cnt_i], seqio[seqi].sam_pos2s[v_cnt_i], seqio[seqi].cigar_p2s[v_cnt_i], h_chars, seqio[seqi].lv_re2s[v_cnt_i]);
						}
					}
#ifdef	PR_SINGLE
					if(seqio[seqi].xa_n2 > 0)
					{
						if((seqio[seqi].xa_n_p2 == 0) && (seqio[seqi].xa_n2 > 0))	fprintf(fp_sam, "\tXA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n2; v_cnt_i++)
							fprintf(fp_sam, "%s,%c%u,%uM,%u;", chr_names[seqio[seqi].chr_res2[v_cnt_i]], seqio[seqi].xa_ds2[v_cnt_i], seqio[seqi].sam_poss2[v_cnt_i], seqio[seqi].read_length2, seqio[seqi].lv_res2[v_cnt_i]);

					}
#endif

#ifdef	FIX_SA

					if(seqio[seqi].xa_n_x2 > 0)
					{
						xa_n_p1_tmp = seqio[seqi].xa_n_p2;

						fprintf(fp_sam, "\tSA:Z:");
						for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x2; v_cnt_i++)
						{
							fprintf(fp_sam, "%s,%u,%c,%s,%u,%d;", chr_names[seqio[seqi].chr_res_s2[v_cnt_i]], seqio[seqi].sam_pos2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].xa_d2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].cigar_p2s[v_cnt_i + xa_n_p1_tmp], seqio[seqi].qualc2, seqio[seqi].lv_re2s[v_cnt_i + xa_n_p1_tmp]);

						}
					}
#endif

#ifdef	RD_FILED

					fprintf(fp_sam, "\tNM:i:%u\tRG:Z:L2\n", seqio[seqi].nm2);

#else

					fprintf(fp_sam, "\tNM:i:%u\n", seqio[seqi].nm2);

#endif
				}


#ifdef FIX_SV

#ifndef	FIX_SA
			
				for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_x2; v_cnt_i++)
				{
					if(seqio[seqi].xa_d2s[v_cnt_i + seqio[seqi].xa_n_p2] == '+')
					{
						flag_tmp = 2209;
						seq_tmp = seqio[seqi].read_seq2;
					}
					else
					{
						flag_tmp = 2193;
						read_tmp = seqio[seqi].read_seq2;
						for(sam_seq_i = 0; sam_seq_i < seqio[seqi].read_length2; sam_seq_i++)
							sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[read_tmp[sam_seq_i]] ^ 0X3];
						sam_seq1[sam_seq_i] = '\0';
						strrev1(sam_seq1);
						seq_tmp = sam_seq1;
					}

					qual_tmp = seqio[seqi].qual2;

					//seqio[seqi].read_seq1
					//read_rev_buffer[seqi]

					fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s\t%c\t%"PRId64"\t%d\t%s\t%s\tNM:i:%u\tRG:Z:L2\n",
					        seqio[seqi].name_other, flag_tmp, chr_names[seqio[seqi].chr_res_s2[v_cnt_i]], seqio[seqi].sam_pos2s[v_cnt_i + seqio[seqi].xa_n_p2],
					        seqio[seqi].qualc2, seqio[seqi].cigar_p2s[v_cnt_i + seqio[seqi].xa_n_p2], tmp_pos_ch, seqio[seqi].sam_pos2s[v_cnt_i + seqio[seqi].xa_n_p2], seqio[seqi].sam_pos2s[v_cnt_i + seqio[seqi].xa_n_p2] - seqio[seqi].pos1, seq_tmp,
					        qual_tmp, seqio[seqi].lv_re2s[v_cnt_i + seqio[seqi].xa_n_p2]
					       );

				}
#endif

#endif
			}

#endif
		}


#ifdef	POST_DEBUG
		printf("print end\n");
#endif

#ifdef	OUTPUT_ARR
		for(seqi = 0; seqi < seqii; seqi++)
		{
			if(seqio[seqi].v_cnt > 0)
			{
				for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p1; v_cnt_i++)
					free(seqio[seqi].cigar_p1s[v_cnt_i]);

				for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n_p2; v_cnt_i++)
					free(seqio[seqi].cigar_p2s[v_cnt_i]);
			}
		}
#endif

	}


	dtime = omp_get_wtime() - dtime;
	fprintf(stderr, "%lf seconds is used in mapping\n", dtime);
	fflush(stdout);

#ifdef	R_W_LOCK
	pthread_rwlock_destroy(&rwlock);
#endif

	if(flag_std == 0)
		fclose(fp_sam);

	kseq_destroy(seq1);
	gzclose(fp1);
	kseq_destroy(seq2);
	gzclose(fp2);

#ifdef	STA_SEED
	fclose(fp_read);
	fclose(fp_seed);
#endif

#ifdef	PTHREAD_USE
	fprintf(stderr, "free memory\n");

	if(g_low != NULL)	free(g_low);
	if(r_low != NULL)	free(r_low);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedm[r_i] != NULL)	free(seedm[r_i]);
	if(seedm != NULL)	free(seedm);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedu[r_i] != NULL)	free(seedu[r_i]);
	if(seedu != NULL)	free(seedu);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedsets[r_i] != NULL)	free(seedsets[r_i]);
	if(seedsets != NULL)	free(seedsets);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_set_off[r_i] != NULL)	free(seed_set_off[r_i]);
	if(seed_set_off != NULL)	free(seed_set_off);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_set_pos[0][0][r_i] != NULL)	free(seed_set_pos[0][0][r_i]);
	if(seed_set_pos[0][0] != NULL)	free(seed_set_pos[0][0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_set_pos[0][1][r_i] != NULL)	free(seed_set_pos[0][1][r_i]);
	if(seed_set_pos[0][1] != NULL)	free(seed_set_pos[0][1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_set_pos[1][0][r_i] != NULL)	free(seed_set_pos[1][0][r_i]);
	if(seed_set_pos[1][0] != NULL)	free(seed_set_pos[1][0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_set_pos[1][1][r_i] != NULL)	free(seed_set_pos[1][1][r_i]);
	if(seed_set_pos[1][1] != NULL)	free(seed_set_pos[1][1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpa1[0][r_i] != NULL)	free(seedpa1[0][r_i]);
	if(seedpa1[0] != NULL)	free(seedpa1[0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpa1[1][r_i] != NULL)	free(seedpa1[1][r_i]);
	if(seedpa1[1] != NULL)	free(seedpa1[1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpa2[0][r_i] != NULL)	free(seedpa2[0][r_i]);
	if(seedpa2[0] != NULL)	free(seedpa2[0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpa2[1][r_i] != NULL)	free(seedpa2[1][r_i]);
	if(seedpa2[1] != NULL)	free(seedpa2[1]);

#ifdef	ALTER_DEBUG
	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_length_arr[r_i] != NULL)	free(seed_length_arr[r_i]);
	if(seed_length_arr != NULL)	free(seed_length_arr);

	if(rep_go != NULL)	free(rep_go);
#endif



#ifdef	REDUCE_ANCHOR
	for(r_i = 0; r_i < thread_n; r_i++)
		if(poses1[r_i] != NULL)	free(poses1[r_i]);
	if(poses1 != NULL)	free(poses1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(poses2[r_i] != NULL)	free(poses2[r_i]);
	if(poses2 != NULL)	free(poses2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ls1[r_i] != NULL)	free(ls1[r_i]);
	if(ls1 != NULL)	free(ls1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ls2[r_i] != NULL)	free(ls2[r_i]);
	if(ls2 != NULL)	free(ls2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(rs1[r_i] != NULL)	free(rs1[r_i]);
	if(rs1 != NULL)	free(rs1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(rs2[r_i] != NULL)	free(rs2[r_i]);
	if(rs2 != NULL)	free(rs2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(rcs1[r_i] != NULL)	free(rcs1[r_i]);
	if(rcs1 != NULL)	free(rcs1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(rcs2[r_i] != NULL)	free(rcs2[r_i]);
	if(rcs2 != NULL)	free(rcs2);

#endif

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_vector_pos1[r_i] != NULL)	free(op_vector_pos1[r_i]);
	if(op_vector_pos1 != NULL)	free(op_vector_pos1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_vector_pos2[r_i] != NULL)	free(op_vector_pos2[r_i]);
	if(op_vector_pos2 != NULL)	free(op_vector_pos2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_vector_pos1[r_i] != NULL)	free(ops_vector_pos1[r_i]);
	if(ops_vector_pos1 != NULL)	free(ops_vector_pos1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_vector_pos2[r_i] != NULL)	free(ops_vector_pos2[r_i]);
	if(ops_vector_pos2 != NULL)	free(ops_vector_pos2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_l1[r_i] != NULL)	free(op_dm_l1[r_i]);
	if(op_dm_l1 != NULL)	free(op_dm_l1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_r1[r_i] != NULL)	free(op_dm_r1[r_i]);
	if(op_dm_r1 != NULL)	free(op_dm_r1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_l2[r_i] != NULL)	free(op_dm_l2[r_i]);
	if(op_dm_l2 != NULL)	free(op_dm_l2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_r2[r_i] != NULL)	free(op_dm_r2[r_i]);
	if(op_dm_r2 != NULL)	free(op_dm_r2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_l1[r_i] != NULL)	free(ops_dm_l1[r_i]);
	if(ops_dm_l1 != NULL)	free(ops_dm_l1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_r1[r_i] != NULL)	free(ops_dm_r1[r_i]);
	if(ops_dm_r1 != NULL)	free(ops_dm_r1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_l2[r_i] != NULL)	free(ops_dm_l2[r_i]);
	if(ops_dm_l2 != NULL)	free(ops_dm_l2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_r2[r_i] != NULL)	free(ops_dm_r2[r_i]);
	if(ops_dm_r2 != NULL)	free(ops_dm_r2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_ex1[r_i] != NULL)	free(op_dm_ex1[r_i]);
	if(op_dm_ex1 != NULL)	free(op_dm_ex1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_ex2[r_i] != NULL)	free(op_dm_ex2[r_i]);
	if(op_dm_ex2 != NULL)	free(op_dm_ex2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_ex1[r_i] != NULL)	free(ops_dm_ex1[r_i]);
	if(ops_dm_ex1 != NULL)	free(ops_dm_ex1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_ex2[r_i] != NULL)	free(ops_dm_ex2[r_i]);
	if(ops_dm_ex2 != NULL)	free(ops_dm_ex2);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < cus_max_output_ali; m++)
			if(op_vector_seq1[r_i][m] != NULL)	free(op_vector_seq1[r_i][m]);
		if(op_vector_seq1[r_i] != NULL)	free(op_vector_seq1[r_i]);
	}
	if(op_vector_seq1 != NULL)	free(op_vector_seq1);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < cus_max_output_ali; m++)
			if(op_vector_seq2[r_i][m] != NULL)	free(op_vector_seq2[r_i][m]);
		if(op_vector_seq2[r_i] != NULL)	free(op_vector_seq2[r_i]);
	}
	if(op_vector_seq2 != NULL)	free(op_vector_seq2);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < cus_max_output_ali; m++)
			if(ops_vector_seq1[r_i][m] != NULL)	free(ops_vector_seq1[r_i][m]);

		if(ops_vector_seq1[r_i] != NULL)	free(ops_vector_seq1[r_i]);
	}
	if(ops_vector_seq1 != NULL)	free(ops_vector_seq1);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < cus_max_output_ali; m++)
			if(ops_vector_seq2[r_i][m] != NULL)	free(ops_vector_seq2[r_i][m]);
		if(ops_vector_seq2[r_i] != NULL)	free(ops_vector_seq2[r_i]);
	}
	if(ops_vector_seq2 != NULL)	free(ops_vector_seq2);

#ifdef ALT_ALL

	for(r_i = 0; r_i < thread_n; r_i++)
		if(chr_res[r_i] != NULL)	free(chr_res[r_i]);
	if(chr_res != NULL)	free(chr_res);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(sam_pos1s[r_i] != NULL)	free(sam_pos1s[r_i]);
	if(sam_pos1s != NULL)	free(sam_pos1s);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(sam_pos2s[r_i] != NULL)	free(sam_pos2s[r_i]);
	if(sam_pos2s != NULL)	free(sam_pos2s);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < CUS_MAX_OUTPUT_ALI2; m++)
			if(cigar_p1s[r_i][m] != NULL)	free(cigar_p1s[r_i][m]);
		if(cigar_p1s[r_i] != NULL)	free(cigar_p1s[r_i]);
	}
	if(cigar_p1s != NULL)	free(cigar_p1s);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < CUS_MAX_OUTPUT_ALI2; m++)
			if(cigar_p2s[r_i][m] != NULL)	free(cigar_p2s[r_i][m]);
		if(cigar_p2s[r_i] != NULL)	free(cigar_p2s[r_i]);
	}
	if(cigar_p2s != NULL)	free(cigar_p2s);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(xa_d1s[r_i] != NULL)	free(xa_d1s[r_i]);
	if(xa_d1s != NULL)	free(xa_d1s);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(xa_d2s[r_i] != NULL)	free(xa_d2s[r_i]);
	if(xa_d2s != NULL)	free(xa_d2s);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(lv_re1s[r_i] != NULL)	free(lv_re1s[r_i]);
	if(lv_re1s != NULL)	free(lv_re1s);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(lv_re2s[r_i] != NULL)	free(lv_re2s[r_i]);
	if(lv_re2s != NULL)	free(lv_re2s);

#endif

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_rc[r_i] != NULL)	free(op_rc[r_i]);
	if(op_rc != NULL)	free(op_rc);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_rc[r_i] != NULL)	free(ops_rc[r_i]);
	if(ops_rc != NULL)	free(ops_rc);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ref_seq_tmp1[r_i] != NULL)	free(ref_seq_tmp1[r_i]);
	if(ref_seq_tmp1 != NULL)	free(ref_seq_tmp1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ref_seq_tmp2[r_i] != NULL)	free(ref_seq_tmp2[r_i]);
	if(ref_seq_tmp2 != NULL)	free(ref_seq_tmp2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(mat_pos1[r_i] != NULL)	free(mat_pos1[r_i]);
	if(mat_pos1 != NULL)	free(mat_pos1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(mat_pos2[r_i] != NULL)	free(mat_pos2[r_i]);
	if(mat_pos2 != NULL)	free(mat_pos2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(mat_rc[r_i] != NULL)	free(mat_rc[r_i]);
	if(mat_rc != NULL)	free(mat_rc);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_no1[r_i] != NULL)	free(seed_no1[r_i]);
	if(seed_no1 != NULL)	free(seed_no1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_no2[r_i] != NULL)	free(seed_no2[r_i]);
	if(seed_no2 != NULL)	free(seed_no2);

#ifdef	PAIR_RANDOM

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < readlen_max / pair_ran_intvp; m++)
			if(seed_k_pos[r_i][m] != NULL)	free(seed_k_pos[r_i][m]);
		if(seed_k_pos[r_i] != NULL)	free(seed_k_pos[r_i]);
	}
	if(seed_k_pos != NULL)	free(seed_k_pos);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedposk[r_i] != NULL)	free(seedposk[r_i]);
	if(seedposk != NULL)	free(seedposk);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpos[0][0][r_i] != NULL)	free(seedpos[0][0][r_i]);
	if(seedpos[0][0] != NULL)	free(seedpos[0][0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpos[0][1][r_i] != NULL)	free(seedpos[0][1][r_i]);
	if(seedpos[0][1] != NULL)	free(seedpos[0][1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpos[1][0][r_i] != NULL)	free(seedpos[1][0][r_i]);
	if(seedpos[1][0] != NULL)	free(seedpos[1][0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpos[1][1][r_i] != NULL)	free(seedpos[1][1][r_i]);
	if(seedpos[1][1] != NULL)	free(seedpos[1][1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_posf[0][r_i] != NULL)	free(seed_posf[0][r_i]);
	if(seed_posf[0] != NULL)	free(seed_posf[0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_posf[1][r_i] != NULL)	free(seed_posf[1][r_i]);
	if(seed_posf[1] != NULL)	free(seed_posf[1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_single_pos[0][r_i] != NULL)	free(seed_single_pos[0][r_i]);
	if(seed_single_pos[0] != NULL)	free(seed_single_pos[0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_single_pos[1][r_i] != NULL)	free(seed_single_pos[1][r_i]);
	if(seed_single_pos[1] != NULL)	free(seed_single_pos[1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_single_ld[0][r_i] != NULL)	free(seed_single_ld[0][r_i]);
	if(seed_single_ld[0] != NULL)	free(seed_single_ld[0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_single_ld[1][r_i] != NULL)	free(seed_single_ld[1][r_i]);
	if(seed_single_ld[1] != NULL)	free(seed_single_ld[1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_single_rd[0][r_i] != NULL)	free(seed_single_rd[0][r_i]);
	if(seed_single_rd[0] != NULL)	free(seed_single_rd[0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_single_rd[1][r_i] != NULL)	free(seed_single_rd[1][r_i]);
	if(seed_single_rd[1] != NULL)	free(seed_single_rd[1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_single_dm[0][r_i] != NULL)	free(seed_single_dm[0][r_i]);
	if(seed_single_dm[0] != NULL)	free(seed_single_dm[0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_single_dm[1][r_i] != NULL)	free(seed_single_dm[1][r_i]);
	if(seed_single_dm[1] != NULL)	free(seed_single_dm[1]);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < (readlen_max / pair_ran_intvp) * RANDOM_RANGE; m++)
			if(seed_single_refs[0][r_i][m] != NULL)	free(seed_single_refs[0][r_i][m]);

		if(seed_single_refs[0][r_i] != NULL)	free(seed_single_refs[0][r_i]);
	}
	if(seed_single_refs[0] != NULL)	free(seed_single_refs[0]);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < (readlen_max / pair_ran_intvp) * RANDOM_RANGE; m++)
			if(seed_single_refs[1][r_i][m] != NULL)	free(seed_single_refs[1][r_i][m]);

		if(seed_single_refs[1][r_i] != NULL)	free(seed_single_refs[1][r_i]);
	}
	if(seed_single_refs[1] != NULL)	free(seed_single_refs[1]);

#ifdef	PR_SINGLE

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpos_mis[0][0][r_i] != NULL)	free(seedpos_mis[0][0][r_i]);
	if(seedpos_mis[0][0] != NULL)	free(seedpos_mis[0][0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpos_mis[0][1][r_i] != NULL)	free(seedpos_mis[0][1][r_i]);
	if(seedpos_mis[0][1] != NULL)	free(seedpos_mis[0][1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpos_mis[1][0][r_i] != NULL)	free(seedpos_mis[1][0][r_i]);
	if(seedpos_mis[1][0] != NULL)	free(seedpos_mis[1][0]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpos_mis[1][1][r_i] != NULL)	free(seedpos_mis[1][1][r_i]);
	if(seedpos_mis[1][1] != NULL)	free(seedpos_mis[1][1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pr_chr_res1[r_i] != NULL)	free(pr_chr_res1[r_i]);
	if(pr_chr_res1 != NULL)	free(pr_chr_res1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pr_sam_pos1[r_i] != NULL)	free(pr_sam_pos1[r_i]);
	if(pr_sam_pos1 != NULL)	free(pr_sam_pos1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pr_xa_d1[r_i] != NULL)	free(pr_xa_d1[r_i]);
	if(pr_xa_d1 != NULL)	free(pr_xa_d1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pr_lv_re1[r_i] != NULL)	free(pr_lv_re1[r_i]);
	if(pr_lv_re1 != NULL)	free(pr_lv_re1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pr_chr_res2[r_i] != NULL)	free(pr_chr_res2[r_i]);
	if(pr_chr_res2 != NULL)	free(pr_chr_res2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pr_sam_pos2[r_i] != NULL)	free(pr_sam_pos2[r_i]);
	if(pr_sam_pos2 != NULL)	free(pr_sam_pos2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pr_xa_d2[r_i] != NULL)	free(pr_xa_d2[r_i]);
	if(pr_xa_d2 != NULL)	free(pr_xa_d2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pr_lv_re2[r_i] != NULL)	free(pr_lv_re2[r_i]);
	if(pr_lv_re2 != NULL)	free(pr_lv_re2);

#endif

#ifdef	PAIR_RANDOM_SEED
	if(seed_r_dup != NULL)	free(seed_r_dup);
	if(random_buffer != NULL)	free(random_buffer);
#endif

#endif

#ifdef	READN_RANDOM_SEED
	if(random_buffer_readn)	free(random_buffer_readn);
#endif

	//thread gloabl value
	if(rc_cov_f != NULL)	free(rc_cov_f);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(cov_a_n[r_i] != NULL)	free(cov_a_n[r_i]);
	if(cov_a_n != NULL)	free(cov_a_n);

	if(cov_a_n_s != NULL)	free(cov_a_n_s);
	if(s_uid_f != NULL)	free(s_uid_f);
	if(min_mis != NULL)	free(min_mis);

	if(de_m_p != NULL)	free(de_m_p);
	if(de_m_p_o != NULL)	free(de_m_p_o);

	if(seed_l != NULL)	free(seed_l);

	if(max_mismatch1 != NULL)	free(max_mismatch1);
	if(max_mismatch2 != NULL)	free(max_mismatch2);

	if(max_mismatch1_single != NULL)	free(max_mismatch1_single);
	if(max_mismatch2_single != NULL)	free(max_mismatch2_single);

#ifdef	PR_SINGLE
	if(pr_o_f != NULL)	free(pr_o_f);
	if(seed_re_r != NULL)	free(seed_re_r);
#endif
	//if(off_start != NULL)	free(off_start);
	if(end_dis1 != NULL)	free(end_dis1);
	if(end_dis2 != NULL)	free(end_dis2);
	if(mat_posi != NULL)	free(mat_posi);

#ifdef	PAIR_SEED_LENGTH_FILT
	if(max_lengtht != NULL)	free(max_lengtht);
#endif

#ifdef	REDUCE_ANCHOR
	for(r_i = 0; r_i < thread_n; r_i++)
		if(orders1[r_i] != NULL)	free(orders1[r_i]);
	if(orders1 != NULL)	free(orders1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(orders2[r_i] != NULL)	free(orders2[r_i]);
	if(orders2 != NULL)	free(orders2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(dms1[r_i] != NULL)	free(dms1[r_i]);
	if(dms1 != NULL)	free(dms1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(dms2[r_i] != NULL)	free(dms2[r_i]);
	if(dms2 != NULL)	free(dms2);
#endif

	if(dm_op != NULL)	free(dm_op);
	if(dm_ops != NULL)	free(dm_ops);
	if(low_mask1 != NULL)	free(low_mask1);
	if(low_mask2 != NULL)	free(low_mask2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pos_add[r_i] != NULL)	free(pos_add[r_i]);
	if(pos_add != NULL)	free(pos_add);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(cigar_m1[r_i] != NULL)	free(cigar_m1[r_i]);
	if(cigar_m1 != NULL)	free(cigar_m1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(cigar_m2[r_i] != NULL)	free(cigar_m2[r_i]);
	if(cigar_m2 != NULL)	free(cigar_m2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ali_ref_seq[r_i] != NULL)	free(ali_ref_seq[r_i]);
	if(ali_ref_seq != NULL)	free(ali_ref_seq);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(read_char[r_i] != NULL)	free(read_char[r_i]);
	if(read_char != NULL)	free(read_char);

	if(pos_ren[0][0] != NULL)	free(pos_ren[0][0]);
	if(pos_ren[0][1] != NULL)	free(pos_ren[0][1]);
	if(pos_ren[1][0] != NULL)	free(pos_ren[1][0]);
	if(pos_ren[1][1] != NULL)	free(pos_ren[1][1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(sub_mask1[r_i] != NULL)	free(sub_mask1[r_i]);
	if(sub_mask1 != NULL)	free(sub_mask1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(sub_mask2[r_i] != NULL)	free(sub_mask2[r_i]);
	if(sub_mask2 != NULL)	free(sub_mask2);

	if(seedpos_misn[0][0] != NULL)	free(seedpos_misn[0][0]);
	if(seedpos_misn[0][1] != NULL)	free(seedpos_misn[0][1]);
	if(seedpos_misn[1][0] != NULL)	free(seedpos_misn[1][0]);
	if(seedpos_misn[1][1] != NULL)	free(seedpos_misn[1][1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ex_d_mask1[r_i] != NULL)	free(ex_d_mask1[r_i]);
	if(ex_d_mask1 != NULL)	free(ex_d_mask1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ex_d_mask2[r_i] != NULL)	free(ex_d_mask2[r_i]);
	if(ex_d_mask2 != NULL)	free(ex_d_mask2);

#ifdef	PAIR_RANDOM

	for(r_i = 0; r_i < thread_n; r_i++)
		if(b[r_i] != NULL)	free(b[r_i]);
	if(b != NULL)	free(b);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ls[r_i] != NULL)	free(ls[r_i]);
	if(ls != NULL)	free(ls);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(read_bit_pr[r_i] != NULL)	free(read_bit_pr[r_i]);
	if(read_bit_pr != NULL)	free(read_bit_pr);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(read_val_pr[r_i] != NULL)	free(read_val_pr[r_i]);
	if(read_val_pr != NULL)	free(read_val_pr);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(kcol[r_i] != NULL)	free(kcol[r_i]);
	if(kcol != NULL)	free(kcol);

	if(seed_posn[0] != NULL)	free(seed_posn[0]);
	if(seed_posn[1] != NULL)	free(seed_posn[1]);
	if(seedn != NULL)	free(seedn);

#endif
	if(read_bit_1 != NULL)	free(read_bit_1);
	if(read_bit_2 != NULL)	free(read_bit_2);
	if(read_val_1 != NULL)	free(read_val_1);
	if(read_val_2 != NULL)	free(read_val_2);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if(read_val1[r_i][0] != NULL)	free(read_val1[r_i][0]);
		if(read_val1[r_i][1] != NULL)	free(read_val1[r_i][1]);
		if(read_val1[r_i] != NULL)	free(read_val1[r_i]);
	}
	if(read_val1 != NULL)	free(read_val1);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if(read_val2[r_i][0] != NULL)	free(read_val2[r_i][0]);
		if(read_val2[r_i][1] != NULL)	free(read_val2[r_i][1]);
		if(read_val2[r_i] != NULL)	free(read_val2[r_i]);
	}
	if(read_val2 != NULL)	free(read_val2);


#ifdef	QUAL_FILT_LV
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if(qual_filt_lv1[r_i][0] != NULL)	free(qual_filt_lv1[r_i][0]);
		if(qual_filt_lv1[r_i][1] != NULL)	free(qual_filt_lv1[r_i][1]);
		if(qual_filt_lv1[r_i] != NULL)	free(qual_filt_lv1[r_i]);
	}
	if(qual_filt_lv1 != NULL)	free(qual_filt_lv1);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if(qual_filt_lv2[r_i][0] != NULL)	free(qual_filt_lv2[r_i][0]);
		if(qual_filt_lv2[r_i][1] != NULL)	free(qual_filt_lv2[r_i][1]);
		if(qual_filt_lv2[r_i] != NULL)	free(qual_filt_lv2[r_i]);
	}
	if(qual_filt_lv2 != NULL)	free(qual_filt_lv2);
#endif

#ifdef	PR_COV_FILTER
	if(cov_filt_f != NULL)	free(cov_filt_f);
#endif

#endif

	if(buffer_ref_seq != NULL)	free(buffer_ref_seq);
	if(buffer_seq != NULL)	free(buffer_seq);
	if(buffer_seqf != NULL)	free(buffer_seqf);
	//if(buffer_edge != NULL)	free(buffer_edge);
	if(buffer_p != NULL)	free(buffer_p);
	if(buffer_pp != NULL)	free(buffer_pp);
	if(buffer_hash_g != NULL)	free(buffer_hash_g);
	if(buffer_kmer_g != NULL)	free(buffer_kmer_g);
	if(buffer_off_g != NULL)	free(buffer_off_g);

#ifdef KSW_ALN_PAIR
	for(r_i = 0; r_i < thread_n; r_i++)
		if(ali_ref_seq2[r_i] != NULL)	free(ali_ref_seq2[r_i]);
	if(ali_ref_seq2 != NULL)	free(ali_ref_seq2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(read_char2[r_i] != NULL)	free(read_char2[r_i]);
	if(read_char2 != NULL)	free(read_char2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_kl1[r_i] != NULL)	free(op_dm_kl1[r_i]);
	if(op_dm_kl1 != NULL)	free(op_dm_kl1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_kr1[r_i] != NULL)	free(op_dm_kr1[r_i]);
	if(op_dm_kr1 != NULL)	free(op_dm_kr1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_kl2[r_i] != NULL)	free(op_dm_kl2[r_i]);
	if(op_dm_kl2 != NULL)	free(op_dm_kl2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_kr2[r_i] != NULL)	free(op_dm_kr2[r_i]);
	if(op_dm_kr2 != NULL)	free(op_dm_kr2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_kl1[r_i] != NULL)	free(ops_dm_kl1[r_i]);
	if(ops_dm_kl1 != NULL)	free(ops_dm_kl1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_kr1[r_i] != NULL)	free(ops_dm_kr1[r_i]);
	if(ops_dm_kr1 != NULL)	free(ops_dm_kr1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_kl2[r_i] != NULL)	free(ops_dm_kl2[r_i]);
	if(ops_dm_kl2 != NULL)	free(ops_dm_kl2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_kr2[r_i] != NULL)	free(ops_dm_kr2[r_i]);
	if(ops_dm_kr2 != NULL)	free(ops_dm_kr2);

	if(mat != NULL)	free(mat);
#endif

#ifdef	MAPPING_QUALITY

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if(mp_subs1[r_i][0] != NULL)	free(mp_subs1[r_i][0]);
		if(mp_subs1[r_i][1] != NULL)	free(mp_subs1[r_i][1]);
		if(mp_subs1[r_i] != NULL)	free(mp_subs1[r_i]);
	}
	if(mp_subs1 != NULL)	free(mp_subs1);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if(mp_subs2[r_i][0] != NULL)	free(mp_subs2[r_i][0]);
		if(mp_subs2[r_i][1] != NULL)	free(mp_subs2[r_i][1]);
		if(mp_subs2[r_i] != NULL)	free(mp_subs2[r_i]);
	}
	if(mp_subs2 != NULL)	free(mp_subs2);

	if(sub_t != NULL)	free(sub_t);

#endif

	return 0;
}
#ifdef	R_W_LOCK
int seed_ali_core(int read_seq_core, uint8_t tid)
#else
int seed_ali_core(uint32_t seqn, uint8_t tid)
#endif
{
	uint8_t rc_i = 0;
	uint8_t un_ii = 0;
	uint8_t lv_ref_length_re1 = 0;
	uint8_t lv_ref_length_re2 = 0;
	uint8_t q_n1 = 0;
	uint8_t q_n2 = 0;
	uint8_t c_m_f = 0;
	uint8_t end1_uc_f = 0;
	uint8_t end2_uc_f = 0;
	uint8_t nuc1_f = 0;
	uint8_t nuc2_f = 0;
	uint8_t b_t_n_r = 0;
	uint8_t re_d = 0;
	uint8_t single_lv_re = 0;
	uint8_t rc_ii = 0;
	uint8_t unmatch[2][2];
	uint16_t max_mismatch = 0;
	//uint8_t flow_f = 0;

	char cigar_p1[MAX_LV_CIGARCOM];
	char cigar_p2[MAX_LV_CIGARCOM];

	uint16_t f_cigar[MAX_LV_CIGAR];
	char sam_seq1[MAX_READLEN] = {};
	char sam_seq2[MAX_READLEN] = {};
	char cigarBuf1[MAX_LV_CIGAR] = {};
	char cigarBuf2[MAX_LV_CIGAR] = {};
	char str_o[MAX_LV_CIGAR];
	char b_cigar[MAX_LV_CIGAR];
	char* cp1 = NULL;
	char* cp2 = NULL;
	char* seq1p = NULL;
	char* seq2p = NULL;
	char* pch = NULL;
	char* saveptr = NULL;
	char* cigar_m_1 = NULL;
	char* cigar_m_2 = NULL;

	uint16_t ref_copy_num = 0;
	uint16_t ref_copy_num1 = 0;
	uint16_t ref_copy_num2 = 0;
	uint16_t ref_copy_num_1 = 0;
	uint16_t ref_copy_num_2 = 0;
	uint16_t f_cigarn = 0;
	uint16_t read_bit_char1 = 0;
	uint16_t read_bit_char2 = 0;
	int16_t rst_i = 0;
	uint16_t s_m_t = 0;
	uint16_t read_b_i = 0;
	uint16_t f_c = 0;
	uint16_t pchl = 0;
	uint16_t f_i = 0;
	uint16_t s_o = 0;
	uint16_t snt = 0;
	uint16_t d_n1 = 0;
	uint16_t i_n1 = 0;
	uint16_t no1_tmp = 0;
	uint16_t no2_tmp = 0;
	uint16_t read_length = 0;
	uint16_t read_length1 = 0;
	uint16_t read_length2 = 0;
	uint16_t read_length_1 = 0;
	uint16_t read_length_2 = 0;
	uint16_t dm_i = 0;

	uint32_t read_length_a1 = 0;
	uint32_t read_length_a2 = 0;
	uint32_t seqi = 0;
	uint32_t v_cnt = 0;
	uint32_t vs_cnt = 0;
	uint32_t ref_copy_num_chars = 0;
	uint32_t ref_copy_num_chars1 = 0;
	uint32_t ref_copy_num_chars2 = 0;
	uint32_t ref_copy_num_chars_1 = 0;
	uint32_t ref_copy_num_chars_2 = 0;
	uint32_t r_i = 0;
	uint32_t psp_i = 0;
	uint32_t d_l1 = 0;
	uint32_t d_l2 = 0;
	uint32_t d_r1 = 0;
	uint32_t d_r2 = 0;
	uint32_t mis_c_n = 0;
	uint32_t mis_c_n_filt = 0;
	uint32_t xa_i = 0;
	uint32_t xa_i_1 = 0;
	uint32_t xa_i_2 = 0;
	uint32_t v_cnt_i = 0;
	uint32_t va_cnt_i = 0;
	uint32_t sam_flag1 = 0;
	uint32_t sam_flag2 = 0;
	uint32_t seed1_i = 0;
	uint32_t seed2_i = 0;
	uint32_t sam_seq_i = 0;
	uint32_t rc_f = 0;
	uint32_t max_pair_score = 0;
	uint32_t seed_posn_filter = 0;

#ifdef QUALS_CHECK
	uint8_t error_f = 0;
#endif

	int	dm_op_p = 0;

	int s_r_o_l1 = 0;
	int	s_r_o_r1 = 0;
	int s_r_o_l2 = 0;
	int	s_r_o_r2 = 0;
	int lv_re1f = 0;
	int lv_re1b = 0;
	int lv_re2f = 0;
	int lv_re2b = 0;
	int chr_re = 0;
	int chr_re1 = 0;
	int chr_re2 = 0;
	int mid = 0;
	int low = 0;
	int high = 0;
	int m_m_n = 0;
	int sn = 0;
	int dm12 = 0;
	int dm1 = 0;
	int dm2 = 0;
	int dm_l1 = 0;
	int dm_r1 = 0;
	int dm_l2 = 0;
	int dm_r2 = 0;
	int dmt1 = 0;
	int dmt2 = 0;
	int lv_dmt1 = 0;
	int lv_dmt2 = 0;
	int ld1 = 0;
	int rd1 = 0;
	int ld2 = 0;
	int rd2 = 0;
	int ksw_re = 0;
	int cmp_re = 0;
	int q_rear_i = 0;
	int q_rear1 = 0;
	int q_rear2 = 0;
	int cache_dml1[MAX_Q_NUM];
	int cache_dml2[MAX_Q_NUM];
	int cache_dmr1[MAX_Q_NUM];
	int cache_dmr2[MAX_Q_NUM];
	int cache_dis1[MAX_Q_NUM];
	int cache_dis2[MAX_Q_NUM];
	int cache_kl1[MAX_Q_NUM];
	int cache_kr1[MAX_Q_NUM];
	int cache_kl2[MAX_Q_NUM];
	int cache_kr2[MAX_Q_NUM];
	int nm_score1 = 0;
	int nm_score2 = 0;

	uint32_t spa_i[2][2];
	uint64_t c_tmp = 0;
	uint64_t xor_tmp = 0;
	uint64_t ref_tmp_ori = 0;
	uint64_t ref_tmp_ori2 = 0;
	uint64_t tran_tmp_p = 0;
	uint64_t tran_tmp = 0;
	int64_t sam_cross = 0;
	int64_t sam_pos1 = 0;
	int64_t sam_pos2  = 0;
	uint64_t cache_end1[MAX_Q_NUM][MAX_REF_SEQ_C];
	uint64_t cache_end2[MAX_Q_NUM][MAX_REF_SEQ_C];
	uint64_t low_mask = 0;
	uint16_t* sub_mask = NULL;
	uint64_t* ex_d_mask = NULL;

#ifdef	SEED_FILTER_LENGTH
	uint16_t seed_no[2][2];
	uint16_t max_read_length[2][2];
	uint16_t max_sets_n[2][2];
	uint16_t off_i = 0;
#endif

	int64_t bit_char_i = 0;
	int64_t bit_char_i_ref = 0;
	cnt_re cnt;

#ifdef	ANCHOR_HASH_ALI
	uint8_t ori_i = 0;
	uint16_t tra_i;
	uint8_t base_re = 0;
	uint8_t anchor_hash_i = 0;

	uint16_t des_i = 0;
	uint16_t print_i;
	uint16_t array_index;
	uint16_t char_i = 0;
	uint16_t buffer_i = 0;
	uint16_t seed_length = 0;
	uint16_t max_seed_length = 0;
	uint16_t lv_k = 0;
	uint16_t lv_k1 = 0;
	uint16_t lv_k2 = 0;
	uint16_t lv_k_1 = 0;
	uint16_t lv_k_2 = 0;

	uint32_t read_k_p = 0;
	uint32_t read_k_t = 0;
	uint32_t anchor_mask = 0;
	uint32_t array_i = 0;
	uint32_t array_i_p = 0;
	uint32_t anchor_ref = 0;
	uint32_t r_b_v = 0;
	uint32_t max_right = 0;

	int left_i = 0;
	int right_i = 0;

	uint64_t* ori = NULL;
	anchor_mask = ((uint32_t )1 << k_anchor_b) - 1;

	uint16_t* anchor_hash_p = NULL;
	uint8_t* anchor_array_p = NULL;
	uint16_t* anchor_point_p = NULL;
	uint16_t* anchor_pos_p = NULL;

#endif

	seed_pa* seed_pr1 = NULL;
	seed_pa* seed_pr2 = NULL;

	seed_pa* seed_pr[2][2];
	seed_pr[0][0] = NULL;
	seed_pr[0][1] = NULL;
	seed_pr[1][0] = NULL;
	seed_pr[1][1] = NULL;

	//for # position
	int16_t pound_pos1_f_forward = 0;
	int16_t pound_pos1_f_reverse = 0;
	int16_t pound_pos1_r_forward = 0;
	int16_t pound_pos1_r_reverse = 0;
	int16_t pound_pos2_f_forward = 0;
	int16_t pound_pos2_f_reverse = 0;
	int16_t pound_pos2_r_forward = 0;
	int16_t pound_pos2_r_reverse = 0;

	int16_t pound_mis = 0;
	int16_t pound_pos_1_f = 0;
	int16_t pound_pos_1_r = 0;
	int16_t pound_pos_2_f = 0;
	int16_t pound_pos_2_r = 0;
	int16_t lv_up_left = 0;
	int16_t lv_up_right = 0;
	int16_t lv_down_right = 0;
	int16_t lv_down_left = 0;
	int16_t dm_cir_1 = 0;
	int16_t dm_cir_2 = 0;
	int16_t dm_cir = 0;
	int16_t dm_cir_min = 0;
	int16_t cache_dm_cir1[MAX_Q_NUM];
	int16_t cache_dm_cir2[MAX_Q_NUM];

	int16_t read_pos_re = 0;
	int16_t read_pos_end_re = 0;
	int16_t read_pos_start_num = 0;
	int16_t read_pos_end_num = 0;
	int16_t rst_i_1 = 0;
	int16_t m_n_f = 0;
	int16_t m_n_b = 0;
	int16_t s_r_o_l = 0;
	int16_t s_r_o_r = 0;
	uint8_t cir_n = 0;

#ifdef	KSW_ALN_PAIR
	int band_with = 33;//100
	int zdrop = 0;//100
	int end_bonus = 5;
	int tle;
	int gtle;
	int gscore;
	int max_off;
#endif

#ifdef	QUAL_FILT
	uint64_t* qual_filt_1 = NULL;
	uint64_t* qual_filt_2 = NULL;
	uint8_t qual_filt_fix = 55;
#endif

#ifdef	QUAL_FILT_LV
	uint8_t* qual_filt_lv_1 = NULL;
	uint8_t* qual_filt_lv_2 = NULL;
	uint8_t* qual_filt_lv_1_o = NULL;
	uint8_t* qual_filt_lv_2_o = NULL;
	uint16_t mis_n1 = 0;
	uint16_t mis_n2 = 0;
#endif

	uint16_t s_offset1 = 0;
	uint16_t s_offset2 = 0;

#ifdef	ALTER_DEBUG
	uint16_t seed_length1 = 0;
	uint16_t seed_length2 = 0;
#endif

#ifdef	MAPPING_QUALITY
	float m_sub_tmp = 0;
#endif

#ifdef UNPIPATH_OFF_K20
	uint64_t pos_l = 0;
	uint64_t posi = 0;
	uint64_t x = 0;
	uint64_t ksw_s = 0;
	uint64_t ksw_e = 0;
	uint64_t base_i = 0;
	uint64_t base_i_off_l = 0;
	uint64_t base_i_off_r = 0;
#else
	uint32_t pos_l = 0;
	uint32_t posi = 0;
	uint32_t x = 0;
	uint32_t ksw_s = 0;
	uint32_t ksw_e = 0;
	uint32_t base_i = 0;
	uint32_t base_i_off_l = 0;
	uint32_t base_i_off_r = 0;
#endif

	uint64_t low_mask_1 = 0;
	uint64_t low_mask_2 = 0;

#ifdef	READN_RANDOM_SEED
	char tmp_char;
	uint16_t readn_re1 = 0;
	uint16_t readn_re2 = 0;
#endif

#ifdef	CIGAR_LEN_ERR
	int cigar_len = 0;
	int cigar_len_tmp = 0;
	int cigar_len_re = 0;
	uint16_t pchl_tmp = 0;
	uint16_t s_o_tmp = 0;
	char* pch_tmp = NULL;
	char* saveptr_tmp = NULL;
	char cigar_tmp[MAX_LV_CIGARCOM];
#endif

	uint8_t qual_flag = 0;

#ifdef	REDUCE_ANCHOR
	int ls = 0;
	int rs = 0;
	uint8_t rcs = 0;
	uint16_t tra1_i = 0;
	uint16_t tra2_i = 0;
	uint16_t anchor_n1 = 0;
	uint16_t anchor_n2 = 0;
	uint8_t other_end_flag = 0;
	uint16_t op_mask[MAX_REDUCE_ANCHOR_NUM];
	uint16_t ops_mask[MAX_REDUCE_ANCHOR_NUM];
#endif

#ifdef	ANCHOR_LV_S
	int16_t s_offset_l = 0;
	int16_t s_offset_r = 0;
#endif

	uint64_t* op_vector_seq1_tmp = NULL;
	uint16_t v_cnt_i_tmp = 0;
#ifdef	FIX_SV
	uint16_t tra_i_n = 0;
	uint8_t sv_add = 0;
#endif

#ifdef FIX_SA
	uint8_t op_rc_tmp = 0;
	uint16_t sv_s_len = 0;
	uint16_t sv_s_len_p = 0;
#endif

#ifndef	R_W_LOCK
	for(seqi = 0; seqi < seqn; seqi++)
#else
	if(1)
#endif
	{

#ifndef	R_W_LOCK

#ifdef PTHREAD_USE
		if ((seqi % thread_n) != tid) continue;
#endif

#else
		seqi = read_seq_core;
#endif

		read_length1 = seqio[seqi].read_length1;
		read_length2 = seqio[seqi].read_length2;

		lv_k1 = (read_length1 * lv_rate) + 1;
		lv_k2 = (read_length2 * lv_rate) + 1;

		if((insert_dis < read_length1) || (insert_dis < read_length2))
		{
			fprintf(stderr, "Input error: wrong -u or -f\n");
			exit(1);
		}

		end_dis1[tid] = insert_dis - read_length2;
		end_dis2[tid] = insert_dis - read_length1;

		read_length_a1 = read_length1 - 1;
		read_length_a2 = read_length2 - 1;

		f_cigarn = (read_length1 > read_length2) ? read_length1: read_length2;

		f_cigar[f_cigarn] = '\0';

		//sprintf(cigar_m1[tid],"%uM\0",read_length1);
		//sprintf(cigar_m2[tid],"%uM\0",read_length2);
		sprintf(cigar_m1[tid],"%uM",read_length1);
		sprintf(cigar_m2[tid],"%uM",read_length2);

		//for bit operation
		read_bit_char1 = (((uint16_t )((read_length_a1 >> 5) + 1)) << 3);
		read_bit_char2 = (((uint16_t )((read_length_a2 >> 5) + 1)) << 3);

		ref_copy_num1 = ((read_length_a1) >> 5) + 3;
		ref_copy_num_chars1 = (ref_copy_num1 << 3);
		lv_ref_length_re1 = (read_length1 & 0X1f);

		for(r_i = 0; r_i <= 32 - lv_ref_length_re1; r_i++)
		{
			ex_d_mask1[tid][r_i] = bit_assi[lv_ref_length_re1 + r_i];
			sub_mask1[tid][r_i] = ref_copy_num1 - 2;
		}

		for(r_i = 33 - lv_ref_length_re1; r_i <= 32; r_i++)
		{
			ex_d_mask1[tid][r_i] = bit_assi[lv_ref_length_re1 + r_i - 32];
			sub_mask1[tid][r_i] = ref_copy_num1 - 1;
		}

		low_mask1[tid] = bit_assi[lv_ref_length_re1];

		//for bit operation
		ref_copy_num2 = ((read_length_a2) >> 5) + 3;
		ref_copy_num_chars2 = (ref_copy_num2 << 3);
		lv_ref_length_re2 = (read_length2 & 0X1f);

		for(r_i = 0; r_i <= 32 - lv_ref_length_re2; r_i++)
		{
			ex_d_mask2[tid][r_i] = bit_assi[lv_ref_length_re2 + r_i];
			sub_mask2[tid][r_i] = ref_copy_num2 - 2;
		}

		for(r_i = 33 - lv_ref_length_re2; r_i <= 32; r_i++)
		{
			ex_d_mask2[tid][r_i] = bit_assi[lv_ref_length_re2 + r_i - 32];
			sub_mask2[tid][r_i] = ref_copy_num2 - 1;
		}

		low_mask2[tid] = bit_assi[lv_ref_length_re2];

		max_pair_score = (uint32_t )(((float )(read_length1 + read_length2)) * max_pair_score_r);
		cov_a_n[tid][0] = ((read_length_a1) >> 6) + 1;//2
		cov_a_n[tid][1] = ((read_length_a2) >> 6) + 1;
		max_mismatch1[tid] = (uint16_t )(((float )read_length1) * mis_match_r);
		max_mismatch2[tid] = (uint16_t )(((float )read_length2) * mis_match_r);

		max_mismatch1_single[tid] = (uint16_t )(((float )read_length1) * mis_match_r_single);
		max_mismatch2_single[tid] = (uint16_t )(((float )read_length2) * mis_match_r_single);

		memset(read_bit1[tid][0], 0, read_bit_char1);
		memset(read_bit1[tid][1], 0, read_bit_char1);
		memset(read_bit2[tid][0], 0, read_bit_char2);
		memset(read_bit2[tid][1], 0, read_bit_char2);

		r_i = 0;
		while ((seqio[seqi].read_seq1)[r_i])
		{
#ifdef	READN_RANDOM_SEED
			tmp_char = (seqio[seqi].read_seq1)[r_i];
			if(tmp_char == 'N')
			{
				readn_re1 = readn_cnt & 0X3ff;
				readn_re2 = readn_cnt & 0Xf;
				c_tmp = ((random_buffer_readn[readn_re1] >> (readn_re2 << 1)) & 0X3);
				readn_cnt++;
			}
			else	c_tmp = charToDna5n[tmp_char];
#else
			c_tmp = charToDna5n[(seqio[seqi].read_seq1)[r_i]];
#endif

#ifdef ALI_LV
			read_val1[tid][0][r_i] = c_tmp;
			read_val1[tid][1][read_length_a1 - r_i] = c_tmp ^ 0X3;
#endif
			read_bit1[tid][0][r_i >> 5] |= (((uint64_t )c_tmp) << ((31 - (r_i & 0X1f)) << 1));
			read_bit1[tid][1][(read_length_a1 - r_i) >> 5] |= (((uint64_t )(c_tmp ^ 0X3)) << ((31 - ((read_length_a1 - r_i) & 0X1f)) << 1));

			r_i++;
		}

		r_i = 0;
		while ((seqio[seqi].read_seq2)[r_i])
		{
#ifdef	READN_RANDOM_SEED
			tmp_char = (seqio[seqi].read_seq2)[r_i];
			if(tmp_char == 'N')
			{
				readn_re1 = readn_cnt & 0X3ff;
				readn_re2 = readn_cnt & 0Xf;
				c_tmp = ((random_buffer_readn[readn_re1] >> (readn_re2 << 1)) & 0X3);
				readn_cnt++;
			}
			else	c_tmp = charToDna5n[tmp_char];
#else
			c_tmp = charToDna5n[(seqio[seqi].read_seq2)[r_i]];
#endif

#ifdef ALI_LV
			read_val2[tid][0][r_i] = c_tmp;
			read_val2[tid][1][read_length_a2 - r_i] = c_tmp ^ 0X3;
#endif
			read_bit2[tid][0][r_i >> 5] |= (((uint64_t )c_tmp) << ((31 - (r_i & 0X1f)) << 1));
			read_bit2[tid][1][(read_length_a2 - r_i) >> 5] |= (((uint64_t )(c_tmp ^ 0X3)) << ((31 - ((read_length_a2 - r_i) & 0X1f)) << 1));

			r_i++;
		}

#ifdef	QUAL_FILT

		memset(qual_filt1[tid][0], 0, read_bit_char1);
		memset(qual_filt1[tid][1], 0, read_bit_char1);
		memset(qual_filt2[tid][0], 0, read_bit_char2);
		memset(qual_filt2[tid][1], 0, read_bit_char2);

		c_tmp = 3;
		for(r_i = 0; r_i < read_length1; r_i++)
		{
			if(seqio[seqi].qual1[r_i] < qual_filt_fix)   //'7': 55 63
			{
				qual_filt1[tid][0][r_i >> 5] |= (((uint64_t )c_tmp) << ((31 - (r_i & 0X1f)) << 1));
				qual_filt1[tid][1][(read_length_a1 - r_i) >> 5] |= (((uint64_t )c_tmp) << ((31 - ((read_length_a1 - r_i) & 0X1f)) << 1));
			}
		}

		for(r_i = 0; r_i < read_length2; r_i++)
		{
			if(seqio[seqi].qual2[r_i] < qual_filt_fix)   //'7': 55 63
			{
				qual_filt2[tid][0][r_i >> 5] |= (((uint64_t )c_tmp) << ((31 - (r_i & 0X1f)) << 1));
				qual_filt2[tid][1][(read_length_a2 - r_i) >> 5] |= (((uint64_t )c_tmp) << ((31 - ((read_length_a2 - r_i) & 0X1f)) << 1));
			}
		}

		for(r_i = 0; r_i < (read_length1 >> 5) + 1; r_i++)
			qual_filt1[tid][0][r_i] = ~qual_filt1[tid][0][r_i];
		for(r_i = 0; r_i < (read_length1 >> 5) + 1; r_i++)
			qual_filt1[tid][1][r_i] = ~qual_filt1[tid][1][r_i];
		for(r_i = 0; r_i < (read_length2 >> 5) + 1; r_i++)
			qual_filt2[tid][0][r_i] = ~qual_filt2[tid][0][r_i];
		for(r_i = 0; r_i < (read_length2 >> 5) + 1; r_i++)
			qual_filt2[tid][1][r_i] = ~qual_filt2[tid][1][r_i];

#ifdef	QUAL_FILT_LV
		for(r_i = 0; r_i < read_length1; r_i++)
		{
			if(seqio[seqi].qual1[r_i] < qual_filt_fix)	qual_filt_lv1[tid][0][r_i] = 0;
			else	qual_filt_lv1[tid][0][r_i] = 3;
		}
		for(r_i = 0; r_i < read_length1; r_i++)
		{
			if(seqio[seqi].qual1[r_i] < qual_filt_fix)	qual_filt_lv1[tid][1][read_length_a1 - r_i] = 0;
			else	qual_filt_lv1[tid][1][read_length_a1 - r_i] = 3;
		}

		for(r_i = 0; r_i < read_length2; r_i++)
		{
			if(seqio[seqi].qual2[r_i] < qual_filt_fix)	qual_filt_lv2[tid][0][r_i] = 0;
			else	qual_filt_lv2[tid][0][r_i] = 3;
		}
		for(r_i = 0; r_i < read_length2; r_i++)
		{
			if(seqio[seqi].qual2[r_i] < qual_filt_fix)	qual_filt_lv2[tid][1][read_length_a2 - r_i] = 0;
			else	qual_filt_lv2[tid][1][read_length_a2 - r_i] = 3;
		}
#endif

#endif



#ifdef	MAPPING_QUALITY
		for(r_i = 0; r_i < read_length1; r_i++)
		{
			m_sub_tmp = mp_sub_bp[seqio[seqi].qual1[r_i]];
			mp_subs1[tid][0][r_i] = m_sub_tmp;
			mp_subs1[tid][1][read_length_a1 - r_i] = m_sub_tmp;
		}

		for(r_i = 0; r_i < read_length2; r_i++)
		{
			m_sub_tmp = mp_sub_bp[seqio[seqi].qual2[r_i]];
			mp_subs2[tid][0][r_i] = m_sub_tmp;
			mp_subs2[tid][1][read_length_a2 - r_i] = m_sub_tmp;
		}

#endif
		pound_pos1_f_forward = read_length1;
		pound_pos1_r_forward = read_length1;
		pound_pos1_f_reverse = 0;
		pound_pos1_r_reverse = 0;

		pound_pos2_f_forward = read_length2;
		pound_pos2_r_forward = read_length2;
		pound_pos2_f_reverse = 0;
		pound_pos2_r_reverse = 0;

		for(r_i = 0; r_i < read_length1; r_i++)
			if(seqio[seqi].qual1[r_i] == '#')
			{
				pound_pos1_f_forward = r_i;
				pound_pos1_r_reverse = read_length1 - r_i;
				break;
			}

		for(r_i = 0; r_i < read_length2; r_i++)
			if(seqio[seqi].qual2[r_i] == '#')
			{
				pound_pos2_f_forward = r_i;
				pound_pos2_r_reverse = read_length2 - r_i;
				break;
			}

		de_m_p[tid] = 1;
		de_m_p_o[tid] = 1;
		seed_l[tid] = seed_l_max;

#ifdef	PAIR_SEED_LENGTH_FILT
		max_lengtht[tid] = 0;
#endif
		rc_f = 0;

#ifdef	PR_COV_FILTER
		cov_filt_f[tid] = 0;
#endif

#ifdef	PR_SINGLE
		seed_re_r[tid] = 0;
#endif

#ifdef	TER_J
		dm_op_p = 0;
		v_cnt_p = 0;
#endif
		dm_op_p = MAX_OP_SCORE;

#ifdef	ALTER_DEBUG
		rep_go[tid] = 1;
#endif

		cir_n = 1;
		while(cir_n <= cir_fix_n)
		{
			dm_op[tid] = MAX_OP_SCORE;
			dm_ops[tid] = MAX_OP_SCORE;
			v_cnt = 0;
			vs_cnt = 0;
			mat_posi[tid] = 0;
			rc_f = 0Xffffffff;
			dm_cir_min = 0Xfff;

			unmatch[0][0] = 0;
			unmatch[0][1] = 0;
			unmatch[1][0] = 0;
			unmatch[1][1] = 0;

#ifdef	PAIR_RANDOM
			if(cir_n == 1)
			{
				pos_ren[0][0][tid] = 0;
				pos_ren[0][1][tid] = 0;
				pos_ren[1][0][tid] = 0;
				pos_ren[1][1][tid] = 0;
			}
#endif

			for(rc_i = 0; rc_i < 2; rc_i++)
			{
				rc_cov_f[tid] = 0;
#ifdef	SEED_FILTER_LENGTH
				max_read_length[rc_i][rc_i] = 0;

#ifdef UNPIPATH_OFF_K20
				if((seed_pr[rc_i][rc_i] = single_seed_reduction_core_filter64(seedpa1[rc_i][tid], read_bit1[tid][rc_i], read_val1[tid][rc_i], &(spa_i[rc_i][rc_i]), &(seed_set_pos[rc_i][rc_i][tid]), &(pos_ren[rc_i][rc_i][tid]), tid, read_length1, &(seed_no[rc_i][rc_i]), &(max_read_length[rc_i][rc_i]))) == NULL)
#else
				if((seed_pr[rc_i][rc_i] = single_seed_reduction_core_filter(seedpa1[rc_i][tid], read_bit1[tid][rc_i], read_val1[tid][rc_i], &(spa_i[rc_i][rc_i]), &(seed_set_pos[rc_i][rc_i][tid]), &(pos_ren[rc_i][rc_i][tid]), tid, read_length1, &(seed_no[rc_i][rc_i]), &(max_read_length[rc_i][rc_i]))) == NULL)
#endif

#else
				if((seed_pr[rc_i][rc_i] = single_seed_reduction_core(seedpa1[rc_i][tid], read_bit1[tid][rc_i], read_val1[tid][rc_i], &(spa_i[rc_i][rc_i]), &(seed_set_pos[rc_i][rc_i][tid]), &(pos_ren[rc_i][rc_i][tid]), tid, read_length1)) == NULL)
#endif
				{
#ifdef	UNMATCH_SINGLE_END
					unmatch[rc_i][rc_i] = 1;
#endif

#ifdef	UNMATCH_SINGLE_END_MIS
					unmatch[rc_i][rc_i] = 1;
#endif
				}

				rc_cov_f[tid] = 1;
#ifdef	SEED_FILTER_LENGTH
				max_read_length[rc_i][1 - rc_i] = 0;
#ifdef UNPIPATH_OFF_K20
				if((seed_pr[rc_i][1 - rc_i] = single_seed_reduction_core_filter64(seedpa2[rc_i][tid], read_bit2[tid][1 - rc_i], read_val2[tid][1 - rc_i], &(spa_i[rc_i][1 - rc_i]), &(seed_set_pos[rc_i][1 - rc_i][tid]), &(pos_ren[rc_i][1 - rc_i][tid]), tid, read_length2, &(seed_no[rc_i][1 - rc_i]), &(max_read_length[rc_i][1 - rc_i]))) == NULL)
#else
				if((seed_pr[rc_i][1 - rc_i] = single_seed_reduction_core_filter(seedpa2[rc_i][tid], read_bit2[tid][1 - rc_i], read_val2[tid][1 - rc_i], &(spa_i[rc_i][1 - rc_i]), &(seed_set_pos[rc_i][1 - rc_i][tid]), &(pos_ren[rc_i][1 - rc_i][tid]), tid, read_length2, &(seed_no[rc_i][1 - rc_i]), &(max_read_length[rc_i][1 - rc_i]))) == NULL)
#endif
#else
				if((seed_pr[rc_i][1 - rc_i] = single_seed_reduction_core(seedpa2[rc_i][tid], read_bit2[tid][1 - rc_i], read_val2[tid][1 - rc_i], &(spa_i[rc_i][1 - rc_i]), &(seed_set_pos[rc_i][1 - rc_i][tid]), &(pos_ren[rc_i][1 - rc_i][tid]), tid, read_length2)) == NULL)
#endif
				{
#ifdef	UNMATCH_SINGLE_END
					unmatch[rc_i][1 - rc_i] = 1;
#endif

#ifdef	UNMATCH_SINGLE_END_MIS
					unmatch[rc_i][1 - rc_i] = 1;
#endif
				}

				if((seed_pr[rc_i][rc_i] != NULL) && (seed_pr[rc_i][1 - rc_i] != NULL))
				{
#ifdef UNPIPATH_OFF_K20
					merge_pair_end64(spa_i[rc_i][0], spa_i[rc_i][1], seed_pr[rc_i][0], seed_pr[rc_i][1], seed_set_pos[rc_i][0][tid], seed_set_pos[rc_i][1][tid], rc_i, tid);
#else
					merge_pair_end(spa_i[rc_i][0], spa_i[rc_i][1], seed_pr[rc_i][0], seed_pr[rc_i][1], seed_set_pos[rc_i][0][tid], seed_set_pos[rc_i][1][tid], rc_i, tid);
#endif
					if(rc_i == 0)	rc_f = mat_posi[tid];
				}
			}

			if(de_m_p[tid] == 1)
			{
				seed_l[tid] -= seed_step;
			}
			else
			{
				no1_tmp = 0Xffff;
				no2_tmp = 0Xffff;

				for(r_i = 0; r_i < mat_posi[tid]; r_i++)
				{
					seed1_i = seed_no1[tid][r_i];
					seed2_i = seed_no2[tid][r_i];

					if(mat_rc[tid][r_i] == 0)
					{
						read_bit_1[tid] = read_bit1[tid][0];
						read_bit_2[tid] = read_bit2[tid][1];
						read_val_1[tid] = read_val1[tid][0];
						read_val_2[tid] = read_val2[tid][1];
#ifdef	QUAL_FILT
						qual_filt_1 = qual_filt1[tid][0];
						qual_filt_2 = qual_filt2[tid][1];
#endif

#ifdef	QUAL_FILT_LV
						qual_filt_lv_1 = qual_filt_lv1[tid][0];
						qual_filt_lv_2 = qual_filt_lv2[tid][1];
						qual_filt_lv_1_o = qual_filt_lv1[tid][1];
						qual_filt_lv_2_o = qual_filt_lv2[tid][0];
#endif
						seed_pr1 = seed_pr[0][0];
						seed_pr2 = seed_pr[0][1];

						read_length_1 = read_length1;
						read_length_2 = read_length2;

						ref_copy_num_1 = ref_copy_num1;
						ref_copy_num_2 = ref_copy_num2;
						ref_copy_num_chars_1 = ref_copy_num_chars1;
						ref_copy_num_chars_2 = ref_copy_num_chars2;

						low_mask_1 = low_mask1[tid];
						low_mask_2 = low_mask2[tid];

						pound_pos_1_f = pound_pos1_f_forward;
						pound_pos_1_r = pound_pos1_r_forward;
						pound_pos_2_f = pound_pos2_f_reverse;
						pound_pos_2_r = pound_pos2_r_reverse;
					}
					if(mat_rc[tid][r_i] == 1)
					{
						read_bit_1[tid] = read_bit2[tid][0];
						read_bit_2[tid] = read_bit1[tid][1];
						read_val_1[tid] = read_val2[tid][0];
						read_val_2[tid] = read_val1[tid][1];
#ifdef	QUAL_FILT
						qual_filt_1 = qual_filt2[tid][0];
						qual_filt_2 = qual_filt1[tid][1];
#endif

#ifdef	QUAL_FILT_LV
						qual_filt_lv_1 = qual_filt_lv2[tid][0];
						qual_filt_lv_2 = qual_filt_lv1[tid][1];
						qual_filt_lv_1_o = qual_filt_lv2[tid][1];
						qual_filt_lv_2_o = qual_filt_lv1[tid][0];
#endif
						seed_pr1 = seed_pr[1][0];
						seed_pr2 = seed_pr[1][1];

						read_length_1 = read_length2;
						read_length_2 = read_length1;

						ref_copy_num_1 = ref_copy_num2;
						ref_copy_num_2 = ref_copy_num1;
						ref_copy_num_chars_1 = ref_copy_num_chars2;
						ref_copy_num_chars_2 = ref_copy_num_chars1;

						low_mask_1 = low_mask2[tid];
						low_mask_2 = low_mask1[tid];

						pound_pos_1_f = pound_pos2_f_forward;
						pound_pos_1_r = pound_pos2_r_forward;
						pound_pos_2_f = pound_pos1_f_reverse;
						pound_pos_2_r = pound_pos1_r_reverse;
					}

#ifdef	PAIR_SEED_LENGTH_FILT
					if(seed_pr1[seed1_i].length + seed_pr2[seed2_i].length < max_lengtht[tid] - length_reducet)	continue;
#endif

					if(r_i == rc_f)
					{
						no1_tmp = 0Xffff;
						no2_tmp = 0Xffff;
					}

					if(seed1_i != no1_tmp)
					{
						//current seed does not cross one unipath
						if(seed_pr1[seed1_i].ui == 1)
						{
							end1_uc_f = 0;
							nuc1_f = 0;

							d_l1 = seed_pr1[seed1_i].ref_pos_off;
							d_r1 = seed_pr1[seed1_i].ref_pos_off_r;
						}
						else
						{
							end1_uc_f = 1;
						}
						q_rear1 = 0;
						q_n1 = 0;

						dmt1 = ali_exl;
						lv_dmt1 = lv_k1;

						s_r_o_l1 = seed_pr1[seed1_i].s_r_o_l;
						s_r_o_r1 = seed_pr1[seed1_i].s_r_o_r;
#ifdef	ALTER_DEBUG
						seed_length1 = seed_pr1[seed1_i].length;
#endif
						no1_tmp = seed1_i;
					}
					if(seed2_i != no2_tmp)
					{
						//current seed does not cross one unipath
						if(seed_pr2[seed2_i].ui == 1)
						{
							end2_uc_f = 0;
							nuc2_f = 0;

							d_l2 = seed_pr2[seed2_i].ref_pos_off;
							d_r2 = seed_pr2[seed2_i].ref_pos_off_r;
						}
						else
						{
							end2_uc_f = 1;
						}
						q_rear2 = 0;
						q_n2 = 0;

						dmt2 = ali_exl;
						lv_dmt2 = lv_k2;

						s_r_o_l2 = seed_pr2[seed2_i].s_r_o_l;
						s_r_o_r2 = seed_pr2[seed2_i].s_r_o_r;
#ifdef	ALTER_DEBUG
						seed_length2 = seed_pr2[seed2_i].length;
#endif
						no2_tmp = seed2_i;
					}

#ifdef ALI_B_V_R
					//alignment

					//for end1
					if((end1_uc_f == 0) && ((dmt1 <= d_l1) && (dmt1 <= d_r1)))   // && (nuc1_f == 0)
					{
						if(nuc1_f == 0)
						{
							//lv_f1 = 0;
							pos_l = mat_pos1[tid][r_i] - max_extension_length - 1;
							re_d = pos_l & 0X1f;
							b_t_n_r = 32 - re_d;

							if(re_d != 0)
							{
								tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);
								memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars_1);

								for(rst_i = 0; rst_i < ref_copy_num_1 - 1; rst_i++)
								{
									tran_tmp = (ref_seq_tmp1[tid][rst_i] & bit_tran_re[re_d]);

									ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
									ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
									tran_tmp_p = tran_tmp;
								}

								ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
								ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

								//clear the lowest n bit
								ref_seq_tmp1[tid][rst_i] &= low_mask_1;
							}
							else
							{
								memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars_1);

								//clear the lowest n bit
								ref_seq_tmp1[tid][ref_copy_num_1 - 1] &= low_mask_1;
							}

							//exact match
							mis_c_n = 0;
#ifdef	QUAL_FILT
							mis_c_n_filt = 0;
#endif
							ref_tmp_ori = ref_seq_tmp1[tid][ref_copy_num_1 - 2];
							ref_seq_tmp1[tid][ref_copy_num_1 - 2] &= low_mask_1;

							for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num_1 - 1; rst_i++, read_b_i++)
							{
								xor_tmp = ref_seq_tmp1[tid][rst_i] ^ read_bit_1[tid][read_b_i];
#ifdef	QUAL_FILT
								mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								xor_tmp &= qual_filt_1[read_b_i];
								mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								if(mis_c_n_filt > max_mismatch1[tid])	break;
#else
								mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								mis_c_n_filt = mis_c_n;
								if(mis_c_n > max_mismatch1[tid])	break;
#endif
							}

							ref_seq_tmp1[tid][ref_copy_num_1 - 2] = ref_tmp_ori;

							if(mis_c_n_filt > max_mismatch1[tid])
							{
								if((cir_n == cir_fix_n) && (local_ksw))
								{
#ifdef	KSW_ALN_PAIR
									for(bit_char_i = s_r_o_l1, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + s_r_o_l1, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									ksw_extend(s_r_o_l1 + 1, read_char[tid], 33 + s_r_o_l1, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r1 - s_r_o_l1, &dm_l1, &tle, &gtle, &gscore, &max_off);

									for(bit_char_i = s_r_o_r1, read_b_i = 0; bit_char_i < read_length_1; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + s_r_o_r1, read_b_i = 0; bit_char_i < read_length_1 + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									ksw_extend(read_length_1 - s_r_o_r1, read_char[tid], 32 + read_length_1 - s_r_o_r1, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r1 - s_r_o_l1, &dm_r1, &tle, &gtle, &gscore, &max_off);

									dm1 = MAX_OP_SCORE_P1 - (dm_l1 + dm_r1); //read length cannot be more than 1024

									ld1 = s_r_o_l1;
									rd1 = s_r_o_r1;
#endif

								}
								else
								{
#ifdef SPLIT_LV
									pound_mis = 0;
									if(pound_pos_1_f >= s_r_o_r1)   //1
									{
#ifdef	POUND_MIS
										for(bit_char_i = pound_pos_1_f, bit_char_i_ref = pound_pos_1_f + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;

#else
										if(pound_pos_1_f < pound_pos_1_r)
										{
											read_pos_start_num = pound_pos_1_f >> 5;
											read_pos_end_num = (pound_pos_1_r - 1) >> 5;
											read_pos_re = (pound_pos_1_f & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = s_r_o_l1;
										lv_down_right = pound_pos_1_f;
										lv_down_left = s_r_o_r1;
									}
									else if(pound_pos_1_r <= s_r_o_l1 + 1)     //5
									{
#ifdef	POUND_MIS
										for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(pound_pos_1_r > 0)
										{
											read_pos_end_num = (pound_pos_1_r - 1) >> 5;
											read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

											if(read_pos_end_num == 0)
											{
												xor_tmp = (read_bit_1[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = read_bit_1[tid][0] ^ ref_seq_tmp1[tid][1];
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif

										lv_up_left = pound_pos_1_r;//
										lv_up_right = s_r_o_l1;
										lv_down_right = read_length_1;
										lv_down_left = s_r_o_r1;
									}
									else if((pound_pos_1_f <= s_r_o_l1 + 1) && (pound_pos_1_r >= s_r_o_r1))     //2
									{
#ifdef	POUND_MIS
										for(bit_char_i = pound_pos_1_f, bit_char_i_ref = pound_pos_1_f + 32; (bit_char_i <= s_r_o_l1) && (bit_char_i_ref <= s_r_o_l1 + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;


										for(bit_char_i = s_r_o_r1, bit_char_i_ref = s_r_o_r1 + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(pound_pos_1_f <= s_r_o_l1)
										{
											read_pos_start_num = pound_pos_1_f >> 5;
											read_pos_end_num = s_r_o_l1 >> 5;
											read_pos_re = (pound_pos_1_f & 0X1f) << 1;
											read_pos_end_re = (s_r_o_l1 & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}

										if(s_r_o_r1 < pound_pos_1_r)
										{
											read_pos_start_num = s_r_o_r1 >> 5;
											read_pos_end_num = (pound_pos_1_r - 1) >> 5;
											read_pos_re = (s_r_o_r1 & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = pound_pos_1_f - 1;
										lv_down_right = read_length_1;
										lv_down_left = pound_pos_1_r;
									}
									else if((pound_pos_1_f > s_r_o_l1 + 1) && (pound_pos_1_f < s_r_o_r1))     //3
									{
#ifdef	POUND_MIS
										for(bit_char_i = s_r_o_r1, bit_char_i_ref = s_r_o_r1 + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(s_r_o_r1 < pound_pos_1_r)
										{
											read_pos_start_num = s_r_o_r1 >> 5;
											read_pos_end_num = (pound_pos_1_r - 1) >> 5;
											read_pos_re = (s_r_o_r1 & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = s_r_o_l1;
										lv_down_right = read_length_1;
										lv_down_left = pound_pos_1_r;
									}
									else     //4
									{
#ifdef	POUND_MIS
										for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < s_r_o_l1) && (bit_char_i_ref < s_r_o_l1 + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(s_r_o_l1 > 0)
										{
											read_pos_end_num = (s_r_o_l1 - 1) >> 5;
											read_pos_end_re = ((s_r_o_l1 - 1) & 0X1f) << 1;

											if(read_pos_end_num == 0)
											{
												xor_tmp = (read_bit_1[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = read_bit_1[tid][0] ^ ref_seq_tmp1[tid][1];
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = -1;
										lv_down_right = read_length_1;
#ifdef	POUND_MODIFY
										lv_down_left = s_r_o_r1;
#else
										lv_down_left = s_r_o_l1;
#endif
									}
#ifdef	QUAL_FILT_LV

#ifdef	QUAL_FILT_LV_MIS
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#ifdef	S_DEBUG
									printf("extension lv1 %u: %u %d %d\n", mat_rc[tid][r_i], lv_dmt1, lv_up_right, lv_up_left);
#endif
									dm_l1 = computeEditDistance_mis(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid], qual_filt_lv_1_o + read_length_a1 - lv_up_right);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r1 = computeEditDistance_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid], qual_filt_lv_1 + lv_down_left);

									dm1 = dm_l1 + dm_r1 + pound_mis;

									dm_cir_1 = dm_l1 + dm_r1;// + pound_mis
#else
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l1 = computeEditDistance_misboth(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid], L_mis[tid], qual_filt_lv_1_o + read_length_a1 - lv_up_right, &mis_n1);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r1 = computeEditDistance_misboth(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid], L_mis[tid], qual_filt_lv_1 + lv_down_left, &mis_n2);

									dm1 = mis_n1 + mis_n2 + pound_mis;

									dm_cir_1 = dm_l1 + dm_r1;// + pound_mis

#endif
									ld1 = s_r_o_l1;
									rd1 = s_r_o_r1;

#else

#ifdef	LAST_CIRCLE_NOPOUND
									pound_mis = 0;
#endif
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l1 = computeEditDistance(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid]);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r1 = computeEditDistance(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid]);

									dm1 = dm_l1 + dm_r1 + pound_mis;

									dm_cir_1 = dm_l1 + dm_r1;// + pound_mis

									ld1 = s_r_o_l1;
									rd1 = s_r_o_r1;
#endif

#endif
								}

							}
							else
							{
								dm1 = mis_c_n;
								dm_cir_1 = mis_c_n;
								ld1 = 0;
								rd1 = 0;
								dm_l1 = 0;
								dm_r1 = 0;
							}

							//these two values could be different
							if((dm_l1 != -1) && (dm_r1 != -1))
							{
								if(dm1 < dmt1)	dmt1 = dm1 + 1;

								if(dm1 < lv_dmt1)	lv_dmt1 = dm1 + 1;

								if(dm1 < max_mismatch1[tid] - 1)	dmt1 = 0;
							}
							else
							{
								dm1 = MAX_EDIT_SCORE;
								dm_cir_1 = MAX_EDIT_SCORE;
							}
							nuc1_f = 1;
						}
					}
					else
					{
						pos_l = mat_pos1[tid][r_i] - max_extension_length - 1;
						re_d = pos_l & 0X1f;
						b_t_n_r = 32 - re_d;

						if(re_d != 0)
						{
							tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);

							memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars_1);

							for(rst_i = 0; rst_i < ref_copy_num_1 - 1; rst_i++)
							{
								tran_tmp = (ref_seq_tmp1[tid][rst_i] & bit_tran_re[re_d]);

								ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
								ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
								tran_tmp_p = tran_tmp;
							}

							ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
							ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

							//clear the lowest n bit
							ref_seq_tmp1[tid][rst_i] &= low_mask_1;
						}
						else
						{
							memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars_1);

							//clear the lowest n bit
							ref_seq_tmp1[tid][ref_copy_num_1 - 1] &= low_mask_1;
						}
						//trim the beginning and end of the current ref seq based on current minimum edit distance dm_t
						s_m_t = sub_mask1[tid][dmt1];
						ref_seq_tmp1[tid][0] &= bit_tran[dmt1];
						ref_seq_tmp1[tid][s_m_t] &= ex_d_mask1[tid][dmt1];

						//traverse and check whether there is an existing seq that is as same as current new ref seq
						c_m_f = 0;
						for(q_rear_i = q_rear1 - 1; q_rear_i >= 0; q_rear_i--)
						{
							ref_tmp_ori = cache_end1[q_rear_i][0];
							cache_end1[q_rear_i][0] &= bit_tran[dmt1];

							ref_tmp_ori2 = cache_end1[q_rear_i][s_m_t];
							cache_end1[q_rear_i][s_m_t] &= ex_d_mask1[tid][dmt1];

							cmp_re = memcmp(cache_end1[q_rear_i], ref_seq_tmp1[tid], (s_m_t + 1) << 3);

							cache_end1[q_rear_i][0] = ref_tmp_ori;
							cache_end1[q_rear_i][s_m_t] = ref_tmp_ori2;

							if(cmp_re == 0)
							{
								//deal with finding an alignment

								//lv_f1 = cache_lvf1[q_rear_i];
								dm1 = cache_dis1[q_rear_i];
								ld1 = cache_dml1[q_rear_i];
								rd1 = cache_dmr1[q_rear_i];

								dm_cir_1 = cache_dm_cir1[q_rear_i];

#ifdef	KSW_ALN_PAIR
								if(local_ksw)
								{
									dm_l1 = cache_kl1[q_rear_i];
									dm_r1 = cache_kr1[q_rear_i];
								}
#endif
								c_m_f = 1;
								break;

							}
						}

						if((q_n1 > MAX_Q_NUM) && (q_rear_i < 0))
						{
							for(q_rear_i = MAX_Q_NUM - 1; q_rear_i >= q_rear1; q_rear_i--)
							{
								ref_tmp_ori = cache_end1[q_rear_i][0];
								cache_end1[q_rear_i][0] &= bit_tran[dmt1];

								ref_tmp_ori2 = cache_end1[q_rear_i][s_m_t];
								cache_end1[q_rear_i][s_m_t] &= ex_d_mask1[tid][dmt1];

								cmp_re = memcmp(cache_end1[q_rear_i], ref_seq_tmp1[tid], (s_m_t + 1) << 3);

								cache_end1[q_rear_i][0] = ref_tmp_ori;
								cache_end1[q_rear_i][s_m_t] = ref_tmp_ori2;

								if(cmp_re == 0)
								{
									//deal with finding an alignment

									//lv_f1 = cache_lvf1[q_rear_i];
									dm1 = cache_dis1[q_rear_i];
									ld1 = cache_dml1[q_rear_i];
									rd1 = cache_dmr1[q_rear_i];

									dm_cir_1 = cache_dm_cir1[q_rear_i];

#ifdef	KSW_ALN_PAIR
									if(local_ksw)
									{
										dm_l1 = cache_kl1[q_rear_i];
										dm_r1 = cache_kr1[q_rear_i];
									}
#endif
									c_m_f = 1;
									break;

								}
							}
						}

						//do not find the seq in cache, exact match or lv and add into cache
						if(c_m_f == 0)
						{
							//exact match

							mis_c_n = 0;
#ifdef	QUAL_FILT
							mis_c_n_filt = 0;
#endif
							ref_tmp_ori = ref_seq_tmp1[tid][ref_copy_num_1 - 2];
							ref_seq_tmp1[tid][ref_copy_num_1 - 2] &= low_mask_1;

							for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num_1 - 1; rst_i++, read_b_i++)
							{
								xor_tmp = ref_seq_tmp1[tid][rst_i] ^ read_bit_1[tid][read_b_i];
#ifdef	QUAL_FILT
								mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								xor_tmp &= qual_filt_1[read_b_i];
								mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								if(mis_c_n_filt > max_mismatch1[tid])	break;
#else
								mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								mis_c_n_filt = mis_c_n;
								if(mis_c_n > max_mismatch1[tid])	break;
#endif
								//mis_c_n += popcount_3(ref_seq_tmp1[rst_i] ^ read_bit_1[read_b_i]);
							}

							ref_seq_tmp1[tid][ref_copy_num_1 - 2] = ref_tmp_ori;

							//lv
							if(mis_c_n_filt > max_mismatch1[tid])
							{
								if((cir_n == cir_fix_n) && (local_ksw))
								{
#ifdef	KSW_ALN_PAIR
									for(bit_char_i = s_r_o_l1, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + s_r_o_l1, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									ksw_extend(s_r_o_l1 + 1, read_char[tid], 33 + s_r_o_l1, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r1 - s_r_o_l1, &dm_l1, &tle, &gtle, &gscore, &max_off);

									for(bit_char_i = s_r_o_r1, read_b_i = 0; bit_char_i < read_length_1; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + s_r_o_r1, read_b_i = 0; bit_char_i < read_length_1 + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									ksw_extend(read_length_1 - s_r_o_r1, read_char[tid], 32 + read_length_1 - s_r_o_r1, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r1 - s_r_o_l1, &dm_r1, &tle, &gtle, &gscore, &max_off);

									dm1 = MAX_OP_SCORE_P1 - (dm_l1 + dm_r1); //read length cannot be more than 1024

									ld1 = s_r_o_l1;
									rd1 = s_r_o_r1;
#endif
								}
								else
								{

#ifdef SPLIT_LV
									pound_mis = 0;
									if(pound_pos_1_f >= s_r_o_r1)   //1
									{
#ifdef	POUND_MIS
										for(bit_char_i = pound_pos_1_f, bit_char_i_ref = pound_pos_1_f + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(pound_pos_1_f < pound_pos_1_r)
										{
											read_pos_start_num = pound_pos_1_f >> 5;
											read_pos_end_num = (pound_pos_1_r - 1) >> 5;
											read_pos_re = (pound_pos_1_f & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = s_r_o_l1;
										lv_down_right = pound_pos_1_f;
										lv_down_left = s_r_o_r1;
									}
									else if(pound_pos_1_r <= s_r_o_l1 + 1)     //5
									{
#ifdef	POUND_MIS
										for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else
										if(pound_pos_1_r > 0)
										{
											read_pos_end_num = (pound_pos_1_r - 1) >> 5;
											read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

											if(read_pos_end_num == 0)
											{
												xor_tmp = (read_bit_1[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = read_bit_1[tid][0] ^ ref_seq_tmp1[tid][1];
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif

										lv_up_left = pound_pos_1_r;//
										lv_up_right = s_r_o_l1;
										lv_down_right = read_length_1;
										lv_down_left = s_r_o_r1;
									}
									else if((pound_pos_1_f <= s_r_o_l1 + 1) && (pound_pos_1_r >= s_r_o_r1))     //2
									{
#ifdef	POUND_MIS
										for(bit_char_i = pound_pos_1_f, bit_char_i_ref = pound_pos_1_f + 32; (bit_char_i <= s_r_o_l1) && (bit_char_i_ref <= s_r_o_l1 + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;


										for(bit_char_i = s_r_o_r1, bit_char_i_ref = s_r_o_r1 + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(pound_pos_1_f <= s_r_o_l1)
										{
											read_pos_start_num = pound_pos_1_f >> 5;
											read_pos_end_num = s_r_o_l1 >> 5;
											read_pos_re = (pound_pos_1_f & 0X1f) << 1;
											read_pos_end_re = (s_r_o_l1 & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}

										if(s_r_o_r1 < pound_pos_1_r)
										{
											read_pos_start_num = s_r_o_r1 >> 5;
											read_pos_end_num = (pound_pos_1_r - 1) >> 5;
											read_pos_re = (s_r_o_r1 & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif

										lv_up_left = 0;
										lv_up_right = pound_pos_1_f - 1;
										lv_down_right = read_length_1;
										lv_down_left = pound_pos_1_r;
									}
									else if((pound_pos_1_f > s_r_o_l1 + 1) && (pound_pos_1_f < s_r_o_r1))     //3
									{
#ifdef	POUND_MIS
										for(bit_char_i = s_r_o_r1, bit_char_i_ref = s_r_o_r1 + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(s_r_o_r1 < pound_pos_1_r)
										{
											read_pos_start_num = s_r_o_r1 >> 5;
											read_pos_end_num = (pound_pos_1_r - 1) >> 5;
											read_pos_re = (s_r_o_r1 & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = s_r_o_l1;
										lv_down_right = read_length_1;
										lv_down_left = pound_pos_1_r;
									}
									else     //pound_pos_1_r > s_r_o_l1) && (pound_pos_1_r < s_r_o_r1)//4
									{
#ifdef	POUND_MIS
										for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < s_r_o_l1) && (bit_char_i_ref < s_r_o_l1 + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(s_r_o_l1 > 0)
										{
											read_pos_end_num = (s_r_o_l1 - 1) >> 5;
											read_pos_end_re = ((s_r_o_l1 - 1) & 0X1f) << 1;

											if(read_pos_end_num == 0)
											{
												xor_tmp = (read_bit_1[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = read_bit_1[tid][0] ^ ref_seq_tmp1[tid][1];
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = -1;
										lv_down_right = read_length_1;
#ifdef	POUND_MODIFY
										lv_down_left = s_r_o_r1;
#else
										lv_down_left = s_r_o_l1;
#endif
									}
#ifdef	QUAL_FILT_LV

#ifdef	QUAL_FILT_LV_MIS
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l1 = computeEditDistance_mis(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid], qual_filt_lv_1_o + read_length_a1 - lv_up_right);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r1 = computeEditDistance_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid], qual_filt_lv_1 + lv_down_left);

									dm1 = dm_l1 + dm_r1 + pound_mis;

									dm_cir_1 = dm_l1 + dm_r1;// + pound_mis
#else
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l1 = computeEditDistance_misboth(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid], L_mis[tid], qual_filt_lv_1_o + read_length_a1 - lv_up_right, &mis_n1);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r1 = computeEditDistance_misboth(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid], L_mis[tid], qual_filt_lv_1 + lv_down_left, &mis_n2);

									dm1 = mis_n1 + mis_n2 + pound_mis;

									dm_cir_1 = dm_l1 + dm_r1;// + pound_mis

#endif
									ld1 = s_r_o_l1;
									rd1 = s_r_o_r1;
#else


#ifdef	LAST_CIRCLE_NOPOUND
									pound_mis = 0;
#endif
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l1 = computeEditDistance(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid]);

									//lv_x_r = (dm_l1 == -1) ? lv_dmt1 : lv_dmt1 - dm_l1;

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r1 = computeEditDistance(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid]);

									dm1 = dm_l1 + dm_r1 + pound_mis;

									dm_cir_1 = dm_l1 + dm_r1;// + pound_mis

									ld1 = s_r_o_l1;
									rd1 = s_r_o_r1;
#endif

#endif
								}

							}
							else
							{
								dm1 = mis_c_n;
								dm_cir_1 = mis_c_n;
								ld1 = 0;
								rd1 = 0;
								dm_l1 = 0;
								dm_r1 = 0;
							}

							//these two values could be different
							if((dm_l1 != -1) && (dm_r1 != -1))
							{
								if(dm1 < dmt1)	dmt1 = dm1 + 1;

								if(dm1 < lv_dmt1)	lv_dmt1 = dm1 + 1;

								if(dm1 < max_mismatch1[tid] - 1)	dmt1 = 0;
							}
							else
							{
								dm1 = MAX_EDIT_SCORE;
								dm_cir_1 = MAX_EDIT_SCORE;
							}

							//add the ref sequence at the end of queue
							memcpy(cache_end1[q_rear1], ref_seq_tmp1[tid], ref_copy_num_chars_1);
							cache_dis1[q_rear1] = dm1;
							cache_dml1[q_rear1] = ld1;
							cache_dmr1[q_rear1] = rd1;

							cache_dm_cir1[q_rear1] = dm_cir_1;

#ifdef	KSW_ALN_PAIR
							if(local_ksw)
							{
								cache_kl1[q_rear1] = dm_l1;
								cache_kr1[q_rear1] = dm_r1;
							}
#endif
							q_rear1 = ((q_rear1 + 1) & 0X1f);
							++q_n1;

							//add edit distance
						}
					}

					//for end2
					if((end2_uc_f == 0) && ((dmt2 <= d_l2) && (dmt2 <= d_r2)))   // && (nuc2_f == 0)
					{
						if(nuc2_f == 0)
						{
							pos_l = mat_pos2[tid][r_i] - max_extension_length - 1;
							re_d = pos_l & 0X1f;
							b_t_n_r = 32 - re_d;

							if(re_d != 0)
							{
								tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);
								memcpy(ref_seq_tmp2[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars_2);

								for(rst_i = 0; rst_i < ref_copy_num_2 - 1; rst_i++)
								{
									tran_tmp = (ref_seq_tmp2[tid][rst_i] & bit_tran_re[re_d]);

									ref_seq_tmp2[tid][rst_i] >>= (b_t_n_r << 1);
									ref_seq_tmp2[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
									tran_tmp_p = tran_tmp;
								}
								ref_seq_tmp2[tid][rst_i] >>= (b_t_n_r << 1);
								ref_seq_tmp2[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

								//clear the lowest n bit
								ref_seq_tmp2[tid][rst_i] &= low_mask_2;

							}
							else
							{
								memcpy(ref_seq_tmp2[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars_2);

								//clear the lowest n bit
								ref_seq_tmp2[tid][ref_copy_num_2 - 1] &= low_mask_2;
							}

							//exact match

							mis_c_n = 0;
#ifdef	QUAL_FILT
							mis_c_n_filt = 0;
#endif
							ref_tmp_ori = ref_seq_tmp2[tid][ref_copy_num_2 - 2];
							ref_seq_tmp2[tid][ref_copy_num_2 - 2] &= low_mask_2;

							for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num_2 - 1; rst_i++, read_b_i++)
							{
								xor_tmp = ref_seq_tmp2[tid][rst_i] ^ read_bit_2[tid][read_b_i];
#ifdef	QUAL_FILT
								mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								xor_tmp &= qual_filt_2[read_b_i];
								mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								if(mis_c_n_filt > max_mismatch2[tid])	break;
#else
								mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								mis_c_n_filt = mis_c_n;
								if(mis_c_n > max_mismatch2[tid])	break;
#endif
							}

							ref_seq_tmp2[tid][ref_copy_num_2 - 2] = ref_tmp_ori;

							//lv
							if(mis_c_n_filt > max_mismatch2[tid])
							{
								if((cir_n == cir_fix_n) && (local_ksw))
								{
#ifdef	KSW_ALN_PAIR
									for(bit_char_i = s_r_o_l2, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + s_r_o_l2, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									ksw_extend(s_r_o_l2 + 1, read_char[tid], 33 + s_r_o_l2, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r2 - s_r_o_l2, &dm_l2, &tle, &gtle, &gscore, &max_off);

									for(bit_char_i = s_r_o_r2, read_b_i = 0; bit_char_i < read_length_2; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + s_r_o_r2, read_b_i = 0; bit_char_i < read_length_2 + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									ksw_extend(read_length_2 - s_r_o_r2, read_char[tid], 32 + read_length_2 - s_r_o_r2, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r2 - s_r_o_l2, &dm_r2, &tle, &gtle, &gscore, &max_off);

									dm2 = MAX_OP_SCORE_P2 - (dm_l2 + dm_r2); //read length cannot be more than 1024

									ld2 = s_r_o_l2;
									rd2 = s_r_o_r2;
#endif

								}
								else
								{


#ifdef SPLIT_LV
									pound_mis = 0;
									if(pound_pos_2_f >= s_r_o_r2)   //1
									{
#ifdef	POUND_MIS
										for(bit_char_i = pound_pos_2_f, bit_char_i_ref = pound_pos_2_f + 32; (bit_char_i < pound_pos_2_r) && (bit_char_i_ref < pound_pos_2_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(pound_pos_2_f < pound_pos_2_r)
										{
											read_pos_start_num = pound_pos_2_f >> 5;
											read_pos_end_num = (pound_pos_2_r - 1) >> 5;
											read_pos_re = (pound_pos_2_f & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_2_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_2[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_2[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = s_r_o_l2;
										lv_down_right = pound_pos_2_f;
										lv_down_left = s_r_o_r2;
									}
									else if(pound_pos_2_r <= s_r_o_l2 + 1)     //5
									{
#ifdef	POUND_MIS
										for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < pound_pos_2_r) && (bit_char_i_ref < pound_pos_2_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else
										if(pound_pos_2_r > 0)
										{
											read_pos_end_num = (pound_pos_2_r - 1) >> 5;
											read_pos_end_re = ((pound_pos_2_r - 1) & 0X1f) << 1;

											if(read_pos_end_num == 0)
											{
												xor_tmp = (read_bit_2[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = read_bit_2[tid][0] ^ ref_seq_tmp2[tid][1];
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif

										lv_up_left = pound_pos_2_r;
										lv_up_right = s_r_o_l2;
										lv_down_right = read_length_2;
										lv_down_left = s_r_o_r2;
									}
									else if((pound_pos_2_f <= s_r_o_l2 + 1) && (pound_pos_2_r >= s_r_o_r2))     //2
									{
#ifdef	POUND_MIS
										for(bit_char_i = pound_pos_2_f, bit_char_i_ref = pound_pos_2_f + 32; (bit_char_i <= s_r_o_l2) && (bit_char_i_ref <= s_r_o_l2 + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;


										for(bit_char_i = s_r_o_r2, bit_char_i_ref = s_r_o_r2 + 32; (bit_char_i < pound_pos_2_r) && (bit_char_i_ref < pound_pos_2_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(pound_pos_2_f <= s_r_o_l2)
										{
											read_pos_start_num = pound_pos_2_f >> 5;
											read_pos_end_num = s_r_o_l2 >> 5;
											read_pos_re = (pound_pos_2_f & 0X1f) << 1;
											read_pos_end_re = (s_r_o_l2 & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_2[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_2[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}

										if(s_r_o_r2 < pound_pos_2_r)
										{
											read_pos_start_num = s_r_o_r2 >> 5;
											read_pos_end_num = (pound_pos_2_r - 1) >> 5;
											read_pos_re = (s_r_o_r2 & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_2_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_2[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_2[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif

										lv_up_left = 0;
										lv_up_right = pound_pos_2_f - 1;
										lv_down_right = read_length_2;
										lv_down_left = pound_pos_2_r;
									}
									else if((pound_pos_2_f > s_r_o_l2 + 1) && (pound_pos_2_f < s_r_o_r2))     //3
									{
#ifdef	POUND_MIS
										for(bit_char_i = s_r_o_r2, bit_char_i_ref = s_r_o_r2 + 32; (bit_char_i < pound_pos_2_r) && (bit_char_i_ref < pound_pos_2_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(s_r_o_r2 < pound_pos_2_r)
										{
											read_pos_start_num = s_r_o_r2 >> 5;
											read_pos_end_num = (pound_pos_2_r - 1) >> 5;
											read_pos_re = (s_r_o_r2 & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_2_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_2[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_2[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = s_r_o_l2;
										lv_down_right = read_length_2;
										lv_down_left = pound_pos_2_r;
									}
									else     //4
									{
#ifdef	POUND_MIS
										for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < s_r_o_l2) && (bit_char_i_ref < s_r_o_l2 + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(s_r_o_l2 > 0)
										{
											read_pos_end_num = (s_r_o_l2 - 1) >> 5;
											read_pos_end_re = ((s_r_o_l2 - 1) & 0X1f) << 1;

											if(read_pos_end_num == 0)
											{
												xor_tmp = (read_bit_2[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = read_bit_2[tid][0] ^ ref_seq_tmp2[tid][1];
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = -1;
										lv_down_right = read_length_2;
#ifdef	POUND_MODIFY
										lv_down_left = s_r_o_r2;
#else
										lv_down_left = s_r_o_l2;
#endif
									}
#ifdef	QUAL_FILT_LV

#ifdef	QUAL_FILT_LV_MIS
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l2 = computeEditDistance_mis(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt2, L[tid], qual_filt_lv_2_o + read_length_a2 - lv_up_right);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r2 = computeEditDistance_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt2, L[tid], qual_filt_lv_2 + lv_down_left);

									dm2 = dm_l2 + dm_r2 + pound_mis;

									dm_cir_2 = dm_l2 + dm_r2;// + pound_mis

#else
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l2 = computeEditDistance_misboth(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt2, L[tid], L_mis[tid], qual_filt_lv_2_o + read_length_a2 - lv_up_right, &mis_n1);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r2 = computeEditDistance_misboth(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt2, L[tid], L_mis[tid], qual_filt_lv_2 + lv_down_left, &mis_n2);

									dm2 = mis_n1 + mis_n2 + pound_mis;

									dm_cir_2 = dm_l2 + dm_r2;// + pound_mis

#endif
									ld2 = s_r_o_l2;
									rd2 = s_r_o_r2;
#else

#ifdef	LAST_CIRCLE_NOPOUND
									pound_mis = 0;
#endif
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l2 = computeEditDistance(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt2, L[tid]);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r2 = computeEditDistance(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt2, L[tid]);

									dm2 = dm_l2 + dm_r2 + pound_mis;

									dm_cir_2 = dm_l2 + dm_r2;// + pound_mis

									ld2 = s_r_o_l2;
									rd2 = s_r_o_r2;
#endif

#endif
								}
							}
							else
							{
								dm2 = mis_c_n;
								dm_cir_2 = mis_c_n;
								ld2 = 0;
								rd2 = 0;
								dm_l2 = 0;
								dm_r2 = 0;
							}

							//these two values could be different
							if((dm_l2 != -1) && (dm_r2 != -1))
							{
								if(dm2 < dmt2)	dmt2 = dm2 + 1;

								if(dm2 < lv_dmt2)	lv_dmt2 = dm2 + 1;

								if(dm2 < max_mismatch2[tid] - 1)	dmt2 = 0;
							}
							else
							{
								dm2 = MAX_EDIT_SCORE;
								dm_cir_2 = MAX_EDIT_SCORE;
							}
							nuc2_f = 1;
						}
					}
					else
					{
						pos_l = mat_pos2[tid][r_i] - max_extension_length - 1;
						re_d = pos_l & 0X1f;
						b_t_n_r = 32 - re_d;

						if(re_d != 0)
						{
							tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);
							memcpy(ref_seq_tmp2[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars_2);

							for(rst_i = 0; rst_i < ref_copy_num_2 - 1; rst_i++)
							{
								tran_tmp = (ref_seq_tmp2[tid][rst_i] & bit_tran_re[re_d]);

								ref_seq_tmp2[tid][rst_i] >>= (b_t_n_r << 1);
								ref_seq_tmp2[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
								tran_tmp_p = tran_tmp;
							}

							ref_seq_tmp2[tid][rst_i] >>= (b_t_n_r << 1);
							ref_seq_tmp2[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

							//clear the lowest n bit
							ref_seq_tmp2[tid][rst_i] &= low_mask_2;
						}
						else
						{
							memcpy(ref_seq_tmp2[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars_2);

							//clear the lowest n bit
							ref_seq_tmp2[tid][ref_copy_num_2 - 1] &= low_mask_2;
						}

						//trim the beginning and end of the current ref seq based on current minimum edit distance dm_t
						s_m_t = sub_mask2[tid][dmt2];
						ref_seq_tmp2[tid][0] &= bit_tran[dmt2];
						ref_seq_tmp2[tid][s_m_t] &= ex_d_mask2[tid][dmt2];

						//traverse and check whether there is an existing seq that is as same as current new ref seq
						c_m_f = 0;
						for(q_rear_i = q_rear2 - 1; q_rear_i >= 0; q_rear_i--)
						{
							ref_tmp_ori = cache_end2[q_rear_i][0];
							cache_end2[q_rear_i][0] &= bit_tran[dmt2];

							ref_tmp_ori2 = cache_end2[q_rear_i][s_m_t];
							cache_end2[q_rear_i][s_m_t] &= ex_d_mask2[tid][dmt2];

							cmp_re = memcmp(cache_end2[q_rear_i], ref_seq_tmp2[tid], (s_m_t + 1) << 3);

							cache_end2[q_rear_i][0] = ref_tmp_ori;
							cache_end2[q_rear_i][s_m_t] = ref_tmp_ori2;

							if(cmp_re == 0)
							{
								//deal with finding an alignment
								dm2 = cache_dis2[q_rear_i];
								ld2 = cache_dml2[q_rear_i];
								rd2 = cache_dmr2[q_rear_i];

								dm_cir_2 = cache_dm_cir2[q_rear_i];
#ifdef	KSW_ALN_PAIR
								if(local_ksw)
								{
									dm_l2 = cache_kl2[q_rear_i];
									dm_r2 = cache_kr2[q_rear_i];
								}
#endif

								c_m_f = 1;
								break;
							}
						}

						if((q_n2 > MAX_Q_NUM) && (q_rear_i < 0))
						{
							for(q_rear_i = MAX_Q_NUM - 1; q_rear_i >= q_rear2; q_rear_i--)
							{
								ref_tmp_ori = cache_end2[q_rear_i][0];
								cache_end2[q_rear_i][0] &= bit_tran[dmt2];

								ref_tmp_ori2 = cache_end2[q_rear_i][s_m_t];
								cache_end2[q_rear_i][s_m_t] &= ex_d_mask2[tid][dmt2];

								cmp_re = memcmp(cache_end2[q_rear_i], ref_seq_tmp2[tid], (s_m_t + 1) << 3);

								cache_end2[q_rear_i][0] = ref_tmp_ori;
								cache_end2[q_rear_i][s_m_t] = ref_tmp_ori2;

								if(cmp_re == 0)
								{
									//deal with finding an alignment
									dm2 = cache_dis2[q_rear_i];
									ld2 = cache_dml2[q_rear_i];
									rd2 = cache_dmr2[q_rear_i];

									dm_cir_2 = cache_dm_cir2[q_rear_i];
#ifdef	KSW_ALN_PAIR
									if(local_ksw)
									{
										dm_l2 = cache_kl2[q_rear_i];
										dm_r2 = cache_kr2[q_rear_i];
									}
#endif
									c_m_f = 1;
									break;
								}
							}
						}

						//do not find the seq in cache, exact match or lv and add into cache
						if(c_m_f == 0)
						{
							//exact match
							mis_c_n = 0;
#ifdef	QUAL_FILT
							mis_c_n_filt = 0;
#endif
							ref_tmp_ori = ref_seq_tmp2[tid][ref_copy_num_2 - 2];
							ref_seq_tmp2[tid][ref_copy_num_2 - 2] &= low_mask_2;

							for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num_2 - 1; rst_i++, read_b_i++)
							{
								xor_tmp = ref_seq_tmp2[tid][rst_i] ^ read_bit_2[tid][read_b_i];
#ifdef	QUAL_FILT
								mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								xor_tmp &= qual_filt_2[read_b_i];
								mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								if(mis_c_n_filt > max_mismatch2[tid])	break;
#else
								mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								mis_c_n_filt = mis_c_n;
								if(mis_c_n > max_mismatch2[tid])	break;
#endif
							}

							ref_seq_tmp2[tid][ref_copy_num_2 - 2] = ref_tmp_ori;

							//lv
							if(mis_c_n_filt > max_mismatch2[tid])
							{
								if((cir_n == cir_fix_n) && (local_ksw))
								{
#ifdef	KSW_ALN_PAIR
									for(bit_char_i = s_r_o_l2, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + s_r_o_l2, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									ksw_extend(s_r_o_l2 + 1, read_char[tid], 33 + s_r_o_l2, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r2 - s_r_o_l2, &dm_l2, &tle, &gtle, &gscore, &max_off);

									for(bit_char_i = s_r_o_r2, read_b_i = 0; bit_char_i < read_length_2; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + s_r_o_r2, read_b_i = 0; bit_char_i < read_length_2 + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									ksw_extend(read_length_2 - s_r_o_r2, read_char[tid], 32 + read_length_2 - s_r_o_r2, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r2 - s_r_o_l2, &dm_r2, &tle, &gtle, &gscore, &max_off);

									dm2 = MAX_OP_SCORE_P2 - (dm_l2 + dm_r2); //read length cannot be more than 1024

									ld2 = s_r_o_l2;
									rd2 = s_r_o_r2;
#endif

								}
								else
								{
#ifdef SPLIT_LV
									pound_mis = 0;
									if(pound_pos_2_f >= s_r_o_r2)   //1
									{
#ifdef	POUND_MIS
										for(bit_char_i = pound_pos_2_f, bit_char_i_ref = pound_pos_2_f + 32; (bit_char_i < pound_pos_2_r) && (bit_char_i_ref < pound_pos_2_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(pound_pos_2_f < pound_pos_2_r)
										{
											read_pos_start_num = pound_pos_2_f >> 5;
											read_pos_end_num = (pound_pos_2_r - 1) >> 5;
											read_pos_re = (pound_pos_2_f & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_2_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_2[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_2[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = s_r_o_l2;
										lv_down_right = pound_pos_2_f;
										lv_down_left = s_r_o_r2;
									}
									else if(pound_pos_2_r <= s_r_o_l2 + 1)     //5
									{
#ifdef	POUND_MIS
										for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < pound_pos_2_r) && (bit_char_i_ref < pound_pos_2_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else
										if(pound_pos_2_r > 0)
										{
											read_pos_end_num = (pound_pos_2_r - 1) >> 5;
											read_pos_end_re = ((pound_pos_2_r - 1) & 0X1f) << 1;

											if(read_pos_end_num == 0)
											{
												xor_tmp = (read_bit_2[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = read_bit_2[tid][0] ^ ref_seq_tmp2[tid][1];
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif

										lv_up_left = pound_pos_2_r;
										lv_up_right = s_r_o_l2;
										lv_down_right = read_length_2;
										lv_down_left = s_r_o_r2;
									}
									else if((pound_pos_2_f <= s_r_o_l2 + 1) && (pound_pos_2_r >= s_r_o_r2))     //2
									{
#ifdef	POUND_MIS
										for(bit_char_i = pound_pos_2_f, bit_char_i_ref = pound_pos_2_f + 32; (bit_char_i <= s_r_o_l2) && (bit_char_i_ref <= s_r_o_l2 + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;


										for(bit_char_i = s_r_o_r2, bit_char_i_ref = s_r_o_r2 + 32; (bit_char_i < pound_pos_2_r) && (bit_char_i_ref < pound_pos_2_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(pound_pos_2_f <= s_r_o_l2)
										{
											read_pos_start_num = pound_pos_2_f >> 5;
											read_pos_end_num = s_r_o_l2 >> 5;
											read_pos_re = (pound_pos_2_f & 0X1f) << 1;
											read_pos_end_re = (s_r_o_l2 & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_2[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_2[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}

										if(s_r_o_r2 < pound_pos_2_r)
										{
											read_pos_start_num = s_r_o_r2 >> 5;
											read_pos_end_num = (pound_pos_2_r - 1) >> 5;
											read_pos_re = (s_r_o_r2 & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_2_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_2[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_2[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif

										lv_up_left = 0;
										lv_up_right = pound_pos_2_f - 1;
										lv_down_right = read_length_2;
										lv_down_left = pound_pos_2_r;
									}
									else if((pound_pos_2_f > s_r_o_l2 + 1) && (pound_pos_2_f < s_r_o_r2))     //3
									{
#ifdef	POUND_MIS
										for(bit_char_i = s_r_o_r2, bit_char_i_ref = s_r_o_r2 + 32; (bit_char_i < pound_pos_2_r) && (bit_char_i_ref < pound_pos_2_r + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(s_r_o_r2 < pound_pos_2_r)
										{
											read_pos_start_num = s_r_o_r2 >> 5;
											read_pos_end_num = (pound_pos_2_r - 1) >> 5;
											read_pos_re = (s_r_o_r2 & 0X1f) << 1;
											read_pos_end_re = ((pound_pos_2_r - 1) & 0X1f) << 1;

											if(read_pos_start_num == read_pos_end_num)
											{
												xor_tmp = ((read_bit_2[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = (read_bit_2[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp2[tid][read_pos_start_num + 1] << read_pos_re);
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = s_r_o_l2;
										lv_down_right = read_length_2;
										lv_down_left = pound_pos_2_r;
									}
									else     //4
									{
#ifdef	POUND_MIS
										for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < s_r_o_l2) && (bit_char_i_ref < s_r_o_l2 + 32); bit_char_i++, bit_char_i_ref++)
											if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp2[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
												pound_mis++;
#else

										if(s_r_o_l2 > 0)
										{
											read_pos_end_num = (s_r_o_l2 - 1) >> 5;
											read_pos_end_re = ((s_r_o_l2 - 1) & 0X1f) << 1;

											if(read_pos_end_num == 0)
											{
												xor_tmp = (read_bit_2[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

											}
											else
											{
												xor_tmp = read_bit_2[tid][0] ^ ref_seq_tmp2[tid][1];
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
												{
													xor_tmp = read_bit_2[tid][rst_i] ^ ref_seq_tmp2[tid][rst_i_1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}

												xor_tmp = (read_bit_2[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp2[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
												pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
											}
										}
#endif
										lv_up_left = 0;
										lv_up_right = -1;
										lv_down_right = read_length_2;
#ifdef	POUND_MODIFY
										lv_down_left = s_r_o_r2;
#else
										lv_down_left = s_r_o_l2;
#endif
									}
#ifdef	QUAL_FILT_LV

#ifdef	QUAL_FILT_LV_MIS
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#ifdef	S_DEBUG
									printf("extension lv2 %u: %u %d %d\n", mat_rc[tid][r_i], lv_dmt2, lv_up_right, lv_up_left);
#endif
									dm_l2 = computeEditDistance_mis(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt2, L[tid], qual_filt_lv_2_o + read_length_a2 - lv_up_right);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r2 = computeEditDistance_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt2, L[tid], qual_filt_lv_2 + lv_down_left);

									dm2 = dm_l2 + dm_r2 + pound_mis;

									dm_cir_2 = dm_l2 + dm_r2;// + pound_mis
#else
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l2 = computeEditDistance_misboth(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt2, L[tid], L_mis[tid], qual_filt_lv_2_o + read_length_a2 - lv_up_right, &mis_n1);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r2 = computeEditDistance_misboth(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt2, L[tid], L_mis[tid], qual_filt_lv_2 + lv_down_left, &mis_n2);

									dm2 = mis_n1 + mis_n2 + pound_mis;

									dm_cir_2 = dm_l2 + dm_r2;// + pound_mis
#endif
									ld2 = s_r_o_l2;
									rd2 = s_r_o_r2;
#else


#ifdef	LAST_CIRCLE_NOPOUND
									pound_mis = 0;
#endif
									for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_l2 = computeEditDistance(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt2, L[tid]);

									for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
										read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
										ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

									dm_r2 = computeEditDistance(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt2, L[tid]);

									dm2 = dm_l2 + dm_r2 + pound_mis;

									dm_cir_2 = dm_l2 + dm_r2;// + pound_mis

									ld2 = s_r_o_l2;
									rd2 = s_r_o_r2;
#endif

#endif
								}
							}
							else
							{
								dm2 = mis_c_n;
								dm_cir_2 = mis_c_n;
								ld2 = 0;
								rd2 = 0;
								dm_l2 = 0;
								dm_r2 = 0;
							}

							//these two values could be different
							if((dm_l2 != -1) && (dm_r2 != -1))
							{
								if(dm2 < dmt2)	dmt2 = dm2 + 1;

								if(dm2 < lv_dmt2)	lv_dmt2 = dm2 + 1;

								if(dm2 < max_mismatch2[tid] - 1)	dmt2 = 0;
							}
							else
							{
								dm2 = MAX_EDIT_SCORE;
								dm_cir_2 = MAX_EDIT_SCORE;
							}

							//add the ref sequence at the end of queue
							memcpy(cache_end2[q_rear2], ref_seq_tmp2[tid], ref_copy_num_chars_2);
							cache_dis2[q_rear2] = dm2;
							//cache_lvf1[q_rear1] = lv_f1;
							cache_dml2[q_rear2] = ld2;
							cache_dmr2[q_rear2] = rd2;

							cache_dm_cir2[q_rear2] = dm_cir_2;

#ifdef	KSW_ALN_PAIR
							if(local_ksw)
							{
								cache_kl2[q_rear2] = dm_l2;
								cache_kr2[q_rear2] = dm_r2;
							}
#endif
							q_rear2 = ((q_rear2 + 1) & 0X1f);
							++q_n2;

							//add edit distance
						}
					}
#endif
					dm_cir = dm_cir_1 + dm_cir_2;
					if(dm_cir < dm_cir_min)	dm_cir_min = dm_cir;

					dm12 = dm1 + dm2;

					if(dm12 < dm_op[tid])
					{
#ifdef	DM_COPY_PAIR
						for(dm_i = 0; dm_i < v_cnt; dm_i++)
						{
							ops_vector_pos1[tid][dm_i] = op_vector_pos1[tid][dm_i];

							ops_dm_l1[tid][dm_i] = op_dm_l1[tid][dm_i];
							ops_dm_r1[tid][dm_i] = op_dm_r1[tid][dm_i];
#ifdef	MAPPING_QUALITY

							memcpy(ops_vector_seq1[tid][dm_i], op_vector_seq1[tid][dm_i], ref_copy_num_chars1);
#else
							if(!((op_dm_l1[tid][dm_i] == 0) && (op_dm_r1[tid][dm_i] == 0)))
								memcpy(ops_vector_seq1[tid][dm_i], op_vector_seq1[tid][dm_i], ref_copy_num_chars1);
#endif
							ops_dm_ex1[tid][dm_i] = op_dm_ex1[tid][dm_i];

							ops_vector_pos2[tid][dm_i] = op_vector_pos2[tid][dm_i];

							ops_dm_l2[tid][dm_i] = op_dm_l2[tid][dm_i];
							ops_dm_r2[tid][dm_i] = op_dm_r2[tid][dm_i];
#ifdef	MAPPING_QUALITY
							memcpy(ops_vector_seq2[tid][dm_i], op_vector_seq2[tid][dm_i], ref_copy_num_chars2);
#else
							if(!((op_dm_l2[tid][dm_i] == 0) && (op_dm_r2[tid][dm_i] == 0)))
								memcpy(ops_vector_seq2[tid][dm_i], op_vector_seq2[tid][dm_i], ref_copy_num_chars2);
#endif
							ops_dm_ex2[tid][dm_i] = op_dm_ex2[tid][dm_i];

#ifdef	KSW_ALN_PAIR
							if(local_ksw)
							{
								ops_dm_kl1[tid][dm_i] = op_dm_kl1[tid][dm_i];
								ops_dm_kr1[tid][dm_i] = op_dm_kr1[tid][dm_i];
								ops_dm_kl2[tid][dm_i] = op_dm_kl2[tid][dm_i];
								ops_dm_kr2[tid][dm_i] = op_dm_kr2[tid][dm_i];
							}
#endif
							ops_rc[tid][dm_i] = op_rc[tid][dm_i];
						}
						vs_cnt = v_cnt;
						dm_ops[tid] = dm_op[tid];
#endif
						v_cnt = 0;
						if(mat_rc[tid][r_i] == 0)
						{
							op_vector_pos1[tid][v_cnt] = mat_pos1[tid][r_i];
							op_vector_pos2[tid][v_cnt] = mat_pos2[tid][r_i];
#ifdef	MAPPING_QUALITY
							memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
							memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#else
							if(!((ld1 == 0) && (rd1 == 0)))
								memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
							if(!((ld2 == 0) && (rd2 == 0)))
								memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#endif
							op_dm_l1[tid][v_cnt] = ld1;
							op_dm_r1[tid][v_cnt] = rd1;
							op_dm_l2[tid][v_cnt] = ld2;
							op_dm_r2[tid][v_cnt] = rd2;

							op_dm_ex1[tid][v_cnt] = dm1;
							op_dm_ex2[tid][v_cnt] = dm2;

#ifdef	KSW_ALN_PAIR
							if(local_ksw)
							{
								op_dm_kl1[tid][v_cnt] = dm_l1;
								op_dm_kr1[tid][v_cnt] = dm_r1;
								op_dm_kl2[tid][v_cnt] = dm_l2;
								op_dm_kr2[tid][v_cnt] = dm_r2;
							}

#endif
						}
						else
						{
							op_vector_pos1[tid][v_cnt] = mat_pos2[tid][r_i];
							op_vector_pos2[tid][v_cnt] = mat_pos1[tid][r_i];
#ifdef	MAPPING_QUALITY
							memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
							memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#else
							if(!((ld2 == 0) && (rd2 == 0)))
								memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
							if(!((ld1 == 0) && (rd1 == 0)))
								memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#endif
							op_dm_l1[tid][v_cnt] = ld2;
							op_dm_r1[tid][v_cnt] = rd2;
							op_dm_l2[tid][v_cnt] = ld1;
							op_dm_r2[tid][v_cnt] = rd1;

							op_dm_ex1[tid][v_cnt] = dm2;
							op_dm_ex2[tid][v_cnt] = dm1;

#ifdef	KSW_ALN_PAIR
							if(local_ksw)
							{
								op_dm_kl1[tid][v_cnt] = dm_l2;
								op_dm_kr1[tid][v_cnt] = dm_r2;
								op_dm_kl2[tid][v_cnt] = dm_l1;
								op_dm_kr2[tid][v_cnt] = dm_r1;
							}
#endif
						}

#ifdef	ALTER_DEBUG
						seed_length_arr[tid][v_cnt].seed_length = seed_length1 + seed_length2;
						seed_length_arr[tid][v_cnt].index = v_cnt;
#endif

						op_rc[tid][v_cnt] = mat_rc[tid][r_i];
						++v_cnt;
						dm_op[tid] = dm12;
					}
					else if(dm12 == dm_op[tid])
					{
						if(v_cnt < cus_max_output_ali)
						{
							if(mat_rc[tid][r_i] == 0)
							{
								op_vector_pos1[tid][v_cnt] = mat_pos1[tid][r_i];
								op_vector_pos2[tid][v_cnt] = mat_pos2[tid][r_i];
#ifdef	MAPPING_QUALITY
								memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
								memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#else
								if(!((ld1 == 0) && (rd1 == 0)))
									memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
								if(!((ld2 == 0) && (rd2 == 0)))
									memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#endif

								op_dm_l1[tid][v_cnt] = ld1;
								op_dm_r1[tid][v_cnt] = rd1;
								op_dm_l2[tid][v_cnt] = ld2;
								op_dm_r2[tid][v_cnt] = rd2;

								op_dm_ex1[tid][v_cnt] = dm1;
								op_dm_ex2[tid][v_cnt] = dm2;

#ifdef	KSW_ALN_PAIR
								if(local_ksw)
								{
									op_dm_kl1[tid][v_cnt] = dm_l1;
									op_dm_kr1[tid][v_cnt] = dm_r1;
									op_dm_kl2[tid][v_cnt] = dm_l2;
									op_dm_kr2[tid][v_cnt] = dm_r2;
								}
#endif
							}
							else
							{
								op_vector_pos1[tid][v_cnt] = mat_pos2[tid][r_i];
								op_vector_pos2[tid][v_cnt] = mat_pos1[tid][r_i];
#ifdef	MAPPING_QUALITY
								memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
								memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#else
								if(!((ld2 == 0) && (rd2 == 0)))
									memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
								if(!((ld1 == 0) && (rd1 == 0)))
									memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#endif
								op_dm_l1[tid][v_cnt] = ld2;
								op_dm_r1[tid][v_cnt] = rd2;
								op_dm_l2[tid][v_cnt] = ld1;
								op_dm_r2[tid][v_cnt] = rd1;

								op_dm_ex1[tid][v_cnt] = dm2;
								op_dm_ex2[tid][v_cnt] = dm1;

#ifdef	KSW_ALN_PAIR
								if(local_ksw)
								{
									op_dm_kl1[tid][v_cnt] = dm_l2;
									op_dm_kr1[tid][v_cnt] = dm_r2;
									op_dm_kl2[tid][v_cnt] = dm_l1;
									op_dm_kr2[tid][v_cnt] = dm_r1;
								}
#endif
							}
#ifdef	ALTER_DEBUG
							seed_length_arr[tid][v_cnt].seed_length = seed_length1 + seed_length2;
							seed_length_arr[tid][v_cnt].index = v_cnt;
#endif
							op_rc[tid][v_cnt] = mat_rc[tid][r_i];
							++v_cnt;
						}
					}
					else if(dm12 < dm_ops[tid])
					{
						vs_cnt = 0;
						if(mat_rc[tid][r_i] == 0)
						{
							ops_vector_pos1[tid][vs_cnt] = mat_pos1[tid][r_i];
							ops_vector_pos2[tid][vs_cnt] = mat_pos2[tid][r_i];
#ifdef	MAPPING_QUALITY
							memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
							memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#else
							if(!((ld1 == 0) && (rd1 == 0)))
								memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
							if(!((ld2 == 0) && (rd2 == 0)))
								memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#endif
							ops_dm_l1[tid][vs_cnt] = ld1;
							ops_dm_r1[tid][vs_cnt] = rd1;
							ops_dm_l2[tid][vs_cnt] = ld2;
							ops_dm_r2[tid][vs_cnt] = rd2;

							ops_dm_ex1[tid][vs_cnt] = dm1;
							ops_dm_ex2[tid][vs_cnt] = dm2;

#ifdef	KSW_ALN_PAIR
							if(local_ksw)
							{
								ops_dm_kl1[tid][vs_cnt] = dm_l1;
								ops_dm_kr1[tid][vs_cnt] = dm_r1;
								ops_dm_kl2[tid][vs_cnt] = dm_l2;
								ops_dm_kr2[tid][vs_cnt] = dm_r2;
							}
#endif
						}
						else
						{
							ops_vector_pos1[tid][vs_cnt] = mat_pos2[tid][r_i];
							ops_vector_pos2[tid][vs_cnt] = mat_pos1[tid][r_i];
#ifdef	MAPPING_QUALITY
							memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
							memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#else
							if(!((ld2 == 0) && (rd2 == 0)))
								memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
							if(!((ld1 == 0) && (rd1 == 0)))
								memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#endif
							ops_dm_l1[tid][vs_cnt] = ld2;
							ops_dm_r1[tid][vs_cnt] = rd2;
							ops_dm_l2[tid][vs_cnt] = ld1;
							ops_dm_r2[tid][vs_cnt] = rd1;

							ops_dm_ex1[tid][vs_cnt] = dm2;
							ops_dm_ex2[tid][vs_cnt] = dm1;

							//for S
#ifdef QUALS_CHECK
							ops_err[tid][vs_cnt] = error_f;
#endif

#ifdef	KSW_ALN_PAIR
							if(local_ksw)
							{
								ops_dm_kl1[tid][vs_cnt] = dm_l2;
								ops_dm_kr1[tid][vs_cnt] = dm_r2;
								ops_dm_kl2[tid][vs_cnt] = dm_l1;
								ops_dm_kr2[tid][vs_cnt] = dm_r1;
							}
#endif
						}
						ops_rc[tid][vs_cnt] = mat_rc[tid][r_i];
						++vs_cnt;
						dm_ops[tid] = dm12;
					}
					else if(dm12 == dm_ops[tid])
					{
						if(vs_cnt < cus_max_output_ali)
						{
							if(mat_rc[tid][r_i] == 0)
							{
								ops_vector_pos1[tid][vs_cnt] = mat_pos1[tid][r_i];
								ops_vector_pos2[tid][vs_cnt] = mat_pos2[tid][r_i];
#ifdef	MAPPING_QUALITY
								memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
								memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#else
								if(!((ld1 == 0) && (rd1 == 0)))
									memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
								if(!((ld2 == 0) && (rd2 == 0)))
									memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#endif
								ops_dm_l1[tid][vs_cnt] = ld1;
								ops_dm_r1[tid][vs_cnt] = rd1;
								ops_dm_l2[tid][vs_cnt] = ld2;
								ops_dm_r2[tid][vs_cnt] = rd2;

								ops_dm_ex1[tid][vs_cnt] = dm1;
								ops_dm_ex2[tid][vs_cnt] = dm2;

#ifdef	KSW_ALN_PAIR
								if(local_ksw)
								{
									ops_dm_kl1[tid][vs_cnt] = dm_l1;
									ops_dm_kr1[tid][vs_cnt] = dm_r1;
									ops_dm_kl2[tid][vs_cnt] = dm_l2;
									ops_dm_kr2[tid][vs_cnt] = dm_r2;
								}
#endif
							}
							else
							{
								ops_vector_pos1[tid][vs_cnt] = mat_pos2[tid][r_i];
								ops_vector_pos2[tid][vs_cnt] = mat_pos1[tid][r_i];
#ifdef	MAPPING_QUALITY
								memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
								memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#else
								if(!((ld2 == 0) && (rd2 == 0)))
									memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
								if(!((ld1 == 0) && (rd1 == 0)))
									memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#endif
								ops_dm_l1[tid][vs_cnt] = ld2;
								ops_dm_r1[tid][vs_cnt] = rd2;
								ops_dm_l2[tid][vs_cnt] = ld1;
								ops_dm_r2[tid][vs_cnt] = rd1;

								ops_dm_ex1[tid][vs_cnt] = dm2;
								ops_dm_ex2[tid][vs_cnt] = dm1;

#ifdef	KSW_ALN_PAIR
								if(local_ksw)
								{
									ops_dm_kl1[tid][vs_cnt] = dm_l2;
									ops_dm_kr1[tid][vs_cnt] = dm_r2;
									ops_dm_kl2[tid][vs_cnt] = dm_l1;
									ops_dm_kr2[tid][vs_cnt] = dm_r1;
								}
#endif
							}
							ops_rc[tid][vs_cnt] = mat_rc[tid][r_i];
							++vs_cnt;
						}
					}
				}

				if(local_ksw)
				{
					if((dm_cir_min > max_pair_score) && (cir_n != cir_fix_n))
					{
						seed_l[tid] -= seed_step;
						cir_n = cir_fix_n;
						continue;
					}
				}
				else
				{
					if(dm_cir_min > max_pair_score)
					{
						seed_l[tid] -= seed_step;
						cir_n++;
#ifdef	CIR_JUMP
						if(last_circle_rate)
						{
							lv_k1 = (read_length1 * last_circle_rate);
							lv_k2 = (read_length2 * last_circle_rate);
							max_pair_score = (read_length1 + read_length2) * last_circle_rate;
						}
#endif
						continue;
					}
				}
#ifdef ALI_OUT
				seqio[seqi].v_cnt = v_cnt;

				if(v_cnt > 0)
				{

#ifdef	ALTER_DEBUG
					if(v_cnt > 1)	qsort(seed_length_arr[tid], v_cnt, sizeof(seed_length_array), compare_seed_length);
#endif
					de_m_p_o[tid] = 0;
#ifdef	PR_SINGLE
					pr_o_f[tid] = 0;
#endif
					cnt.v_cnt = v_cnt;
					cnt.vs_cnt = vs_cnt;

					pair_sam_output(tid, read_length1, read_length2, f_cigarn, &cnt, seqi, lv_k1, lv_k2, pound_pos1_f_forward, pound_pos1_f_reverse, pound_pos1_r_forward, pound_pos1_r_reverse, pound_pos2_f_forward, pound_pos2_f_reverse, pound_pos2_r_forward, pound_pos2_r_reverse, cir_n);
				}
#endif
				break;
			}
			cir_n++;
		}

#ifdef	PAIR_RANDOM

#ifdef	PR_COV_FILTER
		if((((pos_ren[0][0][tid] > 0) && (pos_ren[0][1][tid] > 0)) || ((pos_ren[1][0][tid] > 0) && (pos_ren[1][1][tid] > 0))) && (de_m_p_o[tid] == 1) && (cov_filt_f[tid] == 0))
#else
		if((((pos_ren[0][0][tid] > 0) && (pos_ren[0][1][tid] > 0)) || ((pos_ren[1][0][tid] > 0) && (pos_ren[1][1][tid] > 0))) && (de_m_p_o[tid] == 1))
#endif
		{
#ifdef	PR_SINGLE
			seed_re_r[tid] = 1;
#endif

#ifdef	ALTER_DEBUG
			rep_go[tid] = 0;
#endif
			seed_repetitive(tid, read_length1, read_length2, f_cigarn, ref_copy_num1, ref_copy_num2, ref_copy_num_chars1, ref_copy_num_chars2, &cnt, seqi, lv_k1, lv_k2, pound_pos1_f_forward, pound_pos1_f_reverse, pound_pos1_r_forward, pound_pos1_r_reverse, pound_pos2_f_forward, pound_pos2_f_reverse, pound_pos2_r_forward, pound_pos2_r_reverse);
		}

#endif

		//deal with unmatched reads
#ifdef	UNMATCH_SINGLE_END

		if(de_m_p_o[tid] == 1)
		{
#ifdef	SINGLE_END_NOEXECUTE

			dm_op[tid] = MAX_OP_SCORE;
			dm_ops[tid] = MAX_OP_SCORE;
			v_cnt = 0;
			vs_cnt = 0;

			lv_k1 = (read_length1 * lv_rate_anchor) + 1;//lv_rate
			lv_k2 = (read_length2 * lv_rate_anchor) + 1;//lv_rate

#ifdef	SEED_FILTER_LENGTH
			max_sets_n[0][0] = (max_read_length[0][0] > max_read_length[1][1]) ? max_read_length[0][0]:max_read_length[1][1];
			max_sets_n[0][1] = (max_read_length[0][1] > max_read_length[1][0]) ? max_read_length[0][1]:max_read_length[1][0];

			max_sets_n[1][1] = max_sets_n[0][0];
			max_sets_n[1][0] = max_sets_n[0][1];
#endif

			for(rc_i = 0; rc_i < 2; rc_i++)
			{
				for(un_ii = 0; un_ii < 2; un_ii++)
				{
					if(unmatch[rc_i][un_ii] == 1)	continue;

					seed_pr1 = seed_pr[rc_i][un_ii];

#ifdef	SEED_FILTER_LENGTH
					for(off_i = 0; off_i < spa_i[rc_i][un_ii]; off_i++)
						if(seed_pr1[off_i].length < max_sets_n[rc_i][un_ii] - length_reduce)
						{
							++off_i;
							break;
						}
					spa_i[rc_i][un_ii] = off_i;
#endif

#ifdef	SEED_FILTER_POS
					qsort(seed_pr1, spa_i[rc_i][un_ii], sizeof(seed_pa), compare_seed_filter_posn);
					seed_posn_filter = 0;
#endif

#ifdef	SEED_FILTER
					single_lv_re = 0;
#endif
					if(rc_i == 0)
					{
						if(un_ii == 0)
						{
							read_bit_1[tid] = read_bit1[tid][0];
#ifdef	QUAL_FILT_SINGLE
							qual_filt_1 = qual_filt1[tid][0];
#endif

#ifdef	QUAL_FILT_LV_SINGLE
							qual_filt_lv_1 = qual_filt_lv1[tid][0];
							qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif
							ref_copy_num_chars = ref_copy_num_chars1;
							ref_copy_num = ref_copy_num1;
							low_mask = low_mask1[tid];
							sub_mask = sub_mask1[tid];
							ex_d_mask = ex_d_mask1[tid];
							max_mismatch = max_mismatch1_single[tid];
							lv_k = lv_k1;
							read_length = read_length1;

							pound_pos_1_f = pound_pos1_f_forward;
							pound_pos_1_r = pound_pos1_r_forward;
						}
						else
						{
							read_bit_1[tid] = read_bit2[tid][1];
#ifdef	QUAL_FILT_SINGLE
							qual_filt_1 = qual_filt2[tid][1];
#endif

#ifdef	QUAL_FILT_LV_SINGLE
							qual_filt_lv_1 = qual_filt_lv2[tid][1];
							qual_filt_lv_1_o = qual_filt_lv2[tid][0];
#endif
							ref_copy_num_chars = ref_copy_num_chars2;
							ref_copy_num = ref_copy_num2;
							low_mask = low_mask2[tid];
							sub_mask = sub_mask2[tid];
							ex_d_mask = ex_d_mask2[tid];
							max_mismatch = max_mismatch2_single[tid];
							lv_k = lv_k2;
							read_length = read_length2;

							pound_pos_1_f = pound_pos2_f_reverse;
							pound_pos_1_r = pound_pos2_r_reverse;

						}
					}
					else
					{
						if(un_ii == 0)
						{
							read_bit_1[tid] = read_bit2[tid][0];
#ifdef	QUAL_FILT_SINGLE
							qual_filt_1 = qual_filt2[tid][0];
#endif

#ifdef	QUAL_FILT_LV_SINGLE
							qual_filt_lv_1 = qual_filt_lv2[tid][0];
							qual_filt_lv_1_o = qual_filt_lv2[tid][1];
#endif
							ref_copy_num_chars = ref_copy_num_chars2;
							ref_copy_num = ref_copy_num2;
							low_mask = low_mask2[tid];
							sub_mask = sub_mask2[tid];
							ex_d_mask = ex_d_mask2[tid];
							max_mismatch = max_mismatch2_single[tid];
							lv_k = lv_k2;
							read_length = read_length2;

							pound_pos_1_f = pound_pos2_f_forward;
							pound_pos_1_r = pound_pos2_r_forward;

						}
						else
						{
							read_bit_1[tid] = read_bit1[tid][1];
#ifdef	QUAL_FILT_SINGLE
							qual_filt_1 = qual_filt1[tid][1];
#endif

#ifdef	QUAL_FILT_LV_SINGLE
							qual_filt_lv_1 = qual_filt_lv1[tid][1];
							qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif
							ref_copy_num_chars = ref_copy_num_chars1;
							ref_copy_num = ref_copy_num1;
							low_mask = low_mask1[tid];
							sub_mask = sub_mask1[tid];
							ex_d_mask = ex_d_mask1[tid];
							max_mismatch = max_mismatch1_single[tid];
							lv_k = lv_k1;
							read_length = read_length1;

							pound_pos_1_f = pound_pos1_f_reverse;
							pound_pos_1_r = pound_pos1_r_reverse;
						}
					}

					for(seed1_i = 0; seed1_i < spa_i[rc_i][un_ii]; seed1_i++)
					{
#ifdef	SEED_FILTER
						if(seed_pr1[seed1_i].length > (read_length >> 1))	single_lv_re = 1;
						else
						{
							if(single_lv_re == 1)	continue;
						}
#endif
						if(seed_pr1[seed1_i].ui == 1)
						{
							end1_uc_f = 0;
							nuc1_f = 0;

							d_l1 = seed_pr1[seed1_i].ref_pos_off;
							d_r1 = seed_pr1[seed1_i].ref_pos_off_r;
						}
						else
						{
							end1_uc_f = 1;
						}

						q_rear1 = 0;
						q_n1 = 0;

						dmt1 = ali_exl;
						lv_dmt1 = lv_k;

						s_r_o_l1 = seed_pr1[seed1_i].s_r_o_l;
						s_r_o_r1 = seed_pr1[seed1_i].s_r_o_r;

#ifdef	ALTER_DEBUG_ANCHOR
						seed_length1 = seed_pr1[seed1_i].length;
#endif
						for(psp_i = 0; psp_i < seed_pr1[seed1_i].pos_n; psp_i++)
						{

#ifdef	SEED_FILTER_POS
							if(seed_posn_filter > seed_filter_pos_num)	break;

							seed_posn_filter++;
#endif
							posi = seed_set_pos[rc_i][un_ii][tid][seed_pr1[seed1_i].pos_start + psp_i];
#ifdef ANCHOR_LV_S
							s_offset_l = 0;
							s_offset_r = 0;
#endif
							//for end1
							if((end1_uc_f == 0) && ((dmt1 <= d_l1) && (dmt1 <= d_r1)))   // && (nuc1_f == 0)
							{
								if(nuc1_f == 0)
								{
									pos_l = posi - max_extension_length - 1;
									re_d = pos_l & 0X1f;
									b_t_n_r = 32 - re_d;

									if(re_d != 0)
									{
										tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);
										memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars);

										for(rst_i = 0; rst_i < ref_copy_num - 1; rst_i++)
										{
											tran_tmp = (ref_seq_tmp1[tid][rst_i] & bit_tran_re[re_d]);

											ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
											ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
											tran_tmp_p = tran_tmp;
										}
										ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
										ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

										//clear the lowest n bit
										ref_seq_tmp1[tid][rst_i] &= low_mask;

									}
									else
									{
										memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars);

										//clear the lowest n bit
										ref_seq_tmp1[tid][ref_copy_num - 1] &= low_mask;
									}

									//exact match
									mis_c_n = 0;
#ifdef	QUAL_FILT_SINGLE
									mis_c_n_filt = 0;
#endif
									ref_tmp_ori = ref_seq_tmp1[tid][ref_copy_num - 2];
									ref_seq_tmp1[tid][ref_copy_num - 2] &= low_mask;

									for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num - 1; rst_i++, read_b_i++)
									{
										xor_tmp = ref_seq_tmp1[tid][rst_i] ^ read_bit_1[tid][read_b_i];
#ifdef	QUAL_FILT_SINGLE
										mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
										xor_tmp &= qual_filt_1[read_b_i];
										mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
										if(mis_c_n_filt > max_mismatch)	break;
#else
										mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
										mis_c_n_filt = mis_c_n;
										if(mis_c_n > max_mismatch)	break;
#endif
									}

									ref_seq_tmp1[tid][ref_copy_num - 2] = ref_tmp_ori;

									//dm = mis_c_n;
									//lv
									if(mis_c_n_filt > max_mismatch)
									{
#ifdef SPLIT_LV
										pound_mis = 0;
										if(pound_pos_1_f >= s_r_o_r1)   //1
										{
#ifdef	POUND_MIS
											for(bit_char_i = pound_pos_1_f, bit_char_i_ref = pound_pos_1_f + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(pound_pos_1_f < pound_pos_1_r)
											{
												read_pos_start_num = pound_pos_1_f >> 5;
												read_pos_end_num = (pound_pos_1_r - 1) >> 5;
												read_pos_re = (pound_pos_1_f & 0X1f) << 1;
												read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

												if(read_pos_start_num == read_pos_end_num)
												{
													xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = 0;
											lv_up_right = s_r_o_l1;
											lv_down_right = pound_pos_1_f;
											lv_down_left = s_r_o_r1;
										}
										else if(pound_pos_1_r <= s_r_o_l1 + 1)     //5
										{
#ifdef	POUND_MIS
											for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(pound_pos_1_r > 0)
											{
												read_pos_end_num = (pound_pos_1_r - 1) >> 5;
												read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

												if(read_pos_end_num == 0)
												{
													xor_tmp = (read_bit_1[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = read_bit_1[tid][0] ^ ref_seq_tmp1[tid][1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = pound_pos_1_r;
											lv_up_right = s_r_o_l1;
											lv_down_right = read_length;
											lv_down_left = s_r_o_r1;
										}
										else if((pound_pos_1_f <= s_r_o_l1 + 1) && (pound_pos_1_r >= s_r_o_r1))     //2
										{
#ifdef	POUND_MIS
											for(bit_char_i = pound_pos_1_f, bit_char_i_ref = pound_pos_1_f + 32; (bit_char_i <= s_r_o_l1) && (bit_char_i_ref <= s_r_o_l1 + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;


											for(bit_char_i = s_r_o_r1, bit_char_i_ref = s_r_o_r1 + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(pound_pos_1_f <= s_r_o_l1)
											{
												read_pos_start_num = pound_pos_1_f >> 5;
												read_pos_end_num = s_r_o_l1 >> 5;
												read_pos_re = (pound_pos_1_f & 0X1f) << 1;
												read_pos_end_re = (s_r_o_l1 & 0X1f) << 1;

												if(read_pos_start_num == read_pos_end_num)
												{
													xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}

											if(s_r_o_r1 < pound_pos_1_r)
											{
												read_pos_start_num = s_r_o_r1 >> 5;
												read_pos_end_num = (pound_pos_1_r - 1) >> 5;
												read_pos_re = (s_r_o_r1 & 0X1f) << 1;
												read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

												if(read_pos_start_num == read_pos_end_num)
												{
													xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = 0;
											lv_up_right = pound_pos_1_f - 1;
											lv_down_right = read_length;
											lv_down_left = pound_pos_1_r;
										}
										else if((pound_pos_1_f > s_r_o_l1 + 1) && (pound_pos_1_f < s_r_o_r1))     //3
										{
#ifdef	POUND_MIS
											for(bit_char_i = s_r_o_r1, bit_char_i_ref = s_r_o_r1 + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(s_r_o_r1 < pound_pos_1_r)
											{
												read_pos_start_num = s_r_o_r1 >> 5;
												read_pos_end_num = (pound_pos_1_r - 1) >> 5;
												read_pos_re = (s_r_o_r1 & 0X1f) << 1;
												read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

												if(read_pos_start_num == read_pos_end_num)
												{
													xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = 0;
											lv_up_right = s_r_o_l1;
											lv_down_right = read_length;
											lv_down_left = pound_pos_1_r;
										}
										else     //4
										{
#ifdef	POUND_MIS
											for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < s_r_o_l1) && (bit_char_i_ref < s_r_o_l1 + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(s_r_o_l1 > 0)
											{
												read_pos_end_num = (s_r_o_l1 - 1) >> 5;
												read_pos_end_re = ((s_r_o_l1 - 1) & 0X1f) << 1;

												if(read_pos_end_num == 0)
												{
													xor_tmp = (read_bit_1[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = read_bit_1[tid][0] ^ ref_seq_tmp1[tid][1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = 0;
											lv_up_right = -1;
											lv_down_right = read_length;
#ifdef	POUND_MODIFY
											lv_down_left = s_r_o_r1;
#else
											lv_down_left = s_r_o_l1;
#endif
										}
#ifdef	QUAL_FILT_LV_SINGLE
										for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#ifdef	ANCHOR_LV_S
										dm_l1 = computeEditDistance_mis_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid], qual_filt_lv_1_o + read_length - 1 - lv_up_right, &s_offset_l);
#else
										dm_l1 = computeEditDistance_mis(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid], qual_filt_lv_1_o + read_length - 1 - lv_up_right);
#endif
										for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#ifdef	ANCHOR_LV_S
										dm_r1 = computeEditDistance_mis_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid], qual_filt_lv_1 + lv_down_left, &s_offset_r);
#else
										dm_r1 = computeEditDistance_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid], qual_filt_lv_1 + lv_down_left);
#endif
										dm1 = dm_l1 + dm_r1 + pound_mis;

#else
										for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_l1 = computeEditDistance(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid]);

										for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_r1 = computeEditDistance(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid]);

										dm1 = dm_l1 + dm_r1 + pound_mis;

#endif

										ld1 = s_r_o_l1;
										rd1 = s_r_o_r1;
#endif
									}
									else
									{
										dm1 = mis_c_n;
										ld1 = 0;
										rd1 = 0;
										dm_l1 = 0;
										dm_r1 = 0;
									}

									//these two values could be different
									if((dm_l1 != -1) && (dm_r1 != -1))
									{
										if(dm1 < dmt1)	dmt1 = dm1 + 1;

										if(dm1 < lv_dmt1)	lv_dmt1 = dm1 + 1;

										if(dm1 < max_mismatch - 1)	dmt1 = 0;
									}
									else
									{
#ifdef	ANCHOR_LV_S
										if((s_offset_l + s_offset_r) < ((float )read_length * ANCHOR_LV_S_FLOAT))
											dm1 = lv_dmt1 + pound_mis;
										else	dm1 = MAX_EDIT_SCORE;
#else
										dm1 = MAX_EDIT_SCORE;
#endif
									}
									nuc1_f = 1;
								}
							}
							else
							{
								pos_l = posi - max_extension_length - 1;
								re_d = pos_l & 0X1f;
								b_t_n_r = 32 - re_d;

								if(re_d != 0)
								{
									tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);
									memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars);

									for(rst_i = 0; rst_i < ref_copy_num - 1; rst_i++)
									{
										tran_tmp = (ref_seq_tmp1[tid][rst_i] & bit_tran_re[re_d]);

										ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
										ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
										tran_tmp_p = tran_tmp;
									}
									ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
									ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

									//clear the lowest n bit
									ref_seq_tmp1[tid][rst_i] &= low_mask;

								}
								else
								{
									memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars);

									//clear the lowest n bit
									ref_seq_tmp1[tid][ref_copy_num - 1] &= low_mask;
								}

								//trim the beginning and end of the current ref seq based on current minimum edit distance dm_t
								s_m_t = sub_mask[dmt1];
								ref_seq_tmp1[tid][0] &= bit_tran[dmt1];
								ref_seq_tmp1[tid][s_m_t] &= ex_d_mask[dmt1];

								//traverse and check whether there is an existing seq that is as same as current new ref seq
								c_m_f = 0;
								for(q_rear_i = q_rear1 - 1; q_rear_i >= 0; q_rear_i--)
								{
									ref_tmp_ori = cache_end1[q_rear_i][0];
									cache_end1[q_rear_i][0] &= bit_tran[dmt1];

									ref_tmp_ori2 = cache_end1[q_rear_i][s_m_t];
									cache_end1[q_rear_i][s_m_t] &= ex_d_mask[dmt1];

									cmp_re = memcmp(cache_end1[q_rear_i], ref_seq_tmp1[tid], (s_m_t + 1) << 3);

									cache_end1[q_rear_i][0] = ref_tmp_ori;
									cache_end1[q_rear_i][s_m_t] = ref_tmp_ori2;

									if(cmp_re == 0)
									{
										//deal with finding an alignment
										dm1 = cache_dis1[q_rear_i];
										ld1 = cache_dml1[q_rear_i];
										rd1 = cache_dmr1[q_rear_i];

										c_m_f = 1;
										break;
									}
								}

								if((q_n1 > MAX_Q_NUM) && (q_rear_i < 0))
								{
									for(q_rear_i = MAX_Q_NUM - 1; q_rear_i >= q_rear1; q_rear_i--)
									{
										ref_tmp_ori = cache_end1[q_rear_i][0];
										cache_end1[q_rear_i][0] &= bit_tran[dmt1];

										ref_tmp_ori2 = cache_end1[q_rear_i][s_m_t];
										cache_end1[q_rear_i][s_m_t] &= ex_d_mask[dmt1];

										cmp_re = memcmp(cache_end1[q_rear_i], ref_seq_tmp1[tid], (s_m_t + 1) << 3);

										cache_end1[q_rear_i][0] = ref_tmp_ori;
										cache_end1[q_rear_i][s_m_t] = ref_tmp_ori2;

										if(cmp_re == 0)
										{
											//deal with finding an alignment
											dm1 = cache_dis1[q_rear_i];
											ld1 = cache_dml1[q_rear_i];
											rd1 = cache_dmr1[q_rear_i];

											c_m_f = 1;
											break;
										}
									}
								}
								//do not find the seq in cache, exact match or lv and add into cache
								if(c_m_f == 0)
								{
									//exact match
									mis_c_n = 0;
#ifdef	QUAL_FILT_SINGLE
									mis_c_n_filt = 0;
#endif
									ref_tmp_ori = ref_seq_tmp1[tid][ref_copy_num - 2];
									ref_seq_tmp1[tid][ref_copy_num - 2] &= low_mask;

									for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num - 1; rst_i++, read_b_i++)
									{
										xor_tmp = ref_seq_tmp1[tid][rst_i] ^ read_bit_1[tid][read_b_i];
#ifdef	QUAL_FILT_SINGLE
										mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
										xor_tmp &= qual_filt_1[read_b_i];
										mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
										if(mis_c_n_filt > max_mismatch)	break;
#else
										mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
										mis_c_n_filt = mis_c_n;
										if(mis_c_n > max_mismatch)	break;
#endif
									}

									ref_seq_tmp1[tid][ref_copy_num - 2] = ref_tmp_ori;

									//lv
									if(mis_c_n_filt > max_mismatch)
									{
#ifdef SPLIT_LV
										pound_mis = 0;
										if(pound_pos_1_f >= s_r_o_r1)   //1
										{
#ifdef	POUND_MIS
											for(bit_char_i = pound_pos_1_f, bit_char_i_ref = pound_pos_1_f + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(pound_pos_1_f < pound_pos_1_r)
											{
												read_pos_start_num = pound_pos_1_f >> 5;
												read_pos_end_num = (pound_pos_1_r - 1) >> 5;
												read_pos_re = (pound_pos_1_f & 0X1f) << 1;
												read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

												if(read_pos_start_num == read_pos_end_num)
												{
													xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = 0;
											lv_up_right = s_r_o_l1;
											lv_down_right = pound_pos_1_f;
											lv_down_left = s_r_o_r1;
										}
										else if(pound_pos_1_r <= s_r_o_l1 + 1)     //5
										{
#ifdef	POUND_MIS
											for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(pound_pos_1_r > 0)
											{
												read_pos_end_num = (pound_pos_1_r - 1) >> 5;
												read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

												if(read_pos_end_num == 0)
												{
													xor_tmp = (read_bit_1[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = read_bit_1[tid][0] ^ ref_seq_tmp1[tid][1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = pound_pos_1_r;
											lv_up_right = s_r_o_l1;
											lv_down_right = read_length;
											lv_down_left = s_r_o_r1;
										}
										else if((pound_pos_1_f <= s_r_o_l1 + 1) && (pound_pos_1_r >= s_r_o_r1))     //2
										{
#ifdef	POUND_MIS
											for(bit_char_i = pound_pos_1_f, bit_char_i_ref = pound_pos_1_f + 32; (bit_char_i <= s_r_o_l1) && (bit_char_i_ref <= s_r_o_l1 + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;


											for(bit_char_i = s_r_o_r1, bit_char_i_ref = s_r_o_r1 + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(pound_pos_1_f <= s_r_o_l1)
											{
												read_pos_start_num = pound_pos_1_f >> 5;
												read_pos_end_num = s_r_o_l1 >> 5;
												read_pos_re = (pound_pos_1_f & 0X1f) << 1;
												read_pos_end_re = (s_r_o_l1 & 0X1f) << 1;

												if(read_pos_start_num == read_pos_end_num)
												{
													xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}

											if(s_r_o_r1 < pound_pos_1_r)
											{
												read_pos_start_num = s_r_o_r1 >> 5;
												read_pos_end_num = (pound_pos_1_r - 1) >> 5;
												read_pos_re = (s_r_o_r1 & 0X1f) << 1;
												read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

												if(read_pos_start_num == read_pos_end_num)
												{
													xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = 0;
											lv_up_right = pound_pos_1_f - 1;
											lv_down_right = read_length;
											lv_down_left = pound_pos_1_r;
										}
										else if((pound_pos_1_f > s_r_o_l1 + 1) && (pound_pos_1_f < s_r_o_r1))     //3
										{
#ifdef	POUND_MIS
											for(bit_char_i = s_r_o_r1, bit_char_i_ref = s_r_o_r1 + 32; (bit_char_i < pound_pos_1_r) && (bit_char_i_ref < pound_pos_1_r + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(s_r_o_r1 < pound_pos_1_r)
											{
												read_pos_start_num = s_r_o_r1 >> 5;
												read_pos_end_num = (pound_pos_1_r - 1) >> 5;
												read_pos_re = (s_r_o_r1 & 0X1f) << 1;
												read_pos_end_re = ((pound_pos_1_r - 1) & 0X1f) << 1;

												if(read_pos_start_num == read_pos_end_num)
												{
													xor_tmp = ((read_bit_1[tid][read_pos_start_num] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re)) ^ ((ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re) >> (62 - read_pos_end_re + read_pos_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = (read_bit_1[tid][read_pos_start_num] << read_pos_re) ^ (ref_seq_tmp1[tid][read_pos_start_num + 1] << read_pos_re);
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = read_pos_start_num + 1, rst_i_1 = read_pos_start_num + 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = 0;
											lv_up_right = s_r_o_l1;
											lv_down_right = read_length;
											lv_down_left = pound_pos_1_r;
										}
										else     //4
										{
#ifdef	POUND_MIS
											for(bit_char_i = 0, bit_char_i_ref = 32; (bit_char_i < s_r_o_l1) && (bit_char_i_ref < s_r_o_l1 + 32); bit_char_i++, bit_char_i_ref++)
												if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ref_seq_tmp1[tid][bit_char_i_ref >> 5] >> ((31 - (bit_char_i_ref & 0X1f)) << 1)) & 0X3))
													pound_mis++;
#else

											if(s_r_o_l1 > 0)
											{
												read_pos_end_num = (s_r_o_l1 - 1) >> 5;
												read_pos_end_re = ((s_r_o_l1 - 1) & 0X1f) << 1;

												if(read_pos_end_num == 0)
												{
													xor_tmp = (read_bit_1[tid][0] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

												}
												else
												{
													xor_tmp = read_bit_1[tid][0] ^ ref_seq_tmp1[tid][1];
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));

													for(rst_i = 1, rst_i_1 = 2; rst_i < read_pos_end_num; rst_i++, rst_i_1++)
													{
														xor_tmp = read_bit_1[tid][rst_i] ^ ref_seq_tmp1[tid][rst_i_1];
														pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
													}

													xor_tmp = (read_bit_1[tid][read_pos_end_num] >> (62 - read_pos_end_re)) ^ (ref_seq_tmp1[tid][read_pos_end_num + 1] >> (62 - read_pos_end_re));
													pound_mis += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
												}
											}
#endif
											lv_up_left = 0;
											lv_up_right = -1;
											lv_down_right = read_length;
#ifdef	POUND_MODIFY
											lv_down_left = s_r_o_r1;
#else
											lv_down_left = s_r_o_l1;
#endif
										}
#ifdef	QUAL_FILT_LV_SINGLE
										for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

#ifdef	ANCHOR_LV_S
										dm_l1 = computeEditDistance_mis_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid], qual_filt_lv_1_o + read_length - 1 - lv_up_right, &s_offset_l);
#else
										dm_l1 = computeEditDistance_mis(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid], qual_filt_lv_1_o + read_length - 1 - lv_up_right);
#endif
										for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);
#ifdef	ANCHOR_LV_S
										dm_r1 = computeEditDistance_mis_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid], qual_filt_lv_1 + lv_down_left, &s_offset_r);
#else
										dm_r1 = computeEditDistance_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid], qual_filt_lv_1 + lv_down_left);
#endif
										dm1 = dm_l1 + dm_r1 + pound_mis;

#else
										for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

										dm_l1 = computeEditDistance(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_dmt1, L[tid]);

										for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

										dm_r1 = computeEditDistance(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_dmt1, L[tid]);

										dm1 = dm_l1 + dm_r1 + pound_mis;
#endif
										ld1 = s_r_o_l1;
										rd1 = s_r_o_r1;
#endif
									}
									else
									{
										dm1 = mis_c_n;
										ld1 = 0;
										rd1 = 0;
										dm_l1 = 0;
										dm_r1 = 0;
									}

									//these two values could be different
									if((dm_l1 != -1) && (dm_r1 != -1))
									{
										if(dm1 < dmt1)	dmt1 = dm1 + 1;

										if(dm1 < lv_dmt1)	lv_dmt1 = dm1 + 1;

										if(dm1 < max_mismatch - 1)	dmt1 = 0;
									}
									else
									{
#ifdef	ANCHOR_LV_S
										if((s_offset_l + s_offset_r) < ((float )read_length * ANCHOR_LV_S_FLOAT))
											dm1 = lv_dmt1 + pound_mis;
										else	dm1 = MAX_EDIT_SCORE;
#else
										dm1 = MAX_EDIT_SCORE;
#endif
									}

									//add the ref sequence at the end of queue
									memcpy(cache_end1[q_rear1], ref_seq_tmp1[tid], ref_copy_num_chars);
									cache_dis1[q_rear1] = dm1;
									cache_dml1[q_rear1] = ld1;
									cache_dmr1[q_rear1] = rd1;

									q_rear1 = ((q_rear1 + 1) & 0X1f);
									++q_n1;

									//add edit distance
								}
							}


							/*
								x = posi;
													low = 0;
													high = chr_file_n - 1;

													while ( low <= high )
													{
							mid = (low + high) >> 1;
							if(x < (chr_end_n[mid]))
							{
								high = mid - 1;
							}
							else if(x > (chr_end_n[mid]))
							{
								low = mid + 1;
							}
							else
							{
								chr_re =  mid;
								break;
							}
							chr_re = low;
													}

													posi = x - chr_end_n[chr_re - 1] + 1;

													fprintf(stderr, "to add: %u %u %d %d %d rc: %u\n", chr_re, posi, dm1, dm_op[tid], dm_ops[tid],((rc_i << 1) + un_ii));

													posi = x;
								*/		


							if(dm1 < dm_op[tid])
							{
#ifdef	DM_COPY_SINGLE
								for(dm_i = 0; dm_i < v_cnt; dm_i++)
								{
									ops_vector_pos1[tid][dm_i] = op_vector_pos1[tid][dm_i];

									ops_dm_l1[tid][dm_i] = op_dm_l1[tid][dm_i];
									ops_dm_r1[tid][dm_i] = op_dm_r1[tid][dm_i];

									if(!((op_dm_l1[tid][dm_i] == 0) && (op_dm_r1[tid][dm_i] == 0)))
										memcpy(ops_vector_seq1[tid][dm_i], op_vector_seq1[tid][dm_i], ref_copy_num_chars);

									ops_rc[tid][dm_i] = op_rc[tid][dm_i];
								}
								vs_cnt = v_cnt;
								dm_ops[tid] = dm_op[tid];
#endif

								v_cnt = 0;
								op_vector_pos1[tid][v_cnt] = posi;
								if(!((ld1 == 0) && (rd1 == 0)))
									memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);

								op_dm_l1[tid][v_cnt] = ld1;
								op_dm_r1[tid][v_cnt] = rd1;

#ifdef	ALTER_DEBUG_ANCHOR
								seed_length_arr[tid][v_cnt].seed_length = seed_length1;
								seed_length_arr[tid][v_cnt].index = v_cnt;
#endif
								op_rc[tid][v_cnt] = ((rc_i << 1) + un_ii);
								++v_cnt;
								dm_op[tid] = dm1;
							}
							else if(dm1 == dm_op[tid])
							{
								if(v_cnt < cus_max_output_ali)
								{
									op_vector_pos1[tid][v_cnt] = posi;
									if(!((ld1 == 0) && (rd1 == 0)))
										memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);

									op_dm_l1[tid][v_cnt] = ld1;
									op_dm_r1[tid][v_cnt] = rd1;

#ifdef	ALTER_DEBUG_ANCHOR
									seed_length_arr[tid][v_cnt].seed_length = seed_length1;
									seed_length_arr[tid][v_cnt].index = v_cnt;
#endif
									op_rc[tid][v_cnt] = ((rc_i << 1) + un_ii);
									++v_cnt;
								}
							}
							else if(dm1 < dm_ops[tid])
							{
								vs_cnt = 0;

								ops_vector_pos1[tid][vs_cnt] = posi;
								if(!((ld1 == 0) && (rd1 == 0)))
									memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);

								ops_dm_l1[tid][vs_cnt] = ld1;
								ops_dm_r1[tid][vs_cnt] = rd1;

								ops_rc[tid][vs_cnt] = ((rc_i << 1) + un_ii);

								++vs_cnt;
								dm_ops[tid] = dm1;
							}
							else if(dm1 == dm_ops[tid])
							{
								if(vs_cnt < cus_max_output_ali)
								{
									ops_vector_pos1[tid][vs_cnt] = posi;
									if(!((ld1 == 0) && (rd1 == 0)))
										memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);

									ops_dm_l1[tid][vs_cnt] = ld1;
									ops_dm_r1[tid][vs_cnt] = rd1;

									ops_rc[tid][vs_cnt] = ((rc_i << 1) + un_ii);

									++vs_cnt;
								}
							}
						}
					}
				}
			}

#ifdef	PR_SINGLE
			if(seed_re_r[tid] == 1)
			{
				if(min_mis[tid] <= dm_op[tid])
				{
					if(min_mis[tid] < dm_op[tid])
					{
						dm_op[tid] = min_mis[tid];
						v_cnt = 0;
					}

					for(rc_i = 0; rc_i < 2; rc_i++)
						for(rc_ii = 0; rc_ii < 2; rc_ii++)
							for(v_cnt_i = 0; v_cnt_i < seedpos_misn[rc_i][rc_ii][tid]; v_cnt_i++)
							{
								if(seedpos_mis[rc_i][rc_ii][tid][v_cnt_i] == min_mis[tid])
								{
									if(v_cnt < cus_max_output_ali)
									{
										op_vector_pos1[tid][v_cnt] = seedpos[rc_i][rc_ii][tid][v_cnt_i];
										op_dm_l1[tid][v_cnt] = 0;
										op_dm_r1[tid][v_cnt] = 0;
										op_rc[tid][v_cnt] = ((rc_i << 1) + rc_ii);
#ifdef	ALTER_DEBUG_ANCHOR
										seed_length_arr[tid][v_cnt].seed_length = 0;
										seed_length_arr[tid][v_cnt].index = v_cnt;
#endif
										++v_cnt;
									}
								}
							}
				}
				else if(min_mis[tid] <= dm_ops[tid])
				{
					if(min_mis[tid] < dm_ops[tid])
					{
						dm_ops[tid] = min_mis[tid];
						vs_cnt = 0;
					}

					for(rc_i = 0; rc_i < 2; rc_i++)
						for(rc_ii = 0; rc_ii < 2; rc_ii++)
							for(v_cnt_i = 0; v_cnt_i < seedpos_misn[rc_i][rc_ii][tid]; v_cnt_i++)
							{
								if(seedpos_mis[rc_i][rc_ii][tid][v_cnt_i] == min_mis[tid])
								{
									if(vs_cnt < cus_max_output_ali)
									{
										ops_vector_pos1[tid][vs_cnt] = seedpos[rc_i][rc_ii][tid][v_cnt_i];
										ops_dm_l1[tid][vs_cnt] = 0;
										ops_dm_r1[tid][vs_cnt] = 0;
										ops_rc[tid][vs_cnt] = ((rc_i << 1) + rc_ii);
										++vs_cnt;
									}
								}
							}
				}
			}
#endif

			if(v_cnt > 0)
			{
#ifdef	ALTER_DEBUG_ANCHOR
				if(v_cnt > 1)	qsort(seed_length_arr[tid], v_cnt, sizeof(seed_length_array), compare_seed_length);
#endif

#ifdef	REDUCE_ANCHOR

				tra1_i = 0;
				tra2_i = 0;
				anchor_n1 = 0;
				anchor_n2 = 0;

				memset(op_mask, 0, v_cnt << 1);
				memset(ops_mask, 0, vs_cnt << 1);

				for(va_cnt_i = 0; va_cnt_i < v_cnt; va_cnt_i++)
				{
#ifdef	ALTER_DEBUG_ANCHOR
					if(v_cnt > 1)	tra_i = seed_length_arr[tid][va_cnt_i].index;
					else	tra_i = va_cnt_i;
#else
					tra_i = va_cnt_i;
#endif
					if((op_rc[tid][tra_i] == 0) || (op_rc[tid][tra_i] == 3))
					{
						poses1[tid][anchor_n1] = op_vector_pos1[tid][tra_i];
						ls1[tid][anchor_n1] = op_dm_l1[tid][tra_i];
						rs1[tid][anchor_n1] = op_dm_r1[tid][tra_i];
						rcs1[tid][anchor_n1] = op_rc[tid][tra_i];
						dms1[tid][anchor_n1] = dm_op[tid];
						orders1[tid][anchor_n1] = tra_i;

						anchor_n1++;
					}
					else
					{
						poses2[tid][anchor_n2] = op_vector_pos1[tid][tra_i];
						ls2[tid][anchor_n2] = op_dm_l1[tid][tra_i];
						rs2[tid][anchor_n2] = op_dm_r1[tid][tra_i];
						rcs2[tid][anchor_n2] = op_rc[tid][tra_i];
						dms2[tid][anchor_n2] = dm_op[tid];
						orders2[tid][anchor_n2] = tra_i;

						anchor_n2++;
					}
				}

				for(tra_i = 0; tra_i < vs_cnt; tra_i++)
				{
					if((ops_rc[tid][tra_i] == 0) || (ops_rc[tid][tra_i] == 3))
					{
						poses1[tid][anchor_n1] = ops_vector_pos1[tid][tra_i];
						ls1[tid][anchor_n1] = ops_dm_l1[tid][tra_i];
						rs1[tid][anchor_n1] = ops_dm_r1[tid][tra_i];
						rcs1[tid][anchor_n1] = ops_rc[tid][tra_i];
						dms1[tid][anchor_n1] = dm_ops[tid];
						orders1[tid][anchor_n1] = tra_i + MAX_REDUCE_ANCHOR_NUM;

						anchor_n1++;
					}
					else
					{
						poses2[tid][anchor_n2] = ops_vector_pos1[tid][tra_i];
						ls2[tid][anchor_n2] = ops_dm_l1[tid][tra_i];
						rs2[tid][anchor_n2] = ops_dm_r1[tid][tra_i];
						rcs2[tid][anchor_n2] = ops_rc[tid][tra_i];
						dms2[tid][anchor_n2] = dm_ops[tid];
						orders2[tid][anchor_n2] = tra_i + MAX_REDUCE_ANCHOR_NUM;

						anchor_n2++;
					}
				}
#endif

#ifdef	ANCHOR_HASH_ALI
				for(anchor_hash_i = 0; anchor_hash_i < 4; anchor_hash_i++)
				{
					if(anchor_hash_i == 0)
					{
						read_length = read_length2;
						ori = read_bit2[tid][1];
					}
					else if(anchor_hash_i == 1)
					{
						read_length = read_length1;
						ori = read_bit1[tid][0];
					}
					else if(anchor_hash_i == 2)
					{
						read_length = read_length1;
						ori = read_bit1[tid][1];
					}
					else
					{
						read_length = read_length2;
						ori = read_bit2[tid][0];
					}

					des_i = 0;
					for(ori_i = 0; ori_i < read_length >> 5; ori_i++)
					{
						for(char_i = 0; (char_i <= 64 - k_anchor_b) && (des_i <= read_length - k_anchor); char_i += 2, des_i++)
						{
							rh[tid][des_i].des = ((ori[ori_i] >> (64 - char_i - k_anchor_b)) & anchor_mask);
							rh[tid][des_i].off_set = des_i;
						}
						//note that boundary value
						for(char_i = 2; (char_i < k_anchor_b) && (des_i <= read_length - k_anchor); char_i += 2, des_i++)
						{
							rh[tid][des_i].des = ((ori[ori_i] & anchor_mask_boundary_re[char_i >> 1]) << char_i) | (ori[ori_i + 1] >> (64 - char_i));// & anchor_mask_boundary[char_i >> 1]
							rh[tid][des_i].off_set = des_i;
						}
					}

					qsort(rh[tid], des_i, sizeof(read_h), compare_read_hash);

					memset(anchor_hash[tid][anchor_hash_i], 0, 1 << 17);

					read_k_p = rh[tid][0].des;
					anchor_hash[tid][anchor_hash_i][read_k_p >> k_anchor_back] = 0;
					anchor_array[tid][anchor_hash_i][0] = (read_k_p & anchor_back_mask);
					anchor_point[tid][anchor_hash_i][0] = 0;
					anchor_pos[tid][anchor_hash_i][0] = rh[tid][0].off_set;
					array_i = 1;
					array_i_p = 0;

					for(char_i = 1; char_i < des_i; char_i++)
					{
						anchor_pos[tid][anchor_hash_i][char_i] = rh[tid][char_i].off_set;

						read_k_t = rh[tid][char_i].des;

						if((read_k_t >> k_anchor_back) == (read_k_p >> k_anchor_back))
						{
							if((read_k_t & anchor_back_mask) != (read_k_p & anchor_back_mask))
							{
								anchor_array[tid][anchor_hash_i][array_i] = (read_k_t & anchor_back_mask);
								anchor_point[tid][anchor_hash_i][array_i] = char_i;
								array_i++;
							}
						}
						else
						{
							anchor_hash[tid][anchor_hash_i][read_k_t >> k_anchor_back] = (array_i << 5);
							anchor_hash[tid][anchor_hash_i][read_k_p >> k_anchor_back] |= (array_i - array_i_p);

							array_i_p = array_i;

							anchor_array[tid][anchor_hash_i][array_i] = (read_k_t & anchor_back_mask);
							anchor_point[tid][anchor_hash_i][array_i] = char_i;
							array_i++;
						}
						read_k_p = read_k_t;
					}
					anchor_hash[tid][anchor_hash_i][read_k_p >> k_anchor_back] |= (array_i - array_i_p);
					anchor_point[tid][anchor_hash_i][array_i] = char_i;
				}
#endif

				v_cnt_i = 0;
#ifdef	ALTER_DEBUG_ANCHOR
				if(v_cnt > 1)	v_cnt_i = seed_length_arr[tid][v_cnt_i].index;
#endif

#ifdef	REDUCE_ANCHOR
				op_mask[v_cnt_i] = 1;
#endif
				x = op_vector_pos1[tid][v_cnt_i];
				low = 0;
				high = chr_file_n - 1;

				while ( low <= high )
				{
					mid = (low + high) >> 1;
					if(x < (chr_end_n[mid]))
					{
						high = mid - 1;
					}
					else if(x > (chr_end_n[mid]))
					{
						low = mid + 1;
					}
					else
					{
						chr_re1 =  mid;
						break;
					}
					chr_re1 = low;
				}

				sam_pos1 = op_vector_pos1[tid][v_cnt_i] - chr_end_n[chr_re1 - 1] + 1;

#ifdef	FIX_SA
				op_rc_tmp = op_rc[tid][v_cnt_i];
#endif

				if(op_rc[tid][v_cnt_i] == 0)
				{
#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit1[tid][0];
					read_bit_2[tid] = read_bit2[tid][1];
#else
					strcpy(sam_seq1, seqio[seqi].read_seq1);

					for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
						sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
					sam_seq2[sam_seq_i] = '\0';

					strrev1(sam_seq2);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv1[tid][0];
					qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif
					ksw_s = x + end_dis1[tid] - devi;
					ksw_e = x + insert_dis + devi;

					cigar_m_1 = cigar_m1[tid];
					cigar_m_2 = cigar_m2[tid];
					lv_k_1 = lv_k1;
					lv_k_2 = lv_k2;
					read_length_1 = read_length1;
					read_length_2 = read_length2;

					pound_pos_1_f = pound_pos1_f_forward;
					pound_pos_1_r = pound_pos1_r_forward;
					pound_pos_2_f = pound_pos2_f_reverse;
					pound_pos_2_r = pound_pos2_r_reverse;
				}
				else if(op_rc[tid][v_cnt_i] == 1)
				{
#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit2[tid][1];
					read_bit_2[tid] = read_bit1[tid][0];
#else
					for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
						sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

					sam_seq1[sam_seq_i] = '\0';

					strrev1(sam_seq1);
					strcpy(sam_seq2, seqio[seqi].read_seq1);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv2[tid][1];
					qual_filt_lv_1_o = qual_filt_lv2[tid][0];
#endif
					ksw_s = x - (end_dis1[tid] + devi);
					ksw_e = x - (end_dis1[tid] - devi) + read_length1;

					cigar_m_1 = cigar_m2[tid];
					cigar_m_2 = cigar_m1[tid];
					lv_k_1 = lv_k2;
					lv_k_2 = lv_k1;
					read_length_1 = read_length2;
					read_length_2 = read_length1;

					pound_pos_1_f = pound_pos2_f_reverse;
					pound_pos_1_r = pound_pos2_r_reverse;
					pound_pos_2_f = pound_pos1_f_forward;
					pound_pos_2_r = pound_pos1_r_forward;
				}
				else if(op_rc[tid][v_cnt_i] == 2)
				{
#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit2[tid][0];
					read_bit_2[tid] = read_bit1[tid][1];
#else
					strcpy(sam_seq1, seqio[seqi].read_seq2);
					for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
						sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

					sam_seq2[sam_seq_i] = '\0';

					strrev1(sam_seq2);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv2[tid][0];
					qual_filt_lv_1_o = qual_filt_lv2[tid][1];
#endif

					ksw_s = x + (end_dis2[tid] - devi);
					ksw_e = x + insert_dis + devi;

					cigar_m_1 = cigar_m2[tid];
					cigar_m_2 = cigar_m1[tid];
					lv_k_1 = lv_k2;
					lv_k_2 = lv_k1;
					read_length_1 = read_length2;
					read_length_2 = read_length1;

					pound_pos_1_f = pound_pos2_f_forward;
					pound_pos_1_r = pound_pos2_r_forward;
					pound_pos_2_f = pound_pos1_f_reverse;
					pound_pos_2_r = pound_pos1_r_reverse;
				}
				else
				{
#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit1[tid][1];
					read_bit_2[tid] = read_bit2[tid][0];
#else
					for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
						sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

					sam_seq1[sam_seq_i] = '\0';

					strrev1(sam_seq1);
					strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv1[tid][1];
					qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif
					ksw_s = x - (end_dis2[tid] + devi);
					ksw_e = x - (end_dis2[tid] - devi) + read_length2;

					cigar_m_1 = cigar_m1[tid];
					cigar_m_2 = cigar_m2[tid];
					lv_k_1 = lv_k1;
					lv_k_2 = lv_k2;
					read_length_1 = read_length1;
					read_length_2 = read_length2;

					pound_pos_1_f = pound_pos1_f_reverse;
					pound_pos_1_r = pound_pos1_r_reverse;
					pound_pos_2_f = pound_pos2_f_forward;
					pound_pos_2_r = pound_pos2_r_forward;
				}

#ifdef	ANCHOR_HASH_ALI
				anchor_hash_p = anchor_hash[tid][op_rc[tid][v_cnt_i]];
				anchor_array_p = anchor_array[tid][op_rc[tid][v_cnt_i]];
				anchor_point_p = anchor_point[tid][op_rc[tid][v_cnt_i]];
				anchor_pos_p = anchor_pos[tid][op_rc[tid][v_cnt_i]];
#endif

				d_n1 = 0;
				i_n1 = 0;
				s_offset1 = 0;
				s_offset2 = 0;
				s_r_o_l = op_dm_l1[tid][v_cnt_i];
				s_r_o_r = op_dm_r1[tid][v_cnt_i];

				if((s_r_o_l == 0) && (s_r_o_r == 0))
				{
					strcpy(cigar_p1, cigar_m_1);
				}
				else     //indel
				{
#ifdef	OUTPUT_DEBUG
					if(pound_pos_1_f >= s_r_o_r)   //1
					{
						lv_up_left = 0;
						lv_up_right = s_r_o_l;
						lv_down_right = pound_pos_1_f;
						lv_down_left = s_r_o_r;
						m_n_f = 0;
						m_n_b = read_length_1 - pound_pos_1_f;
						m_m_n = s_r_o_r - s_r_o_l - 1;
					}
					else if(pound_pos_1_r <= s_r_o_l + 1)     //5
					{
						lv_up_left = pound_pos_1_r;//
						lv_up_right = s_r_o_l;
						lv_down_right = read_length_1;
						lv_down_left = s_r_o_r;
						m_n_f = pound_pos_1_r;
						m_n_b = 0;
						m_m_n = s_r_o_r - s_r_o_l - 1;
					}
					else if((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
					{
						lv_up_left = 0;
						lv_up_right = pound_pos_1_f - 1;
						lv_down_right = read_length_1;
						lv_down_left = pound_pos_1_r;
						m_n_f = 0;
						m_n_b = 0;
						m_m_n = pound_pos_1_r - pound_pos_1_f;
					}
					else if((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
					{
						lv_up_left = 0;
						lv_up_right = s_r_o_l;
						lv_down_right = read_length_1;
						lv_down_left = pound_pos_1_r;
						m_n_f = 0;
						m_n_b = 0;
						m_m_n = read_length_1 - s_r_o_l - 1;
					}
					else     //4
					{
						lv_up_left = 0;
						lv_up_right = -1;
						lv_down_right = read_length_1;
						lv_down_left = s_r_o_r;
						m_n_f = 0;
						m_n_b = 0;
						m_m_n = s_r_o_r;
					}
#ifdef	QUAL_FILT_SINGLE_OUT
					for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length_1 - 1- lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]

					for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
					for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
					for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#ifdef	CHAR_CP
					for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]

#endif
#endif

					//deal with front and back lv cigar
					strncpy(str_o, cigarBuf1, f_cigarn);
					s_o = 0;
					f_c = 0;
					pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

					while (pch != NULL)
					{
						pchl = strlen(pch);
						f_cigar[f_cigarn - f_c - 2] = atoi(pch);
						s_o += (pchl + 1);
						f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

						f_c += 2;

						if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
						if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

						pch = strtok_r(NULL, "DMIS", &saveptr);
					}

					strncpy(b_cigar, cigarBuf2, f_cigarn);
					pch = strtok(cigarBuf2,"DMIS");

					if(pch != NULL)
						pchl = strlen(pch);

					snt = 0;

#ifdef	CIGAR_S_MODIFY
					if(lv_up_left)
					{
						sn = sprintf(cigar_p1 + snt, "%uS", lv_up_left);
						snt += sn;
					}
#else
					if(m_n_f)
					{
						if(f_c)
						{
							if(f_cigar[f_cigarn + 1 - f_c] == 'M')
							{
								f_cigar[f_cigarn - f_c] += m_n_f;
							}
							else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
							{
								f_cigar[f_cigarn - f_c] += m_n_f;
							}
							else
							{
								sn = sprintf(cigar_p1 + snt, "%uM", m_n_f);
								snt += sn;
							}
						}
						else	m_m_n += m_n_f;
					}
#endif

					if((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))   //(op_dm_l1[tid][v_cnt_i] != -1) && (op_dm_r1[tid][v_cnt_i] != read_length1)
					{
						if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
						{
							f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s", b_cigar + pchl + 1);
							snt += sn;
						}
						else if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s",b_cigar);
							snt += sn;
						}
						else if(b_cigar[pchl] == 'M')
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}

							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							snt += sn;
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
							snt += sn;
						}
					}
					else if((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
					{
						if(b_cigar[pchl] == 'M')
						{
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							snt += sn;
						}
						else
						{
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
							snt += sn;
						}
					}
					else if((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
					{
						if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
							snt += sn;
						}
					}
					else
					{
						sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
						snt += sn;
					}

#ifdef	CIGAR_S_MODIFY
					if(lv_down_right < read_length_1)
					{
						sn = sprintf(cigar_p1 + snt, "%uS", read_length_1 - lv_down_right);
						snt += sn;
					}
#else
					if(m_n_b)
					{
						if(cigar_p1[snt - 1] == 'M')
						{
							for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
							{
								if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
								m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
							}
							sn = sprintf(cigar_p1 + bit_char_i + 1, "%uM", m_n_b);
							snt = bit_char_i + 1 + sn;
						}
						else if(cigar_p1[snt - 1] == 'S')
						{
							for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
							{
								if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
								m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
							}
							sn = sprintf(cigar_p1 + bit_char_i + 1, "%uS", m_n_b);
							snt = bit_char_i + 1 + sn;
						}
						else
						{
							sn = sprintf(cigar_p1 + snt, "%uM", m_n_b);
							snt += sn;
						}
					}
#endif

#ifdef	CIGAR_LEN_ERR
					cigar_len = 0;
					s_o_tmp = 0;
					strncpy(cigar_tmp, cigar_p1, snt);
					cigar_tmp[snt] = '\0';
					pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

					while (pch_tmp != NULL)
					{
						pchl_tmp = strlen(pch_tmp);
						s_o_tmp += (pchl_tmp + 1);

						if(cigar_p1[s_o_tmp - 1] != 'D')
						{
							cigar_len_tmp = atoi(pch_tmp);
							cigar_len += cigar_len_tmp;
						}
						pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
					}

					if(read_length_1 != cigar_len)
					{
						if(read_length_1 < cigar_len)
						{
							cigar_len_re = cigar_len_tmp - (cigar_len - read_length_1);
							if(cigar_len_re > 0)	sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
							else if(cigar_len_re == 0)	sprintf(cigar_p1 + snt - sn, "\0");
							else	strcpy(cigar_p1, cigar_m_1);
						}
						else
						{
							cigar_len_re = cigar_len_tmp + (read_length_1 - cigar_len);
							sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
						}
					}
#endif
				}
#ifdef	NO_S_OFF
				s_offset1 = 0;
#endif
				sam_pos1 = sam_pos1 + i_n1 - d_n1 + s_offset1;

				//use ksw to find the other end alignment
				ksw_re = 0;
#ifdef	REDUCE_ANCHOR
				other_end_flag = 0;
#endif

#ifdef	ANCHOR_HASH_ALI

				buffer_i = 0;
				r_b_v = 0;

				for(base_i = ksw_s - 1; base_i < ksw_e - k_anchor; base_i += anchor_seed_d)
				{
					if(base_i + k_anchor - 1 < r_b_v)	continue;

					base_re = (base_i & 0X1f);
					if(base_re <= k_anchor_re)
					{
						anchor_ref = ((buffer_ref_seq[base_i >> 5] >> ((k_anchor_re - base_re) << 1)) & anchor_mask);
					}
					else
					{
						anchor_ref = (((buffer_ref_seq[base_i >> 5] & anchor_mask_boundary[32 - base_re]) << ((base_re - k_anchor_re) << 1)) | (buffer_ref_seq[(base_i >> 5) + 1] >> ((32 + k_anchor_re - base_re) << 1)));
					}

					max_right = 0;
					for(tra_i = 0; tra_i < (anchor_hash_p[anchor_ref >> 4] & 0X1f); tra_i++)
					{
						array_index = (anchor_hash_p[anchor_ref >> 4] >> 5) + tra_i;
						if(anchor_array_p[array_index] == (anchor_ref & 0Xf))
						{
							max_seed_length = 0;
							for(print_i = anchor_point_p[array_index]; print_i < anchor_point_p[array_index + 1]; print_i++)
							{
								//extension on both sides
								for(left_i = anchor_pos_p[print_i] - 1, base_i_off_l = base_i - 1; (left_i >= 0) && (base_i_off_l >= ksw_s - 1); left_i--, base_i_off_l--)
								{
#ifdef	CHAR_CP
									if(((read_bit_2[tid][left_i >> 5] >> ((31 - (left_i  & 0X1f)) << 1)) & 0X3) != ((buffer_ref_seq[base_i_off_l >> 5] >> ((31 - (base_i_off_l & 0X1f)) << 1)) & 0X3))
										break;
#else
									if(sam_seq2[left_i] != Dna5Tochar[((buffer_ref_seq[base_i_off_l >> 5] >> ((31 - (base_i_off_l & 0X1f)) << 1)) & 0X3)])
										break;
#endif
								}

								for(right_i = anchor_pos_p[print_i] + k_anchor, base_i_off_r = base_i + k_anchor; (right_i < read_length_2) && (base_i_off_r < ksw_e); right_i++, base_i_off_r++)
								{
#ifdef	CHAR_CP
									if(((read_bit_2[tid][right_i >> 5] >> ((31 - (right_i  & 0X1f)) << 1)) & 0X3) != ((buffer_ref_seq[base_i_off_r >> 5] >> ((31 - (base_i_off_r & 0X1f)) << 1)) & 0X3))
										break;
#else
									if(sam_seq2[right_i] != Dna5Tochar[((buffer_ref_seq[base_i_off_r >> 5] >> ((31 - (base_i_off_r & 0X1f)) << 1)) & 0X3)])
										break;
#endif
								}

								seed_length = right_i - left_i - 1;
								if(seed_length > max_seed_length)
								{
									max_seed_length = seed_length;
									max_right = base_i_off_r;
								}

								anchor_seed_buffer[tid][buffer_i].read_left_off = left_i;
								anchor_seed_buffer[tid][buffer_i].read_right_off = right_i;
								anchor_seed_buffer[tid][buffer_i].ref_left_off = base_i_off_l;
								anchor_seed_buffer[tid][buffer_i].ref_right_off = base_i_off_r;
								anchor_seed_buffer[tid][buffer_i].seed_length = seed_length;
								buffer_i++;
							}
							break;
						}
					}
					r_b_v = max_right;
				}

				if((buffer_i > 0) && (max_seed_length > anchor_seed_length_thr))
				{
					qsort(anchor_seed_buffer[tid], buffer_i, sizeof(anchor_seed), comepare_anchor_seed);

					//LV
					left_i = anchor_seed_buffer[tid][0].read_left_off;
					right_i = anchor_seed_buffer[tid][0].read_right_off;
					base_i_off_l = anchor_seed_buffer[tid][0].ref_left_off;
					base_i_off_r = anchor_seed_buffer[tid][0].ref_right_off;

					d_n1 = 0;
					i_n1 = 0;
					if((left_i == -1) && (right_i == read_length_2))
					{
						strcpy(cigar_p2, cigar_m_2);
					}
					else     //indel
					{
#ifdef	LV_CCIGAR

#ifdef	CHAR_CP
						for(bit_char_i = left_i, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
							read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);
						//33
						for(bit_char_i = base_i_off_l, read_b_i = 0; bit_char_i > base_i_off_l - left_i - 33; bit_char_i--, read_b_i++)
							ali_ref_seq[tid][read_b_i] = ((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
						for(bit_char_i = left_i, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
							read_char[tid][read_b_i] = sam_seq2[bit_char_i];
						//33
						for(bit_char_i = base_i_off_l, read_b_i = 0; bit_char_i > base_i_off_l - left_i - 33; bit_char_i--, read_b_i++)
							ali_ref_seq[tid][read_b_i] = Dna5Tochar[((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
						dm_l2 = computeEditDistanceWithCigar_s_nm_left(ali_ref_seq[tid], left_i + 33, read_char[tid], left_i + 1, lv_k_2, cigarBuf1, f_cigarn, L[tid], &nm_score1, &s_offset2);//, 0

#ifdef	CHAR_CP
						for(bit_char_i = right_i, read_b_i = 0; bit_char_i < read_length_2; bit_char_i++, read_b_i++)
							read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

						for(bit_char_i = base_i_off_r, read_b_i = 0; bit_char_i < base_i_off_r + read_length_2 - right_i + 32; bit_char_i++, read_b_i++)
							ali_ref_seq[tid][read_b_i] = ((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
						for(bit_char_i = right_i, read_b_i = 0; bit_char_i < read_length_2; bit_char_i++, read_b_i++)
							read_char[tid][read_b_i] = sam_seq2[bit_char_i];

						for(bit_char_i = base_i_off_r, read_b_i = 0; bit_char_i < base_i_off_r + read_length_2 - right_i + 32; bit_char_i++, read_b_i++)
							ali_ref_seq[tid][read_b_i] = Dna5Tochar[((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
						dm_r2 = computeEditDistanceWithCigar_s_nm(ali_ref_seq[tid], read_length_2 - right_i + 32, read_char[tid], read_length_2 - right_i, lv_k_2, cigarBuf2, f_cigarn, L[tid], &nm_score2);//, 0
#endif

						//deal with cigar of LV
						//deal with front and back lv cigar
						m_m_n = anchor_seed_buffer[tid][0].seed_length;

						strncpy(str_o, cigarBuf1, f_cigarn);
						s_o = 0;
						f_c = 0;

						pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

						while (pch != NULL)
						{
							pchl = strlen(pch);
							f_cigar[f_cigarn - f_c - 2] = atoi(pch);
							s_o += (pchl + 1);
							f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];
							f_c += 2;

							if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
							if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

							pch = strtok_r(NULL, "DMIS", &saveptr);
						}

						strncpy(b_cigar, cigarBuf2, f_cigarn);
						pch = strtok (cigarBuf2,"DMIS");

						if(pch != NULL)
							pchl = strlen(pch);

						snt = 0;
						if((left_i != -1) && (right_i != read_length_2))
						{
							if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
							{
								f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
								snt += sn;
							}
							else if(f_cigar[f_cigarn - 1] == 'M')
							{
								f_cigar[f_cigarn - 2] += m_m_n;

								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
								snt += sn;
							}
							else if(b_cigar[pchl] == 'M')
							{
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
								snt += sn;
							}
							else
							{
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
								snt += sn;
							}
						}
						else if(left_i == -1)
						{
							if(b_cigar[pchl] == 'M')
								sn = sprintf(cigar_p2, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							else	sn = sprintf(cigar_p2, "%uM%s", m_m_n, b_cigar);

							snt += sn;
						}
						else
						{
							if(f_cigar[f_cigarn - 1] == 'M')
							{
								f_cigar[f_cigarn - 2] += m_m_n;
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
							}
							else
							{
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
								snt += sn;
							}
						}
#ifdef	CIGAR_LEN_ERR

#ifdef FIX_SA
						sv_s_len_p = 0;
#endif
						cigar_len = 0;
						s_o_tmp = 0;
						strncpy(cigar_tmp, cigar_p2, snt);
						cigar_tmp[snt] = '\0';
						pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

						while (pch_tmp != NULL)
						{
							pchl_tmp = strlen(pch_tmp);
							s_o_tmp += (pchl_tmp + 1);

							
							if(cigar_p2[s_o_tmp - 1] != 'D')
							{
								cigar_len_tmp = atoi(pch_tmp);
								cigar_len += cigar_len_tmp;
#ifdef FIX_SA
								if(cigar_p2[s_o_tmp - 1] == 'S')
									sv_s_len_p += cigar_len_tmp;
#endif
								if(cigar_p2[s_o_tmp - 1] == 'I')
									dm_l2 += cigar_len_tmp;
							}else{
								cigar_len_tmp = atoi(pch_tmp);
								dm_l2 += cigar_len_tmp;
							}
							pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
						}

						if(read_length_2 != cigar_len)
						{
							if(read_length_2 < cigar_len)
							{
								cigar_len_re = cigar_len_tmp - (cigar_len - read_length_2);
								if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
								else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
								else	strcpy(cigar_p2, cigar_m_2);
							}
							else
							{
								cigar_len_re = cigar_len_tmp + (read_length_2 - cigar_len);
								sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
							}
						}
#endif
					}
#ifdef	NO_S_OFF
					s_offset2 = 0;
#endif
					sam_pos2 = base_i_off_l - left_i + i_n1 - d_n1 - chr_end_n[chr_re1 - 1] + 2 + s_offset2;
					chr_re2 = chr_re1;
					ksw_re = 1;
				}
				else
				{

#ifdef	REDUCE_ANCHOR

					if((op_rc[tid][v_cnt_i] == 0) || (op_rc[tid][v_cnt_i] == 3))
					{
						other_end_flag = 0;
						while(tra2_i < anchor_n2)
						{
							rcs = rcs2[tid][tra2_i];
							rs = rs2[tid][tra2_i];
							ls = ls2[tid][tra2_i];
							sam_pos2 = poses2[tid][tra2_i];
							dm_l2 = dms2[tid][tra2_i];
							dm_r2 = 0;
							if(orders2[tid][tra2_i] >= MAX_REDUCE_ANCHOR_NUM)
							{
								v_cnt_i_tmp = orders2[tid][tra2_i] - MAX_REDUCE_ANCHOR_NUM;
								if(ops_mask[v_cnt_i_tmp] == 0)
								{
									op_vector_seq1_tmp = ops_vector_seq1[tid][v_cnt_i_tmp];
									ops_mask[v_cnt_i_tmp] = 1;
									other_end_flag = 1;
									tra2_i++;
									break;
								}
							}
							else
							{
								v_cnt_i_tmp = orders2[tid][tra2_i];
								if(op_mask[v_cnt_i_tmp] == 0)
								{
									op_vector_seq1_tmp = op_vector_seq1[tid][v_cnt_i_tmp];
									op_mask[v_cnt_i_tmp] = 1;
									other_end_flag = 1;
									tra2_i++;
									break;
								}
							}
							tra2_i++;
						}
					}
					else
					{
						other_end_flag = 0;
						while(tra1_i < anchor_n1)
						{
							rcs = rcs1[tid][tra1_i];
							rs = rs1[tid][tra1_i];
							ls = ls1[tid][tra1_i];
							sam_pos2 = poses1[tid][tra1_i];
							dm_l2 = dms1[tid][tra1_i];
							dm_r2 = 0;
							if(orders1[tid][tra1_i] >= MAX_REDUCE_ANCHOR_NUM)
							{
								v_cnt_i_tmp = orders1[tid][tra1_i] - MAX_REDUCE_ANCHOR_NUM;
								op_vector_seq1_tmp = ops_vector_seq1[tid][v_cnt_i_tmp];
								if(ops_mask[v_cnt_i_tmp] == 1)
								{
									ops_mask[v_cnt_i_tmp] = 1;
									tra1_i++;
									other_end_flag = 1;
									break;
								}
							}
							else
							{
								v_cnt_i_tmp = orders1[tid][tra1_i];
								op_vector_seq1_tmp = op_vector_seq1[tid][v_cnt_i_tmp];
								if(op_mask[v_cnt_i_tmp] == 1)
								{
									op_mask[v_cnt_i_tmp] = 1;
									tra1_i++;
									other_end_flag = 1;
									break;
								}
							}
							tra1_i++;
						}
					}

					if(other_end_flag)
					{
						x = sam_pos2;
						low = 0;
						high = chr_file_n - 1;

						while ( low <= high )
						{
							mid = (low + high) >> 1;
							if(x < (chr_end_n[mid]))
							{
								high = mid - 1;
							}
							else if(x > (chr_end_n[mid]))
							{
								low = mid + 1;
							}
							else
							{
								chr_re2 =  mid;
								break;
							}
							chr_re2 = low;
						}
						sam_pos2 = x - chr_end_n[chr_re2 - 1] + 1;

						if(rcs == 0)
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit1[tid][0];
							read_bit_2[tid] = read_bit2[tid][1];
#else
							strcpy(sam_seq1, seqio[seqi].read_seq1);

							/*
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
							*/
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv1[tid][0];
							qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif

							cigar_m_1 = cigar_m1[tid];
							cigar_m_2 = cigar_m2[tid];
							lv_k_1 = lv_k1;
							lv_k_2 = lv_k2;
							read_length_1 = read_length1;
							read_length_2 = read_length2;

							pound_pos_1_f = pound_pos1_f_forward;
							pound_pos_1_r = pound_pos1_r_forward;
							pound_pos_2_f = pound_pos2_f_reverse;
							pound_pos_2_r = pound_pos2_r_reverse;
						}
						else if(rcs == 1)
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit2[tid][1];
							read_bit_2[tid] = read_bit1[tid][0];
#else
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

							sam_seq1[sam_seq_i] = '\0';

							strrev1(sam_seq1);
							//strcpy(sam_seq2, seqio[seqi].read_seq1);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv2[tid][1];
							qual_filt_lv_1_o = qual_filt_lv2[tid][0];
#endif

							cigar_m_1 = cigar_m2[tid];
							cigar_m_2 = cigar_m1[tid];
							lv_k_1 = lv_k2;
							lv_k_2 = lv_k1;
							read_length_1 = read_length2;
							read_length_2 = read_length1;

							pound_pos_1_f = pound_pos2_f_reverse;
							pound_pos_1_r = pound_pos2_r_reverse;
							pound_pos_2_f = pound_pos1_f_forward;
							pound_pos_2_r = pound_pos1_r_forward;
						}
						else if(rcs == 2)
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit2[tid][0];
							read_bit_2[tid] = read_bit1[tid][1];
#else
							strcpy(sam_seq1, seqio[seqi].read_seq2);
							/*
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
							*/
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv2[tid][0];
							qual_filt_lv_1_o = qual_filt_lv2[tid][1];
#endif

							cigar_m_1 = cigar_m2[tid];
							cigar_m_2 = cigar_m1[tid];
							lv_k_1 = lv_k2;
							lv_k_2 = lv_k1;
							read_length_1 = read_length2;
							read_length_2 = read_length1;

							pound_pos_1_f = pound_pos2_f_forward;
							pound_pos_1_r = pound_pos2_r_forward;
							pound_pos_2_f = pound_pos1_f_reverse;
							pound_pos_2_r = pound_pos1_r_reverse;
						}
						else
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit1[tid][1];
							read_bit_2[tid] = read_bit2[tid][0];
#else
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

							sam_seq1[sam_seq_i] = '\0';

							strrev1(sam_seq1);
							//strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv1[tid][1];
							qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif

							cigar_m_1 = cigar_m1[tid];
							cigar_m_2 = cigar_m2[tid];
							lv_k_1 = lv_k1;
							lv_k_2 = lv_k2;
							read_length_1 = read_length1;
							read_length_2 = read_length2;

							pound_pos_1_f = pound_pos1_f_reverse;
							pound_pos_1_r = pound_pos1_r_reverse;
							pound_pos_2_f = pound_pos2_f_forward;
							pound_pos_2_r = pound_pos2_r_forward;
						}

						d_n1 = 0;
						i_n1 = 0;
						s_offset1 = 0;
						s_offset2 = 0;
						s_r_o_l = ls;
						s_r_o_r = rs;

						if((s_r_o_l == 0) && (s_r_o_r == 0))
						{
							strcpy(cigar_p2, cigar_m_2);
						}
						else     //indel
						{
#ifdef	OUTPUT_DEBUG
							if(pound_pos_1_f >= s_r_o_r)   //1
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = pound_pos_1_f;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = read_length_1 - pound_pos_1_f;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if(pound_pos_1_r <= s_r_o_l + 1)     //5
							{
								lv_up_left = pound_pos_1_r;//
								lv_up_right = s_r_o_l;
								lv_down_right = read_length_1;
								lv_down_left = s_r_o_r;
								m_n_f = pound_pos_1_r;
								m_n_b = 0;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
							{
								lv_up_left = 0;
								lv_up_right = pound_pos_1_f - 1;
								lv_down_right = read_length_1;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = pound_pos_1_r - pound_pos_1_f;
							}
							else if((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = read_length_1;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = read_length_1 - s_r_o_l - 1;
							}
							else     //4
							{
								lv_up_left = 0;
								lv_up_right = -1;
								lv_down_right = read_length_1;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = s_r_o_r;
							}
#ifdef	QUAL_FILT_SINGLE_OUT
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length_1 - 1- lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]

							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#ifdef	CHAR_CP
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]

#endif
#endif

							//deal with front and back lv cigar
							strncpy(str_o, cigarBuf1, f_cigarn);
							s_o = 0;
							f_c = 0;
							pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

							while (pch != NULL)
							{
								pchl = strlen(pch);
								f_cigar[f_cigarn - f_c - 2] = atoi(pch);
								s_o += (pchl + 1);
								f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

								f_c += 2;

								if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
								if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

								pch = strtok_r(NULL, "DMIS", &saveptr);
							}

							strncpy(b_cigar, cigarBuf2, f_cigarn);
							pch = strtok(cigarBuf2,"DMIS");

							if(pch != NULL)
								pchl = strlen(pch);

							snt = 0;
#ifdef	CIGAR_S_MODIFY
							if(lv_up_left)
							{
								sn = sprintf(cigar_p2 + snt, "%uS", lv_up_left);
								snt += sn;
							}
#else
							if(m_n_f)
							{
								if(f_c)
								{
									if(f_cigar[f_cigarn + 1 - f_c] == 'M')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else
									{
										sn = sprintf(cigar_p2 + snt, "%uM", m_n_f);
										snt += sn;
									}
								}
								else	m_m_n += m_n_f;
							}
#endif
							if((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))   //(op_dm_l1[tid][v_cnt_i] != -1) && (op_dm_r1[tid][v_cnt_i] != read_length1)
							{
								if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
								{
									f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
									snt += sn;
								}
								else if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
									snt += sn;
								}
								else if(b_cigar[pchl] == 'M')
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}

									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
							{
								if(b_cigar[pchl] == 'M')
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
							{
								if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
									snt += sn;
								}
							}
							else
							{
								sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
								snt += sn;
							}
#ifdef	CIGAR_S_MODIFY
							if(lv_down_right < read_length_1)
							{
								sn = sprintf(cigar_p2 + snt, "%uS", read_length_1 - lv_down_right);
								snt += sn;
							}
#else
							if(m_n_b)
							{
								if(cigar_p2[snt - 1] == 'M')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uM", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else if(cigar_p2[snt - 1] == 'S')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uS", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM", m_n_b);
									snt += sn;
								}
							}
#endif

#ifdef	CIGAR_LEN_ERR

#ifdef FIX_SA
							sv_s_len_p = 0;
#endif
							cigar_len = 0;
							s_o_tmp = 0;
							strncpy(cigar_tmp, cigar_p2, snt);
							cigar_tmp[snt] = '\0';
							pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

							while (pch_tmp != NULL)
							{
								pchl_tmp = strlen(pch_tmp);
								s_o_tmp += (pchl_tmp + 1);

								if(cigar_p2[s_o_tmp - 1] != 'D')
								{
									cigar_len_tmp = atoi(pch_tmp);
									cigar_len += cigar_len_tmp;
#ifdef FIX_SA
									if(cigar_p2[s_o_tmp - 1] == 'S')
										sv_s_len_p += cigar_len_tmp;
#endif
								}
								pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
							}

							if(read_length_1 != cigar_len)
							{
								if(read_length_1 < cigar_len)
								{
									cigar_len_re = cigar_len_tmp - (cigar_len - read_length_1);
									if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
									else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
									else	strcpy(cigar_p2, cigar_m_1);
								}
								else
								{
									cigar_len_re = cigar_len_tmp + (read_length_1 - cigar_len);
									sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
								}
							}
#endif
						}
#ifdef	NO_S_OFF
						s_offset1 = 0;
#endif
						sam_pos2 = sam_pos2 + i_n1 - d_n1 + s_offset1;

						ksw_re = 1;
					}
					else
					{
						sam_pos2 = sam_pos1;
#ifdef	PICARD_BUG
						strcpy(cigar_p2, cigar_p1);
#else
						strcpy(cigar_p2, "*");
#endif
						chr_re2 = chr_re1;
						ksw_re = 0;
					}
#else
					sam_pos2 = sam_pos1;
#ifdef	PICARD_BUG
					strcpy(cigar_p2, cigar_p1);
#else
					strcpy(cigar_p2, "*");
#endif
					chr_re2 = chr_re1;
					ksw_re = 0;

#endif
				}
#endif

				if(op_rc[tid][v_cnt_i] == 0)
				{
#ifdef	FIX_SV
					sv_add = 2;
#endif

#ifdef	CHAR_CP
					strcpy(sam_seq1, seqio[seqi].read_seq1);

#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 2)
							strcpy(sam_seq2, seqio[seqi].read_seq2);
						else
						{
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
						}

						//if((rcs == 0) || (rcs == 3))	printf("Error: %u\n", rcs);

					}
					else
					{
						for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
							sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
						sam_seq2[sam_seq_i] = '\0';

						strrev1(sam_seq2);
					}
#else
					for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
						sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
					sam_seq2[sam_seq_i] = '\0';

					strrev1(sam_seq2);
#endif

#endif
					seq1p = sam_seq1;
					if(ksw_re == 1)
					{
						sam_cross = sam_pos2 + read_length2 - sam_pos1;
						if((sam_cross > insert_dis + devi) || (sam_cross < insert_dis - devi) || (chr_re1 != chr_re2))
						{
							if((other_end_flag) && (rcs == 2))
							{
								sam_flag1 = 65;
								sam_flag2 = 129;
							}
							else
							{
								sam_flag1 = 97;
								sam_flag2 = 145;
							}
						}
						else
						{
							if((other_end_flag) && (rcs == 2))
							{
								sam_flag1 = 67;
								sam_flag2 = 131;
							}
							else
							{
								sam_flag1 = 99;
								sam_flag2 = 147;
							}
						}

						seq2p = sam_seq2;
					}
					else
					{
						sam_flag1 = 73;
						sam_flag2 = 133;
						sam_cross = 0;
						seq2p = seqio[seqi].read_seq2;

						qual_flag = 2;
#ifdef	SAMTOOLS_BUG
						sam_pos2 = 0;
						chr_re2 = chr_file_n;
#endif

					}
					cp1 = cigar_p1;
					cp2 = cigar_p2;

					seqio[seqi].nm1 = dm_op[tid];
					seqio[seqi].nm2 = dm_l2 + dm_r2;
				}
				else if(op_rc[tid][v_cnt_i] == 1)
				{
#ifdef	FIX_SV
					sv_add = 1;
#endif

#ifdef	CHAR_CP
					for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
						sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

					sam_seq1[sam_seq_i] = '\0';
					strrev1(sam_seq1);
#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 0)
							strcpy(sam_seq2, seqio[seqi].read_seq1);
						else
						{
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];
							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
						}

						//if((rcs == 1) || (rcs == 2))	printf("Error: %u\n", rcs);

					}
					else
					{
						strcpy(sam_seq2, seqio[seqi].read_seq1);
					}
#else
					strcpy(sam_seq2, seqio[seqi].read_seq1);
#endif

#endif
					seq2p = sam_seq1;
					if(ksw_re == 1)
					{
						chr_re1 = chr_re1 ^ chr_re2;
						chr_re2 = chr_re1 ^ chr_re2;
						chr_re1 = chr_re1 ^ chr_re2;

						sam_pos1 = sam_pos1 ^ sam_pos2;
						sam_pos2 = sam_pos1 ^ sam_pos2;
						sam_pos1 = sam_pos1 ^ sam_pos2;

						sam_cross = sam_pos2 + read_length2 - sam_pos1;

						if((sam_cross > insert_dis + devi) || (sam_cross < insert_dis - devi) || (chr_re1 != chr_re2))
						{
							if((other_end_flag) && (rcs != 0))
							{
								sam_flag1 = 113;
								sam_flag2 = 177;
							}
							else
							{
								sam_flag1 = 97;
								sam_flag2 = 145;
							}
						}
						else
						{
							if((other_end_flag) && (rcs != 0))
							{
								sam_flag1 = 115;
								sam_flag2 = 179;
							}
							else
							{
								sam_flag1 = 99;
								sam_flag2 = 147;
							}
						}
						seq1p = sam_seq2;
					}
					else
					{
						sam_flag1 = 117;
						sam_flag2 = 153;
						sam_cross = 0;

						seq1p = seqio[seqi].read_seq1;

						qual_flag = 1;
#ifdef	SAMTOOLS_BUG
						sam_pos2 = sam_pos1;
						sam_pos1 = 0;
						chr_re2 = chr_re1;
						chr_re1 = chr_file_n;
#endif
					}
					cp1 = cigar_p2;
					cp2 = cigar_p1;

					seqio[seqi].nm1 = dm_l2 + dm_r2;
					seqio[seqi].nm2 = dm_op[tid];
				}
				else if(op_rc[tid][v_cnt_i] == 2)
				{
#ifdef	FIX_SV
					sv_add = 1;
#endif

#ifdef	CHAR_CP
					strcpy(sam_seq1, seqio[seqi].read_seq2);

#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 0)
							strcpy(sam_seq2, seqio[seqi].read_seq1);
						else
						{
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];
							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
						}

						//if((rcs == 1) || (rcs == 2))	printf("Error: %u\n", rcs);

					}
					else
					{
						for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
							sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];
						sam_seq2[sam_seq_i] = '\0';
						strrev1(sam_seq2);
					}
#else
					for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
						sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];
					sam_seq2[sam_seq_i] = '\0';
					strrev1(sam_seq2);
#endif

#endif
					seq2p = sam_seq1;
					if(ksw_re == 1)
					{
						chr_re1 = chr_re1 ^ chr_re2;
						chr_re2 = chr_re1 ^ chr_re2;
						chr_re1 = chr_re1 ^ chr_re2;

						sam_pos1 = sam_pos1 ^ sam_pos2;
						sam_pos2 = sam_pos1 ^ sam_pos2;
						sam_pos1 = sam_pos1 ^ sam_pos2;

						sam_cross = sam_pos2 - read_length1 - sam_pos1;
						if((sam_cross > insert_dis + devi) || (sam_cross < insert_dis - devi) || (chr_re1 != chr_re2))
						{
							if((other_end_flag) && (rcs == 0))
							{
								sam_flag1 = 65;
								sam_flag2 = 129;
							}
							else
							{
								sam_flag1 = 81;
								sam_flag2 = 161;
							}
						}
						else
						{
							if((other_end_flag) && (rcs == 0))
							{
								sam_flag1 = 67;
								sam_flag2 = 131;
							}
							else
							{
								sam_flag1 = 83;
								sam_flag2 = 163;
							}
						}

						seq1p = sam_seq2;
					}
					else
					{
						sam_flag1 = 69;
						sam_flag2 = 137;
						sam_cross = 0;

						seq1p = seqio[seqi].read_seq1;

						qual_flag = 1;
#ifdef	SAMTOOLS_BUG
						sam_pos2 = sam_pos1;
						sam_pos1 = 0;
						chr_re2 = chr_re1;
						chr_re1 = chr_file_n;
#endif
					}
					cp1 = cigar_p2;
					cp2 = cigar_p1;

					seqio[seqi].nm1 = dm_l2 + dm_r2;
					seqio[seqi].nm2 = dm_op[tid];
				}
				else
				{
#ifdef	FIX_SV
					sv_add = 2;
#endif

#ifdef	CHAR_CP
					for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
						sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];
					sam_seq1[sam_seq_i] = '\0';
					strrev1(sam_seq1);

#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 2)
							strcpy(sam_seq2, seqio[seqi].read_seq2);
						else
						{
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
							sam_seq2[sam_seq_i] = '\0';
							strrev1(sam_seq2);
						}

						//if((rcs == 0) || (rcs == 3))	printf("Error: %u\n", rcs);

					}
					else
					{
						strcpy(sam_seq2, seqio[seqi].read_seq2);
					}
#else
					strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#endif
					seq1p = sam_seq1;
					seq2p = sam_seq2;
					if(ksw_re == 1)
					{
						sam_cross = sam_pos2 - read_length1 - sam_pos1;
						if((sam_cross > insert_dis + devi) || (sam_cross < insert_dis - devi) || (chr_re1 != chr_re2))
						{
							if((other_end_flag) && (rcs != 2))
							{
								sam_flag1 = 113;
								sam_flag2 = 177;
							}
							else
							{
								sam_flag1 = 81;
								sam_flag2 = 161;
							}
						}
						else
						{

							if((other_end_flag) && (rcs != 2))
							{
								sam_flag1 = 115;
								sam_flag2 = 179;
							}
							else
							{
								sam_flag1 = 83;
								sam_flag2 = 163;
							}
						}
					}
					else
					{
						sam_flag1 = 89;
						sam_flag2 = 181;
						sam_cross = 0;

						qual_flag = 2;
#ifdef	SAMTOOLS_BUG
						sam_pos2 = 0;
						chr_re2 = chr_file_n;
#endif
					}
					cp1 = cigar_p1;
					cp2 = cigar_p2;

					seqio[seqi].nm1 = dm_op[tid];
					seqio[seqi].nm2 = dm_l2 + dm_r2;
				}

#ifdef	OUTPUT_ARR

				seqio[seqi].flag1 = sam_flag1;
				seqio[seqi].flag2 = sam_flag2;
				//seqio[seqi].chr_re = chr_re;
				seqio[seqi].chr_re1 = chr_re1;
				seqio[seqi].chr_re2 = chr_re2;

				if(sam_pos1 <= 0)
				{
					sam_pos1 = 1;
				}
				else
				{
					if((sam_pos1 + read_length1 - 1) > (chr_end_n[chr_re1] - chr_end_n[chr_re1 - 1]))
						sam_pos1 = chr_end_n[chr_re1] - chr_end_n[chr_re1 - 1] - 1 - read_length1;
				}
				if(sam_pos2 <= 0)
				{
					sam_pos2 = 1;
				}
				else
				{
					if((sam_pos2 + read_length2 - 1) > (chr_end_n[chr_re2] - chr_end_n[chr_re2 - 1]))
						sam_pos2 = chr_end_n[chr_re2] - chr_end_n[chr_re2 - 1] - 1 - read_length2;
				}

				seqio[seqi].pos1 = sam_pos1;
				seqio[seqi].pos2 = sam_pos2;
				seqio[seqi].cross = sam_cross;

				strcpy(pr_cigar1_buffer[seqi], cp1);
				seqio[seqi].cigar1 = pr_cigar1_buffer[seqi];
				strcpy(pr_cigar2_buffer[seqi], cp2);
				seqio[seqi].cigar2 = pr_cigar2_buffer[seqi];

				if((sam_flag1 == 99) || (sam_flag1 == 117) || (sam_flag1 == 97))
				{
					strcpy(read_rev_buffer[seqi], seq2p);
					read_rev_buffer[seqi][read_length2] = '\0';

					seqio[seqi].seq2 = read_rev_buffer[seqi];
					seqio[seqi].seq1 = seqio[seqi].read_seq1;

					strrev1(qual2_buffer[seqi]);
				}
				else if((sam_flag1 == 83) || (sam_flag1 == 89) || (sam_flag1 == 81))
				{
					strcpy(read_rev_buffer[seqi], seq1p);
					read_rev_buffer[seqi][read_length1] = '\0';

					seqio[seqi].seq1 = read_rev_buffer[seqi];
					seqio[seqi].seq2 = seqio[seqi].read_seq2;

					strrev1(qual1_buffer[seqi]);
				}
				else if((sam_flag1 == 113) || (sam_flag1 == 115))
				{
					strcpy(read_rev_buffer[seqi], seq1p);
					read_rev_buffer[seqi][read_length1] = '\0';

					strcpy(read_rev_buffer_1[seqi], seq2p);
					read_rev_buffer_1[seqi][read_length2] = '\0';

					seqio[seqi].seq1 = read_rev_buffer[seqi];
					seqio[seqi].seq2 = read_rev_buffer_1[seqi];
				}
				else
				{
					seqio[seqi].seq1 = seqio[seqi].read_seq1;
					seqio[seqi].seq2 = seqio[seqi].read_seq2;
				}
#endif
			}

			lv_re1f = dm_op[tid];
			xa_i_1 = 0;
			xa_i_2 = 0;
			for(va_cnt_i = 1; va_cnt_i < v_cnt; va_cnt_i++)
			{
#ifdef	ALTER_DEBUG_ANCHOR
				v_cnt_i = seed_length_arr[tid][va_cnt_i].index;
#else
				v_cnt_i = va_cnt_i;
#endif

#ifdef	REDUCE_ANCHOR
				if(op_mask[v_cnt_i] == 1)	continue;

				op_mask[v_cnt_i] = 1;
#endif
				x = op_vector_pos1[tid][v_cnt_i];
				low = 0;
				high = chr_file_n - 1;

				while ( low <= high )
				{
					mid = (low + high) >> 1;
					if(x < (chr_end_n[mid]))
					{
						high = mid - 1;
					}
					else if(x > (chr_end_n[mid]))
					{
						low = mid + 1;
					}
					else
					{
						chr_re =  mid;
						break;
					}
					chr_re = low;
				}

				if(chr_re)
					sam_pos1 = op_vector_pos1[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;
				else
				{
					sam_pos1 = op_vector_pos1[tid][v_cnt_i];
					chr_re = 1;
				}
				//chr_res[tid][xa_i] = chr_re;

				if(op_rc[tid][v_cnt_i] == 0)
				{
#ifdef OUPUT_SINGLE_REPEAT

#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit1[tid][0];
					read_bit_2[tid] = read_bit2[tid][1];
#else
					strcpy(sam_seq1, seqio[seqi].read_seq1);
					for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
						sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

					sam_seq2[sam_seq_i] = '\0';

					strrev1(sam_seq2);
#endif


#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv1[tid][0];
					qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif
					ksw_s = x + end_dis1[tid] - devi;
					ksw_e = x + insert_dis + devi;

					cigar_m_1 = cigar_m1[tid];
					cigar_m_2 = cigar_m2[tid];
					lv_k_1 = lv_k1;
					lv_k_2 = lv_k2;
					read_length_1 = read_length1;
					read_length_2 = read_length2;

					pound_pos_1_f = pound_pos1_f_forward;
					pound_pos_1_r = pound_pos1_r_forward;
					pound_pos_2_f = pound_pos2_f_reverse;
					pound_pos_2_r = pound_pos2_r_reverse;
				}
				else if(op_rc[tid][v_cnt_i] == 1)
				{
#ifdef OUPUT_SINGLE_REPEAT

#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit2[tid][1];
					read_bit_2[tid] = read_bit1[tid][0];
#else
					for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
						sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

					sam_seq1[sam_seq_i] = '\0';

					strrev1(sam_seq1);
					strcpy(sam_seq2, seqio[seqi].read_seq1);
#endif

#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv2[tid][1];
					qual_filt_lv_1_o = qual_filt_lv2[tid][0];
#endif

					ksw_s = x - (end_dis1[tid] + devi);
					ksw_e = x - (end_dis1[tid] - devi) + read_length1;

					cigar_m_1 = cigar_m2[tid];
					cigar_m_2 = cigar_m1[tid];
					lv_k_1 = lv_k2;
					lv_k_2 = lv_k1;
					read_length_1 = read_length2;
					read_length_2 = read_length1;

					pound_pos_1_f = pound_pos2_f_reverse;
					pound_pos_1_r = pound_pos2_r_reverse;
					pound_pos_2_f = pound_pos1_f_forward;
					pound_pos_2_r = pound_pos1_r_forward;
				}
				else if(op_rc[tid][v_cnt_i] == 2)
				{
#ifdef OUPUT_SINGLE_REPEAT

#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit2[tid][0];
					read_bit_2[tid] = read_bit1[tid][1];
#else
					strcpy(sam_seq1, seqio[seqi].read_seq2);
					for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
						sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

					sam_seq2[sam_seq_i] = '\0';
					strrev1(sam_seq2);
#endif

#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv2[tid][0];
					qual_filt_lv_1_o = qual_filt_lv2[tid][1];
#endif

					ksw_s = x + (end_dis2[tid] - devi);
					ksw_e = x + insert_dis + devi;

					cigar_m_1 = cigar_m2[tid];
					cigar_m_2 = cigar_m1[tid];
					lv_k_1 = lv_k2;
					lv_k_2 = lv_k1;
					read_length_1 = read_length2;
					read_length_2 = read_length1;

					pound_pos_1_f = pound_pos2_f_forward;
					pound_pos_1_r = pound_pos2_r_forward;
					pound_pos_2_f = pound_pos1_f_reverse;
					pound_pos_2_r = pound_pos1_r_reverse;
				}
				else
				{
#ifdef OUPUT_SINGLE_REPEAT

#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit1[tid][1];
					read_bit_2[tid] = read_bit2[tid][0];
#else
					for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
						sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

					sam_seq1[sam_seq_i] = '\0';

					strrev1(sam_seq1);
					strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv1[tid][1];
					qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif

					ksw_s = x - (end_dis2[tid] + devi);
					ksw_e = x - (end_dis2[tid] - devi) + read_length2;

					cigar_m_1 = cigar_m1[tid];
					cigar_m_2 = cigar_m2[tid];
					lv_k_1 = lv_k1;
					lv_k_2 = lv_k2;
					read_length_1 = read_length1;
					read_length_2 = read_length2;

					pound_pos_1_f = pound_pos1_f_reverse;
					pound_pos_1_r = pound_pos1_r_reverse;
					pound_pos_2_f = pound_pos2_f_forward;
					pound_pos_2_r = pound_pos2_r_forward;
				}

#ifdef	ANCHOR_HASH_ALI
				anchor_hash_p = anchor_hash[tid][op_rc[tid][v_cnt_i]];
				anchor_array_p = anchor_array[tid][op_rc[tid][v_cnt_i]];
				anchor_point_p = anchor_point[tid][op_rc[tid][v_cnt_i]];
				anchor_pos_p = anchor_pos[tid][op_rc[tid][v_cnt_i]];
#endif
				d_n1 = 0;
				i_n1 = 0;
				s_offset1 = 0;
				s_offset2 = 0;
				s_r_o_l = op_dm_l1[tid][v_cnt_i];
				s_r_o_r = op_dm_r1[tid][v_cnt_i];

				if((s_r_o_l == 0) && (s_r_o_r == 0))
				{
					strcpy(cigar_p1, cigar_m_1);
				}
				else     //indel
				{
#ifdef	OUTPUT_DEBUG
					if(pound_pos_1_f >= s_r_o_r)   //1
					{
						lv_up_left = 0;
						lv_up_right = s_r_o_l;
						lv_down_right = pound_pos_1_f;
						lv_down_left = s_r_o_r;
						m_n_f = 0;
						m_n_b = read_length_1 - pound_pos_1_f;
						m_m_n = s_r_o_r - s_r_o_l - 1;
					}
					else if(pound_pos_1_r <= s_r_o_l + 1)     //5
					{
						lv_up_left = pound_pos_1_r;//
						lv_up_right = s_r_o_l;
						lv_down_right = read_length_1;
						lv_down_left = s_r_o_r;
						m_n_f = pound_pos_1_r;
						m_n_b = 0;
						m_m_n = s_r_o_r - s_r_o_l - 1;
					}
					else if((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
					{
						lv_up_left = 0;
						lv_up_right = pound_pos_1_f - 1;
						lv_down_right = read_length_1;
						lv_down_left = pound_pos_1_r;
						m_n_f = 0;
						m_n_b = 0;
						m_m_n = pound_pos_1_r - pound_pos_1_f;
					}
					else if((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
					{
						lv_up_left = 0;
						lv_up_right = s_r_o_l;
						lv_down_right = read_length_1;
						lv_down_left = pound_pos_1_r;
						m_n_f = 0;
						m_n_b = 0;
						m_m_n = read_length_1 - s_r_o_l - 1;
					}
					else     //4
					{
						lv_up_left = 0;
						lv_up_right = -1;
						lv_down_right = read_length_1;
						lv_down_left = s_r_o_r;
						m_n_f = 0;
						m_n_b = 0;
						m_m_n = s_r_o_r;
					}

#ifdef	QUAL_FILT_SINGLE_OUT
					for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length_1 - 1- lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]

					for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
					for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]

#ifdef	CHAR_CP
					for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#endif

					//deal with front and back lv cigar
					strncpy(str_o, cigarBuf1, f_cigarn);
					s_o = 0;
					f_c = 0;

					pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

					while (pch != NULL)
					{
						pchl = strlen(pch);
						f_cigar[f_cigarn - f_c - 2] = atoi(pch);
						s_o += (pchl + 1);
						f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

						f_c += 2;

						if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
						if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

						pch = strtok_r(NULL, "DMIS", &saveptr);
					}

					strncpy(b_cigar, cigarBuf2, f_cigarn);
					pch = strtok(cigarBuf2,"DMIS");

					if(pch != NULL)
						pchl = strlen(pch);

					snt = 0;
#ifdef	CIGAR_S_MODIFY
					if(lv_up_left)
					{
						sn = sprintf(cigar_p1 + snt, "%uS", lv_up_left);
						snt += sn;
					}
#else
					if(m_n_f)
					{
						if(f_c)
						{
							if(f_cigar[f_cigarn + 1 - f_c] == 'M')
							{
								f_cigar[f_cigarn - f_c] += m_n_f;
							}
							else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
							{
								f_cigar[f_cigarn - f_c] += m_n_f;
							}
							else
							{
								sn = sprintf(cigar_p1 + snt, "%uM", m_n_f);
								snt += sn;
							}
						}
						else	m_m_n += m_n_f;

					}
#endif
					if((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))   //(op_dm_l1[tid][v_cnt_i] != -1) && (op_dm_r1[tid][v_cnt_i] != read_length1)
					{
						if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
						{
							f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s", b_cigar + pchl + 1);
							snt += sn;
						}
						else if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s",b_cigar);
							snt += sn;
						}
						else if(b_cigar[pchl] == 'M')
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}

							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							snt += sn;
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
							snt += sn;
						}
					}
					else if((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
					{
						if(b_cigar[pchl] == 'M')
						{
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							snt += sn;
						}
						else
						{
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
							snt += sn;
						}
					}
					else if((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
					{
						if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
							snt += sn;
						}
					}
					else
					{
						sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
						snt += sn;
					}
#ifdef	CIGAR_S_MODIFY
					if(lv_down_right < read_length_1)
					{
						sn = sprintf(cigar_p1 + snt, "%uS", read_length_1 - lv_down_right);
						snt += sn;
					}
#else
					if(m_n_b)
					{
						if(cigar_p1[snt - 1] == 'M')
						{
							for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
							{
								if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
								m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
							}
							sn = sprintf(cigar_p1 + bit_char_i + 1, "%uM", m_n_b);
							snt = bit_char_i + 1 + sn;
						}
						else if(cigar_p1[snt - 1] == 'S')
						{
							for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
							{
								if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
								m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
							}
							sn = sprintf(cigar_p1 + bit_char_i + 1, "%uS", m_n_b);
							snt = bit_char_i + 1 + sn;
						}
						else
						{
							sn = sprintf(cigar_p1 + snt, "%uM", m_n_b);
							snt += sn;
						}
					}
#endif

#ifdef	CIGAR_LEN_ERR
					cigar_len = 0;
					s_o_tmp = 0;
					strncpy(cigar_tmp, cigar_p1, snt);
					cigar_tmp[snt] = '\0';
					pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

					while (pch_tmp != NULL)
					{
						pchl_tmp = strlen(pch_tmp);
						s_o_tmp += (pchl_tmp + 1);

						if(cigar_p1[s_o_tmp - 1] != 'D')
						{
							cigar_len_tmp = atoi(pch_tmp);
							cigar_len += cigar_len_tmp;
						}
						pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
					}

					if(read_length_1 != cigar_len)
					{
						if(read_length_1 < cigar_len)
						{
							cigar_len_re = cigar_len_tmp - (cigar_len - read_length_1);
							if(cigar_len_re > 0)	sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
							else if(cigar_len_re == 0)	sprintf(cigar_p1 + snt - sn, "\0");
							else	strcpy(cigar_p1, cigar_m_1);
						}
						else
						{
							cigar_len_re = cigar_len_tmp + (read_length_1 - cigar_len);
							sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
						}
					}
#endif

				}
#ifdef	NO_S_OFF
				s_offset1 = 0;
#endif
				sam_pos1 = sam_pos1 + i_n1 - d_n1 + s_offset1;

				ksw_re = 0;
#ifdef	REDUCE_ANCHOR
				other_end_flag = 0;
#endif

#ifdef	ANCHOR_HASH_ALI

				buffer_i = 0;
				r_b_v = 0;

				for(base_i = ksw_s - 1; base_i < ksw_e - k_anchor; base_i += anchor_seed_d)
				{
					if(base_i + k_anchor - 1 < r_b_v)	continue;

					base_re = (base_i & 0X1f);
					if(base_re <= k_anchor_re)
					{
						anchor_ref = ((buffer_ref_seq[base_i >> 5] >> ((k_anchor_re - base_re) << 1)) & anchor_mask);
					}
					else
					{
						anchor_ref = (((buffer_ref_seq[base_i >> 5] & anchor_mask_boundary[32 - base_re]) << ((base_re - k_anchor_re) << 1)) | (buffer_ref_seq[(base_i >> 5) + 1] >> ((32 + k_anchor_re - base_re) << 1)));
					}

					max_right = 0;
					for(tra_i = 0; tra_i < (anchor_hash_p[anchor_ref >> 4] & 0X1f); tra_i++)
					{
						array_index = (anchor_hash_p[anchor_ref >> 4] >> 5) + tra_i;

						if(anchor_array_p[array_index] == (anchor_ref & 0Xf))
						{
							max_seed_length = 0;
							for(print_i = anchor_point_p[array_index]; print_i < anchor_point_p[array_index + 1]; print_i++)
							{
								//extension on both sides
								for(left_i = anchor_pos_p[print_i] - 1, base_i_off_l = base_i - 1; (left_i >= 0) && (base_i_off_l >= ksw_s - 1); left_i--, base_i_off_l--)
								{
#ifdef	CHAR_CP
									if(((read_bit_2[tid][left_i >> 5] >> ((31 - (left_i  & 0X1f)) << 1)) & 0X3) != ((buffer_ref_seq[base_i_off_l >> 5] >> ((31 - (base_i_off_l & 0X1f)) << 1)) & 0X3))
										break;
#else
									if(sam_seq2[left_i] != Dna5Tochar[((buffer_ref_seq[base_i_off_l >> 5] >> ((31 - (base_i_off_l & 0X1f)) << 1)) & 0X3)])
										break;
#endif
								}

								for(right_i = anchor_pos_p[print_i] + k_anchor, base_i_off_r = base_i + k_anchor; (right_i < read_length_2) && (base_i_off_r < ksw_e); right_i++, base_i_off_r++)
								{
#ifdef	CHAR_CP
									if(((read_bit_2[tid][right_i >> 5] >> ((31 - (right_i  & 0X1f)) << 1)) & 0X3) != ((buffer_ref_seq[base_i_off_r >> 5] >> ((31 - (base_i_off_r & 0X1f)) << 1)) & 0X3))
										break;
#else
									if(sam_seq2[right_i] != Dna5Tochar[((buffer_ref_seq[base_i_off_r >> 5] >> ((31 - (base_i_off_r & 0X1f)) << 1)) & 0X3)])
										break;
#endif
								}


								seed_length = right_i - left_i - 1;
								if(seed_length > max_seed_length)
								{
									max_seed_length = seed_length;
									max_right = base_i_off_r;
								}

								anchor_seed_buffer[tid][buffer_i].read_left_off = left_i;
								anchor_seed_buffer[tid][buffer_i].read_right_off = right_i;
								anchor_seed_buffer[tid][buffer_i].ref_left_off = base_i_off_l;
								anchor_seed_buffer[tid][buffer_i].ref_right_off = base_i_off_r;
								anchor_seed_buffer[tid][buffer_i].seed_length = seed_length;
								buffer_i++;
							}

							break;
						}
					}

					r_b_v = max_right;
				}

				if((buffer_i > 0) && (max_seed_length > anchor_seed_length_thr))
				{
					qsort(anchor_seed_buffer[tid], buffer_i, sizeof(anchor_seed), comepare_anchor_seed);

					//LV

					left_i = anchor_seed_buffer[tid][0].read_left_off;
					right_i = anchor_seed_buffer[tid][0].read_right_off;
					base_i_off_l = anchor_seed_buffer[tid][0].ref_left_off;
					base_i_off_r = anchor_seed_buffer[tid][0].ref_right_off;

					d_n1 = 0;
					i_n1 = 0;
					if((left_i == -1) && (right_i == read_length_2))
					{
						strcpy(cigar_p2, cigar_m_2);
					}
					else     //indel
					{
#ifdef	LV_CCIGAR

#ifdef	CHAR_CP
						for(bit_char_i = left_i, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
							read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

						for(bit_char_i = base_i_off_l, read_b_i = 0; bit_char_i > base_i_off_l - left_i - 33; bit_char_i--, read_b_i++)
							ali_ref_seq[tid][read_b_i] = ((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
						for(bit_char_i = left_i, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
							read_char[tid][read_b_i] = sam_seq2[bit_char_i];

						for(bit_char_i = base_i_off_l, read_b_i = 0; bit_char_i > base_i_off_l - left_i - 33; bit_char_i--, read_b_i++)
							ali_ref_seq[tid][read_b_i] = Dna5Tochar[((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
						dm1 = computeEditDistanceWithCigar_s_left(ali_ref_seq[tid], left_i + 33, read_char[tid], left_i + 1, lv_k_2, cigarBuf1, f_cigarn, L[tid], &s_offset2);//, 0
#ifdef	CHAR_CP
						for(bit_char_i = right_i, read_b_i = 0; bit_char_i < read_length_2; bit_char_i++, read_b_i++)
							read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

						for(bit_char_i = base_i_off_r, read_b_i = 0; bit_char_i < base_i_off_r + read_length_2 - right_i + 32; bit_char_i++, read_b_i++)
							ali_ref_seq[tid][read_b_i] = ((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
						for(bit_char_i = right_i, read_b_i = 0; bit_char_i < read_length_2; bit_char_i++, read_b_i++)
							read_char[tid][read_b_i] = sam_seq2[bit_char_i];

						for(bit_char_i = base_i_off_r, read_b_i = 0; bit_char_i < base_i_off_r + read_length_2 - right_i + 32; bit_char_i++, read_b_i++)
							ali_ref_seq[tid][read_b_i] = Dna5Tochar[((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
						dm2 = computeEditDistanceWithCigar_s(ali_ref_seq[tid], read_length_2 - right_i + 32, read_char[tid], read_length_2 - right_i, lv_k_2, cigarBuf2, f_cigarn, L[tid]);//, 0
#endif

						//deal with cigar of LV
						//deal with front and back lv cigar
						m_m_n = anchor_seed_buffer[tid][0].seed_length;

						strncpy(str_o, cigarBuf1, f_cigarn);
						s_o = 0;
						f_c = 0;

						pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

						while (pch != NULL)
						{
							pchl = strlen(pch);
							f_cigar[f_cigarn - f_c - 2] = atoi(pch);
							s_o += (pchl + 1);
							f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];
							f_c += 2;

							if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
							if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

							pch = strtok_r(NULL, "DMIS", &saveptr);
						}

						strncpy(b_cigar, cigarBuf2, f_cigarn);
						pch = strtok (cigarBuf2,"DMIS");

						if(pch != NULL)
							pchl = strlen(pch);

						snt = 0;
						if((left_i != -1) && (right_i != read_length_2))
						{
							if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
							{
								f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
								snt += sn;
							}
							else if(f_cigar[f_cigarn - 1] == 'M')
							{
								f_cigar[f_cigarn - 2] += m_m_n;

								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
								snt += sn;
							}
							else if(b_cigar[pchl] == 'M')
							{
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}

								sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
								snt += sn;
							}
							else
							{
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
								snt += sn;
							}
						}
						else if(left_i == -1)
						{
							if(b_cigar[pchl] == 'M')
								sn = sprintf(cigar_p2, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							else	sn = sprintf(cigar_p2, "%uM%s", m_m_n, b_cigar);

							snt += sn;
						}
						else
						{
							if(f_cigar[f_cigarn - 1] == 'M')
							{
								f_cigar[f_cigarn - 2] += m_m_n;
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
							}
							else
							{
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
								snt += sn;
							}
						}
#ifdef	CIGAR_LEN_ERR
						cigar_len = 0;
						s_o_tmp = 0;
						strncpy(cigar_tmp, cigar_p2, snt);
						cigar_tmp[snt] = '\0';
						pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

						while (pch_tmp != NULL)
						{
							pchl_tmp = strlen(pch_tmp);
							s_o_tmp += (pchl_tmp + 1);

							if(cigar_p2[s_o_tmp - 1] != 'D')
							{
								cigar_len_tmp = atoi(pch_tmp);
								cigar_len += cigar_len_tmp;
							}

							pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
						}

						if(read_length_2 != cigar_len)
						{
							if(read_length_2 < cigar_len)
							{
								cigar_len_re = cigar_len_tmp - (cigar_len - read_length_2);
								if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
								else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
								else	strcpy(cigar_p2, cigar_m_2);
							}
							else
							{
								cigar_len_re = cigar_len_tmp + (read_length_2 - cigar_len);
								sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
							}
						}
#endif

					}
#ifdef	NO_S_OFF
					s_offset2 = 0;
#endif
					sam_pos2 = base_i_off_l - left_i + i_n1 - d_n1 - chr_end_n[chr_re - 1] + 2 + s_offset2;

					ksw_re = 1;
					lv_re2f = dm1 + dm2;
				}
				else
				{
#ifdef	REDUCE_ANCHOR

					if((op_rc[tid][v_cnt_i] == 0) || (op_rc[tid][v_cnt_i] == 3))
					{
						other_end_flag = 0;
						while(tra2_i < anchor_n2)
						{
							rcs = rcs2[tid][tra2_i];
							rs = rs2[tid][tra2_i];
							ls = ls2[tid][tra2_i];
							sam_pos2 = poses2[tid][tra2_i];
							lv_re2f = dms2[tid][tra2_i];
							if(orders2[tid][tra2_i] >= MAX_REDUCE_ANCHOR_NUM)
							{
								v_cnt_i_tmp = orders2[tid][tra2_i] - MAX_REDUCE_ANCHOR_NUM;
								op_vector_seq1_tmp = ops_vector_seq1[tid][v_cnt_i_tmp];
								if(ops_mask[v_cnt_i_tmp] == 0)
								{
									ops_mask[v_cnt_i_tmp] = 1;
									other_end_flag = 1;
									tra2_i++;
									break;
								}
							}
							else
							{
								v_cnt_i_tmp = orders2[tid][tra2_i];
								op_vector_seq1_tmp = op_vector_seq1[tid][v_cnt_i_tmp];
								if(op_mask[v_cnt_i_tmp] == 0)
								{
									op_mask[v_cnt_i_tmp] = 1;
									other_end_flag = 1;
									tra2_i++;
									break;
								}
							}
							tra2_i++;
						}
					}
					else
					{
						other_end_flag = 0;
						while(tra1_i < anchor_n1)
						{
							rcs = rcs1[tid][tra1_i];
							rs = rs1[tid][tra1_i];
							ls = ls1[tid][tra1_i];
							sam_pos2 = poses1[tid][tra1_i];
							lv_re2f = dms1[tid][tra1_i];
							if(orders1[tid][tra1_i] >= MAX_REDUCE_ANCHOR_NUM)
							{
								v_cnt_i_tmp = orders1[tid][tra1_i] - MAX_REDUCE_ANCHOR_NUM;
								op_vector_seq1_tmp = ops_vector_seq1[tid][v_cnt_i_tmp];
								if(ops_mask[v_cnt_i_tmp] == 0)
								{
									ops_mask[v_cnt_i_tmp] = 1;
									tra1_i++;
									other_end_flag = 1;
									break;
								}
							}
							else
							{
								v_cnt_i_tmp = orders1[tid][tra1_i];
								op_vector_seq1_tmp = op_vector_seq1[tid][v_cnt_i_tmp];
								if(op_mask[v_cnt_i_tmp] == 0)
								{
									op_mask[v_cnt_i_tmp] = 1;
									tra1_i++;
									other_end_flag = 1;
									break;
								}
							}
							tra1_i++;
						}
					}

					if(other_end_flag)
					{
						x = sam_pos2;
						low = 0;
						high = chr_file_n - 1;

						while ( low <= high )
						{
							mid = (low + high) >> 1;
							if(x < (chr_end_n[mid]))
							{
								high = mid - 1;
							}
							else if(x > (chr_end_n[mid]))
							{
								low = mid + 1;
							}
							else
							{
								chr_re =  mid;
								break;
							}
							chr_re = low;
						}

						sam_pos2 = x - chr_end_n[chr_re - 1] + 1;

						if(rcs == 0)
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit1[tid][0];
							read_bit_2[tid] = read_bit2[tid][1];
#else
							strcpy(sam_seq1, seqio[seqi].read_seq1);

							/*
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
							*/
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv1[tid][0];
							qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif

							cigar_m_1 = cigar_m1[tid];
							cigar_m_2 = cigar_m2[tid];
							lv_k_1 = lv_k1;
							lv_k_2 = lv_k2;
							read_length_1 = read_length1;
							read_length_2 = read_length2;

							pound_pos_1_f = pound_pos1_f_forward;
							pound_pos_1_r = pound_pos1_r_forward;
							pound_pos_2_f = pound_pos2_f_reverse;
							pound_pos_2_r = pound_pos2_r_reverse;
						}
						else if(rcs == 1)
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit2[tid][1];
							read_bit_2[tid] = read_bit1[tid][0];
#else
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

							sam_seq1[sam_seq_i] = '\0';

							strrev1(sam_seq1);
							//strcpy(sam_seq2, seqio[seqi].read_seq1);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv2[tid][1];
							qual_filt_lv_1_o = qual_filt_lv2[tid][0];
#endif

							cigar_m_1 = cigar_m2[tid];
							cigar_m_2 = cigar_m1[tid];
							lv_k_1 = lv_k2;
							lv_k_2 = lv_k1;
							read_length_1 = read_length2;
							read_length_2 = read_length1;

							pound_pos_1_f = pound_pos2_f_reverse;
							pound_pos_1_r = pound_pos2_r_reverse;
							pound_pos_2_f = pound_pos1_f_forward;
							pound_pos_2_r = pound_pos1_r_forward;
						}
						else if(rcs == 2)
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit2[tid][0];
							read_bit_2[tid] = read_bit1[tid][1];
#else
							strcpy(sam_seq1, seqio[seqi].read_seq2);
							/*
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
							*/
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv2[tid][0];
							qual_filt_lv_1_o = qual_filt_lv2[tid][1];
#endif

							cigar_m_1 = cigar_m2[tid];
							cigar_m_2 = cigar_m1[tid];
							lv_k_1 = lv_k2;
							lv_k_2 = lv_k1;
							read_length_1 = read_length2;
							read_length_2 = read_length1;

							pound_pos_1_f = pound_pos2_f_forward;
							pound_pos_1_r = pound_pos2_r_forward;
							pound_pos_2_f = pound_pos1_f_reverse;
							pound_pos_2_r = pound_pos1_r_reverse;
						}
						else
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit1[tid][1];
							read_bit_2[tid] = read_bit2[tid][0];
#else
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

							sam_seq1[sam_seq_i] = '\0';

							strrev1(sam_seq1);
							//strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv1[tid][1];
							qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif

							cigar_m_1 = cigar_m1[tid];
							cigar_m_2 = cigar_m2[tid];
							lv_k_1 = lv_k1;
							lv_k_2 = lv_k2;
							read_length_1 = read_length1;
							read_length_2 = read_length2;

							pound_pos_1_f = pound_pos1_f_reverse;
							pound_pos_1_r = pound_pos1_r_reverse;
							pound_pos_2_f = pound_pos2_f_forward;
							pound_pos_2_r = pound_pos2_r_forward;
						}

						d_n1 = 0;
						i_n1 = 0;
						s_offset1 = 0;
						s_offset2 = 0;
						s_r_o_l = ls;
						s_r_o_r = rs;

						if((s_r_o_l == 0) && (s_r_o_r == 0))
						{
							strcpy(cigar_p2, cigar_m_2);
						}
						else     //indel
						{
#ifdef	OUTPUT_DEBUG
							if(pound_pos_1_f >= s_r_o_r)   //1
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = pound_pos_1_f;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = read_length_1 - pound_pos_1_f;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if(pound_pos_1_r <= s_r_o_l + 1)     //5
							{
								lv_up_left = pound_pos_1_r;//
								lv_up_right = s_r_o_l;
								lv_down_right = read_length_1;
								lv_down_left = s_r_o_r;
								m_n_f = pound_pos_1_r;
								m_n_b = 0;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
							{
								lv_up_left = 0;
								lv_up_right = pound_pos_1_f - 1;
								lv_down_right = read_length_1;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = pound_pos_1_r - pound_pos_1_f;
							}
							else if((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = read_length_1;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = read_length_1 - s_r_o_l - 1;
							}
							else     //4
							{
								lv_up_left = 0;
								lv_up_right = -1;
								lv_down_right = read_length_1;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = s_r_o_r;
							}
#ifdef	QUAL_FILT_SINGLE_OUT
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length_1 - 1- lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]

							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#ifdef	CHAR_CP
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]

#endif
#endif

							//deal with front and back lv cigar
							strncpy(str_o, cigarBuf1, f_cigarn);
							s_o = 0;
							f_c = 0;
							pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

							while (pch != NULL)
							{
								pchl = strlen(pch);
								f_cigar[f_cigarn - f_c - 2] = atoi(pch);
								s_o += (pchl + 1);
								f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

								f_c += 2;

								if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
								if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

								pch = strtok_r(NULL, "DMIS", &saveptr);
							}

							strncpy(b_cigar, cigarBuf2, f_cigarn);
							pch = strtok(cigarBuf2,"DMIS");

							if(pch != NULL)
								pchl = strlen(pch);

							snt = 0;
#ifdef	CIGAR_S_MODIFY
							if(lv_up_left)
							{
								sn = sprintf(cigar_p2 + snt, "%uS", lv_up_left);
								snt += sn;
							}
#else
							if(m_n_f)
							{
								if(f_c)
								{
									if(f_cigar[f_cigarn + 1 - f_c] == 'M')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else
									{
										sn = sprintf(cigar_p2 + snt, "%uM", m_n_f);
										snt += sn;
									}
								}
								else	m_m_n += m_n_f;
							}
#endif
							if((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))   //(op_dm_l1[tid][v_cnt_i] != -1) && (op_dm_r1[tid][v_cnt_i] != read_length1)
							{
								if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
								{
									f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
									snt += sn;
								}
								else if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
									snt += sn;
								}
								else if(b_cigar[pchl] == 'M')
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}

									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
							{
								if(b_cigar[pchl] == 'M')
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
							{
								if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
									snt += sn;
								}
							}
							else
							{
								sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
								snt += sn;
							}
#ifdef	CIGAR_S_MODIFY
							if(lv_down_right < read_length_1)
							{
								sn = sprintf(cigar_p2 + snt, "%uS", read_length_1 - lv_down_right);
								snt += sn;
							}
#else
							if(m_n_b)
							{
								if(cigar_p2[snt - 1] == 'M')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uM", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else if(cigar_p2[snt - 1] == 'S')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uS", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM", m_n_b);
									snt += sn;
								}
							}
#endif

#ifdef	CIGAR_LEN_ERR
							cigar_len = 0;
							s_o_tmp = 0;
							strncpy(cigar_tmp, cigar_p2, snt);
							cigar_tmp[snt] = '\0';
							pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

							while (pch_tmp != NULL)
							{
								pchl_tmp = strlen(pch_tmp);
								s_o_tmp += (pchl_tmp + 1);

								if(cigar_p2[s_o_tmp - 1] != 'D')
								{
									cigar_len_tmp = atoi(pch_tmp);
									cigar_len += cigar_len_tmp;
								}
								pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
							}

							if(read_length_1 != cigar_len)
							{
								if(read_length_1 < cigar_len)
								{
									cigar_len_re = cigar_len_tmp - (cigar_len - read_length_1);
									if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
									else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
									else	strcpy(cigar_p2, cigar_m_1);
								}
								else
								{
									cigar_len_re = cigar_len_tmp + (read_length_1 - cigar_len);
									sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
								}
							}
#endif
						}
#ifdef	NO_S_OFF
						s_offset1 = 0;
#endif
						sam_pos2 = sam_pos2 + i_n1 - d_n1 + s_offset1;

						ksw_re = 1;
					}
					else
					{
						sam_pos2 = sam_pos1;
#ifdef	PICARD_BUG
						strcpy(cigar_p2, cigar_p1);
#else
						strcpy(cigar_p2, "*");
#endif
						ksw_re = 0;
						lv_re2f = 0;

					}
#else
					sam_pos2 = sam_pos1;
#ifdef	PICARD_BUG
					strcpy(cigar_p2, cigar_p1);
#else
					strcpy(cigar_p2, "*");
#endif
					ksw_re = 0;
					lv_re2f = 0;
#endif
				}
#endif

				if(op_rc[tid][v_cnt_i] == 0)
				{
					cp1 = cigar_p1;
					cp2 = cigar_p2;

					xa_d1s[tid][xa_i_1] = '+';
#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 2)	xa_d2s[tid][xa_i_2] = '+';
						else	xa_d2s[tid][xa_i_2] = '-';
					}
					else	xa_d2s[tid][xa_i_2] = '-';
#else
					xa_d2s[tid][xa_i_2] = '-';
#endif
					lv_re1b = lv_re1f;
					lv_re2b = lv_re2f;
				}
				else if(op_rc[tid][v_cnt_i] == 1)
				{
					if(ksw_re == 1)
					{
						sam_pos1 = sam_pos1 ^ sam_pos2;
						sam_pos2 = sam_pos1 ^ sam_pos2;
						sam_pos1 = sam_pos1 ^ sam_pos2;
					}

					cp1 = cigar_p2;
					cp2 = cigar_p1;

#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 0)	xa_d1s[tid][xa_i_1] = '+';
						else	xa_d1s[tid][xa_i_1] = '-';
					}
					else	xa_d1s[tid][xa_i_1] = '+';
#else
					xa_d1s[tid][xa_i_1] = '+';
#endif
					xa_d2s[tid][xa_i_2] = '-';

					lv_re1b = lv_re2f;
					lv_re2b = lv_re1f;
				}
				else if(op_rc[tid][v_cnt_i] == 2)
				{
					if(ksw_re == 1)
					{
						sam_pos1 = sam_pos1 ^ sam_pos2;
						sam_pos2 = sam_pos1 ^ sam_pos2;
						sam_pos1 = sam_pos1 ^ sam_pos2;
					}

					cp1 = cigar_p2;
					cp2 = cigar_p1;

#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 0)	xa_d1s[tid][xa_i_1] = '+';
						else	xa_d1s[tid][xa_i_1] = '-';
					}
					else	xa_d1s[tid][xa_i_1] = '-';
#else
					xa_d1s[tid][xa_i_1] = '-';
#endif
					xa_d2s[tid][xa_i_2] = '+';

					lv_re1b = lv_re2f;
					lv_re2b = lv_re1f;

				}
				else
				{
					cp1 = cigar_p1;
					cp2 = cigar_p2;

					xa_d1s[tid][xa_i_1] = '-';
#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 2)	xa_d2s[tid][xa_i_2] = '+';
						else	xa_d2s[tid][xa_i_2] = '-';
					}
					else	xa_d2s[tid][xa_i_2] = '+';
#else
					xa_d2s[tid][xa_i_2] = '+';
#endif
					lv_re1b = lv_re1f;
					lv_re2b = lv_re2f;
				}

				if(sam_pos1 <= 0) sam_pos1 = 1;
				if(sam_pos2 <= 0) sam_pos2 = 1;

				/*
				if((sam_flag1 == 117) || (sam_flag1 == 69))
				{
				#ifdef	PICARD_BUG
					//strcpy(cigar_p1s[tid][xa_i], cp2);
				#else
					//strcpy(cigar_p1s[tid][xa_i], "*");
				#endif
					strcpy(cigar_p2s[tid][xa_i_2], cp2);
					lv_re2s[tid][xa_i_2] = lv_re2b;
					chr_res[tid][xa_i_2] = chr_re;
					sam_pos2s[tid][xa_i_2] = (uint32_t )sam_pos2;
					xa_i_2++;
				}
				else if((sam_flag1 == 73) || (sam_flag1 == 89))
				{
					strcpy(cigar_p1s[tid][xa_i_1], cp1);
					lv_re1s[tid][xa_i_1] = lv_re1b;
					chr_res[tid][xa_i_1] = chr_re;
					sam_pos1s[tid][xa_i_1] = (uint32_t )sam_pos1;
					xa_i_1++;
				#ifdef	PICARD_BUG
					//strcpy(cigar_p2s[tid][xa_i], cp1);
				#else
					//strcpy(cigar_p2s[tid][xa_i], "*");
				#endif
				}
				else
				{
					strcpy(cigar_p1s[tid][xa_i_1], cp1);
					strcpy(cigar_p2s[tid][xa_i_2], cp2);

					lv_re1s[tid][xa_i_1] = lv_re1b;
					lv_re2s[tid][xa_i_2] = lv_re2b;

					chr_res[tid][xa_i_1] = chr_re;

					sam_pos1s[tid][xa_i_1] = (uint32_t )sam_pos1;
					sam_pos2s[tid][xa_i_2] = (uint32_t )sam_pos2;

					xa_i_1++;
					xa_i_2++;
				}
				*/

				strcpy(cigar_p1s[tid][xa_i_1], cp1);
				strcpy(cigar_p2s[tid][xa_i_2], cp2);

				lv_re1s[tid][xa_i_1] = lv_re1b;
				lv_re2s[tid][xa_i_2] = lv_re2b;

				chr_res[tid][xa_i_1] = chr_re;

				sam_pos1s[tid][xa_i_1] = (uint32_t )sam_pos1;
				sam_pos2s[tid][xa_i_2] = (uint32_t )sam_pos2;

				xa_i_1++;
				xa_i_2++;

			}

			lv_re1f = dm_ops[tid];
#ifdef	ALT_ALL

#ifdef	ALL_ALL_SINGLE
			for(v_cnt_i = 0; v_cnt_i < vs_cnt; v_cnt_i++)
#else
			for(v_cnt_i = 0; (v_cnt_i < vs_cnt) && (dm_op[tid] + 3 > dm_ops[tid]); v_cnt_i++)
#endif
#else
			for(v_cnt_i = 0; (v_cnt_i < vs_cnt) && (v_cnt < 2) && (dm_op[tid] + 3 > dm_ops[tid]); v_cnt_i++)
#endif
			{
#ifdef	REDUCE_ANCHOR
				if(ops_mask[v_cnt_i] == 1)	continue;

				ops_mask[v_cnt_i] = 1;
#endif
				x = ops_vector_pos1[tid][v_cnt_i];
				low = 0;
				high = chr_file_n - 1;

				while ( low <= high )
				{
					mid = (low + high) >> 1;
					if(x < (chr_end_n[mid]))
					{
						high = mid - 1;
					}
					else if(x > (chr_end_n[mid]))
					{
						low = mid + 1;
					}
					else
					{
						chr_re =  mid;
						break;
					}
					chr_re = low;
				}
				if(chr_re)
					sam_pos1 = ops_vector_pos1[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;
				else
				{
					sam_pos1 = ops_vector_pos1[tid][v_cnt_i];
					chr_re = 1;
				}
				//chr_res[tid][xa_i] = chr_re;

				if(ops_rc[tid][v_cnt_i] == 0)
				{
#ifdef OUPUT_SINGLE_REPEAT2

#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit1[tid][0];
					read_bit_2[tid] = read_bit2[tid][1];
#else
					strcpy(sam_seq1, seqio[seqi].read_seq1);
					for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
						sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

					sam_seq2[sam_seq_i] = '\0';

					strrev1(sam_seq2);
#endif

#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv1[tid][0];
					qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif

					ksw_s = x + end_dis1[tid] - devi;
					ksw_e = x + insert_dis + devi;

					cigar_m_1 = cigar_m1[tid];
					cigar_m_2 = cigar_m2[tid];
					lv_k_1 = lv_k1;
					lv_k_2 = lv_k2;
					read_length_1 = read_length1;
					read_length_2 = read_length2;

					pound_pos_1_f = pound_pos1_f_forward;
					pound_pos_1_r = pound_pos1_r_forward;
					pound_pos_2_f = pound_pos2_f_reverse;
					pound_pos_2_r = pound_pos2_r_reverse;
				}
				else if(ops_rc[tid][v_cnt_i] == 1)
				{
#ifdef OUPUT_SINGLE_REPEAT2

#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit2[tid][1];
					read_bit_2[tid] = read_bit1[tid][0];
#else
					for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
						sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

					sam_seq1[sam_seq_i] = '\0';

					strrev1(sam_seq1);
					strcpy(sam_seq2, seqio[seqi].read_seq1);
#endif

#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv2[tid][1];
					qual_filt_lv_1_o = qual_filt_lv2[tid][0];
#endif

					ksw_s = x - (end_dis1[tid] + devi);
					ksw_e = x - (end_dis1[tid] - devi) + read_length1;

					cigar_m_1 = cigar_m2[tid];
					cigar_m_2 = cigar_m1[tid];
					lv_k_1 = lv_k2;
					lv_k_2 = lv_k1;
					read_length_1 = read_length2;
					read_length_2 = read_length1;

					pound_pos_1_f = pound_pos2_f_reverse;
					pound_pos_1_r = pound_pos2_r_reverse;
					pound_pos_2_f = pound_pos1_f_forward;
					pound_pos_2_r = pound_pos1_r_forward;
				}
				else if(ops_rc[tid][v_cnt_i] == 2)
				{
#ifdef OUPUT_SINGLE_REPEAT2

#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit2[tid][0];
					read_bit_2[tid] = read_bit1[tid][1];
#else
					strcpy(sam_seq1, seqio[seqi].read_seq2);
					for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
						sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

					sam_seq2[sam_seq_i] = '\0';

					strrev1(sam_seq2);
#endif

#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv2[tid][0];
					qual_filt_lv_1_o = qual_filt_lv2[tid][1];
#endif

					ksw_s = x + (end_dis2[tid] - devi);
					ksw_e = x + insert_dis + devi;

					cigar_m_1 = cigar_m2[tid];
					cigar_m_2 = cigar_m1[tid];
					lv_k_1 = lv_k2;
					lv_k_2 = lv_k1;
					read_length_1 = read_length2;
					read_length_2 = read_length1;

					pound_pos_1_f = pound_pos2_f_forward;
					pound_pos_1_r = pound_pos2_r_forward;
					pound_pos_2_f = pound_pos1_f_reverse;
					pound_pos_2_r = pound_pos1_r_reverse;
				}
				else
				{
#ifdef OUPUT_SINGLE_REPEAT2

#ifdef	CHAR_CP
					read_bit_1[tid] = read_bit1[tid][1];
					read_bit_2[tid] = read_bit2[tid][0];
#else
					for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
						sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

					sam_seq1[sam_seq_i] = '\0';

					strrev1(sam_seq1);
					strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#endif

#ifdef	QUAL_FILT_SINGLE_OUT
					qual_filt_lv_1 = qual_filt_lv1[tid][1];
					qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif

					ksw_s = x - (end_dis2[tid] + devi);
					ksw_e = x - (end_dis2[tid] - devi) + read_length2;

					cigar_m_1 = cigar_m1[tid];
					cigar_m_2 = cigar_m2[tid];
					lv_k_1 = lv_k1;
					lv_k_2 = lv_k2;
					read_length_1 = read_length1;
					read_length_2 = read_length2;

					pound_pos_1_f = pound_pos1_f_reverse;
					pound_pos_1_r = pound_pos1_r_reverse;
					pound_pos_2_f = pound_pos2_f_forward;
					pound_pos_2_r = pound_pos2_r_forward;
				}

#ifdef	ANCHOR_HASH_ALI
				anchor_hash_p = anchor_hash[tid][ops_rc[tid][v_cnt_i]];
				anchor_array_p = anchor_array[tid][ops_rc[tid][v_cnt_i]];
				anchor_point_p = anchor_point[tid][ops_rc[tid][v_cnt_i]];
				anchor_pos_p = anchor_pos[tid][ops_rc[tid][v_cnt_i]];
#endif
				d_n1 = 0;
				i_n1 = 0;
				s_offset1 = 0;
				s_offset2 = 0;
				s_r_o_l = ops_dm_l1[tid][v_cnt_i];
				s_r_o_r = ops_dm_r1[tid][v_cnt_i];

				if((s_r_o_l == 0) && (s_r_o_r == 0))
				{
					strcpy(cigar_p1, cigar_m_1);
				}
				else     //indel
				{
#ifdef	OUTPUT_DEBUG
					if(pound_pos_1_f >= s_r_o_r)   //1
					{
						lv_up_left = 0;
						lv_up_right = s_r_o_l;
						lv_down_right = pound_pos_1_f;
						lv_down_left = s_r_o_r;
						m_n_f = 0;
						m_n_b = read_length_1 - pound_pos_1_f;
						m_m_n = s_r_o_r - s_r_o_l - 1;
					}
					else if(pound_pos_1_r <= s_r_o_l + 1)     //5
					{
						lv_up_left = pound_pos_1_r;//
						lv_up_right = s_r_o_l;
						lv_down_right = read_length_1;
						lv_down_left = s_r_o_r;
						m_n_f = pound_pos_1_r;
						m_n_b = 0;
						m_m_n = s_r_o_r - s_r_o_l - 1;
					}
					else if((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
					{
						lv_up_left = 0;
						lv_up_right = pound_pos_1_f - 1;
						lv_down_right = read_length_1;
						lv_down_left = pound_pos_1_r;
						m_n_f = 0;
						m_n_b = 0;
						m_m_n = pound_pos_1_r - pound_pos_1_f;
					}
					else if((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
					{
						lv_up_left = 0;
						lv_up_right = s_r_o_l;
						lv_down_right = read_length_1;
						lv_down_left = pound_pos_1_r;
						m_n_f = 0;
						m_n_b = 0;
						m_m_n = read_length_1 - s_r_o_l - 1;
					}
					else     //4
					{
						lv_up_left = 0;
						lv_up_right = -1;
						lv_down_right = read_length_1;
						lv_down_left = s_r_o_r;
						m_n_f = 0;
						m_n_b = 0;
						m_m_n = s_r_o_r;
					}

#ifdef	QUAL_FILT_SINGLE_OUT
					for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length_1 - 1- lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]

					for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
					for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#ifdef	CHAR_CP
					for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#endif
					//deal with front and back lv cigar
					strncpy(str_o, cigarBuf1, f_cigarn);
					s_o = 0;
					f_c = 0;

					pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

					while (pch != NULL)
					{
						pchl = strlen(pch);
						f_cigar[f_cigarn - f_c - 2] = atoi(pch);
						s_o += (pchl + 1);
						f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

						f_c += 2;

						if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
						if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

						pch = strtok_r(NULL, "DMIS", &saveptr);
					}

					strncpy(b_cigar, cigarBuf2, f_cigarn);
					pch = strtok(cigarBuf2,"DMIS");

					if(pch != NULL)
						pchl = strlen(pch);

					snt = 0;
#ifdef	CIGAR_S_MODIFY
					if(lv_up_left)
					{
						sn = sprintf(cigar_p1 + snt, "%uS", lv_up_left);
						snt += sn;
					}
#else
					if(m_n_f)
					{
						if(f_c)
						{
							if(f_cigar[f_cigarn + 1 - f_c] == 'M')
							{
								f_cigar[f_cigarn - f_c] += m_n_f;
							}
							else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
							{
								f_cigar[f_cigarn - f_c] += m_n_f;
							}
							else
							{
								sn = sprintf(cigar_p1 + snt, "%uM", m_n_f);
								snt += sn;
							}
						}
						else	m_m_n += m_n_f;

					}
#endif
					if((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))   //(op_dm_l1[tid][v_cnt_i] != -1) && (op_dm_r1[tid][v_cnt_i] != read_length1)
					{
						if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
						{
							f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s", b_cigar + pchl + 1);
							snt += sn;
						}
						else if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s",b_cigar);
							snt += sn;
						}
						else if(b_cigar[pchl] == 'M')
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}

							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							snt += sn;
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
							snt += sn;
						}
					}
					else if((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
					{
						if(b_cigar[pchl] == 'M')
						{
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							snt += sn;
						}
						else
						{
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
							snt += sn;
						}
					}
					else if((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
					{
						if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
							snt += sn;
						}
					}
					else
					{
						sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
						snt += sn;
					}
#ifdef	CIGAR_S_MODIFY
					if(lv_down_right < read_length_1)
					{
						sn = sprintf(cigar_p1 + snt, "%uS", read_length_1 - lv_down_right);
						snt += sn;
					}
#else
					if(m_n_b)
					{
						if(cigar_p1[snt - 1] == 'M')
						{
							for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
							{
								if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
								m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
							}
							sn = sprintf(cigar_p1 + bit_char_i + 1, "%uM", m_n_b);
							snt = bit_char_i + 1 +sn;
						}
						else if(cigar_p1[snt - 1] == 'S')
						{
							for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
							{
								if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
								m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
							}
							sn = sprintf(cigar_p1 + bit_char_i + 1, "%uS", m_n_b);
							snt = bit_char_i + 1 + sn;
						}
						else
						{
							sn = sprintf(cigar_p1 + snt, "%uM", m_n_b);
							snt += sn;
						}
					}
#endif

#ifdef	CIGAR_LEN_ERR
					cigar_len = 0;
					s_o_tmp = 0;
					strncpy(cigar_tmp, cigar_p1, snt);
					cigar_tmp[snt] = '\0';
					pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

					while (pch_tmp != NULL)
					{
						pchl_tmp = strlen(pch_tmp);
						s_o_tmp += (pchl_tmp + 1);

						if(cigar_p1[s_o_tmp - 1] != 'D')
						{
							cigar_len_tmp = atoi(pch_tmp);
							cigar_len += cigar_len_tmp;
						}

						pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
					}

					if(read_length_1 != cigar_len)
					{
						if(read_length_1 < cigar_len)
						{
							cigar_len_re = cigar_len_tmp - (cigar_len - read_length_1);
							if(cigar_len_re > 0)	sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
							else if(cigar_len_re == 0)	sprintf(cigar_p1 + snt - sn, "\0");
							else	strcpy(cigar_p1, cigar_m_1);
						}
						else
						{
							cigar_len_re = cigar_len_tmp + (read_length_1 - cigar_len);
							sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
						}
					}
#endif

				}
#ifdef	NO_S_OFF
				s_offset1 = 0;
#endif
				sam_pos1 = sam_pos1 + i_n1 - d_n1 + s_offset1;

				ksw_re = 0;
#ifdef	REDUCE_ANCHOR
				other_end_flag = 0;
#endif

#ifdef	ANCHOR_HASH_ALI

				buffer_i = 0;
				r_b_v = 0;

				for(base_i = ksw_s - 1; base_i < ksw_e - k_anchor; base_i += anchor_seed_d)
				{
					if(base_i + k_anchor - 1 < r_b_v)	continue;

					base_re = (base_i & 0X1f);
					if(base_re <= k_anchor_re)
					{
						anchor_ref = ((buffer_ref_seq[base_i >> 5] >> ((k_anchor_re - base_re) << 1)) & anchor_mask);
					}
					else
					{
						anchor_ref = (((buffer_ref_seq[base_i >> 5] & anchor_mask_boundary[32 - base_re]) << ((base_re - k_anchor_re) << 1)) | (buffer_ref_seq[(base_i >> 5) + 1] >> ((32 + k_anchor_re - base_re) << 1)));
					}

					max_right = 0;
					for(tra_i = 0; tra_i < (anchor_hash_p[anchor_ref >> 4] & 0X1f); tra_i++)
					{
						array_index = (anchor_hash_p[anchor_ref >> 4] >> 5) + tra_i;
						if(anchor_array_p[array_index] == (anchor_ref & 0Xf))
						{

							max_seed_length = 0;
							for(print_i = anchor_point_p[array_index]; print_i < anchor_point_p[array_index + 1]; print_i++)
							{
								//extension on both sides
								for(left_i = anchor_pos_p[print_i] - 1, base_i_off_l = base_i - 1; (left_i >= 0) && (base_i_off_l >= ksw_s - 1); left_i--, base_i_off_l--)
								{
#ifdef	CHAR_CP
									if(((read_bit_2[tid][left_i >> 5] >> ((31 - (left_i  & 0X1f)) << 1)) & 0X3) != ((buffer_ref_seq[base_i_off_l >> 5] >> ((31 - (base_i_off_l & 0X1f)) << 1)) & 0X3))
										break;
#else
									if(sam_seq2[left_i] != Dna5Tochar[((buffer_ref_seq[base_i_off_l >> 5] >> ((31 - (base_i_off_l & 0X1f)) << 1)) & 0X3)])
										break;
#endif
								}

								for(right_i = anchor_pos_p[print_i] + k_anchor, base_i_off_r = base_i + k_anchor; (right_i < read_length_2) && (base_i_off_r < ksw_e); right_i++, base_i_off_r++)
								{
#ifdef	CHAR_CP
									if(((read_bit_2[tid][right_i >> 5] >> ((31 - (right_i  & 0X1f)) << 1)) & 0X3) != ((buffer_ref_seq[base_i_off_r >> 5] >> ((31 - (base_i_off_r & 0X1f)) << 1)) & 0X3))
										break;
#else
									if(sam_seq2[right_i] != Dna5Tochar[((buffer_ref_seq[base_i_off_r >> 5] >> ((31 - (base_i_off_r & 0X1f)) << 1)) & 0X3)])
										break;
#endif
								}

								seed_length = right_i - left_i - 1;
								if(seed_length > max_seed_length)
								{
									max_seed_length = seed_length;
									max_right = base_i_off_r;
								}

								anchor_seed_buffer[tid][buffer_i].read_left_off = left_i;
								anchor_seed_buffer[tid][buffer_i].read_right_off = right_i;
								anchor_seed_buffer[tid][buffer_i].ref_left_off = base_i_off_l;
								anchor_seed_buffer[tid][buffer_i].ref_right_off = base_i_off_r;
								anchor_seed_buffer[tid][buffer_i].seed_length = seed_length;
								buffer_i++;
							}

							break;
						}
					}

					r_b_v = max_right;
				}

				if((buffer_i > 0) && (max_seed_length > anchor_seed_length_thr))
				{
					qsort(anchor_seed_buffer[tid], buffer_i, sizeof(anchor_seed), comepare_anchor_seed);

					//LV
					left_i = anchor_seed_buffer[tid][0].read_left_off;
					right_i = anchor_seed_buffer[tid][0].read_right_off;
					base_i_off_l = anchor_seed_buffer[tid][0].ref_left_off;
					base_i_off_r = anchor_seed_buffer[tid][0].ref_right_off;

					d_n1 = 0;
					i_n1 = 0;
					if((left_i == -1) && (right_i == read_length_2))
					{
						strcpy(cigar_p2, cigar_m_2);
					}
					else     //indel
					{
#ifdef	LV_CCIGAR

#ifdef	CHAR_CP
						for(bit_char_i = left_i, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
							read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

						for(bit_char_i = base_i_off_l, read_b_i = 0; bit_char_i > base_i_off_l - left_i - 33; bit_char_i--, read_b_i++)
							ali_ref_seq[tid][read_b_i] = ((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
						for(bit_char_i = left_i, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
							read_char[tid][read_b_i] = sam_seq2[bit_char_i];

						for(bit_char_i = base_i_off_l, read_b_i = 0; bit_char_i > base_i_off_l - left_i - 33; bit_char_i--, read_b_i++)
							ali_ref_seq[tid][read_b_i] = Dna5Tochar[((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
						dm1 = computeEditDistanceWithCigar_s_left(ali_ref_seq[tid], left_i + 33, read_char[tid], left_i + 1, lv_k_2, cigarBuf1, f_cigarn, L[tid], &s_offset2);//, 0
#ifdef	CHAR_CP
						for(bit_char_i = right_i, read_b_i = 0; bit_char_i < read_length_2; bit_char_i++, read_b_i++)
							read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

						for(bit_char_i = base_i_off_r, read_b_i = 0; bit_char_i < base_i_off_r + read_length_2 - right_i + 32; bit_char_i++, read_b_i++)
							ali_ref_seq[tid][read_b_i] = ((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else

						for(bit_char_i = right_i, read_b_i = 0; bit_char_i < read_length_2; bit_char_i++, read_b_i++)
							read_char[tid][read_b_i] = sam_seq2[bit_char_i];

						for(bit_char_i = base_i_off_r, read_b_i = 0; bit_char_i < base_i_off_r + read_length_2 - right_i + 32; bit_char_i++, read_b_i++)
							ali_ref_seq[tid][read_b_i] = Dna5Tochar[((buffer_ref_seq[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
						dm2 = computeEditDistanceWithCigar_s(ali_ref_seq[tid], read_length_2 - right_i + 32, read_char[tid], read_length_2 - right_i, lv_k_2, cigarBuf2, f_cigarn, L[tid]);//, 0
#endif

						//deal with cigar of LV
						//deal with front and back lv cigar
						m_m_n = anchor_seed_buffer[tid][0].seed_length;

						strncpy(str_o, cigarBuf1, f_cigarn);
						s_o = 0;
						f_c = 0;

						pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

						while (pch != NULL)
						{
							pchl = strlen(pch);
							f_cigar[f_cigarn - f_c - 2] = atoi(pch);
							s_o += (pchl + 1);
							f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];
							f_c += 2;

							if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
							if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

							pch = strtok_r(NULL, "DMIS", &saveptr);
						}

						strncpy(b_cigar, cigarBuf2, f_cigarn);
						pch = strtok (cigarBuf2,"DMIS");

						if(pch != NULL)
							pchl = strlen(pch);

						snt = 0;
						if((left_i != -1) && (right_i != read_length_2))
						{
							if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
							{
								f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
								snt += sn;
							}
							else if(f_cigar[f_cigarn - 1] == 'M')
							{
								f_cigar[f_cigarn - 2] += m_m_n;

								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
								snt += sn;
							}
							else if(b_cigar[pchl] == 'M')
							{
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}

								sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
								snt += sn;
							}
							else
							{
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
								snt += sn;
							}
						}
						else if(left_i == -1)
						{
							if(b_cigar[pchl] == 'M')
								sn = sprintf(cigar_p2, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							else	sn = sprintf(cigar_p2, "%uM%s", m_m_n, b_cigar);

							snt += sn;
						}
						else
						{
							if(f_cigar[f_cigarn - 1] == 'M')
							{
								f_cigar[f_cigarn - 2] += m_m_n;
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
							}
							else
							{
								for(f_i = 0; f_i < f_c; f_i += 2)
								{
									sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
									snt += sn;
								}
								sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
								snt += sn;
							}
						}
#ifdef	CIGAR_LEN_ERR
						cigar_len = 0;
						s_o_tmp = 0;
						strncpy(cigar_tmp, cigar_p2, snt);
						cigar_tmp[snt] = '\0';
						pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

						while (pch_tmp != NULL)
						{
							pchl_tmp = strlen(pch_tmp);
							s_o_tmp += (pchl_tmp + 1);

							if(cigar_p2[s_o_tmp - 1] != 'D')
							{
								cigar_len_tmp = atoi(pch_tmp);
								cigar_len += cigar_len_tmp;
							}

							pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
						}

						if(read_length_2 != cigar_len)
						{
							if(read_length_2 < cigar_len)
							{
								cigar_len_re = cigar_len_tmp - (cigar_len - read_length_2);
								if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
								else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
								else	strcpy(cigar_p2, cigar_m_2);
							}
							else
							{
								cigar_len_re = cigar_len_tmp + (read_length_2 - cigar_len);
								sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
							}
						}
#endif
					}
#ifdef	NO_S_OFF
					s_offset2 = 0;
#endif
					sam_pos2 = base_i_off_l - left_i + i_n1 - d_n1 - chr_end_n[chr_re - 1] + 2 + s_offset2;

					ksw_re = 1;
					lv_re2f = dm1 + dm2;
				}
				else
				{
#ifdef	REDUCE_ANCHOR

					if((ops_rc[tid][v_cnt_i] == 0) || (ops_rc[tid][v_cnt_i] == 3))
					{
						other_end_flag = 0;
						while(tra2_i < anchor_n2)
						{
							rcs = rcs2[tid][tra2_i];
							rs = rs2[tid][tra2_i];
							ls = ls2[tid][tra2_i];
							sam_pos2 = poses2[tid][tra2_i];
							lv_re2f = dms2[tid][tra2_i];
							if(orders2[tid][tra2_i] >= MAX_REDUCE_ANCHOR_NUM)
							{
								v_cnt_i_tmp = orders2[tid][tra2_i] - MAX_REDUCE_ANCHOR_NUM;
								op_vector_seq1_tmp = ops_vector_seq1[tid][v_cnt_i_tmp];
								if(ops_mask[v_cnt_i_tmp] == 0)
								{
									ops_mask[v_cnt_i_tmp] = 1;
									tra2_i++;
									other_end_flag = 1;
									break;
								}
							}
							else
							{
								v_cnt_i_tmp = orders2[tid][tra2_i];
								op_vector_seq1_tmp = op_vector_seq1[tid][v_cnt_i_tmp];
								if(op_mask[v_cnt_i_tmp] == 0)
								{
									op_mask[v_cnt_i_tmp] = 1;
									tra2_i++;
									other_end_flag = 1;
									break;
								}
							}
							tra2_i++;
						}
					}
					else
					{
						other_end_flag = 0;
						while(tra1_i < anchor_n1)
						{
							rcs = rcs1[tid][tra1_i];
							rs = rs1[tid][tra1_i];
							ls = ls1[tid][tra1_i];
							sam_pos2 = poses1[tid][tra1_i];
							lv_re2f = dms1[tid][tra1_i];
							if(orders1[tid][tra1_i] >= MAX_REDUCE_ANCHOR_NUM)
							{
								v_cnt_i_tmp = orders1[tid][tra1_i] - MAX_REDUCE_ANCHOR_NUM;
								op_vector_seq1_tmp = ops_vector_seq1[tid][v_cnt_i_tmp];
								if(ops_mask[v_cnt_i_tmp] == 0)
								{
									ops_mask[v_cnt_i_tmp] = 1;
									tra1_i++;
									other_end_flag = 1;
									break;
								}

							}
							else
							{
								v_cnt_i_tmp = orders1[tid][tra1_i];
								op_vector_seq1_tmp = op_vector_seq1[tid][v_cnt_i_tmp];
								if(op_mask[v_cnt_i_tmp] == 0)
								{
									op_mask[v_cnt_i_tmp] = 1;
									tra1_i++;
									other_end_flag = 1;
									break;
								}
							}
							tra1_i++;
						}
					}

					if(other_end_flag)
					{
						x = sam_pos2;
						low = 0;
						high = chr_file_n - 1;

						while ( low <= high )
						{
							mid = (low + high) >> 1;
							if(x < (chr_end_n[mid]))
							{
								high = mid - 1;
							}
							else if(x > (chr_end_n[mid]))
							{
								low = mid + 1;
							}
							else
							{
								chr_re =  mid;
								break;
							}
							chr_re = low;
						}

						sam_pos2 = x - chr_end_n[chr_re - 1] + 1;

						if(rcs == 0)
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit1[tid][0];
							read_bit_2[tid] = read_bit2[tid][1];
#else
							strcpy(sam_seq1, seqio[seqi].read_seq1);

							/*
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
							*/
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv1[tid][0];
							qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif

							cigar_m_1 = cigar_m1[tid];
							cigar_m_2 = cigar_m2[tid];
							lv_k_1 = lv_k1;
							lv_k_2 = lv_k2;
							read_length_1 = read_length1;
							read_length_2 = read_length2;

							pound_pos_1_f = pound_pos1_f_forward;
							pound_pos_1_r = pound_pos1_r_forward;
							pound_pos_2_f = pound_pos2_f_reverse;
							pound_pos_2_r = pound_pos2_r_reverse;
						}
						else if(rcs == 1)
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit2[tid][1];
							read_bit_2[tid] = read_bit1[tid][0];
#else
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

							sam_seq1[sam_seq_i] = '\0';

							strrev1(sam_seq1);
							//strcpy(sam_seq2, seqio[seqi].read_seq1);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv2[tid][1];
							qual_filt_lv_1_o = qual_filt_lv2[tid][0];
#endif

							cigar_m_1 = cigar_m2[tid];
							cigar_m_2 = cigar_m1[tid];
							lv_k_1 = lv_k2;
							lv_k_2 = lv_k1;
							read_length_1 = read_length2;
							read_length_2 = read_length1;

							pound_pos_1_f = pound_pos2_f_reverse;
							pound_pos_1_r = pound_pos2_r_reverse;
							pound_pos_2_f = pound_pos1_f_forward;
							pound_pos_2_r = pound_pos1_r_forward;
						}
						else if(rcs == 2)
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit2[tid][0];
							read_bit_2[tid] = read_bit1[tid][1];
#else
							strcpy(sam_seq1, seqio[seqi].read_seq2);
							/*
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
							*/
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv2[tid][0];
							qual_filt_lv_1_o = qual_filt_lv2[tid][1];
#endif

							cigar_m_1 = cigar_m2[tid];
							cigar_m_2 = cigar_m1[tid];
							lv_k_1 = lv_k2;
							lv_k_2 = lv_k1;
							read_length_1 = read_length2;
							read_length_2 = read_length1;

							pound_pos_1_f = pound_pos2_f_forward;
							pound_pos_1_r = pound_pos2_r_forward;
							pound_pos_2_f = pound_pos1_f_reverse;
							pound_pos_2_r = pound_pos1_r_reverse;
						}
						else
						{
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit1[tid][1];
							read_bit_2[tid] = read_bit2[tid][0];
#else
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

							sam_seq1[sam_seq_i] = '\0';

							strrev1(sam_seq1);
							//strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv1[tid][1];
							qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif

							cigar_m_1 = cigar_m1[tid];
							cigar_m_2 = cigar_m2[tid];
							lv_k_1 = lv_k1;
							lv_k_2 = lv_k2;
							read_length_1 = read_length1;
							read_length_2 = read_length2;

							pound_pos_1_f = pound_pos1_f_reverse;
							pound_pos_1_r = pound_pos1_r_reverse;
							pound_pos_2_f = pound_pos2_f_forward;
							pound_pos_2_r = pound_pos2_r_forward;
						}

						d_n1 = 0;
						i_n1 = 0;
						s_offset1 = 0;
						s_offset2 = 0;
						s_r_o_l = ls;
						s_r_o_r = rs;

						if((s_r_o_l == 0) && (s_r_o_r == 0))
						{
							strcpy(cigar_p2, cigar_m_2);
						}
						else     //indel
						{
#ifdef	OUTPUT_DEBUG
							if(pound_pos_1_f >= s_r_o_r)   //1
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = pound_pos_1_f;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = read_length_1 - pound_pos_1_f;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if(pound_pos_1_r <= s_r_o_l + 1)     //5
							{
								lv_up_left = pound_pos_1_r;//
								lv_up_right = s_r_o_l;
								lv_down_right = read_length_1;
								lv_down_left = s_r_o_r;
								m_n_f = pound_pos_1_r;
								m_n_b = 0;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
							{
								lv_up_left = 0;
								lv_up_right = pound_pos_1_f - 1;
								lv_down_right = read_length_1;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = pound_pos_1_r - pound_pos_1_f;
							}
							else if((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = read_length_1;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = read_length_1 - s_r_o_l - 1;
							}
							else     //4
							{
								lv_up_left = 0;
								lv_up_right = -1;
								lv_down_right = read_length_1;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = s_r_o_r;
							}
#ifdef	QUAL_FILT_SINGLE_OUT
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length_1 - 1- lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]

							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#ifdef	CHAR_CP
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]

#endif
#endif

							//deal with front and back lv cigar
							strncpy(str_o, cigarBuf1, f_cigarn);
							s_o = 0;
							f_c = 0;
							pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

							while (pch != NULL)
							{
								pchl = strlen(pch);
								f_cigar[f_cigarn - f_c - 2] = atoi(pch);
								s_o += (pchl + 1);
								f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

								f_c += 2;

								if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
								if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

								pch = strtok_r(NULL, "DMIS", &saveptr);
							}

							strncpy(b_cigar, cigarBuf2, f_cigarn);
							pch = strtok(cigarBuf2,"DMIS");

							if(pch != NULL)
								pchl = strlen(pch);

							snt = 0;
#ifdef	CIGAR_S_MODIFY
							if(lv_up_left)
							{
								sn = sprintf(cigar_p2 + snt, "%uS", lv_up_left);
								snt += sn;
							}
#else
							if(m_n_f)
							{
								if(f_c)
								{
									if(f_cigar[f_cigarn + 1 - f_c] == 'M')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else
									{
										sn = sprintf(cigar_p2 + snt, "%uM", m_n_f);
										snt += sn;
									}
								}
								else	m_m_n += m_n_f;
							}
#endif
							if((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))   //(op_dm_l1[tid][v_cnt_i] != -1) && (op_dm_r1[tid][v_cnt_i] != read_length1)
							{
								if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
								{
									f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
									snt += sn;
								}
								else if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
									snt += sn;
								}
								else if(b_cigar[pchl] == 'M')
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}

									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
							{
								if(b_cigar[pchl] == 'M')
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
							{
								if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
									snt += sn;
								}
							}
							else
							{
								sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
								snt += sn;
							}
#ifdef	CIGAR_S_MODIFY
							if(lv_down_right < read_length_1)
							{
								sn = sprintf(cigar_p2 + snt, "%uS", read_length_1 - lv_down_right);
								snt += sn;
							}
#else
							if(m_n_b)
							{
								if(cigar_p2[snt - 1] == 'M')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uM", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else if(cigar_p2[snt - 1] == 'S')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uS", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM", m_n_b);
									snt += sn;
								}
							}
#endif

#ifdef	CIGAR_LEN_ERR
							cigar_len = 0;
							s_o_tmp = 0;
							strncpy(cigar_tmp, cigar_p2, snt);
							cigar_tmp[snt] = '\0';
							pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

							while (pch_tmp != NULL)
							{
								pchl_tmp = strlen(pch_tmp);
								s_o_tmp += (pchl_tmp + 1);

								if(cigar_p2[s_o_tmp - 1] != 'D')
								{
									cigar_len_tmp = atoi(pch_tmp);
									cigar_len += cigar_len_tmp;
								}
								pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
							}

							if(read_length_1 != cigar_len)
							{
								if(read_length_1 < cigar_len)
								{
									cigar_len_re = cigar_len_tmp - (cigar_len - read_length_1);
									if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
									else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
									else	strcpy(cigar_p2, cigar_m_1);
								}
								else
								{
									cigar_len_re = cigar_len_tmp + (read_length_1 - cigar_len);
									sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
								}
							}
#endif
						}
#ifdef	NO_S_OFF
						s_offset1 = 0;
#endif
						sam_pos2 = sam_pos2 + i_n1 - d_n1 + s_offset1;

						ksw_re = 1;
					}
					else
					{
						sam_pos2 = sam_pos1;
#ifdef	PICARD_BUG
						strcpy(cigar_p2, cigar_p1);
#else
						strcpy(cigar_p2, "*");
#endif
						ksw_re = 0;
						lv_re2f = 0;

					}
#else

					sam_pos2 = sam_pos1;
#ifdef	PICARD_BUG
					strcpy(cigar_p2, cigar_p1);
#else
					strcpy(cigar_p2, "*");
#endif
					ksw_re = 0;
					lv_re2f = 0;
#endif
				}

#endif

				if(ops_rc[tid][v_cnt_i] == 0)
				{
					cp1 = cigar_p1;
					cp2 = cigar_p2;

					xa_d1s[tid][xa_i_1] = '+';
#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 2)	xa_d2s[tid][xa_i_2] = '+';
						else	xa_d2s[tid][xa_i_2] = '-';
					}
					else	xa_d2s[tid][xa_i_2] = '-';
#else
					xa_d2s[tid][xa_i_2] = '-';
#endif

					lv_re1b = lv_re1f;
					lv_re2b = lv_re2f;
				}
				else if(ops_rc[tid][v_cnt_i] == 1)
				{
					if(ksw_re == 1)
					{
						sam_pos1 = sam_pos1 ^ sam_pos2;
						sam_pos2 = sam_pos1 ^ sam_pos2;
						sam_pos1 = sam_pos1 ^ sam_pos2;
					}

					cp1 = cigar_p2;
					cp2 = cigar_p1;
#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 0)	xa_d1s[tid][xa_i_1] = '+';
						else	xa_d1s[tid][xa_i_1] = '-';
					}
					else	xa_d1s[tid][xa_i_1] = '+';
#else
					xa_d1s[tid][xa_i_1] = '+';
#endif
					xa_d2s[tid][xa_i_2] = '-';

					lv_re1b = lv_re2f;
					lv_re2b = lv_re1f;
				}
				else if(ops_rc[tid][v_cnt_i] == 2)
				{
					if(ksw_re == 1)
					{
						sam_pos1 = sam_pos1 ^ sam_pos2;
						sam_pos2 = sam_pos1 ^ sam_pos2;
						sam_pos1 = sam_pos1 ^ sam_pos2;
					}

					cp1 = cigar_p2;
					cp2 = cigar_p1;
#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 0)	xa_d1s[tid][xa_i_1] = '+';
						else	xa_d1s[tid][xa_i_1] = '-';
					}
					else	xa_d1s[tid][xa_i_1] = '-';
#else
					xa_d1s[tid][xa_i_1] = '-';
#endif
					xa_d2s[tid][xa_i_2] = '+';

					lv_re1b = lv_re2f;
					lv_re2b = lv_re1f;
				}
				else
				{
					cp1 = cigar_p1;
					cp2 = cigar_p2;

					xa_d1s[tid][xa_i_1] = '-';
#ifdef	REDUCE_ANCHOR
					if(other_end_flag)
					{
						if(rcs == 2)	xa_d2s[tid][xa_i_2] = '+';
						else	xa_d2s[tid][xa_i_2] = '-';
					}
					else	xa_d2s[tid][xa_i_2] = '+';
#else
					xa_d2s[tid][xa_i_2] = '+';
#endif
					lv_re1b = lv_re1f;
					lv_re2b = lv_re2f;
				}

				if(sam_pos1 <= 0) sam_pos1 = 1;
				if(sam_pos2 <= 0) sam_pos2 = 1;

				/*
				if((sam_flag1 == 117) || (sam_flag1 == 69))
				{
				#ifdef	PICARD_BUG
					//strcpy(cigar_p1s[tid][xa_i], cp2);
				#else
					//strcpy(cigar_p1s[tid][xa_i], "*");
				#endif
					strcpy(cigar_p2s[tid][xa_i_2], cp2);
					lv_re2s[tid][xa_i_2] = lv_re2b;
					chr_res[tid][xa_i_2] = chr_re;
					sam_pos2s[tid][xa_i_2] = (uint32_t )sam_pos2;
					xa_i_2++;
				}
				else if((sam_flag1 == 73) || (sam_flag1 == 89))
				{
					strcpy(cigar_p1s[tid][xa_i_1], cp1);
					lv_re1s[tid][xa_i_1] = lv_re1b;
					chr_res[tid][xa_i_1] = chr_re;
					sam_pos1s[tid][xa_i_1] = (uint32_t )sam_pos1;
					xa_i_1++;
				#ifdef	PICARD_BUG
					//strcpy(cigar_p2s[tid][xa_i], cp1);
				#else
					//strcpy(cigar_p2s[tid][xa_i], "*");
				#endif
				}
				else
				{
					strcpy(cigar_p1s[tid][xa_i_1], cp1);
					strcpy(cigar_p2s[tid][xa_i_2], cp2);

					lv_re1s[tid][xa_i_1] = lv_re1b;
					lv_re2s[tid][xa_i_2] = lv_re2b;

					chr_res[tid][xa_i_1] = chr_re;

					sam_pos1s[tid][xa_i_1] = (uint32_t )sam_pos1;
					sam_pos2s[tid][xa_i_2] = (uint32_t )sam_pos2;

					xa_i_1++;
					xa_i_2++;
				}
				*/

				strcpy(cigar_p1s[tid][xa_i_1], cp1);
				strcpy(cigar_p2s[tid][xa_i_2], cp2);

				lv_re1s[tid][xa_i_1] = lv_re1b;
				lv_re2s[tid][xa_i_2] = lv_re2b;

				chr_res[tid][xa_i_1] = chr_re;

				sam_pos1s[tid][xa_i_1] = (uint32_t )sam_pos1;
				sam_pos2s[tid][xa_i_2] = (uint32_t )sam_pos2;

				xa_i_1++;
				xa_i_2++;
			}
#else
			v_cnt = 0;
#endif

			seqio[seqi].v_cnt = v_cnt;
			if(v_cnt > 0)
			{
				if(qual_flag == 1)
					xa_i_1 = 0;

				if(qual_flag == 2)
					xa_i_2 = 0;

				//seqio[seqi].xa_n = xa_i;
				seqio[seqi].xa_n_p1 = xa_i_1;
				seqio[seqi].xa_n_p2 = xa_i_2;
				seqio[seqi].xa_n1 = 0;
				seqio[seqi].xa_n2 = 0;

				if((v_cnt == 1) || ((xa_i_1 == 0) && (xa_i_2 == 0)))
				{
					if(xa_i_1)
						seqio[seqi].qualc1 = 20;
					else
						seqio[seqi].qualc1 = 60;

					if(xa_i_2)
						seqio[seqi].qualc2 = 20;
					else
						seqio[seqi].qualc2 = 60;
				}
				else
				{
					seqio[seqi].qualc1 = 0;
					seqio[seqi].qualc2 = 0;
				}

				if(qual_flag == 1)
					seqio[seqi].qualc1 = 0;

				if(qual_flag == 2)
					seqio[seqi].qualc2 = 0;

				memcpy(chr_res_buffer[seqi], chr_res[tid], (xa_i_1 > xa_i_2 ? xa_i_1:xa_i_2) << 2);
				seqio[seqi].chr_res = chr_res_buffer[seqi];

#ifdef	FIX_SV

				seqio[seqi].xa_n_x1 = 0;
				seqio[seqi].xa_n_x2 = 0;

				tra_i_n = 0;
				if(sv_add == 1)
				{
					for(tra_i = tra1_i; (tra_i_n < 1) && (tra_i < anchor_n1) && ((xa_i_1 + tra_i_n) < CUS_MAX_OUTPUT_ALI2); tra_i++)
					{
						if(orders1[tid][tra_i] >= MAX_REDUCE_ANCHOR_NUM)
						{
							v_cnt_i_tmp = orders1[tid][tra_i] - MAX_REDUCE_ANCHOR_NUM;
							if(ops_mask[v_cnt_i_tmp] == 1)	continue;
							op_vector_seq1_tmp = ops_vector_seq1[tid][v_cnt_i_tmp];
						}
						else
						{
							v_cnt_i_tmp = orders1[tid][tra_i];
							if(op_mask[v_cnt_i_tmp] == 1)	continue;
							op_vector_seq1_tmp = op_vector_seq1[tid][v_cnt_i_tmp];
						}

						sam_pos2 = poses1[tid][tra_i];
						x = sam_pos2;
						low = 0;
						high = chr_file_n - 1;

						while ( low <= high )
						{
							mid = (low + high) >> 1;
							if(x < (chr_end_n[mid]))
							{
								high = mid - 1;
							}
							else if(x > (chr_end_n[mid]))
							{
								low = mid + 1;
							}
							else
							{
								chr_re =  mid;
								break;
							}
							chr_re = low;
						}

						sam_pos2 = x - chr_end_n[chr_re - 1] + 1;

						rcs = rcs1[tid][tra_i];
						rs = rs1[tid][tra_i];
						ls = ls1[tid][tra_i];
						chr_res_buffer1[seqi][tra_i_n] = chr_re;

						if(rcs == 0)
						{
							xa_d1s[tid][xa_i_1 + tra_i_n] = '+';
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit1[tid][0];
							read_bit_2[tid] = read_bit2[tid][1];
#else
							strcpy(sam_seq1, seqio[seqi].read_seq1);

							/*
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
							*/
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv1[tid][0];
							qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif

							//cigar_m_1 = cigar_m1[tid];
							//cigar_m_2 = cigar_m2[tid];
							lv_k_1 = lv_k1;
							lv_k_2 = lv_k2;
							//read_length_1 = read_length1;
							//read_length_2 = read_length2;

							pound_pos_1_f = pound_pos1_f_forward;
							pound_pos_1_r = pound_pos1_r_forward;
							pound_pos_2_f = pound_pos2_f_reverse;
							pound_pos_2_r = pound_pos2_r_reverse;
						}
						else
						{
							xa_d1s[tid][xa_i_1 + tra_i_n] = '-';
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit1[tid][1];
							read_bit_2[tid] = read_bit2[tid][0];
#else
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

							sam_seq1[sam_seq_i] = '\0';

							strrev1(sam_seq1);
							//strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv1[tid][1];
							qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif

							//cigar_m_1 = cigar_m1[tid];
							//cigar_m_2 = cigar_m2[tid];
							lv_k_1 = lv_k1;
							lv_k_2 = lv_k2;
							//read_length_1 = read_length1;
							//read_length_2 = read_length2;

							pound_pos_1_f = pound_pos1_f_reverse;
							pound_pos_1_r = pound_pos1_r_reverse;
							pound_pos_2_f = pound_pos2_f_forward;
							pound_pos_2_r = pound_pos2_r_forward;
						}

						d_n1 = 0;
						i_n1 = 0;
						s_offset1 = 0;
						s_r_o_l = ls;
						s_r_o_r = rs;

						if((s_r_o_l == 0) && (s_r_o_r == 0))
						{
							strcpy(cigar_p2, cigar_m1[tid]);
						}
						else     //indel
						{
#ifdef	OUTPUT_DEBUG
							if(pound_pos_1_f >= s_r_o_r)   //1
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = pound_pos_1_f;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = read_length1 - pound_pos_1_f;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if(pound_pos_1_r <= s_r_o_l + 1)     //5
							{
								lv_up_left = pound_pos_1_r;//
								lv_up_right = s_r_o_l;
								lv_down_right = read_length1;
								lv_down_left = s_r_o_r;
								m_n_f = pound_pos_1_r;
								m_n_b = 0;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
							{
								lv_up_left = 0;
								lv_up_right = pound_pos_1_f - 1;
								lv_down_right = read_length1;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = pound_pos_1_r - pound_pos_1_f;
							}
							else if((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = read_length1;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = read_length1 - s_r_o_l - 1;
							}
							else//4
							{
								lv_up_left = 0;
								lv_up_right = -1;
								lv_down_right = read_length1;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = s_r_o_r;
							}
#ifdef	QUAL_FILT_SINGLE_OUT
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length1 - 1- lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]

							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#ifdef	CHAR_CP
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]

#endif
#endif

							//deal with front and back lv cigar
							strncpy(str_o, cigarBuf1, f_cigarn);
							s_o = 0;
							f_c = 0;
							pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

							while (pch != NULL)
							{
								pchl = strlen(pch);
								f_cigar[f_cigarn - f_c - 2] = atoi(pch);
								s_o += (pchl + 1);
								if(str_o[s_o - 1] == 'S')	f_cigar[f_cigarn - f_c - 1] = 'S';
								else	f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

								f_c += 2;

								if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
								if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

								pch = strtok_r(NULL, "DMIS", &saveptr);
							}

							strncpy(b_cigar, cigarBuf2, f_cigarn);
							pch = strtok(cigarBuf2,"DMIS");

							if(pch != NULL)
								pchl = strlen(pch);

							snt = 0;
#ifdef	CIGAR_S_MODIFY
							if(lv_up_left)
							{
								sn = sprintf(cigar_p2 + snt, "%uS", lv_up_left);
								snt += sn;
							}
#else
							if(m_n_f)
							{
								if(f_c)
								{
									if(f_cigar[f_cigarn + 1 - f_c] == 'M')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else
									{
										sn = sprintf(cigar_p2 + snt, "%uM", m_n_f);
										snt += sn;
									}
								}
								else	m_m_n += m_n_f;
							}
#endif
							if((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))   //(op_dm_l1[tid][v_cnt_i] != -1) && (op_dm_r1[tid][v_cnt_i] != read_length1)
							{
								if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
								{
									f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
									snt += sn;
								}
								else if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
									snt += sn;
								}
								else if(b_cigar[pchl] == 'M')
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}

									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
							{
								if(b_cigar[pchl] == 'M')
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
							{
								if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
									snt += sn;
								}
							}
							else
							{
								sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
								snt += sn;
							}
#ifdef	CIGAR_S_MODIFY
							if(lv_down_right < read_length1)
							{
								sn = sprintf(cigar_p2 + snt, "%uS", read_length1 - lv_down_right);
								snt += sn;
							}
#else
							if(m_n_b)
							{
								if(cigar_p2[snt - 1] == 'M')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uM", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else if(cigar_p2[snt - 1] == 'S')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uS", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM", m_n_b);
									snt += sn;
								}
							}
#endif

#ifdef	CIGAR_LEN_ERR

#ifdef FIX_SA
							sv_s_len = 0;
#endif
							cigar_len = 0;
							s_o_tmp = 0;
							strncpy(cigar_tmp, cigar_p2, snt);
							cigar_tmp[snt] = '\0';
							pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

							while (pch_tmp != NULL)
							{
								pchl_tmp = strlen(pch_tmp);
								s_o_tmp += (pchl_tmp + 1);

								if(cigar_p2[s_o_tmp - 1] != 'D')
								{
									cigar_len_tmp = atoi(pch_tmp);
									cigar_len += cigar_len_tmp;
#ifdef FIX_SA
									if(cigar_p2[s_o_tmp - 1] == 'S')
										sv_s_len += cigar_len_tmp;
#endif
								}
								pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
							}

							if(read_length1 != cigar_len)
							{
								if(read_length1 < cigar_len)
								{
									cigar_len_re = cigar_len_tmp - (cigar_len - read_length1);
									if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
									else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
									else	strcpy(cigar_p2, cigar_m1[tid]);
								}
								else
								{
									cigar_len_re = cigar_len_tmp + (read_length1 - cigar_len);
									sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
								}
							}
#endif
						}
#ifdef	NO_S_OFF
						s_offset1 = 0;
#endif
						sam_pos2 = sam_pos2 + i_n1 - d_n1 + s_offset1;
						if(sam_pos2 == seqio[seqi].pos1)	continue;
						if(sam_pos2 <= 0)	sam_pos2 = 1;

						if(sv_s_len_p > sv_s_len)
						{
							chr_res_buffer1[seqi][tra_i_n] = seqio[seqi].chr_re1;
							seqio[seqi].chr_re1 = chr_re;

							sam_pos1s[tid][xa_i_1 + tra_i_n] = seqio[seqi].pos1;

							if((sam_pos2 + read_length1 - 1) > (chr_end_n[chr_re] - chr_end_n[chr_re - 1]))
								sam_pos2 = chr_end_n[chr_re] - chr_end_n[chr_re - 1] - 1 - read_length1;

							seqio[seqi].pos1 = sam_pos2;
							
							lv_re1s[tid][xa_i_1 + tra_i_n] = seqio[seqi].nm1;

							seqio[seqi].nm1 = dms1[tid][tra_i];
							
							if(seqio[seqi].flag1 & 0X10)
								xa_d1s[tid][xa_i_1 + tra_i_n] = '-';
							else xa_d1s[tid][xa_i_1 + tra_i_n] = '+';

							if(rcs == 0)
							{
								if(op_rc_tmp == 1)	//1+ 2-
								{
									seqio[seqi].flag1 = 97;
									seqio[seqi].flag2 = 145;
									//xa_d1s[tid][xa_i_1 + tra_i_n] = '+';
									if(seqio[seqi].pos2 > sam_pos2)
										sam_cross = seqio[seqi].pos2 + read_length2 - sam_pos2;
									else
										sam_cross = seqio[seqi].pos2 - sam_pos2 - read_length1;
									
									seqio[seqi].seq1 = seqio[seqi].read_seq1;

									for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
										sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
									sam_seq1[sam_seq_i] = '\0';
									strrev1(sam_seq1);

									strcpy(read_rev_buffer[seqi], sam_seq1);
									read_rev_buffer[seqi][read_length2] = '\0';
									seqio[seqi].seq2 = read_rev_buffer[seqi];
								}
								else 				//1+ 2+
								{
									seqio[seqi].flag1 = 65;
									seqio[seqi].flag2 = 129;
									//xa_d1s[tid][xa_i_1 + tra_i_n] = '-';
									sam_cross = seqio[seqi].pos2 - sam_pos2;

									seqio[seqi].seq1 = seqio[seqi].read_seq1;
									seqio[seqi].seq2 = seqio[seqi].read_seq2;
									//strrev1(qual2_buffer[seqi]);
								}
							}
							else
							{
								if(op_rc_tmp == 1)	//1- 2-
								{
									seqio[seqi].flag1 = 113;
									seqio[seqi].flag2 = 177;
									//xa_d1s[tid][xa_i_1 + tra_i_n] = '+';
									sam_cross = seqio[seqi].pos2 - sam_pos2;

									for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
										sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];
									sam_seq1[sam_seq_i] = '\0';
									strrev1(sam_seq1);

									strcpy(read_rev_buffer[seqi], sam_seq1);
									read_rev_buffer[seqi][read_length1] = '\0';
									seqio[seqi].seq1 = read_rev_buffer[seqi];

									for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
										sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
									sam_seq1[sam_seq_i] = '\0';
									strrev1(sam_seq1);

									strcpy(read_rev_buffer_1[seqi], sam_seq1);
									read_rev_buffer_1[seqi][read_length2] = '\0';
									seqio[seqi].seq2 = read_rev_buffer_1[seqi];

									//strrev1(qual1_buffer[seqi]);
								}
								else 				//1- 2+
								{
									seqio[seqi].flag1 = 81;
									seqio[seqi].flag2 = 161;
									//xa_d1s[tid][xa_i_1 + tra_i_n] = '-';
									if(sam_pos2 > seqio[seqi].pos2)
										sam_cross = seqio[seqi].pos2 - read_length1 - sam_pos2;
									else
										sam_cross = seqio[seqi].pos2 + read_length2 - sam_pos2;
									
									for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
										sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];
									sam_seq1[sam_seq_i] = '\0';
									strrev1(sam_seq1);

									strcpy(read_rev_buffer[seqi], sam_seq1);
									read_rev_buffer[seqi][read_length1] = '\0';
									seqio[seqi].seq1 = read_rev_buffer[seqi];

									seqio[seqi].seq2 = seqio[seqi].read_seq2;
								}
							}
							seqio[seqi].cross = sam_cross;

							strcpy(cigar_p1s[tid][xa_i_1 + tra_i_n], pr_cigar1_buffer[seqi]);
							strcpy(pr_cigar1_buffer[seqi], cigar_p2);
						}
						else
						{
							sam_pos1s[tid][xa_i_1 + tra_i_n] = (uint32_t )sam_pos2;

							lv_re1s[tid][xa_i_1 + tra_i_n] = dms1[tid][tra_i];
							strcpy(cigar_p1s[tid][xa_i_1 + tra_i_n], cigar_p2);
						}
						tra_i_n++;
					}

					seqio[seqi].chr_res_s1 = chr_res_buffer1[seqi];
					seqio[seqi].xa_n_x1 = tra_i_n;
					//xa_i_1 += seqio[seqi].xa_n_x1;
					xa_i_1 += tra_i_n;

				}
				else
				{
					for(tra_i = tra2_i; (tra_i_n < 1) && (tra_i < anchor_n2) && ((xa_i_2 + tra_i_n) < CUS_MAX_OUTPUT_ALI2); tra_i++)
					{
						if(orders2[tid][tra_i] >= MAX_REDUCE_ANCHOR_NUM)
						{
							v_cnt_i_tmp = orders2[tid][tra_i] - MAX_REDUCE_ANCHOR_NUM;
							if(ops_mask[v_cnt_i_tmp] == 1)	continue;
							op_vector_seq1_tmp = ops_vector_seq1[tid][v_cnt_i_tmp];
						}
						else
						{
							v_cnt_i_tmp = orders2[tid][tra_i];
							if(op_mask[v_cnt_i_tmp] == 1)	continue;
							op_vector_seq1_tmp = op_vector_seq1[tid][v_cnt_i_tmp];
						}

						sam_pos2 = poses2[tid][tra_i];

						x = sam_pos2;
						low = 0;
						high = chr_file_n - 1;

						while ( low <= high )
						{
							mid = (low + high) >> 1;
							if(x < (chr_end_n[mid]))
							{
								high = mid - 1;
							}
							else if(x > (chr_end_n[mid]))
							{
								low = mid + 1;
							}
							else
							{
								chr_re =  mid;
								break;
							}
							chr_re = low;
						}

						sam_pos2 = x - chr_end_n[chr_re - 1] + 1;

						rcs = rcs2[tid][tra_i];
						rs = rs2[tid][tra_i];
						ls = ls2[tid][tra_i];
						chr_res_buffer2[seqi][tra_i_n] = chr_re;

						if(rcs == 1)
						{
							xa_d2s[tid][xa_i_2 + tra_i_n] = '-';
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit2[tid][1];
							read_bit_2[tid] = read_bit1[tid][0];
#else
							for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
								sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

							sam_seq1[sam_seq_i] = '\0';

							strrev1(sam_seq1);
							//strcpy(sam_seq2, seqio[seqi].read_seq1);
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv2[tid][1];
							qual_filt_lv_1_o = qual_filt_lv2[tid][0];
#endif

							//cigar_m_1 = cigar_m2[tid];
							//cigar_m_2 = cigar_m1[tid];
							lv_k_1 = lv_k2;
							lv_k_2 = lv_k1;
							//read_length_1 = read_length2;
							//read_length_2 = read_length1;

							pound_pos_1_f = pound_pos2_f_reverse;
							pound_pos_1_r = pound_pos2_r_reverse;
							pound_pos_2_f = pound_pos1_f_forward;
							pound_pos_2_r = pound_pos1_r_forward;
						}
						else
						{
							xa_d2s[tid][xa_i_2 + tra_i_n] = '+';
#ifdef	CHAR_CP
							read_bit_1[tid] = read_bit2[tid][0];
							read_bit_2[tid] = read_bit1[tid][1];
#else
							strcpy(sam_seq1, seqio[seqi].read_seq2);
							/*
							for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
								sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

							sam_seq2[sam_seq_i] = '\0';

							strrev1(sam_seq2);
							*/
#endif

#ifdef	QUAL_FILT_SINGLE_OUT
							qual_filt_lv_1 = qual_filt_lv2[tid][0];
							qual_filt_lv_1_o = qual_filt_lv2[tid][1];
#endif

							//cigar_m_1 = cigar_m2[tid];
							//cigar_m_2 = cigar_m1[tid];
							lv_k_1 = lv_k2;
							lv_k_2 = lv_k1;
							//read_length_1 = read_length2;
							//read_length_2 = read_length1;

							pound_pos_1_f = pound_pos2_f_forward;
							pound_pos_1_r = pound_pos2_r_forward;
							pound_pos_2_f = pound_pos1_f_reverse;
							pound_pos_2_r = pound_pos1_r_reverse;
						}

						d_n1 = 0;
						i_n1 = 0;
						s_offset1 = 0;
						s_r_o_l = ls;
						s_r_o_r = rs;

						if((s_r_o_l == 0) && (s_r_o_r == 0))
						{
							strcpy(cigar_p2, cigar_m2[tid]);
						}
						else     //indel
						{
#ifdef	OUTPUT_DEBUG
							if(pound_pos_1_f >= s_r_o_r)   //1
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = pound_pos_1_f;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = read_length2 - pound_pos_1_f;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if(pound_pos_1_r <= s_r_o_l + 1)     //5
							{
								lv_up_left = pound_pos_1_r;//
								lv_up_right = s_r_o_l;
								lv_down_right = read_length2;
								lv_down_left = s_r_o_r;
								m_n_f = pound_pos_1_r;
								m_n_b = 0;
								m_m_n = s_r_o_r - s_r_o_l - 1;
							}
							else if((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
							{
								lv_up_left = 0;
								lv_up_right = pound_pos_1_f - 1;
								lv_down_right = read_length2;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = pound_pos_1_r - pound_pos_1_f;
							}
							else if((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
							{
								lv_up_left = 0;
								lv_up_right = s_r_o_l;
								lv_down_right = read_length2;
								lv_down_left = pound_pos_1_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = read_length2 - s_r_o_l - 1;
							}
							else//4
							{
								lv_up_left = 0;
								lv_up_right = -1;
								lv_down_right = read_length2;
								lv_down_left = s_r_o_r;
								m_n_f = 0;
								m_n_b = 0;
								m_m_n = s_r_o_r;
							}
#ifdef	QUAL_FILT_SINGLE_OUT
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length2 - 1- lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]

							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

							computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
							for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k_1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#ifdef	CHAR_CP
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i  & 0X1f)) << 1)) & 0X3);

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = ((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
							for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
								read_char[tid][read_b_i] = sam_seq1[bit_char_i];

							for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
								ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1_tmp[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
							computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k_1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]

#endif
#endif

							//deal with front and back lv cigar
							strncpy(str_o, cigarBuf1, f_cigarn);
							s_o = 0;
							f_c = 0;
							pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

							while (pch != NULL)
							{
								pchl = strlen(pch);
								f_cigar[f_cigarn - f_c - 2] = atoi(pch);
								s_o += (pchl + 1);
								if(str_o[s_o - 1] == 'S')	f_cigar[f_cigarn - f_c - 1] = 'S';
								else	f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

								f_c += 2;

								if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
								if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

								pch = strtok_r(NULL, "DMIS", &saveptr);
							}

							strncpy(b_cigar, cigarBuf2, f_cigarn);
							pch = strtok(cigarBuf2,"DMIS");

							if(pch != NULL)
								pchl = strlen(pch);

							snt = 0;
#ifdef	CIGAR_S_MODIFY
							if(lv_up_left)
							{
								sn = sprintf(cigar_p2 + snt, "%uS", lv_up_left);
								snt += sn;
							}
#else
							if(m_n_f)
							{
								if(f_c)
								{
									if(f_cigar[f_cigarn + 1 - f_c] == 'M')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
									{
										f_cigar[f_cigarn - f_c] += m_n_f;
									}
									else
									{
										sn = sprintf(cigar_p2 + snt, "%uM", m_n_f);
										snt += sn;
									}
								}
								else	m_m_n += m_n_f;
							}
#endif
							if((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))   //(op_dm_l1[tid][v_cnt_i] != -1) && (op_dm_r1[tid][v_cnt_i] != read_length1)
							{
								if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
								{
									f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
									snt += sn;
								}
								else if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
									snt += sn;
								}
								else if(b_cigar[pchl] == 'M')
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}

									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
							{
								if(b_cigar[pchl] == 'M')
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
									snt += sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
									snt += sn;
								}
							}
							else if((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
							{
								if(f_cigar[f_cigarn - 1] == 'M')
								{
									f_cigar[f_cigarn - 2] += m_m_n;
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
								}
								else
								{
									for(f_i = 0; f_i < f_c; f_i += 2)
									{
										sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
										snt += sn;
									}
									sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
									snt += sn;
								}
							}
							else
							{
								sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
								snt += sn;
							}
#ifdef	CIGAR_S_MODIFY
							if(lv_down_right < read_length2)
							{
								sn = sprintf(cigar_p2 + snt, "%uH", read_length2 - lv_down_right);
								snt += sn;
							}
#else
							if(m_n_b)
							{
								if(cigar_p2[snt - 1] == 'M')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uM", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else if(cigar_p2[snt - 1] == 'S')
								{
									for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
									{
										if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
										m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
									}
									sn = sprintf(cigar_p2 + bit_char_i + 1, "%uS", m_n_b);
									snt = bit_char_i + 1 + sn;
								}
								else
								{
									sn = sprintf(cigar_p2 + snt, "%uM", m_n_b);
									snt += sn;
								}
							}
#endif

#ifdef	CIGAR_LEN_ERR

#ifdef FIX_SA
							sv_s_len = 0;
#endif
							cigar_len = 0;
							s_o_tmp = 0;
							strncpy(cigar_tmp, cigar_p2, snt);
							cigar_tmp[snt] = '\0';
							pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

							while (pch_tmp != NULL)
							{
								pchl_tmp = strlen(pch_tmp);
								s_o_tmp += (pchl_tmp + 1);

								if(cigar_p2[s_o_tmp - 1] != 'D')
								{
									cigar_len_tmp = atoi(pch_tmp);
									cigar_len += cigar_len_tmp;
#ifdef FIX_SA
									if(cigar_p2[s_o_tmp - 1] == 'S')
										sv_s_len += cigar_len_tmp;
#endif
								}
								pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
							}

							if(read_length2 != cigar_len)
							{
								if(read_length2 < cigar_len)
								{
									cigar_len_re = cigar_len_tmp - (cigar_len - read_length2);
									if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
									else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
									else	strcpy(cigar_p2, cigar_m2[tid]);
								}
								else
								{
									cigar_len_re = cigar_len_tmp + (read_length2 - cigar_len);
									sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
								}
							}
#endif
						}
#ifdef	NO_S_OFF
						s_offset1 = 0;
#endif
						sam_pos2 = sam_pos2 + i_n1 - d_n1 + s_offset1;
						if(sam_pos2 == seqio[seqi].pos2)	continue;
						if(sam_pos2 <= 0)	sam_pos2 = 1;

						if(sv_s_len_p > sv_s_len)
						{
							chr_res_buffer2[seqi][tra_i_n] = seqio[seqi].chr_re2;
							seqio[seqi].chr_re2 = chr_re;

							sam_pos2s[tid][xa_i_2 + tra_i_n] = seqio[seqi].pos2;

							if((sam_pos2 + read_length2 - 1) > (chr_end_n[chr_re] - chr_end_n[chr_re - 1]))
								sam_pos2 = chr_end_n[chr_re] - chr_end_n[chr_re - 1] - 1 - read_length2;

							seqio[seqi].pos2 = sam_pos2;
							
							lv_re2s[tid][xa_i_2 + tra_i_n] = seqio[seqi].nm2;
							seqio[seqi].nm2 = dms2[tid][tra_i];
							
							if(seqio[seqi].flag2 & 0X10)
								xa_d2s[tid][xa_i_1 + tra_i_n] = '-';
							else xa_d2s[tid][xa_i_1 + tra_i_n] = '+';


							if(op_rc_tmp == 0)
							{
								if(rcs == 1)		//1+ 2-
								{
									seqio[seqi].flag1 = 97;
									seqio[seqi].flag2 = 145;
									//xa_d2s[tid][xa_i_2 + tra_i_n] = '-';
									if(sam_pos2 > seqio[seqi].pos1)
										sam_cross = sam_pos2 + read_length2 - seqio[seqi].pos1;
									else
										sam_cross = sam_pos2 - seqio[seqi].pos1 - read_length1 ;
									
									seqio[seqi].seq1 = seqio[seqi].read_seq1;

									for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
										sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
									sam_seq1[sam_seq_i] = '\0';
									strrev1(sam_seq1);

									strcpy(read_rev_buffer[seqi], sam_seq1);
									read_rev_buffer[seqi][read_length2] = '\0';
									seqio[seqi].seq2 = read_rev_buffer[seqi];
								}
								else 				//1+ 2+
								{
									seqio[seqi].flag1 = 65;
									seqio[seqi].flag2 = 129;
									//xa_d2s[tid][xa_i_2 + tra_i_n] = '+';
									sam_cross = sam_pos2 - seqio[seqi].pos1;

									seqio[seqi].seq1 = seqio[seqi].read_seq1;
									seqio[seqi].seq2 = seqio[seqi].read_seq2;
								}
							}
							else
							{
								if(rcs == 1)		//1- 2-
								{
									seqio[seqi].flag1 = 113;
									seqio[seqi].flag2 = 177;
									//xa_d2s[tid][xa_i_2 + tra_i_n] = '-';
									sam_cross = sam_pos2 - seqio[seqi].pos1;

									for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
										sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];
									sam_seq1[sam_seq_i] = '\0';
									strrev1(sam_seq1);

									strcpy(read_rev_buffer[seqi], sam_seq1);
									read_rev_buffer[seqi][read_length1] = '\0';
									seqio[seqi].seq1 = read_rev_buffer[seqi];

									for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
										sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];
									sam_seq1[sam_seq_i] = '\0';
									strrev1(sam_seq1);

									strcpy(read_rev_buffer_1[seqi], sam_seq1);
									read_rev_buffer_1[seqi][read_length2] = '\0';
									seqio[seqi].seq2 = read_rev_buffer_1[seqi];
								}
								else 				//1- 2+
								{
									seqio[seqi].flag1 = 81;
									seqio[seqi].flag2 = 161;
									//xa_d2s[tid][xa_i_2 + tra_i_n] = '+';
									if(seqio[seqi].pos1 > sam_pos2)
										sam_cross = sam_pos2 - read_length1 - seqio[seqi].pos1;
									else
										sam_cross = sam_pos2 + read_length2 - seqio[seqi].pos1;
									
									for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
										sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];
									sam_seq1[sam_seq_i] = '\0';
									strrev1(sam_seq1);

									strcpy(read_rev_buffer[seqi], sam_seq1);
									read_rev_buffer[seqi][read_length1] = '\0';
									seqio[seqi].seq1 = read_rev_buffer[seqi];

									seqio[seqi].seq2 = seqio[seqi].read_seq2;
								}
							}
							seqio[seqi].cross = sam_cross;

							strcpy(cigar_p2s[tid][xa_i_2 + tra_i_n], pr_cigar2_buffer[seqi]);
							strcpy(pr_cigar2_buffer[seqi], cigar_p2);
						}
						else
						{
							sam_pos2s[tid][xa_i_2 + tra_i_n] = (uint32_t )sam_pos2;

							lv_re2s[tid][xa_i_2 + tra_i_n] = dms2[tid][tra_i];
							strcpy(cigar_p2s[tid][xa_i_2 + tra_i_n], cigar_p2);

						}
						tra_i_n++;
					}

					seqio[seqi].chr_res_s2 = chr_res_buffer2[seqi];
					seqio[seqi].xa_n_x2 = tra_i_n;

					//xa_i_2 += seqio[seqi].xa_n_x2;
					xa_i_2 += tra_i_n;
				}

				if(xa_i_1 > 0)
				{
					memcpy(xa_d1s_buffer[seqi], xa_d1s[tid], xa_i_1);
					seqio[seqi].xa_d1s = xa_d1s_buffer[seqi];

					memcpy(sam_pos1s_buffer[seqi], sam_pos1s[tid], xa_i_1 << 2);
					seqio[seqi].sam_pos1s = sam_pos1s_buffer[seqi];

					memcpy(lv_re1s_buffer[seqi], lv_re1s[tid], xa_i_1 << 2);
					seqio[seqi].lv_re1s = lv_re1s_buffer[seqi];

					for(v_cnt_i = 0; v_cnt_i < xa_i_1; v_cnt_i++)
					{
						pos_l = strlen(cigar_p1s[tid][v_cnt_i]) + 1;
						seqio[seqi].cigar_p1s[v_cnt_i] = (char* )malloc(pos_l);

						memcpy(seqio[seqi].cigar_p1s[v_cnt_i], cigar_p1s[tid][v_cnt_i], pos_l);
					}
				}

				if(xa_i_2 > 0)
				{
					//memcpy(chr_res_buffer[seqi], chr_res[tid], xa_i_2 << 2);
					//seqio[seqi].chr_res = chr_res_buffer[seqi];

					memcpy(xa_d2s_buffer[seqi], xa_d2s[tid], xa_i_2);
					seqio[seqi].xa_d2s = xa_d2s_buffer[seqi];

					memcpy(sam_pos2s_buffer[seqi], sam_pos2s[tid], xa_i_2 << 2);
					seqio[seqi].sam_pos2s = sam_pos2s_buffer[seqi];

					memcpy(lv_re2s_buffer[seqi], lv_re2s[tid], xa_i_2 << 2);
					seqio[seqi].lv_re2s = lv_re2s_buffer[seqi];

					for(v_cnt_i = 0; v_cnt_i < xa_i_2; v_cnt_i++)
					{
						pos_l = strlen(cigar_p2s[tid][v_cnt_i]) + 1;
						seqio[seqi].cigar_p2s[v_cnt_i] = (char* )malloc(pos_l);

						memcpy(seqio[seqi].cigar_p2s[v_cnt_i], cigar_p2s[tid][v_cnt_i], pos_l);
					}
				}
#else

				if(xa_i_1 > 0)
				{
					memcpy(xa_d1s_buffer[seqi], xa_d1s[tid], xa_i_1);
					seqio[seqi].xa_d1s = xa_d1s_buffer[seqi];

					memcpy(sam_pos1s_buffer[seqi], sam_pos1s[tid], xa_i_1 << 2);
					seqio[seqi].sam_pos1s = sam_pos1s_buffer[seqi];

					memcpy(lv_re1s_buffer[seqi], lv_re1s[tid], xa_i_1 << 2);
					seqio[seqi].lv_re1s = lv_re1s_buffer[seqi];

					for(v_cnt_i = 0; v_cnt_i < xa_i_1; v_cnt_i++)
					{
						pos_l = strlen(cigar_p1s[tid][v_cnt_i]) + 1;
						seqio[seqi].cigar_p1s[v_cnt_i] = (char* )malloc(pos_l);

						memcpy(seqio[seqi].cigar_p1s[v_cnt_i], cigar_p1s[tid][v_cnt_i], pos_l);
					}
				}

				if(xa_i_2 > 0)
				{
					//memcpy(chr_res_buffer[seqi], chr_res[tid], xa_i_2 << 2);
					//seqio[seqi].chr_res = chr_res_buffer[seqi];

					memcpy(xa_d2s_buffer[seqi], xa_d2s[tid], xa_i_2);
					seqio[seqi].xa_d2s = xa_d2s_buffer[seqi];

					memcpy(sam_pos2s_buffer[seqi], sam_pos2s[tid], xa_i_2 << 2);
					seqio[seqi].sam_pos2s = sam_pos2s_buffer[seqi];

					memcpy(lv_re2s_buffer[seqi], lv_re2s[tid], xa_i_2 << 2);
					seqio[seqi].lv_re2s = lv_re2s_buffer[seqi];

					for(v_cnt_i = 0; v_cnt_i < xa_i_2; v_cnt_i++)
					{
						pos_l = strlen(cigar_p2s[tid][v_cnt_i]) + 1;
						seqio[seqi].cigar_p2s[v_cnt_i] = (char* )malloc(pos_l);

						memcpy(seqio[seqi].cigar_p2s[v_cnt_i], cigar_p2s[tid][v_cnt_i], pos_l);
					}
				}
#endif

			}
			else
			{
				//seqio[seqi].xa_n = 0;
				seqio[seqi].xa_n_p1 = 0;
				seqio[seqi].xa_n_p2 = 0;
#ifdef	FIX_SV
				seqio[seqi].xa_n_x1 = 0;
				seqio[seqi].xa_n_x2 = 0;
#endif
				seqio[seqi].xa_n1 = 0;
				seqio[seqi].xa_n2 = 0;
				seqio[seqi].flag1 = 77;
				seqio[seqi].flag2 = 141;
				//seqio[seqi].chr_re = chr_file_n;
				seqio[seqi].chr_re1 = chr_file_n;
				seqio[seqi].chr_re2 = chr_file_n;
				seqio[seqi].pos1 = 0;
				seqio[seqi].pos2 = 0;
				seqio[seqi].qualc1 = 0;
				seqio[seqi].qualc2 = 0;
				strcpy(pr_cigar1_buffer[seqi], "*");
				seqio[seqi].cigar1 = pr_cigar1_buffer[seqi];
				strcpy(pr_cigar2_buffer[seqi], "*");
				seqio[seqi].cigar2 = pr_cigar2_buffer[seqi];

				seqio[seqi].seq1 = seqio[seqi].read_seq1;
				seqio[seqi].seq2 = seqio[seqi].read_seq2;
			}
		}
#else
		seqio[seqi].v_cnt = 0;
#endif

	}

	return 0;
}
#ifdef UNPIPATH_OFF_K20
seed_pa_single* single_seed_reduction_core_single64(seed_pa_single* seedpa, uint64_t* read_bit, uint8_t* read_val, uint64_t** seed_set_pos, uint8_t tid, uint16_t read_length, uint16_t* max_seed_length, uint8_t rc_i)
#else
seed_pa_single* single_seed_reduction_core_single(seed_pa_single* seedpa, uint64_t* read_bit, uint8_t* read_val, uint32_t** seed_set_pos, uint8_t tid, uint16_t read_length, uint16_t* max_seed_length, uint8_t rc_i)
#endif
{

#ifdef UNPIPATH_OFF_K20
	uint64_t kmer_pos_uni = 0;
	uint64_t* pos_tmp_p = NULL;
#else
	uint32_t kmer_pos_uni = 0;
	uint32_t* pos_tmp_p = NULL;
#endif

	uint8_t pop_i = 0;
	uint8_t re_d = 0;
	uint8_t b_t_n_r = 0;

	uint16_t left_i = 1;
	uint16_t right_i = 1;
	uint16_t read_off = 0;
	uint16_t read_pos = 0;
	uint16_t length = 0;
	uint16_t rt_pos = 0;
	uint16_t mem_i = 0;
	uint16_t r_b_v = 0;
	uint16_t su_i = 0;
	uint16_t seed_set_i = 0;
	uint16_t cov_i = 0;

	uint32_t seed_i = 0;
	uint32_t uni_offset_s_l = 0;
	uint32_t uni_offset_s_r = 0;
	uint32_t ref_pos_n = 0;
	uint32_t uni_id = 0;
	uint32_t seed_hash = 0;
	uint32_t seed_kmer = 0;
	uint32_t uni_id_tmp = 0Xffffffff;
	uint32_t spa_i = 0;
	uint32_t pos_n = 0;
	uint32_t set_n = 0;
	uint32_t b_p_i = 0;
	uint32_t id_i = 0;
	uint32_t off_n = 0;
	uint32_t off_i = 0;
	uint32_t set_pos_n = 0;
	uint32_t pos_start_off = 0;
	uint32_t set_i = 0;
	uint32_t set_i_pre = 0;
	uint32_t seed_same_n = 0;
	uint32_t max_sets_n = 0;
	uint32_t interval_n = 0;
	uint32_t bn = 0;
	uint32_t sn = 0;
	int r_p1 = 0;
	int set_c_r = 0;

	uint64_t kmer_bit = 0;
	int64_t posp = 0;
	int64_t reduce_pn = 0;
	int64_t reduce_p = 0;
	int64_t reduce_sn = 0;
	int64_t reduce_s = 0;
	int64_t b_r_s = 0;
	int64_t seed_binary_r = 0;
	int64_t seed_id_r = 0;
	int64_t ref_off1 = 0;

	seed_m* seed_tmp = NULL;

#ifdef	SINGLE_PAIR
	uint16_t cov_num_front = 0;
	uint16_t cov_num_re = 0;
	uint16_t length_front = 0;
	uint16_t length_back = 0;
#endif

	mem_i = 0;
	r_b_v = 0;
	for(read_off = 0; read_off <= read_length - k_t; read_off += seed_l[tid])
	{
		if(read_off + k_t - 1 <= r_b_v)	continue;

		re_d = (read_off & 0X1f);

		b_t_n_r = 32 - re_d;

#ifdef HASH_KMER_READ_J
		if(re_d <= re_b)
		{
			kmer_bit = ((read_bit[read_off >> 5] & bit_tran_re[re_d]) >> ((re_b - re_d) << 1));
		}
		else
		{
			kmer_bit = (((read_bit[read_off >> 5] & bit_tran_re[re_d]) << ((re_d - re_b) << 1)) | (read_bit[(read_off >> 5) + 1] >> ((re_2bt - re_d) << 1)));
		}
#else
		tran_tmp_p = (read_bit[read_off >> 5] & bit_tran_re[re_d]);

		//or use this method to deal with: & bit_tran_re[b_t_n_r]
		kmer_bit = (((read_bit[(read_off >> 5) + 1] >> (b_t_n_r << 1)) & bit_tran_re[b_t_n_r]) | (tran_tmp_p << (re_d << 1)));

		kmer_bit >>= re_bt;
#endif

		seed_kmer = (kmer_bit & bit_tran[k_r]);

		seed_hash = (kmer_bit >> (k_r << 1));

		//find the kmer
#ifdef UNPIPATH_OFF_K20
		seed_binary_r = binsearch_offset64(seed_kmer, buffer_kmer_g, buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], buffer_hash_g[seed_hash]);
#else
		seed_binary_r = binsearch_offset(seed_kmer, buffer_kmer_g, buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], buffer_hash_g[seed_hash]);
#endif
		if(seed_binary_r == -1)
			continue;

		//binary search on unipath offset to get the unipathID
		kmer_pos_uni = buffer_off_g[seed_binary_r];//this kmer's offset on unipath
#ifdef UNPIPATH_OFF_K20
		seed_id_r = binsearch_interval_unipath64(kmer_pos_uni, buffer_seqf, result_seqf);
#else
		seed_id_r = binsearch_interval_unipath(kmer_pos_uni, buffer_seqf, result_seqf);
#endif
		ref_pos_n = buffer_pp[seed_id_r + 1] - buffer_pp[seed_id_r];

		if(ref_pos_n > pos_n_max)
		{
			continue;
		}

		uni_offset_s_l = kmer_pos_uni - buffer_seqf[seed_id_r];
		uni_offset_s_r = buffer_seqf[seed_id_r + 1] - (kmer_pos_uni + k_t);

		for(left_i = 1; (left_i <= uni_offset_s_l) && (left_i <= read_off); left_i++)
		{
#ifdef UNI_SEQ64
			if(((buffer_seq[(kmer_pos_uni - left_i) >> 5] >> ((31 - ((kmer_pos_uni - left_i) & 0X1f)) << 1)) & 0X3)
			        != ((read_bit[(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3)
			  )	break;
#else
			if(((buffer_seq[(kmer_pos_uni - left_i) >> 2] >> (((kmer_pos_uni - left_i) & 0X3) << 1)) & 0X3)
			        != ((read_bit[(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3)
			  )	break;
#endif
		}

		for(right_i = 1; (right_i <= uni_offset_s_r) && (right_i <= read_length - read_off - k_t); right_i++)
		{
#ifdef UNI_SEQ64
			if(((buffer_seq[(kmer_pos_uni + k_t - 1 + right_i) >> 5] >> ((31 - ((kmer_pos_uni + k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
			        != ((read_bit[(read_off + k_t - 1 + right_i) >> 5] >> ((31 - ((read_off + k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
			  )	break;
#else
			if(((buffer_seq[(kmer_pos_uni + k_t - 1 + right_i) >> 2] >> (((kmer_pos_uni + k_t - 1 + right_i) & 0X3) << 1)) & 0X3)
			        != ((read_bit[(read_off + k_t - 1 + right_i) >> 5] >> ((31 - ((read_off + k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
			  )	break;
#endif
		}

		read_pos = read_off + 1 - left_i;
		length = k_t + left_i + right_i - 2;
		rt_pos = read_off + k_t + right_i - 2;

		seedm[tid][mem_i].uni_id = seed_id_r;
		seedm[tid][mem_i].length = length;
		seedm[tid][mem_i].read_pos = read_pos;
		seedm[tid][mem_i].ref_pos_off = uni_offset_s_l + 1 - left_i;
		seedm[tid][mem_i].ref_pos_n = ref_pos_n;

		if((left_i - 1 < read_off) || (right_i - 1 < read_length - read_off - k_t))
		{
			seedm[tid][mem_i].ui = 0;

		}
		else
		{
			seedm[tid][mem_i].ui = 1;
			seedm[tid][mem_i].ref_pos_off_r = uni_offset_s_r + 1 - right_i;
		}

		seedm[tid][mem_i].s_r_o_l = read_off - left_i;
		seedm[tid][mem_i].s_r_o_r = read_off + k_t + right_i - 1;

		memset(seedm[tid][mem_i].cov, 0Xff, (cov_a_n_s[tid]) << 3);

		(seedm[tid][mem_i].cov)[read_pos >> 6] = bv_64[(64 - (read_pos & 0X3f))];

		(seedm[tid][mem_i].cov)[rt_pos >> 6] ^= bv_64[(63 - (rt_pos & 0X3f))];

		for(cov_i = (read_pos >> 6); cov_i > 0; cov_i--)
			(seedm[tid][mem_i].cov)[cov_i - 1] = 0;
		for(cov_i = (rt_pos >> 6) + 1; cov_i < cov_a_n_s[tid]; cov_i++)
			(seedm[tid][mem_i].cov)[cov_i] = 0;

		seedm[tid][mem_i].tid = tid;

		r_b_v = read_off + k_t + right_i - 1;

		++mem_i;

	}

	//end MEM

	if(mem_i == 0)	return 0;

	//merge seeds in the same unipath
	s_uid_f[tid] = 0;
	qsort(seedm[tid], mem_i, sizeof(seed_m), compare_uniid);

	if(s_uid_f[tid] == 1)
	{
		su_i = 0;
		for(seed_i = 0; seed_i < mem_i; seed_i++)
		{
			if(uni_id_tmp != seedm[tid][seed_i].uni_id)
			{
				seedu[tid][su_i].uni_id = seedm[tid][seed_i].uni_id;
				seedu[tid][su_i].read_pos = seedm[tid][seed_i].read_pos;
				seedu[tid][su_i].ref_pos_off = seedm[tid][seed_i].ref_pos_off;
				seedu[tid][su_i].ref_pos_off_r = seedm[tid][seed_i].ref_pos_off_r;
				seedu[tid][su_i].ui = seedm[tid][seed_i].ui;
				seedu[tid][su_i].s_r_o_l = seedm[tid][seed_i].s_r_o_l;
				seedu[tid][su_i].s_r_o_r = seedm[tid][seed_i].s_r_o_r;
				seedu[tid][su_i].ref_pos_n = seedm[tid][seed_i].ref_pos_n;
				memcpy(seedu[tid][su_i].cov, seedm[tid][seed_i].cov, (cov_a_n_s[tid]) << 3);

				seedu[tid][su_i].tid = tid;

				//cal length
				if(su_i > 0)
				{
					seedu[tid][su_i - 1].length = 0;
					for(cov_i = 0; cov_i < cov_a_n_s[tid]; cov_i++)
						seedu[tid][su_i - 1].length += popcount_3((seedu[tid][su_i - 1].cov)[cov_i]);
				}
				++su_i;
			}
			else
			{
				for(cov_i = 0; cov_i < cov_a_n_s[tid]; cov_i++)
					(seedu[tid][su_i - 1].cov)[cov_i] |= (seedm[tid][seed_i].cov)[cov_i];
			}

			uni_id_tmp = seedm[tid][seed_i].uni_id;
		}
		seedu[tid][su_i - 1].length = 0;
		for(cov_i = 0; cov_i < cov_a_n_s[tid]; cov_i++)
			seedu[tid][su_i - 1].length += popcount_3((seedu[tid][su_i - 1].cov)[cov_i]);

		//end merge of same unipath
		seed_tmp = seedu[tid];
	}
	else
	{
		su_i = mem_i;

		seed_tmp = seedm[tid];
	}

	//add print
	for ( seed_i = 0; seed_i < su_i; ++seed_i)
	{
		fprintf(stderr, "read_pos = %u, length = %u\n", seedu[tid][seed_i].read_pos, seedu[tid][seed_i].length);
	}

	//uniqueness
	qsort(seed_tmp, su_i, sizeof(seed_m), compare_posn);

	uni_id = seed_tmp[0].uni_id;
	ref_off1 = seed_tmp[0].ref_pos_off;
	r_p1 = seed_tmp[0].read_pos;

	pos_tmp_p = buffer_p + buffer_pp[uni_id];
	pos_n = buffer_pp[uni_id + 1] - buffer_pp[uni_id];

	//fingerprint
	for(id_i = 0; id_i < pos_n; id_i++)
	{
		seedsets[tid][id_i].seed_set = pos_tmp_p[id_i] + ref_off1 - r_p1;
		memcpy(seedsets[tid][id_i].cov, seed_tmp[0].cov, (cov_a_n_s[tid]) << 3);

		seedsets[tid][id_i].ui = seed_tmp[0].ui;
		seedsets[tid][id_i].ref_pos_off = ref_off1;
		seedsets[tid][id_i].ref_pos_off_r = seed_tmp[0].ref_pos_off_r;

		seedsets[tid][id_i].s_r_o_l = seed_tmp[0].s_r_o_l;
		seedsets[tid][id_i].s_r_o_r = seed_tmp[0].s_r_o_r;

		seedsets[tid][id_i].tid = tid;
	}

	set_n += pos_n;
	seed_set_off[tid][0] = 0;
	seed_set_off[tid][1] = pos_n;
	off_n += 2;

#ifdef MERGE_B_R_VS
	for(seed_i = 1; seed_i < su_i; seed_i++)
	{
		uni_id = seed_tmp[seed_i].uni_id;

		ref_off1 = seed_tmp[seed_i].ref_pos_off;
		r_p1 = seed_tmp[seed_i].read_pos;

		reduce_pn = buffer_pp[uni_id + 1];

		bn = reduce_pn - buffer_pp[uni_id];

		for(seed_set_i = 0; seed_set_i < seed_i; seed_set_i++)
		{
			reduce_sn = seed_set_off[tid][seed_set_i + 1];
			sn = reduce_sn - seed_set_off[tid][seed_set_i];

			if(bn <= sn)
			{
				reduce_s = seed_set_off[tid][seed_set_i];

				for(b_p_i = buffer_pp[uni_id]; b_p_i < reduce_pn; b_p_i++)
				{
					posp = buffer_p[b_p_i] + ref_off1 - r_p1;
#ifdef UNPIPATH_OFF_K20
					b_r_s = binsearch_seed_set_reduce64(posp, seedsets[tid], reduce_sn - reduce_s, reduce_s, tid);
#else
					b_r_s = binsearch_seed_set_reduce(posp, seedsets[tid], reduce_sn - reduce_s, reduce_s, tid);
#endif
					if(b_r_s != -1)
					{
						for(cov_i = 0; cov_i < cov_a_n_s[tid]; cov_i++)
							(seedsets[tid][b_r_s].cov)[cov_i] |= (seed_tmp[seed_i].cov)[cov_i];

						pos_add[tid][b_p_i - buffer_pp[uni_id]] = 1;

						if(b_r_s == reduce_sn - 1)
							break;

						reduce_s = b_r_s + 1;
					}
					else
					{
						if(g_low[tid] == reduce_sn - reduce_s)
							break;
						else	reduce_s += g_low[tid];
					}
				}
			}
			else
			{
				reduce_p = buffer_pp[uni_id];

				for(b_p_i = seed_set_off[tid][seed_set_i]; b_p_i < reduce_sn; b_p_i++)
				{
					posp = seedsets[tid][b_p_i].seed_set + r_p1 - ref_off1;
#ifdef UNPIPATH_OFF_K20
					b_r_s = binsearch_seed_pos_reduce64(posp, buffer_p, reduce_pn - reduce_p, reduce_p, tid);
#else
					b_r_s = binsearch_seed_pos_reduce(posp, buffer_p, reduce_pn - reduce_p, reduce_p, tid);
#endif
					if(b_r_s != -1)
					{
						for(cov_i = 0; cov_i < cov_a_n_s[tid]; cov_i++)
							(seedsets[tid][b_p_i].cov)[cov_i] |= (seed_tmp[seed_i].cov)[cov_i];

						pos_add[tid][b_r_s - buffer_pp[uni_id]] = 1;

						if(b_r_s == reduce_pn - 1)
							break;

						reduce_p = b_r_s + 1;

					}
					else
					{
						if(g_low[tid] == reduce_pn - reduce_p)
							break;
						else	reduce_p += g_low[tid];

					}
				}
			}
		}

		for(b_p_i = 0; b_p_i < bn; b_p_i++)
		{
			if(pos_add[tid][b_p_i] == 0)
			{
				seedsets[tid][set_n].seed_set = buffer_p[buffer_pp[uni_id] + b_p_i] + ref_off1 - r_p1;

				memcpy(seedsets[tid][set_n].cov, seed_tmp[seed_i].cov, (cov_a_n_s[tid]) << 3);

				seedsets[tid][set_n].ui = seed_tmp[seed_i].ui;
				seedsets[tid][set_n].ref_pos_off = ref_off1;
				seedsets[tid][set_n].ref_pos_off_r = seed_tmp[seed_i].ref_pos_off_r;

				seedsets[tid][set_n].s_r_o_l = seed_tmp[seed_i].s_r_o_l;
				seedsets[tid][set_n].s_r_o_r = seed_tmp[seed_i].s_r_o_r;

				seedsets[tid][set_n].tid = tid;

				++set_n;
			}
		}

		memset(pos_add[tid], 0, bn);

		seed_set_off[tid][seed_i + 1] = set_n;
		++off_n;
	}
#endif

#ifdef	SINGLE_PAIR
	cov_num_front = cov_num_front_single[tid];
	cov_num_re = cov_num_re_single[tid];
#endif

	spa_i = spa_i_single[tid];
	set_pos_n = set_pos_n_single[tid];
	pos_start_off = set_pos_n_single[tid];
	for(off_i = 0; off_i < off_n - 1; off_i++)
	{
		interval_n = seed_set_off[tid][off_i + 1] - seed_set_off[tid][off_i];

		if(interval_n == 0)	continue;

		qsort(seedsets[tid] + seed_set_off[tid][off_i], interval_n, sizeof(seed_sets), compare_sets_s);

		//copy elements from seed_tmp[off_i] as for other elements such as quality
		set_i_pre = seed_set_off[tid][off_i];
		seedpa[spa_i].rc = rc_i;
		seedpa[spa_i].pos_start = set_i_pre + pos_start_off;
		memcpy(seedpa[spa_i].cov, seedsets[tid][set_i_pre].cov, (cov_a_n_s[tid]) << 3);

		seedpa[spa_i].ref_pos_off = seedsets[tid][set_i_pre].ref_pos_off;
		seedpa[spa_i].ref_pos_off_r = seedsets[tid][set_i_pre].ref_pos_off_r;
		seedpa[spa_i].ui = seedsets[tid][set_i_pre].ui;

		seedpa[spa_i].s_r_o_l = seedsets[tid][set_i_pre].s_r_o_l;
		seedpa[spa_i].s_r_o_r = seedsets[tid][set_i_pre].s_r_o_r;

#ifdef	SINGLE_PAIR
		length_front = 0;
		length_back = 0;
		for(pop_i = 0; pop_i < cov_num_front; pop_i++)
			length_front += popcount_3((seedpa[spa_i].cov)[pop_i]);

		length_front += popcount_3(((seedpa[spa_i].cov)[cov_num_front]) >> cov_num_re);
		length_back += popcount_3(((seedpa[spa_i].cov)[cov_num_front]) & bv_64[cov_num_re]);

		for(pop_i = cov_num_front + 1; pop_i < cov_a_n_s[tid]; pop_i++)
			length_back += popcount_3((seedpa[spa_i].cov)[pop_i]);

		seedpa[spa_i].length = length_front + length_back;
		if(seedpa[spa_i].length > max_sets_n)	max_sets_n = seedpa[spa_i].length;

		if((length_front > k_t) && (length_back > k_t))	seedpa[spa_i].pair_flag = 1;
		else	seedpa[spa_i].pair_flag = 0;
#else
		//cal number of 1 (length)
		seedpa[spa_i].length = 0;
		for(pop_i = 0; pop_i < cov_a_n_s[tid]; pop_i++)
			seedpa[spa_i].length += popcount_3((seedpa[spa_i].cov)[pop_i]);
		if(seedpa[spa_i].length > max_sets_n)	max_sets_n = seedpa[spa_i].length;

#endif
		++spa_i;

		seed_same_n = 1;

		//copy pos to seed_set_pos
		(*seed_set_pos)[set_pos_n++] = seedsets[tid][set_i_pre].seed_set; //+ ref_off - r_p

		for(set_i = set_i_pre + 1; set_i < seed_set_off[tid][off_i + 1]; set_i++)
		{
			set_c_r = memcmp(seedsets[tid][set_i_pre].cov, seedsets[tid][set_i].cov, (cov_a_n_s[tid]) << 3);
			if(set_c_r != 0)
			{
				seedpa[spa_i - 1].pos_n = seed_same_n;

				//copy elements from seed_tmp[off_i] as for other elements such as quality
				seedpa[spa_i].rc = rc_i;
				seedpa[spa_i].pos_start = set_i + pos_start_off;
				memcpy(seedpa[spa_i].cov, seedsets[tid][set_i].cov, (cov_a_n_s[tid]) << 3);

				seedpa[spa_i].ref_pos_off = seedsets[tid][set_i].ref_pos_off;
				seedpa[spa_i].ref_pos_off_r = seedsets[tid][set_i].ref_pos_off_r;
				seedpa[spa_i].ui = seedsets[tid][set_i].ui;

				seedpa[spa_i].s_r_o_l = seedsets[tid][set_i].s_r_o_l;
				seedpa[spa_i].s_r_o_r = seedsets[tid][set_i].s_r_o_r;

#ifdef	SINGLE_PAIR
				length_front = 0;
				length_back = 0;
				for(pop_i = 0; pop_i < cov_num_front; pop_i++)
					length_front += popcount_3((seedpa[spa_i].cov)[pop_i]);

				length_front += popcount_3(((seedpa[spa_i].cov)[cov_num_front]) >> cov_num_re);
				length_back += popcount_3(((seedpa[spa_i].cov)[cov_num_front]) & bv_64[cov_num_re]);

				for(pop_i = cov_num_front + 1; pop_i < cov_a_n_s[tid]; pop_i++)
					length_back += popcount_3((seedpa[spa_i].cov)[pop_i]);

				seedpa[spa_i].length = length_front + length_back;
				if(seedpa[spa_i].length > max_sets_n)	max_sets_n = seedpa[spa_i].length;

				if((length_front > k_t) && (length_back > k_t))	seedpa[spa_i].pair_flag = 1;
				else	seedpa[spa_i].pair_flag = 0;
#else
				//cal number of 1 (length)
				seedpa[spa_i].length = 0;
				for(pop_i = 0; pop_i < cov_a_n_s[tid]; pop_i++)
					seedpa[spa_i].length += popcount_3((seedpa[spa_i].cov)[pop_i]);
				if(seedpa[spa_i].length > max_sets_n)	max_sets_n = seedpa[spa_i].length;
#endif
				++spa_i;

				seed_same_n = 1;
			}
			else	++seed_same_n;

			//copy pos to seed_set_pos
			(*seed_set_pos)[set_pos_n++] = seedsets[tid][set_i].seed_set;// + ref_off - r_p

			set_i_pre = set_i;
		}

		seedpa[spa_i - 1].pos_n = seed_same_n;
	}

	set_pos_n_single[tid] = set_pos_n;

	//qsort(seedpa, spa_i, sizeof(seed_pa), compare_plen);
	spa_i_single[tid] = spa_i;

	(*max_seed_length) = max_sets_n;

	return seedpa;
}


#ifdef UNPIPATH_OFF_K20
seed_pa* single_seed_reduction_core_filter64(seed_pa* seedpa, uint64_t* read_bit, uint8_t* read_val, uint32_t* spa1_i, uint64_t** seed_set_pos, uint16_t* pos_ren, uint8_t tid, uint16_t read_length, uint16_t* seed_num, uint16_t* max_seed_length)
#else
seed_pa* single_seed_reduction_core_filter(seed_pa* seedpa, uint64_t* read_bit, uint8_t* read_val, uint32_t* spa1_i, uint32_t** seed_set_pos, uint16_t* pos_ren, uint8_t tid, uint16_t read_length, uint16_t* seed_num, uint16_t* max_seed_length)
#endif
{
	uint8_t pop_i = 0;
	uint8_t re_d = 0;
	uint8_t b_t_n_r = 0;
	uint8_t cov_a_t = 0;
	uint16_t left_i = 1;
	uint16_t right_i = 1;
	uint16_t read_off = 0;
	uint16_t read_pos = 0;
	uint16_t length = 0;
	uint16_t rt_pos = 0;
	uint16_t mem_i = 0;
	uint16_t r_b_v = 0;
	uint16_t su_i = 0;
	uint16_t seed_set_i = 0;
	uint16_t cov_i = 0;

	uint32_t seed_i = 0;
	uint32_t uni_offset_s_l = 0;
	uint32_t uni_offset_s_r = 0;
	uint32_t ref_pos_n = 0;

	uint32_t seed_hash = 0;
	uint32_t seed_kmer = 0;
	uint32_t uni_id_tmp = 0Xffffffff;
	uint32_t spa_i = 0;
	uint32_t pos_n = 0;
	uint32_t set_n = 0;
	uint32_t id_i = 0;
	uint32_t off_n = 0;
	uint32_t off_i = 0;
	uint32_t set_pos_n = 0;
	uint32_t set_i = 0;
	uint32_t set_i_pre = 0;
	uint32_t seed_same_n = 0;
	uint32_t max_sets_n = 0;
	uint32_t interval_n = 0;
	uint32_t bn = 0;
	uint32_t sn = 0;

#ifdef UNPIPATH_OFF_K20
	uint64_t* pos_tmp_p = NULL;
	uint64_t kmer_pos_uni = 0;
	uint64_t b_p_i = 0;
	uint64_t uni_id = 0;
#else
	uint32_t* pos_tmp_p = NULL;
	uint32_t kmer_pos_uni = 0;
	uint32_t b_p_i = 0;
	uint32_t uni_id = 0;
#endif

	int r_p1 = 0;
	int set_c_r = 0;

	uint64_t kmer_bit = 0;
	int64_t posp = 0;
	int64_t reduce_pn = 0;
	int64_t reduce_p = 0;
	int64_t reduce_sn = 0;
	int64_t reduce_s = 0;
	int64_t b_r_s = 0;
	int64_t seed_binary_r = 0;
	int64_t seed_id_r = 0;
	int64_t ref_off1 = 0;

	cov_a_t = cov_a_n[tid][rc_cov_f[tid]];
	seed_m* seed_tmp = NULL;

	mem_i = 0;
	r_b_v = 0;
	for(read_off = 0; read_off <= read_length - k_t; read_off += seed_l[tid])
	{
		if(read_off + k_t - 1 <= r_b_v)	continue;

		re_d = (read_off & 0X1f);

		b_t_n_r = 32 - re_d;

#ifdef HASH_KMER_READ_J
		if(re_d <= re_b)
		{
			kmer_bit = ((read_bit[read_off >> 5] & bit_tran_re[re_d]) >> ((re_b - re_d) << 1));
		}
		else
		{
			kmer_bit = (((read_bit[read_off >> 5] & bit_tran_re[re_d]) << ((re_d - re_b) << 1)) | (read_bit[(read_off >> 5) + 1] >> ((re_2bt - re_d) << 1)));
		}
#else
		tran_tmp_p = (read_bit[read_off >> 5] & bit_tran_re[re_d]);

		//or use this method to deal with: & bit_tran_re[b_t_n_r]
		kmer_bit = (((read_bit[(read_off >> 5) + 1] >> (b_t_n_r << 1)) & bit_tran_re[b_t_n_r]) | (tran_tmp_p << (re_d << 1)));

		kmer_bit >>= re_bt;
#endif

		seed_kmer = (kmer_bit & bit_tran[k_r]);

		seed_hash = (kmer_bit >> (k_r << 1));

		//find the kmer
#ifdef UNPIPATH_OFF_K20
		seed_binary_r = binsearch_offset64(seed_kmer, buffer_kmer_g, buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], buffer_hash_g[seed_hash]);
#else
		seed_binary_r = binsearch_offset(seed_kmer, buffer_kmer_g, buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], buffer_hash_g[seed_hash]);
#endif
		if(seed_binary_r == -1)
			continue;

		//binary search on unipath offset to get the unipathID
		kmer_pos_uni = buffer_off_g[seed_binary_r];//this kmer's offset on unipath
#ifdef UNPIPATH_OFF_K20
		seed_id_r = binsearch_interval_unipath64(kmer_pos_uni, buffer_seqf, result_seqf);
#else
		seed_id_r = binsearch_interval_unipath(kmer_pos_uni, buffer_seqf, result_seqf);
#endif
		ref_pos_n = buffer_pp[seed_id_r + 1] - buffer_pp[seed_id_r];

		if(ref_pos_n > pos_n_max)
		{
#ifdef	PAIR_RANDOM
			if(seed_l[tid] == pair_ran_intvp)	(*pos_ren)++;
#endif
			continue;
		}

		uni_offset_s_l = kmer_pos_uni - buffer_seqf[seed_id_r];
		uni_offset_s_r = buffer_seqf[seed_id_r + 1] - (kmer_pos_uni + k_t);

		for(left_i = 1; (left_i <= uni_offset_s_l) && (left_i <= read_off); left_i++)
		{
#ifdef UNI_SEQ64
			if(((buffer_seq[(kmer_pos_uni - left_i) >> 5] >> ((31 - ((kmer_pos_uni - left_i) & 0X1f)) << 1)) & 0X3)
			        != ((read_bit[(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3)
			  )	break;
#else
			if(((buffer_seq[(kmer_pos_uni - left_i) >> 2] >> (((kmer_pos_uni - left_i) & 0X3) << 1)) & 0X3)
			        != ((read_bit[(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3)
			  )	break;
#endif
		}

		for(right_i = 1; (right_i <= uni_offset_s_r) && (right_i <= read_length - read_off - k_t); right_i++)
		{
#ifdef UNI_SEQ64
			if(((buffer_seq[(kmer_pos_uni + k_t - 1 + right_i) >> 5] >> ((31 - ((kmer_pos_uni + k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
			        != ((read_bit[(read_off + k_t - 1 + right_i) >> 5] >> ((31 - ((read_off + k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
			  )	break;
#else
			if(((buffer_seq[(kmer_pos_uni + k_t - 1 + right_i) >> 2] >> (((kmer_pos_uni + k_t - 1 + right_i) & 0X3) << 1)) & 0X3)
			        != ((read_bit[(read_off + k_t - 1 + right_i) >> 5] >> ((31 - ((read_off + k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
			  )	break;
#endif
		}

		read_pos = read_off + 1 - left_i;
		length = k_t + left_i + right_i - 2;
		rt_pos = read_off + k_t + right_i - 2;

		seedm[tid][mem_i].uni_id = seed_id_r;
		seedm[tid][mem_i].length = length;
		seedm[tid][mem_i].read_pos = read_pos;
		seedm[tid][mem_i].ref_pos_off = uni_offset_s_l + 1 - left_i;
		seedm[tid][mem_i].ref_pos_n = ref_pos_n;

		if((left_i - 1 < read_off) || (right_i - 1 < read_length - read_off - k_t))   //cross the unipath
		{
			seedm[tid][mem_i].ui = 0;

		}
		else
		{
			seedm[tid][mem_i].ui = 1;
			seedm[tid][mem_i].ref_pos_off_r = uni_offset_s_r + 1 - right_i;
		}

		seedm[tid][mem_i].s_r_o_l = read_off - left_i;
		seedm[tid][mem_i].s_r_o_r = read_off + k_t + right_i - 1;

		memset(seedm[tid][mem_i].cov, 0Xff, (cov_a_t) << 3);

		(seedm[tid][mem_i].cov)[read_pos >> 6] = bv_64[(64 - (read_pos & 0X3f))];

		(seedm[tid][mem_i].cov)[rt_pos >> 6] ^= bv_64[(63 - (rt_pos & 0X3f))];

		for(cov_i = (read_pos >> 6); cov_i > 0; cov_i--)
			(seedm[tid][mem_i].cov)[cov_i - 1] = 0;
		for(cov_i = (rt_pos >> 6) + 1; cov_i < cov_a_t; cov_i++)
			(seedm[tid][mem_i].cov)[cov_i] = 0;

		seedm[tid][mem_i].tid = tid;

		r_b_v = read_off + k_t + right_i - 1;

		++mem_i;
	}

	//end MEM

	if(mem_i == 0)	return 0;

	//merge seeds in the same unipath
	s_uid_f[tid] = 0;
	qsort(seedm[tid], mem_i, sizeof(seed_m), compare_uniid);

	if(s_uid_f[tid] == 1)
	{
		su_i = 0;
		for(seed_i = 0; seed_i < mem_i; seed_i++)
		{
			if(uni_id_tmp != seedm[tid][seed_i].uni_id)
			{
				seedu[tid][su_i].uni_id = seedm[tid][seed_i].uni_id;
				seedu[tid][su_i].read_pos = seedm[tid][seed_i].read_pos;//choose the one with more ref positions
				seedu[tid][su_i].ref_pos_off = seedm[tid][seed_i].ref_pos_off;
				seedu[tid][su_i].ref_pos_off_r = seedm[tid][seed_i].ref_pos_off_r;
				seedu[tid][su_i].ui = seedm[tid][seed_i].ui;
				seedu[tid][su_i].s_r_o_l = seedm[tid][seed_i].s_r_o_l;
				seedu[tid][su_i].s_r_o_r = seedm[tid][seed_i].s_r_o_r;
				seedu[tid][su_i].ref_pos_n = seedm[tid][seed_i].ref_pos_n;
				memcpy(seedu[tid][su_i].cov, seedm[tid][seed_i].cov, (cov_a_t) << 3);

				seedu[tid][su_i].tid = tid;

				//cal length
				if(su_i > 0)
				{
					seedu[tid][su_i - 1].length = 0;
					for(cov_i = 0; cov_i < cov_a_t; cov_i++)
						seedu[tid][su_i - 1].length += popcount_3((seedu[tid][su_i - 1].cov)[cov_i]);
				}
				++su_i;
			}
			else
			{
				for(cov_i = 0; cov_i < cov_a_t; cov_i++)
					(seedu[tid][su_i - 1].cov)[cov_i] |= (seedm[tid][seed_i].cov)[cov_i];
			}

			uni_id_tmp = seedm[tid][seed_i].uni_id;
		}

		seedu[tid][su_i - 1].length = 0;
		for(cov_i = 0; cov_i < cov_a_t; cov_i++)
			seedu[tid][su_i - 1].length += popcount_3((seedu[tid][su_i - 1].cov)[cov_i]);

		//end merge of same unipath
		seed_tmp = seedu[tid];
	}
	else
	{
		su_i = mem_i;
		seed_tmp = seedm[tid];
	}

	//uniqueness
	qsort(seed_tmp, su_i, sizeof(seed_m), compare_posn);

	uni_id = seed_tmp[0].uni_id;
	ref_off1 = seed_tmp[0].ref_pos_off;
	r_p1 = seed_tmp[0].read_pos;

	pos_tmp_p = buffer_p + buffer_pp[uni_id];//119843814
	pos_n = buffer_pp[uni_id + 1] - buffer_pp[uni_id];

	//fingerprint
	for(id_i = 0; id_i < pos_n; id_i++)
	{
		seedsets[tid][id_i].seed_set = pos_tmp_p[id_i] + ref_off1 - r_p1;
		memcpy(seedsets[tid][id_i].cov, seed_tmp[0].cov, (cov_a_t) << 3);

		seedsets[tid][id_i].ui = seed_tmp[0].ui;
		seedsets[tid][id_i].ref_pos_off = ref_off1;
		seedsets[tid][id_i].ref_pos_off_r = seed_tmp[0].ref_pos_off_r;

		seedsets[tid][id_i].s_r_o_l = seed_tmp[0].s_r_o_l;
		seedsets[tid][id_i].s_r_o_r = seed_tmp[0].s_r_o_r;

		seedsets[tid][id_i].tid = tid;
	}

	set_n += pos_n;
	seed_set_off[tid][0] = 0;
	seed_set_off[tid][1] = pos_n;
	off_n += 2;

#ifdef MERGE_B_R_VS
	for(seed_i = 1; seed_i < su_i; seed_i++)
	{
		uni_id = seed_tmp[seed_i].uni_id;

		ref_off1 = seed_tmp[seed_i].ref_pos_off;
		r_p1 = seed_tmp[seed_i].read_pos;

		reduce_pn = buffer_pp[uni_id + 1];

		bn = reduce_pn - buffer_pp[uni_id];

		for(seed_set_i = 0; seed_set_i < seed_i; seed_set_i++)
		{
			reduce_sn = seed_set_off[tid][seed_set_i + 1];
			sn = reduce_sn - seed_set_off[tid][seed_set_i];

			if(bn <= sn)
			{
				reduce_s = seed_set_off[tid][seed_set_i];

				for(b_p_i = buffer_pp[uni_id]; b_p_i < reduce_pn; b_p_i++)
				{
					posp = buffer_p[b_p_i] + ref_off1 - r_p1;
#ifdef UNPIPATH_OFF_K20
					b_r_s = binsearch_seed_set_reduce64(posp, seedsets[tid], reduce_sn - reduce_s, reduce_s, tid);
#else
					b_r_s = binsearch_seed_set_reduce(posp, seedsets[tid], reduce_sn - reduce_s, reduce_s, tid);
#endif
					if(b_r_s != -1)
					{
						for(cov_i = 0; cov_i < cov_a_t; cov_i++)
							(seedsets[tid][b_r_s].cov)[cov_i] |= (seed_tmp[seed_i].cov)[cov_i];

						pos_add[tid][b_p_i - buffer_pp[uni_id]] = 1;

						if(b_r_s == reduce_sn - 1)
							break;

						reduce_s = b_r_s + 1;
					}
					else
					{
						if(g_low[tid] == reduce_sn - reduce_s)
							break;
						else	reduce_s += g_low[tid];
					}
				}
			}
			else
			{
				reduce_p = buffer_pp[uni_id];

				for(b_p_i = seed_set_off[tid][seed_set_i]; b_p_i < reduce_sn; b_p_i++)
				{
					posp = seedsets[tid][b_p_i].seed_set + r_p1 - ref_off1;
#ifdef UNPIPATH_OFF_K20
					b_r_s = binsearch_seed_pos_reduce64(posp, buffer_p, reduce_pn - reduce_p, reduce_p, tid);
#else
					b_r_s = binsearch_seed_pos_reduce(posp, buffer_p, reduce_pn - reduce_p, reduce_p, tid);
#endif
					if(b_r_s != -1)
					{
						for(cov_i = 0; cov_i < cov_a_t; cov_i++)
							(seedsets[tid][b_p_i].cov)[cov_i] |= (seed_tmp[seed_i].cov)[cov_i];

						pos_add[tid][b_r_s - buffer_pp[uni_id]] = 1;

						if(b_r_s == reduce_pn - 1)
							break;

						reduce_p = b_r_s + 1;
					}
					else
					{
						if(g_low[tid] == reduce_pn - reduce_p)
							break;
						else	reduce_p += g_low[tid];
					}
				}
			}
		}

		for(b_p_i = 0; b_p_i < bn; b_p_i++)
		{
			if(pos_add[tid][b_p_i] == 0)
			{
				seedsets[tid][set_n].seed_set = buffer_p[buffer_pp[uni_id] + b_p_i] + ref_off1 - r_p1;

				memcpy(seedsets[tid][set_n].cov, seed_tmp[seed_i].cov, (cov_a_t) << 3);

				seedsets[tid][set_n].ui = seed_tmp[seed_i].ui;
				seedsets[tid][set_n].ref_pos_off = ref_off1;
				seedsets[tid][set_n].ref_pos_off_r = seed_tmp[seed_i].ref_pos_off_r;

				seedsets[tid][set_n].s_r_o_l = seed_tmp[seed_i].s_r_o_l;
				seedsets[tid][set_n].s_r_o_r = seed_tmp[seed_i].s_r_o_r;

				seedsets[tid][set_n].tid = tid;

				++set_n;
			}
		}

		memset(pos_add[tid], 0, bn);

		seed_set_off[tid][seed_i + 1] = set_n;
		++off_n;
	}
#endif

	for(off_i = 0; off_i < off_n - 1; off_i++)
	{
		interval_n = seed_set_off[tid][off_i + 1] - seed_set_off[tid][off_i];

		if(interval_n == 0)	continue;

		qsort(seedsets[tid] + seed_set_off[tid][off_i], interval_n, sizeof(seed_sets), compare_sets);

		//copy elements from seed_tmp[off_i] as for other elements such as quality
		set_i_pre = seed_set_off[tid][off_i];
		seedpa[spa_i].pos_start = set_i_pre;
		memcpy(seedpa[spa_i].cov, seedsets[tid][set_i_pre].cov, (cov_a_t) << 3);

		seedpa[spa_i].ref_pos_off = seedsets[tid][set_i_pre].ref_pos_off;
		seedpa[spa_i].ref_pos_off_r = seedsets[tid][set_i_pre].ref_pos_off_r;
		seedpa[spa_i].ui = seedsets[tid][set_i_pre].ui;

		seedpa[spa_i].s_r_o_l = seedsets[tid][set_i_pre].s_r_o_l;
		seedpa[spa_i].s_r_o_r = seedsets[tid][set_i_pre].s_r_o_r;

		//cal number of 1 (length)
		seedpa[spa_i].length = 0;
		for(pop_i = 0; pop_i < cov_a_t; pop_i++)
			seedpa[spa_i].length += popcount_3((seedpa[spa_i].cov)[pop_i]);
		if(seedpa[spa_i].length > max_sets_n)	max_sets_n = seedpa[spa_i].length;

		++spa_i;

		seed_same_n = 1;

		//copy pos to seed_set_pos
		(*seed_set_pos)[set_pos_n++] = seedsets[tid][set_i_pre].seed_set; //+ ref_off - r_p

		for(set_i = set_i_pre + 1; set_i < seed_set_off[tid][off_i + 1]; set_i++)
		{
			set_c_r = memcmp(seedsets[tid][set_i_pre].cov, seedsets[tid][set_i].cov, (cov_a_t) << 3);
			if(set_c_r != 0)
			{
				seedpa[spa_i - 1].pos_n = seed_same_n;

				//copy elements from seed_tmp[off_i] as for other elements such as quality
				seedpa[spa_i].pos_start = set_i;
				memcpy(seedpa[spa_i].cov, seedsets[tid][set_i].cov, (cov_a_t) << 3);

				seedpa[spa_i].ref_pos_off = seedsets[tid][set_i].ref_pos_off;
				seedpa[spa_i].ref_pos_off_r = seedsets[tid][set_i].ref_pos_off_r;
				seedpa[spa_i].ui = seedsets[tid][set_i].ui;

				seedpa[spa_i].s_r_o_l = seedsets[tid][set_i].s_r_o_l;
				seedpa[spa_i].s_r_o_r = seedsets[tid][set_i].s_r_o_r;

				//cal number of 1 (length)
				seedpa[spa_i].length = 0;
				for(pop_i = 0; pop_i < cov_a_t; pop_i++)
					seedpa[spa_i].length += popcount_3((seedpa[spa_i].cov)[pop_i]);
				if(seedpa[spa_i].length > max_sets_n)	max_sets_n = seedpa[spa_i].length;

				++spa_i;

				seed_same_n = 1;
			}
			else	++seed_same_n;

			//copy pos to seed_set_pos
			(*seed_set_pos)[set_pos_n++] = seedsets[tid][set_i].seed_set;// + ref_off - r_p

			set_i_pre = set_i;
		}

		seedpa[spa_i - 1].pos_n = seed_same_n;
	}

#ifdef	SINGLE_COV_FILT
	qsort(seedpa, spa_i, sizeof(seed_pa), compare_plen);

	for(off_i = 0; off_i < spa_i; off_i++)
		if(seedpa[off_i].length < max_sets_n - length_reduce)
		{
			++off_i;
			break;
		}

	(*spa1_i) = off_i;
#else
	(*spa1_i) = spa_i;
#endif

	(*seed_num) = spa_i;
	(*max_seed_length) = max_sets_n;

#ifdef	PR_COV_FILTER

	if(seed_l[tid] == COV_FIT_SI)
	{
		if(max_sets_n > (read_length >> 1))	cov_filt_f[tid] = 1;
	}

#endif

	return seedpa;
}

#ifdef UNPIPATH_OFF_K20
void merge_pair_end64(uint32_t spa1_i, uint32_t spa2_i, seed_pa* seed_pr1, seed_pa* seed_pr2, uint64_t* seed_set_pos1, uint64_t* seed_set_pos2, uint8_t rc_i, uint8_t tid)
#else
void merge_pair_end(uint32_t spa1_i, uint32_t spa2_i, seed_pa* seed_pr1, seed_pa* seed_pr2, uint32_t* seed_set_pos1, uint32_t* seed_set_pos2, uint8_t rc_i, uint8_t tid)
#endif
{

#ifdef UNPIPATH_OFF_K20
	uint64_t posi = 0;
	uint64_t posp = 0;
#else
	uint32_t posi = 0;
	uint32_t posp = 0;
#endif

	//uint32_t pospi = 0;
	uint32_t seed1_i = 0;
	uint32_t seed2_i = 0;
	uint32_t psp_i = 0;
	uint32_t mat_posi1 = 0;

	int64_t b_s_p_r_i = 0;
	int64_t b_s_p_r = 0;
	int64_t reduce_r = 0;
	int64_t reduce_rn = 0;

	uint16_t end_dis = 0;

	if(rc_i == 0)	end_dis = end_dis1[tid];
	else	end_dis = end_dis2[tid];

	mat_posi1 = mat_posi[tid];

	//match each side of pair-end read
	for(seed1_i = 0; seed1_i < spa1_i; seed1_i++)
	{
		for(seed2_i = 0; seed2_i < spa2_i; seed2_i++)
		{
#ifdef PAIR_B_R_V
			if(seed_pr1[seed1_i].pos_n <= seed_pr2[seed2_i].pos_n)
			{
				reduce_r = seed_pr2[seed2_i].pos_start;
				reduce_rn = seed_pr2[seed2_i].pos_n + seed_pr2[seed2_i].pos_start;

				for(psp_i = 0; psp_i < seed_pr1[seed1_i].pos_n; psp_i++)
				{
					posi = seed_set_pos1[seed_pr1[seed1_i].pos_start + psp_i];

					posp = posi + end_dis;

					//can add if( != 0Xffffffff) in following binary search function
#ifdef UNPIPATH_OFF_K20
					b_s_p_r = binsearch_pair_pos_reduce64(posp, seed_set_pos2, reduce_rn - reduce_r, reduce_r, tid);
#else
					b_s_p_r = binsearch_pair_pos_reduce(posp, seed_set_pos2, reduce_rn - reduce_r, reduce_r, tid);
#endif
					if(b_s_p_r != -1)
					{
						de_m_p[tid] = 0;

						if(mat_posi1 < CUS_SEED_SET)
						{
							mat_pos1[tid][mat_posi1] = posi;
							mat_pos2[tid][mat_posi1] = seed_set_pos2[b_s_p_r];
							seed_no1[tid][mat_posi1] = seed1_i;
							seed_no2[tid][mat_posi1] = seed2_i;
							mat_rc[tid][mat_posi1] = rc_i;
							++mat_posi1;
						}
#ifdef	PAIR_SEED_LENGTH_FILT
						if(seed_pr1[seed1_i].length + seed_pr2[seed2_i].length > max_lengtht[tid])
							max_lengtht[tid] = seed_pr1[seed1_i].length + seed_pr2[seed2_i].length;
#endif
						if(b_s_p_r == reduce_rn - 1)
							break;

#ifdef	PAIR_TRAVERSE
						for(b_s_p_r_i = b_s_p_r + 1; (seed_set_pos2[b_s_p_r_i] < posp + devi) && (b_s_p_r_i < reduce_rn); b_s_p_r_i++)
						{
							if(mat_posi1 < CUS_SEED_SET)
							{
								mat_pos1[tid][mat_posi1] = posi;
								mat_pos2[tid][mat_posi1] = seed_set_pos2[b_s_p_r_i];
								seed_no1[tid][mat_posi1] = seed1_i;
								seed_no2[tid][mat_posi1] = seed2_i;
								mat_rc[tid][mat_posi1] = rc_i;
								++mat_posi1;
							}
						}

						for(b_s_p_r_i = b_s_p_r - 1; (seed_set_pos2[b_s_p_r_i] > posp - devi) && (b_s_p_r_i >= reduce_r); b_s_p_r_i--)
						{
							if(mat_posi1 < CUS_SEED_SET)
							{
								mat_pos1[tid][mat_posi1] = posi;
								mat_pos2[tid][mat_posi1] = seed_set_pos2[b_s_p_r_i];
								seed_no1[tid][mat_posi1] = seed1_i;
								seed_no2[tid][mat_posi1] = seed2_i;
								mat_rc[tid][mat_posi1] = rc_i;
								++mat_posi1;
							}
						}
#endif
						reduce_r = b_s_p_r + 1;
					}
					else
					{
						if(r_low[tid] == reduce_rn - reduce_r)
							break;
						else	reduce_r += r_low[tid];
					}
				}
			}
			else
			{
				reduce_r = seed_pr1[seed1_i].pos_start;
				reduce_rn = seed_pr1[seed1_i].pos_n + seed_pr1[seed1_i].pos_start;

				for(psp_i = 0; psp_i < seed_pr2[seed2_i].pos_n; psp_i++)
				{
					posi = seed_set_pos2[seed_pr2[seed2_i].pos_start + psp_i];
					posp = posi - end_dis;

					//can add if( != 0Xffffffff) in following binary search function
#ifdef UNPIPATH_OFF_K20
					b_s_p_r = binsearch_pair_pos_reduce64(posp, seed_set_pos1, reduce_rn - reduce_r, reduce_r, tid);
#else
					b_s_p_r = binsearch_pair_pos_reduce(posp, seed_set_pos1, reduce_rn - reduce_r, reduce_r, tid);
#endif
					if(b_s_p_r != -1)
					{
						de_m_p[tid] = 0;
						if(mat_posi1 < CUS_SEED_SET)
						{
							mat_pos1[tid][mat_posi1] = seed_set_pos1[b_s_p_r];
							mat_pos2[tid][mat_posi1] = posi;
							seed_no1[tid][mat_posi1] = seed1_i;
							seed_no2[tid][mat_posi1] = seed2_i;
							mat_rc[tid][mat_posi1] = rc_i;
							++mat_posi1;
						}

#ifdef	PAIR_SEED_LENGTH_FILT
						if(seed_pr1[seed1_i].length + seed_pr2[seed2_i].length > max_lengtht[tid])
							max_lengtht[tid] = seed_pr1[seed1_i].length + seed_pr2[seed2_i].length;
#endif
						if(b_s_p_r == reduce_rn - 1)
							break;

#ifdef	PAIR_TRAVERSE
						for(b_s_p_r_i = b_s_p_r + 1; (seed_set_pos1[b_s_p_r_i] < posp + devi) && (b_s_p_r_i < reduce_rn); b_s_p_r_i++)
						{
							if(mat_posi1 < CUS_SEED_SET)
							{
#ifdef	MERGE_MODIFY
								mat_pos1[tid][mat_posi1] = seed_set_pos1[b_s_p_r_i];
								mat_pos2[tid][mat_posi1] = posi;
#else
								mat_pos1[tid][mat_posi1] = posi;
								mat_pos2[tid][mat_posi1] = seed_set_pos1[b_s_p_r_i];
#endif
								seed_no1[tid][mat_posi1] = seed1_i;
								seed_no2[tid][mat_posi1] = seed2_i;
								mat_rc[tid][mat_posi1] = rc_i;
								++mat_posi1;
							}
						}

						for(b_s_p_r_i = b_s_p_r - 1; (seed_set_pos1[b_s_p_r_i] > posp - devi) && (b_s_p_r_i >= reduce_r); b_s_p_r_i--)
						{
							if(mat_posi1 < CUS_SEED_SET)
							{
#ifdef	MERGE_MODIFY
								mat_pos1[tid][mat_posi1] = seed_set_pos1[b_s_p_r_i];
								mat_pos2[tid][mat_posi1] = posi;
#else
								mat_pos1[tid][mat_posi1] = posi;
								mat_pos2[tid][mat_posi1] = seed_set_pos1[b_s_p_r_i];
#endif
								seed_no1[tid][mat_posi1] = seed1_i;
								seed_no2[tid][mat_posi1] = seed2_i;
								mat_rc[tid][mat_posi1] = rc_i;
								++mat_posi1;
							}
						}
#endif
						reduce_r = b_s_p_r + 1;
					}
					else
					{
						if(r_low[tid] == reduce_rn - reduce_r)
							break;
						else	reduce_r += r_low[tid];
					}
				}
			}
#endif
		}
	}
	mat_posi[tid] = mat_posi1;
}

void pair_sam_output(uint8_t tid, uint16_t read_length1, uint16_t read_length2, uint16_t f_cigarn, cnt_re* cnt, uint32_t seqi, uint16_t lv_k1, uint16_t lv_k2, int16_t pound_pos1_f_forward, int16_t pound_pos1_f_reverse, int16_t pound_pos1_r_forward, int16_t pound_pos1_r_reverse, int16_t pound_pos2_f_forward, int16_t pound_pos2_f_reverse, int16_t pound_pos2_r_forward, int16_t pound_pos2_r_reverse, uint8_t cir_n)
{
	//alignment in optimal and suboptimal vectors and output to sam file
	uint8_t sam_flag1 = 0;
	uint8_t sam_flag2 = 0;

	uint16_t v_cnt_out = 0;
	uint16_t vs_cnt_out = 0;
	uint16_t d_n1 = 0;
	uint16_t d_n2 = 0;
	uint16_t i_n1 = 0;
	uint16_t i_n2 = 0;
	uint16_t f_c = 0;
	uint16_t pchl = 0;
	uint16_t f_i = 0;
	uint16_t s_o = 0;
	uint16_t sn = 0;
	uint16_t snt = 0;
	uint16_t read_b_i = 0;
	uint16_t sam_seq_i = 0;
	uint16_t xa_i = 0;
	uint16_t xa_i1 = 0;
	uint16_t xa_i2 = 0;
	uint16_t v_cnt_i = 0;
	uint16_t va_cnt_i = 0;
	uint32_t v_cnt = 0;
	uint32_t vs_cnt = 0;

	int low = 0;
	int high = 0;
	int mid = 0;
	int chr_re = 0;
	int bit_char_i = 0;
	int lv_re1f = 0;
	int lv_re1b = 0;
	int lv_re2f = 0;
	int lv_re2b = 0;
	int lv_re1 = 0;
	int lv_re2 = 0;
	int sam_qual1 = 0;
	int sam_qual2 = 0;
	int m_m_n = 0;
	int16_t m_n_f = 0;
	int16_t m_n_b = 0;

	int16_t s_r_o_l = 0;
	int16_t s_r_o_r = 0;
	int16_t lv_up_left = 0;
	int16_t lv_up_right = 0;
	int16_t lv_down_right = 0;
	int16_t lv_down_left = 0;
	int16_t pound_pos_1_f = 0;
	int16_t	pound_pos_1_r = 0;
	int16_t	pound_pos_2_f = 0;
	int16_t	pound_pos_2_r = 0;

	int64_t sam_pos1 = 0;
	int64_t sam_pos2  = 0;
	int64_t sam_cross = 0;
	int64_t sam_pos1_pr = 0;
	int64_t sam_pos2_pr = 0;
	uint8_t rc_i = 0;
	uint8_t rc_ii = 0;
	char sam_seq1[MAX_READLEN] = {};
	char sam_seq2[MAX_READLEN] = {};
	char cigarBuf1[MAX_LV_CIGAR] = {};
	char cigarBuf2[MAX_LV_CIGAR] = {};
	uint16_t f_cigar[MAX_LV_CIGAR];
	char str_o[MAX_LV_CIGAR];
	char b_cigar[MAX_LV_CIGAR];
	char* pch = NULL;
	char* saveptr = NULL;

#ifdef	KSW_ALN_PAIR
	int band_with = 33;//100
	int x_score, y_score, op_score, len_score, op_score1, len_score1, op_score2, len_score2;
	int k_start1 = 0;
	int k_start2 = 0;
	int k_middle = 0;
	uint32_t* cigar1 = NULL;
	uint32_t* cigar2 = NULL;
	int n_cigar1 = 0;
	int n_cigar2 = 0;
	int nm_score = 0;
#endif

#ifdef	QUAL_FILT_LV_OUT
	uint8_t* qual_filt_lv_1 = NULL;
	uint8_t* qual_filt_lv_1_o = NULL;
	uint8_t* qual_filt_lv_2 = NULL;
	uint8_t* qual_filt_lv_2_o = NULL;
#endif

	char cigar_p1[MAX_LV_CIGARCOM];
	char cigar_p2[MAX_LV_CIGARCOM];

	uint16_t s_offset1 = 0;
	uint16_t s_offset2 = 0;

#ifdef	MAPPING_QUALITY
	float log_tmp = 0;
	float* mp_subs_1 = NULL;
	float* mp_subs_1_o = NULL;
	float* mp_subs_2 = NULL;
	float* mp_subs_2_o = NULL;
	uint8_t mp_flag = 0;
#endif

	uint8_t ops_flag = 1;

#ifdef UNPIPATH_OFF_K20
	uint64_t x = 0;
#else
	uint32_t x = 0;
#endif

	uint64_t chr_re_tmp = 0;

	v_cnt = cnt->v_cnt;
	vs_cnt = cnt->vs_cnt;

	v_cnt_out = (cus_ali_n < v_cnt) ? cus_ali_n:v_cnt;

	v_cnt_i = 0;

#ifdef	CIGAR_LEN_ERR
	int cigar_len = 0;
	int cigar_len_tmp = 0;
	int cigar_len_re = 0;
	uint16_t pchl_tmp = 0;
	uint16_t s_o_tmp = 0;
	char* pch_tmp = NULL;
	char* saveptr_tmp = NULL;
	char cigar_tmp[MAX_LV_CIGARCOM];
#endif

#ifdef	ALTER_DEBUG
	if((v_cnt > 1) && (rep_go[tid]))	v_cnt_i = seed_length_arr[tid][v_cnt_i].index;
#endif

#ifdef	MAPPING_QUALITY
	if(v_cnt > 1)
	{
		sam_qual1 = 0;
		sam_qual2 = 0;
		mp_flag = 0;
	}
	else if(vs_cnt == 0)
	{
		sam_qual1 = 60;
		sam_qual2 = 60;
		mp_flag = 0;
	}
	else
	{
		mp_flag = 1;
	}
#endif

	x = op_vector_pos1[tid][v_cnt_i];
	low = 0;
	high = chr_file_n - 1;

	while ( low <= high )
	{
		mid = (low + high) >> 1;
		if(x < (chr_end_n[mid]))
		{
			high = mid - 1;
		}
		else if(x > (chr_end_n[mid]))
		{
			low = mid + 1;
		}
		else
		{
			chr_re =  mid;
			break;
		}
		chr_re = low;
	}

	sam_pos1 = op_vector_pos1[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;
	sam_pos2 = op_vector_pos2[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;

	if(op_rc[tid][v_cnt_i] == 0)
	{
		sam_flag1 = 99;
		sam_flag2 = 147;

#ifdef	CHAR_CP
		read_bit_1[tid] = read_bit1[tid][0];
		read_bit_2[tid] = read_bit2[tid][1];
#else
		for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
			sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

		sam_seq2[sam_seq_i] = '\0';
		strrev1(sam_seq2);

		strcpy(sam_seq1, seqio[seqi].read_seq1);
#endif


#ifdef	QUAL_FILT_LV_OUT
		qual_filt_lv_1 = qual_filt_lv1[tid][0];
		qual_filt_lv_1_o = qual_filt_lv1[tid][1];
		qual_filt_lv_2 = qual_filt_lv2[tid][1];
		qual_filt_lv_2_o = qual_filt_lv2[tid][0];
#endif
		sam_cross = sam_pos2 + read_length2 - sam_pos1;

		pound_pos_1_f = pound_pos1_f_forward;
		pound_pos_1_r = pound_pos1_r_forward;
		pound_pos_2_f = pound_pos2_f_reverse;
		pound_pos_2_r = pound_pos2_r_reverse;
#ifdef	MAPPING_QUALITY
		if(mp_flag)
		{
			mp_subs_1 = mp_subs1[tid][0];
			mp_subs_1_o = mp_subs1[tid][1];
			mp_subs_2 = mp_subs2[tid][1];
			mp_subs_2_o = mp_subs2[tid][0];
		}
#endif
	}
	else
	{
		sam_flag1 = 83;
		sam_flag2 = 163;

#ifdef	CHAR_CP
		read_bit_1[tid] = read_bit1[tid][1];
		read_bit_2[tid] = read_bit2[tid][0];
#else
		for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
			sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

		sam_seq1[sam_seq_i] = '\0';

		strrev1(sam_seq1);
		strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif


#ifdef	QUAL_FILT_LV_OUT
		qual_filt_lv_1 = qual_filt_lv1[tid][1];
		qual_filt_lv_1_o = qual_filt_lv1[tid][0];
		qual_filt_lv_2 = qual_filt_lv2[tid][0];
		qual_filt_lv_2_o = qual_filt_lv2[tid][1];
#endif

		sam_cross = sam_pos2 - read_length1 - sam_pos1;

		pound_pos_1_f = pound_pos1_f_reverse;
		pound_pos_1_r = pound_pos1_r_reverse;
		pound_pos_2_f = pound_pos2_f_forward;
		pound_pos_2_r = pound_pos2_r_forward;
#ifdef	MAPPING_QUALITY
		if(mp_flag)
		{
			mp_subs_1 = mp_subs1[tid][1];
			mp_subs_1_o = mp_subs1[tid][0];
			mp_subs_2 = mp_subs2[tid][0];
			mp_subs_2_o = mp_subs2[tid][1];
		}
#endif
	}

	d_n1 = 0;
	i_n1 = 0;
	s_offset1 = 0;
	s_offset2 = 0;
	s_r_o_l = op_dm_l1[tid][v_cnt_i];
	s_r_o_r = op_dm_r1[tid][v_cnt_i];

#ifdef	MAPPING_QUALITY
	if(mp_flag)	sub_t[tid] = 0;
#endif

	if((s_r_o_l == 0) && (s_r_o_r == 0))
	{
		strcpy(cigar_p1, cigar_m1[tid]);

#ifdef	MAPPING_QUALITY

		if(mp_flag)
		{
			for (bit_char_i = 0, read_b_i = 32; bit_char_i < read_length1; bit_char_i++, read_b_i++)
				if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((op_vector_seq1[tid][v_cnt_i][read_b_i >> 5] >> ((31 - (read_b_i & 0X1f)) << 1)) & 0X3))
				{
					sub_t[tid] += mp_subs_1[bit_char_i];
				}
		}
#endif

	}
	else     //indel
	{
		if((cir_n == cir_fix_n) && (local_ksw))
		{
#ifdef	KSW_ALN_PAIR

#ifdef	CHAR_CP
			for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
				read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
			for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
				read_char[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif
			for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i > -1
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			ksw_global(op_dm_kl1[tid][v_cnt_i], read_char[tid], op_dm_kl1[tid][v_cnt_i], ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar1, &cigar1);

#ifdef	CHAR_CP
			for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
				read_char2[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
			for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
				read_char2[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif

			for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length + 64
				ali_ref_seq2[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			ksw_global(op_dm_kr1[tid][v_cnt_i], read_char2[tid], op_dm_kr1[tid][v_cnt_i], ali_ref_seq2[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar2, &cigar2);

			m_m_n = s_r_o_r - s_r_o_l - 1;

			nm_score = 0;
			snt = 0;
			if (n_cigar1)
			{
				op_score1  = cigar1[0]&0xf;
				len_score1 = cigar1[0]>>4;
			}
			else	op_score1 = 3;

			if (n_cigar2)
			{
				op_score2  = cigar2[0]&0xf;
				len_score2 = cigar2[0]>>4;
			}
			else	op_score2 = 3;

			if (s_r_o_l >= op_dm_kl1[tid][v_cnt_i])
			{
				sn = sprintf(cigar_p1 + snt, "%dS", s_r_o_l + 1 - op_dm_kl1[tid][v_cnt_i]);
				snt += sn;
			}

			if ((s_r_o_l != -1) && (s_r_o_r != read_length1))
			{
				if ((op_score1 == 0) && (op_score2 == 0))
				{
					k_start1 = 0;
					k_start2 = 1;
					k_middle = len_score1 + len_score2 + m_m_n;
				}
				else if (op_score1 == 0)
				{
					k_start1 = 0;
					k_start2 = 0;
					k_middle = len_score1 + m_m_n;
				}
				else if (op_score2 == 0)
				{
					k_start1 = -1;
					k_start2 = 1;
					k_middle = len_score2 + m_m_n;
				}
				else
				{
					k_start1 = -1;
					k_start2 = 0;
					k_middle = m_m_n;
				}

				x_score = y_score = 0;
				for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
				{
					op_score  = cigar1[bit_char_i]&0xf;
					len_score = cigar1[bit_char_i]>>4;

					sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
					snt += sn;

					if (op_score == 0)
					{
						// match
						for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
							if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
							{
#ifdef	MAPPING_QUALITY
								if(mp_flag)	sub_t[tid] += mp_subs_1_o[read_length1 + n_cigar1 - s_r_o_l - x_score - read_b_i - 2];
#endif
								++nm_score;
							}
						x_score += len_score;
						y_score += len_score;
					}
					else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
					else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
				}

				sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
				snt += sn;

				x_score = y_score = 0;
				for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
				{
					op_score  = cigar2[bit_char_i]&0xf;
					len_score = cigar2[bit_char_i]>>4;

					sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
					snt += sn;

					if (op_score == 0)
					{
						// match
						for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
							if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
							{
#ifdef	MAPPING_QUALITY
								if(mp_flag)	sub_t[tid] += mp_subs_1[s_r_o_r + x_score + read_b_i];
#endif
								++nm_score;
							}
						x_score += len_score;
						y_score += len_score;
					}
					else if (op_score == 1) x_score += len_score, nm_score += len_score;
					else if (op_score == 2) y_score += len_score, nm_score += len_score;
				}
			}
			else if (s_r_o_l == -1)
			{
				if (op_score2 == 0)
				{
					k_start2 = 1;
					k_middle = len_score2 + m_m_n;
				}
				else
				{
					k_start2 = 0;
					k_middle = m_m_n;
				}
				sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
				snt += sn;

				x_score = y_score = 0;
				for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
				{
					op_score  = cigar2[bit_char_i]&0xf;
					len_score = cigar2[bit_char_i]>>4;

					sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
					snt += sn;

					if (op_score == 0)
					{
						// match
						for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
							if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
							{
#ifdef	MAPPING_QUALITY
								if(mp_flag)	sub_t[tid] += mp_subs_1[s_r_o_r + x_score + read_b_i];
#endif
								++nm_score;
							}
						x_score += len_score;
						y_score += len_score;
					}
					else if (op_score == 1) x_score += len_score, nm_score += len_score;
					else if (op_score == 2) y_score += len_score, nm_score += len_score;
				}
			}
			else
			{
				if (op_score1 == 0)
				{
					k_start1 = 0;
					k_middle = len_score1 + m_m_n;
				}
				else
				{
					k_start1 = -1;
					k_middle = m_m_n;
				}
				x_score = y_score = 0;
				for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
				{
					op_score  = cigar1[bit_char_i]&0xf;
					len_score = cigar1[bit_char_i]>>4;

					sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
					snt += sn;

					if (op_score == 0)
					{
						// match
						for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
							if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
							{
#ifdef	MAPPING_QUALITY
								if(mp_flag)	sub_t[tid] += mp_subs_1_o[read_length1 + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
								++nm_score;
							}
						x_score += len_score;
						y_score += len_score;
					}
					else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
					else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
				}

				sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
				snt += sn;
			}

			if (read_length1 - s_r_o_r > op_dm_kr1[tid][v_cnt_i])
			{
				sn = sprintf(cigar_p1 + snt, "%dS", read_length1 - s_r_o_r - op_dm_kr1[tid][v_cnt_i]);
				snt += sn;
			}
			//sn = sprintf(cigar_p1 + snt, "\0");
			//snt += sn;

			if (n_cigar1)	free(cigar1);
			if (n_cigar2)	free(cigar2);

			op_dm_ex1[tid][v_cnt_i] = nm_score;
#endif
		}
		else
		{
			if(pound_pos_1_f >= s_r_o_r)   //1
			{
				lv_up_left = 0;
				lv_up_right = s_r_o_l;
				lv_down_right = pound_pos_1_f;
				lv_down_left = s_r_o_r;
				m_n_f = 0;
				m_n_b = read_length1 - pound_pos_1_f;
				m_m_n = s_r_o_r - s_r_o_l - 1;
			}
			else if(pound_pos_1_r <= s_r_o_l + 1)     //5
			{
				lv_up_left = pound_pos_1_r;//
				lv_up_right = s_r_o_l;
				lv_down_right = read_length1;
				lv_down_left = s_r_o_r;
				m_n_f = pound_pos_1_r;
				m_n_b = 0;
				m_m_n = s_r_o_r - s_r_o_l - 1;
			}
			else if((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
			{
				lv_up_left = 0;
				lv_up_right = pound_pos_1_f - 1;
				lv_down_right = read_length1;
				lv_down_left = pound_pos_1_r;
				m_n_f = 0;
				m_n_b = 0;
				m_m_n = pound_pos_1_r - pound_pos_1_f;
			}
			else if((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
			{
				lv_up_left = 0;
				lv_up_right = s_r_o_l;
				lv_down_right = read_length1;
				lv_down_left = pound_pos_1_r;
				m_n_f = 0;
				m_n_b = 0;
				m_m_n = read_length1 - s_r_o_l - 1;
			}
			else     //4
			{
				lv_up_left = 0;
				lv_up_right = -1;
				lv_down_right = read_length1;
				lv_down_left = s_r_o_r;
				m_n_f = 0;
				m_n_b = 0;
				m_m_n = s_r_o_r;
			}


#ifdef	QUAL_FILT_LV_OUT


#ifdef	CHAR_CP
			for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
				read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
			for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
				read_char[tid][read_b_i] = sam_seq1[bit_char_i];

			for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
				ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY
			if(mp_flag)
				computeEditDistanceWithCigar_s_mis_left_mp(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length1 - 1 - lv_up_right, &s_offset1, mp_subs_1_o + read_length1 - 1 - lv_up_right, &(sub_t[tid]));//, 0, op_dm_sl1[tid][v_cnt_i]
			else
				computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length1 - 1 - lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]
#else
			computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length1 - 1 - lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]
#endif


#ifdef CHAR_CP
			for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
				read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
			for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
				read_char[tid][read_b_i] = sam_seq1[bit_char_i];

			for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
				ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif


#ifdef	MAPPING_QUALITY
			if(mp_flag)
				computeEditDistanceWithCigar_s_mis_mp(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + lv_down_left, mp_subs_1 + lv_down_left, &(sub_t[tid]));//, 0, op_dm_sr1[tid][v_cnt_i]
			else	computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]
#else
			computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#else

#ifdef	CHAR_CP
			for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
				read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
			for(bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
				read_char[tid][read_b_i] = sam_seq1[bit_char_i];

			for(bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
				ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY
			if(mp_flag)
				computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid], mp_subs_1_o + read_length1 - 1 - lv_up_right, &(sub_t[tid]));//, 0, op_dm_sl1[tid][v_cnt_i]
			else
				computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#else
			computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#endif

#ifdef	CHAR_CP
			for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
				read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
			for(bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
				read_char[tid][read_b_i] = sam_seq1[bit_char_i];

			for(bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
				ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif


#ifdef	MAPPING_QUALITY
			if(mp_flag)
				computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid], mp_subs_1 + lv_down_left, &(sub_t[tid]));//, 0, op_dm_sr1[tid][v_cnt_i]
			else	computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#else
			computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#endif

			//deal with front and back lv cigar
			strncpy(str_o, cigarBuf1, f_cigarn);
			s_o = 0;
			f_c = 0;

			pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

			while (pch != NULL)
			{
				pchl = strlen(pch);
				f_cigar[f_cigarn - f_c - 2] = atoi(pch);
				s_o += (pchl + 1);
				f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

				f_c += 2;

				if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
				if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

				pch = strtok_r(NULL, "DMIS", &saveptr);
			}

			strncpy(b_cigar, cigarBuf2, f_cigarn);
			pch = strtok(cigarBuf2,"DMIS");

			if(pch != NULL)
				pchl = strlen(pch);

			snt = 0;

#ifdef	CIGAR_S_MODIFY
			if(lv_up_left)
			{
				sn = sprintf(cigar_p1 + snt, "%uS", lv_up_left);
				snt += sn;
			}
#else
			if(m_n_f)
			{
				if(f_c)
				{
					if(f_cigar[f_cigarn + 1 - f_c] == 'M')
					{
						f_cigar[f_cigarn - f_c] += m_n_f;
					}
					else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
					{
						f_cigar[f_cigarn - f_c] += m_n_f;
					}
					else
					{
						sn = sprintf(cigar_p1 + snt, "%uM", m_n_f);
						snt += sn;
					}
				}
				else	m_m_n += m_n_f;

			}
#endif
			if((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))   //(op_dm_l1[tid][v_cnt_i] != -1) && (op_dm_r1[tid][v_cnt_i] != read_length1)
			{
				if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
				{
					f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

					for(f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
					sn = sprintf(cigar_p1 + snt, "%s", b_cigar + pchl + 1);
					snt += sn;
				}
				else if(f_cigar[f_cigarn - 1] == 'M')
				{
					f_cigar[f_cigarn - 2] += m_m_n;
					for(f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
					sn = sprintf(cigar_p1 + snt, "%s",b_cigar);
					snt += sn;
				}
				else if(b_cigar[pchl] == 'M')
				{
					for(f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}

					sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
					snt += sn;
				}
				else
				{
					for(f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
					sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
					snt += sn;
				}
			}
			else if((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
			{
				if(b_cigar[pchl] == 'M')
				{
					sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
					snt += sn;
				}
				else
				{
					sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
					snt += sn;
				}
			}
			else if((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
			{
				if(f_cigar[f_cigarn - 1] == 'M')
				{
					f_cigar[f_cigarn - 2] += m_m_n;
					for(f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
				}
				else
				{
					for(f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
					sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
					snt += sn;
				}
			}
			else
			{
				sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
				snt += sn;
			}
#ifdef	CIGAR_S_MODIFY
			if(lv_down_right < read_length1)
			{
				sn = sprintf(cigar_p1 + snt, "%uS", read_length1 - lv_down_right);
				snt += sn;
			}
#else
			if(m_n_b)
			{
				if(cigar_p1[snt - 1] == 'M')
				{
					for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
					{
						if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
						m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
					}
					sn = sprintf(cigar_p1 + bit_char_i + 1, "%uM", m_n_b);
					snt = bit_char_i + 1 + sn;
				}
				else if(cigar_p1[snt - 1] == 'S')
				{
					for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
					{
						if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
						m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
					}
					sn = sprintf(cigar_p1 + bit_char_i + 1, "%uS", m_n_b);
					snt = bit_char_i + 1 + sn;
				}
				else
				{
					sn = sprintf(cigar_p1 + snt, "%uM", m_n_b);
					snt += sn;
				}
			}
#endif

			//sprintf(cigar_p1 + snt, "\0");
		}

#ifdef	CIGAR_LEN_ERR
		cigar_len = 0;
		s_o_tmp = 0;
		strncpy(cigar_tmp, cigar_p1, snt);
		cigar_tmp[snt] = '\0';
		pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

		while (pch_tmp != NULL)
		{
			pchl_tmp = strlen(pch_tmp);
			s_o_tmp += (pchl_tmp + 1);

			if(cigar_p1[s_o_tmp - 1] != 'D')
			{
				cigar_len_tmp = atoi(pch_tmp);
				cigar_len += cigar_len_tmp;
			}
			pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
		}

		if(read_length1 != cigar_len)
		{
			if(read_length1 < cigar_len)
			{
				cigar_len_re = cigar_len_tmp - (cigar_len - read_length1);
				if(cigar_len_re > 0)	sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
				else if(cigar_len_re == 0)	sprintf(cigar_p1 + snt - sn, "\0");
				else	strcpy(cigar_p1, cigar_m1[tid]);
			}
			else
			{
				cigar_len_re = cigar_len_tmp + (read_length1 - cigar_len);
				sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
			}
		}
#endif
	}

	d_n2 = 0;
	i_n2 = 0;
	s_r_o_l = op_dm_l2[tid][v_cnt_i];
	s_r_o_r = op_dm_r2[tid][v_cnt_i];

	if((s_r_o_l == 0) && (s_r_o_r == 0))
	{
		strcpy(cigar_p2, cigar_m2[tid]);
#ifdef	MAPPING_QUALITY

		if(mp_flag)
		{
			for (bit_char_i = 0, read_b_i = 32; bit_char_i < read_length2; bit_char_i++, read_b_i++)
				if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((op_vector_seq2[tid][v_cnt_i][read_b_i >> 5] >> ((31 - (read_b_i & 0X1f)) << 1)) & 0X3))
				{
					sub_t[tid] += mp_subs_2[bit_char_i];
				}
		}

#endif
	}
	else     //indel
	{
		if((cir_n == cir_fix_n) && (local_ksw))
		{
#ifdef	KSW_ALN_PAIR

#ifdef	CHAR_CP
			for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl2[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
				read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
			for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl2[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
				read_char[tid][read_b_i] = charToDna5n[sam_seq2[bit_char_i]];
#endif

			for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl2[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i > -1
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			ksw_global(op_dm_kl2[tid][v_cnt_i], read_char[tid], op_dm_kl2[tid][v_cnt_i], ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar1, &cigar1);
#ifdef	CHAR_CP
			for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr2[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
				read_char2[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
			for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr2[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
				read_char2[tid][read_b_i] = charToDna5n[sam_seq2[bit_char_i]];
#endif

			for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr2[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length + 64
				ali_ref_seq2[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			ksw_global(op_dm_kr2[tid][v_cnt_i], read_char2[tid], op_dm_kr2[tid][v_cnt_i], ali_ref_seq2[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar2, &cigar2);

			m_m_n = s_r_o_r - s_r_o_l - 1;

			nm_score = 0;
			snt = 0;
			if (n_cigar1)
			{
				op_score1  = cigar1[0]&0xf;
				len_score1 = cigar1[0]>>4;
			}
			else	op_score1 = 3;

			if (n_cigar2)
			{
				op_score2  = cigar2[0]&0xf;
				len_score2 = cigar2[0]>>4;
			}
			else	op_score2 = 3;


			if (s_r_o_l >= op_dm_kl2[tid][v_cnt_i])
			{
				sn = sprintf(cigar_p2 + snt, "%dS", s_r_o_l + 1 - op_dm_kl2[tid][v_cnt_i]);
				snt += sn;
			}

			if ((s_r_o_l != -1) && (s_r_o_r != read_length2))
			{
				if ((op_score1 == 0) && (op_score2 == 0))
				{
					k_start1 = 0;
					k_start2 = 1;
					k_middle = len_score1 + len_score2 + m_m_n;
				}
				else if (op_score1 == 0)
				{
					k_start1 = 0;
					k_start2 = 0;
					k_middle = len_score1 + m_m_n;
				}
				else if (op_score2 == 0)
				{
					k_start1 = -1;
					k_start2 = 1;
					k_middle = len_score2 + m_m_n;
				}
				else
				{
					k_start1 = -1;
					k_start2 = 0;
					k_middle = m_m_n;
				}

				x_score = y_score = 0;
				for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
				{
					op_score  = cigar1[bit_char_i]&0xf;
					len_score = cigar1[bit_char_i]>>4;

					sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
					snt += sn;

					if (op_score == 0)
					{
						// match
						for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
							if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
							{
#ifdef	MAPPING_QUALITY
								if(mp_flag)	sub_t[tid] += mp_subs_2_o[read_length2 + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
								++nm_score;
							}
						x_score += len_score;
						y_score += len_score;
					}
					else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
					else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
				}

				sn = sprintf(cigar_p2 + snt, "%dM", k_middle);
				snt += sn;

				x_score = y_score = 0;
				for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
				{
					op_score  = cigar2[bit_char_i]&0xf;
					len_score = cigar2[bit_char_i]>>4;

					sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
					snt += sn;

					if (op_score == 0)
					{
						// match
						for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
							if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
							{
#ifdef	MAPPING_QUALITY
								if(mp_flag)	sub_t[tid] += mp_subs_2[s_r_o_r + x_score + read_b_i];
#endif
								++nm_score;
							}
						x_score += len_score;
						y_score += len_score;
					}
					else if (op_score == 1) x_score += len_score, nm_score += len_score;
					else if (op_score == 2) y_score += len_score, nm_score += len_score;
				}
			}
			else if (s_r_o_l == -1)
			{
				if (op_score2 == 0)
				{
					k_start2 = 1;
					k_middle = len_score2 + m_m_n;
				}
				else
				{
					k_start2 = 0;
					k_middle = m_m_n;
				}
				sn = sprintf(cigar_p2 + snt, "%dM", k_middle);
				snt += sn;

				x_score = y_score = 0;
				for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
				{
					op_score  = cigar2[bit_char_i]&0xf;
					len_score = cigar2[bit_char_i]>>4;

					sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
					snt += sn;

					if (op_score == 0)
					{
						// match
						for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
							if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
							{
#ifdef	MAPPING_QUALITY
								if(mp_flag)	sub_t[tid] += mp_subs_2[s_r_o_r + x_score + read_b_i];
#endif
								++nm_score;
							}
						x_score += len_score;
						y_score += len_score;
					}
					else if (op_score == 1) x_score += len_score, nm_score += len_score;
					else if (op_score == 2) y_score += len_score, nm_score += len_score;
				}
			}
			else
			{
				if (op_score1 == 0)
				{
					k_start1 = 0;
					k_middle = len_score1 + m_m_n;
				}
				else
				{
					k_start1 = -1;
					k_middle = m_m_n;
				}
				x_score = y_score = 0;
				for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
				{
					op_score  = cigar1[bit_char_i]&0xf;
					len_score = cigar1[bit_char_i]>>4;

					sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
					snt += sn;

					if (op_score == 0)
					{
						// match
						for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
							if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
							{
#ifdef	MAPPING_QUALITY
								if(mp_flag)	sub_t[tid] += mp_subs_2_o[read_length2 + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
								++nm_score;
							}
						x_score += len_score;
						y_score += len_score;
					}
					else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
					else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
				}

				sn = sprintf(cigar_p2 + snt, "%dM", k_middle);
				snt += sn;
			}

			if (read_length2 - s_r_o_r > op_dm_kr2[tid][v_cnt_i])
			{
				sn = sprintf(cigar_p2 + snt, "%dS", read_length2 - s_r_o_r - op_dm_kr2[tid][v_cnt_i]);
				snt += sn;
			}
			//sn = sprintf(cigar_p2 + snt, "\0");
			//snt += sn;

			if (n_cigar1)	free(cigar1);
			if (n_cigar2)	free(cigar2);

			op_dm_ex2[tid][v_cnt_i] = nm_score;
#endif

		}
		else
		{
#ifdef	OUTPUT_DEBUG

			if (pound_pos_2_f >= s_r_o_r)   //1
			{
				lv_up_left = 0;
				lv_up_right = s_r_o_l;
				lv_down_right = pound_pos_2_f;
				lv_down_left = s_r_o_r;
				m_n_f = 0;
				m_n_b = read_length2 - pound_pos_2_f;
				m_m_n = s_r_o_r - s_r_o_l - 1;
			}
			else if (pound_pos_2_r <= s_r_o_l + 1)     //5
			{
				lv_up_left = pound_pos_2_r;//21
				lv_up_right = s_r_o_l;//20
				lv_down_right = read_length2;
				lv_down_left = s_r_o_r;
				m_n_f = pound_pos_2_r;
				m_n_b = 0;
				m_m_n = s_r_o_r - s_r_o_l - 1;
			}
			else if ((pound_pos_2_f <= s_r_o_l + 1) && (pound_pos_2_r >= s_r_o_r))     //2
			{
				lv_up_left = 0;
				lv_up_right = pound_pos_2_f - 1;
				lv_down_right = read_length2;
				lv_down_left = pound_pos_2_r;
				m_n_f = 0;
				m_n_b = 0;
				m_m_n = pound_pos_2_r - pound_pos_2_f;
			}
			else if ((pound_pos_2_f > s_r_o_l + 1) && (pound_pos_2_f < s_r_o_r))     //3
			{
				lv_up_left = 0;
				lv_up_right = s_r_o_l;
				lv_down_right = read_length2;
				lv_down_left = pound_pos_2_r;
				m_n_f = 0;
				m_n_b = 0;
				m_m_n = read_length2 - s_r_o_l - 1;
			}
			else     //4
			{
				lv_up_left = 0;
				lv_up_right = -1;
				lv_down_right = read_length2;
				lv_down_left = s_r_o_r;
				m_n_f = 0;
				m_n_b = 0;
				m_m_n = s_r_o_r;
			}

#ifdef	QUAL_FILT_LV_OUT

#ifdef	CHAR_CP
			for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
				read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
			for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
				read_char[tid][read_b_i] = sam_seq2[bit_char_i];

			for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
				ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY
			if(mp_flag)
				computeEditDistanceWithCigar_s_mis_left_mp(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_2_o + read_length2 - 1 - lv_up_right, &s_offset2, mp_subs_2_o + read_length2 - 1 - lv_up_right, &(sub_t[tid]));//, 0, op_dm_sl1[tid][v_cnt_i]
			else
				computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_2_o + read_length2 - 1 - lv_up_right, &s_offset2);//, 0, op_dm_sl1[tid][v_cnt_i]
#else
			computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_2_o + read_length2 - 1 - lv_up_right, &s_offset2);//, 0, op_dm_sl1[tid][v_cnt_i]
#endif

#ifdef	CHAR_CP
			for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
				read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
			for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
				read_char[tid][read_b_i] = sam_seq2[bit_char_i];

			for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
				ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

#ifdef	MAPPING_QUALITY
			if(mp_flag)
				computeEditDistanceWithCigar_s_mis_mp(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_2 + lv_down_left, mp_subs_2 + lv_down_left, &(sub_t[tid]));//, 0, op_dm_sr1[tid][v_cnt_i]
			else	computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_2 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]
#else
			computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_2 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#else

#ifdef CHAR_CP
			for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
				read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
			for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
				read_char[tid][read_b_i] = sam_seq2[bit_char_i];

			for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
				ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY
			if(mp_flag)
				computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid], mp_subs_2_o + read_length2 - 1 - lv_up_right, &(sub_t[tid]));//, 0, op_dm_sl1[tid][v_cnt_i]
			else
				computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]

#else
			computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#endif

#ifdef	CHAR_CP
			for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
				read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

			for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
				ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
			for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
				read_char[tid][read_b_i] = sam_seq2[bit_char_i];

			for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
				ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY
			if(mp_flag)
				computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid], mp_subs_2 + lv_down_left, &(sub_t[tid]));//, 0, op_dm_sr1[tid][v_cnt_i]
			else	computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#else
			computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#endif

#endif

			//deal with front and back lv cigar
			strncpy(str_o, cigarBuf1, f_cigarn);
			s_o = 0;
			f_c = 0;

			pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

			while (pch != NULL)
			{
				pchl = strlen(pch);
				f_cigar[f_cigarn - f_c - 2] = atoi(pch);
				s_o += (pchl + 1);
				f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

				f_c += 2;

				if (str_o[s_o - 1] == 'D')	d_n2 += atoi(pch);
				if (str_o[s_o - 1] == 'I')	i_n2 += atoi(pch);

				pch = strtok_r(NULL, "DMIS", &saveptr);
			}

			strncpy(b_cigar, cigarBuf2, f_cigarn);
			pch = strtok(cigarBuf2,"DMIS");

			if (pch != NULL)
				pchl = strlen(pch);

			snt = 0;

#ifdef	CIGAR_S_MODIFY
			if(lv_up_left)
			{
				sn = sprintf(cigar_p2 + snt, "%uS", lv_up_left);
				snt += sn;
			}
#else
			if (m_n_f)
			{
				if(f_c)
				{
					if (f_cigar[f_cigarn + 1 - f_c] == 'M')
					{
						f_cigar[f_cigarn - f_c] += m_n_f;
					}
					else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
					{
						f_cigar[f_cigarn - f_c] += m_n_f;
					}
					else
					{
						sn = sprintf(cigar_p2 + snt, "%uM", m_n_f);
						snt += sn;
					}
				}
				else	m_m_n += m_n_f;
			}
#endif

			if ((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))
			{
				if ((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
				{
					f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

					for (f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
					sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
					snt += sn;
				}
				else if (f_cigar[f_cigarn - 1] == 'M')
				{
					f_cigar[f_cigarn - 2] += m_m_n;
					for (f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
					sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
					snt += sn;
				}
				else if (b_cigar[pchl] == 'M')
				{
					for (f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}

					sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
					snt += sn;
				}
				else
				{
					for (f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
					sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
					snt += sn;
				}
			}
			else if ((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
			{
				if (b_cigar[pchl] == 'M')
				{
					sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
					snt += sn;
				}
				else
				{
					sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
					snt += sn;
				}
			}
			else if ((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
			{
				if (f_cigar[f_cigarn - 1] == 'M')
				{
					f_cigar[f_cigarn - 2] += m_m_n;
					for (f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
				}
				else
				{
					for (f_i = 0; f_i < f_c; f_i += 2)
					{
						sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
						snt += sn;
					}
					sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
					snt += sn;
				}
			}
			else
			{
				sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
				snt += sn;
			}
#ifdef	CIGAR_S_MODIFY
			if(lv_down_right < read_length2)
			{
				sn = sprintf(cigar_p2 + snt, "%uS", read_length2 - lv_down_right);
				snt += sn;
			}
#else
			if (m_n_b)
			{
				if (cigar_p2[snt - 1] == 'M')
				{
					for (bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
					{
						if ((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
						m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
					}
					sn = sprintf(cigar_p2 + bit_char_i + 1, "%uM", m_n_b);
					snt = bit_char_i + 1 + sn;
				}
				else if(cigar_p2[snt - 1] == 'S')
				{
					for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
					{
						if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
						m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
					}
					sn = sprintf(cigar_p2 + bit_char_i + 1, "%uS", m_n_b);
					snt = bit_char_i + 1 + sn;
				}
				else
				{
					sn = sprintf(cigar_p2 + snt, "%uM", m_n_b);
					snt += sn;
				}
			}
#endif
			//sprintf(cigar_p2 + snt, "\0");
		}


#ifdef	CIGAR_LEN_ERR
		cigar_len = 0;
		s_o_tmp = 0;
		strncpy(cigar_tmp, cigar_p2, snt);
		cigar_tmp[snt] = '\0';
		pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

		while (pch_tmp != NULL)
		{
			pchl_tmp = strlen(pch_tmp);
			s_o_tmp += (pchl_tmp + 1);

			if(cigar_p2[s_o_tmp - 1] != 'D')
			{
				cigar_len_tmp = atoi(pch_tmp);
				cigar_len += cigar_len_tmp;
			}

			pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
		}

		if(read_length2 != cigar_len)
		{
			if(read_length2 < cigar_len)
			{
				cigar_len_re = cigar_len_tmp - (cigar_len - read_length2);
				if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
				else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
				else	strcpy(cigar_p2, cigar_m2[tid]);
			}
			else
			{
				cigar_len_re = cigar_len_tmp + (read_length2 - cigar_len);
				sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
			}
		}
#endif

	}
#ifdef	MAPPING_QUALITY
	if(mp_flag)
	{
		sam_qual1 = sub_t[tid];
	}
#endif


	sam_pos1_pr = sam_pos1 + i_n1 - d_n1 + s_offset1;
	sam_pos2_pr = sam_pos2 + i_n2 - d_n2 + s_offset2;

	chr_re_tmp = chr_end_n[chr_re] - chr_end_n[chr_re - 1];

	if(sam_pos1_pr <= 0)
	{
		sam_pos1_pr = 1;
	}
	else
	{
		
		if((sam_pos1_pr + read_length1 - 1) > chr_re_tmp)
			sam_pos1_pr = chr_re_tmp - 1 - read_length1;
		
	}
	if(sam_pos2_pr <= 0)
	{
		sam_pos2_pr = 1;
	}
	else
	{
		if((sam_pos2_pr + read_length2 - 1) > chr_re_tmp)
			sam_pos2_pr = chr_re_tmp - 1 - read_length2;
	}

#ifdef	OUTPUT_ARR
	seqio[seqi].flag1 = sam_flag1;
	seqio[seqi].flag2 = sam_flag2;
	//seqio[seqi].chr_re = chr_re;
	seqio[seqi].chr_re1 = chr_re;
	seqio[seqi].chr_re2 = chr_re;
	seqio[seqi].pos1 = sam_pos1_pr;
	seqio[seqi].pos2 = sam_pos2_pr;
	seqio[seqi].cross = sam_cross;
	seqio[seqi].nm1 = op_dm_ex1[tid][v_cnt_i];
	seqio[seqi].nm2 = op_dm_ex2[tid][v_cnt_i];

	strcpy(pr_cigar1_buffer[seqi], cigar_p1);
	seqio[seqi].cigar1 = pr_cigar1_buffer[seqi];
	strcpy(pr_cigar2_buffer[seqi], cigar_p2);
	seqio[seqi].cigar2 = pr_cigar2_buffer[seqi];

	if(sam_flag1 == 99)
	{
#ifdef	CHAR_CP
		for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
			sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

		sam_seq2[sam_seq_i] = '\0';

		strrev1(sam_seq2);
#endif
		strcpy(read_rev_buffer[seqi], sam_seq2);
		read_rev_buffer[seqi][read_length2] = '\0';

		seqio[seqi].seq2 = read_rev_buffer[seqi];
		seqio[seqi].seq1 = seqio[seqi].read_seq1;

		strrev1(qual2_buffer[seqi]);
		seqio[seqi].qual1 = qual1_buffer[seqi];
		seqio[seqi].qual2 = qual2_buffer[seqi];
	}
	else
	{
#ifdef	CHAR_CP
		for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
			sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

		sam_seq1[sam_seq_i] = '\0';

		strrev1(sam_seq1);
#endif
		strcpy(read_rev_buffer[seqi], sam_seq1);
		read_rev_buffer[seqi][read_length1] = '\0';

		seqio[seqi].seq1 = read_rev_buffer[seqi];
		seqio[seqi].seq2 = seqio[seqi].read_seq2;

		strrev1(qual1_buffer[seqi]);
		seqio[seqi].qual1 = qual1_buffer[seqi];
		seqio[seqi].qual2 = qual2_buffer[seqi];
	}
#endif

	//}

	xa_i = 0;

	for(va_cnt_i = 1; va_cnt_i < v_cnt_out; va_cnt_i++)
	{
#ifdef	ALTER_DEBUG
		if(rep_go[tid])	v_cnt_i = seed_length_arr[tid][va_cnt_i].index;
		else	v_cnt_i = va_cnt_i;
#else
		v_cnt_i = va_cnt_i;
#endif
		x = op_vector_pos1[tid][v_cnt_i];
		low = 0;
		high = chr_file_n - 1;

		while ( low <= high )
		{
			mid = (low + high) >> 1;
			if(x < (chr_end_n[mid]))
			{
				high = mid - 1;
			}
			else if(x > (chr_end_n[mid]))
			{
				low = mid + 1;
			}
			else
			{
				chr_re =  mid;
				break;
			}
			chr_re = low;
		}

		sam_pos1 = op_vector_pos1[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;
		sam_pos2 = op_vector_pos2[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;

		chr_res[tid][xa_i] = chr_re;

		if(op_rc[tid][v_cnt_i] == 0)
		{
#ifdef OUPUT_REPEAT
			sam_flag1 = 99;
			sam_flag2 = 147;

#ifdef	CHAR_CP
			read_bit_1[tid] = read_bit1[tid][0];
			read_bit_2[tid] = read_bit2[tid][1];
#else
			for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
				sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

			sam_seq2[sam_seq_i] = '\0';
			strrev1(sam_seq2);

			strcpy(sam_seq1, seqio[seqi].read_seq1);
#endif

#endif


#ifdef	QUAL_FILT_LV_OUT
			qual_filt_lv_1 = qual_filt_lv1[tid][0];
			qual_filt_lv_1_o = qual_filt_lv1[tid][1];
			qual_filt_lv_2 = qual_filt_lv2[tid][1];
			qual_filt_lv_2_o = qual_filt_lv2[tid][0];
#endif

			sam_cross = sam_pos2 + read_length2 - sam_pos1;

			xa_d1s[tid][xa_i] = '+';
			xa_d2s[tid][xa_i] = '-';

			pound_pos_1_f = pound_pos1_f_forward;
			pound_pos_1_r = pound_pos1_r_forward;
			pound_pos_2_f = pound_pos2_f_reverse;
			pound_pos_2_r = pound_pos2_r_reverse;
		}
		else
		{
#ifdef OUPUT_REPEAT
			sam_flag1 = 83;
			sam_flag2 = 163;

#ifdef	CHAR_CP
			read_bit_1[tid] = read_bit1[tid][1];
			read_bit_2[tid] = read_bit2[tid][0];
#else
			for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
				sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

			sam_seq1[sam_seq_i] = '\0';

			strrev1(sam_seq1);
			strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#endif


#ifdef	QUAL_FILT_LV_OUT
			qual_filt_lv_1 = qual_filt_lv1[tid][1];
			qual_filt_lv_1_o = qual_filt_lv1[tid][0];
			qual_filt_lv_2 = qual_filt_lv2[tid][0];
			qual_filt_lv_2_o = qual_filt_lv2[tid][1];
#endif

			sam_cross = sam_pos2 - read_length1 - sam_pos1;

			xa_d1s[tid][xa_i] = '-';
			xa_d2s[tid][xa_i] = '+';

			pound_pos_1_f = pound_pos1_f_reverse;
			pound_pos_1_r = pound_pos1_r_reverse;
			pound_pos_2_f = pound_pos2_f_forward;
			pound_pos_2_r = pound_pos2_r_forward;
		}

		d_n1 = 0;
		i_n1 = 0;
		s_offset1 = 0;
		s_offset2 = 0;
		s_r_o_l = op_dm_l1[tid][v_cnt_i];
		s_r_o_r = op_dm_r1[tid][v_cnt_i];
		if((s_r_o_l == 0) && (s_r_o_r == 0))
		{
			strcpy(cigar_p1, cigar_m1[tid]);

			lv_re1 = op_dm_ex1[tid][v_cnt_i];
		}
		else     //indel
		{
			if((cir_n == cir_fix_n) && (local_ksw))
			{
#ifdef	KSW_ALN_PAIR

#ifdef	CHAR_CP
				for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
					read_char[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif

				for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i > -1
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				ksw_global(op_dm_kl1[tid][v_cnt_i], read_char[tid], op_dm_kl1[tid][v_cnt_i], ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar1, &cigar1);

#ifdef	CHAR_CP
				for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
					read_char2[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
					read_char2[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif

				for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length + 64
					ali_ref_seq2[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				ksw_global(op_dm_kr1[tid][v_cnt_i], read_char2[tid], op_dm_kr1[tid][v_cnt_i], ali_ref_seq2[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar2, &cigar2);

				m_m_n = s_r_o_r - s_r_o_l - 1;

				nm_score = 0;
				snt = 0;
				if (n_cigar1)
				{
					op_score1  = cigar1[0]&0xf;
					len_score1 = cigar1[0]>>4;
				}
				else	op_score1 = 3;

				if (n_cigar2)
				{
					op_score2  = cigar2[0]&0xf;
					len_score2 = cigar2[0]>>4;
				}
				else	op_score2 = 3;

				if (s_r_o_l >= op_dm_kl1[tid][v_cnt_i])
				{
					sn = sprintf(cigar_p1 + snt, "%dS", s_r_o_l + 1 - op_dm_kl1[tid][v_cnt_i]);
					snt += sn;
				}

				if ((s_r_o_l != -1) && (s_r_o_r != read_length1))
				{
					if ((op_score1 == 0) && (op_score2 == 0))
					{
						k_start1 = 0;
						k_start2 = 1;
						k_middle = len_score1 + len_score2 + m_m_n;
					}
					else if (op_score1 == 0)
					{
						k_start1 = 0;
						k_start2 = 0;
						k_middle = len_score1 + m_m_n;
					}
					else if (op_score2 == 0)
					{
						k_start1 = -1;
						k_start2 = 1;
						k_middle = len_score2 + m_m_n;
					}
					else
					{
						k_start1 = -1;
						k_start2 = 0;
						k_middle = m_m_n;
					}

					x_score = y_score = 0;
					for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
					{
						op_score  = cigar1[bit_char_i]&0xf;
						len_score = cigar1[bit_char_i]>>4;

						sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i]) ++nm_score;
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
					}

					sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
					snt += sn;

					x_score = y_score = 0;
					for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
					{
						op_score  = cigar2[bit_char_i]&0xf;
						len_score = cigar2[bit_char_i]>>4;

						sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i]) ++nm_score;
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score;
					}
				}
				else if (s_r_o_l == -1)
				{
					if (op_score2 == 0)
					{
						k_start2 = 1;
						k_middle = len_score2 + m_m_n;
					}
					else
					{
						k_start2 = 0;
						k_middle = m_m_n;
					}
					sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
					snt += sn;

					x_score = y_score = 0;
					for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
					{
						op_score  = cigar2[bit_char_i]&0xf;
						len_score = cigar2[bit_char_i]>>4;

						sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i]) ++nm_score;
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score;
					}
				}
				else
				{
					if (op_score1 == 0)
					{
						k_start1 = 0;
						k_middle = len_score1 + m_m_n;
					}
					else
					{
						k_start1 = -1;
						k_middle = m_m_n;
					}
					x_score = y_score = 0;
					for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
					{
						op_score  = cigar1[bit_char_i]&0xf;
						len_score = cigar1[bit_char_i]>>4;

						sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i]) ++nm_score;
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
					}

					sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
					snt += sn;
				}

				if (read_length1 - s_r_o_r > op_dm_kr1[tid][v_cnt_i])
				{
					sn = sprintf(cigar_p1 + snt, "%dS", read_length1 - s_r_o_r - op_dm_kr1[tid][v_cnt_i]);
					snt += sn;
				}
				//sn = sprintf(cigar_p1 + snt, "\0");
				//snt += sn;

				if (n_cigar1)	free(cigar1);
				if (n_cigar2)	free(cigar2);

				lv_re1 = nm_score;
#endif

			}
			else
			{
#ifdef	OUTPUT_DEBUG
				if (pound_pos_1_f >= s_r_o_r)   //1
				{
					lv_up_left = 0;
					lv_up_right = s_r_o_l;
					lv_down_right = pound_pos_1_f;
					lv_down_left = s_r_o_r;
					m_n_f = 0;
					m_n_b = read_length1 - pound_pos_1_f;
					m_m_n = s_r_o_r - s_r_o_l - 1;
				}
				else if (pound_pos_1_r <= s_r_o_l + 1)     //5
				{
					lv_up_left = pound_pos_1_r;//
					lv_up_right = s_r_o_l;
					lv_down_right = read_length1;
					lv_down_left = s_r_o_r;
					m_n_f = pound_pos_1_r;
					m_n_b = 0;
					m_m_n = s_r_o_r - s_r_o_l - 1;
				}
				else if ((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
				{
					lv_up_left = 0;
					lv_up_right = pound_pos_1_f - 1;
					lv_down_right = read_length1;
					lv_down_left = pound_pos_1_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = pound_pos_1_r - pound_pos_1_f;
				}
				else if ((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
				{
					lv_up_left = 0;
					lv_up_right = s_r_o_l;
					lv_down_right = read_length1;
					lv_down_left = pound_pos_1_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = read_length1 - s_r_o_l - 1;
				}
				else     //4
				{
					lv_up_left = 0;
					lv_up_right = -1;
					lv_down_right = read_length1;
					lv_down_left = s_r_o_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = s_r_o_r;
				}

#ifdef	QUAL_FILT_LV_OUT

#ifdef	CHAR_CP
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = sam_seq1[bit_char_i];

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

				lv_re1f = computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length1 - 1 - lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]

#ifdef	CHAR_CP
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = sam_seq1[bit_char_i];

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif


				lv_re1b = computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = sam_seq1[bit_char_i];

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

				lv_re1f = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#ifdef	CHAR_CP
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = sam_seq1[bit_char_i];

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

				lv_re1b = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#endif
				lv_re1 = lv_re1f + lv_re1b;

				//deal with front and back lv cigar
				strncpy(str_o, cigarBuf1, f_cigarn);
				s_o = 0;
				f_c = 0;

				pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

				while (pch != NULL)
				{
					pchl = strlen(pch);
					f_cigar[f_cigarn - f_c - 2] = atoi(pch);
					s_o += (pchl + 1);
					f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

					f_c += 2;

					if (str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
					if (str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

					pch = strtok_r(NULL, "DMIS", &saveptr);
				}

				strncpy(b_cigar, cigarBuf2, f_cigarn);
				pch = strtok(cigarBuf2,"DMIS");

				if (pch != NULL)
					pchl = strlen(pch);

				snt = 0;
#ifdef	CIGAR_S_MODIFY
				if(lv_up_left)
				{
					sn = sprintf(cigar_p1 + snt, "%uS", lv_up_left);
					snt += sn;
				}
#else
				if (m_n_f)
				{
					if(f_c)
					{
						if (f_cigar[f_cigarn + 1 - f_c] == 'M')
						{
							f_cigar[f_cigarn - f_c] += m_n_f;
						}
						else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
						{
							f_cigar[f_cigarn - f_c] += m_n_f;
						}
						else
						{
							sn = sprintf(cigar_p1 + snt, "%uM", m_n_f);
							snt += sn;
						}
					}
					else	m_m_n += m_n_f;

				}
#endif
				if ((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))
				{
					if ((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
					{
						f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p1 + snt, "%s", b_cigar + pchl + 1);
						snt += sn;
					}
					else if (f_cigar[f_cigarn - 1] == 'M')
					{
						f_cigar[f_cigarn - 2] += m_m_n;
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p1 + snt, "%s",b_cigar);
						snt += sn;
					}
					else if (b_cigar[pchl] == 'M')
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}

						sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						snt += sn;
					}
					else
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
				}
				else if ((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
				{
					if (b_cigar[pchl] == 'M')
					{
						sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						snt += sn;
					}
					else
					{
						sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
				}
				else if ((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
				{
					if (f_cigar[f_cigarn - 1] == 'M')
					{
						f_cigar[f_cigarn - 2] += m_m_n;
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
					}
					else
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
						snt += sn;
					}
				}
				else
				{
					sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
					snt += sn;
				}
#ifdef	CIGAR_S_MODIFY
				if(lv_down_right < read_length1)
				{
					sn = sprintf(cigar_p1 + snt, "%uS", read_length1 - lv_down_right);
					snt += sn;
				}
#else
				if (m_n_b)
				{
					if (cigar_p1[snt - 1] == 'M')
					{
						for (bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
						{
							if ((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
							m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
						}
						sn = sprintf(cigar_p1 + bit_char_i + 1, "%uM", m_n_b);
						snt = bit_char_i + 1 + sn;
					}
					else if(cigar_p1[snt - 1] == 'S')
					{
						for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
						{
							if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
							m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
						}
						sn = sprintf(cigar_p1 + bit_char_i + 1, "%uS", m_n_b);
						snt = bit_char_i + 1 + sn;
					}
					else
					{
						sn = sprintf(cigar_p1 + snt, "%uM", m_n_b);
						snt += sn;
					}
				}
#endif

				//sprintf(cigar_p1 + snt, "\0");
			}
#ifdef	CIGAR_LEN_ERR
			cigar_len = 0;
			s_o_tmp = 0;
			strncpy(cigar_tmp, cigar_p1, snt);
			cigar_tmp[snt] = '\0';
			pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

			while (pch_tmp != NULL)
			{
				pchl_tmp = strlen(pch_tmp);
				s_o_tmp += (pchl_tmp + 1);

				if(cigar_p1[s_o_tmp - 1] != 'D')
				{
					cigar_len_tmp = atoi(pch_tmp);
					cigar_len += cigar_len_tmp;
				}

				pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
			}

			if(read_length1 != cigar_len)
			{
				if(read_length1 < cigar_len)
				{
					cigar_len_re = cigar_len_tmp - (cigar_len - read_length1);
					if(cigar_len_re > 0)	sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
					else if(cigar_len_re == 0)	sprintf(cigar_p1 + snt - sn, "\0");
					else	strcpy(cigar_p1, cigar_m1[tid]);
				}
				else
				{
					cigar_len_re = cigar_len_tmp + (read_length1 - cigar_len);
					sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
				}
			}
#endif
		}

		d_n2 = 0;
		i_n2 = 0;
		s_r_o_l = op_dm_l2[tid][v_cnt_i];
		s_r_o_r = op_dm_r2[tid][v_cnt_i];

		if((s_r_o_l == 0) && (s_r_o_r == 0))
		{
			strcpy(cigar_p2, cigar_m2[tid]);

			lv_re2 = op_dm_ex2[tid][v_cnt_i];
		}
		else     //indel
		{
			if((cir_n == cir_fix_n) && (local_ksw))
			{
#ifdef	KSW_ALN_PAIR

#ifdef	CHAR_CP
				for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl2[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl2[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
					read_char[tid][read_b_i] = charToDna5n[sam_seq2[bit_char_i]];
#endif

				for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl2[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i > -1
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				ksw_global(op_dm_kl2[tid][v_cnt_i], read_char[tid], op_dm_kl2[tid][v_cnt_i], ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar1, &cigar1);

#ifdef CHAR_CP
				for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr2[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
					read_char2[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr2[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
					read_char2[tid][read_b_i] = charToDna5n[sam_seq2[bit_char_i]];
#endif

				for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr2[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length + 64
					ali_ref_seq2[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				ksw_global(op_dm_kr2[tid][v_cnt_i], read_char2[tid], op_dm_kr2[tid][v_cnt_i], ali_ref_seq2[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar2, &cigar2);

				m_m_n = s_r_o_r - s_r_o_l - 1;

				nm_score = 0;
				snt = 0;
				if (n_cigar1)
				{
					op_score1  = cigar1[0]&0xf;
					len_score1 = cigar1[0]>>4;
				}
				else	op_score1 = 3;

				if (n_cigar2)
				{
					op_score2  = cigar2[0]&0xf;
					len_score2 = cigar2[0]>>4;
				}
				else	op_score2 = 3;


				if (s_r_o_l >= op_dm_kl2[tid][v_cnt_i])
				{
					sn = sprintf(cigar_p2 + snt, "%dS", s_r_o_l + 1 - op_dm_kl2[tid][v_cnt_i]);
					snt += sn;
				}

				if ((s_r_o_l != -1) && (s_r_o_r != read_length2))
				{
					if ((op_score1 == 0) && (op_score2 == 0))
					{
						k_start1 = 0;
						k_start2 = 1;
						k_middle = len_score1 + len_score2 + m_m_n;
					}
					else if (op_score1 == 0)
					{
						k_start1 = 0;
						k_start2 = 0;
						k_middle = len_score1 + m_m_n;
					}
					else if (op_score2 == 0)
					{
						k_start1 = -1;
						k_start2 = 1;
						k_middle = len_score2 + m_m_n;
					}
					else
					{
						k_start1 = -1;
						k_start2 = 0;
						k_middle = m_m_n;
					}

					x_score = y_score = 0;
					for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
					{
						op_score  = cigar1[bit_char_i]&0xf;
						len_score = cigar1[bit_char_i]>>4;

						sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i]) ++nm_score;
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
					}

					sn = sprintf(cigar_p2 + snt, "%dM", k_middle);
					snt += sn;

					x_score = y_score = 0;
					for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
					{
						op_score  = cigar2[bit_char_i]&0xf;
						len_score = cigar2[bit_char_i]>>4;

						sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i]) ++nm_score;
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score;
					}
				}
				else if (s_r_o_l == -1)
				{
					if (op_score2 == 0)
					{
						k_start2 = 1;
						k_middle = len_score2 + m_m_n;
					}
					else
					{
						k_start2 = 0;
						k_middle = m_m_n;
					}
					sn = sprintf(cigar_p2 + snt, "%dM", k_middle);
					snt += sn;

					x_score = y_score = 0;
					for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
					{
						op_score  = cigar2[bit_char_i]&0xf;
						len_score = cigar2[bit_char_i]>>4;

						sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i]) ++nm_score;
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score;
					}
				}
				else
				{
					if (op_score1 == 0)
					{
						k_start1 = 0;
						k_middle = len_score1 + m_m_n;
					}
					else
					{
						k_start1 = -1;
						k_middle = m_m_n;
					}
					x_score = y_score = 0;
					for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
					{
						op_score  = cigar1[bit_char_i]&0xf;
						len_score = cigar1[bit_char_i]>>4;

						sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i]) ++nm_score;
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
					}

					sn = sprintf(cigar_p2 + snt, "%dM", k_middle);
					snt += sn;
				}

				if (read_length2 - s_r_o_r > op_dm_kr2[tid][v_cnt_i])
				{
					sn = sprintf(cigar_p2 + snt, "%dS", read_length2 - s_r_o_r - op_dm_kr2[tid][v_cnt_i]);
					snt += sn;
				}
				//sn = sprintf(cigar_p2 + snt, "\0");
				//snt += sn;

				if (n_cigar1)	free(cigar1);
				if (n_cigar2)	free(cigar2);

				lv_re2 = nm_score;
#endif

			}
			else
			{
#ifdef	OUTPUT_DEBUG
				if (pound_pos_2_f >= s_r_o_r)   //1
				{
					lv_up_left = 0;
					lv_up_right = s_r_o_l;
					lv_down_right = pound_pos_2_f;
					lv_down_left = s_r_o_r;
					m_n_f = 0;
					m_n_b = read_length2 - pound_pos_2_f;
					m_m_n = s_r_o_r - s_r_o_l - 1;
				}
				else if (pound_pos_2_r <= s_r_o_l + 1)     //5
				{
					lv_up_left = pound_pos_2_r;//
					lv_up_right = s_r_o_l;
					lv_down_right = read_length2;
					lv_down_left = s_r_o_r;
					m_n_f = pound_pos_2_r;
					m_n_b = 0;
					m_m_n = s_r_o_r - s_r_o_l - 1;
				}
				else if ((pound_pos_2_f <= s_r_o_l + 1) && (pound_pos_2_r >= s_r_o_r))     //2
				{
					lv_up_left = 0;
					lv_up_right = pound_pos_2_f - 1;
					lv_down_right = read_length2;
					lv_down_left = pound_pos_2_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = pound_pos_2_r - pound_pos_2_f;
				}
				else if ((pound_pos_2_f > s_r_o_l + 1) && (pound_pos_2_f < s_r_o_r))     //3
				{
					lv_up_left = 0;
					lv_up_right = s_r_o_l;
					lv_down_right = read_length2;
					lv_down_left = pound_pos_2_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = read_length2 - s_r_o_l - 1;
				}
				else     //4
				{
					lv_up_left = 0;
					lv_up_right = -1;
					lv_down_right = read_length2;
					lv_down_left = s_r_o_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = s_r_o_r;
				}

#ifdef	QUAL_FILT_LV_OUT

#ifdef	CHAR_CP
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = sam_seq2[bit_char_i];

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

				lv_re2f = computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_2_o + read_length2 - 1 - lv_up_right, &s_offset2);//, 0, op_dm_sl1[tid][v_cnt_i]

#ifdef CHAR_CP
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = sam_seq2[bit_char_i];

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

				lv_re2b = computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_2 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]

#else

#ifdef	CHAR_CP
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = sam_seq2[bit_char_i];

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

				lv_re2f = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#ifdef	CHAR_CP
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = sam_seq2[bit_char_i];

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

				lv_re2b = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#endif
				lv_re2 = lv_re2f + lv_re2b;

				//deal with front and back lv cigar
				strncpy(str_o, cigarBuf1, f_cigarn);
				s_o = 0;
				f_c = 0;

				pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

				while (pch != NULL)
				{
					pchl = strlen(pch);
					f_cigar[f_cigarn - f_c - 2] = atoi(pch);
					s_o += (pchl + 1);
					f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

					f_c += 2;

					if (str_o[s_o - 1] == 'D')	d_n2 += atoi(pch);
					if (str_o[s_o - 1] == 'I')	i_n2 += atoi(pch);

					pch = strtok_r(NULL, "DMIS", &saveptr);
				}

				strncpy(b_cigar, cigarBuf2, f_cigarn);
				pch = strtok(cigarBuf2,"DMIS");

				if (pch != NULL)
					pchl = strlen(pch);

				snt = 0;
#ifdef	CIGAR_S_MODIFY
				if(lv_up_left)
				{
					sn = sprintf(cigar_p2 + snt, "%uS", lv_up_left);
					snt += sn;
				}
#else
				if (m_n_f)
				{
					if(f_c)
					{
						if (f_cigar[f_cigarn + 1 - f_c] == 'M')
						{
							f_cigar[f_cigarn - f_c] += m_n_f;
						}
						else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
						{
							f_cigar[f_cigarn - f_c] += m_n_f;
						}
						else
						{
							sn = sprintf(cigar_p2 + snt, "%uM", m_n_f);
							snt += sn;
						}
					}
					else	m_m_n += m_n_f;

				}
#endif
				if ((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))
				{
					if ((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
					{
						f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
						snt += sn;
					}
					else if (f_cigar[f_cigarn - 1] == 'M')
					{
						f_cigar[f_cigarn - 2] += m_m_n;
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
						snt += sn;
					}
					else if (b_cigar[pchl] == 'M')
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}

						sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						snt += sn;
					}
					else
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
				}
				else if ((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
				{
					if (b_cigar[pchl] == 'M')
					{
						sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						snt += sn;
					}
					else
					{
						sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
				}
				else if ((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
				{
					if (f_cigar[f_cigarn - 1] == 'M')
					{
						f_cigar[f_cigarn - 2] += m_m_n;
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
					}
					else
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
						snt += sn;
					}
				}
				else
				{
					sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
					snt += sn;
				}
#ifdef	CIGAR_S_MODIFY
				if(lv_down_right < read_length2)
				{
					sn = sprintf(cigar_p2 + snt, "%uS", read_length2 - lv_down_right);
					snt += sn;
				}
#else
				if (m_n_b)
				{
					if (cigar_p2[snt - 1] == 'M')
					{
						for (bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
						{
							if ((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
							m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
						}
						sn = sprintf(cigar_p2 + bit_char_i + 1, "%uM", m_n_b);
						snt = bit_char_i + 1 + sn;
					}
					else if(cigar_p2[snt - 1] == 'S')
					{
						for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
						{
							if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
							m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
						}
						sn = sprintf(cigar_p2 + bit_char_i + 1, "%uS", m_n_b);
						snt = bit_char_i + 1 + sn;
					}
					else
					{
						sn = sprintf(cigar_p2 + snt, "%uM", m_n_b);
						snt += sn;
					}
				}
#endif

				//sprintf(cigar_p2 + snt, "\0");
			}
#ifdef	CIGAR_LEN_ERR
			cigar_len = 0;
			s_o_tmp = 0;
			strncpy(cigar_tmp, cigar_p2, snt);
			cigar_tmp[snt] = '\0';
			pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

			while (pch_tmp != NULL)
			{
				pchl_tmp = strlen(pch_tmp);
				s_o_tmp += (pchl_tmp + 1);

				if(cigar_p2[s_o_tmp - 1] != 'D')
				{
					cigar_len_tmp = atoi(pch_tmp);
					cigar_len += cigar_len_tmp;
				}

				pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
			}

			if(read_length2 != cigar_len)
			{
				if(read_length2 < cigar_len)
				{
					cigar_len_re = cigar_len_tmp - (cigar_len - read_length2);
					if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
					else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
					else	strcpy(cigar_p2, cigar_m2[tid]);
				}
				else
				{
					cigar_len_re = cigar_len_tmp + (read_length2 - cigar_len);
					sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
				}
			}
#endif

		}

		sam_pos1s[tid][xa_i] = (uint32_t )(sam_pos1 + i_n1 - d_n1 + s_offset1);
		sam_pos2s[tid][xa_i] = (uint32_t )(sam_pos2 + i_n2 - d_n2 + s_offset2);

		strcpy(cigar_p1s[tid][xa_i], cigar_p1);
		strcpy(cigar_p2s[tid][xa_i], cigar_p2);

		lv_re1s[tid][xa_i] = lv_re1;
		lv_re2s[tid][xa_i] = lv_re2;
		++xa_i;
	}

	//suboptimal
	vs_cnt_out = (cus_ali_n < vs_cnt) ? cus_ali_n:vs_cnt;

#ifdef	ALL_ALL
	for(v_cnt_i = 0; v_cnt_i < vs_cnt_out; v_cnt_i++)
#else
	if((mgn_flag) && (v_cnt_out > 1))
		ops_flag = 0;

	for(v_cnt_i = 0; (v_cnt_i < vs_cnt_out) && (dm_op[tid] + 3 > dm_ops[tid]) && ops_flag; v_cnt_i++)// && (v_cnt_out < 2)
#endif
	{
		x = ops_vector_pos1[tid][v_cnt_i];
		low = 0;
		high = chr_file_n - 1;

		while ( low <= high )
		{
			mid = (low + high) >> 1;
			if(x < (chr_end_n[mid]))
			{
				high = mid - 1;
			}
			else if(x > (chr_end_n[mid]))
			{
				low = mid + 1;
			}
			else
			{
				chr_re =  mid;
				break;
			}
			chr_re = low;
		}

		sam_pos1 = ops_vector_pos1[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;
		sam_pos2 = ops_vector_pos2[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;

		chr_res[tid][xa_i] = chr_re;

		if(ops_rc[tid][v_cnt_i] == 0)
		{
#ifdef OUPUT_REPEAT
			sam_flag1 = 99;
			sam_flag2 = 147;

#ifdef	CHAR_CP
			read_bit_1[tid] = read_bit1[tid][0];
			read_bit_2[tid] = read_bit2[tid][1];
#else

			for(sam_seq_i = 0; sam_seq_i < read_length2; sam_seq_i++)
				sam_seq2[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[sam_seq_i]] ^ 0X3];

			sam_seq2[sam_seq_i] = '\0';

			strrev1(sam_seq2);
			strcpy(sam_seq1, seqio[seqi].read_seq1);
#endif

#endif


#ifdef	QUAL_FILT_LV_OUT
			qual_filt_lv_1 = qual_filt_lv1[tid][0];
			qual_filt_lv_1_o = qual_filt_lv1[tid][1];
			qual_filt_lv_2 = qual_filt_lv2[tid][1];
			qual_filt_lv_2_o = qual_filt_lv2[tid][0];
#endif
			sam_cross = sam_pos2 + read_length2 - sam_pos1;

			xa_d1s[tid][xa_i] = '+';
			xa_d2s[tid][xa_i] = '-';

			pound_pos_1_f = pound_pos1_f_forward;
			pound_pos_1_r = pound_pos1_r_forward;
			pound_pos_2_f = pound_pos2_f_reverse;
			pound_pos_2_r = pound_pos2_r_reverse;

#ifdef	MAPPING_QUALITY
			if((mp_flag) && (v_cnt_i == 0))
			{
				mp_subs_1 = mp_subs1[tid][0];
				mp_subs_1_o = mp_subs1[tid][1];
				mp_subs_2 = mp_subs2[tid][1];
				mp_subs_2_o = mp_subs2[tid][0];
			}
#endif
		}
		else
		{
#ifdef OUPUT_REPEAT
			sam_flag1 = 83;
			sam_flag2 = 163;

#ifdef	CHAR_CP
			read_bit_1[tid] = read_bit1[tid][1];
			read_bit_2[tid] = read_bit2[tid][0];
#else
			for(sam_seq_i = 0; sam_seq_i < read_length1; sam_seq_i++)
				sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

			sam_seq1[sam_seq_i] = '\0';

			strrev1(sam_seq1);
			strcpy(sam_seq2, seqio[seqi].read_seq2);
#endif

#endif


#ifdef	QUAL_FILT_LV_OUT
			qual_filt_lv_1 = qual_filt_lv1[tid][1];
			qual_filt_lv_1_o = qual_filt_lv1[tid][0];
			qual_filt_lv_2 = qual_filt_lv2[tid][0];
			qual_filt_lv_2_o = qual_filt_lv2[tid][1];
#endif
			sam_cross = sam_pos2 - read_length1 - sam_pos1;

			xa_d1s[tid][xa_i] = '-';
			xa_d2s[tid][xa_i] = '+';

			pound_pos_1_f = pound_pos1_f_reverse;
			pound_pos_1_r = pound_pos1_r_reverse;
			pound_pos_2_f = pound_pos2_f_forward;
			pound_pos_2_r = pound_pos2_r_forward;

#ifdef	MAPPING_QUALITY
			if((mp_flag) && (v_cnt_i == 0))
			{
				mp_subs_1 = mp_subs1[tid][1];
				mp_subs_1_o = mp_subs1[tid][0];
				mp_subs_2 = mp_subs2[tid][0];
				mp_subs_2_o = mp_subs2[tid][1];
			}
#endif
		}

		d_n1 = 0;
		i_n1 = 0;
		s_offset1 = 0;
		s_offset2 = 0;
		s_r_o_l = ops_dm_l1[tid][v_cnt_i];
		s_r_o_r = ops_dm_r1[tid][v_cnt_i];

#ifdef	MAPPING_QUALITY
		if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] = 0;
#endif

		if((s_r_o_l == 0) && (s_r_o_r == 0))
		{
			strcpy(cigar_p1, cigar_m1[tid]);

			lv_re1 = ops_dm_ex1[tid][v_cnt_i];
#ifdef	MAPPING_QUALITY

			if((mp_flag) && (v_cnt_i == 0))
			{
				for (bit_char_i = 0, read_b_i = 32; bit_char_i < read_length1; bit_char_i++, read_b_i++)
					if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ops_vector_seq1[tid][v_cnt_i][read_b_i >> 5] >> ((31 - (read_b_i & 0X1f)) << 1)) & 0X3))
					{
						sub_t[tid] += mp_subs_1[bit_char_i];
					}
			}

#endif
		}
		else     //indel
		{
			if((cir_n == cir_fix_n) && (local_ksw))
			{
#ifdef	KSW_ALN_PAIR

#ifdef CHAR_CP
				for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < ops_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < ops_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
					read_char[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif

				for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; read_b_i < ops_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i > -1
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				ksw_global(ops_dm_kl1[tid][v_cnt_i], read_char[tid], ops_dm_kl1[tid][v_cnt_i], ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar1, &cigar1);

#ifdef	CHAR_CP
				for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < ops_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
					read_char2[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < ops_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
					read_char2[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif

				for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; read_b_i < ops_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length + 64
					ali_ref_seq2[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				ksw_global(ops_dm_kr1[tid][v_cnt_i], read_char2[tid], ops_dm_kr1[tid][v_cnt_i], ali_ref_seq2[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar2, &cigar2);

				m_m_n = s_r_o_r - s_r_o_l - 1;

				nm_score = 0;
				snt = 0;
				if (n_cigar1)
				{
					op_score1  = cigar1[0]&0xf;
					len_score1 = cigar1[0]>>4;
				}
				else	op_score1 = 3;

				if (n_cigar2)
				{
					op_score2  = cigar2[0]&0xf;
					len_score2 = cigar2[0]>>4;
				}
				else	op_score2 = 3;

				if (s_r_o_l >= ops_dm_kl1[tid][v_cnt_i])
				{
					sn = sprintf(cigar_p1 + snt, "%dS", s_r_o_l + 1 - ops_dm_kl1[tid][v_cnt_i]);
					snt += sn;
				}

				if ((s_r_o_l != -1) && (s_r_o_r != read_length1))
				{
					if ((op_score1 == 0) && (op_score2 == 0))
					{
						k_start1 = 0;
						k_start2 = 1;
						k_middle = len_score1 + len_score2 + m_m_n;
					}
					else if (op_score1 == 0)
					{
						k_start1 = 0;
						k_start2 = 0;
						k_middle = len_score1 + m_m_n;
					}
					else if (op_score2 == 0)
					{
						k_start1 = -1;
						k_start2 = 1;
						k_middle = len_score2 + m_m_n;
					}
					else
					{
						k_start1 = -1;
						k_start2 = 0;
						k_middle = m_m_n;
					}

					x_score = y_score = 0;
					for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
					{
						op_score  = cigar1[bit_char_i]&0xf;
						len_score = cigar1[bit_char_i]>>4;

						sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
								{
#ifdef	MAPPING_QUALITY
									if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_1_o[read_length1 + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
									++nm_score;
								}
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
					}

					sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
					snt += sn;

					x_score = y_score = 0;
					for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
					{
						op_score  = cigar2[bit_char_i]&0xf;
						len_score = cigar2[bit_char_i]>>4;

						sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
								{
#ifdef	MAPPING_QUALITY
									if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_1[s_r_o_r + x_score + read_b_i];
#endif
									++nm_score;
								}
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score;
					}
				}
				else if (s_r_o_l == -1)
				{
					if (op_score2 == 0)
					{
						k_start2 = 1;
						k_middle = len_score2 + m_m_n;
					}
					else
					{
						k_start2 = 0;
						k_middle = m_m_n;
					}
					sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
					snt += sn;

					x_score = y_score = 0;
					for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
					{
						op_score  = cigar2[bit_char_i]&0xf;
						len_score = cigar2[bit_char_i]>>4;

						sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
								{
#ifdef	MAPPING_QUALITY
									if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_1[s_r_o_r + x_score + read_b_i];
#endif
									++nm_score;
								}
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score;
					}
				}
				else
				{
					if (op_score1 == 0)
					{
						k_start1 = 0;
						k_middle = len_score1 + m_m_n;
					}
					else
					{
						k_start1 = -1;
						k_middle = m_m_n;
					}
					x_score = y_score = 0;
					for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
					{
						op_score  = cigar1[bit_char_i]&0xf;
						len_score = cigar1[bit_char_i]>>4;

						sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
								{
#ifdef	MAPPING_QUALITY
									if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_1_o[read_length1 + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
									++nm_score;
								}
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
					}

					sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
					snt += sn;
				}

				if (read_length1 - s_r_o_r > ops_dm_kr1[tid][v_cnt_i])
				{
					sn = sprintf(cigar_p1 + snt, "%dS", read_length1 - s_r_o_r - ops_dm_kr1[tid][v_cnt_i]);
					snt += sn;
				}
				//sn = sprintf(cigar_p1 + snt, "\0");
				//snt += sn;

				if (n_cigar1)	free(cigar1);
				if (n_cigar2)	free(cigar2);

				lv_re1 = nm_score;
#endif

			}
			else
			{
#ifdef	OUTPUT_DEBUG
				if (pound_pos_1_f >= s_r_o_r)   //1
				{
					lv_up_left = 0;
					lv_up_right = s_r_o_l;
					lv_down_right = pound_pos_1_f;
					lv_down_left = s_r_o_r;
					m_n_f = 0;
					m_n_b = read_length1 - pound_pos_1_f;
					m_m_n = s_r_o_r - s_r_o_l - 1;
				}
				else if (pound_pos_1_r <= s_r_o_l + 1)     //5
				{
					lv_up_left = pound_pos_1_r;//
					lv_up_right = s_r_o_l;
					lv_down_right = read_length1;
					lv_down_left = s_r_o_r;
					m_n_f = pound_pos_1_r;
					m_n_b = 0;
					m_m_n = s_r_o_r - s_r_o_l - 1;
				}
				else if ((pound_pos_1_f <= s_r_o_l + 1) && (pound_pos_1_r >= s_r_o_r))     //2
				{
					lv_up_left = 0;
					lv_up_right = pound_pos_1_f - 1;
					lv_down_right = read_length1;
					lv_down_left = pound_pos_1_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = pound_pos_1_r - pound_pos_1_f;
				}
				else if ((pound_pos_1_f > s_r_o_l + 1) && (pound_pos_1_f < s_r_o_r))     //3
				{
					lv_up_left = 0;
					lv_up_right = s_r_o_l;
					lv_down_right = read_length1;
					lv_down_left = pound_pos_1_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = read_length1 - s_r_o_l - 1;
				}
				else     //4
				{
					lv_up_left = 0;
					lv_up_right = -1;
					lv_down_right = read_length1;
					lv_down_left = s_r_o_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = s_r_o_r;
				}

#ifdef	QUAL_FILT_LV_OUT

#ifdef	CHAR_CP
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = sam_seq1[bit_char_i];

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY
				if((mp_flag) && (v_cnt_i == 0))
					lv_re1f = computeEditDistanceWithCigar_s_mis_left_mp(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length1 - 1 - lv_up_right, &s_offset1, mp_subs_1_o + read_length1 - 1 - lv_up_right, &(sub_t[tid]));//, 0, op_dm_sl1[tid][v_cnt_i]
				else
					lv_re1f = computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length1 - 1 - lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]
#else
				lv_re1f = computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length1 - 1 - lv_up_right, &s_offset1);//, 0, op_dm_sl1[tid][v_cnt_i]
#endif

#ifdef	CHAR_CP
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = sam_seq1[bit_char_i];

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

#ifdef	MAPPING_QUALITY
				if((mp_flag) && (v_cnt_i == 0))
					lv_re1b = computeEditDistanceWithCigar_s_mis_mp(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + lv_down_left, mp_subs_1 + lv_down_left, &(sub_t[tid]));//, 0, op_dm_sr1[tid][v_cnt_i]
				else
					lv_re1b = computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]
#else
				lv_re1b = computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#else

#ifdef	CHAR_CP
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = sam_seq1[bit_char_i];

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

#ifdef	MAPPING_QUALITY
				if((mp_flag) && (v_cnt_i == 0))
					lv_re1f = computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid], mp_subs_1_o + read_length1 - 1 - lv_up_right, &(sub_t[tid]));//, 0, op_dm_sl1[tid][v_cnt_i]
				else
					lv_re1f = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#else
				lv_re1f = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k1, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#endif

#ifdef	CHAR_CP
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = sam_seq1[bit_char_i];

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

#ifdef	MAPPING_QUALITY
				if((mp_flag) && (v_cnt_i == 0))
					lv_re1b = computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid], mp_subs_1 + lv_down_left, &(sub_t[tid]));//, 0, op_dm_sr1[tid][v_cnt_i]
				else
					lv_re1b = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#else
				lv_re1b = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k1, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#endif

#endif

				lv_re1 = lv_re1f + lv_re1b;

				//deal with front and back lv cigar
				strncpy(str_o, cigarBuf1, f_cigarn);
				s_o = 0;
				f_c = 0;

				pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

				while (pch != NULL)
				{
					pchl = strlen(pch);
					f_cigar[f_cigarn - f_c - 2] = atoi(pch);
					s_o += (pchl + 1);
					f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

					f_c += 2;

					if (str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
					if (str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

					pch = strtok_r(NULL, "DMIS", &saveptr);
				}

				strncpy(b_cigar, cigarBuf2, f_cigarn);
				pch = strtok(cigarBuf2,"DMIS");

				if (pch != NULL)
					pchl = strlen(pch);

				snt = 0;
#ifdef	CIGAR_S_MODIFY
				if(lv_up_left)
				{
					sn = sprintf(cigar_p1 + snt, "%uS", lv_up_left);
					snt += sn;
				}
#else
				if (m_n_f)
				{
					if(f_c)
					{
						if (f_cigar[f_cigarn + 1 - f_c] == 'M')
						{
							f_cigar[f_cigarn - f_c] += m_n_f;
						}
						else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
						{
							f_cigar[f_cigarn - f_c] += m_n_f;
						}
						else
						{
							sn = sprintf(cigar_p1 + snt, "%uM", m_n_f);
							snt += sn;
						}
					}
					else	m_m_n += m_n_f;

				}
#endif
				if ((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))
				{
					if ((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
					{
						f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p1 + snt, "%s", b_cigar + pchl + 1);
						snt += sn;
					}
					else if (f_cigar[f_cigarn - 1] == 'M')
					{
						f_cigar[f_cigarn - 2] += m_m_n;
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p1 + snt, "%s",b_cigar);
						snt += sn;
					}
					else if (b_cigar[pchl] == 'M')
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}

						sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						snt += sn;
					}
					else
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
				}
				else if ((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
				{
					if (b_cigar[pchl] == 'M')
					{
						sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						snt += sn;
					}
					else
					{
						sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
				}
				else if ((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
				{
					if (f_cigar[f_cigarn - 1] == 'M')
					{
						f_cigar[f_cigarn - 2] += m_m_n;
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
					}
					else
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
						snt += sn;
					}
				}
				else
				{
					sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
					snt += sn;
				}
#ifdef	CIGAR_S_MODIFY
				if(lv_down_right < read_length1)
				{
					sn = sprintf(cigar_p1 + snt, "%uS", read_length1 - lv_down_right);
					snt += sn;
				}
#else
				if (m_n_b)
				{
					if (cigar_p1[snt - 1] == 'M')
					{
						for (bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
						{
							if ((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
							m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
						}
						sn = sprintf(cigar_p1 + bit_char_i + 1, "%uM", m_n_b);
						snt = bit_char_i + 1 + sn;
					}
					else if(cigar_p1[snt - 1] == 'S')
					{
						for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
						{
							if((cigar_p1[bit_char_i] > 64) && (cigar_p1[bit_char_i] < 91))	break;
							m_n_b += (cigar_p1[bit_char_i] - '0') * carry_ten[f_i];
						}
						sn = sprintf(cigar_p1 + bit_char_i + 1, "%uS", m_n_b);
						snt = bit_char_i + 1 + sn;
					}
					else
					{
						sn = sprintf(cigar_p1 + snt, "%uM", m_n_b);
						snt += sn;
					}
				}
#endif

				//sprintf(cigar_p1 + snt, "\0");
			}

#ifdef	CIGAR_LEN_ERR
			cigar_len = 0;
			s_o_tmp = 0;
			strncpy(cigar_tmp, cigar_p1, snt);
			cigar_tmp[snt] = '\0';
			pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

			while (pch_tmp != NULL)
			{
				pchl_tmp = strlen(pch_tmp);
				s_o_tmp += (pchl_tmp + 1);

				if(cigar_p1[s_o_tmp - 1] != 'D')
				{
					cigar_len_tmp = atoi(pch_tmp);
					cigar_len += cigar_len_tmp;
				}

				pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
			}

			if(read_length1 != cigar_len)
			{
				if(read_length1 < cigar_len)
				{
					cigar_len_re = cigar_len_tmp - (cigar_len - read_length1);
					if(cigar_len_re > 0)	sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
					else if(cigar_len_re == 0)	sprintf(cigar_p1 + snt - sn, "\0");
					else	strcpy(cigar_p1, cigar_m1[tid]);
				}
				else
				{
					cigar_len_re = cigar_len_tmp + (read_length1 - cigar_len);
					sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
				}
			}
#endif

		}

		d_n2 = 0;
		i_n2 = 0;
		s_r_o_l = ops_dm_l2[tid][v_cnt_i];
		s_r_o_r = ops_dm_r2[tid][v_cnt_i];

		if((s_r_o_l == 0) && (s_r_o_r == 0))
		{
			strcpy(cigar_p2, cigar_m2[tid]);
			lv_re2 = ops_dm_ex2[tid][v_cnt_i];
#ifdef	MAPPING_QUALITY

			if((mp_flag) && (v_cnt_i == 0))
			{
				for (bit_char_i = 0, read_b_i = 32; bit_char_i < read_length2; bit_char_i++, read_b_i++)
					if(((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ops_vector_seq2[tid][v_cnt_i][read_b_i >> 5] >> ((31 - (read_b_i & 0X1f)) << 1)) & 0X3))
					{
						sub_t[tid] += mp_subs_2[bit_char_i];
					}
			}

#endif

		}
		else     //indel
		{
			if((cir_n == cir_fix_n) && (local_ksw))
			{
#ifdef	KSW_ALN_PAIR

#ifdef CHAR_CP
				for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < ops_dm_kl2[tid][v_cnt_i]; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
				for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < ops_dm_kl2[tid][v_cnt_i]; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = charToDna5n[sam_seq2[bit_char_i]];
#endif
				for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; read_b_i < ops_dm_kl2[tid][v_cnt_i]; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				ksw_global(ops_dm_kl2[tid][v_cnt_i], read_char[tid], ops_dm_kl2[tid][v_cnt_i], ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar1, &cigar1);

#ifdef CHAR_CP
				for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < ops_dm_kr2[tid][v_cnt_i]; bit_char_i++, read_b_i++)
					read_char2[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < ops_dm_kr2[tid][v_cnt_i]; bit_char_i++, read_b_i++)
					read_char2[tid][read_b_i] = charToDna5n[sam_seq2[bit_char_i]];
#endif

				for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; read_b_i < ops_dm_kr2[tid][v_cnt_i]; bit_char_i++, read_b_i++)
					ali_ref_seq2[tid][read_b_i] = ((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				ksw_global(ops_dm_kr2[tid][v_cnt_i], read_char2[tid], ops_dm_kr2[tid][v_cnt_i], ali_ref_seq2[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar2, &cigar2);

				m_m_n = s_r_o_r - s_r_o_l - 1;

				nm_score = 0;
				snt = 0;
				if (n_cigar1)
				{
					op_score1  = cigar1[0]&0xf;
					len_score1 = cigar1[0]>>4;
				}
				else	op_score1 = 3;

				if (n_cigar2)
				{
					op_score2  = cigar2[0]&0xf;
					len_score2 = cigar2[0]>>4;
				}
				else	op_score2 = 3;


				if (s_r_o_l >= ops_dm_kl2[tid][v_cnt_i])
				{
					sn = sprintf(cigar_p2 + snt, "%dS", s_r_o_l + 1 - ops_dm_kl2[tid][v_cnt_i]);
					snt += sn;
				}

				if ((s_r_o_l != -1) && (s_r_o_r != read_length2))
				{
					if ((op_score1 == 0) && (op_score2 == 0))
					{
						k_start1 = 0;
						k_start2 = 1;
						k_middle = len_score1 + len_score2 + m_m_n;
					}
					else if (op_score1 == 0)
					{
						k_start1 = 0;
						k_start2 = 0;
						k_middle = len_score1 + m_m_n;
					}
					else if (op_score2 == 0)
					{
						k_start1 = -1;
						k_start2 = 1;
						k_middle = len_score2 + m_m_n;
					}
					else
					{
						k_start1 = -1;
						k_start2 = 0;
						k_middle = m_m_n;
					}

					x_score = y_score = 0;
					for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
					{
						op_score  = cigar1[bit_char_i]&0xf;
						len_score = cigar1[bit_char_i]>>4;

						sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
								{
#ifdef	MAPPING_QUALITY
									if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_2_o[read_length2 + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
									++nm_score;
								}
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
					}

					sn = sprintf(cigar_p2 + snt, "%dM", k_middle);
					snt += sn;

					x_score = y_score = 0;
					for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
					{
						op_score  = cigar2[bit_char_i]&0xf;
						len_score = cigar2[bit_char_i]>>4;

						sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
								{
#ifdef	MAPPING_QUALITY
									if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_2[s_r_o_r + x_score + read_b_i];
#endif
									++nm_score;
								}
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score;
					}
				}
				else if (s_r_o_l == -1)
				{
					if (op_score2 == 0)
					{
						k_start2 = 1;
						k_middle = len_score2 + m_m_n;
					}
					else
					{
						k_start2 = 0;
						k_middle = m_m_n;
					}
					sn = sprintf(cigar_p2 + snt, "%dM", k_middle);
					snt += sn;

					x_score = y_score = 0;
					for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
					{
						op_score  = cigar2[bit_char_i]&0xf;
						len_score = cigar2[bit_char_i]>>4;

						sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
								{
#ifdef	MAPPING_QUALITY
									if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_2[s_r_o_r + x_score + read_b_i];
#endif
									++nm_score;
								}
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score;
					}
				}
				else
				{
					if (op_score1 == 0)
					{
						k_start1 = 0;
						k_middle = len_score1 + m_m_n;
					}
					else
					{
						k_start1 = -1;
						k_middle = m_m_n;
					}
					x_score = y_score = 0;
					for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
					{
						op_score  = cigar1[bit_char_i]&0xf;
						len_score = cigar1[bit_char_i]>>4;

						sn = sprintf(cigar_p2 + snt, "%d%c", len_score, ksw_cigars[op_score]);
						snt += sn;

						if (op_score == 0)
						{
							// match
							for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
								if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
								{
#ifdef	MAPPING_QUALITY
									if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_2_o[read_length2 + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
									++nm_score;
								}
							x_score += len_score;
							y_score += len_score;
						}
						else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
						else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
					}

					sn = sprintf(cigar_p2 + snt, "%dM", k_middle);
					snt += sn;
				}

				if (read_length2 - s_r_o_r > ops_dm_kr2[tid][v_cnt_i])
				{
					sn = sprintf(cigar_p2 + snt, "%dS", read_length2 - s_r_o_r - ops_dm_kr2[tid][v_cnt_i]);
					snt += sn;
				}
				//sn = sprintf(cigar_p2 + snt, "\0");
				//snt += sn;

				if (n_cigar1)	free(cigar1);
				if (n_cigar2)	free(cigar2);

				lv_re2 = nm_score;
#endif

			}
			else
			{
#ifdef	OUTPUT_DEBUG
				if (pound_pos_2_f >= s_r_o_r)   //1
				{
					lv_up_left = 0;
					lv_up_right = s_r_o_l;
					lv_down_right = pound_pos_2_f;
					lv_down_left = s_r_o_r;
					m_n_f = 0;
					m_n_b = read_length2 - pound_pos_2_f;
					m_m_n = s_r_o_r - s_r_o_l - 1;
				}
				else if (pound_pos_2_r <= s_r_o_l + 1)     //5
				{
					lv_up_left = pound_pos_2_r;//
					lv_up_right = s_r_o_l;
					lv_down_right = read_length2;
					lv_down_left = s_r_o_r;
					m_n_f = pound_pos_2_r;
					m_n_b = 0;
					m_m_n = s_r_o_r - s_r_o_l - 1;
				}
				else if ((pound_pos_2_f <= s_r_o_l + 1) && (pound_pos_2_r >= s_r_o_r))     //2
				{
					lv_up_left = 0;
					lv_up_right = pound_pos_2_f - 1;
					lv_down_right = read_length2;
					lv_down_left = pound_pos_2_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = pound_pos_2_r - pound_pos_2_f;
				}
				else if ((pound_pos_2_f > s_r_o_l + 1) && (pound_pos_2_f < s_r_o_r))     //3
				{
					lv_up_left = 0;
					lv_up_right = s_r_o_l;
					lv_down_right = read_length2;
					lv_down_left = pound_pos_2_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = read_length2 - s_r_o_l - 1;
				}
				else     //4
				{
					lv_up_left = 0;
					lv_up_right = -1;
					lv_down_right = read_length2;
					lv_down_left = s_r_o_r;
					m_n_f = 0;
					m_n_b = 0;
					m_m_n = s_r_o_r;
				}

#ifdef	QUAL_FILT_LV_OUT

#ifdef CHAR_CP
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = sam_seq2[bit_char_i];

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY
				if((mp_flag) && (v_cnt_i == 0))
					lv_re2f = computeEditDistanceWithCigar_s_mis_left_mp(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_2_o + read_length2 - 1 - lv_up_right, &s_offset2, mp_subs_2_o + read_length2 - 1 - lv_up_right, &(sub_t[tid]));//, 0, op_dm_sl1[tid][v_cnt_i]
				else
					lv_re2f = computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_2_o + read_length2 - 1 - lv_up_right, &s_offset2);//, 0, op_dm_sl1[tid][v_cnt_i]
#else
				lv_re2f = computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_2_o + read_length2 - 1 - lv_up_right, &s_offset2);//, 0, op_dm_sl1[tid][v_cnt_i]
#endif

#ifdef CHAR_CP
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = sam_seq2[bit_char_i];

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY
				if((mp_flag) && (v_cnt_i == 0))
					lv_re2b = computeEditDistanceWithCigar_s_mis_mp(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_2 + lv_down_left, mp_subs_2 + lv_down_left, &(sub_t[tid]));//, 0, op_dm_sr1[tid][v_cnt_i]
				else
					lv_re2b = computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_2 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]
#else
				lv_re2b = computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_2 + lv_down_left);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#else

#ifdef	CHAR_CP
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_up_right, read_b_i = 0; bit_char_i >= lv_up_left; bit_char_i--, read_b_i++)
					read_char[tid][read_b_i] = sam_seq2[bit_char_i];

				for (bit_char_i = 32 + lv_up_right, read_b_i = 0; bit_char_i > lv_up_left - 1; bit_char_i--, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

#ifdef	MAPPING_QUALITY
				if((mp_flag) && (v_cnt_i == 0))
					lv_re2f = computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid], mp_subs_2_o + read_length2 - 1 - lv_up_right, &(sub_t[tid]));//, 0, op_dm_sl1[tid][v_cnt_i]
				else
					lv_re2f = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#else
				lv_re2f = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + lv_up_right - lv_up_left, read_char[tid], lv_up_right + 1 - lv_up_left, lv_k2, cigarBuf1, f_cigarn, L[tid]);//, 0, op_dm_sl1[tid][v_cnt_i]
#endif

#ifdef	CHAR_CP
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = ((read_bit_2[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = ((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
				for (bit_char_i = lv_down_left, read_b_i = 0; bit_char_i < lv_down_right; bit_char_i++, read_b_i++)
					read_char[tid][read_b_i] = sam_seq2[bit_char_i];

				for (bit_char_i = 32 + lv_down_left, read_b_i = 0; bit_char_i < lv_down_right + 64; bit_char_i++, read_b_i++)
					ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq2[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];

#endif

#ifdef	MAPPING_QUALITY
				if((mp_flag) && (v_cnt_i == 0))
					lv_re2b = computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid], mp_subs_2 + lv_down_left, &(sub_t[tid]));//, 0, op_dm_sr1[tid][v_cnt_i]
				else
					lv_re2b = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#else
				lv_re2b = computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + lv_down_right - lv_down_left, read_char[tid], lv_down_right - lv_down_left, lv_k2, cigarBuf2, f_cigarn, L[tid]);//, 0, op_dm_sr1[tid][v_cnt_i]
#endif

#endif

#endif

				lv_re2 = lv_re2f + lv_re2b;

				//deal with front and back lv cigar
				strncpy(str_o, cigarBuf1, f_cigarn);
				s_o = 0;
				f_c = 0;

				pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

				while (pch != NULL)
				{
					pchl = strlen(pch);
					f_cigar[f_cigarn - f_c - 2] = atoi(pch);
					s_o += (pchl + 1);
					f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

					f_c += 2;

					if (str_o[s_o - 1] == 'D')	d_n2 += atoi(pch);
					if (str_o[s_o - 1] == 'I')	i_n2 += atoi(pch);

					pch = strtok_r(NULL, "DMIS", &saveptr);
				}

				strncpy(b_cigar, cigarBuf2, f_cigarn);
				pch = strtok(cigarBuf2,"DMIS");

				if (pch != NULL)
					pchl = strlen(pch);

				snt = 0;
#ifdef	CIGAR_S_MODIFY
				if(lv_up_left)
				{
					sn = sprintf(cigar_p2 + snt, "%uS", lv_up_left);
					snt += sn;
				}
#else
				if (m_n_f)
				{
					if(f_c)
					{
						if (f_cigar[f_cigarn + 1 - f_c] == 'M')
						{
							f_cigar[f_cigarn - f_c] += m_n_f;
						}
						else if(f_cigar[f_cigarn + 1 - f_c] == 'S')
						{
							f_cigar[f_cigarn - f_c] += m_n_f;
						}
						else
						{
							sn = sprintf(cigar_p2 + snt, "%uM", m_n_f);
							snt += sn;
						}
					}
					else	m_m_n += m_n_f;

				}
#endif
				if ((lv_up_right >= lv_up_left) && (lv_down_right > lv_down_left))
				{
					if ((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
					{
						f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p2 + snt, "%s", b_cigar + pchl + 1);
						snt += sn;
					}
					else if (f_cigar[f_cigarn - 1] == 'M')
					{
						f_cigar[f_cigarn - 2] += m_m_n;
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p2 + snt, "%s",b_cigar);
						snt += sn;
					}
					else if (b_cigar[pchl] == 'M')
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}

						sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						snt += sn;
					}
					else
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
				}
				else if ((lv_up_right < lv_up_left) && (lv_down_right > lv_down_left))     //op_dm_l1[tid][v_cnt_i] == -1
				{
					if (b_cigar[pchl] == 'M')
					{
						sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						snt += sn;
					}
					else
					{
						sn = sprintf(cigar_p2 + snt, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
				}
				else if ((lv_down_right <= lv_down_left) && (lv_up_right >= lv_up_left))
				{
					if (f_cigar[f_cigarn - 1] == 'M')
					{
						f_cigar[f_cigarn - 2] += m_m_n;
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
					}
					else
					{
						for (f_i = 0; f_i < f_c; f_i += 2)
						{
							sn = sprintf(cigar_p2 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
							snt += sn;
						}
						sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
						snt += sn;
					}
				}
				else
				{
					sn = sprintf(cigar_p2 + snt, "%uM", m_m_n);
					snt += sn;
				}
#ifdef	CIGAR_S_MODIFY
				if(lv_down_right < read_length2)
				{
					sn = sprintf(cigar_p2 + snt, "%uS", read_length2 - lv_down_right);
					snt += sn;
				}
#else
				if (m_n_b)
				{
					if (cigar_p2[snt - 1] == 'M')
					{
						for (bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
						{
							if ((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
							m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
						}
						sn = sprintf(cigar_p2 + bit_char_i + 1, "%uM", m_n_b);
						snt = bit_char_i + 1 + sn;
					}
					else if(cigar_p2[snt - 1] == 'S')
					{
						for(bit_char_i = snt - 2, f_i = 0; bit_char_i > -1; bit_char_i--, f_i++)
						{
							if((cigar_p2[bit_char_i] > 64) && (cigar_p2[bit_char_i] < 91))	break;
							m_n_b += (cigar_p2[bit_char_i] - '0') * carry_ten[f_i];
						}
						sn = sprintf(cigar_p2 + bit_char_i + 1, "%uS", m_n_b);
						snt = bit_char_i + 1 + sn;
					}
					else
					{
						sn = sprintf(cigar_p2 + snt, "%uM", m_n_b);
						snt += sn;
					}
				}
#endif
				//sprintf(cigar_p2 + snt, "\0");
			}

#ifdef	CIGAR_LEN_ERR
			cigar_len = 0;
			s_o_tmp = 0;
			strncpy(cigar_tmp, cigar_p2, snt);
			cigar_tmp[snt] = '\0';
			pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

			while (pch_tmp != NULL)
			{
				pchl_tmp = strlen(pch_tmp);
				s_o_tmp += (pchl_tmp + 1);

				if(cigar_p2[s_o_tmp - 1] != 'D')
				{
					cigar_len_tmp = atoi(pch_tmp);
					cigar_len += cigar_len_tmp;
				}

				pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
			}

			if(read_length2 != cigar_len)
			{
				if(read_length2 < cigar_len)
				{
					cigar_len_re = cigar_len_tmp - (cigar_len - read_length2);
					if(cigar_len_re > 0)	sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
					else if(cigar_len_re == 0)	sprintf(cigar_p2 + snt - sn, "\0");
					else	strcpy(cigar_p2, cigar_m2[tid]);
				}
				else
				{
					cigar_len_re = cigar_len_tmp + (read_length2 - cigar_len);
					sprintf(cigar_p2 + snt - sn, "%u%c", cigar_len_re, cigar_p2[snt - 1]);
				}
			}
#endif
		}

		sam_pos1s[tid][xa_i] = (uint32_t )(sam_pos1 + i_n1 - d_n1 + s_offset1);
		sam_pos2s[tid][xa_i] = (uint32_t )(sam_pos2 + i_n2 - d_n2 + s_offset2);

		strcpy(cigar_p1s[tid][xa_i], cigar_p1);
		strcpy(cigar_p2s[tid][xa_i], cigar_p2);

		lv_re1s[tid][xa_i] = lv_re1;
		lv_re2s[tid][xa_i] = lv_re2;

		++xa_i;
	}

#ifdef	MAPPING_QUALITY
	if(mp_flag)
	{
		if(xa_i)
		{
			log_tmp = log(vs_cnt);

			sam_qual1 = (sam_qual1 - sub_t[tid] - log_tmp) * 3.434;
			sam_qual1 = (sam_qual1 > 2 ? sam_qual1 : 2);
#ifdef	QUAL_20
			sam_qual1 = (sam_qual1 < 20 ? 20 : sam_qual1);
#endif
			sam_qual1 = (sam_qual1 > 60 ? 60 : sam_qual1);
			sam_qual2 = sam_qual1;

#ifdef	MP_DEBUG
			printf("2: %d %d\n\n", sam_qual1, sam_qual2);
#endif
		}
		else
		{
			sam_qual1 = 60;
			sam_qual2 = 60;
		}

	}
#endif

	seqio[seqi].qualc1 = sam_qual1;
	seqio[seqi].qualc2 = sam_qual2;

#ifdef	PR_SINGLE

	//plug-in single positions output
	xa_i1 = 0;
	xa_i2 = 0;
	if(pr_o_f[tid] == 1)
	{
		for(rc_i = 0; rc_i < 2; rc_i++)
			for(rc_ii = 0; rc_ii < 2; rc_ii++)
				for(v_cnt_i = 0; v_cnt_i < seedpos_misn[rc_i][rc_ii][tid]; v_cnt_i++)
				{
					if(seedpos_mis[rc_i][rc_ii][tid][v_cnt_i] == min_mis[tid])
					{
						x = seedpos[rc_i][rc_ii][tid][v_cnt_i];
						low = 0;
						high = chr_file_n - 1;

						while ( low <= high )
						{
							mid = (low + high) >> 1;
							if(x < (chr_end_n[mid]))
							{
								high = mid - 1;
							}
							else if(x > (chr_end_n[mid]))
							{
								low = mid + 1;
							}
							else
							{
								chr_re =  mid;
								break;
							}
							chr_re = low;
						}

						if(rc_i == rc_ii)
						{
							if(xa_i1 < pr_single_outputn)
							{
								pr_chr_res1[tid][xa_i1] = chr_re;
								pr_sam_pos1[tid][xa_i1] = x - chr_end_n[chr_re - 1] + 1;
								pr_xa_d1[tid][xa_i1] = pr_dir[rc_ii];
								pr_lv_re1[tid][xa_i1] = min_mis[tid];

								xa_i1++;
							}
						}
						else
						{
							if(xa_i2 < pr_single_outputn)
							{
								pr_chr_res2[tid][xa_i2] = chr_re;
								pr_sam_pos2[tid][xa_i2] = x - chr_end_n[chr_re - 1] + 1;
								pr_xa_d2[tid][xa_i2] = pr_dir[rc_ii];
								pr_lv_re2[tid][xa_i2] = min_mis[tid];

								xa_i2++;
							}
						}
					}
				}
	}

#endif

	seqio[seqi].v_cnt = v_cnt_out;

	//seqio[seqi].xa_n = xa_i;
	seqio[seqi].xa_n_p1 = xa_i;
	seqio[seqi].xa_n_p2 = xa_i;

#ifdef	FIX_SV
	seqio[seqi].xa_n_x1 = 0;
	seqio[seqi].xa_n_x2 = 0;
#endif

	if(xa_i > 0)
	{
		memcpy(chr_res_buffer[seqi], chr_res[tid], xa_i << 2);
		seqio[seqi].chr_res = chr_res_buffer[seqi];

		memcpy(xa_d1s_buffer[seqi], xa_d1s[tid], xa_i);
		seqio[seqi].xa_d1s = xa_d1s_buffer[seqi];

		memcpy(sam_pos1s_buffer[seqi], sam_pos1s[tid], xa_i << 2);
		seqio[seqi].sam_pos1s = sam_pos1s_buffer[seqi];

		memcpy(lv_re1s_buffer[seqi], lv_re1s[tid], xa_i << 2);
		seqio[seqi].lv_re1s = lv_re1s_buffer[seqi];

		for(v_cnt_i = 0; v_cnt_i < xa_i; v_cnt_i++)
		{
			f_i = strlen(cigar_p1s[tid][v_cnt_i]) + 1;
			seqio[seqi].cigar_p1s[v_cnt_i] = (char* )malloc(f_i);

			memcpy(seqio[seqi].cigar_p1s[v_cnt_i], cigar_p1s[tid][v_cnt_i], f_i);

			f_i = strlen(cigar_p2s[tid][v_cnt_i]) + 1;
			seqio[seqi].cigar_p2s[v_cnt_i] = (char* )malloc(f_i);

			memcpy(seqio[seqi].cigar_p2s[v_cnt_i], cigar_p2s[tid][v_cnt_i], f_i);
		}
	}
	seqio[seqi].xa_n1 = xa_i1;

#ifdef	PR_SINGLE
	if(pr_o_f[tid] == 1)
	{
		if(xa_i1 > 0)
		{
			memcpy(chr_res1_buffer[seqi], pr_chr_res1[tid], xa_i1);
			seqio[seqi].chr_res1 = chr_res1_buffer[seqi];

			memcpy(xa_ds1_buffer[seqi], pr_xa_d1[tid], xa_i1);
			seqio[seqi].xa_ds1 = xa_ds1_buffer[seqi];

			memcpy(sam_poss1_buffer[seqi], pr_sam_pos1[tid], xa_i1 << 2);
			seqio[seqi].sam_poss1 = sam_poss1_buffer[seqi];

			memcpy(lv_res1_buffer[seqi], pr_lv_re1[tid], xa_i1);
			seqio[seqi].lv_res1 = lv_res1_buffer[seqi];
		}
	}

#endif

	if(xa_i > 0)
	{
		memcpy(xa_d2s_buffer[seqi], xa_d2s[tid], xa_i);
		seqio[seqi].xa_d2s = xa_d2s_buffer[seqi];

		memcpy(sam_pos2s_buffer[seqi], sam_pos2s[tid], xa_i << 2);
		seqio[seqi].sam_pos2s = sam_pos2s_buffer[seqi];

		memcpy(lv_re2s_buffer[seqi], lv_re2s[tid], xa_i << 2);
		seqio[seqi].lv_re2s = lv_re2s_buffer[seqi];
	}
	seqio[seqi].xa_n2 = xa_i2;

#ifdef	PR_SINGLE
	if(pr_o_f[tid] == 1)
	{
		if(xa_i2 > 0)
		{
			memcpy(chr_res2_buffer[seqi], pr_chr_res2[tid], xa_i2);
			seqio[seqi].chr_res2 = chr_res2_buffer[seqi];

			memcpy(xa_ds2_buffer[seqi], pr_xa_d2[tid], xa_i2);
			seqio[seqi].xa_ds2 = xa_ds2_buffer[seqi];

			memcpy(sam_poss2_buffer[seqi], pr_sam_pos2[tid], xa_i2 << 2);
			seqio[seqi].sam_poss2 = sam_poss2_buffer[seqi];

			memcpy(lv_res2_buffer[seqi], pr_lv_re2[tid], xa_i2);
			seqio[seqi].lv_res2 = lv_res2_buffer[seqi];
		}
	}
#endif
}

#ifdef	PAIR_RANDOM

void Adjust(uint32_t* ls, uint32_t s, uint8_t tid)
{
	uint32_t i, t;

	t = ((s + seedn[tid]) >> 1);
	while(t > 0)
	{
		if(b[tid][s] > b[tid][ls[t]])
		{
			i = s;
			s = ls[t];
			ls[t] = i;
		}
		t >>= 1;
	}
	ls[0] = s;
}

void CreateLoserTree(uint32_t* ls, uint8_t tid)
{
	int64_t i;
	b[tid][seedn[tid]] = MINKEY;

	for(i = 0; i < seedn[tid]; ++i)
	{
		ls[i] = seedn[tid];
	}

	for(i = seedn[tid] - 1; i >= 0; --i)
	{
		Adjust(ls, i, tid);
	}
}

#ifdef UNPIPATH_OFF_K20
void seed_repetitive_single64(uint64_t* read_bit, uint32_t* pos_n, uint64_t** seedpos, uint8_t tid, uint16_t read_length)
#else
void seed_repetitive_single(uint64_t* read_bit, uint32_t* pos_n, uint32_t** seedpos, uint8_t tid, uint16_t read_length)
#endif
{
	uint8_t re_d = 0;
	uint8_t b_t_n_r = 0;
	uint32_t k_i, q, posk_n;
	uint32_t ref_pos_n = 0;
	uint32_t uni_offset_s_l = 0;
	uint16_t read_off = 0;
	uint32_t posn = 0;
	uint32_t seed_hash = 0;
	uint32_t seed_kmer = 0;
	uint32_t pos_r_ir = 0;
	uint32_t r_ru = 0;
	uint32_t ran_p = 0;
	uint64_t kmer_bit = 0;
	int64_t seed_id_r = 0;
	int64_t seed_binary_r = 0;
	float r_r = 0;

#ifdef UNPIPATH_OFF_K20
	uint64_t kmer_pos_uni = 0;
#else
	uint32_t kmer_pos_uni = 0;
#endif

	seedn[tid] = 0;
	//off_start[tid] = 0;
	seed_l[tid] = pair_ran_intvp;
	for(read_off = 0; read_off <= read_length - k_t; read_off += seed_l[tid])
	{
		//if(read_off + k_t - 1 <= r_b_v)	continue;
		re_d = (read_off & 0X1f);

		b_t_n_r = 32 - re_d;

#ifdef HASH_KMER_READ_J

		if(re_d <= re_b)
		{
			kmer_bit = ((read_bit[read_off >> 5] & bit_tran_re[re_d]) >> ((re_b - re_d) << 1));
		}
		else
		{
			kmer_bit = (((read_bit[read_off >> 5] & bit_tran_re[re_d]) << ((re_d - re_b) << 1)) | (read_bit[(read_off >> 5) + 1] >> ((re_2bt - re_d) << 1)));
		}
#else
		tran_tmp_p = (read_bit[read_off >> 5] & bit_tran_re[re_d]);

		//or use this method to deal with: & bit_tran_re[b_t_n_r]
		kmer_bit = (((read_bit[(read_off >> 5) + 1] >> (b_t_n_r << 1)) & bit_tran_re[b_t_n_r]) | (tran_tmp_p << (re_d << 1)));

		kmer_bit >>= re_bt;
#endif

		seed_kmer = (kmer_bit & bit_tran[k_r]);

		seed_hash = (kmer_bit >> (k_r << 1));

		//find the kmer
#ifdef UNPIPATH_OFF_K20
		seed_binary_r = binsearch_offset64(seed_kmer, buffer_kmer_g, buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], buffer_hash_g[seed_hash]);
#else
		seed_binary_r = binsearch_offset(seed_kmer, buffer_kmer_g, buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], buffer_hash_g[seed_hash]);
#endif
		if(seed_binary_r == -1)
			continue;

		//binary search on unipath offset to get the unipathID
		kmer_pos_uni = buffer_off_g[seed_binary_r];//this kmer's offset on unipath
#ifdef UNPIPATH_OFF_K20
		seed_id_r = binsearch_interval_unipath64(kmer_pos_uni, buffer_seqf, result_seqf);
#else
		seed_id_r = binsearch_interval_unipath(kmer_pos_uni, buffer_seqf, result_seqf);
#endif

		ref_pos_n = buffer_pp[seed_id_r + 1] - buffer_pp[seed_id_r];

		uni_offset_s_l = kmer_pos_uni - buffer_seqf[seed_id_r];

		if(ref_pos_n > pos_n_max)
		{
			if(ref_pos_n <= RANDOM_RANGE)
			{
				for(pos_r_ir = 0; pos_r_ir < ref_pos_n; pos_r_ir++)
					seed_k_pos[tid][seedn[tid]][pos_r_ir] = buffer_p[buffer_pp[seed_id_r] + pos_r_ir] + uni_offset_s_l - read_off;
				seed_k_pos[tid][seedn[tid]][pos_r_ir] = MAXKEY;
				seedn[tid]++;
			}
			else
			{
#ifdef	PAIR_RANDOM_SEED
				ran_p = (uint32_t )rand();
				r_r = ((float)ref_pos_n / RANDOM_RANGE_MAX);

				if(r_r >= 2)
				{
					r_ru = (uint32_t )r_r;
					for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE; pos_r_ir++, ran_p++)
						seed_k_pos[tid][seedn[tid]][pos_r_ir] = buffer_p[buffer_pp[seed_id_r] + (pos_r_ir * r_ru) + (random_buffer[ran_p & 0Xfff] % r_ru)] + uni_offset_s_l - read_off;
				}
				else
				{
					for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE; pos_r_ir++)
					{
						seed_k_pos[tid][seedn[tid]][pos_r_ir] = buffer_p[buffer_pp[seed_id_r] + (uint32_t )(seed_r_dup[pos_r_ir] * r_r)] + uni_offset_s_l - read_off;
					}
				}

				seed_k_pos[tid][seedn[tid]][pos_r_ir] = MAXKEY;
				seedn[tid]++;
#else
				r_r = ((float)ref_pos_n / RANDOM_RANGE);
				for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE; pos_r_ir++)
					seed_k_pos[tid][seedn[tid]][pos_r_ir] = buffer_p[buffer_pp[seed_id_r] + (uint32_t )(pos_r_ir * r_r)] + uni_offset_s_l - read_off;
				seed_k_pos[tid][seedn[tid]][pos_r_ir] = MAXKEY;
				seedn[tid]++;
#endif
			}
		}
	}

	posk_n = 0;

	for(k_i = 0; k_i < seedn[tid]; ++k_i)
	{
		b[tid][k_i] = seed_k_pos[tid][k_i][kcol[tid][k_i]++];
	}

	CreateLoserTree(ls[tid], tid);

	while(b[tid][ls[tid][0]] != MAXKEY)
	{
		q = ls[tid][0];

		seedposk[tid][posk_n++] = b[tid][q];

		if(kcol[tid][q] < MAX_COL)
		{
			b[tid][q] = seed_k_pos[tid][q][kcol[tid][q]++];
			Adjust(ls[tid],q, tid);
		}
	}

	memset(kcol[tid], 0, seedn[tid] << 2);

	posn = 1;
	(*seedpos)[0] = seedposk[tid][0];
	for(k_i = 1; k_i < posk_n; k_i++)
	{
		if(seedposk[tid][k_i] != seedposk[tid][k_i - 1])
			(*seedpos)[posn++] = seedposk[tid][k_i];
	}

	(*pos_n) = posn;
}

void seed_repetitive(uint8_t tid, uint16_t read_length1, uint16_t read_length2, uint16_t f_cigarn, uint16_t ref_copy_num1, uint16_t ref_copy_num2, uint32_t ref_copy_num_chars1, uint32_t ref_copy_num_chars2, cnt_re* cnt, uint32_t seqi, uint16_t lv_k1, uint16_t lv_k2, int16_t pound_pos1_f_forward, int16_t pound_pos1_f_reverse, int16_t pound_pos1_r_forward, int16_t pound_pos1_r_reverse, int16_t pound_pos2_f_forward, int16_t pound_pos2_f_reverse, int16_t pound_pos2_r_forward, int16_t pound_pos2_r_reverse)
{
	uint8_t rc_i = 0;
	uint8_t rc_ii = 0;
	uint8_t re_d = 0;
	uint8_t b_t_n_r = 0;
	uint8_t q_n1 = 0;
	uint8_t c_m_f = 0;
	uint8_t max_mismatch = 0;
	uint16_t rst_i = 0;
	uint16_t s_m_t = 0;
	uint16_t read_b_i = 0;
	uint16_t ref_copy_num = 0;
	uint16_t end_dis[2];
	uint16_t read_length = 0;

	int bit_char_i = 0;
	int dm1 = 0;
	int dmt1 = 0;
	int lv_dmt1 = 0;
	int dm12 = 0;
	int dm_l1 = 0;
	int dm_r1 = 0;
	int ld1 = 0;
	int rd1 = 0;
	int ld2 = 0;
	int rd2 = 0;
	int s_r_o_l1 = 0;
	int	s_r_o_r1 = 0;
	int cmp_re = 0;
	int q_rear_i = 0;
	int q_rear1 = 0;
	int cache_dml1[MAX_Q_NUM];
	int cache_dmr1[MAX_Q_NUM];
	int cache_dis1[MAX_Q_NUM];

	uint32_t ref_copy_num_chars = 0;
	uint32_t ref_copy_num_chars_1 = 0;
	uint32_t ref_copy_num_chars_2 = 0;
	uint32_t pos1_t = 0;
	uint32_t pos2_t = 0;
	uint32_t seed_nt = 0;
	uint32_t seed_nt1 = 0;
	uint32_t seed_nt2 = 0;
	uint32_t inner_n = 0;
	uint32_t seed_sn = 0;
	uint32_t seednn[2];
	uint32_t mis_c_n = 0;
	uint32_t mis_c_n_filt = 0;


	int64_t b_r_s = 0;
	int64_t b_r_s_i = 0;
	uint64_t tran_tmp_p = 0;
	uint64_t tran_tmp = 0;
	uint64_t xor_tmp = 0;
	uint64_t ref_tmp_ori = 0;
	uint64_t ref_tmp_ori2 = 0;
	uint64_t cache_end1[MAX_Q_NUM][MAX_REF_SEQ_C];
	uint64_t low_mask = 0;
	uint16_t* sub_mask = NULL;
	uint64_t* ex_d_mask = NULL;

	dm_op[tid] = MAX_OP_SCORE;
	dm_ops[tid] = MAX_OP_SCORE;
	cnt->v_cnt = 0;
	cnt->vs_cnt = 0;

#ifdef	QUAL_FILT_REP
	uint64_t* qual_filt_1 = NULL;
#endif

#ifdef	PR_SINGLE
	min_mis[tid] = 0Xff;
	seedpos_misn[0][0][tid] = 0;
	seedpos_misn[0][1][tid] = 0;
	seedpos_misn[1][0][tid] = 0;
	seedpos_misn[1][1][tid] = 0;
#endif

#ifdef UNPIPATH_OFF_K20
	uint64_t pos_l = 0;
	uint64_t pos_t1 = 0;
	uint64_t pos_t2 = 0;
	uint64_t posi = 0;
	uint64_t up_pos = 0;
	uint64_t down_pos = 0;
	uint64_t* seed_pt = NULL;
#else
	uint32_t pos_l = 0;
	uint32_t pos_t1 = 0;
	uint32_t pos_t2 = 0;
	uint32_t posi = 0;
	uint32_t up_pos = 0;
	uint32_t down_pos = 0;
	uint32_t* seed_pt = NULL;
#endif

	end_dis[0] = end_dis1[tid];
	end_dis[1] = end_dis2[tid];

	for(rc_i = 0; rc_i < 2; rc_i++)
	{
		//check whether can make pair
		if((pos_ren[rc_i][0][tid] == 0) || (pos_ren[rc_i][1][tid] == 0))	continue;
#ifdef UNPIPATH_OFF_K20
		seed_repetitive_single64(read_bit1[tid][rc_i], &(seed_posn[rc_i][tid]), &(seedpos[rc_i][rc_i][tid]), tid, read_length1);
		seed_repetitive_single64(read_bit2[tid][1 - rc_i], &(seed_posn[1 - rc_i][tid]), &(seedpos[rc_i][1 - rc_i][tid]), tid, read_length2);

#else
		seed_repetitive_single(read_bit1[tid][rc_i], &(seed_posn[rc_i][tid]), &(seedpos[rc_i][rc_i][tid]), tid, read_length1);
		seed_repetitive_single(read_bit2[tid][1 - rc_i], &(seed_posn[1 - rc_i][tid]), &(seedpos[rc_i][1 - rc_i][tid]), tid, read_length2);

#endif

		seed_nt1 = seed_posn[0][tid];
		seed_nt2 = seed_posn[1][tid];

		if((seed_nt1 == 0) || (seed_nt2 == 0))	continue;

#ifdef	PR_SINGLE
		//do exact match to get mismatch number on every position in seedpos
		for(rc_ii = 0; rc_ii < 2; rc_ii++)
		{
			if(rc_i == rc_ii)
			{
				ref_copy_num_chars = ref_copy_num_chars1;
				ref_copy_num = ref_copy_num1;
				low_mask = low_mask1[tid];
				sub_mask = sub_mask1[tid];
				ex_d_mask = ex_d_mask1[tid];
				max_mismatch = max_mismatch1[tid];
#ifdef	QUAL_FILT_REP
				if(rc_i == 0)	qual_filt_1 = qual_filt1[tid][0];
				else	qual_filt_1 = qual_filt1[tid][1];
#endif
			}
			else
			{
				ref_copy_num_chars = ref_copy_num_chars2;
				ref_copy_num = ref_copy_num2;
				low_mask = low_mask2[tid];
				sub_mask = sub_mask2[tid];
				ex_d_mask = ex_d_mask2[tid];
				max_mismatch = max_mismatch2[tid];
#ifdef	QUAL_FILT_REP
				if(rc_i == 0)	qual_filt_1 = qual_filt2[tid][1];
				else	qual_filt_1 = qual_filt2[tid][0];
#endif
			}

			read_bit_1[tid] = read_bit_pr[tid][(rc_i << 1) + rc_ii];

			for(pos1_t = 0; pos1_t < seed_posn[rc_ii][tid]; pos1_t++)
			{
				pos_l = seedpos[rc_i][rc_ii][tid][pos1_t] - max_extension_length - 1;
				re_d = pos_l & 0X1f;
				b_t_n_r = 32 - re_d;

				if(re_d != 0)
				{
					tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);
					memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars);

					for(rst_i = 0; rst_i < ref_copy_num - 1; rst_i++)
					{
						tran_tmp = (ref_seq_tmp1[tid][rst_i] & bit_tran_re[re_d]);

						ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
						ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
						tran_tmp_p = tran_tmp;
					}
					ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
					ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

					//clear the lowest n bit
					ref_seq_tmp1[tid][rst_i] &= low_mask;

				}
				else
				{
					memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars);

					//clear the lowest n bit
					ref_seq_tmp1[tid][ref_copy_num - 1] &= low_mask;
				}

				//trim the beginning and end of the current ref seq based on current minimum edit distance dm_t
				s_m_t = sub_mask[0];
				ref_seq_tmp1[tid][0] &= bit_tran[0];
				ref_seq_tmp1[tid][s_m_t] &= ex_d_mask[0];

				mis_c_n = 0;
#ifdef	QUAL_FILT_REP
				mis_c_n_filt = 0;
#endif
				//ref_tmp_ori = ref_seq_tmp1[ref_copy_num - 2];
				ref_seq_tmp1[tid][ref_copy_num - 2] &= low_mask;

				for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num - 1; rst_i++, read_b_i++)
				{
					xor_tmp = ref_seq_tmp1[tid][rst_i] ^ read_bit_1[tid][read_b_i];
#ifdef	QUAL_FILT_REP
					mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
					xor_tmp &= qual_filt_1[read_b_i];
					mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
					if(mis_c_n_filt > max_mismatch)	break;
#else
					mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
					mis_c_n_filt = mis_c_n;
					if(mis_c_n > max_mismatch)	break;
#endif
				}

				if(mis_c_n_filt > max_mismatch)
				{
					seedpos_mis[rc_i][rc_ii][tid][pos1_t] = MAX_SPM;
				}
				else
				{
					seedpos_mis[rc_i][rc_ii][tid][pos1_t] = mis_c_n;
					if(mis_c_n < min_mis[tid])	min_mis[tid] = mis_c_n;
				}
			}
			seedpos_misn[rc_i][rc_ii][tid] = seed_posn[rc_ii][tid];
		}
#endif

		memset(seed_posf[0][tid], 0, seed_nt1);
		memset(seed_posf[1][tid], 0, seed_nt2);

		for(pos1_t = 0, pos2_t = 0; (pos1_t < seed_nt1) && (pos2_t < seed_nt2); pos1_t++)
		{
			up_pos = seedpos[rc_i][0][tid][pos1_t] + end_dis[rc_i] - devi;
			down_pos = seedpos[rc_i][0][tid][pos1_t] + end_dis[rc_i] + devi;
			seed_pt = seedpos[rc_i][1][tid];

			inner_n = 0;
			for( ; (seed_pt[pos2_t] < down_pos) && (pos2_t < seed_nt2); pos2_t++)
				if(seed_pt[pos2_t] > up_pos)
				{
					seed_posf[1][tid][pos2_t] = 1;
					inner_n++;
				}
			if(inner_n > 0)	seed_posf[0][tid][pos1_t] = 1;
		}

		for(pos1_t = 0, pos2_t = 0; (pos1_t < seed_nt1) && (pos2_t < seed_nt2); pos2_t++)
		{
			up_pos = seedpos[rc_i][1][tid][pos2_t] - end_dis[rc_i] - devi;
			down_pos = seedpos[rc_i][1][tid][pos2_t] - end_dis[rc_i] + devi;
			seed_pt = seedpos[rc_i][0][tid];

			inner_n = 0;
			for( ; (seed_pt[pos1_t] < down_pos) && (pos1_t < seed_nt2); pos1_t++)
				if(seed_pt[pos1_t] > up_pos)
				{
					seed_posf[0][tid][pos1_t] = 1;
					inner_n++;
				}
			if(inner_n > 0)	seed_posf[1][tid][pos2_t] = 1;
		}

		for(rc_ii = 0; rc_ii < 2; rc_ii++)
		{
			if(rc_i == rc_ii)
			{
				ref_copy_num_chars = ref_copy_num_chars1;
				ref_copy_num = ref_copy_num1;
				low_mask = low_mask1[tid];
				sub_mask = sub_mask1[tid];
				ex_d_mask = ex_d_mask1[tid];
				lv_dmt1 = lv_k1;
				read_length = read_length1;
				max_mismatch = max_mismatch1[tid];
#ifdef	QUAL_FILT_REP
				if(rc_i == 0)	qual_filt_1 = qual_filt1[tid][0];
				else	qual_filt_1 = qual_filt1[tid][1];
#endif
			}
			else
			{
				ref_copy_num_chars = ref_copy_num_chars2;
				ref_copy_num = ref_copy_num2;
				low_mask = low_mask2[tid];
				sub_mask = sub_mask2[tid];
				ex_d_mask = ex_d_mask2[tid];
				lv_dmt1 = lv_k2;
				read_length = read_length2;
				max_mismatch = max_mismatch2[tid];
#ifdef	QUAL_FILT_REP
				if(rc_i == 0)	qual_filt_1 = qual_filt2[tid][1];
				else	qual_filt_1 = qual_filt2[tid][0];
#endif
			}

			read_bit_1[tid] = read_bit_pr[tid][(rc_i << 1) + rc_ii];
			read_val_1[tid] = read_val_pr[tid][(rc_i << 1) + rc_ii];

			q_rear1 = 0;
			q_n1 = 0;
			dmt1 = ali_exl;

			s_r_o_l1 = 0;
			s_r_o_r1 = 0;
			seed_sn = 0;
			for(pos1_t = 0; pos1_t < seed_posn[rc_ii][tid]; pos1_t++)
			{
				if(seed_posf[rc_ii][tid][pos1_t] == 1)
				{

#ifdef	PR_SINGLE
					if(seedpos_mis[rc_i][rc_ii][tid][pos1_t] == MAX_SPM)
					{
#endif
						//lv_f1 = 0;
						pos_l = seedpos[rc_i][rc_ii][tid][pos1_t] - max_extension_length - 1;
						re_d = pos_l & 0X1f;
						b_t_n_r = 32 - re_d;

						if(re_d != 0)
						{
							tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);
							memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars);

							for(rst_i = 0; rst_i < ref_copy_num - 1; rst_i++)
							{
								tran_tmp = (ref_seq_tmp1[tid][rst_i] & bit_tran_re[re_d]);

								ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
								ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
								tran_tmp_p = tran_tmp;
							}
							ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
							ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

							//clear the lowest n bit
							ref_seq_tmp1[tid][rst_i] &= low_mask;

						}
						else
						{
							memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars);

							//clear the lowest n bit
							ref_seq_tmp1[tid][ref_copy_num - 1] &= low_mask;
						}

						//trim the beginning and end of the current ref seq based on current minimum edit distance dm_t
						s_m_t = sub_mask[dmt1];
						ref_seq_tmp1[tid][0] &= bit_tran[dmt1];
						ref_seq_tmp1[tid][s_m_t] &= ex_d_mask[dmt1];

						//traverse and check whether there is an existing seq that is as same as current new ref seq
						c_m_f = 0;

						for(q_rear_i = q_rear1 - 1; q_rear_i >= 0; q_rear_i--)
						{
							ref_tmp_ori = cache_end1[q_rear_i][0];
							cache_end1[q_rear_i][0] &= bit_tran[dmt1];

							ref_tmp_ori2 = cache_end1[q_rear_i][s_m_t];
							cache_end1[q_rear_i][s_m_t] &= ex_d_mask[dmt1];

							cmp_re = memcmp(cache_end1[q_rear_i], ref_seq_tmp1[tid], (s_m_t + 1) << 3);

							cache_end1[q_rear_i][0] = ref_tmp_ori;
							cache_end1[q_rear_i][s_m_t] = ref_tmp_ori2;

							if(cmp_re == 0)
							{
								//deal with finding an alignment

								dm1 = cache_dis1[q_rear_i];
								ld1 = cache_dml1[q_rear_i];
								rd1 = cache_dmr1[q_rear_i];

								c_m_f = 1;
								break;

							}
						}

						if((q_n1 > MAX_Q_NUM) && (q_rear_i < 0))
						{
							for(q_rear_i = MAX_Q_NUM - 1; q_rear_i >= q_rear1; q_rear_i--)
							{
								ref_tmp_ori = cache_end1[q_rear_i][0];
								cache_end1[q_rear_i][0] &= bit_tran[dmt1];

								ref_tmp_ori2 = cache_end1[q_rear_i][s_m_t];
								cache_end1[q_rear_i][s_m_t] &= ex_d_mask[dmt1];

								cmp_re = memcmp(cache_end1[q_rear_i], ref_seq_tmp1[tid], (s_m_t + 1) << 3);

								cache_end1[q_rear_i][0] = ref_tmp_ori;
								cache_end1[q_rear_i][s_m_t] = ref_tmp_ori2;

								if(cmp_re == 0)
								{
									//deal with finding an alignment

									dm1 = cache_dis1[q_rear_i];
									ld1 = cache_dml1[q_rear_i];
									rd1 = cache_dmr1[q_rear_i];

									c_m_f = 1;
									break;

								}
							}
						}

						//do not find the seq in cache, exact match or lv and add into cache
						if(c_m_f == 0)
						{
#ifdef	PR_SINGLE
							mis_c_n = max_mismatch + 1;
#else
							//exact match
							mis_c_n = 0;
#ifdef	QUAL_FILT_REP
							mis_c_n_filt = 0;
#endif
							ref_tmp_ori = ref_seq_tmp1[tid][ref_copy_num - 2];
							ref_seq_tmp1[tid][ref_copy_num - 2] &= low_mask;

							for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num - 1; rst_i++, read_b_i++)
							{
								xor_tmp = ref_seq_tmp1[tid][rst_i] ^ read_bit_1[tid][read_b_i];
#ifdef	QUAL_FILT_REP
								mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								xor_tmp &= qual_filt_1[read_b_i];
								mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
								if(mis_c_n_filt > max_mismatch)	break;
#else
							mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
							mis_c_n_filt = mis_c_n;
							if(mis_c_n > max_mismatch)	break;
#endif
							}

							ref_seq_tmp1[tid][ref_copy_num - 2] = ref_tmp_ori;
#endif
							//lv
							if(mis_c_n_filt > max_mismatch)
							{
#ifdef SPLIT_LV
								//split lv
								for(bit_char_i = s_r_o_l1, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
									read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - bit_char_i) << 1)) & 0X3);

								for(bit_char_i = 32 + s_r_o_l1, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
									ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - bit_char_i) << 1)) & 0X3);
#ifdef LV_CIGAR
								dm_l1 = computeEditDistanceWithCigar(ali_ref_seq[tid], 33 + s_r_o_l1, read_char[tid], s_r_o_l1 + 1, dmt1, cigarBuf, f_cigarn, 0, 0, L[tid]);
#else
								dm_l1 = computeEditDistance(ali_ref_seq[tid], 33 + s_r_o_l1, read_char[tid], s_r_o_l1 + 1, lv_dmt1, L[tid]);
#endif
								for(bit_char_i = s_r_o_r1, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
									read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - bit_char_i) << 1)) & 0X3);

								for(bit_char_i = 32 + s_r_o_r1, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
									ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - bit_char_i) << 1)) & 0X3);
#ifdef	LV_CIGAR
								dm_r1 = computeEditDistanceWithCigar(ali_ref_seq[tid], 32 + read_length - s_r_o_r1, read_char[tid], read_length - s_r_o_r1, dmt1, cigarBuf, f_cigarn, 0, 0, L[tid]);
#else
								dm_r1 = computeEditDistance(ali_ref_seq[tid], 32 + read_length - s_r_o_r1, read_char[tid], read_length - s_r_o_r1, lv_dmt1, L[tid]);
#endif
								dm1 = dm_l1 + dm_r1;

								ld1 = s_r_o_l1;
								rd1 = s_r_o_r1;
#endif
							}
							else
							{
								dm1 = mis_c_n;
								ld1 = 0;
								rd1 = 0;
								dm_l1 = 0;
								dm_r1 = 0;
							}

							//these two values could be different

							if((dm_l1 != -1) && (dm_r1 != -1))
							{
								if(dm1 < dmt1)	dmt1 = dm1 + 1;

								if(dm1 < lv_dmt1)	lv_dmt1 = dm1 + 1;

								if(dm1 < max_mismatch - 1)	dmt1 = 0;
							}
							else	dm1 = MAX_EDIT_SCORE;

							//add the ref sequence at the end of queue
							memcpy(cache_end1[q_rear1], ref_seq_tmp1[tid], ref_copy_num_chars);
							cache_dis1[q_rear1] = dm1;
							cache_dml1[q_rear1] = ld1;
							cache_dmr1[q_rear1] = rd1;

							q_rear1 = ((q_rear1 + 1) & 0X1f);
							++q_n1;

							//add edit distance
						}

						if(dm1 != MAX_EDIT_SCORE)
						{
							seed_single_pos[rc_ii][tid][seed_sn] = seedpos[rc_i][rc_ii][tid][pos1_t];
							seed_single_ld[rc_ii][tid][seed_sn] = ld1;
							seed_single_rd[rc_ii][tid][seed_sn] = rd1;
							seed_single_dm[rc_ii][tid][seed_sn] = dm1;
							memcpy(seed_single_refs[rc_ii][tid][seed_sn], ref_seq_tmp1[tid], ref_copy_num_chars);
							seed_sn++;
						}
#ifdef	PR_SINGLE
					}
					else
					{
						seed_single_pos[rc_ii][tid][seed_sn] = seedpos[rc_i][rc_ii][tid][pos1_t];
						seed_single_ld[rc_ii][tid][seed_sn] = 0;
						seed_single_rd[rc_ii][tid][seed_sn] = 0;
						seed_single_dm[rc_ii][tid][seed_sn] = seedpos_mis[rc_i][rc_ii][tid][pos1_t];
						seed_sn++;
					}
#endif
				}
			}
			seednn[rc_ii] = seed_sn;
		}

		if((seednn[0] == 0) || (seednn[1] == 0)) continue;

		if(rc_i == 0)
		{
			ref_copy_num_chars_1 = ref_copy_num_chars1;
			ref_copy_num_chars_2 = ref_copy_num_chars2;
		}
		else
		{
			ref_copy_num_chars_1 = ref_copy_num_chars2;
			ref_copy_num_chars_2 = ref_copy_num_chars1;
		}


		if(seednn[0] <= seednn[1])
		{
			for(pos1_t = 0; pos1_t < seednn[0]; pos1_t++)
			{
				pos_t1 = seed_single_pos[0][tid][pos1_t];
				posi = pos_t1 + end_dis[rc_i];
#ifdef UNPIPATH_OFF_K20
				b_r_s = binsearch_seed_pos_ss64(posi, seed_single_pos[1][tid], seednn[1]);
#else
				b_r_s = binsearch_seed_pos_ss(posi, seed_single_pos[1][tid], seednn[1]);
#endif
				if(b_r_s != -1)
				{
					dm12 = seed_single_dm[0][tid][pos1_t] + seed_single_dm[1][tid][b_r_s];
					pos_t2 = seed_single_pos[1][tid][b_r_s];
					ld1 = seed_single_ld[0][tid][pos1_t];
					ld2 = seed_single_ld[1][tid][b_r_s];
					rd1 = seed_single_rd[0][tid][pos1_t];
					rd2 = seed_single_rd[1][tid][b_r_s];

					memcpy(ref_seq_tmp1[tid], seed_single_refs[0][tid][pos1_t], ref_copy_num_chars_1);
					memcpy(ref_seq_tmp2[tid], seed_single_refs[1][tid][b_r_s], ref_copy_num_chars_2);
#ifdef UNPIPATH_OFF_K20
					cnt = pair_op_add64(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][pos1_t], seed_single_dm[1][tid][b_r_s], tid, ref_copy_num_chars1,  ref_copy_num_chars2, cnt);
#else
					cnt = pair_op_add(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][pos1_t], seed_single_dm[1][tid][b_r_s], tid, ref_copy_num_chars1,  ref_copy_num_chars2, cnt);
#endif

#ifdef	PAIR_TRAVERSE
					for(b_r_s_i = b_r_s + 1; (seed_single_pos[1][tid][b_r_s_i] < posi + devi) && (b_r_s_i < seednn[1]); b_r_s_i++)
					{
						dm12 = seed_single_dm[0][tid][pos1_t] + seed_single_dm[1][tid][b_r_s_i];
						pos_t2 = seed_single_pos[1][tid][b_r_s_i];
						ld2 = seed_single_ld[1][tid][b_r_s_i];
						rd2 = seed_single_rd[1][tid][b_r_s_i];
						memcpy(ref_seq_tmp2[tid], seed_single_refs[1][tid][b_r_s_i], ref_copy_num_chars_2);
#ifdef UNPIPATH_OFF_K20
						cnt = pair_op_add64(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][pos1_t], seed_single_dm[1][tid][b_r_s_i], tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#else
						cnt = pair_op_add(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][pos1_t], seed_single_dm[1][tid][b_r_s_i], tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#endif
					}

					for(b_r_s_i = b_r_s - 1; (seed_single_pos[1][tid][b_r_s_i] > posi - devi) && (b_r_s_i >= 0); b_r_s_i--)
					{
						dm12 = seed_single_dm[0][tid][pos1_t] + seed_single_dm[1][tid][b_r_s_i];
						pos_t2 = seed_single_pos[1][tid][b_r_s_i];
						ld2 = seed_single_ld[1][tid][b_r_s_i];
						rd2 = seed_single_rd[1][tid][b_r_s_i];
						memcpy(ref_seq_tmp2[tid], seed_single_refs[1][tid][b_r_s_i], ref_copy_num_chars_2);
#ifdef UNPIPATH_OFF_K20
						cnt = pair_op_add64(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][pos1_t], seed_single_dm[1][tid][b_r_s_i], tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#else
						cnt = pair_op_add(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][pos1_t], seed_single_dm[1][tid][b_r_s_i], tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#endif
					}
#endif
				}
			}
		}
		else
		{
			for(pos2_t = 0; pos2_t < seednn[1]; pos2_t++)
			{
				pos_t2 = seed_single_pos[1][tid][pos2_t];
				posi = pos_t2 - end_dis[rc_i];
#ifdef UNPIPATH_OFF_K20
				b_r_s = binsearch_seed_pos_ss64(posi, seed_single_pos[0][tid], seednn[0]);
#else
				b_r_s = binsearch_seed_pos_ss(posi, seed_single_pos[0][tid], seednn[0]);
#endif
				if(b_r_s != -1)
				{
					dm12 = seed_single_dm[1][tid][pos2_t] + seed_single_dm[0][tid][b_r_s];
					pos_t1 = seed_single_pos[0][tid][b_r_s];
					ld1 = seed_single_ld[0][tid][b_r_s];
					ld2 = seed_single_ld[1][tid][pos2_t];
					rd1 = seed_single_rd[0][tid][b_r_s];
					rd2 = seed_single_rd[1][tid][pos2_t];

					memcpy(ref_seq_tmp1[tid], seed_single_refs[0][tid][b_r_s], ref_copy_num_chars_1);
					memcpy(ref_seq_tmp2[tid], seed_single_refs[1][tid][pos2_t], ref_copy_num_chars_2);
#ifdef UNPIPATH_OFF_K20
					cnt = pair_op_add64(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][b_r_s], seed_single_dm[1][tid][pos2_t],tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#else
					cnt = pair_op_add(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][b_r_s], seed_single_dm[1][tid][pos2_t],tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#endif

#ifdef	PAIR_TRAVERSE
					for(b_r_s_i = b_r_s + 1; (seed_single_pos[0][tid][b_r_s_i] < posi + devi) && (b_r_s_i < seednn[0]); b_r_s_i++)
					{
						dm12 = seed_single_dm[1][tid][pos2_t] + seed_single_dm[0][tid][b_r_s_i];
						pos_t1 = seed_single_pos[0][tid][b_r_s_i];
						ld1 = seed_single_ld[0][tid][b_r_s_i];
						rd1 = seed_single_rd[0][tid][b_r_s_i];

						memcpy(ref_seq_tmp1[tid], seed_single_refs[0][tid][b_r_s_i], ref_copy_num_chars_1);
#ifdef UNPIPATH_OFF_K20
						cnt = pair_op_add64(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][b_r_s_i], seed_single_dm[1][tid][pos2_t], tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#else
						cnt = pair_op_add(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][b_r_s_i], seed_single_dm[1][tid][pos2_t], tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#endif
					}

					for(b_r_s_i = b_r_s - 1; (seed_single_pos[0][tid][b_r_s_i] > posi - devi) && (b_r_s_i >= 0); b_r_s_i--)
					{
						dm12 = seed_single_dm[1][tid][pos2_t] + seed_single_dm[0][tid][b_r_s_i];
						pos_t1 = seed_single_pos[0][tid][b_r_s_i];
						ld1 = seed_single_ld[0][tid][b_r_s_i];
						rd1 = seed_single_rd[0][tid][b_r_s_i];

						memcpy(ref_seq_tmp1[tid], seed_single_refs[0][tid][b_r_s_i], ref_copy_num_chars_1);
#ifdef UNPIPATH_OFF_K20
						cnt = pair_op_add64(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][b_r_s_i], seed_single_dm[1][tid][pos2_t], tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#else
						cnt = pair_op_add(pos_t1, pos_t2, ld1, rd1, ld2, rd2, dm12, rc_i, seed_single_dm[0][tid][b_r_s_i], seed_single_dm[1][tid][pos2_t], tid, ref_copy_num_chars1, ref_copy_num_chars2, cnt);
#endif
					}
#endif
				}
			}
		}

	}
	if(cnt->v_cnt > 0)
	{
		de_m_p_o[tid] = 0;
#ifdef ALI_OUT
#ifdef	PR_SINGLE
		pr_o_f[tid] = 1;
#endif
		pair_sam_output(tid, read_length1, read_length2, f_cigarn, cnt, seqi, lv_k1, lv_k2, pound_pos1_f_forward, pound_pos1_f_reverse, pound_pos1_r_forward, pound_pos1_r_reverse, pound_pos2_f_forward, pound_pos2_f_reverse, pound_pos2_r_forward, pound_pos2_r_reverse, cir_fix_n - 1);
#endif
	}
}
#ifdef UNPIPATH_OFF_K20
cnt_re* pair_op_add64(uint64_t pos_t1, uint64_t pos_t2, int ld1, int rd1, int ld2, int rd2, int dm12, uint8_t rc_i, int dm1, int dm2, uint8_t tid, uint32_t ref_copy_num_chars1, uint32_t ref_copy_num_chars2, cnt_re* cnt)
#else
cnt_re* pair_op_add(uint32_t pos_t1, uint32_t pos_t2, int ld1, int rd1, int ld2, int rd2, int dm12, uint8_t rc_i, int dm1, int dm2, uint8_t tid, uint32_t ref_copy_num_chars1, uint32_t ref_copy_num_chars2, cnt_re* cnt)
#endif
{
	uint16_t dm_i = 0;
	uint32_t v_cnt = 0;
	uint32_t vs_cnt = 0;
	v_cnt = cnt->v_cnt;
	vs_cnt = cnt->vs_cnt;

	if(dm12 < dm_op[tid])
	{
#ifdef	DM_COPY_REP
		for(dm_i = 0; dm_i < v_cnt; dm_i++)
		{
			ops_vector_pos1[tid][dm_i] = op_vector_pos1[tid][dm_i];
			ops_dm_l1[tid][dm_i] = op_dm_l1[tid][dm_i];
			ops_dm_r1[tid][dm_i] = op_dm_r1[tid][dm_i];

#ifdef	MAPPING_QUALITY
			memcpy(ops_vector_seq1[tid][dm_i], op_vector_seq1[tid][dm_i], ref_copy_num_chars1);
#else
			if(!((op_dm_l1[tid][dm_i] == 0) && (op_dm_r1[tid][dm_i] == 0)))
				memcpy(ops_vector_seq1[tid][dm_i], op_vector_seq1[tid][dm_i], ref_copy_num_chars1);
#endif
			ops_dm_ex1[tid][dm_i] = op_dm_ex1[tid][dm_i];

			ops_vector_pos2[tid][dm_i] = op_vector_pos2[tid][dm_i];
			ops_dm_l2[tid][dm_i] = op_dm_l2[tid][dm_i];
			ops_dm_r2[tid][dm_i] = op_dm_r2[tid][dm_i];
			if(!((op_dm_l2[tid][dm_i] == 0) && (op_dm_r2[tid][dm_i] == 0)))
				memcpy(ops_vector_seq2[tid][dm_i], op_vector_seq2[tid][dm_i], ref_copy_num_chars2);
			ops_dm_ex2[tid][dm_i] = op_dm_ex2[tid][dm_i];

		}
		vs_cnt = v_cnt;
		dm_ops[tid] = dm_op[tid];
#endif
		v_cnt = 0;
		if(rc_i == 0)
		{
			op_vector_pos1[tid][v_cnt] = pos_t1;
			op_vector_pos2[tid][v_cnt] = pos_t2;
#ifdef	MAPPING_QUALITY
			memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
			memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#else
			if(!((ld1 == 0) && (rd1 == 0)))
				memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
			if(!((ld2 == 0) && (rd2 == 0)))
				memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#endif
			op_dm_l1[tid][v_cnt] = ld1;
			op_dm_r1[tid][v_cnt] = rd1;
			op_dm_l2[tid][v_cnt] = ld2;
			op_dm_r2[tid][v_cnt] = rd2;

			op_dm_ex1[tid][v_cnt] = dm1;
			op_dm_ex2[tid][v_cnt] = dm2;
		}
		else
		{
			op_vector_pos1[tid][v_cnt] = pos_t2;
			op_vector_pos2[tid][v_cnt] = pos_t1;
#ifdef	MAPPING_QUALITY
			memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
			memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#else
			if(!((ld2 == 0) && (rd2 == 0)))
				memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
			if(!((ld1 == 0) && (rd1 == 0)))
				memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#endif
			op_dm_l1[tid][v_cnt] = ld2;
			op_dm_r1[tid][v_cnt] = rd2;
			op_dm_l2[tid][v_cnt] = ld1;
			op_dm_r2[tid][v_cnt] = rd1;

			op_dm_ex1[tid][v_cnt] = dm2;
			op_dm_ex2[tid][v_cnt] = dm1;
		}

		op_rc[tid][v_cnt] = rc_i;
		++v_cnt;
		dm_op[tid] = dm12;
	}
	else if(dm12 == dm_op[tid])
	{
		if(v_cnt < cus_max_output_ali)
		{
			if(rc_i == 0)
			{
				op_vector_pos1[tid][v_cnt] = pos_t1;
				op_vector_pos2[tid][v_cnt] = pos_t2;
#ifdef	MAPPING_QUALITY
				memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
				memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#else
				if(!((ld1 == 0) && (rd1 == 0)))
					memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
				if(!((ld2 == 0) && (rd2 == 0)))
					memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#endif
				op_dm_l1[tid][v_cnt] = ld1;
				op_dm_r1[tid][v_cnt] = rd1;
				op_dm_l2[tid][v_cnt] = ld2;
				op_dm_r2[tid][v_cnt] = rd2;

				op_dm_ex1[tid][v_cnt] = dm1;
				op_dm_ex2[tid][v_cnt] = dm2;
			}
			else
			{
				op_vector_pos1[tid][v_cnt] = pos_t2;
				op_vector_pos2[tid][v_cnt] = pos_t1;
#ifdef	MAPPING_QUALITY
				memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
				memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#else
				if(!((ld2 == 0) && (rd2 == 0)))
					memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
				if(!((ld1 == 0) && (rd1 == 0)))
					memcpy(op_vector_seq2[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#endif
				op_dm_l1[tid][v_cnt] = ld2;
				op_dm_r1[tid][v_cnt] = rd2;
				op_dm_l2[tid][v_cnt] = ld1;
				op_dm_r2[tid][v_cnt] = rd1;

				op_dm_ex1[tid][v_cnt] = dm2;
				op_dm_ex2[tid][v_cnt] = dm1;
			}
			op_rc[tid][v_cnt] = rc_i;
			++v_cnt;
		}
	}
	else if(dm12 < dm_ops[tid])
	{
		vs_cnt = 0;
		if(rc_i == 0)
		{
			ops_vector_pos1[tid][vs_cnt] = pos_t1;
			ops_vector_pos2[tid][vs_cnt] = pos_t2;
#ifdef	MAPPING_QUALITY
			memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
			memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#else
			if(!((ld1 == 0) && (rd1 == 0)))
				memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
			if(!((ld2 == 0) && (rd2 == 0)))
				memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#endif
			ops_dm_l1[tid][vs_cnt] = ld1;
			ops_dm_r1[tid][vs_cnt] = rd1;
			ops_dm_l2[tid][vs_cnt] = ld2;
			ops_dm_r2[tid][vs_cnt] = rd2;

			ops_dm_ex1[tid][vs_cnt] = dm1;
			ops_dm_ex2[tid][vs_cnt] = dm2;
		}
		else
		{
			ops_vector_pos1[tid][vs_cnt] = pos_t2;
			ops_vector_pos2[tid][vs_cnt] = pos_t1;
#ifdef	MAPPING_QUALITY
			memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
			memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#else
			if(!((ld2 == 0) && (rd2 == 0)))
				memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
			if(!((ld1 == 0) && (rd1 == 0)))
				memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#endif
			ops_dm_l1[tid][vs_cnt] = ld2;
			ops_dm_r1[tid][vs_cnt] = rd2;
			ops_dm_l2[tid][vs_cnt] = ld1;
			ops_dm_r2[tid][vs_cnt] = rd1;

			ops_dm_ex1[tid][vs_cnt] = dm2;
			ops_dm_ex2[tid][vs_cnt] = dm1;
		}
		ops_rc[tid][vs_cnt] = rc_i;
		++vs_cnt;
		dm_ops[tid] = dm12;
	}
	else if(dm12 == dm_ops[tid])
	{
		if(vs_cnt < cus_max_output_ali)
		{
			if(rc_i == 0)
			{
				ops_vector_pos1[tid][vs_cnt] = pos_t1;
				ops_vector_pos2[tid][vs_cnt] = pos_t2;

#ifdef	MAPPING_QUALITY
				memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
				memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#else
				if(!((ld1 == 0) && (rd1 == 0)))
					memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars1);
				if(!((ld2 == 0) && (rd2 == 0)))
					memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars2);
#endif
				ops_dm_l1[tid][vs_cnt] = ld1;
				ops_dm_r1[tid][vs_cnt] = rd1;
				ops_dm_l2[tid][vs_cnt] = ld2;
				ops_dm_r2[tid][vs_cnt] = rd2;

				ops_dm_ex1[tid][vs_cnt] = dm1;
				ops_dm_ex2[tid][vs_cnt] = dm2;
			}
			else
			{
				ops_vector_pos1[tid][vs_cnt] = pos_t2;
				ops_vector_pos2[tid][vs_cnt] = pos_t1;

#ifdef	MAPPING_QUALITY
				memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
				memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#else
				if(!((ld2 == 0) && (rd2 == 0)))
					memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp2[tid], ref_copy_num_chars1);
				if(!((ld1 == 0) && (rd1 == 0)))
					memcpy(ops_vector_seq2[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars2);
#endif
				ops_dm_l1[tid][vs_cnt] = ld2;
				ops_dm_r1[tid][vs_cnt] = rd2;
				ops_dm_l2[tid][vs_cnt] = ld1;
				ops_dm_r2[tid][vs_cnt] = rd1;

				ops_dm_ex1[tid][vs_cnt] = dm2;
				ops_dm_ex2[tid][vs_cnt] = dm1;
			}
			ops_rc[tid][vs_cnt] = rc_i;
			++vs_cnt;
		}
	}
	cnt->v_cnt = v_cnt;
	cnt->vs_cnt = vs_cnt;
	return cnt;
}

void seed_repetitive_single_end(uint64_t* read_bit, uint8_t tid, uint16_t read_length)
{
	uint8_t re_d = 0;
	uint8_t b_t_n_r = 0;
	uint32_t ref_pos_n = 0;
	uint32_t uni_offset_s_l = 0;
	uint32_t uni_offset_s_r = 0;
	uint16_t read_off = 0;
	uint16_t left_i = 1;
	uint16_t right_i = 1;
	uint16_t r_b_v = 0;
	uint32_t seed_hash = 0;
	uint32_t seed_kmer = 0;
	uint32_t pos_r_ir = 0;
	uint32_t r_ru = 0;
	uint32_t ran_p = 0;
	uint64_t kmer_bit = 0;
	int64_t seed_id_r = 0;
	int64_t seed_binary_r = 0;
	float r_r = 0;

#ifdef UNPIPATH_OFF_K20
	uint64_t kmer_pos_uni = 0;
#else
	uint32_t kmer_pos_uni = 0;
#endif

	r_b_v = 0;
	seedn[tid] = 0;
	seed_l[tid] = pair_ran_intvp;
	for(read_off = 0; read_off <= read_length - k_t; read_off += seed_l[tid])
	{
		if(read_off + k_t - 1 <= r_b_v)	continue;

		re_d = (read_off & 0X1f);

		b_t_n_r = 32 - re_d;

#ifdef HASH_KMER_READ_J
		if(re_d <= re_b)
		{
			kmer_bit = ((read_bit[read_off >> 5] & bit_tran_re[re_d]) >> ((re_b - re_d) << 1));
		}
		else
		{
			kmer_bit = (((read_bit[read_off >> 5] & bit_tran_re[re_d]) << ((re_d - re_b) << 1)) | (read_bit[(read_off >> 5) + 1] >> ((re_2bt - re_d) << 1)));
		}
#else
		tran_tmp_p = (read_bit[read_off >> 5] & bit_tran_re[re_d]);

		//or use this method to deal with: & bit_tran_re[b_t_n_r]
		kmer_bit = (((read_bit[(read_off >> 5) + 1] >> (b_t_n_r << 1)) & bit_tran_re[b_t_n_r]) | (tran_tmp_p << (re_d << 1)));

		kmer_bit >>= re_bt;
#endif

		seed_kmer = (kmer_bit & bit_tran[k_r]);

		seed_hash = (kmer_bit >> (k_r << 1));

		//find the kmer
#ifdef UNPIPATH_OFF_K20
		seed_binary_r = binsearch_offset64(seed_kmer, buffer_kmer_g, buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], buffer_hash_g[seed_hash]);
#else
		seed_binary_r = binsearch_offset(seed_kmer, buffer_kmer_g, buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], buffer_hash_g[seed_hash]);
#endif
		if(seed_binary_r == -1)
			continue;

		//binary search on unipath offset to get the unipathID
		kmer_pos_uni = buffer_off_g[seed_binary_r];//this kmer's offset on unipath
#ifdef UNPIPATH_OFF_K20
		seed_id_r = binsearch_interval_unipath64(kmer_pos_uni, buffer_seqf, result_seqf);
#else
		seed_id_r = binsearch_interval_unipath(kmer_pos_uni, buffer_seqf, result_seqf);
#endif

		ref_pos_n = buffer_pp[seed_id_r + 1] - buffer_pp[seed_id_r];

		uni_offset_s_l = kmer_pos_uni - buffer_seqf[seed_id_r];
		uni_offset_s_r = buffer_seqf[seed_id_r + 1] - (kmer_pos_uni + k_t);

		if(ref_pos_n > pos_n_max)
		{
			if(ref_pos_n <= RANDOM_RANGE)
			{
				for(pos_r_ir = 0; pos_r_ir < ref_pos_n; pos_r_ir++)
					seed_k_pos[tid][seedn[tid]][pos_r_ir] = buffer_p[buffer_pp[seed_id_r] + pos_r_ir] + uni_offset_s_l - read_off;
				seed_k_pos[tid][seedn[tid]][pos_r_ir] = MAXKEY;
				++seedn[tid];
			}
			else
			{
#ifdef	PAIR_RANDOM_SEED
				ran_p = (uint32_t )rand();
				r_r = ((float)ref_pos_n / RANDOM_RANGE_MAX);

				if(r_r >= 2)
				{
					r_ru = (uint32_t )r_r;
					for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE; pos_r_ir++, ran_p++)
						seed_k_pos[tid][seedn[tid]][pos_r_ir] = buffer_p[buffer_pp[seed_id_r] + (pos_r_ir * r_ru) + (random_buffer[ran_p & 0Xfff] % r_ru)] + uni_offset_s_l - read_off;
				}
				else
				{
					for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE; pos_r_ir++)
					{
						seed_k_pos[tid][seedn[tid]][pos_r_ir] = buffer_p[buffer_pp[seed_id_r] + (uint32_t )(seed_r_dup[pos_r_ir] * r_r)] + uni_offset_s_l - read_off;
					}
				}

				seed_k_pos[tid][seedn[tid]][pos_r_ir] = MAXKEY;
				++seedn[tid];

#else
				r_r = ((float)ref_pos_n / RANDOM_RANGE);
				for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE; pos_r_ir++)
					seed_k_pos[tid][seedn[tid]][pos_r_ir] = buffer_p[buffer_pp[seed_id_r] + (uint32_t )(pos_r_ir * r_r)] + uni_offset_s_l - read_off;
				seed_k_pos[tid][seedn[tid]][pos_r_ir] = MAXKEY;
				++seedn[tid];
#endif
			}

			for(left_i = 1; (left_i <= uni_offset_s_l) && (left_i <= read_off); left_i++)
			{
#ifdef UNI_SEQ64
				if(((buffer_seq[(kmer_pos_uni - left_i) >> 5] >> ((31 - ((kmer_pos_uni - left_i) & 0X1f)) << 1)) & 0X3)
				        != ((read_bit[(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3)
				  )	break;
#else
				if(((buffer_seq[(kmer_pos_uni - left_i) >> 2] >> (((kmer_pos_uni - left_i) & 0X3) << 1)) & 0X3)
				        != ((read_bit[(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3)
				  )	break;
#endif
			}

			for(right_i = 1; (right_i <= uni_offset_s_r) && (right_i <= read_length - read_off - k_t); right_i++)
			{
#ifdef UNI_SEQ64
				if(((buffer_seq[(kmer_pos_uni + k_t - 1 + right_i) >> 5] >> ((31 - ((kmer_pos_uni + k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
				        != ((read_bit[(read_off + k_t - 1 + right_i) >> 5] >> ((31 - ((read_off + k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
				  )	break;
#else
				if(((buffer_seq[(kmer_pos_uni + k_t - 1 + right_i) >> 2] >> (((kmer_pos_uni + k_t - 1 + right_i) & 0X3) << 1)) & 0X3)
				        != ((read_bit[(read_off + k_t - 1 + right_i) >> 5] >> ((31 - ((read_off + k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
				  )	break;
#endif
			}

			left_dis[tid][seedn[tid] - 1] = read_off - left_i;
			right_dis[tid][seedn[tid] - 1] = read_off + k_t + right_i - 1;

			r_b_v = read_off + k_t + right_i - 1;
		}
	}
}
#endif


int seed_ali_single_end(int argc, char *argv[])
{
	fprintf(stderr, "single end reads mapping\n\n");

	double t = 0.00;
	clock_t start = 0, end = 0;
	uint16_t v_cnt_i = 0;
	uint32_t r_i = 0;
	uint32_t seqi = 0;
	uint32_t seqii = 0;
	int64_t kr1 = 0;
	int32_t m = 0;
	uint16_t read_length_tmp = 0;


#ifdef	PAIR_RANDOM_SEED
	uint32_t r_dup_i = 0;
	int pos_r_nr = 0;
	uint32_t pos_r_ir = 0;
#endif

	char h_chars[6];

	uint32_t read_in = 0X10000;
	uint8_t readlen_name = 255;

	start = clock();

	load_index_file();

	end = clock();

	t = (double)(end - start) / CLOCKS_PER_SEC;

	seed_l_max = seed_step * cir_fix_n;

	fprintf(stderr, "%lf seconds is used for loading index\n", t);

	fp1 = gzopen(read_fastq1, "r");
	seq1 = kseq_init(fp1);

	if(fp1 == NULL)
	{
		fprintf(stderr, "wrong input file route or name: %s\n", read_fastq1);
		exit(1);
	}

	if(flag_std == 0)
	{
		fp_sam = fopen(sam_result, "w");
		if(fp_sam == NULL)
		{
			fprintf(stderr, "wrong output file route or name: %s\n", sam_result);
			exit(1);
		}
	}
	

	k_r = k_t - k;
	re_b = 32 - k_t;
	re_bt = (re_b << 1);
	re_2bt = 64 - k_t;

	fprintf(stderr, "number of threads running: %u\nmax read length: %u\nseed positions filter: %u\nthe filter number of seed positions: %u\nthe minimum seed interval: %u\nthe number of seed circulation: %u\n", thread_n, readlen_max, pos_n_max, seed_filter_pos_num_singlen, seed_step, cir_fix_n);

	pair_ran_intvp = seed_l_max;

#ifdef	SEED_FILTER_POS
	seed_filter_pos_num_single = seed_filter_pos_num_singlen >> 1;
#endif

	//should allocate at beginning
	seed_num = ((readlen_max - k_t) / seed_l_l) + 1;
	new_seed_cnt = seed_num * pos_n_max;

	uint16_t cigar_max_n = MAX_LV_CIGARCOM;

	start = clock();

#ifdef	PTHREAD_USE

	g_low = (int64_t* )calloc(thread_n, 8);
	r_low = (int64_t* )calloc(thread_n, 8);

	seedm = (seed_m** )calloc(thread_n, sizeof(seed_m* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedm[r_i] = (seed_m* )calloc(seed_num, sizeof(seed_m));

	seedu = (seed_m** )calloc(thread_n, sizeof(seed_m* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedu[r_i] = (seed_m* )calloc(seed_num, sizeof(seed_m));


	seedsets = (seed_sets** )calloc(thread_n, sizeof(seed_sets* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedsets[r_i] = (seed_sets* )calloc(new_seed_cnt, sizeof(seed_sets));

	seed_set_off = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_off[r_i] = (uint32_t* )calloc(new_seed_cnt, sizeof(uint32_t ));

#ifdef	SINGLE_PAIR
	cov_num_front_single = (uint16_t* )calloc(thread_n, sizeof(uint16_t ));
	cov_num_re_single = (uint16_t* )calloc(thread_n, sizeof(uint16_t ));
#endif

#ifdef UNPIPATH_OFF_K20
	seed_set_pos_single = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos_single[r_i] = (uint64_t* )calloc(new_seed_cnt << 1, 8);
#else
	seed_set_pos_single = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_set_pos_single[r_i] = (uint32_t* )calloc(new_seed_cnt << 1, 4);
#endif

	set_pos_n_single = (uint32_t* )calloc(thread_n, 4);
	spa_i_single = (uint32_t* )calloc(thread_n, 4);

	seedpa1_single = (seed_pa_single** )calloc(thread_n, sizeof(seed_pa_single* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seedpa1_single[r_i] = (seed_pa_single* )calloc(new_seed_cnt << 1, sizeof(seed_pa_single ));

	pos_add = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		pos_add[r_i] = (uint8_t* )calloc(pos_n_max, 1);

	fprintf(stderr, "Load seed reduction allocation\n");

#ifdef UNPIPATH_OFF_K20
	op_vector_pos1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_vector_pos1[r_i] = (uint64_t* )calloc(cus_max_output_ali, 8);

	ops_vector_pos1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_vector_pos1[r_i] = (uint64_t* )calloc(cus_max_output_ali, 4);
#else
	op_vector_pos1 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_vector_pos1[r_i] = (uint32_t* )calloc(cus_max_output_ali, 4);

	ops_vector_pos1 = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_vector_pos1[r_i] = (uint32_t* )calloc(cus_max_output_ali, 4);
#endif

	op_dm_l1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_l1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_r1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_r1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_l1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_l1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_r1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_r1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_vector_seq1 = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		op_vector_seq1[r_i] = (uint64_t** )calloc(cus_max_output_ali, sizeof(uint64_t* ));
		for(m = 0; m < cus_max_output_ali; m++)
			if((op_vector_seq1[r_i][m] = (uint64_t* )calloc((((readlen_max - 1) >> 5) + 3), 8)) == NULL)
				exit(1);
	}

	ops_vector_seq1 = (uint64_t*** )calloc(thread_n, sizeof(uint64_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		ops_vector_seq1[r_i] = (uint64_t** )calloc(cus_max_output_ali, sizeof(uint64_t* ));
		for(m = 0; m < cus_max_output_ali; m++)
			if((ops_vector_seq1[r_i][m] = (uint64_t* )calloc((((readlen_max - 1) >> 5) + 3), 8)) == NULL)
				exit(1);
	}

	fprintf(stderr, "Load alignment allocation\n");

#ifdef ALT_ALL
	chr_res = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		chr_res[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

	sam_pos1s = (uint32_t** )calloc(thread_n, sizeof(uint32_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		if((sam_pos1s[r_i] = (uint32_t* )calloc(CUS_MAX_OUTPUT_ALI2, sizeof(uint32_t ))) == NULL)
			exit(1);

	cigar_p1s = (char*** )calloc(thread_n, sizeof(char** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		cigar_p1s[r_i] = (char** )calloc(CUS_MAX_OUTPUT_ALI2, sizeof(char* ));
		for(m = 0; m < CUS_MAX_OUTPUT_ALI2; m++)
			if((cigar_p1s[r_i][m] = (char* )calloc(cigar_max_n, 1)) == NULL)
				exit(1);
	}

	xa_d1s = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		xa_d1s[r_i] = (char* )calloc(CUS_MAX_OUTPUT_ALI2, 1);

	lv_re1s = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		lv_re1s[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, sizeof(int));

#endif

	fprintf(stderr, "Load output allocation\n");

	op_rc = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_rc[r_i] = (uint8_t* )calloc(cus_max_output_ali, 1);

	ops_rc = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_rc[r_i] = (uint8_t* )calloc(cus_max_output_ali, 1);

	ref_seq_tmp1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ref_seq_tmp1[r_i] = (uint64_t* )calloc((((readlen_max - 1) >> 5) + 3), 8);

	cov_a_n_s = (uint8_t* )calloc(thread_n, 1);
	s_uid_f = (uint8_t* )calloc(thread_n, 1);
	seed_l = (int16_t* )calloc(thread_n, 2);//cannot be 0
	max_mismatch = (uint8_t* )calloc(thread_n, 1);
	max_mismatch_p = (uint8_t* )calloc(thread_n, 1);
	dm_op = (int* )calloc(thread_n, 4);
	dm_ops = (int* )calloc(thread_n, 4);
	low_mask = (uint64_t* )calloc(thread_n, 8);

	cigar_m = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		cigar_m[r_i] = (char* )calloc(CIGARMN, 1);

	ali_ref_seq = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ali_ref_seq[r_i] = (char* )calloc(readlen_max + 64, 1);

	read_char = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		read_char[r_i] = (char* )calloc(readlen_max + 32, 1);

	pos_ren[0][0] = (uint16_t* )calloc(thread_n, sizeof(uint16_t ));
	pos_ren[1][1] = (uint16_t* )calloc(thread_n, sizeof(uint16_t ));

	sub_mask = (uint16_t** )calloc(thread_n, sizeof(uint16_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		sub_mask[r_i] = (uint16_t* )calloc(33, 2);

	ex_d_mask = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ex_d_mask[r_i] = (uint64_t* )calloc(33, 8);

	read_bit_1 = (uint64_t** )calloc(thread_n, sizeof(uint64_t* ));
	read_val_1 = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));

	read_val1 = (uint8_t*** )calloc(thread_n, sizeof(uint8_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		read_val1[r_i] = (uint8_t** )calloc(2, sizeof(uint8_t* ));
		read_val1[r_i][0] = (uint8_t* )calloc(readlen_max, 1);
		read_val1[r_i][1] = (uint8_t* )calloc(readlen_max, 1);
	}


	qual_arr_single = (float** )calloc(thread_n, sizeof(float* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		qual_arr_single[r_i] = (float* )calloc(readlen_max, sizeof(float ));
#endif

#ifdef	PAIR_RANDOM_SEED
	seed_r_dup = (uint32_t* )calloc(RANDOM_RANGE, 4);
	srand((unsigned)time(0));
	pos_r_nr = RANDOM_RANGE;

	r_dup_i = 0;
	for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE_MAX; pos_r_ir++)
		if((rand() % (RANDOM_RANGE_MAX - pos_r_ir)) < pos_r_nr)
		{
			seed_r_dup[r_dup_i++] = pos_r_ir;
			pos_r_nr--;
		}

	random_buffer = (uint32_t* )malloc((RAN_CIR) << 2);
	for(pos_r_ir = 0; pos_r_ir < RAN_CIR; pos_r_ir++)
		random_buffer[pos_r_ir] = (uint32_t )rand();
#endif

#ifdef	READN_RANDOM_SEED
	random_buffer_readn = (uint32_t* )calloc(RANDOM_RANGE_READN, 4);

	for(pos_r_ir = 0; pos_r_ir < RANDOM_RANGE_READN; pos_r_ir++)
		random_buffer_readn[pos_r_ir] = (uint32_t )rand();
#endif

	InitiateLVCompute(L, thread_n);

	//for pthread read input
	fprintf(stderr, "Load seq input memory\n");

	seqio = (seq_io* )calloc(read_in, sizeof(seq_io));

	char** read_seq1_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		read_seq1_buffer[r_i] = (char* )calloc(readlen_max, 1);

	char** name_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		name_buffer[r_i] = (char* )calloc(readlen_name, 1);

	qual1_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		qual1_buffer[r_i] = (char* )calloc(readlen_max, 1);

#ifdef	OUTPUT_ARR
	read_rev_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		read_rev_buffer[r_i] = (char* )calloc(readlen_max, 1);

	pr_cigar1_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		pr_cigar1_buffer[r_i] = (char* )calloc(cigar_max_n, 1);

	chr_res_buffer = (int** )calloc(read_in, sizeof(int* ));
	for(r_i = 0; r_i < read_in; r_i++)
		chr_res_buffer[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

	xa_d1s_buffer = (char** )calloc(read_in, sizeof(char* ));
	for(r_i = 0; r_i < read_in; r_i++)
		xa_d1s_buffer[r_i] = (uint8_t* )calloc(CUS_MAX_OUTPUT_ALI2, 1);

	sam_pos1s_buffer = (uint32_t** )calloc(read_in, sizeof(uint32_t* ));
	for(r_i = 0; r_i < read_in; r_i++)
		sam_pos1s_buffer[r_i] = (uint32_t* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

	lv_re1s_buffer = (int** )calloc(read_in, sizeof(int* ));
	for(r_i = 0; r_i < read_in; r_i++)
		lv_re1s_buffer[r_i] = (int* )calloc(CUS_MAX_OUTPUT_ALI2, 4);

#endif

#ifdef	KSW_ALN
	ali_ref_seq2 = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ali_ref_seq2[r_i] = (char* )calloc(readlen_max + 64, 1);

	read_char2 = (char** )calloc(thread_n, sizeof(char* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		read_char2[r_i] = (char* )calloc(readlen_max + 32, 1);

	op_dm_kl1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_kl1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	op_dm_kr1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		op_dm_kr1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_kl1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_kl1[r_i] = (int* )calloc(cus_max_output_ali, 4);

	ops_dm_kr1 = (int** )calloc(thread_n, sizeof(int* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		ops_dm_kr1[r_i] = (int* )calloc(cus_max_output_ali, 4);


	int8_t l, k1;
	mat = (int8_t*)calloc(25, sizeof(int8_t));
	for (l = k1 = 0; l < 4; ++l)
	{
		for (m = 0; m < 4; ++m) mat[k1++] = l == m ? 1 : -4;	/* weight_match : -weight_mismatch */
		mat[k1++] = -1; // ambiguous base
	}
	for (m = 0; m < 5; ++m) mat[k1++] = 0;
#endif

#ifdef	QUAL_FILT_LV_MIS_SINGLE_END
	qual_filt_lv1 = (uint8_t*** )calloc(thread_n, sizeof(uint8_t** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		qual_filt_lv1[r_i] = (uint8_t** )calloc(2, sizeof(uint8_t* ));
		qual_filt_lv1[r_i][0] = (uint8_t* )calloc(readlen_max, 1);
		qual_filt_lv1[r_i][1] = (uint8_t* )calloc(readlen_max, 1);
	}
#endif

#ifdef	ALTER_DEBUG_SINGLE_END
	seed_length_arr = (seed_length_array** )calloc(thread_n, sizeof(seed_length_array* ));
	for(r_i = 0; r_i < thread_n; r_i++)
		seed_length_arr[r_i] = (seed_length_array* )calloc(cus_max_output_ali, sizeof(seed_length_array));

	rep_go = (uint8_t* )calloc(thread_n, 1);

#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
	mp_subs1 = (float*** )calloc(thread_n, sizeof(float** ));
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		mp_subs1[r_i] = (float** )calloc(2, sizeof(float* ));
		mp_subs1[r_i][0] = (float* )calloc(readlen_max, sizeof(float));
		mp_subs1[r_i][1] = (float* )calloc(readlen_max, sizeof(float));
	}
	sub_t = (float* )calloc(thread_n, sizeof(float ));
#endif

	end = clock();

	t = (double)(end - start) / CLOCKS_PER_SEC;

	fprintf(stderr, "%lf seconds is used for allocating memory\n", t);

	fprintf(stderr, "begin reading fastq single-end reads and doing seed reduction and alignment\n");

#ifdef	R_W_LOCK
	pthread_rwlock_init(&rwlock, NULL);
	int seqii_i = 0;
#endif

	double dtime = omp_get_wtime(); //value in seconds

	while (kr1 >= 0)
	{
		for(seqii = 0; (seqii < read_in) && ((kr1 = kseq_read(seq1)) > 0); seqii++)   //(kr1 >= 0)
		{
			strncpy(name_buffer[seqii], seq1->name.s, seq1->name.l);

			if((seq1->name.s)[seq1->name.l - 2] == '/')
				name_buffer[seqii][seq1->name.l - 2] = '\0';
			else	name_buffer[seqii][seq1->name.l] = '\0';

			if(seq1->seq.l > readlen_max)
			{
				seqio[seqii].read_length1 = readlen_max;
				seqio[seqii].length_h1 = seq1->seq.l - readlen_max;
			}
			else
			{
				seqio[seqii].read_length1 = seq1->seq.l;
				seqio[seqii].length_h1 = 0;
			}

			read_length_tmp = seqio[seqii].read_length1;
			strncpy(read_seq1_buffer[seqii], seq1->seq.s, read_length_tmp);
			read_seq1_buffer[seqii][read_length_tmp] = '\0';

			strncpy(qual1_buffer[seqii], seq1->qual.s, read_length_tmp);
			qual1_buffer[seqii][read_length_tmp] = '\0';

			seqio[seqii].read_seq1 = read_seq1_buffer[seqii];
			seqio[seqii].name = name_buffer[seqii];
			seqio[seqii].qual1 = qual1_buffer[seqii];
		}

		if(thread_n <= 1)
		{
#ifdef	R_W_LOCK
			for(seqii_i = 0; seqii_i < seqii; seqii_i++)
				seed_ali_core_single_end(seqii_i, 0);
#else
			seed_ali_core_single_end(seqii, 0);
#endif
		}
		else
		{
			pthread_t* tid;
			thread_ali_t* data;
			pthread_attr_t attr;

			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_ali_t* )calloc(thread_n, sizeof(thread_ali_t ));
			tid = (pthread_t* )calloc(thread_n, sizeof(pthread_t ));

			for(r_i = 0; r_i < thread_n; ++r_i)
			{
				data[r_i].tid = r_i;
				data[r_i].seqn = seqii;
				pthread_create(&tid[r_i], &attr, worker_single_end, data + r_i);
			}
			for(r_i = 0; r_i < thread_n; ++r_i) pthread_join(tid[r_i], 0);

			free(data);
			free(tid);
		}

		//output
#ifdef	OUTPUT_ARR
		for(seqi = 0; seqi < seqii; seqi++)
		{
			if(seqio[seqi].length_h1 == 0)
			{
				fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s\t*\t0\t0\t%s\t%s",
				        seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re], seqio[seqi].pos1,
				        seqio[seqi].qualc1, seqio[seqi].cigar1, seqio[seqi].seq1,
				        seqio[seqi].qual1
				       );

				if(seqio[seqi].xa_n > 0)
				{
					fprintf(fp_sam, "\tXA:Z:");
					for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n; v_cnt_i++)
						fprintf(fp_sam, "%s,%c%u,%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d1s[v_cnt_i], seqio[seqi].sam_pos1s[v_cnt_i], seqio[seqi].cigar_p1s[v_cnt_i], seqio[seqi].lv_re1s[v_cnt_i]);
				}
				fprintf(fp_sam, "\tNM:i:%u\n",seqio[seqi].nm1);
			}
			else
			{
				sprintf(h_chars, "%uH", seqio[seqi].length_h1);

				fprintf(fp_sam, "%s\t%u\t%s\t%"PRId64"\t%d\t%s%s\t*\t0\t0\t%s\t%s",
				        seqio[seqi].name, seqio[seqi].flag1, chr_names[seqio[seqi].chr_re], seqio[seqi].pos1,
				        seqio[seqi].qualc1, seqio[seqi].cigar1, h_chars, seqio[seqi].seq1,
				        seqio[seqi].qual1
				       );

				if(seqio[seqi].xa_n > 0)
				{
					fprintf(fp_sam, "\tXA:Z:");
					for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n; v_cnt_i++)
					{
						fprintf(fp_sam, "%s,%c%u,%s%s,%d;",chr_names[seqio[seqi].chr_res[v_cnt_i]], seqio[seqi].xa_d1s[v_cnt_i], seqio[seqi].sam_pos1s[v_cnt_i], seqio[seqi].cigar_p1s[v_cnt_i], h_chars, seqio[seqi].lv_re1s[v_cnt_i]);
					}
				}
				fprintf(fp_sam, "\tNM:i:%u\n",seqio[seqi].nm1);
			}
		}
#endif

#ifdef	OUTPUT_ARR
		for(seqi = 0; seqi < seqii; seqi++)
		{
			if(seqio[seqi].xa_n > 0)   //
			{
				for(v_cnt_i = 0; v_cnt_i < seqio[seqi].xa_n; v_cnt_i++)
				{
					free(seqio[seqi].cigar_p1s[v_cnt_i]);
				}
			}
		}
#endif

	}

	dtime = omp_get_wtime() - dtime;

	fprintf(stderr, "%lf seconds is used \n", dtime);
	fflush(stdout);

#ifdef	R_W_LOCK
	pthread_rwlock_destroy(&rwlock);
#endif

#ifdef	PTHREAD_USE

	if(g_low != NULL)	free(g_low);
	if(r_low != NULL)	free(r_low);


	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedm[r_i] != NULL)	free(seedm[r_i]);
	if(seedm != NULL)	free(seedm);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedu[r_i] != NULL)	free(seedu[r_i]);
	if(seedu != NULL)	free(seedu);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedsets[r_i] != NULL)	free(seedsets[r_i]);
	if(seedsets != NULL)	free(seedsets);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_set_off[r_i] != NULL)	free(seed_set_off[r_i]);
	if(seed_set_off != NULL)	free(seed_set_off);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_set_pos_single[r_i] != NULL)	free(seed_set_pos_single[r_i]);
	if(seed_set_pos_single != NULL)	free(seed_set_pos_single);

	if(set_pos_n_single != NULL)	free(set_pos_n_single);
	if(spa_i_single != NULL)	free(spa_i_single);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(seedpa1_single[r_i] != NULL)	free(seedpa1_single[r_i]);
	if(seedpa1_single != NULL)	free(seedpa1_single);

#ifdef	SINGLE_PAIR
	if(cov_num_front_single != NULL)	free(cov_num_front_single);
	if(cov_num_re_single != NULL)	free(cov_num_re_single);
#endif

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pos_add[r_i] != NULL)	free(pos_add[r_i]);
	if(pos_add != NULL)	free(pos_add);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_vector_pos1[r_i] != NULL)	free(op_vector_pos1[r_i]);
	if(op_vector_pos1 != NULL)	free(op_vector_pos1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_vector_pos1[r_i] != NULL)	free(ops_vector_pos1[r_i]);
	if(ops_vector_pos1 != NULL)	free(ops_vector_pos1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_l1[r_i] != NULL)	free(op_dm_l1[r_i]);
	if(op_dm_l1 != NULL)	free(op_dm_l1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_r1[r_i] != NULL)	free(op_dm_r1[r_i]);
	if(op_dm_r1 != NULL)	free(op_dm_r1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_l1[r_i] != NULL)	free(ops_dm_l1[r_i]);
	if(ops_dm_l1 != NULL)	free(ops_dm_l1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_r1[r_i] != NULL)	free(ops_dm_r1[r_i]);
	if(ops_dm_r1 != NULL)	free(ops_dm_r1);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < cus_max_output_ali; m++)
			if(op_vector_seq1[r_i][m] != NULL)	free(op_vector_seq1[r_i][m]);
		if(op_vector_seq1[r_i] != NULL)	free(op_vector_seq1[r_i]);
	}
	if(op_vector_seq1 != NULL)	free(op_vector_seq1);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < cus_max_output_ali; m++)
			if(ops_vector_seq1[r_i][m] != NULL)	free(ops_vector_seq1[r_i][m]);
		if(ops_vector_seq1[r_i] != NULL)	free(ops_vector_seq1[r_i]);
	}
	if(ops_vector_seq1 != NULL)	free(ops_vector_seq1);

#ifdef ALT_ALL

	for(r_i = 0; r_i < thread_n; r_i++)
		if(chr_res[r_i] != NULL)	free(chr_res[r_i]);
	if(chr_res != NULL)	free(chr_res);


	for(r_i = 0; r_i < thread_n; r_i++)
		if(sam_pos1s[r_i] != NULL)	free(sam_pos1s[r_i]);
	if(sam_pos1s != NULL)	free(sam_pos1s);


	for(r_i = 0; r_i < thread_n; r_i++)
	{
		for(m = 0; m < CUS_MAX_OUTPUT_ALI2; m++)
			if(cigar_p1s[r_i][m] != NULL)	free(cigar_p1s[r_i][m]);
		if(cigar_p1s[r_i] != NULL)	free(cigar_p1s[r_i]);
	}
	if(cigar_p1s != NULL)	free(cigar_p1s);


	for(r_i = 0; r_i < thread_n; r_i++)
		if(xa_d1s[r_i] != NULL)	free(xa_d1s[r_i]);
	if(xa_d1s != NULL)	free(xa_d1s);


	for(r_i = 0; r_i < thread_n; r_i++)
		if(lv_re1s[r_i] != NULL)	free(lv_re1s[r_i]);
	if(lv_re1s != NULL)	free(lv_re1s);

#endif

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_rc[r_i] != NULL)	free(op_rc[r_i]);
	if(op_rc != NULL)	free(op_rc);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_rc[r_i] != NULL)	free(ops_rc[r_i]);
	if(ops_rc != NULL)	free(ops_rc);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ref_seq_tmp1[r_i] != NULL)	free(ref_seq_tmp1[r_i]);
	if(ref_seq_tmp1 != NULL)	free(ref_seq_tmp1);

	if(cov_a_n_s != NULL)	free(cov_a_n_s);
	if(s_uid_f != NULL)	free(s_uid_f);
	if(max_mismatch != NULL)	free(max_mismatch);
	if(max_mismatch_p != NULL)	free(max_mismatch_p);
	if(dm_op != NULL)	free(dm_op);
	if(dm_ops != NULL)	free(dm_ops);
	if(low_mask != NULL)	free(low_mask);


	for(r_i = 0; r_i < thread_n; r_i++)
		if(cigar_m[r_i] != NULL)	free(cigar_m[r_i]);
	if(cigar_m != NULL)	free(cigar_m);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ali_ref_seq[r_i] != NULL)	free(ali_ref_seq[r_i]);
	if(ali_ref_seq != NULL)	free(ali_ref_seq);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(read_char[r_i] != NULL)	free(read_char[r_i]);
	if(read_char != NULL)	free(read_char);

	if(pos_ren[0][0] != NULL)	free(pos_ren[0][0]);
	if(pos_ren[1][1] != NULL)	free(pos_ren[1][1]);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(sub_mask[r_i] != NULL)	free(sub_mask[r_i]);
	if(sub_mask != NULL)	free(sub_mask);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ex_d_mask[r_i] != NULL)	free(ex_d_mask[r_i]);
	if(ex_d_mask != NULL)	free(ex_d_mask);

	if(read_bit_1 != NULL)	free(read_bit_1);
	if(read_val_1 != NULL)	free(read_val_1);

	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if(read_val1[r_i][0] != NULL)	free(read_val1[r_i][0]);
		if(read_val1[r_i][1] != NULL)	free(read_val1[r_i][1]);
		if(read_val1[r_i] != NULL)	free(read_val1[r_i]);
	}
	if(read_val1 != NULL)	free(read_val1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(qual_arr_single[r_i] != NULL)	free(qual_arr_single[r_i]);
	if(qual_arr_single != NULL)	free(qual_arr_single);

#endif

#ifdef	PAIR_RANDOM_SEED

	if(seed_r_dup != NULL)	free(seed_r_dup);
	if(random_buffer != NULL)	free(random_buffer);

#endif

#ifdef	READN_RANDOM_SEED
	if(random_buffer_readn)	free(random_buffer_readn);
#endif

	if(seqio != NULL)	free(seqio);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(read_seq1_buffer[r_i] != NULL)	free(read_seq1_buffer[r_i]);
	if(read_seq1_buffer != NULL)	free(read_seq1_buffer);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(name_buffer[r_i] != NULL)	free(name_buffer[r_i]);
	if(name_buffer != NULL)	free(name_buffer);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(qual1_buffer[r_i] != NULL)	free(qual1_buffer[r_i]);
	if(qual1_buffer != NULL)	free(qual1_buffer);


#ifdef	OUTPUT_ARR
	for(r_i = 0; r_i < thread_n; r_i++)
		if(read_rev_buffer[r_i] != NULL)	free(read_rev_buffer[r_i]);
	if(read_rev_buffer != NULL)	free(read_rev_buffer);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(pr_cigar1_buffer[r_i] != NULL)	free(pr_cigar1_buffer[r_i]);
	if(pr_cigar1_buffer != NULL)	free(pr_cigar1_buffer);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(chr_res_buffer[r_i] != NULL)	free(chr_res_buffer[r_i]);
	if(chr_res_buffer != NULL)	free(chr_res_buffer);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(xa_d1s_buffer[r_i] != NULL)	free(xa_d1s_buffer[r_i]);
	if(xa_d1s_buffer != NULL)	free(xa_d1s_buffer);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(sam_pos1s_buffer[r_i] != NULL)	free(sam_pos1s_buffer[r_i]);
	if(sam_pos1s_buffer != NULL)	free(sam_pos1s_buffer);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(lv_re1s_buffer[r_i] != NULL)	free(lv_re1s_buffer[r_i]);
	if(lv_re1s_buffer != NULL)	free(lv_re1s_buffer);

#endif

#ifdef	KSW_ALN
	for(r_i = 0; r_i < thread_n; r_i++)
		if(ali_ref_seq2[r_i] != NULL)	free(ali_ref_seq2[r_i]);
	if(ali_ref_seq2 != NULL)	free(ali_ref_seq2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(read_char2[r_i] != NULL)	free(read_char2[r_i]);
	if(read_char2 != NULL)	free(read_char2);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_kl1[r_i] != NULL)	free(op_dm_kl1[r_i]);
	if(op_dm_kl1 != NULL)	free(op_dm_kl1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(op_dm_kr1[r_i] != NULL)	free(op_dm_kr1[r_i]);
	if(op_dm_kr1 != NULL)	free(op_dm_kr1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_kl1[r_i] != NULL)	free(ops_dm_kl1[r_i]);
	if(ops_dm_kl1 != NULL)	free(ops_dm_kl1);

	for(r_i = 0; r_i < thread_n; r_i++)
		if(ops_dm_kr1[r_i] != NULL)	free(ops_dm_kr1[r_i]);
	if(ops_dm_kr1 != NULL)	free(ops_dm_kr1);

	free(mat);
#endif

#ifdef	QUAL_FILT_LV_MIS_SINGLE_END
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if(qual_filt_lv1[r_i][0] != NULL)	free(qual_filt_lv1[r_i][0]);
		if(qual_filt_lv1[r_i][1] != NULL)	free(qual_filt_lv1[r_i][1]);
		if(qual_filt_lv1[r_i] != NULL)	free(qual_filt_lv1[r_i]);
	}
	if(qual_filt_lv1 != NULL)	free(qual_filt_lv1);
#endif

#ifdef	ALTER_DEBUG_SINGLE_END
	for(r_i = 0; r_i < thread_n; r_i++)
		if(seed_length_arr[r_i] != NULL)	free(seed_length_arr[r_i]);
	if(seed_length_arr != NULL)	free(seed_length_arr);

	if(rep_go != NULL)	free(rep_go);
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
	for(r_i = 0; r_i < thread_n; r_i++)
	{
		if(mp_subs1[r_i][0] != NULL)	free(mp_subs1[r_i][0]);
		if(mp_subs1[r_i][1] != NULL)	free(mp_subs1[r_i][1]);
		if(mp_subs1[r_i] != NULL)	free(mp_subs1[r_i]);
	}

	if(sub_t != NULL)	free(sub_t);
#endif

	fclose(fp_sam);

	kseq_destroy(seq1);
	gzclose(fp1);

	return 0;
}
#ifdef	R_W_LOCK
int seed_ali_core_single_end(int read_seq_core, uint8_t tid)
#else
int seed_ali_core_single_end(uint32_t seqn, uint8_t tid)
#endif
{
	uint16_t dm_i = 0;
	uint8_t extension_stop = 0;
	char cigar_p1[MAX_LV_CIGAR];
	uint8_t rc_i = 0;
	uint8_t lv_ref_length_re = 0;
	uint8_t q_n1 = 0;
	uint8_t c_m_f = 0;
	uint8_t end1_uc_f = 0;
	uint8_t nuc1_f = 0;
	uint8_t b_t_n_r = 0;
	uint8_t re_d = 0;
	uint16_t f_cigar[MAX_LV_CIGAR];
	char sam_seq1[MAX_READLEN + 1] = {};
	char cigarBuf1[MAX_LV_CIGAR] = {};
	char cigarBuf2[MAX_LV_CIGAR] = {};
	char str_o[MAX_LV_CIGAR];
	char b_cigar[MAX_LV_CIGAR];
	char* pch = NULL;
	char* saveptr = NULL;

	uint16_t ref_copy_num = 0;
	uint16_t f_cigarn = 0;
	uint16_t read_bit_char = 0;
	uint16_t rst_i = 0;
	uint16_t s_m_t = 0;
	uint16_t read_b_i = 0;
	uint16_t f_c = 0;
	uint16_t pchl = 0;
	uint16_t f_i = 0;
	uint16_t s_o = 0;
	uint16_t snt = 0;
	uint16_t d_n1 = 0;
	uint16_t i_n1 = 0;
	uint16_t read_length = 0;
	uint16_t lv_k = 0;

	uint32_t read_length_a = 0;
	uint32_t seqi = 0;
	uint32_t v_cnt = 0;
	uint32_t vs_cnt = 0;
	uint32_t ref_copy_num_chars = 0;
	uint32_t r_i = 0;
	uint32_t psp_i = 0;
	uint32_t d_l1 = 0;
	uint32_t d_r1 = 0;
	uint32_t mis_c_n = 0;
	uint32_t xa_i = 0;
	uint32_t v_cnt_i = 0;
	uint32_t va_cnt_i = 0;
	uint32_t sam_flag1 = 0;
	uint32_t seed1_i = 0;
	uint32_t sam_seq_i = 0;
	uint32_t max_single_score = 0;
	uint32_t seed_posn_filter_single = 0;

	int s_r_o_l1 = 0;
	int	s_r_o_r1 = 0;
	int lv_re1f = 0;
	int sam_qual1 = 0;
	int chr_re = 0;
	int mid = 0;
	int low = 0;
	int high = 0;
	int m_m_n = 0;
	int sn = 0;
	int bit_char_i = 0;
	int dm1 = 0;
	int dm_l1 = 0;
	int dm_r1 = 0;
	int dmt1 = 0;
	int lv_dmt1 = 0;
	int ld1 = 0;
	int rd1 = 0;
	int cmp_re = 0;
	int q_rear_i = 0;
	int q_rear1 = 0;
	int cache_dml1[MAX_Q_NUM];
	int cache_dmr1[MAX_Q_NUM];
	int cache_dis1[MAX_Q_NUM];
	int cache_kl1[MAX_Q_NUM];
	int cache_kr1[MAX_Q_NUM];

	uint16_t off_i = 0;
	uint16_t max_sets_n = 0;
	uint16_t max_seed_length[2];
	uint64_t c_tmp = 0;
	uint64_t xor_tmp = 0;
	uint64_t ref_tmp_ori = 0;
	uint64_t ref_tmp_ori2 = 0;
	uint64_t tran_tmp_p = 0;
	uint64_t tran_tmp = 0;
	int64_t sam_pos1 = 0;
	uint64_t cache_end1[MAX_Q_NUM][MAX_REF_SEQ_C];

	seed_pa_single* seed_pr1 = NULL;

	uint16_t s_offset1 = 0;
	int16_t s_r_o_l = 0;
	int16_t s_r_o_r = 0;

#ifdef	KSW_ALN
	int band_with = 33;//100
	int zdrop = 0;//100
	int end_bonus = 5;
	int tle;
	int gtle;
	int gscore;
	int max_off;
	int x_score, y_score, op_score, len_score, op_score1, len_score1, op_score2, len_score2;
	int k_start1 = 0;
	int k_start2 = 0;
	int k_middle = 0;
	uint32_t* cigar1 = NULL;
	uint32_t* cigar2 = NULL;
	int n_cigar1 = 0;
	int n_cigar2 = 0;
	int nm_score = 0;
#endif

#ifdef	QUAL_FILT_SINGLE_END
	uint64_t* qual_filt_1 = NULL;
	uint8_t qual_filt_fix_single_end = 55;//
#endif

#ifdef	QUAL_FILT_LV_MIS_SINGLE_END
	uint8_t* qual_filt_lv_1 = NULL;
	uint8_t* qual_filt_lv_1_o = NULL;
#endif

#ifdef	ALTER_DEBUG_SINGLE_END
	uint16_t seed_length1 = 0;
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
	float log_tmp = 0;
	float* mp_subs_1 = NULL;
	float* mp_subs_1_o = NULL;
	float m_sub_tmp = 0;
#endif

	uint8_t mp_flag = 0;
	uint8_t cir_n = 0;

	uint16_t mis_c_n_filt = 0;
	int16_t dm_cir_min = 0;

#ifdef UNPIPATH_OFF_K20
	uint64_t posi = 0;
	uint64_t pos_l = 0;
	uint64_t x = 0;
#else
	uint32_t posi = 0;
	uint32_t pos_l = 0;
	uint32_t x = 0;
#endif

#ifdef	READN_RANDOM_SEED
	char tmp_char;
	uint16_t readn_re1 = 0;
	uint16_t readn_re2 = 0;
#endif

#ifdef	CIGAR_LEN_ERR
	int cigar_len = 0;
	int cigar_len_tmp = 0;
	int cigar_len_re = 0;
	uint16_t pchl_tmp = 0;
	uint16_t s_o_tmp = 0;
	char* pch_tmp = NULL;
	char* saveptr_tmp = NULL;
	char cigar_tmp[MAX_LV_CIGARCOM];
#endif

#ifndef	R_W_LOCK
	for(seqi = 0; seqi < seqn; seqi++)
#else
	if(1)
#endif
	{
#ifndef	R_W_LOCK

#ifdef PTHREAD_USE
		if ((seqi % thread_n) != tid) continue;
#endif

#else
		seqi = read_seq_core;
#endif
		read_length = seqio[seqi].read_length1;

		lv_k = (read_length * max_single_score_r) + 1;

		read_length_a = read_length - 1;

#ifdef	SINGLE_PAIR
		cov_num_front_single[tid] = (read_length_a >> 7);
		cov_num_re_single[tid] = 63 - ((read_length_a >> 1) & 0X3f);
#endif

		f_cigarn = read_length;

		f_cigar[f_cigarn] = '\0';

		//sprintf(cigar_m[tid],"%uM\0",read_length);
		sprintf(cigar_m[tid],"%uM",read_length);

		//for bit operation
		read_bit_char = (((uint16_t )((read_length_a >> 5) + 1)) << 3);

		ref_copy_num = ((read_length_a) >> 5) + 3;
		ref_copy_num_chars = (ref_copy_num << 3);
		lv_ref_length_re = (read_length & 0X1f);

		for(r_i = 0; r_i <= 32 - lv_ref_length_re; r_i++)
		{
			ex_d_mask[tid][r_i] = bit_assi[lv_ref_length_re + r_i];
			sub_mask[tid][r_i] = ref_copy_num - 2;
		}

		for(r_i = 33 - lv_ref_length_re; r_i <= 32; r_i++)
		{
			ex_d_mask[tid][r_i] = bit_assi[lv_ref_length_re + r_i - 32];
			sub_mask[tid][r_i] = ref_copy_num - 1;
		}

		low_mask[tid] = bit_assi[lv_ref_length_re];

#ifdef	QUAL_FILT_SINGLE_END
		memset(qual_filt1[tid][0], 0, read_bit_char);
		memset(qual_filt1[tid][1], 0, read_bit_char);

		c_tmp = 3;
		for(r_i = 0; r_i < read_length; r_i++)
		{
			if(seqio[seqi].qual1[r_i] < qual_filt_fix_single_end)   //'7': 55 63
			{
				qual_filt1[tid][0][r_i >> 5] |= (((uint64_t )c_tmp) << ((31 - (r_i & 0X1f)) << 1));
				qual_filt1[tid][1][(read_length_a - r_i) >> 5] |= (((uint64_t )c_tmp) << ((31 - ((read_length_a - r_i) & 0X1f)) << 1));
			}
		}

		for(r_i = 0; r_i < (read_length >> 5) + 1; r_i++)
			qual_filt1[tid][0][r_i] = ~qual_filt1[tid][0][r_i];
		for(r_i = 0; r_i < (read_length >> 5) + 1; r_i++)
			qual_filt1[tid][1][r_i] = ~qual_filt1[tid][1][r_i];

#endif

#ifdef	QUAL_FILT_LV_MIS_SINGLE_END
		for(r_i = 0; r_i < read_length; r_i++)
		{
			if(seqio[seqi].qual1[r_i] < qual_filt_fix_single_end)	qual_filt_lv1[tid][0][r_i] = 0;
			else	qual_filt_lv1[tid][0][r_i] = 3;
		}
		for(r_i = 0; r_i < read_length; r_i++)
		{
			if(seqio[seqi].qual1[r_i] < qual_filt_fix_single_end)	qual_filt_lv1[tid][1][read_length_a - r_i] = 0;
			else	qual_filt_lv1[tid][1][read_length_a - r_i] = 3;
		}
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
		for(r_i = 0; r_i < read_length; r_i++)
		{
			m_sub_tmp = mp_sub_bp[seqio[seqi].qual1[r_i]];
			mp_subs1[tid][0][r_i] = m_sub_tmp;
			mp_subs1[tid][1][read_length_a - r_i] = m_sub_tmp;
		}
#endif

		max_single_score = (uint32_t )(((float )read_length) * max_single_score_r);
		cov_a_n_s[tid] = ((read_length_a) >> 6) + 1;
		max_mismatch[tid] = (uint8_t )(((float )read_length) * mis_match_r);
		max_mismatch_p[tid] = max_mismatch[tid] - 1;

		memset(read_bit1[tid][0], 0, read_bit_char);
		memset(read_bit1[tid][1], 0, read_bit_char);

		r_i = 0;
		while ((seqio[seqi].read_seq1)[r_i])
		{
#ifdef	READN_RANDOM_SEED
			tmp_char = (seqio[seqi].read_seq1)[r_i];
			if(tmp_char == 'N')
			{
				readn_re1 = readn_cnt & 0X3ff;
				readn_re2 = readn_cnt & 0Xf;
				c_tmp = ((random_buffer_readn[readn_re1] >> (readn_re2 << 1)) & 0X3);
				readn_cnt++;
			}
			else	c_tmp = charToDna5n[tmp_char];
#else
			c_tmp = charToDna5n[(seqio[seqi].read_seq1)[r_i]];
#endif

#ifdef ALI_LV
			read_val1[tid][0][r_i] = c_tmp;
			read_val1[tid][1][read_length_a - r_i] = c_tmp ^ 0X3;
#endif

			read_bit1[tid][0][r_i >> 5] |= (((uint64_t )c_tmp) << ((31 - (r_i & 0X1f)) << 1));
			read_bit1[tid][1][(read_length_a - r_i) >> 5] |= (((uint64_t )(c_tmp ^ 0X3)) << ((31 - ((read_length_a - r_i) & 0X1f)) << 1));

			r_i++;
		}

		seed_l[tid] = seed_l_max;

#ifdef	ALTER_DEBUG_SINGLE_END
		rep_go[tid] = 1;
#endif

		cir_n = 1;
		while(cir_n <= cir_fix_n)
		{
			cir_cnt++;

			dm_op[tid] = MAX_OP_SCORE;
			dm_ops[tid] = MAX_OP_SCORE;
			v_cnt = 0;
			vs_cnt = 0;
			dm_cir_min = 0Xfff;

			spa_i_single[tid] = 0;
			set_pos_n_single[tid] = 0;

			for(rc_i = 0; rc_i < 2; rc_i++)
			{
				max_seed_length[rc_i] = 0;
#ifdef UNPIPATH_OFF_K20
				single_seed_reduction_core_single64(seedpa1_single[tid], read_bit1[tid][rc_i], read_val1[tid][rc_i], &(seed_set_pos_single[tid]), tid, read_length, &(max_seed_length[rc_i]), rc_i);
#else
				single_seed_reduction_core_single(seedpa1_single[tid], read_bit1[tid][rc_i], read_val1[tid][rc_i], &(seed_set_pos_single[tid]), tid, read_length, &(max_seed_length[rc_i]), rc_i);
#endif
			}

			if(spa_i_single[tid] == 0)
			{
				seed_l[tid] -= seed_step;
			}
			else
			{
				seed_posn_filter_single = 0;
				extension_stop = 0;

				seed_pr1 = seedpa1_single[tid];
				qsort(seed_pr1, spa_i_single[tid], sizeof(seed_pa_single), compare_plen_single);

				max_sets_n = (max_seed_length[0] > max_seed_length[1]) ? max_seed_length[0]:max_seed_length[1];
				for(off_i = 0; off_i < spa_i_single[tid]; off_i++)
					if(seed_pr1[off_i].length < max_sets_n >> 1)
					{
						++off_i;
						break;
					}

				for(seed1_i = 0; (seed1_i < off_i) && (extension_stop == 0); seed1_i++)   //
				{
					rc_i = seed_pr1[seed1_i].rc;
					read_bit_1[tid] = read_bit1[tid][rc_i];

#ifdef	QUAL_FILT_SINGLE_END
					qual_filt_1 = qual_filt1[tid][rc_i];
#endif

#ifdef	QUAL_FILT_LV_MIS_SINGLE_END
					qual_filt_lv_1 = qual_filt_lv1[tid][rc_i];
					qual_filt_lv_1_o = qual_filt_lv1[tid][1 - rc_i];
#endif
					if(seed_pr1[seed1_i].ui == 1)
					{
						end1_uc_f = 0;
						nuc1_f = 0;

						d_l1 = seed_pr1[seed1_i].ref_pos_off;
						d_r1 = seed_pr1[seed1_i].ref_pos_off_r;
					}
					else
					{
						end1_uc_f = 1;
					}

					q_rear1 = 0;
					q_n1 = 0;

					dmt1 = ali_exl;
					lv_dmt1 = lv_k;

					s_r_o_l1 = seed_pr1[seed1_i].s_r_o_l;
					s_r_o_r1 = seed_pr1[seed1_i].s_r_o_r;
#ifdef	ALTER_DEBUG_SINGLE_END
					seed_length1 = seed_pr1[seed1_i].length;
#endif

					for(psp_i = 0; psp_i < seed_pr1[seed1_i].pos_n; psp_i++)
					{
						if(seed_posn_filter_single > seed_filter_pos_num_singlen)	extension_stop = 1;//break;
						seed_posn_filter_single++;

						posi = seed_set_pos_single[tid][seed_pr1[seed1_i].pos_start + psp_i];

						extension_cnt++;

						//for end1
						if((end1_uc_f == 0) && ((dmt1 <= d_l1) && (dmt1 <= d_r1)))   // && (nuc1_f == 0)
						{
							if(nuc1_f == 0)
							{
								pos_l = posi - max_extension_length - 1;
								re_d = pos_l & 0X1f;
								b_t_n_r = 32 - re_d;

								if(re_d != 0)
								{
									tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);
									memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars);

									for(rst_i = 0; rst_i < ref_copy_num - 1; rst_i++)
									{
										tran_tmp = (ref_seq_tmp1[tid][rst_i] & bit_tran_re[re_d]);

										ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
										ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
										tran_tmp_p = tran_tmp;
									}
									ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
									ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

									//clear the lowest n bit
									ref_seq_tmp1[tid][rst_i] &= low_mask[tid];

								}
								else
								{
									memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars);

									//clear the lowest n bit
									ref_seq_tmp1[tid][ref_copy_num - 1] &= low_mask[tid];
								}

								//exact match
								mis_c_n = 0;
#ifdef	QUAL_FILT_SINGLE_END
								mis_c_n_filt = 0;
#endif
								ref_tmp_ori = ref_seq_tmp1[tid][ref_copy_num - 2];
								ref_seq_tmp1[tid][ref_copy_num - 2] &= low_mask[tid];

								for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num - 1; rst_i++, read_b_i++)
								{
									xor_tmp = ref_seq_tmp1[tid][rst_i] ^ read_bit_1[tid][read_b_i];
#ifdef	QUAL_FILT_SINGLE_END
									mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
									xor_tmp &= qual_filt_1[read_b_i];
									mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
									if(mis_c_n_filt > max_mismatch[tid])	break;
#else
									mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
									mis_c_n_filt = mis_c_n;
									if(mis_c_n > max_mismatch[tid])	break;
#endif
								}

								ref_seq_tmp1[tid][ref_copy_num - 2] = ref_tmp_ori;

								//lv
								if(mis_c_n_filt > max_mismatch[tid])
								{
									if((cir_n == cir_fix_n) && (local_ksw))
									{
#ifdef	KSW_ALN
										for(bit_char_i = s_r_o_l1, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_l1, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										ksw_extend(s_r_o_l1 + 1, read_char[tid], 33 + s_r_o_l1, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r1 - s_r_o_l1, &dm_l1, &tle, &gtle, &gscore, &max_off);

										for(bit_char_i = s_r_o_r1, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_r1, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										ksw_extend(read_length - s_r_o_r1, read_char[tid], 32 + read_length - s_r_o_r1, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r1 - s_r_o_l1, &dm_r1, &tle, &gtle, &gscore, &max_off);

										dm1 = MAX_OP_SCORE - (dm_l1 + dm_r1);

										ld1 = s_r_o_l1;
										rd1 = s_r_o_r1;
#endif
									}
									else
									{
#ifdef SPLIT_LV_SINGLE

#ifdef	QUAL_FILT_LV_MIS_SINGLE_END
										for(bit_char_i = s_r_o_l1, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_l1, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_l1 = computeEditDistance_mis(ali_ref_seq[tid], 33 + s_r_o_l1, read_char[tid], s_r_o_l1 + 1, lv_dmt1, L[tid], qual_filt_lv_1_o + read_length_a - s_r_o_l1);

										for(bit_char_i = s_r_o_r1, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_r1, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_r1 = computeEditDistance_mis(ali_ref_seq[tid], 32 + read_length - s_r_o_r1, read_char[tid], read_length - s_r_o_r1, lv_dmt1, L[tid], qual_filt_lv_1 + s_r_o_r1);

										dm1 = dm_l1 + dm_r1;

#else
										//split lv
										for(bit_char_i = s_r_o_l1, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_l1, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_l1 = computeEditDistance(ali_ref_seq[tid], 33 + s_r_o_l1, read_char[tid], s_r_o_l1 + 1, lv_dmt1, L[tid]);

										for(bit_char_i = s_r_o_r1, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_r1, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_r1 = computeEditDistance(ali_ref_seq[tid], 32 + read_length - s_r_o_r1, read_char[tid], read_length - s_r_o_r1, lv_dmt1, L[tid]);

										dm1 = dm_l1 + dm_r1;
#endif
										ld1 = s_r_o_l1;
										rd1 = s_r_o_r1;
#endif
									}
								}
								else
								{
									dm1 = mis_c_n;
									ld1 = 0;
									rd1 = 0;
									dm_l1 = 0;
									dm_r1 = 0;
								}

								//these two values could be different
								if((dm_l1 != -1) && (dm_r1 != -1))   //need to be modified
								{
									if(dm1 < dmt1)	dmt1 = dm1 + 1;

									if(dm1 < lv_dmt1)	lv_dmt1 = dm1 + 1;

									if(dm1 < max_mismatch_p[tid])
									{
										dmt1 = 0;
									}
								}
								else
								{

									dm1 = MAX_EDIT_SCORE;
								}
								nuc1_f = 1;
							}
						}
						else
						{
							pos_l = posi - max_extension_length - 1;
							re_d = pos_l & 0X1f;
							b_t_n_r = 32 - re_d;

							if(re_d != 0)
							{
								tran_tmp_p = (buffer_ref_seq[pos_l >> 5] & bit_tran_re[re_d]);
								memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5) + 1, ref_copy_num_chars);

								for(rst_i = 0; rst_i < ref_copy_num - 1; rst_i++)
								{
									tran_tmp = (ref_seq_tmp1[tid][rst_i] & bit_tran_re[re_d]);

									ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
									ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));
									tran_tmp_p = tran_tmp;
								}
								ref_seq_tmp1[tid][rst_i] >>= (b_t_n_r << 1);
								ref_seq_tmp1[tid][rst_i] |= (tran_tmp_p << (re_d << 1));

								//clear the lowest n bit
								ref_seq_tmp1[tid][rst_i] &= low_mask[tid];

							}
							else
							{
								memcpy(ref_seq_tmp1[tid], buffer_ref_seq + (pos_l >> 5), ref_copy_num_chars);

								//clear the lowest n bit
								ref_seq_tmp1[tid][ref_copy_num - 1] &= low_mask[tid];
							}

							//trim the beginning and end of the current ref seq based on current minimum edit distance dm_t
							s_m_t = sub_mask[tid][dmt1];
							ref_seq_tmp1[tid][0] &= bit_tran[dmt1];
							ref_seq_tmp1[tid][s_m_t] &= ex_d_mask[tid][dmt1];

							//traverse and check whether there is an existing seq that is as same as current new ref seq
							c_m_f = 0;
							for(q_rear_i = q_rear1 - 1; q_rear_i >= 0; q_rear_i--)
							{
								ref_tmp_ori = cache_end1[q_rear_i][0];
								cache_end1[q_rear_i][0] &= bit_tran[dmt1];

								ref_tmp_ori2 = cache_end1[q_rear_i][s_m_t];
								cache_end1[q_rear_i][s_m_t] &= ex_d_mask[tid][dmt1];

								cmp_re = memcmp(cache_end1[q_rear_i], ref_seq_tmp1[tid], (s_m_t + 1) << 3);

								cache_end1[q_rear_i][0] = ref_tmp_ori;
								cache_end1[q_rear_i][s_m_t] = ref_tmp_ori2;

								if(cmp_re == 0)
								{
									//deal with finding an alignment
									dm1 = cache_dis1[q_rear_i];
									ld1 = cache_dml1[q_rear_i];
									rd1 = cache_dmr1[q_rear_i];
#ifdef	KSW_ALN
									dm_l1 = cache_kl1[q_rear_i];
									dm_r1 = cache_kr1[q_rear_i];
#endif
									c_m_f = 1;
									break;

								}
							}

							if((q_n1 > MAX_Q_NUM) && (q_rear_i < 0))
							{
								for(q_rear_i = MAX_Q_NUM - 1; q_rear_i >= q_rear1; q_rear_i--)
								{
									ref_tmp_ori = cache_end1[q_rear_i][0];
									cache_end1[q_rear_i][0] &= bit_tran[dmt1];

									ref_tmp_ori2 = cache_end1[q_rear_i][s_m_t];
									cache_end1[q_rear_i][s_m_t] &= ex_d_mask[tid][dmt1];

									cmp_re = memcmp(cache_end1[q_rear_i], ref_seq_tmp1[tid], (s_m_t + 1) << 3);

									cache_end1[q_rear_i][0] = ref_tmp_ori;
									cache_end1[q_rear_i][s_m_t] = ref_tmp_ori2;

									if(cmp_re == 0)
									{
										//deal with finding an alignment

										dm1 = cache_dis1[q_rear_i];
										ld1 = cache_dml1[q_rear_i];
										rd1 = cache_dmr1[q_rear_i];
#ifdef	KSW_ALN
										dm_l1 = cache_kl1[q_rear_i];
										dm_r1 = cache_kr1[q_rear_i];
#endif
										c_m_f = 1;
										break;

									}
								}
							}
							//do not find the seq in cache, exact match or lv and add into cache
							if(c_m_f == 0)
							{
								//exact match

								mis_c_n = 0;
#ifdef	QUAL_FILT_SINGLE_END
								mis_c_n_filt = 0;
#endif
								ref_tmp_ori = ref_seq_tmp1[tid][ref_copy_num - 2];
								ref_seq_tmp1[tid][ref_copy_num - 2] &= low_mask[tid];

								for(rst_i = 1, read_b_i = 0; rst_i < ref_copy_num - 1; rst_i++, read_b_i++)
								{
									xor_tmp = ref_seq_tmp1[tid][rst_i] ^ read_bit_1[tid][read_b_i];
#ifdef	QUAL_FILT_SINGLE_END
									mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
									xor_tmp &= qual_filt_1[read_b_i];
									mis_c_n_filt += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
									if(mis_c_n_filt > max_mismatch[tid])	break;
#else
									mis_c_n += popcount_3((((xor_tmp & low_bit_mask) << 1) | (xor_tmp & high_bit_mask)));
									mis_c_n_filt = mis_c_n;

									if(mis_c_n > max_mismatch[tid])	break;
#endif
								}

								ref_seq_tmp1[tid][ref_copy_num - 2] = ref_tmp_ori;

								//lv
								if(mis_c_n_filt > max_mismatch[tid])
								{
									if((cir_n == cir_fix_n) && (local_ksw))
									{
#ifdef	KSW_ALN
										for(bit_char_i = s_r_o_l1, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_l1, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										ksw_extend(s_r_o_l1 + 1, read_char[tid], 33 + s_r_o_l1, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r1 - s_r_o_l1, &dm_l1, &tle, &gtle, &gscore, &max_off);

										for(bit_char_i = s_r_o_r1, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_r1, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										ksw_extend(read_length - s_r_o_r1, read_char[tid], 32 + read_length - s_r_o_r1, ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, end_bonus, zdrop, s_r_o_r1 - s_r_o_l1, &dm_r1, &tle, &gtle, &gscore, &max_off);

										dm1 = MAX_OP_SCORE - (dm_l1 + dm_r1);

										ld1 = s_r_o_l1;
										rd1 = s_r_o_r1;
#endif


									}
									else
									{
#ifdef SPLIT_LV_SINGLE


#ifdef	QUAL_FILT_LV_MIS_SINGLE_END
										for(bit_char_i = s_r_o_l1, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_l1, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_l1 = computeEditDistance_mis(ali_ref_seq[tid], 33 + s_r_o_l1, read_char[tid], s_r_o_l1 + 1, lv_dmt1, L[tid], qual_filt_lv_1_o + read_length_a - s_r_o_l1);

										for(bit_char_i = s_r_o_r1, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_r1, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_r1 = computeEditDistance_mis(ali_ref_seq[tid], 32 + read_length - s_r_o_r1, read_char[tid], read_length - s_r_o_r1, lv_dmt1, L[tid], qual_filt_lv_1 + s_r_o_r1);

										dm1 = dm_l1 + dm_r1;

#else
										//split lv
										for(bit_char_i = s_r_o_l1, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_l1, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_l1 = computeEditDistance(ali_ref_seq[tid], 33 + s_r_o_l1, read_char[tid], s_r_o_l1 + 1, lv_dmt1, L[tid]);

										for(bit_char_i = s_r_o_r1, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
											read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										for(bit_char_i = 32 + s_r_o_r1, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
											ali_ref_seq[tid][read_b_i] = ((ref_seq_tmp1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

										dm_r1 = computeEditDistance(ali_ref_seq[tid], 32 + read_length - s_r_o_r1, read_char[tid], read_length - s_r_o_r1, lv_dmt1, L[tid]);

										dm1 = dm_l1 + dm_r1;
#endif
										ld1 = s_r_o_l1;
										rd1 = s_r_o_r1;
#endif
									}
								}
								else
								{
									dm1 = mis_c_n;
									ld1 = 0;
									rd1 = 0;
									dm_l1 = 0;
									dm_r1 = 0;
								}

								//these two values could be different
								if((dm_l1 != -1) && (dm_r1 != -1))
								{
									if(dm1 < dmt1)	dmt1 = dm1 + 1;

									if(dm1 < lv_dmt1)	lv_dmt1 = dm1 + 1;

									if(dm1 < max_mismatch_p[tid])
									{
										dmt1 = 0;
									}
								}
								else
								{
									dm1 = MAX_EDIT_SCORE;
								}

								//add the ref sequence at the end of queue
								memcpy(cache_end1[q_rear1], ref_seq_tmp1[tid], ref_copy_num_chars);
								cache_dis1[q_rear1] = dm1;
								cache_dml1[q_rear1] = ld1;
								cache_dmr1[q_rear1] = rd1;
#ifdef	KSW_ALN
								cache_kl1[q_rear1] = dm_l1;
								cache_kr1[q_rear1] = dm_r1;
#endif
								q_rear1 = ((q_rear1 + 1) & 0X1f);
								++q_n1;

								//add edit distance
							}
						}

						if(dm1 < dm_cir_min)	dm_cir_min = dm1;

						if(dm1 < dm_op[tid])
						{

#ifdef	DM_COPY_SINGLE
							for(dm_i = 0; dm_i < v_cnt; dm_i++)
							{
								ops_vector_pos1[tid][dm_i] = op_vector_pos1[tid][dm_i];

								ops_dm_l1[tid][dm_i] = op_dm_l1[tid][dm_i];
								ops_dm_r1[tid][dm_i] = op_dm_r1[tid][dm_i];

								ops_dm_kl1[tid][dm_i] = op_dm_kl1[tid][dm_i];
								ops_dm_kr1[tid][dm_i] = op_dm_kr1[tid][dm_i];

								ops_rc[tid][dm_i] = op_rc[tid][dm_i];
#ifdef	MAPPING_QUALITY_SINGLE_END
								memcpy(ops_vector_seq1[tid][dm_i], op_vector_seq1[tid][dm_i], ref_copy_num_chars);
#else
								if(!((op_dm_l1[tid][dm_i] == 0) && (op_dm_r1[tid][dm_i] == 0)))
									memcpy(ops_vector_seq1[tid][dm_i], op_vector_seq1[tid][dm_i], ref_copy_num_chars);
#endif
							}
							vs_cnt = v_cnt;
							dm_ops[tid] = dm_op[tid];
#endif

							v_cnt = 0;
							op_vector_pos1[tid][v_cnt] = posi;
#ifdef	MAPPING_QUALITY_SINGLE_END
							memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);
#else
							if(!((ld1 == 0) && (rd1 == 0)))
								memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);
#endif
							op_dm_l1[tid][v_cnt] = ld1;
							op_dm_r1[tid][v_cnt] = rd1;
#ifdef	KSW_ALN
							op_dm_kl1[tid][v_cnt] = dm_l1;
							op_dm_kr1[tid][v_cnt] = dm_r1;
#endif

#ifdef	ALTER_DEBUG_SINGLE_END
							seed_length_arr[tid][v_cnt].seed_length = seed_length1;
							seed_length_arr[tid][v_cnt].index = v_cnt;
#endif
							op_rc[tid][v_cnt] = ((rc_i << 1) + rc_i);
							++v_cnt;
							dm_op[tid] = dm1;
						}
						else if(dm1 == dm_op[tid])
						{
							if(v_cnt < cus_max_output_ali)
							{
								op_vector_pos1[tid][v_cnt] = posi;
#ifdef	MAPPING_QUALITY_SINGLE_END
								memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);
#else
								if(!((ld1 == 0) && (rd1 == 0)))
									memcpy(op_vector_seq1[tid][v_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);
#endif
								op_dm_l1[tid][v_cnt] = ld1;
								op_dm_r1[tid][v_cnt] = rd1;
#ifdef	KSW_ALN
								op_dm_kl1[tid][v_cnt] = dm_l1;
								op_dm_kr1[tid][v_cnt] = dm_r1;
#endif

#ifdef	ALTER_DEBUG_SINGLE_END
								seed_length_arr[tid][v_cnt].seed_length = seed_length1;
								seed_length_arr[tid][v_cnt].index = v_cnt;
#endif
								op_rc[tid][v_cnt] = ((rc_i << 1) + rc_i);
								++v_cnt;
							}
						}
						else if(dm1 < dm_ops[tid])
						{
							vs_cnt = 0;

							ops_vector_pos1[tid][vs_cnt] = posi;
#ifdef	MAPPING_QUALITY_SINGLE_END
							memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);
#else
							if(!((ld1 == 0) && (rd1 == 0)))
								memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);
#endif
							ops_dm_l1[tid][vs_cnt] = ld1;
							ops_dm_r1[tid][vs_cnt] = rd1;
#ifdef	KSW_ALN
							ops_dm_kl1[tid][vs_cnt] = dm_l1;
							ops_dm_kr1[tid][vs_cnt] = dm_r1;
#endif
							ops_rc[tid][vs_cnt] = ((rc_i << 1) + rc_i);

							++vs_cnt;
							dm_ops[tid] = dm1;
						}
						else if(dm1 == dm_ops[tid])
						{
							if(vs_cnt < cus_max_output_ali)
							{
								ops_vector_pos1[tid][vs_cnt] = posi;
#ifdef	MAPPING_QUALITY_SINGLE_END
								memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);
#else
								if(!((ld1 == 0) && (rd1 == 0)))
									memcpy(ops_vector_seq1[tid][vs_cnt], ref_seq_tmp1[tid], ref_copy_num_chars);
#endif
								ops_dm_l1[tid][vs_cnt] = ld1;
								ops_dm_r1[tid][vs_cnt] = rd1;
#ifdef	KSW_ALN
								ops_dm_kl1[tid][vs_cnt] = dm_l1;
								ops_dm_kr1[tid][vs_cnt] = dm_r1;
#endif
								ops_rc[tid][vs_cnt] = ((rc_i << 1) + rc_i);

								++vs_cnt;
							}
						}

					}
				}

#ifdef	CIR_CONTINUE_SINGLE_END
				if(local_ksw)
				{
					if((dm_cir_min > max_single_score) && (cir_n != cir_fix_n))
					{
						seed_l[tid] -= seed_step;
						cir_n = cir_fix_n;
						continue;
					}
				}
				else
				{
					if(dm_cir_min > max_single_score)
					{
						seed_l[tid] -= seed_step;
						cir_n++;
						if(last_circle_rate)
						{
							lv_k = (read_length * last_circle_rate);
							max_single_score = lv_k;
						}
						continue;
					}
				}
#else
				if((dm_op[tid] > max_single_score) && (cir_n < cir_fix_n))
				{
					seed_l[tid] -= seed_step;
					cir_n++;
					continue;
				}
#endif
				break;
			}

			cir_n++;
		}

		//output
		if(v_cnt > 0)
		{
#ifdef	ALTER_DEBUG_SINGLE_END
			if(v_cnt > 1)	qsort(seed_length_arr[tid], v_cnt, sizeof(seed_length_array), compare_seed_length);
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
			if(v_cnt > 1)
			{
				sam_qual1 = 0;
				mp_flag = 0;
			}
			else if(vs_cnt == 0)
			{
				sam_qual1 = 60;
				mp_flag = 0;
			}
			else
			{
				mp_flag = 1;
			}
#endif
			v_cnt_i = 0;

#ifdef	ALTER_DEBUG_SINGLE_END
			if(v_cnt > 1)	v_cnt_i = seed_length_arr[tid][v_cnt_i].index;
#endif
			x = op_vector_pos1[tid][v_cnt_i];
			low = 0;
			high = chr_file_n - 1;

			while ( low <= high )
			{
				mid = (low + high) >> 1;
				if(x < (chr_end_n[mid]))
				{
					high = mid - 1;
				}
				else if(x > (chr_end_n[mid]))
				{
					low = mid + 1;
				}
				else
				{
					chr_re =  mid;
					break;
				}
				chr_re = low;
			}

			sam_pos1 = op_vector_pos1[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;

			if(op_rc[tid][v_cnt_i] == 0)
			{
#ifdef	CHAR_CP_SINGLE_END
				read_bit_1[tid] = read_bit1[tid][0];
#else
				strcpy(sam_seq1, seqio[seqi].read_seq1);
#endif
				sam_flag1 = 0;

#ifdef	QUAL_FILT_LV_OUT_SINGLE_END
				qual_filt_lv_1 = qual_filt_lv1[tid][0];
				qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
				if(mp_flag)
				{
					mp_subs_1 = mp_subs1[tid][0];
					mp_subs_1_o = mp_subs1[tid][1];
				}
#endif
			}
			else
			{
#ifdef	CHAR_CP_SINGLE_END
				read_bit_1[tid] = read_bit1[tid][1];
#else
				for(sam_seq_i = 0; sam_seq_i < read_length; sam_seq_i++)
					sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

				sam_seq1[sam_seq_i] = '\0';

				strrev1(sam_seq1);
#endif
				sam_flag1 = 16;

#ifdef	QUAL_FILT_LV_OUT_SINGLE_END
				qual_filt_lv_1 = qual_filt_lv1[tid][1];
				qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
				if(mp_flag)
				{
					mp_subs_1 = mp_subs1[tid][1];
					mp_subs_1_o = mp_subs1[tid][0];
				}
#endif
			}
#ifdef	MAPPING_QUALITY_SINGLE_END
			if(mp_flag)	sub_t[tid] = 0;
#endif
			d_n1 = 0;
			i_n1 = 0;
			s_offset1 = 0;
			s_r_o_l = op_dm_l1[tid][v_cnt_i];
			s_r_o_r = op_dm_r1[tid][v_cnt_i];
			if((s_r_o_l == 0) && (s_r_o_r == 0))
			{
				strcpy(cigar_p1, cigar_m[tid]);
#ifdef	MAPPING_QUALITY_SINGLE_END

				if(mp_flag)
				{
					for (bit_char_i = 0, read_b_i = 32; bit_char_i < read_length; bit_char_i++, read_b_i++)
						if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((op_vector_seq1[tid][v_cnt_i][read_b_i >> 5] >> ((31 - (read_b_i & 0X1f)) << 1)) & 0X3))
						{
							sub_t[tid] += mp_subs_1[bit_char_i];
						}
				}
#endif
			}
			else     //indel
			{
				if((cir_n == cir_fix_n) && (local_ksw))
				{
#ifdef	KSW_ALN

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
					for(bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
						read_char[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif
					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i > -1
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					ksw_global(op_dm_kl1[tid][v_cnt_i], read_char[tid], op_dm_kl1[tid][v_cnt_i], ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar1, &cigar1);

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
						read_char2[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
					for(bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
						read_char2[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif
					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length + 64
						ali_ref_seq2[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					ksw_global(op_dm_kr1[tid][v_cnt_i], read_char2[tid], op_dm_kr1[tid][v_cnt_i], ali_ref_seq2[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar2, &cigar2);

					m_m_n = s_r_o_r - s_r_o_l - 1;

					nm_score = 0;
					snt = 0;
					if(n_cigar1)
					{
						op_score1  = cigar1[0]&0xf;
						len_score1 = cigar1[0]>>4;
					}
					else	op_score1 = 3;

					if(n_cigar2)
					{
						op_score2  = cigar2[0]&0xf;
						len_score2 = cigar2[0]>>4;
					}
					else	op_score2 = 3;


					if(s_r_o_l >= op_dm_kl1[tid][v_cnt_i])
					{
						sn = sprintf(cigar_p1 + snt, "%dS", s_r_o_l + 1 - op_dm_kl1[tid][v_cnt_i]);
						snt += sn;
					}

					if((s_r_o_l != -1) && (s_r_o_r != read_length))
					{
						if((op_score1 == 0) && (op_score2 == 0))
						{
							k_start1 = 0;
							k_start2 = 1;
							k_middle = len_score1 + len_score2 + m_m_n;
						}
						else if(op_score1 == 0)
						{
							k_start1 = 0;
							k_start2 = 0;
							k_middle = len_score1 + m_m_n;
						}
						else if(op_score2 == 0)
						{
							k_start1 = -1;
							k_start2 = 1;
							k_middle = len_score2 + m_m_n;
						}
						else
						{
							k_start1 = -1;
							k_start2 = 0;
							k_middle = m_m_n;
						}

						x_score = y_score = 0;
						for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
						{
							op_score  = cigar1[bit_char_i]&0xf;
							len_score = cigar1[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
									{
#ifdef	MAPPING_QUALITY_SINGLE_END
										if(mp_flag)	sub_t[tid] += mp_subs_1_o[read_length + n_cigar1 - s_r_o_l - x_score - read_b_i - 2];

#endif
										++nm_score;
									}
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
						}

						sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
						snt += sn;

						x_score = y_score = 0;
						for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
						{
							op_score  = cigar2[bit_char_i]&0xf;
							len_score = cigar2[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
									{
#ifdef	MAPPING_QUALITY_SINGLE_END
										if(mp_flag)	sub_t[tid] += mp_subs_1[s_r_o_r + x_score + read_b_i];
#endif
										++nm_score;
									}
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score;
						}
					}
					else if(s_r_o_l == -1)
					{
						if(op_score2 == 0)
						{
							k_start2 = 1;
							k_middle = len_score2 + m_m_n;
						}
						else
						{
							k_start2 = 0;
							k_middle = m_m_n;
						}
						sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
						snt += sn;

						x_score = y_score = 0;
						for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
						{
							op_score  = cigar2[bit_char_i]&0xf;
							len_score = cigar2[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
									{
#ifdef	MAPPING_QUALITY_SINGLE_END
										if(mp_flag)	sub_t[tid] += mp_subs_1[s_r_o_r + x_score + read_b_i];
#endif
										++nm_score;
									}
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score;
						}
					}
					else
					{
						if(op_score1 == 0)
						{
							k_start1 = 0;
							k_middle = len_score1 + m_m_n;
						}
						else
						{
							k_start1 = -1;
							k_middle = m_m_n;
						}
						x_score = y_score = 0;
						for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
						{
							op_score  = cigar1[bit_char_i]&0xf;
							len_score = cigar1[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
									{
#ifdef	MAPPING_QUALITY_SINGLE_END
										if(mp_flag)	sub_t[tid] += mp_subs_1_o[read_length + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
										++nm_score;
									}
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
						}

						sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
						snt += sn;
					}

					if(read_length - s_r_o_r > op_dm_kr1[tid][v_cnt_i])
					{
						sn = sprintf(cigar_p1 + snt, "%dS", read_length -s_r_o_r - op_dm_kr1[tid][v_cnt_i]);
						snt += sn;
					}
					//sn = sprintf(cigar_p1 + snt, "\0");
					//snt += sn;

					if(n_cigar1)	free(cigar1);
					if(n_cigar2)	free(cigar2);
#endif

				}
				else
				{
					nm_score = dm_op[tid];

#ifdef	QUAL_FILT_LV_OUT_SINGLE_END

#ifdef	CHAR_CP_SINGLE_END
					for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
					if(mp_flag)	computeEditDistanceWithCigar_s_mis_left_mp(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length - 1 - s_r_o_l, &s_offset1, mp_subs_1_o + read_length - 1 - s_r_o_l, &(sub_t[tid]));
					else	computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length - 1 - s_r_o_l, &s_offset1);
#else
					computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length - 1 - s_r_o_l, &s_offset1);
#endif

#ifdef	CHAR_CP_SINGLE_END
					for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
					if(mp_flag)	computeEditDistanceWithCigar_s_mis_mp(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - op_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + s_r_o_r, mp_subs_1 + s_r_o_r, &(sub_t[tid]));
					else	computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - op_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + s_r_o_r);
#else
					computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - op_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + s_r_o_r);
#endif

#else

#ifdef	CHAR_CP_SINGLE_END
					for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
					for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
					if(mp_flag)	computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid], mp_subs_1_o + read_length - 1 - s_r_o_l, &(sub_t[tid]));
					else	computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid]);
#else
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid]);
#endif

#ifdef	CHAR_CP_SINGLE_END
					for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
					if(mp_flag)	computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - op_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid], mp_subs_1 + s_r_o_r, &(sub_t[tid]));
					else	computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - op_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid]);
#else
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - op_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid]);
#endif

#endif
					//deal with front and back lv cigar
					m_m_n = s_r_o_r - s_r_o_l - 1;

					strncpy(str_o, cigarBuf1, f_cigarn);
					s_o = 0;
					f_c = 0;

					pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

					while (pch != NULL)
					{
						pchl = strlen(pch);
						f_cigar[f_cigarn - f_c - 2] = atoi(pch);
						s_o += (pchl + 1);
						f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];
						f_c += 2;

						if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
						if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

						pch = strtok_r(NULL, "DMIS", &saveptr);
					}

					strncpy(b_cigar, cigarBuf2, f_cigarn);
					pch = strtok (cigarBuf2,"DMIS");

					if(pch != NULL)
						pchl = strlen(pch);

					snt = 0;
					if((s_r_o_l != -1) && (s_r_o_r != read_length))
					{
						if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
						{
							f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s", b_cigar + pchl + 1);
							snt += sn;
						}
						else if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;

							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s",b_cigar);
							snt += sn;
						}
						else if(b_cigar[pchl] == 'M')
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							snt += sn;
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
							snt += sn;
						}

					}
					else if(s_r_o_l == -1)
					{
						if(b_cigar[pchl] == 'M')
							sn = sprintf(cigar_p1, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						else	sn = sprintf(cigar_p1, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
					else
					{
						if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
							snt += sn;
						}
					}
					//sn = sprintf(cigar_p1 + snt, "\0");
					//snt += sn;
				}

#ifdef	CIGAR_LEN_ERR
				cigar_len = 0;
				s_o_tmp = 0;
				cigar_len_tmp = 0;
				strncpy(cigar_tmp, cigar_p1, snt);
				cigar_tmp[snt] = '\0';
				pch_tmp = strtok_r(cigar_tmp,"DMIS", &saveptr_tmp);

				while (pch_tmp != NULL)
				{
					pchl_tmp = strlen(pch_tmp);
					s_o_tmp += (pchl_tmp + 1);

					if(cigar_p1[s_o_tmp - 1] != 'D')
					{
						cigar_len_tmp = atoi(pch_tmp);
						cigar_len += cigar_len_tmp;
					}

					pch_tmp = strtok_r(NULL, "DMIS", &saveptr_tmp);
				}

				if(read_length != cigar_len)
				{
					if(read_length < cigar_len)
					{
						cigar_len_re = cigar_len_tmp - (cigar_len - read_length);
						if(cigar_len_re > 0)	sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
						else if(cigar_len_re == 0)	sprintf(cigar_p1 + snt - sn, "\0");
						else	strcpy(cigar_p1, cigar_m[tid]);
					}
					else
					{
						cigar_len_re = cigar_len_tmp + (read_length - cigar_len);
						sprintf(cigar_p1 + snt - sn, "%u%c", cigar_len_re, cigar_p1[snt - 1]);
					}
				}
#endif

#ifdef	NO_S_OFF
				s_offset1 = 0;
#endif
				sam_pos1 = sam_pos1 + i_n1 - d_n1 + s_offset1;
			}

#ifdef	MAPPING_QUALITY_SINGLE_END
			if(mp_flag)
			{
				sam_qual1 = sub_t[tid];
			}
#endif

#ifdef	OUTPUT_ARR
			seqio[seqi].flag1 = sam_flag1;
			seqio[seqi].chr_re = chr_re;
			seqio[seqi].pos1 = sam_pos1;

			strcpy(pr_cigar1_buffer[seqi], cigar_p1);
			seqio[seqi].cigar1 = pr_cigar1_buffer[seqi];

			seqio[seqi].nm1 = nm_score;

			if(sam_flag1 == 0)
			{
				seqio[seqi].seq1 = seqio[seqi].read_seq1;
			}
			else
			{
#ifdef	CHAR_CP_SINGLE_END
				for(sam_seq_i = 0; sam_seq_i < read_length; sam_seq_i++)
					sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

				sam_seq1[sam_seq_i] = '\0';
				strrev1(sam_seq1);
#endif
				strcpy(read_rev_buffer[seqi], sam_seq1);
				read_rev_buffer[seqi][read_length] = '\0';
				seqio[seqi].seq1 = read_rev_buffer[seqi];
				strrev1(qual1_buffer[seqi]);
			}
#endif
		}

		xa_i = 0;
		for(va_cnt_i = 1; va_cnt_i < v_cnt; va_cnt_i++)
		{
#ifdef	ALTER_DEBUG_SINGLE_END
			v_cnt_i = seed_length_arr[tid][va_cnt_i].index;
#else
			v_cnt_i = va_cnt_i;
#endif
			x = op_vector_pos1[tid][v_cnt_i];
			low = 0;
			high = chr_file_n - 1;

			while ( low <= high )
			{
				mid = (low + high) >> 1;
				if(x < (chr_end_n[mid]))
				{
					high = mid - 1;
				}
				else if(x > (chr_end_n[mid]))
				{
					low = mid + 1;
				}
				else
				{
					chr_re =  mid;
					break;
				}
				chr_re = low;
			}

			sam_pos1 = op_vector_pos1[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;

			chr_res[tid][xa_i] = chr_re;

			if(op_rc[tid][v_cnt_i] == 0)
			{
#ifdef	CHAR_CP_SINGLE_END
				read_bit_1[tid] = read_bit1[tid][0];
#else
				strcpy(sam_seq1, seqio[seqi].read_seq1);
#endif

#ifdef	QUAL_FILT_LV_OUT_SINGLE_END
				qual_filt_lv_1 = qual_filt_lv1[tid][0];
				qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif
			}
			else
			{
#ifdef	CHAR_CP_SINGLE_END
				read_bit_1[tid] = read_bit1[tid][1];
#else
				for(sam_seq_i = 0; sam_seq_i < read_length; sam_seq_i++)
					sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

				sam_seq1[sam_seq_i] = '\0';

				strrev1(sam_seq1);
#endif

#ifdef	QUAL_FILT_LV_OUT_SINGLE_END
				qual_filt_lv_1 = qual_filt_lv1[tid][1];
				qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif

			}

			d_n1 = 0;
			i_n1 = 0;
			lv_re1f = dm_op[tid];
			s_offset1 = 0;
			s_r_o_l = op_dm_l1[tid][v_cnt_i];
			s_r_o_r = op_dm_r1[tid][v_cnt_i];
			if((s_r_o_l == 0) && (s_r_o_r == 0))
			{
				strcpy(cigar_p1, cigar_m[tid]);
			}
			else     //indel
			{
				if((cir_n == cir_fix_n) && (local_ksw))
				{
#ifdef	KSW_ALN

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
					for(bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
						read_char[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif
					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; read_b_i < op_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i > -1
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					ksw_global(op_dm_kl1[tid][v_cnt_i], read_char[tid], op_dm_kl1[tid][v_cnt_i], ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar1, &cigar1);

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
						read_char2[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
					for(bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
						read_char2[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif
					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length + 64
						ali_ref_seq2[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					ksw_global(op_dm_kr1[tid][v_cnt_i], read_char2[tid], op_dm_kr1[tid][v_cnt_i], ali_ref_seq2[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar2, &cigar2);

					m_m_n = s_r_o_r - s_r_o_l - 1;

					nm_score = 0;
					snt = 0;
					if(n_cigar1)
					{
						op_score1  = cigar1[0]&0xf;
						len_score1 = cigar1[0]>>4;
					}
					else	op_score1 = 3;

					if(n_cigar2)
					{
						op_score2  = cigar2[0]&0xf;
						len_score2 = cigar2[0]>>4;
					}
					else	op_score2 = 3;

					if(s_r_o_l >= op_dm_kl1[tid][v_cnt_i])
					{
						sn = sprintf(cigar_p1 + snt, "%dS", s_r_o_l + 1 - op_dm_kl1[tid][v_cnt_i]);
						snt += sn;
					}

					if((s_r_o_l != -1) && (s_r_o_r != read_length))
					{
						if((op_score1 == 0) && (op_score2 == 0))
						{
							k_start1 = 0;
							k_start2 = 1;
							k_middle = len_score1 + len_score2 + m_m_n;
						}
						else if(op_score1 == 0)
						{
							k_start1 = 0;
							k_start2 = 0;
							k_middle = len_score1 + m_m_n;
						}
						else if(op_score2 == 0)
						{
							k_start1 = -1;
							k_start2 = 1;
							k_middle = len_score2 + m_m_n;
						}
						else
						{
							k_start1 = -1;
							k_start2 = 0;
							k_middle = m_m_n;
						}

						x_score = y_score = 0;
						for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
						{
							op_score  = cigar1[bit_char_i]&0xf;
							len_score = cigar1[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i]) ++nm_score;
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
						}

						sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
						snt += sn;

						x_score = y_score = 0;
						for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
						{
							op_score  = cigar2[bit_char_i]&0xf;
							len_score = cigar2[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i]) ++nm_score;
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score;
						}

					}
					else if(s_r_o_l == -1)
					{
						if(op_score2 == 0)
						{
							k_start2 = 1;
							k_middle = len_score2 + m_m_n;
						}
						else
						{
							k_start2 = 0;
							k_middle = m_m_n;
						}
						sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
						snt += sn;

						x_score = y_score = 0;
						for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
						{
							op_score  = cigar2[bit_char_i]&0xf;
							len_score = cigar2[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i]) ++nm_score;
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score;
						}
					}
					else
					{
						if(op_score1 == 0)
						{
							k_start1 = 0;
							k_middle = len_score1 + m_m_n;
						}
						else
						{
							k_start1 = -1;
							k_middle = m_m_n;
						}
						x_score = y_score = 0;
						for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
						{
							op_score  = cigar1[bit_char_i]&0xf;
							len_score = cigar1[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i]) ++nm_score;
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
						}

						sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
						snt += sn;
					}

					if(read_length - s_r_o_r > op_dm_kr1[tid][v_cnt_i])
					{
						sn = sprintf(cigar_p1 + snt, "%dS", read_length - s_r_o_r - op_dm_kr1[tid][v_cnt_i]);
						snt += sn;
					}
					//sn = sprintf(cigar_p1 + snt, "\0");
					//snt += sn;

					lv_re1f = nm_score;

					if(n_cigar1)	free(cigar1);
					if(n_cigar2)	free(cigar2);
#endif
				}
				else
				{
#ifdef	QUAL_FILT_LV_OUT_SINGLE_END

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length - 1 - s_r_o_l, &s_offset1);

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
					for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - op_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + s_r_o_r);

#else

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid]);

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((op_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - op_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid]);
#endif
					//deal with front and back lv cigar
					m_m_n = s_r_o_r - s_r_o_l - 1;

					strncpy(str_o, cigarBuf1, f_cigarn);
					s_o = 0;
					f_c = 0;

					pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

					while (pch != NULL)
					{
						pchl = strlen(pch);
						f_cigar[f_cigarn - f_c - 2] = atoi(pch);
						s_o += (pchl + 1);
						f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];
						f_c += 2;

						if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
						if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

						pch = strtok_r(NULL, "DMIS", &saveptr);
					}

					strncpy(b_cigar, cigarBuf2, f_cigarn);

					pch = strtok (cigarBuf2,"DMIS");

					if(pch != NULL)
						pchl = strlen(pch);

					snt = 0;
					if((s_r_o_l != -1) && (s_r_o_r != read_length))
					{
						if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
						{
							f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s", b_cigar + pchl + 1);
							snt += sn;
						}
						else if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s",b_cigar);
							snt += sn;
						}
						else if(b_cigar[pchl] == 'M')
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							snt += sn;
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
							snt += sn;
						}
					}
					else if(s_r_o_l == -1)
					{
						if(b_cigar[pchl] == 'M')
							sn = sprintf(cigar_p1, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						else	sn = sprintf(cigar_p1, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
					else
					{
						if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
							snt += sn;
						}
					}
					//sn = sprintf(cigar_p1 + snt, "\0");
					//snt += sn;
				}
				sam_pos1 = sam_pos1 + i_n1 - d_n1;
			}

			if(op_rc[tid][v_cnt_i] == 0)
			{
				xa_d1s[tid][xa_i] = '+';
			}
			else
			{
				xa_d1s[tid][xa_i] = '-';
			}
#ifdef	NO_S_OFF
			s_offset1 = 0;
#endif
			sam_pos1s[tid][xa_i] = (uint32_t )sam_pos1 + s_offset1;

			strcpy(cigar_p1s[tid][xa_i], cigar_p1);

			lv_re1s[tid][xa_i] = lv_re1f;
			xa_i++;
		}

		lv_re1f = dm_ops[tid];

#ifdef	ALT_ALL
		for(v_cnt_i = 0; (v_cnt_i < vs_cnt); v_cnt_i++)
#else
		for(v_cnt_i = 0; (v_cnt_i < vs_cnt) && (v_cnt < 2) && (dm_op[tid] + 3 > dm_ops[tid]); v_cnt_i++)
#endif
		{
			x = ops_vector_pos1[tid][v_cnt_i];
			low = 0;
			high = chr_file_n - 1;

			while ( low <= high )
			{
				mid = (low + high) >> 1;
				if(x < (chr_end_n[mid]))
				{
					high = mid - 1;
				}
				else if(x > (chr_end_n[mid]))
				{
					low = mid + 1;
				}
				else
				{
					chr_re =  mid;
					break;
				}
				chr_re = low;
			}

			sam_pos1 = ops_vector_pos1[tid][v_cnt_i] - chr_end_n[chr_re - 1] + 1;

			chr_res[tid][xa_i] = chr_re;

			if(ops_rc[tid][v_cnt_i] == 0)
			{
#ifdef	CHAR_CP_SINGLE_END
				read_bit_1[tid] = read_bit1[tid][0];
#else
				strcpy(sam_seq1, seqio[seqi].read_seq1);
#endif

#ifdef	QUAL_FILT_LV_OUT_SINGLE_END
				qual_filt_lv_1 = qual_filt_lv1[tid][0];
				qual_filt_lv_1_o = qual_filt_lv1[tid][1];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
				if((mp_flag) && (v_cnt_i == 0))
				{
					mp_subs_1 = mp_subs1[tid][0];
					mp_subs_1_o = mp_subs1[tid][1];
				}
#endif
			}
			else
			{
#ifdef	CHAR_CP_SINGLE_END
				read_bit_1[tid] = read_bit1[tid][1];
#else
				for(sam_seq_i = 0; sam_seq_i < read_length; sam_seq_i++)
					sam_seq1[sam_seq_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[sam_seq_i]] ^ 0X3];

				sam_seq1[sam_seq_i] = '\0';
				strrev1(sam_seq1);
#endif

#ifdef	QUAL_FILT_LV_OUT_SINGLE_END
				qual_filt_lv_1 = qual_filt_lv1[tid][1];
				qual_filt_lv_1_o = qual_filt_lv1[tid][0];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
				if((mp_flag) && (v_cnt_i == 0))
				{
					mp_subs_1 = mp_subs1[tid][1];
					mp_subs_1_o = mp_subs1[tid][0];
				}
#endif
			}

			d_n1 = 0;
			i_n1 = 0;
			s_offset1 = 0;
			s_r_o_l = ops_dm_l1[tid][v_cnt_i];
			s_r_o_r = ops_dm_r1[tid][v_cnt_i];

#ifdef	MAPPING_QUALITY_SINGLE_END
			if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] = 0;
#endif
			if(s_r_o_l && (s_r_o_r == 0))
			{
				strcpy(cigar_p1, cigar_m[tid]);

#ifdef	MAPPING_QUALITY_SINGLE_END

				if((mp_flag) && (v_cnt_i == 0))
				{
					for (bit_char_i = 0, read_b_i = 32; bit_char_i < read_length; bit_char_i++, read_b_i++)
						if(((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3) != ((ops_vector_seq1[tid][v_cnt_i][read_b_i >> 5] >> ((31 - (read_b_i & 0X1f)) << 1)) & 0X3))
						{
							sub_t[tid] += mp_subs_1[bit_char_i];
						}
				}
#endif

			}
			else     //indel
			{
				if((cir_n == cir_fix_n) && (local_ksw))
				{
#ifdef	KSW_ALN

#ifdef CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < ops_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
					for(bit_char_i = s_r_o_l, read_b_i = 0; read_b_i < ops_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i >= 0
						read_char[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif
					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; read_b_i < ops_dm_kl1[tid][v_cnt_i]; bit_char_i--, read_b_i++)//bit_char_i > -1
						ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					ksw_global(ops_dm_kl1[tid][v_cnt_i], read_char[tid], ops_dm_kl1[tid][v_cnt_i], ali_ref_seq[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar1, &cigar1);

#ifdef CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < ops_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
						read_char2[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
#else
					for(bit_char_i = s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length
						read_char2[tid][read_b_i] = charToDna5n[sam_seq1[bit_char_i]];
#endif
					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; read_b_i < op_dm_kr1[tid][v_cnt_i]; bit_char_i++, read_b_i++)//bit_char_i < read_length + 64
						ali_ref_seq2[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					ksw_global(op_dm_kr1[tid][v_cnt_i], read_char2[tid], op_dm_kr1[tid][v_cnt_i], ali_ref_seq2[tid], 5, mat, gapo_score, gape_score, band_with, &n_cigar2, &cigar2);

					m_m_n = s_r_o_r - s_r_o_l - 1;

					nm_score = 0;
					snt = 0;
					if(n_cigar1)
					{
						op_score1  = cigar1[0]&0xf;
						len_score1 = cigar1[0]>>4;
					}
					else	op_score1 = 3;

					if(n_cigar2)
					{
						op_score2  = cigar2[0]&0xf;
						len_score2 = cigar2[0]>>4;
					}
					else	op_score2 = 3;

					if(s_r_o_l >= ops_dm_kl1[tid][v_cnt_i])
					{
						sn = sprintf(cigar_p1 + snt, "%dS", s_r_o_l + 1 - ops_dm_kl1[tid][v_cnt_i]);
						snt += sn;
					}

					if((s_r_o_l != -1) && (s_r_o_r != read_length))
					{
						if((op_score1 == 0) && (op_score2 == 0))
						{
							k_start1 = 0;
							k_start2 = 1;
							k_middle = len_score1 + len_score2 + m_m_n;
						}
						else if(op_score1 == 0)
						{
							k_start1 = 0;
							k_start2 = 0;
							k_middle = len_score1 + m_m_n;
						}
						else if(op_score2 == 0)
						{
							k_start1 = -1;
							k_start2 = 1;
							k_middle = len_score2 + m_m_n;
						}
						else
						{
							k_start1 = -1;
							k_start2 = 0;
							k_middle = m_m_n;
						}

						x_score = y_score = 0;
						for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
						{
							op_score  = cigar1[bit_char_i]&0xf;
							len_score = cigar1[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
									{
#ifdef	MAPPING_QUALITY_SINGLE_END
										if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_1_o[read_length + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
										++nm_score;
									}
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
						}

						sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
						snt += sn;

						x_score = y_score = 0;
						for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
						{
							op_score  = cigar2[bit_char_i]&0xf;
							len_score = cigar2[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
									{
#ifdef	MAPPING_QUALITY_SINGLE_END
										if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_1[s_r_o_r + x_score + read_b_i];
#endif
										++nm_score;
									}
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score;
						}

					}
					else if(s_r_o_l == -1)
					{
						if(op_score2 == 0)
						{
							k_start2 = 1;
							k_middle = len_score2 + m_m_n;
						}
						else
						{
							k_start2 = 0;
							k_middle = m_m_n;
						}
						sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
						snt += sn;

						x_score = y_score = 0;
						for (bit_char_i = k_start2; bit_char_i < n_cigar2; bit_char_i++)
						{
							op_score  = cigar2[bit_char_i]&0xf;
							len_score = cigar2[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char2[tid][x_score + read_b_i] != ali_ref_seq2[tid][y_score + read_b_i])
									{
#ifdef	MAPPING_QUALITY_SINGLE_END
										if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_1[s_r_o_r + x_score + read_b_i];
#endif
										++nm_score;
									}
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score;
						}
					}
					else
					{
						if(op_score1 == 0)
						{
							k_start1 = 0;
							k_middle = len_score1 + m_m_n;
						}
						else
						{
							k_start1 = -1;
							k_middle = m_m_n;
						}
						x_score = y_score = 0;
						for (bit_char_i = n_cigar1 - 1; bit_char_i > k_start1; bit_char_i--)
						{
							op_score  = cigar1[bit_char_i]&0xf;
							len_score = cigar1[bit_char_i]>>4;

							sn = sprintf(cigar_p1 + snt, "%d%c", len_score, ksw_cigars[op_score]);
							snt += sn;

							if (op_score == 0)   // match
							{
								for (read_b_i = 0; read_b_i < len_score; ++read_b_i)
									if (read_char[tid][n_cigar1 - 1 - x_score - read_b_i] != ali_ref_seq[tid][n_cigar1 - 1 - y_score - read_b_i])
									{
#ifdef	MAPPING_QUALITY_SINGLE_END
										if((mp_flag) && (v_cnt_i == 0))	sub_t[tid] += mp_subs_1_o[read_length + n_cigar1 - s_r_o_l - 2 - x_score - read_b_i];
#endif
										++nm_score;
									}
								x_score += len_score;
								y_score += len_score;
							}
							else if (op_score == 1) x_score += len_score, nm_score += len_score, i_n1 += len_score;
							else if (op_score == 2) y_score += len_score, nm_score += len_score, d_n1 += len_score;
						}

						sn = sprintf(cigar_p1 + snt, "%dM", k_middle);
						snt += sn;
					}

					if(read_length - s_r_o_r > ops_dm_kr1[tid][v_cnt_i])
					{
						sn = sprintf(cigar_p1 + snt, "%dS", read_length - s_r_o_r - ops_dm_kr1[tid][v_cnt_i]);
						snt += sn;
					}
					//sn = sprintf(cigar_p1 + snt, "\0");
					//snt += sn;

					lv_re1f = nm_score;

					if(n_cigar1)	free(cigar1);
					if(n_cigar2)	free(cigar2);
#endif
				}
				else
				{
#ifdef	QUAL_FILT_LV_OUT_SINGLE_END
#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
					if((mp_flag) && (v_cnt_i == 0))
						computeEditDistanceWithCigar_s_mis_left_mp(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length - 1 - s_r_o_l, &s_offset1, mp_subs_1_o + read_length - 1 - s_r_o_l, &(sub_t[tid]));
					else	computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length - 1 - s_r_o_l, &s_offset1);
#else
					computeEditDistanceWithCigar_s_mis_left(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid], qual_filt_lv_1_o + read_length - 1 - s_r_o_l, &s_offset1);
#endif

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
					if((mp_flag) && (v_cnt_i == 0))
						computeEditDistanceWithCigar_s_mis_mp(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - ops_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + s_r_o_r, mp_subs_1 + s_r_o_r, &(sub_t[tid]));
					else	computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - ops_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + s_r_o_r);
#else
					computeEditDistanceWithCigar_s_mis(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - ops_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid], qual_filt_lv_1 + s_r_o_r);
#endif

#else

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for (bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > - 1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_l, read_b_i = 0; bit_char_i > -1; bit_char_i--, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
					if((mp_flag) && (v_cnt_i == 0))
						computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid], mp_subs_1_o + read_length - 1 - s_r_o_l, &(sub_t[tid]));
					else	computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid]);
#else
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 33 + s_r_o_l, read_char[tid], s_r_o_l + 1, lv_k, cigarBuf1, f_cigarn, L[tid]);
#endif

#ifdef	CHAR_CP_SINGLE_END
					for (bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = ((read_bit_1[tid][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for (bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = ((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

#else
					for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_length; bit_char_i++, read_b_i++)
						read_char[tid][read_b_i] = sam_seq1[bit_char_i];

					for(bit_char_i = 32 + s_r_o_r, read_b_i = 0; bit_char_i < read_length + 64; bit_char_i++, read_b_i++)
						ali_ref_seq[tid][read_b_i] = Dna5Tochar[((ops_vector_seq1[tid][v_cnt_i][bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3)];
#endif

#ifdef	MAPPING_QUALITY_SINGLE_END
					if((mp_flag) && (v_cnt_i == 0))
						computeEditDistanceWithCigar_s_mp(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - ops_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid], mp_subs_1 + s_r_o_r, &(sub_t[tid]));
					else	computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - ops_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid]);
#else
					computeEditDistanceWithCigar_s(ali_ref_seq[tid], 32 + read_length - s_r_o_r, read_char[tid], read_length - ops_dm_r1[tid][v_cnt_i], lv_k, cigarBuf2, f_cigarn, L[tid]);
#endif

#endif
					//deal with front and back lv cigar
					m_m_n = s_r_o_r - s_r_o_l - 1;

					strncpy(str_o, cigarBuf1, f_cigarn);
					s_o = 0;
					f_c = 0;

					pch = strtok_r(cigarBuf1,"DMIS", &saveptr);

					while (pch != NULL)
					{
						pchl = strlen(pch);
						f_cigar[f_cigarn - f_c - 2] = atoi(pch);
						s_o += (pchl + 1);
						f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];
						f_c += 2;

						if(str_o[s_o - 1] == 'D')	d_n1 += atoi(pch);
						if(str_o[s_o - 1] == 'I')	i_n1 += atoi(pch);

						pch = strtok_r(NULL, "DMIS", &saveptr);
					}

					strncpy(b_cigar, cigarBuf2, f_cigarn);
					pch = strtok (cigarBuf2,"DMIS");

					if(pch != NULL)
						pchl = strlen(pch);

					snt = 0;
					if((s_r_o_l != -1) && (s_r_o_r != read_length))
					{
						if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
						{
							f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s", b_cigar + pchl + 1);
							snt += sn;
						}
						else if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%s",b_cigar);
							snt += sn;
						}
						else if(b_cigar[pchl] == 'M')
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}

							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
							snt += sn;
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM%s", m_m_n, b_cigar);
							snt += sn;
						}
					}
					else if(s_r_o_l == -1)
					{
						if(b_cigar[pchl] == 'M')
							sn = sprintf(cigar_p1, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
						else	sn = sprintf(cigar_p1, "%uM%s", m_m_n, b_cigar);
						snt += sn;
					}
					else
					{
						if(f_cigar[f_cigarn - 1] == 'M')
						{
							f_cigar[f_cigarn - 2] += m_m_n;
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
						}
						else
						{
							for(f_i = 0; f_i < f_c; f_i += 2)
							{
								sn = sprintf(cigar_p1 + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
								snt += sn;
							}
							sn = sprintf(cigar_p1 + snt, "%uM", m_m_n);
							snt += sn;
						}
					}
					//sn = sprintf(cigar_p1 + snt, "\0");
					//snt += sn;
				}
				sam_pos1 = sam_pos1 + i_n1 - d_n1;
			}

			if(ops_rc[tid][v_cnt_i] == 0)
			{
				xa_d1s[tid][xa_i] = '+';
			}
			else
			{
				xa_d1s[tid][xa_i] = '-';
			}
#ifdef	NO_S_OFF
			s_offset1 = 0;
#endif
			sam_pos1s[tid][xa_i] = (uint32_t )sam_pos1 + s_offset1;

			strcpy(cigar_p1s[tid][xa_i], cigar_p1);

			lv_re1s[tid][xa_i] = lv_re1f;
			xa_i++;
		}

		seqio[seqi].v_cnt = v_cnt;
		if(v_cnt > 0)
		{
#ifdef	MAPPING_QUALITY
			if(mp_flag)
			{
				if(xa_i)
				{
					log_tmp = log(vs_cnt);

					sam_qual1 = (sam_qual1 - sub_t[tid] - log_tmp) * 3.434;
					sam_qual1 = (sam_qual1 > 2 ? sam_qual1 : 2);
#ifdef	QUAL_20
					sam_qual1 = (sam_qual1 < 20 ? 20 : sam_qual1);
#endif
					sam_qual1 = (sam_qual1 > 60 ? 60 : sam_qual1);
				}
				else
				{
					sam_qual1 = 60;
				}

			}
#endif
			seqio[seqi].qualc1 = sam_qual1;

			seqio[seqi].xa_n = xa_i;

			if(xa_i > 0)
			{
				memcpy(chr_res_buffer[seqi], chr_res[tid], xa_i << 2);
				seqio[seqi].chr_res = chr_res_buffer[seqi];

				memcpy(xa_d1s_buffer[seqi], xa_d1s[tid], xa_i);
				seqio[seqi].xa_d1s = xa_d1s_buffer[seqi];

				memcpy(sam_pos1s_buffer[seqi], sam_pos1s[tid], xa_i << 2);
				seqio[seqi].sam_pos1s = sam_pos1s_buffer[seqi];

				memcpy(lv_re1s_buffer[seqi], lv_re1s[tid], xa_i << 2);
				seqio[seqi].lv_re1s = lv_re1s_buffer[seqi];

				for(v_cnt_i = 0; v_cnt_i < xa_i; v_cnt_i++)
				{
					pos_l = strlen(cigar_p1s[tid][v_cnt_i]) + 1;
					if((seqio[seqi].cigar_p1s[v_cnt_i] = (char* )calloc(pos_l, 1)) == NULL)	exit(1);//

					memcpy(seqio[seqi].cigar_p1s[v_cnt_i], cigar_p1s[tid][v_cnt_i], pos_l - 1);
				}
			}
		}
		else
		{

#ifdef	OUTPUT_ARR
			seqio[seqi].xa_n = 0;
			seqio[seqi].flag1 = 4;
			seqio[seqi].chr_re = chr_file_n;
			seqio[seqi].pos1 = 0;
			seqio[seqi].qualc1 = 0;

			strcpy(pr_cigar1_buffer[seqi], "*");
			seqio[seqi].cigar1 = pr_cigar1_buffer[seqi];

			seqio[seqi].seq1 = seqio[seqi].read_seq1;
#endif
		}
	}
}
