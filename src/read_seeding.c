/*************************************************************************
	> File Name: read_seeding.c
	> Author: Yadong Liu
	> Mail: 
 ************************************************************************/

#include <stdio.h>
#include <time.h>
#include <ctype.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>

#include "read_seeding.h"
#include "bit_operation.h"
#include "binarys_qsort.h"
#include "load_unipath_size.h"
#include "graph.h"
#include "bseq.h"
#include "aln_2pass.h"
#include "desalt_index.h"
#include "ktime.h"

//variable extern
uni_seed*** uniseed = NULL;
int8_t* mata_D = NULL;
int8_t* mata_R = NULL;
uint8_t thread_n;
uint8_t k_t;
uint8_t BASE_true;
uint8_t seed_k_t;
uint32_t new_seed_cnt;
uint32_t max_exon_num_per_read;
int readlen_max;
uint16_t Eindel ; //less than 0.01% of intron length less than 20bp
uint32_t max_intron_length;
int waitingLen;

char *sam_path = NULL;
char *anno_load_script = NULL;
FILE *fp_sam = NULL;

//varibale in this file
uint8_t k_first_level = 14;
uint32_t map2ref_cnt = 0;
uint32_t *map2ref_cnt_arr = NULL;
int THREAD_READ_I;
int TOTAL_READ_COUNTs;
char command[1024];
uint8_t read_type = 5;
int *read_len = NULL;
int POS_N_MAX = 0;

static void init_map_param(param_map *opt)
{ 	
	opt->gap_open_D = 2; //4
	opt->gap_ex_D = 1; //2
	opt->gap_open2_D = 24; //24
	opt->gap_ex2_D = 0; //1
	opt->match_D = 1;
	opt->mismatch_D = 2;

	opt->gap_open_R = 2; //2
	opt->gap_ex_R = 1; //1
	opt->gap_open2_R = 32; //32.  if we think the shortest intron is 20bp, then gap_open2_R = 20*gap_ex_R + gap_open_R = 44;
	//if there is a deletion, length is l. if l > (gap_open2_R - gap_open_R)/gap_ex_R , then the deletion is an intron
	opt->gap_ex2_R = 0; //0
	opt->match_R = 1; //1
	opt->mismatch_R = 2; //2

	opt->bw = 2000;
	opt->zdrop_D = 400;
	opt->zdrop_R = 200;
	opt->noncan = 9;
	opt->end_bonus = -1;
	opt->max_extend_gap = 2000;
	opt->max_extend_left_right = 2000;

	opt->error_overall = 0.134;
	opt->error_ins = 0.36;
	opt->error_del = 0.31;

	opt->secondary_ratio = 0.9;
	opt->min_chain_score = 30;
    opt->strand_diff = 20;
    opt->batch_size = 655350;
	opt->max_read_join_gap = 2000;
	opt->max_intron_length = SPLICDISTANCE;
	// opt->max_sw_mat = opt->max_read_join_gap * opt->max_intron_length;
	// opt->max_sw_mat = 4000000000;
	opt->max_sw_mat = 1000000000;
	opt->Eindel = 20;
	opt->thread_n = 4;
	opt->top_n = TOP_NUM_ALN; //2, 3
	opt->k_t = 22;
	opt->e_shift = 5;
	opt->hash_kmer = 8; //8
	opt->seed_k_t = 15;
	opt->seed_step = 5;
	opt->pos_n_max = 50;
	opt->with_qual = 1;
    opt->with_gtf = 0;
    opt->transcript_strand = 0;
	opt->max_exon_num_per_read = 500; //after statics: human_max = 363, mouse_max = 347, dm_max = 82 
	opt->readlen_max = MAX_READLEN;
}

void init_error(param_map *opt)
{
    if (opt->read_type == 1) //ccs
    {
		opt->error_overall = 0.017;
		opt->error_ins = 0.05;
		opt->error_del = 0.2;
    }
	else if (opt->read_type == 2) //clr
	{
		opt->error_overall = 0.15;
		opt->error_ins = 0.42;
		opt->error_del = 0.22;
	}
	else if (opt->read_type == 3) //on1d
	{
		opt->error_overall = 0.21;
		opt->error_ins = 0.15;
		opt->error_del = 0.37;
	}
	else if (opt->read_type == 4) //ont2d
	{
		opt->error_overall = 0.14;
		opt->error_ins = 0.23;
		opt->error_del = 0.35;
	}
}

static void ksw_gen_mat_D(param_map *opt)
{
	int8_t l,k,m;
	mata_D = (int8_t*)calloc(25, 1);
	mata_R = (int8_t*)calloc(25, 1);
	for (l = k = 0; l < 4; ++l) {
		for (m = 0; m < 4; ++m) {
			mata_D[k] = l == m ? opt->match_D : -(opt->mismatch_D);	/* weight_match : -weight_mismatch */
			mata_R[k] = l == m ? opt->match_R : -(opt->mismatch_R);	/* weight_match : -weight_mismatch */
			k++;
		}
		mata_D[k] = 0; // ambiguous base
		mata_R[k] = 0;
		k++;
	}
	for (m = 0; m < 5; ++m) {
		mata_D[k] = 0;
		mata_D[k] = 0;
		k++;
	}
}


static inline void expand_seed(vertex_u **vertexU, uint32_t *vertexNum, uint32_t *uniseed_length, uint8_t tid)
{
	uint8_t r_c;
	uint32_t i, m;
	uint32_t su_i = 0;

	for (r_c = 0; r_c < 2; ++r_c)
	{	
        if(vertexNum[r_c] == 0)
        {
            uniseed_length[r_c] = 0;
            continue;
        }
		su_i = 0;
		
		for (i = 0; i < vertexNum[r_c]; ++i)
		{
			for (m = 0; m < vertexU[r_c][i].pos_n; ++m)//record the reference position
			{
				uniseed[tid][r_c][su_i].seed_id = i;
				uniseed[tid][r_c][su_i].read_begin = vertexU[r_c][i].read_pos;
				uniseed[tid][r_c][su_i].read_end = vertexU[r_c][i].read_pos + vertexU[r_c][i].length1 - 1;
				uniseed[tid][r_c][su_i].ref_begin =  buffer_p[m + buffer_pp[vertexU[r_c][i].uid]] + vertexU[r_c][i].uni_pos_off - 1;
				uniseed[tid][r_c][su_i].ref_end = uniseed[tid][r_c][su_i].ref_begin + vertexU[r_c][i].length2 - 1;
				uniseed[tid][r_c][su_i].cov = vertexU[r_c][i].cov;
				su_i++;
			}
		}
		uniseed_length[r_c] = su_i;

		//if (su_i > 1)
        //{
        //    //merge anchor in different unipath
        //    qsort(uniseed[tid][r_c], su_i, sizeof(uni_seed), compare_uniseed);
        //    //fprintf(stderr, "after qsort--------------\n");
        //    //show_uniseed(uniseed[tid][r_c], su_i);
        //    int qs, qe, Qs;
        //    int64_t ts, te, Ts;
        //    int cov = 0;
        //    int new_cnt = 0;
        //    qs = uniseed[tid][r_c][0].read_begin;
        //    qe = uniseed[tid][r_c][0].read_end;
        //    ts = uniseed[tid][r_c][0].ref_begin;
        //    te = uniseed[tid][r_c][0].ref_end;
        //    cov = uniseed[tid][r_c][0].cov;
        //    for (i = 1; i < su_i; i++)
        //    {
        //        Qs = uniseed[tid][r_c][i].read_begin;
        //        Ts = uniseed[tid][r_c][i].ref_begin;
        //        if((Qs > qs) && (Ts > ts) && (Qs - qe < waitingLen) && abs((Qs - qe) - (Ts - te)) < Eindel)
        //        {
        //            if (Qs > qe)
        //                cov += uniseed[tid][r_c][i].cov;
        //            else
        //                cov += uniseed[tid][r_c][i].cov - (qe - Qs + 1);
        //            qe = uniseed[tid][r_c][i].read_end;
        //            te = uniseed[tid][r_c][i].ref_end;
        //        }
        //        else
        //        {
        //            uniseed[tid][r_c][new_cnt].read_begin = qs;
        //            uniseed[tid][r_c][new_cnt].read_end = qe;
        //            uniseed[tid][r_c][new_cnt].ref_begin = ts;
        //            uniseed[tid][r_c][new_cnt].ref_end = te;

        //            uniseed[tid][r_c][new_cnt].seed_id = uniseed[tid][r_c][i - 1].seed_id;
        //            uniseed[tid][r_c][new_cnt].cov = cov;

        //            qs = uniseed[tid][r_c][i].read_begin;
        //            qe = uniseed[tid][r_c][i].read_end;
        //            ts = uniseed[tid][r_c][i].ref_begin;
        //            te = uniseed[tid][r_c][i].ref_end;
        //            cov = uniseed[tid][r_c][i].cov;

        //            new_cnt++;
        //        }
        //    }
        //    uniseed[tid][r_c][new_cnt].read_begin = qs;
        //    uniseed[tid][r_c][new_cnt].read_end = qe;
        //    uniseed[tid][r_c][new_cnt].ref_begin = ts;
        //    uniseed[tid][r_c][new_cnt].ref_end = te;

        //    uniseed[tid][r_c][new_cnt].seed_id = uniseed[tid][r_c][i - 1].seed_id;
        //    uniseed[tid][r_c][new_cnt].cov = cov;
        //    new_cnt++;

        //    uniseed_length[r_c] = new_cnt;
        //}
	}
}


//13993
static void get_skeleton_anchor(float dist_max, uint32_t dist_max_index, PATH_t *dist_path, uint32_t vertexNum, uint8_t rc_i, uint8_t tid, uint8_t *out_degree, dpSkeleton_t *dp_skeleton);
static void get_skeleton_anchor_FW_RV(PATH_t *dist_path0, PATH_t *dist_path1, float dist_max0, uint32_t dist_max_index0, float dist_max1, uint32_t dist_max_index1, uint32_t *vertexNum, uint8_t tid, uint8_t **out_degree, dpSkeleton_t *dp_skeleton);

int single_seed_reduction_core_single64(uint64_t (*read_bit)[((MAX_READLEN - 1) >> 5) + 1], uint32_t read_length, uint32_t seqi, uint8_t tid, dpSkeleton_t *dp_skeleton)
{
	int32_t j;
	uint32_t read_off = 0;
	uint32_t r_b_v = 0;
	uint8_t re_d = 0;
	uint64_t kmer_bit = 0;
	uint32_t seed_hash = 0;
	uint32_t seed_kmer = 0;
	int64_t seed_id_r = 0;
	uint32_t ref_pos_n = 0;
	uint32_t uni_offset_s_l = 0;
	uint32_t uni_offset_s_r = 0;

	uint32_t left_i = 1;
	uint32_t right_i = 1;
	uint32_t k_way;
	uint64_t hit_i;
	uint64_t hit_binary0, hit_binary1;
	uint32_t mem_i;
	uint32_t max_right_i = 0;


	uint32_t read_pos = 0;
	uint32_t mem_length = 0;
	int64_t seed_binary[2] = {0,0};

	uint32_t memid[2] = {0,0};
	uint32_t uniseed_length[2] = {0,0}; 

	int result = 0;
	
	uint8_t r_i = 0;
	uint32_t su_i = 0;

	float dist_0 = 0;
	float dist_1 = 0;
	uint64_t uni_id_temp;

	int t = (k_first_level > seed_k_t)? (k_first_level - seed_k_t) : 0;
	uint32_t hash_upper_bound = 0;
	uint32_t hash_low_bound = 0;

	uint64_t kmer_pos_uni = 0;

	for ( r_i = 0; r_i < 2; ++r_i)
	{
		k_way = 1;
		mem_i = 0;
		r_b_v = 0;
        read_off = 0;
		for(read_off = 0; read_off <= read_length - seed_k_t; read_off += seed_step) //every seed_l[tid] we choose a seed
		{
			if(read_off + seed_k_t - 1 <= r_b_v)
			{
				continue;
			}

			re_d = (read_off & 0X1f);/////re_d = read_off & 00011111 = read_off while read_off < 31 , 32 loop

			if(re_d <= re_b)  //re_b = 32 - seed_k_t
			{
				kmer_bit = ((read_bit[r_i][read_off >> 5] & bit_tran_re[re_d]) >> ((re_b - re_d) << 1));
			}
			else
			{
				kmer_bit = (((read_bit[r_i][read_off >> 5] & bit_tran_re[re_d]) << ((re_d - re_b) << 1)) | (read_bit[r_i][(read_off >> 5) + 1] >> ((re_2bt - re_d) << 1)));
			}

			if (seed_k_t == k_first_level) //k = 14
			{
				seed_hash = kmer_bit; //

				seed_binary[0] = 0;
				seed_binary[0] += buffer_hash_g[seed_hash];
				seed_binary[1] = 0;
                seed_binary[1] += buffer_hash_g[seed_hash + 1] - 1;
				if (seed_binary[1] < seed_binary[0])
					continue;
			}
			else
			{
				seed_kmer = (kmer_bit & bit_tran[k_r]); 

				seed_hash = (kmer_bit >> (k_r << 1));
				
				result = binsearch_range(seed_kmer, buffer_kmer_g + buffer_hash_g[seed_hash], buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], seed_binary, seed_offset<<1);

				if (result == -1)
				{
					continue;
				}
				seed_binary[0] += buffer_hash_g[seed_hash];
				seed_binary[1] += buffer_hash_g[seed_hash];
			}

			max_right_i = 1;
			hit_binary0 = seed_binary[0];
			hit_binary1 = seed_binary[1];

			if ((hit_binary1 - hit_binary0 + 1) > uni_pos_n_max)
			{
				continue;
			}
			if (mem_i + hit_binary1 - hit_binary0 + 1 >= new_seed_cnt)
			{
				break;
			}
            
            int ref_cnt = 0;
			for (hit_i = hit_binary0; hit_i <= hit_binary1; ++hit_i)
			{
				kmer_pos_uni = buffer_off_g[hit_i];//this kmer's offset on unipath seq
				//find the UID of this kmer
				seed_id_r = binsearch_interval_unipath64(kmer_pos_uni, buffer_seqf, result_seqf);

				ref_pos_n = buffer_pp[seed_id_r + 1] - buffer_pp[seed_id_r];
                ref_cnt += ref_pos_n;

				//if (ref_cnt > POS_N_MAX)
                if(ref_pos_n > pos_n_max)
				{
                    continue;
					//break;
				}
				uni_offset_s_l = kmer_pos_uni - buffer_seqf[seed_id_r];
				uni_offset_s_r = buffer_seqf[seed_id_r + 1] - (kmer_pos_uni + seed_k_t);
				//extend the kmer to a exact match 
 				for(left_i = 1; (left_i <= uni_offset_s_l) && (left_i <= read_off); left_i++)
 				{
					if(((buffer_seq[(kmer_pos_uni - left_i) >> 5] >> ((31 - ((kmer_pos_uni - left_i) & 0X1f)) << 1)) & 0X3)
					        != ((read_bit[r_i][(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3)
					        ) break;
				}

				for(right_i = 1; (right_i <= uni_offset_s_r) && (right_i <= read_length - read_off - seed_k_t); right_i++)
				{
					if(((buffer_seq[(kmer_pos_uni + seed_k_t - 1 + right_i) >> 5] >> ((31 - ((kmer_pos_uni + seed_k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
					        != ((read_bit[r_i][(read_off + seed_k_t - 1 + right_i) >> 5] >> ((31 - ((read_off + seed_k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
					        ) break;
 				}

				read_pos = read_off + 1 - left_i;
				mem_length = seed_k_t + left_i + right_i - 2;

				vertexm[tid][r_i][mem_i].uid = seed_id_r;
				vertexm[tid][r_i][mem_i].seed_id = k_way;
				vertexm[tid][r_i][mem_i].read_pos = read_pos;
				vertexm[tid][r_i][mem_i].uni_pos_off = uni_offset_s_l + 1 - left_i;
				vertexm[tid][r_i][mem_i].length = mem_length;
				vertexm[tid][r_i][mem_i].pos_n = ref_pos_n;

				if (right_i > max_right_i)
				{
					max_right_i = right_i;
				}
				
				++mem_i;
			}

			k_way++;
			r_b_v = read_off + seed_k_t + max_right_i - 1;
		}

		//merge seeds in the same unipath
		su_i = 0;
		uint32_t s1, e1;
		if (mem_i == 0)
		{
			memid[r_i] = 0;
			continue;
		}
		if (mem_i == 1)
		{
			memid[r_i] = mem_i;
			vertexu[tid][r_i][0].uid = vertexm[tid][r_i][0].uid;
			vertexu[tid][r_i][0].read_pos = vertexm[tid][r_i][0].read_pos;
			vertexu[tid][r_i][0].uni_pos_off = vertexm[tid][r_i][0].uni_pos_off; //whether to set uni_pos_off to the leftest position
			vertexu[tid][r_i][0].pos_n = vertexm[tid][r_i][0].pos_n;
			vertexu[tid][r_i][0].length1 = vertexm[tid][r_i][0].length;
			vertexu[tid][r_i][0].length2 = vertexm[tid][r_i][0].length;
			vertexu[tid][r_i][0].cov = vertexm[tid][r_i][0].length;
			continue;
		}
		qsort(vertexm[tid][r_i], mem_i, sizeof(vertex_m), compare_uniid);

		uni_id_temp = vertexm[tid][r_i][0].uid;
		j = 0;
		uint32_t cov = 0;
		while (j < mem_i)
		{
			s1 = j;
			cov = vertexm[tid][r_i][s1].length;
			j++;
			while ((uni_id_temp == vertexm[tid][r_i][j].uid) && (vertexm[tid][r_i][j].uni_pos_off > vertexm[tid][r_i][j-1].uni_pos_off) && (j < mem_i))
			{
				int diff = (int)(vertexm[tid][r_i][j].read_pos - vertexm[tid][r_i][j-1].read_pos - vertexm[tid][r_i][j-1].length);
				if (diff > waitingLen)
					break;
				if (abs((vertexm[tid][r_i][j].uni_pos_off - vertexm[tid][r_i][j-1].uni_pos_off) - (vertexm[tid][r_i][j].read_pos - vertexm[tid][r_i][j-1].read_pos)) < Eindel)
				{
					cov += (diff > 0)? vertexm[tid][r_i][j].length : (diff + vertexm[tid][r_i][j].length);
					++j;
				}
				else
					break;
			}
			e1 = j - 1;
			vertexu[tid][r_i][su_i].uid = vertexm[tid][r_i][s1].uid;
			vertexu[tid][r_i][su_i].read_pos = vertexm[tid][r_i][s1].read_pos;
			vertexu[tid][r_i][su_i].uni_pos_off = vertexm[tid][r_i][s1].uni_pos_off; //whether to set uni_pos_off to the leftest position
			vertexu[tid][r_i][su_i].pos_n = vertexm[tid][r_i][s1].pos_n;
			vertexu[tid][r_i][su_i].cov = cov;
			cov = 0;
			
			//cal length
			if (s1 == e1)
			{
				vertexu[tid][r_i][su_i].length1 = vertexm[tid][r_i][s1].length;
				vertexu[tid][r_i][su_i].length2 = vertexm[tid][r_i][s1].length;
			}else{
				vertexu[tid][r_i][su_i].length1 = vertexm[tid][r_i][e1].read_pos + vertexm[tid][r_i][e1].length - vertexm[tid][r_i][s1].read_pos;
				vertexu[tid][r_i][su_i].length2 = vertexm[tid][r_i][e1].uni_pos_off + vertexm[tid][r_i][e1].length - vertexm[tid][r_i][s1].uni_pos_off;
			}
			
			uni_id_temp = vertexm[tid][r_i][j].uid;
			++su_i;
		}
		memid[r_i] = su_i;
	}

    if ((memid[0] == 0) && (memid[1] == 0))
    {
		dp_skeleton->multi_n = 0;
        return 1;
    }
    expand_seed(vertexu[tid], memid, uniseed_length, tid);

    PATH_t* path1;
    PATH_t* path2;
	uint32_t max_index1 = 0;
	uint32_t max_index2 = 0;
    path1 = (PATH_t* )malloc(uniseed_length[0]*sizeof(PATH_t));
    path2 = (PATH_t* )malloc(uniseed_length[1]*sizeof(PATH_t));

	uint32_t num = (uniseed_length[0] < uniseed_length[1]) ? uniseed_length[1]: uniseed_length[0];
	uint8_t **out_degree = (uint8_t** )calloc(2, sizeof(uint8_t* ));
	for(j = 0; j < 2; ++j)
	{
		out_degree[j] = (uint8_t* )calloc(num, 1);
	}

    if (uniseed_length[0] != 0)
    {
    	dist_0 = creatGraph(uniseed[tid][0], uniseed_length[0], path1, out_degree[0], &max_index1, tid, max_read_join_gap);
    }else
    {
    	dist_0 = 0;
    }
    if (uniseed_length[1] != 0)
    {
    	dist_1 = creatGraph(uniseed[tid][1], uniseed_length[1], path2, out_degree[1], &max_index2, tid, max_read_join_gap);
    }else
    {
    	dist_1 = 0;
    }

#ifdef Annoation
    fprintf(stderr, "dist_0 = %f, dist_1 = %f\n", dist_0, dist_1);
#endif
	float dist_tmp; 
    uint8_t temp_strand; //1 - 0 +  2 +/-(both)

	int thre = min_chain_score;
	if (dist_0 >= dist_1 && dist_0 > thre)
	{
		dist_tmp = dist_0 * secondary_ratio;
		if (dist_1 > dist_tmp)
			temp_strand = 2;
		else
			temp_strand = 0;
	}
	else if(dist_0 < dist_1 && dist_1 > thre)
	{
		dist_tmp = dist_1 * secondary_ratio;
		if(dist_0 > dist_tmp)
			temp_strand = 2;
		else
			temp_strand = 1;
	}
	else
	{
		dp_skeleton->multi_n = 0;
        if(path1 != NULL)	free(path1);
        if(path2 != NULL)	free(path2);
        for ( j = 0; j < 2; ++j)
        {
            if (out_degree[j] != NULL)	free(out_degree[j]);
        }
        if(out_degree != NULL)	free(out_degree);
		return 0;
	}


    if (temp_strand == 1) //-
    {
    	get_skeleton_anchor(dist_1, max_index2, path2, uniseed_length[1], temp_strand, tid, out_degree[1], dp_skeleton);
    }
	else if (temp_strand == 0)
    {
    	get_skeleton_anchor(dist_0, max_index1, path1, uniseed_length[0], temp_strand, tid, out_degree[0], dp_skeleton);
    }
	else  //both strand
	{
		get_skeleton_anchor_FW_RV(path1, path2, dist_0, max_index1, dist_1, max_index2, uniseed_length, tid, out_degree, dp_skeleton);
	}

	if(path1 != NULL)	free(path1);
	if(path2 != NULL)	free(path2);
	for ( j = 0; j < 2; ++j)
	{
		if (out_degree[j] != NULL)	free(out_degree[j]);
	}
	if(out_degree != NULL)	free(out_degree);
    return 0;
}


static uint32_t adjust_anchor(TARGET_t *target, QUERY_t *query, uint32_t cnt, int *sig)
{
	int32_t i;
    uint32_t ts, te, qs, qe;
	uint32_t cnt_back = 0;
	*sig = 0;
    
    ts = target[0].ts;
	te = target[0].te;
	qs = query[0].qs;
	qe = query[0].qe;
    assert(cnt > 0);
    if(cnt < 2)
    {
		if (qe - qs >= BASE_true)
			//*sig += 1;
            *sig = 1;
        return 1;
    }
	else
    {
        for(i = 1; i < cnt; ++i)
        {
        	//join the gap less than min_gap or have overlap
            if(((query[i].qe > qs) && (target[i].te > ts)) || ((abs)((int)(qs - query[i].qe) - (int)(ts - target[i].te)) < Eindel) && ((int)(qs - query[i].qe) < waitingLen))
            {
                qs = query[i].qs;
                ts = target[i].ts;
            }
            else
            {
				query[cnt_back].qs = qs;
				query[cnt_back].qe = qe;
				target[cnt_back].ts = ts;
				target[cnt_back].te = te;

				if (qe - qs >= BASE_true)
					//*sig += 1;
                    *sig = 1;

				ts = target[i].ts;
				te = target[i].te;
				qs = query[i].qs;
				qe = query[i].qe;   
                cnt_back++;
            }
        }
        //record the last one
		query[cnt_back].qs = qs;
		query[cnt_back].qe = qe;
		target[cnt_back].ts = ts;
		target[cnt_back].te = te;

		if (qe - qs >= BASE_true)
			//*sig += 1;
            *sig = 1;
		cnt_back++;
    }
    return cnt_back;
}

static uint32_t find_min_index(PATH_t *dist_path, uint32_t max_index)
{
	int32_t j;
	int32_t min_index;
	min_index = max_index;
	j = dist_path[max_index].pre_node;
	while(j != -1)
	{
		min_index = j;
		j = dist_path[j].pre_node;
	}

	return min_index;
}

static uint32_t connect_anchor(TARGET_t *target, QUERY_t *query, PATH_t *dist_path, uint32_t max_index, uint8_t rc_i, uint8_t tid)
{
	uint32_t cnt = 0;
	int32_t j;

	//map the anchor to anchor_map2ref array
	target[cnt].ts = uniseed[tid][rc_i][max_index].ref_begin;
	target[cnt].te = uniseed[tid][rc_i][max_index].ref_end;
	query[cnt].qs = uniseed[tid][rc_i][max_index].read_begin;
	query[cnt].qe = uniseed[tid][rc_i][max_index].read_end;
	cnt++;

	j = dist_path[max_index].pre_node;
	while(j != -1)
	{
		target[cnt].ts = uniseed[tid][rc_i][j].ref_begin;
		target[cnt].te = uniseed[tid][rc_i][j].ref_end;
		query[cnt].qs = uniseed[tid][rc_i][j].read_begin;
		query[cnt].qe = uniseed[tid][rc_i][j].read_end;
		cnt++;

		j = dist_path[j].pre_node;
	}

	return cnt;
}


static void get_skeleton_anchor_FW_RV(PATH_t *dist_path0, PATH_t *dist_path1, float dist_max0, uint32_t dist_max_index0, float dist_max1, uint32_t dist_max_index1, uint32_t *vertexNum, uint8_t tid, uint8_t **out_degree, dpSkeleton_t *dp_skeleton)
{
	int32_t i, j;
	float dist_thre;
	uint32_t num;
	uint32_t cnt = 0;
	uint32_t cnt_back = 0;
	uint32_t max_index;
	uint8_t multi_cnt = 1; //multiple map
	uint32_t *TOP_N;
	uint8_t top = 0;
	uint8_t top_0 = 0;
	uint8_t top_1 = 0;

	num = (vertexNum[0] > vertexNum[1])? vertexNum[0] : vertexNum[1];
	TARGET_t *target = (TARGET_t* )calloc(num, sizeof(TARGET_t));
	QUERY_t *query = (QUERY_t* )calloc(num, sizeof(QUERY_t));
	TOP_N = (uint32_t* )calloc(top_n<<1, 4);

	//+
	dist_thre = dist_max0 * secondary_ratio;
	num = vertexNum[0];
	TOP_N[0] = dist_max_index0;
	uint32_t dist_min_index0 = find_min_index(dist_path0, dist_max_index0);

	top_0++;
	for (i = num - 1; i >= 0; --i)
	{
		if ((out_degree[0][i] == 0) && (dist_path0[i].dist > dist_thre) && (i != dist_max_index0))
		{
			if ((top_0 < top_n) && (find_min_index(dist_path0, i) != dist_min_index0))
			{
				TOP_N[top_0++] = i;
			}
			else if (top_0 >= top_n)
				break;
		}
	}
	//-
	dist_thre = dist_max1 * secondary_ratio;
	num = vertexNum[1];
	top = top_0;
	TOP_N[top++] = dist_max_index1;
	uint32_t dist_min_index1 = find_min_index(dist_path1, dist_max_index1);

	top_1++;
	for (i = num - 1; i >= 0; --i)
	{
		if ((out_degree[1][i] == 0) && (dist_path1[i].dist > dist_thre) && (i != dist_max_index1))
		{
			if ((top_1 < top_n) && (find_min_index(dist_path1, i) != dist_min_index1))
			{
				top_1++;
				TOP_N[top++] = i;
			}
			else if (top_1 >= top_n)
				break;
		}
	}

	int sig = 0;
	uint8_t real_multi_n = 0;
	multi_cnt = top;
	anchor_t *Anchor = (anchor_t* )calloc(multi_cnt, sizeof(anchor_t));
	for (i = 0; i < top_0; ++i)
	{
		max_index = TOP_N[i];
		cnt = connect_anchor(target, query, dist_path0, max_index, 0, tid);

#ifdef Annoation
		fprintf(stderr, "top = %d\n", i);
		fprintf(stderr, "before adjust--------------cnt = %u, FOR, head = %d\n", cnt, max_index);
		for (j = 0; j < cnt; ++j)
		{
			fprintf(stderr, "%u-%u-%u-%u\n", query[j].qs, query[j].qe, target[j].ts, target[j].te);
		}
#endif

		cnt_back = adjust_anchor(target, query, cnt, &sig);
		if(sig)
		{
            map2ref_cnt_arr[tid] += cnt_back;
		}
		max_exon_num_per_read = (max_exon_num_per_read < cnt_back)? cnt_back : max_exon_num_per_read;
#ifdef Annoation
		fprintf(stderr, "after adjust--------------cnt_back = %u\n", cnt_back);
		for (j = 0; j < cnt_back; ++j)
		{
			fprintf(stderr, "%u-%u-%u-%u\n", query[j].qs, query[j].qe, target[j].ts, target[j].te);
		}
#endif

		//write the query read information to file 
		/*
			multiple_alignments_cnt strand anchor_count anchor1 anchor2 ... anchorN strand anchor_count anchor1 ......
		*/
		Anchor[real_multi_n].strand = 0;
		Anchor[real_multi_n].anchor_n = cnt_back;
		Anchor[real_multi_n].primary = sig;
		Anchor[real_multi_n].anchor_pos = (uint32_t** )calloc(cnt_back, sizeof(uint32_t* ));
		for(j = 0; j < cnt_back; ++j)
		{
			Anchor[real_multi_n].anchor_pos[j] = (uint32_t* )calloc(4, 4);
		}

		for (j = cnt_back - 1; j >=0; --j)
		{
			// fprintf(fp_tff, "%u\t%u\t%u\t%u\t", target[j].ts, target[j].te, query[j].qs, query[j].qe);
			Anchor[real_multi_n].anchor_pos[j][0] = target[j].ts;
			Anchor[real_multi_n].anchor_pos[j][1] = target[j].te;
			Anchor[real_multi_n].anchor_pos[j][2] = query[j].qs;
			Anchor[real_multi_n].anchor_pos[j][3] = query[j].qe;
		}
		real_multi_n++;
	}

	for (i = top_0; i < top; ++i)
	{
		max_index = TOP_N[i];
		cnt = connect_anchor(target, query, dist_path1, max_index, 1, tid);

#ifdef Annoation
		fprintf(stderr, "top = %d\n", i);
		fprintf(stderr, "before adjust--------------cnt = %u, rev, head = %d\n", cnt, max_index);
		for (j = 0; j < cnt; ++j)
		{
			fprintf(stderr, "%u-%u-%u-%u\n", query[j].qs, query[j].qe, target[j].ts, target[j].te);
		}
#endif

		cnt_back = adjust_anchor(target, query, cnt, &sig);
		if(sig)
		{
            map2ref_cnt_arr[tid] += cnt_back;
		}
		max_exon_num_per_read = (max_exon_num_per_read < cnt_back)? cnt_back : max_exon_num_per_read;
#ifdef Annoation
		fprintf(stderr, "after adjust--------------cnt_back = %u\n", cnt_back);
		for (j = 0; j < cnt_back; ++j)
		{
			fprintf(stderr, "%u-%u-%u-%u\n", query[j].qs, query[j].qe, target[j].ts, target[j].te);
		}
#endif
		//write the query read information to file 
		/*
			multiple_alignments_cnt strand anchor_count anchor1 anchor2 ... anchorN strand anchor_count anchor1 ......
		*/

		Anchor[real_multi_n].strand = 1;
		Anchor[real_multi_n].anchor_n = cnt_back;
		Anchor[real_multi_n].primary = sig;
		Anchor[real_multi_n].anchor_pos = (uint32_t** )calloc(cnt_back, sizeof(uint32_t* ));
		for(j = 0; j < cnt_back; ++j)
		{
			Anchor[real_multi_n].anchor_pos[j] = (uint32_t* )calloc(4, 4);
		}
		// fprintf(fp_tff, "%u\t%u\t", 1, cnt_back);
		for (j = cnt_back - 1; j >=0; --j)
		{
			// fprintf(fp_tff, "%u\t%u\t%u\t%u\t", target[j].ts, target[j].te, query[j].qs, query[j].qe);
			Anchor[real_multi_n].anchor_pos[j][0] = target[j].ts;
			Anchor[real_multi_n].anchor_pos[j][1] = target[j].te;
			Anchor[real_multi_n].anchor_pos[j][2] = query[j].qs;
			Anchor[real_multi_n].anchor_pos[j][3] = query[j].qe;
		}
		real_multi_n++;
	}
	dp_skeleton->multi_n = real_multi_n;
	dp_skeleton->point = Anchor;
	free(target);
	free(query);
	free(TOP_N);
}

static void get_skeleton_anchor(float dist_max, uint32_t dist_max_index, PATH_t *dist_path, uint32_t vertexNum, uint8_t rc_i, uint8_t tid, uint8_t *out_degree, dpSkeleton_t *dp_skeleton)
{
	int32_t i;
	int32_t j;
	int32_t max_index = 0;
	uint32_t cnt = 0;
	uint32_t cnt_back = 0;
	uint8_t multi_cnt = 1; //multiple map

	float dist_thre = dist_max * secondary_ratio;
	// uint8_t N = 5;
	uint32_t *TOP_N;
	uint8_t top = 0;

	TOP_N = (uint32_t* )calloc(top_n, 4);
	TARGET_t *target = (TARGET_t* )calloc(vertexNum, sizeof(TARGET_t));
	QUERY_t *query = (QUERY_t* )calloc(vertexNum, sizeof(QUERY_t));	
	
	TOP_N[0] = dist_max_index;
	//find min index
	uint32_t dist_min_index = find_min_index(dist_path, dist_max_index);
	top++;

	for (i = vertexNum - 1; i >= 0; --i)
	{
		if ((out_degree[i] == 0) && (dist_path[i].dist > dist_thre) && (i != dist_max_index))
		{
			if ((top < top_n) && (find_min_index(dist_path, i) != dist_min_index))
				TOP_N[top++] = i;
			else if (top >= top_n)
				break;
		}
	}
	
	int sig = 0;
	uint8_t real_multi_n = 0;
	multi_cnt = top;
	//record the pos
	anchor_t *Anchor = (anchor_t* )calloc(multi_cnt, sizeof(anchor_t));
	for (i = 0; i < multi_cnt; ++i)
	{
		max_index = TOP_N[i];
		cnt = connect_anchor(target, query, dist_path, max_index, rc_i, tid);
#ifdef Annoation
        fprintf(stderr, "top = %d\n", i);
        fprintf(stderr, "before adjust--------------cnt = %u, strand = %d, head = %d\n", cnt, rc_i, max_index);
        for (j = 0; j < cnt; ++j)
        {
        	fprintf(stderr, "%u-%u-%u-%u\n", query[j].qs, query[j].qe, target[j].ts, target[j].te);
        }
#endif

		cnt_back = adjust_anchor(target, query, cnt, &sig);
		if (sig)
			map2ref_cnt_arr[tid] += cnt_back;
		max_exon_num_per_read = (max_exon_num_per_read < cnt_back)? cnt_back : max_exon_num_per_read;
#ifdef Annoation
		fprintf(stderr, "after adjust--------------cnt_back = %u\n", cnt_back);
		for (j = 0; j < cnt_back; ++j)
		{
			fprintf(stderr, "%u-%u-%u-%u\n", query[j].qs, query[j].qe, target[j].ts, target[j].te);
		}
#endif

		//write the query read information to file 
		/*
			multiple_alignments_cnt strand anchor_count anchor1 anchor2 ... anchorN strand anchor_count anchor1 ......
		*/

		// fprintf(fp_tff, "%u\t%u\t", rc_i, cnt_back);
		Anchor[real_multi_n].strand = rc_i;
		Anchor[real_multi_n].anchor_n = cnt_back;
		Anchor[real_multi_n].primary = sig;
		Anchor[real_multi_n].anchor_pos = (uint32_t** )calloc(cnt_back, sizeof(uint32_t* ));
		for(j = 0; j < cnt_back; ++j)
		{
			Anchor[real_multi_n].anchor_pos[j] = (uint32_t* )calloc(4, 4);
		}
		for (j = cnt_back - 1; j >=0; --j)
		{
			Anchor[real_multi_n].anchor_pos[j][0] = target[j].ts;
			Anchor[real_multi_n].anchor_pos[j][1] = target[j].te;
			Anchor[real_multi_n].anchor_pos[j][2] = query[j].qs;
			Anchor[real_multi_n].anchor_pos[j][3] = query[j].qe;
		}
		real_multi_n++;
	}
	dp_skeleton->multi_n = real_multi_n;
	dp_skeleton->point = Anchor;

	free(target);
	free(query);
	free(TOP_N);
}

int seeding_core(int read_seq_core, uint8_t tid, dpSkeleton_t *dp_skeleton)
{
	uint32_t read_length = 0;
	uint32_t read_length_a = 0;
	uint16_t read_bit_char = 0;

	uint32_t r_i = 0;
	uint32_t seqi = 0;
	uint8_t rc_i;
	uint8_t c_tmp = 0;
	char tmp_char;
	
	seqi = read_seq_core;
	read_length = query_info[seqi].read_length;
    
    readlen_max = (read_length > readlen_max)? read_length : readlen_max;
    
    read_len[tid] += read_length;

	if (read_length < 30)
	{
		dp_skeleton->multi_n = 0;
		return 1;
	}

	read_length_a = read_length - 1;
	read_bit_char = (((uint16_t )((read_length_a >> 5) + 1)) << 3);

	memset(read_bit1[tid][0], 0, read_bit_char);
	memset(read_bit1[tid][1], 0, read_bit_char);

	r_i = 0;
	while ((query_info[seqi].read_seq)[r_i])
	{
		tmp_char = (query_info[seqi].read_seq)[r_i];
		if (tmp_char == 'N')
		{
			//random
			tmp_char = "ACGT"[rand()%4];
		}
		c_tmp = charToDna5n[(uint8_t)tmp_char];

		read_bit1[tid][0][r_i >> 5] |= (((uint64_t )c_tmp) << ((31 - (r_i & 0X1f)) << 1));
		read_bit1[tid][1][(read_length_a - r_i) >> 5] |= (((uint64_t )(c_tmp ^ 0X3)) << ((31 - ((read_length_a - r_i) & 0X1f)) << 1));

		r_i++;
	}

	single_seed_reduction_core_single64(read_bit1[tid], read_length, seqi, tid, dp_skeleton);

	return 0;
}

static void *seeding_core_thread(void *aux)
{
	thread_aln_t *d = (thread_aln_t* )aux;
	int _read_lines;

	while(1)
	{
		pthread_rwlock_wrlock(&rwlock);
		_read_lines = THREAD_READ_I++;
		pthread_rwlock_unlock(&rwlock);

		if (_read_lines < d->seqn)
		{
			seeding_core(_read_lines, d->tid, &d->dp_skeleton[_read_lines]);
		}
		else break;
	}
	return 0;
}

void upper_string(char *string)
{
    while(*string)
    {
        if ( *string >= 'a' && *string <= 'z' )
        {
            *string = *string - 32;
        }
        string++;
    }
}

static void print_sam_head(void)
{
	int r_i;
    for (r_i = 1; r_i < chr_file_n; ++r_i)
    {
        fprintf(fp_sam, "@SQ\tSN:%s\tLN:%u\n", chr_names[r_i], (uint32_t)(chr_end_n[r_i] - chr_end_n[r_i - 1]));
    }

	//write command
    fprintf(fp_sam, "%s\n", command);
}


int load_fasta_1pass(bseq_file_t *bf)
{
	uint32_t read_in = batch_size;
	uint32_t seqii = read_in;
	uint32_t r_i = 0;
	uint32_t r_ii = 0;
	int32_t r_iii = 0;
	uint32_t primary;

	k_r = seed_k_t - k_first_level;
    re_b = 32 - seed_k_t;
    re_bt = (re_b << 1);
	re_2bt = 64 - seed_k_t;

	fp_sam = fopen(sam_path, "w");
	if (fp_sam == NULL)
	{
		fprintf(stderr, "[Wrong] Failed to pen file %s!!!\n", sam_path);
		exit(0);
	}
	print_sam_head();

	fp_tff = fopen(temp_anchor_dir, "w");
	if (fp_tff == NULL)
	{
		fprintf(stderr, "[Wrong] Failed to open file %s!!!\n", temp_anchor_dir);
		exit(0);
	}

	fp_temp = fopen(temp_binary_pos, "wb");
	
	query_info = (READ_t* )calloc(read_in, sizeof(READ_t));

	map2ref_cnt_arr = (uint32_t* )calloc(thread_n, 4); //record the number of anchor in each thread

	vertexm = (vertex_m*** )calloc(thread_n, sizeof(vertex_m** ));
	for ( r_i = 0; r_i < thread_n; ++r_i)
	{
		vertexm[r_i] = (vertex_m** )calloc(2, sizeof(vertex_m* ));
		for(r_ii = 0; r_ii < 2; ++r_ii)
		{
			vertexm[r_i][r_ii] = (vertex_m* )calloc(new_seed_cnt, sizeof(vertex_m));
		}
	}
	if (vertexm == NULL)
	{
		fprintf(stderr, "memory wrong, vertexm\n" );
	}

	vertexu = (vertex_u*** )calloc(thread_n, sizeof(vertex_u** ));
	for ( r_i = 0; r_i < thread_n; ++r_i)
	{
		vertexu[r_i] = (vertex_u** )calloc(2, sizeof(vertex_u* ));
		for(r_ii = 0; r_ii < 2; ++r_ii)
		{
			vertexu[r_i][r_ii] = (vertex_u* )calloc(new_seed_cnt, sizeof(vertex_u));
		}
	}
	if (vertexu == NULL)
	{
		fprintf(stderr, "memory wrong, vertexu\n" );
	}

	uniseed = (uni_seed*** )calloc(thread_n, sizeof(uni_seed** ));
	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		uniseed[r_i] = (uni_seed** )calloc(2, sizeof(uni_seed* ));
		for(r_ii = 0; r_ii < 2; ++r_ii)
		{
			uniseed[r_i][r_ii] = (uni_seed* )calloc(new_seed_cnt*pos_n_max, sizeof(uni_seed));
		}
	}
	if (uniseed == NULL)
	{
		fprintf(stderr, "memory wrong, uniseed\n" );
	}
    
    read_len = (int* )calloc(thread_n, sizeof(int));

	pthread_rwlock_init(&rwlock, NULL);
	thread_aln_t* aux;
	aux = (thread_aln_t* )calloc(thread_n, sizeof(thread_aln_t));
    double t_s;
    while(seqii == read_in)
    {
        t_s = realtime();
    	seqii = bseq_read(bf, read_in, query_info);

        TOTAL_READ_COUNTs += seqii;

		dpSkeleton_t* dp_skeleton;
			
		dp_skeleton = (dpSkeleton_t *)calloc(seqii, sizeof(dpSkeleton_t));

		THREAD_READ_I = 0;
		if (thread_n <= 1)
		{
			uint32_t seqii_i;
			for(seqii_i = 0; seqii_i < seqii; seqii_i++)
			{
				seeding_core(seqii_i, 0, &dp_skeleton[seqii_i]);
			}
		}
		else
		{
			pthread_t* tid;
			pthread_attr_t attr;

			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

			tid = (pthread_t* )calloc(thread_n, sizeof(pthread_t));

			for(r_i = 0; r_i < thread_n; ++r_i)
			{
				aux[r_i].tid = r_i;
				aux[r_i].seqn = seqii;
				aux[r_i].dp_skeleton = dp_skeleton;

				int res = pthread_create(&tid[r_i], &attr, seeding_core_thread, aux + r_i);
				if(res != 0)
				{
					fprintf(stderr, "create pthread error");
					exit(1);
				}
			}

			for(r_i = 0; r_i < thread_n; ++r_i)	pthread_join(tid[r_i], 0);
			free(tid);
		}
		//write all skeleton to temp file
		int cnt_tmp = 0;
		for (r_i = 0; r_i < thread_n; ++r_i)
		{
			map2ref_cnt += map2ref_cnt_arr[r_i];
			cnt_tmp += map2ref_cnt_arr[r_i];
		}
		TARGET_t* anchor_map2ref = (TARGET_t* )malloc(cnt_tmp*sizeof(TARGET_t));
		uint32_t write_size = 0;
		for(r_i = 0; r_i < seqii; ++r_i)
		{
			fprintf(fp_tff, "%u\t", dp_skeleton[r_i].multi_n);
			for(r_ii = 0; r_ii < dp_skeleton[r_i].multi_n; ++r_ii)
			{
                primary = dp_skeleton[r_i].point[r_ii].primary;

				fprintf(fp_tff, "%u\t%u\t%u\t", dp_skeleton[r_i].point[r_ii].strand, primary, dp_skeleton[r_i].point[r_ii].anchor_n);

				for(r_iii = dp_skeleton[r_i].point[r_ii].anchor_n - 1; r_iii >= 0 ; --r_iii)
				{
					//pos in reference
					if (primary)
					{
						anchor_map2ref[write_size].ts = dp_skeleton[r_i].point[r_ii].anchor_pos[r_iii][0];
						anchor_map2ref[write_size].te = dp_skeleton[r_i].point[r_ii].anchor_pos[r_iii][1];
						write_size++;
					}
					fprintf(fp_tff, "%u\t%u\t%u\t%u\t", dp_skeleton[r_i].point[r_ii].anchor_pos[r_iii][0], dp_skeleton[r_i].point[r_ii].anchor_pos[r_iii][1], dp_skeleton[r_i].point[r_ii].anchor_pos[r_iii][2], dp_skeleton[r_i].point[r_ii].anchor_pos[r_iii][3]);
				}
			}
			fprintf(fp_tff, "\n");
		}
		//wirte position on reference genome to file for merging and filtering
		assert(write_size==cnt_tmp);
		fwrite(anchor_map2ref, sizeof(TARGET_t), write_size, fp_temp);
		//reset map2ref_cnt_arr
		for (r_i = 0; r_i < thread_n; ++r_i)
		{
			map2ref_cnt_arr[r_i] = 0;
		}
		free(anchor_map2ref);
		//free
		for(r_i = 0; r_i < seqii; ++r_i)
		{
			for (r_ii = 0; r_ii < dp_skeleton[r_i].multi_n; ++r_ii)
			{
				if (dp_skeleton[r_i].point[r_ii].anchor_pos != NULL)	free(dp_skeleton[r_i].point[r_ii].anchor_pos);
			}
			if(dp_skeleton[r_i].point != NULL)	free(dp_skeleton[r_i].point);
		}
		if (dp_skeleton != NULL)	free(dp_skeleton);
		//free
		for (r_i = 0; r_i < read_in; ++r_i)
		{
			if (query_info[r_i].name != NULL)
			{
				free(query_info[r_i].name);
				query_info[r_i].name = NULL;
			}
			if (query_info[r_i].read_seq != NULL)
			{
				free(query_info[r_i].read_seq);
				query_info[r_i].read_seq = NULL;
			}
		}

        int t_len = 0;
        for (r_i = 0; r_i < thread_n; r_i++)
        {
            t_len += read_len[r_i];
            read_len[r_i] = 0;
        }
        
        fprintf(stderr, "[Skeleton-generation] Generating skeletons of %d reads, total %d bases in %f seconds\n", seqii, t_len, realtime() - t_s); 

		//seqii = 0;
    }

	pthread_rwlock_destroy(&rwlock);
	free(aux);
    // free memory
	if (map2ref_cnt_arr != NULL)	free(map2ref_cnt_arr);
    free(read_len);

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		for (r_ii = 0; r_ii < 2; ++r_ii)
		{
			if (vertexm[r_i][r_ii] != NULL)	free (vertexm[r_i][r_ii]);
		}
		if (vertexm[r_i] != NULL)	free(vertexm[r_i]);
	}
	if (vertexm != NULL) free(vertexm);

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		for (r_ii = 0; r_ii < 2; ++r_ii)
		{
			if (vertexu[r_i][r_ii] != NULL)	free (vertexu[r_i][r_ii]);
		}
		if (vertexu[r_i] != NULL)	free(vertexu[r_i]);
	}
	if (vertexu != NULL) free(vertexu);

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		for (r_ii = 0; r_ii < 2; ++r_ii)
		{
			if (uniseed[r_i][r_ii] != NULL)	free (uniseed[r_i][r_ii]);
		}
		if (uniseed[r_i] != NULL)	free(uniseed[r_i]);
	}
	if (uniseed != NULL) free(uniseed);

	for (r_i = 0; r_i < read_in; ++r_i)
	{
		if (query_info[r_i].read_seq != NULL)	free(query_info[r_i].read_seq);
		if (query_info[r_i].name != NULL)	free(query_info[r_i].name);
	}
	if(query_info != NULL)	free(query_info);

	fclose(fp_tff);
	fclose(fp_temp);
	return 0;
}

void init_memory(param_map *opt, char *index_dir)
{
	ksw_gen_mat_D(opt);

	load_index_file(index_dir);

	initGraph();
}

void del_deBGAmemory()
{
	//if (buffer_ref_seq)	free(buffer_ref_seq);   //free after the program finish
	free(buffer_seqf);
	free(buffer_seq);
	free(buffer_pp);
	free(buffer_p);
	free(buffer_hash_g);
	free(buffer_kmer_g);
	free(buffer_off_g);
}

void del_finalmemory()
{
	free(buffer_ref_seq);   //free after the program finish
 	fclose(fp_sam);
	free(mata_D);
	free(mata_R);
}

static int aln_usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:\tde Brijn Graph-based 3rd RNA sequence alignment\n");
	fprintf(stderr, "Usage:\t\tdeSALT aln [options] -f <temporary file> <index_route> <read.fa/fq>\n\n");


    fprintf(stderr, "    -f <temporary file>           The temporary file for storing alignment skeletons in first pass.\n");
    fprintf(stderr, "                                  If users run two deSALT program in the same time, -f option is necessary.\n");
    fprintf(stderr, "    <index_route>                 The path of RdBG index.\n");
    fprintf(stderr, "    <read.fq/fa>                  The input reads in fasta or fastq format.\n\n");

	fprintf(stderr, "Algorithm options:\n\n");
	fprintf(stderr, "    -t --thread           [INT]    Number of threads. [4]\n");
	fprintf(stderr, "    -k --index-kmer       [INT]    K-mer length of RdBG-index. [%u]\n", INDEX_KMER);
	fprintf(stderr, "    -l --seeding-lmer     [INT]    K-mer length of seeding process (no long than RdBG-index). [%u]\n", SEEDING_KMER);
	fprintf(stderr, "    -a --local-hash-kmer  [INT]    K-mer length of local hash process. [%u]\n", LOCAL_HASH_KMER);
	fprintf(stderr, "    -s --seed-step        [INT]    The interval of seeding. [%u]\n", SEED_STEP);
    fprintf(stderr, "    -B --batch-size       [INT]    The number of reads to be processed in one loop. [%u]\n", BATCH_SIZE);
	fprintf(stderr, "    -n --max-uni-pos      [INT]    Maximum allowed number of hits per seed. [%u]\n", MAX_UNI_POS);
	fprintf(stderr, "    -L --max-readlen      [INT]    Maximum allowed read length. [%u]\n", MAX_READLEN);
	fprintf(stderr, "    -i --min-frag-dis     [INT]    Maximum allowed distance of two fragment can be connect. [%u]\n", MIN_FRAG_DIS);
	fprintf(stderr, "    -I --max-intron-len   [INT]    maximum allowed intron length. [%u]\n", SPLICDISTANCE);
	fprintf(stderr, "    -c --min-chain-score  [INT]    minimal skeleton score(match bases minus gap penalty). [%u]\n", MIN_CHAIN_SCORE);
    fprintf(stderr, "    -d --strand-diff      [INT]    The minimal difference of dp score by two strand to make sure the transcript strand. [%d]\n", STRAND_DIFF);
	fprintf(stderr, "    -g --max-read-gap     [INT]    Maximum allowed gap in read when chaining. [%u]\n", MAX_READ_JOIN_GAP);
	fprintf(stderr, "    -p --secondary-ratio  [FLOAT]  Min secondary-to-primary score ratio. [%.2f]\n", SECONDARY_TO_PRIMARY);
    fprintf(stderr, "    -e --e-shift          [INT]    The number of downstream (upstream) exons will be processed when left (right) extension. [%u]\n", E_SHIFT);
    fprintf(stderr, "    -T --trans-strand              Find splicing site according to transcript strand\n");
    fprintf(stderr, "    -G --gtf              [STR]    Provided annotation information for precise intron donor and acceptor sites.\n");
    fprintf(stderr, "                                   Convert GTF file(now support GTF format only) to fixed format of deSALT by Annotation_Load.py \n");
	fprintf(stderr, "    -x --read-type        [STR]    Specifiy the type of reads and set multiple paramters unless overriden.\n");
	fprintf(stderr, "                                   [null] default parameters.\n");
	fprintf(stderr, "                                   ccs (PacBio SMRT CCS reads): error rate 1%%\n");
	fprintf(stderr, "                                   clr (PacBio SMRT CLR reads): error rate 15%%\n");
	fprintf(stderr, "                                   ont1d (Oxford Nanopore 1D reads): error rate > 20%%\n");
	fprintf(stderr, "                                   ont2d (Oxford Nanopore 2D reads): error rate > 12%%\n\n");

	fprintf(stderr, "Scoring options\n\n");
	fprintf(stderr, "    -O --open-pen         [INT(,INT)]\n");
	fprintf(stderr, "                                   Gap open penealty. [%u,%u]\n", GAP_OPEN, GAP_OPEN2);
	fprintf(stderr, "    -E --ext-pen          [INT(,INT)]\n");
	fprintf(stderr, "                                   Gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2}. [%u,%u]\n", GAP_EXT, GAP_EXT2);
	fprintf(stderr, "    -m --match-score      [INT]    Match score for SW-alginment. [%u]\n", MATCH_SCORE);
	fprintf(stderr, "    -M --mis-score        [INT]    Mismatch score for SW-alignment. [%u]\n", MISMATCH_SCORE);
	fprintf(stderr, "    -z --zdrop            [INT(,INT)]\n");
	fprintf(stderr, "                                   Z-drop score for splice/non-splice alignment. [%u]\n", ZDROP_SCORE);
	fprintf(stderr, "    -w --band-width       [INT]    Bandwidth used in chaining and DP-based alignment. [%u]\n", BANDWIDTH);
    fprintf(stderr, "    -R --noncan           [INT]    Penalty score for non-canonical splice junction sites. [%u]\n\n", NONCAN);

	fprintf(stderr, "Output options\n\n");
	fprintf(stderr, "    -N --top-num-aln      [INT]    Max allowed number of secondary alignment. [%u]\n", TOP_NUM_ALN);
	fprintf(stderr, "    -Q --without-qual              Don't output base quality in SAM\n");
	fprintf(stderr, "    -f --temp-file-perfix [STR]    Route of temporary files after the first-pass alignment. [%s]\n", TEMP_FILE_PERFIRX);
	fprintf(stderr, "                                   If you run more than one deSALT program in the same time, \n");
	fprintf(stderr, "                                   you must point out different routes of temporary files for each program!!!\n");
	fprintf(stderr, "                                   If no, every deSALT program will write temporary data to the same file which \n");
	fprintf(stderr, "                                   will cause crash of program in 2-pass alignment due to inconsistent temporary data.\n");
	fprintf(stderr, "    -o --output           [STR]    Output file (SAM format). [%s]\n", OUTPUT);
	
	return 1;
}

int help_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	deSALT (Third generation RNA sequence alignment)\n");

	fprintf(stderr, "Usage:		deSALT <command> [options]\n\n");
	fprintf(stderr, "Command: \n");
	fprintf(stderr, "		index		index reference sequence\n");
	fprintf(stderr, "		aln		align long RNA sequence to reference\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:	deSALT index <ref.fa> <index_route>\n");
	fprintf(stderr, "		build deBGA index file using default k-mer length of deBGA. You can get more deBGA information from https://github.com/HongzheGuo/deBGA");
	fprintf(stderr, "\n\n");
	fprintf(stderr, "Usage:	deSALT aln [options] -f <temporary file> <index_route> <read.fa/fq>\n\n");

    fprintf(stderr, "    -f <temporary file>           The temporary file for storing alignment skeletons in first pass.\n");
    fprintf(stderr, "                                  If users run two deSALT program in the same time, -f option is necessary.\n");
    fprintf(stderr, "    <index_route>                 The path of RdBG index.\n");
    fprintf(stderr, "    <read.fq/fa>                  The input reads in fasta or fastq format.\n\n");

	fprintf(stderr, "Algorithm options:\n\n");
	fprintf(stderr, "    -t --thread           [INT]    Number of threads. [4]\n");
	fprintf(stderr, "    -k --index-kmer       [INT]    K-mer length of RdBG-index. [%u]\n", INDEX_KMER);
	fprintf(stderr, "    -l --seeding-lmer     [INT]    K-mer length of seeding process (no long than RdBG-index). [%u]\n", SEEDING_KMER);
	fprintf(stderr, "    -a --local-hash-kmer  [INT]    K-mer length of local hash process. [%u]\n", LOCAL_HASH_KMER);
	fprintf(stderr, "    -s --seed-step        [INT]    The interval of seeding. [%u]\n", SEED_STEP);
    fprintf(stderr, "    -B --batch-size       [INT]    The number of reads to be processed in one loop. [%u]\n", BATCH_SIZE);
	fprintf(stderr, "    -n --max-uni-pos      [INT]    Maximum allowed number of hits per seed. [%u]\n", MAX_UNI_POS);
	fprintf(stderr, "    -L --max-readlen      [INT]    Maximum allowed read length. [%u]\n", MAX_READLEN);
	fprintf(stderr, "    -i --min-frag-dis     [INT]    Maximum allowed distance of two fragment can be connect. [%u]\n", MIN_FRAG_DIS);
	fprintf(stderr, "    -I --max-intron-len   [INT]    maximum allowed intron length. [%u]\n", SPLICDISTANCE);
	fprintf(stderr, "    -c --min-chain-score  [INT]    minimal skeleton score(match bases minus gap penalty). [%u]\n", MIN_CHAIN_SCORE);
    fprintf(stderr, "    -d --strand-diff      [INT]    The minimal difference of dp score by two strand to make sure the transcript strand. [%d]\n", STRAND_DIFF);
	fprintf(stderr, "    -g --max-read-gap     [INT]    Maximum allowed gap in read when chaining. [%u]\n", MAX_READ_JOIN_GAP);
	fprintf(stderr, "    -p --secondary-ratio  [FLOAT]  Min secondary-to-primary score ratio. [%.2f]\n", SECONDARY_TO_PRIMARY);
    fprintf(stderr, "    -e --e-shift          [INT]    The number of downstream (upstream) exons will be processed when left (right) extension. [%u]\n", E_SHIFT);
    fprintf(stderr, "    -T --trans-strand              Find splicing sites according to transcript strand\n");
    fprintf(stderr, "    -G --gtf              [STR]    Provided annotation information for precise intron donor and acceptor sites.\n");
    fprintf(stderr, "                                   Convert GTF file(now support GTF format only) to fixed format of deSALT by Annotation_Load.py \n");
	fprintf(stderr, "    -x --read-type        [STR]    Specifiy the type of reads and set multiple paramters unless overriden.\n");
	fprintf(stderr, "                                   [null] default parameters.\n");
	fprintf(stderr, "                                   ccs (PacBio SMRT CCS reads): error rate 1%%\n");
	fprintf(stderr, "                                   clr (PacBio SMRT CLR reads): error rate 15%%\n");
	fprintf(stderr, "                                   ont1d (Oxford Nanopore 1D reads): error rate > 20%%\n");
	fprintf(stderr, "                                   ont2d (Oxford Nanopore 2D reads): error rate > 12%%\n\n");

	fprintf(stderr, "Scoring options\n\n");
	fprintf(stderr, "    -O --open-pen         [INT(,INT)]\n");
	fprintf(stderr, "                                   Gap open penealty. [%u,%u]\n", GAP_OPEN, GAP_OPEN2);
	fprintf(stderr, "    -E --ext-pen          [INT(,INT)]\n");
	fprintf(stderr, "                                   Gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2}. [%u,%u]\n", GAP_EXT, GAP_EXT2);
	fprintf(stderr, "    -m --match-score      [INT]    Match score for SW-alginment. [%u]\n", MATCH_SCORE);
	fprintf(stderr, "    -M --mis-score        [INT]    Mismatch score for SW-alignment. [%u]\n", MISMATCH_SCORE);
	fprintf(stderr, "    -z --zdrop            [INT(,INT)]\n");
	fprintf(stderr, "                                   Z-drop score for splice/non-splice alignment. [%u]\n", ZDROP_SCORE);
	fprintf(stderr, "    -w --band-width       [INT]    Bandwidth used in chaining and DP-based alignment. [%u]\n", BANDWIDTH);
    fprintf(stderr, "    -R --noncan           [INT]    Penalty score for non-canonical splice junction sites. [%u]\n\n", NONCAN);

	fprintf(stderr, "Output options\n\n");
	fprintf(stderr, "    -N --top-num-aln      [INT]    Max allowed number of secondary alignment. [%u]\n", TOP_NUM_ALN);
	fprintf(stderr, "    -Q --without-qual              Don't output base quality in SAM\n");
	fprintf(stderr, "    -f --temp-file-perfix [STR]    Route of temporary files after the first-pass alignment. [%s]\n", TEMP_FILE_PERFIRX);
	fprintf(stderr, "                                   If you run more than one deSALT program in the same time, \n");
	fprintf(stderr, "                                   you must point out different routes of temporary files for each program!!!\n");
	fprintf(stderr, "                                   If no, every deSALT program will write temporary data to the same file which \n");
	fprintf(stderr, "                                   will cause crash of program in 2-pass alignment due to inconsistent temporary data.\n");
	fprintf(stderr, "    -o --output           [STR]    Output file (SAM format). [%s]\n", OUTPUT);

	return 1;
}

static const char *short_option = "k:l:a:t:s:B:n:N:L:c:d:g:O:E:m:M:w:i:I:R:z:p:e:f:QTG:o:hx:";

static struct option long_option[] = {
	{"index-kmer", required_argument, NULL, 'k'},
	{"seeding-lmer", required_argument, NULL, 'l'},
	{"local-hash-kmer", required_argument, NULL, 'a'},
	{"thread", required_argument, NULL, 't'},
	{"seed-step", required_argument, NULL, 's'},
    {"batch-size", required_argument, NULL, 'S'},
	{"max-uni-pos", required_argument, NULL, 'n'},
	{"top-num-aln", required_argument, NULL, 'N'},
	{"max-readlen", required_argument, NULL, 'L'},
	{"max-exon", required_argument, NULL, 'r'},
	{"min-chain-score", required_argument, NULL, 'c'},
    {"strand-diff", required_argument, NULL, 'd'},
    {"trans_strand", no_argument, NULL, 'T'},
	{"max-read-gap", required_argument, NULL, 'g'},
	{"open-pen", required_argument, NULL, 'O'},
	{"ext-pen", required_argument, NULL, 'E'},
	{"match-score", required_argument, NULL, 'm'},
	{"mis-score", required_argument, NULL, 'M'},
	{"band-width", required_argument, NULL, 'w'},
	{"min-frag-dis", required_argument, NULL, 'i'},
	{"max-intron-len", required_argument, NULL, 'I'},
	{"zdrop", required_argument, NULL, 'z'},
    {"noncan", required_argument, NULL, 'R'},
	{"secondary-to-primary", required_argument, NULL, 'p'},
	{"e-shift", required_argument, NULL, 'e'},
	{"temp-file-perfix", required_argument, NULL, 'f'},
	{"without-qual", no_argument, NULL, 'Q'},
    {"gtf", required_argument, NULL, 'G'},
	{"output", required_argument, NULL, 'o'},
	{"help", no_argument, NULL, 'h'},
	{"read-type", required_argument, NULL, 'x'},
	{0,0,0,0}
};

int desalt_aln(int argc, char *argv[], const char *version)
{ 
	fprintf(stderr, "[Main] deSALT - De Bruijn graph-based Spliced Aligner for Long Transcriptome reads\n");
	param_map *opt = (param_map* )calloc(1, sizeof(param_map));
	init_map_param(opt);
	int c;
	char *p;

    sprintf(command, "@PG\tID:deSALT\tPN:deSALT\tVN:%s\tCL:%s", version, argv[0]);
    for (c = 1; c < argc; ++c) sprintf(command+strlen(command), " %s", argv[c]);

	while((c = getopt_long(argc, argv, short_option, long_option, NULL)) != -1)
	{
		switch(c)
		{
			case 'k': opt->k_t = atoi(optarg); break;
			case 'l': opt->seed_k_t = atoi(optarg); break;
			case 'a': opt->hash_kmer = atoi(optarg); break;
			case 't': opt->thread_n = atoi(optarg); break;
			case 's': opt->seed_step = atoi(optarg); break;
            case 'B': opt->batch_size = atoi(optarg); break;
			case 'n': opt->pos_n_max = atoi(optarg); break;
			case 'N': opt->top_n = atoi(optarg); break;
			case 'L': opt->readlen_max = atoi(optarg); break;
			case 'r': opt->max_exon_num_per_read = atoi(optarg); break;
			case 'c': opt->min_chain_score = atoi(optarg); break;
            case 'd': opt->strand_diff = atoi(optarg); break;
			case 'g': opt->max_read_join_gap = atoi(optarg); break;
			case 'O': opt->gap_open_D = opt->gap_open_R = opt->gap_open2_D = opt->gap_open2_R = strtol(optarg, &p, 10);
						if (*p != 0 && ispunct(*p) && isdigit(p[1])) opt->gap_open2_D = opt->gap_open2_R = strtol(p+1, &p, 10); break;
			case 'E': opt->gap_ex_D = opt->gap_ex_R = opt->gap_ex2_D = opt->gap_ex2_R = strtol(optarg, &p, 10);
						if (*p != 0 && ispunct(*p) && isdigit(p[1])) opt->gap_ex2_D = opt->gap_ex2_R = strtol(p+1, &p, 10); break;
			case 'm': opt->match_D = opt->match_R = atoi(optarg); break;
			case 'M': opt->mismatch_D = opt->mismatch_R = atoi(optarg); break;
			case 'w': opt->bw = atoi(optarg); break;
			case 'i': opt->Eindel = atoi(optarg); break;
			case 'I': opt->max_intron_length = atoi(optarg); break;
			case 'z': opt->zdrop_D = opt->zdrop_R = strtol(optarg, &p, 10);
						if (*p != 0 && ispunct(*p) && isdigit(p[1])) opt->zdrop_R = strtol(p+1, &p, 10); break;
			case 'p': opt->secondary_ratio = atof(optarg); break;
            case 'R': opt->noncan = atoi(optarg); break;
			case 'e': opt->e_shift = atoi(optarg); break;
			case 'f': opt->temp_file_perfix = strdup(optarg); break;
			case 'Q': opt->with_qual = 0; break;
            case 'T': opt->transcript_strand = 1; break;
            case 'G': opt->anno_path = strdup(optarg); opt->with_gtf = 1; break;
			case 'o': opt->sam_path = strdup(optarg); break;
			case 'h': return aln_usage(); break;
			case 'x': if (strcmp(optarg, "ccs") == 0) opt->read_type = 1;
					  else if (strcmp(optarg, "clr") == 0) opt->read_type = 2;
					  else if (strcmp(optarg, "ont1d") == 0) opt->read_type = 3;
					  else if (strcmp(optarg, "ont2d") == 0) opt->read_type = 4;
					  else {
						  fprintf(stderr, "[main:] Unkown parameter: %s\n", optarg);
						  return aln_usage();
					  } 
					  init_error(opt); break;
			default: return aln_usage(); break;
		}
	} 

	if (argc - optind < 3)
		return aln_usage();
    if ((opt->seed_k_t < 14) || (opt->seed_k_t > opt->k_t))
    {
        fprintf(stderr, "Input error: -l cannot be less than 14 or more than %d\n", opt->k_t);
        exit(1);
    }
    if ((opt->thread_n < 1) || (opt->thread_n > 48))
    {
        fprintf(stderr, "Input error: -t cannot be less than 1 or more than 48\n");
        exit(1);
    }
    if ((opt->hash_kmer < 6) || (opt->hash_kmer > 10))
    {
        fprintf(stderr, "Input error: -a cannot be less than 6 or more than 10\n");
        exit(1);
    }
    if ((opt->seed_step < 1) || (opt->seed_step > 10))
    {
        fprintf(stderr, "Input error: -s cannot be less than 1 or more than 10\n");
        exit(1);
    }
    if ((opt->pos_n_max < 30) || (opt->pos_n_max > 5000))
    {
        fprintf(stderr, "Input error: -n cannot be less than 30 or more than 5000\n");
        exit(1);
    }
    if ((opt->top_n < 1) || (opt->top_n > 10))
    {
        fprintf(stderr, "Input error: -N cannot be less than 1 or more than 10\n");
        exit(1);
    } 
    if ((opt->Eindel < 15) || (opt->Eindel > 30))
    {
        fprintf(stderr, "Input error: -i cannot be less than 15 or more than 30\n");
        exit(1);
    }
    if ((opt->secondary_ratio < 0.7) || (opt->secondary_ratio >= 1.0))
    {
        fprintf(stderr, "Input error: -p cannot be less than 0.7 or more than 1.0\n");
        exit(1);
    }
	if ((opt->e_shift > 10) || (opt->e_shift < 2))
	{
		fprintf(stderr, "Input error: -e cannot be less than 2 or more than 10\n");
        exit(1);
	}

    
	fprintf(stderr, "[Param-INFO] deSALT parameters:index-kmer:%d\tseed-lmer:%d\thash-kmer:%d\tthread:%d\tstrand_diff:%d\tidentify junction:%s\n", opt->k_t, opt->seed_k_t, opt->hash_kmer, opt->thread_n, opt->strand_diff, (opt->transcript_strand)? "transcript strand": "both_strand");
	
	char *index_dir;
	char *read_fastq;
	index_dir = strdup(argv[optind + 1]);
	if (index_dir[strlen(index_dir) - 1] != '/') strcat(index_dir, "/");
	read_fastq = strdup(argv[optind + 2]);

	if (opt->sam_path == NULL)
		sam_path = "./aln.sam";
	else
		sam_path = opt->sam_path;
	memset(temp_anchor_dir, 0, 1024);
	memset(temp_binary_pos, 0, 1024);
    if (opt->temp_file_perfix == NULL)
    {
        strcpy(temp_anchor_dir, "./skeletons.lines");

        strcpy(temp_binary_pos, "./skeletons.pos");
    }
    else
    {
        strcpy(temp_anchor_dir, opt->temp_file_perfix);
        strcat(temp_anchor_dir, "1pass_anchor.lines");

        strcpy(temp_binary_pos, opt->temp_file_perfix);
        strcat(temp_binary_pos, "1pass_anchor.pos");
    }

	//variable in this file
	thread_n = opt->thread_n;
	k_t = opt->k_t;
	seed_k_t = opt->seed_k_t;
	Eindel = opt->Eindel;
	top_n = opt->top_n;
    batch_size = opt->batch_size;
	read_type = opt->read_type;
    secondary_ratio = opt->secondary_ratio;
	readlen_max = opt->readlen_max;
	seed_step = opt->seed_step;
	max_exon_num_per_read = opt->max_exon_num_per_read;
	pos_n_max = opt->pos_n_max;
	max_intron_length = opt->max_intron_length;
	max_read_join_gap = opt->max_read_join_gap;
	// opt->max_sw_mat = opt->max_read_join_gap * opt->max_intron_length;
	min_chain_score = opt->min_chain_score;
	seed_offset = k_t - seed_k_t;
    
    seed_num = 10000; // extract at moset 10000 seed every read

    // variable uni_pos_n_max and POS_N_MAX got from experience, which considering the speed and accuracy.
    //uni_pos_n_max
    if (seed_offset < 3)
        uni_pos_n_max = pow(4, seed_offset);
    else if (seed_offset < 6)
        uni_pos_n_max = pow(2, seed_offset);
    else
        uni_pos_n_max = 64; //64

    POS_N_MAX = 25;
    if ((seed_k_t == 15) || (seed_k_t == 16))
        POS_N_MAX = pos_n_max;
    else if (seed_k_t == 14)
        POS_N_MAX = 35;

    new_seed_cnt = seed_num * (uni_pos_n_max + 1);
	//waitlength
	float els = 0.05;
	float error = 0.2;

	if (read_type == 1) //ccs
		error = 0.02;
	else if (read_type == 2) //clr
		error = 0.15; 
	else if (read_type == 3) //ont1d
		error = 0.25;
	else if (read_type == 4)
		error = 0.13;
	float t = log10(els)/log10(1 - pow(1 - error, seed_k_t));
	float q = 1/error - seed_k_t * pow(1 - error, seed_k_t)/(1 - pow(1 - error, seed_k_t));
	waitingLen = (int)(t * q);
	BASE_true = seed_k_t + 1/error;
    if (read_type == 1)
    {
        opt->gap_open_D = opt->gap_open_R = 6;
        opt->gap_open2_D = opt->gap_open2_R = 24;
        opt->mismatch_D = opt->mismatch_R = 4;
        opt->noncan = 5; 
        BASE_true = seed_k_t + 7;
    }
    
    bseq_file_t *bf;
    bf = bseq_open(read_fastq);
    if(bf == 0)
    {
        fprintf(stderr, "[Waring] Wrong input file route or name: %s \n", read_fastq);
        exit(1);
    }

	fprintf(stderr, "[Phase-INFO] Loading Index and Reads\n");
	init_memory(opt, index_dir);

	fprintf(stderr, "[Phase-INFO] Seeding and Chaining Phase (first-pass)\n");
    double tt1 = realtime();
	load_fasta_1pass(bf);
    fprintf(stderr, "[Phase-INFO] Total %d reads were processed in %.3f seconds (first-pass)\n", TOTAL_READ_COUNTs, realtime() - tt1);
    
    //reset batch size for accerlation of 2pass alignment
    if (TOTAL_READ_COUNTs < opt->batch_size)
    {
        opt->batch_size = TOTAL_READ_COUNTs / 2 + 2; 
    }
    else if (TOTAL_READ_COUNTs < opt->batch_size * 5)
    {
        opt->batch_size = TOTAL_READ_COUNTs / 5 + 2;
    }

	/*
	free memory of first-pass
	*/
	del_deBGAmemory();
	delGraph();
	bseq_close(bf);

	fprintf(stderr, "[Phase-INFO] Refined Alignment Phase (second-pass)\n");
	int Total_mapped_reads = 0;
    double tt2 = realtime();
    load_fasta_2pass(map2ref_cnt, opt, read_fastq, &Total_mapped_reads);
    fprintf(stderr, "[Phase-INFO] Total %d reads were mapped to genome in %f seconds (second-pass)\n", Total_mapped_reads, realtime() - tt2);

	del_finalmemory();
	char cmd[2048];
	sprintf(cmd, "rm %s %s", temp_anchor_dir, temp_binary_pos);
	system(cmd);
	free(opt);

	fprintf(stderr, "[Phase-INFO] Finishing Alignment\n");
    fprintf(stderr, "[Phase-INFO] Command: %s\n", command);

	return 0;
}

