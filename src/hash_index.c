/*************************************************************************
	> File Name: hash_index.c
	> Author: 
	> Mail: 
	> Created Time: 2017年11月23日 星期四 22时01分22秒
 ************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include "hash_index.h"
#include "load_unipath_size.h"

//#define DEBUG
#define HASH_LEN 10000

hash_table* Table = NULL;
hit_seed** hitseed = NULL;
float RATIO_THRE = 0.1;
uint8_t hash_kmer = 8;

/*
 *  需要注意的地方：check_realign() 中thre 设置成10,而不是20. 因为如果是20的话，后面anchor之间只有一个短exon的时候默认是存在的，但是如果不存在，thre=20时候无法修正
 */

void initHashTable(uint8_t h_k)
{
    int i;
    uint8_t r_i;
    hash_kmer = h_k;
    bucket_num = 1 << (hash_kmer << 1);
    value_num = 50;  // for the same hash-key, we suppose there have at most 50 copies in the target, or we think it is repeat sequence
    // hit_num = readlen_max;
    hit_num = HASH_LEN;

    Table = (hash_table* )calloc(thread_n, sizeof(hash_table));
    if (Table == NULL) 
    {
        fprintf(stderr, "memory allocate for hash table wrong!!!\n");
        exit(0);
    }
    for(r_i = 0; r_i < thread_n; ++r_i)
    {
        Table[r_i].bucket = (hash_entry* )malloc(bucket_num * sizeof(hash_entry));
        for (i = 0; i < bucket_num; ++i)
        {
            Table[r_i].bucket[i].value_index = 0;
            Table[r_i].bucket[i].value = (int* )calloc(value_num ,4);
        }
    }
    //init hitseed
    hitseed = (hit_seed** )calloc(thread_n, sizeof(hit_seed* ));
    for(r_i = 0; r_i < thread_n; ++r_i)
    {
        hitseed[r_i] = (hit_seed* )calloc(hit_num , sizeof(hit_seed)); //the threshold should be change
    }
    
}

void freeHashTable()
{
    uint8_t r_i;
    int i;
    for(r_i = 0; r_i < thread_n; ++r_i)
    {
        for (i = 0; i < bucket_num; ++i)
        {
            if (Table[r_i].bucket[i].value != NULL) free(Table[r_i].bucket[i].value);
        }
        if(Table[r_i].bucket != NULL)  free(Table[r_i].bucket);
    }
    if (Table != NULL)  free(Table);

    for(r_i = 0; r_i < thread_n; ++r_i)
    {
        if (hitseed[r_i] != NULL)   free(hitseed[r_i]);
    }
    if (hitseed != NULL)    free(hitseed);
}


int create_hash_table(uint8_t* target, uint32_t target_len, uint8_t tid)
{
    int i, j;
    int len = target_len - hash_kmer;
    int key = 0;
    int key_part = 0; //the [1...hash_kmer-1] bit key
    int index = 0;
    hash_table *table = &Table[tid];
    
    //init the table
    for (i = 0; i < bucket_num; ++i)
    {
        table->bucket[i].value_index = 0;
    }

    //the first kmer i = 0
    for (j = 0; j < hash_kmer - 1; ++j)
    {
        key_part += target[hash_kmer - 1 - j]*(1 << (j << 1));
    }
    key = key_part + target[0]*(1 << ((hash_kmer - 1) << 1));
    index = table->bucket[key].value_index;
    table->bucket[key].value[index] = 0;//position in target of the kmer
    table->bucket[key].value_index += 1;
    //other kmers
    for (i = 1; i <= len; ++i)
    {
        key = (key_part << 2) + target[hash_kmer - 1 + i];
        if (table->bucket[key].value_index >= value_num)
        {
            key_part = key - target[i]*(1 << ((hash_kmer - 1) << 1));
            continue;
        }
        table->bucket[key].value[table->bucket[key].value_index] = i; //position in target of the kmer
        table->bucket[key].value_index += 1;
        //cal the kmer's [1...hash_kmer-1] bit key as the next kmer's top hash_kmer-2 bit
        key_part = key - target[i]*(1 << ((hash_kmer - 1) << 1));
    }
    return 1;
}

static int compare_hitseed(const void * a, const void * b)  //can remove
{
    hit_seed* sm1 = (hit_seed *)a;
    hit_seed* sm2 = (hit_seed *)b;
    
    if(sm1->target_start > sm2->target_start)
        return 1;
    if(sm1->target_start < sm2->target_start)
        return -1;
    else
    {
        //according to read_pos and uni_pos_off detected if there is an inversion
        if(sm1->query_start > sm2->query_start)
            return 1;
        if(sm1->query_start < sm2->query_start)
            return -1;
        else
            return 0;         
    }
}

int hash_query(uint8_t* query, uint32_t query_len, uint8_t tid, int *smble)
{
    int i, m;
    int read_off;
    int key = 0;
    int KEY = 0;
    int key_part = 0;
    int value_idx;
    int hits = 0;
    int8_t seed_step = 1;
    int8_t j;
    hash_table *table = &Table[tid];

    *smble = 0;

    uint8_t tmp_hash = hash_kmer;
    uint8_t sig = 0;
    int hash_sub = 0;
    if (query_len < 50)
        tmp_hash = 6;

    // tmp_hash = 8;
    if (tmp_hash < hash_kmer)
    {
        sig = 1;
        hash_sub = (hash_kmer - tmp_hash) << 1;
    }

    //fprintf(stderr, "****tmp_hash = %d*******\n", tmp_hash);
    int len = query_len - tmp_hash;
    int SCORE_THRE = tmp_hash + 5;
    //the first kmer
    for (i = 0; i < tmp_hash - seed_step; ++i)
    {
        key_part += query[tmp_hash - 1 - i]*(1 << (i << 1));
    }
    key = key_part;
    for (i = tmp_hash - seed_step; i < tmp_hash; ++i)
    {
        key += query[tmp_hash - 1 - i]*(1 << (i << 1));
    }

    KEY = key;
    if (sig)
    {
        KEY <<= hash_sub;
    }


    int t_k;
    for (m = 0; m < (1<<hash_sub); ++m)
    {
        t_k = KEY + m;
        value_idx = table->bucket[t_k].value_index;
        if (value_idx)
        {
            //found
            for (i = 0; i < value_idx; ++i)
            {
                hitseed[tid][hits].query_start = 0;
                hitseed[tid][hits].target_start = table->bucket[t_k].value[i];
                hitseed[tid][hits].length = tmp_hash;
                hitseed[tid][hits].query_subtract_target = hitseed[tid][hits].query_start - hitseed[tid][hits].target_start;
                hits++;
            }
        }
        //first kmer end
    }

    int break_d = 0;
    //other kmers
    for (read_off = seed_step; read_off <= len; read_off += seed_step)
    {
        key = (key_part << (seed_step << 1));
        for (i = 0; i < seed_step; ++i)
        {
            key += query[read_off + tmp_hash - 1 - i]*(1 << (i << 1));
        }

        KEY = key;
        if (sig)
        {
            KEY <<= hash_sub;
        }
        for (m = 0; m < (1<<hash_sub); ++m)
        {
            t_k = KEY + m;
            value_idx = table->bucket[t_k].value_index;
            //query the key through hash table
            if (value_idx == 0) //can not found the kmer in hash table
            {
                continue;
            }
            if (hits + value_idx > hit_num)
            {
                break_d = 1;
                break;
            }
            for (i = 0; i < value_idx; ++i)
            {
                // hitseed[tid][hits].idx = hits;
                hitseed[tid][hits].query_start = read_off;
                hitseed[tid][hits].target_start = table->bucket[t_k].value[i];
                hitseed[tid][hits].length = tmp_hash;
                hitseed[tid][hits].query_subtract_target = read_off - hitseed[tid][hits].target_start;
                hits++;
            }
        }
        if (break_d)
            break;
        //cal the next kmer key
        for (j = 0; j < seed_step; ++j)
        {
            key -= query[read_off + j]*(1 << ((tmp_hash - 1 - j) << 1));
        }
        key_part = key;
    }

    int qs;
    int qe;
    int ts;
    int te;
    int hits_new = 0;
#ifdef DEBUG
    //fprintf(stderr, "before co-merge----------------\n");
    //for (i = 0; i < hits; ++i)
    //{
    //    fprintf(stderr, "idx = %u, qs = %u, ts = %u, len = %u\n", i, hitseed[tid][i].query_start, hitseed[tid][i].target_start, hitseed[tid][i].length);
    //}
#endif

    if (hits == 0)
        return 0;

    qs = hitseed[tid][0].query_start;
    qe = qs + tmp_hash - 1;
    ts = hitseed[tid][0].target_start;
    te = ts + tmp_hash - 1;
    int max_l_seed_id = 0;
    int max_l_seed_l = 0;
    for (i = 1; i < hits; ++i)
    {
        if((hitseed[tid][i].target_start > ts) && abs((hitseed[tid][i].query_start - qs) - (hitseed[tid][i].target_start - ts)) < tmp_hash)   //hash_kmer as a param  for user defined
        {
            qe = hitseed[tid][i].query_start + tmp_hash - 1;
            te = hitseed[tid][i].target_start + tmp_hash - 1;
        }else
        {
            hitseed[tid][hits_new].query_start = qs;
            hitseed[tid][hits_new].target_start = ts;
            hitseed[tid][hits_new].length = qe - qs + 1;
            
            if (max_l_seed_l < hitseed[tid][hits_new].length)
            {
                max_l_seed_id = hits_new;
                max_l_seed_l = hitseed[tid][hits_new].length;
            }

            hits_new++;
            qs = hitseed[tid][i].query_start;
            qe = qs + tmp_hash - 1;
            ts = hitseed[tid][i].target_start;
            te = ts + tmp_hash - 1;
        }
    }
    //the last 
    hitseed[tid][hits_new].query_start = qs;
    hitseed[tid][hits_new].target_start = ts;
    hitseed[tid][hits_new].length = qe - qs + 1;
    if (max_l_seed_l < hitseed[tid][hits_new].length)
    {
        max_l_seed_id = hits_new;
        max_l_seed_l = hitseed[tid][hits_new].length;
    }
    hits_new++;

#ifdef DEBUG
    fprintf(stderr, "middle co-merge----------------\n");
    for (i = 0; i < hits_new; ++i)
    {
        fprintf(stderr, "idx = %u, qs = %u, ts = %u, len = %u\n", i, hitseed[tid][i].query_start, hitseed[tid][i].target_start, hitseed[tid][i].length);
    }
#endif
    
    //continue merge by the max anchor
    qs = hitseed[tid][max_l_seed_id].query_start;
    qe = qs + max_l_seed_l - 1;
    ts = hitseed[tid][max_l_seed_id].target_start;
    i = 0;
    m = hits_new;
    hits_new = 0;
    //find merged hit up and down stream
    for (i = max_l_seed_id - 1; i >= 0; --i)
    {
        if (abs(hitseed[tid][i].target_start - ts) > waitingLen)
        {
            continue;
        }
        if(abs((hitseed[tid][i].query_start - qs) - (hitseed[tid][i].target_start - ts)) < tmp_hash && (hitseed[tid][i].target_start != ts))   //hash_kmer as a param  for user defined 
        {
            qs = (qs < hitseed[tid][i].query_start)? qs : hitseed[tid][i].query_start;
            ts = (ts < hitseed[tid][i].target_start)? ts : hitseed[tid][i].target_start;
            // max_l_seed_l += hitseed[tid][i].length;
        }
        else
        {
            hitseed[tid][hits_new].query_start = hitseed[tid][i].query_start;
            hitseed[tid][hits_new].target_start = hitseed[tid][i].target_start;
            hitseed[tid][hits_new].length = hitseed[tid][i].length;
            hits_new ++;
        }
    }

    for (i = max_l_seed_id + 1; i < m; ++i)
    {
        if (abs(hitseed[tid][i].target_start - ts) > waitingLen)
        {
            continue;
        }

        if(abs((hitseed[tid][i].query_start - qs) - (hitseed[tid][i].target_start - ts)) < tmp_hash && (hitseed[tid][i].target_start != ts))   //hash_kmer as a param  for user defined 
        {
            // qs = (qs < hitseed[tid][i].query_start)? qs : hitseed[tid][i].query_start;
            // ts = (ts < hitseed[tid][i].target_start)? ts : hitseed[tid][i].target_start;
            qe = hitseed[tid][i].query_start + hitseed[tid][i].length - 1;
            // max_l_seed_l += hitseed[tid][i].length;
        }
        else
        {
            hitseed[tid][hits_new].query_start = hitseed[tid][i].query_start;
            hitseed[tid][hits_new].target_start = hitseed[tid][i].target_start;
            hitseed[tid][hits_new].length = hitseed[tid][i].length;
            hits_new ++;
        }
    }
    max_l_seed_l = qe - qs + 1;
    hitseed[tid][hits_new].query_start = qs; 
    hitseed[tid][hits_new].target_start = ts;
    hitseed[tid][hits_new].length = max_l_seed_l;
    hits_new ++;
#ifdef DEBUG
    fprintf(stderr, "after co-merge----------------\n");
    for (i = 0; i < hits_new; ++i)
    {
        fprintf(stderr, "idx = %u, qs = %u, ts = %u, len = %u\n", i, hitseed[tid][i].query_start, hitseed[tid][i].target_start, hitseed[tid][i].length);
    }
    fprintf(stderr, "max_l_seed_l = %d, query_len = %d\n", max_l_seed_l, query_len);
#endif

    if(max_l_seed_l >= 25 && max_l_seed_l/(float)query_len > 0.5)
        *smble = 1;

    if(max_l_seed_l >= SCORE_THRE || max_l_seed_l/(float)query_len > RATIO_THRE)
        return 1;
    else
        return 0;
    // return max_cov/(float)query_len;
}

static inline void get_refseq(uint8_t *ref, uint32_t len, uint32_t start)
{
    //check 
	uint32_t m;

    for (m = 0; m < len; ++m) 
    {
        ref[m] = (buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3;
    }
}

int local_hash_process(uint8_t *target, uint32_t target_len, int* cov_score, TARGET_t *anchor_map2ref, uint32_t key1, uint32_t key2, uint8_t tid, uint8_t signal)
{
    uint32_t i;
    uint32_t s_s;
	uint32_t s_len;
    uint8_t *query = NULL;
    int max_l = 0;
    int LL = target_len + target_len*0.15*0.3;

    int sig = 0;

    for (i = key1; i < key2; ++i)
    {
        int l = anchor_map2ref[i].te - anchor_map2ref[i].ts + 1;
        if (max_l < l)
            max_l = l;
    }

    query = (uint8_t *)calloc(max_l, 1);

    create_hash_table(target, target_len, tid);

    for (i = key1; i < key2; ++i)
	{
        sig = 0;
		s_s = anchor_map2ref[i].ts;
		s_len = anchor_map2ref[i].te - anchor_map2ref[i].ts + 1;
		get_refseq(query, s_len, s_s);
#ifdef DEBUG
        fprintf(stderr, "exon %d\n", i - key1);
#endif
        if (signal == 1)
        {
            if (s_len < LL)
            {
                if (s_len < 30)// if there short exon, stitch direct
                {
                    cov_score[i - key1] = 1;
                }
                else
                {        
                    cov_score[i - key1] = hash_query(query, s_len, tid, &sig); //later make qseq to be target for everyone using
                }
            }
            else
            {
#ifdef DEBUG
                fprintf(stderr, "not process!\n");
#endif
                cov_score[i - key1] = 0;
            }
            //if (sig)
            //    LL -= s_len;
        }
        else
        {
            cov_score[i - key1] = hash_query(query, s_len, tid, &sig); //later make qseq to be target for everyone using
        }
		//if ((s_len < LL) && (signal == 1))
        ////if (s_len > hash_kmer)
		//	cov_score[i - key1] = hash_query(query, s_len, tid, &sig); //later make qseq to be target for everyone using
		//else
		//	cov_score[i - key1] = 0;
        
	}
    
    free(query);
    return 0;
}
