/*************************************************************************
	> File Name: hash_index.h
	> Author: 
	> Mail: 
	> Created Time: 2017年11月23日 星期四 22时01分17秒
 ************************************************************************/

#ifndef _HASH_INDEX_H
#define _HASH_INDEX_H

#include "read_seeding.h"

typedef struct _HashEntry
{
    int value_index; //order of kmer appear
    int *value;//kmer position in reference
}hash_entry;

typedef struct _HashTable
{
    hash_entry *bucket;
}hash_table;

//target - query array
typedef struct _HitSeed
{
    // int32_t idx;
    int query_start;
    int target_start;
    int query_subtract_target; //can remove
    int length;
}hit_seed;

int bucket_num;
int value_num;
int hit_num;

void initHashTable(uint8_t h_k);
void freeHashTable();
//int create_hash_table(uint8_t* target, uint32_t target_len, uint8_t hash_kmer);

//float local_hash_process(uint8_t *target, uint8_t *query, uint32_t target_len, uint32_t query_len, uint8_t hash_kmer, uint8_t tid);
int local_hash_process(uint8_t *target, uint32_t target_len, int* cov_score, TARGET_t *anchor_map2ref, uint32_t key1, uint32_t key2, uint8_t tid, uint8_t signal);

#endif
