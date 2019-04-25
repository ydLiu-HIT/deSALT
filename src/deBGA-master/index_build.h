
#ifndef INDEX_BUILD_H_
#define INDEX_BUILD_H_

#include <stdio.h>
#include <stdint.h>

#include "load_input.h"

//#define	DEBUG_NOT_FOUND
//#define	DEBUG_INPUT
#define	DEBUG_COMBINE

#define	LAST_DEBUG

#define	DEBUG_64BIT

#define REF_64

#define HASH_MAX (1 << 28)
#define KF_LENGTH_MAX 14
#define KS_LENGTH_MAX 20
#define KT_LENGTH_MAX 28
//200
#define REF_FASTA_LINE  500
#define KMER_UP	28
#define KMER_DOWN	20
#define UNI_SEQ_WRI_ARR	1024

const char charN = 'N';
const char charn = 'n';
const char charZ = 'Z';
const char* identifier = ">";
const char* char_N = "N";
const char* char_n = "n";
const char* chars_N = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
const char* suff = ".fa";

uint64_t sta_num_write = 0;
uint64_t ref_seq_n = 0;

//
#ifdef	UNPIPATH_OFF_K20
typedef struct kmer_pos
{
    uint32_t kp[4];
} k_p;
#else
typedef struct kmer_pos
{
    uint32_t kp[3];
} k_p;
#endif

typedef struct pos_uni
{
    uint32_t uniid;
    uint32_t pos;
} p_u;

typedef struct kmer_us
{
    uint64_t kmer;
    uint32_t unioff;
} k_u;

int index_build(int , char * []);
void upper_string(char * );
void load_reffile_kmer();
void load_reffile_kmer_fa();
int compare (const void * , const void * );
int compare_pu (const void * , const void * );
int compare_us (const void * , const void * );
uint32_t file_kmer_qsort();
void build_pos_unipath();

#ifdef	UNPIPATH_OFF_K20
uint64_t uni_pos_combine64(uint64_t [], uint64_t [], uint64_t [], uint64_t , int64_t , uint64_t** );
#else
uint32_t uni_pos_combine(uint32_t [], uint32_t [], uint32_t [], uint32_t , uint32_t , uint32_t** );
#endif
int64_t binsearch_offset_index_com64(uint64_t , uint64_t [], uint64_t , uint64_t );
int64_t binsearch_offset_index64(uint32_t , uint32_t [], uint64_t , uint64_t );
int64_t binsearch_offset_index(uint32_t , uint32_t [], uint32_t , uint32_t );

uint8_t one_number_uint(uint8_t );
uint8_t node_indentity(uint8_t, uint8_t);
uint8_t edgeout_node(uint8_t );


#endif /* INDEX_BUILD_H_ */
