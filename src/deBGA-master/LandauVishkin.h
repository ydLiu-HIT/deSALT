
// Computes the edit distance between two strings and returns a CIGAR string for the edits.
#include <stdlib.h>
#include "bit_operation.h"

//#define	LV_INI

#define	LV_MP

#define bool int
#define _ASSERT assert
#define _uint64 uint64_t
#define CountTrailingZeroes(x, ans) do{ans = __builtin_ctzll(x);} while(0)
#define __min(a, b) (a<b?a:b)
#define MAX_K 720
//500
//410
//820
//65
//129
#define true 1
#define false 0

//short L[MAX_K+1][2*MAX_K +1];
short L[32][MAX_K+1][2*MAX_K +1];
short L_mis[32][MAX_K+1][2*MAX_K +1];

/*
typedef enum {
    COMPACT_CIGAR_STRING = 0,
    EXPANDED_CIGAR_STRING = 1,
    COMPACT_CIGAR_BINARY = 2,
} CigarFormat;
*/

int InitiateLVCompute(short L[32][MAX_K+1][2*MAX_K+1], short thread_n);

int computeEditDistance(
        const char* text, int textLen,
        const char* pattern, int patternLen,
        int k, short L[MAX_K+1][2*MAX_K+1]);

int computeEditDistanceWithCigar(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, bool useM,
    short L[MAX_K+1][2*MAX_K+1]);//CigarFormat format,

int computeEditDistance_s(
        const char* text, int textLen,
        const char* pattern, int patternLen,
        int k, short L[MAX_K+1][2*MAX_K+1], int* s_position);

int computeEditDistanceWithCigar_s(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen,
    short L[MAX_K+1][2*MAX_K+1]);// bool useM, CigarFormat format,

int computeEditDistanceWithCigar_s_nm(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen,
    short L[MAX_K+1][2*MAX_K+1], int* nm_score);

int computeEditDistance_mis(
        const char* text, int textLen,
        const char* pattern, int patternLen,
        int k, short L[MAX_K+1][2*MAX_K+1], uint8_t* quality);
	
int computeEditDistance_mis_s(
        const char* text, int textLen, 
		const char* pattern, int patternLen, 
		int k, short L[MAX_K+1][2*MAX_K+1], uint8_t* quality, int16_t* s_position);
	
int computeEditDistance_misboth(
        const char* text, int textLen,
        const char* pattern, int patternLen,
        int k, short L[MAX_K+1][2*MAX_K+1], short L_mis[MAX_K+1][2*MAX_K+1], uint8_t* quality, uint16_t* mis_n);

int computeEditDistanceWithCigar_s_mis(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen,
    short L[MAX_K+1][2*MAX_K+1], uint8_t* quality);

int computeEditDistanceWithCigar_s_mis_left(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen,
    short L[MAX_K+1][2*MAX_K+1], uint8_t* quality, uint16_t* s_offset);
	
int computeEditDistanceWithCigar_s_left(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen,
    short L[MAX_K+1][2*MAX_K+1], uint16_t* s_offset);
	
	
//for mapping quality
int computeEditDistanceWithCigar_s_mp(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], float* mp_subs, float* sub_t);

int computeEditDistanceWithCigar_s_mis_mp(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], uint8_t* quality, float* mp_subs, float* sub_t);
	
int computeEditDistanceWithCigar_s_mis_left_mp(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], uint8_t* quality, uint16_t* s_offset, float* mp_subs, float* sub_t);

int computeEditDistanceWithCigar_s_nm_left(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], int* nm_score, uint16_t* s_offset);