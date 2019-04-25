
#ifndef BIT_OPERATION_H_
#define BIT_OPERATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

extern char ksw_cigars[4];
extern uint16_t carry_ten[4];
extern uint32_t anchor_mask_boundary[11];
extern uint32_t anchor_mask_boundary_re[10];
extern uint64_t high_bit_mask;
extern uint64_t low_bit_mask;
extern uint64_t high_bit_mask_lv;
extern uint64_t low_bit_mask_lv;
extern char Dna5Tochar[5];
extern uint8_t charToDna5_N2[];
extern uint8_t charToDna5[];
//modification
extern char charTochar[];

extern uint8_t charToDna5n[];
extern const uint64_t m1; //binary: 0101...
extern const uint64_t m2; //binary: 00110011..
extern const uint64_t m4; //binary:  4 zeros,  4 ones ...
extern const uint64_t m8; //binary:  8 zeros,  8 ones ...
extern const uint64_t m16; //binary: 16 zeros, 16 ones ...
extern const uint64_t m32; //binary: 32 zeros, 32 ones
extern const uint64_t hff; //binary: all ones
extern const uint64_t h01; //the sum of 256 to the power of 0,1,2,3...
extern uint64_t bv_64[65];
extern uint64_t bit_tran[33];
extern uint64_t bit_tran_re[33];
extern uint64_t bit_assi[33];
extern const int8_t nt_table[128];
extern float mp_sub_bp[];		
		
//the first two always have same meaning but sometimes not
//char* uint8_t uint8_t uint32_t
//kmer_i must be single

#define kmer_address(s, k, kmer_i, add)\
	add = 0;\
	for(kmer_i = 0; kmer_i < k; kmer_i++)\
		(add) += (((uint64_t )charToDna5[(uint8_t)s[k - 1 - kmer_i]]) << (kmer_i << 1));

//modification
/*
#define kmer_address(s, k, kmer_i, add)\
	add = 0;\
	for(kmer_i = 0; kmer_i < k; kmer_i++)\
		(add) += (((uint64_t )charToDna54[(uint8_t)s[k - 1 - kmer_i]]) << (kmer_i << 1));
*/
		
#define val_address(s, k, kmer_i, add, offset)\
	add = 0;\
	for(kmer_i = 0; kmer_i < k; kmer_i++)\
		(add) += (((uint64_t )s[offset - kmer_i]) << (kmer_i << 1));

		
//k cannot be > 16
#define kmer_bit32(s, k, kmer_i, kmer_int)\
	kmer_int = 0;\
	for(kmer_i = 0; kmer_i < k; kmer_i++)\
		(kmer_int) |= (((uint32_t )charToDna5[(uint8_t)s[k - 1 - kmer_i]]) << (kmer_i << 1));

//copy the characters from index to end to 32-bit (right-alignment)
#define kmer_bit32_index(s, k, kmer_i, kmer_int, index)\
    kmer_int = 0;\
	for(kmer_i = 0; kmer_i < k - index; kmer_i++)\
		(kmer_int) |= (((uint32_t )charToDna5[(uint8_t)s[k - 1 - kmer_i]]) << (kmer_i << 1));

//kmer_i runs from right to left
//copy the characters from index to end to 32-bit[a] (right-alignment), k is the number

//kmer_bit32a(kmer_s, k_t - f, kmer_i, write_buff, 3)
//a << 2
#define kmer_bit32a(s, k, kmer_i, kmer_int, a)\
    memset(kmer_int, 0, 12);\
	for(kmer_i = 0; kmer_i < k; kmer_i++)\
		(kmer_int[a - 1 - (kmer_i >> 4)]) |= (((uint32_t )charToDna5[(uint8_t)s[k - 1 - kmer_i]]) << ((kmer_i & 0Xf) << 1));


#define kmer_bit64(uint64, k, kmer_i, kmer)\
	uint64 = 0;\
	for(kmer_i = 0; kmer_i < k; kmer_i++)\
		(uint64) |= (((uint64_t )charToDna5[(uint8_t )kmer[kmer_i]]) << ((k - 1 - kmer_i) << 1));


//copy the last a 32bits and get the last k chars, then get first f chars

//bit32a_bit64((k_p_file[k_p_i].kp), 3, k_t - f, k_t - k, kmer_i, kmer_l, kmer_f, kmer_s, input, output)
				
#define bit32a_bit64(kmer_int, a, l, s, kmer_i, uintl, uintf, uints, input, output)\
    uintl = 0;\
    uintf = 0;\
    uints = 0;\
    for(kmer_i = 0; kmer_i < s; kmer_i++)\
        (uints) |= (((uint64_t )((kmer_int[a - 1 - (kmer_i >> 4)] >> ((kmer_i & 0Xf) << 1)) & 0X3)) << (kmer_i << 1));\
    for(kmer_i = 0; kmer_i < l; kmer_i++)\
        (uintl) |= (((uint64_t )((kmer_int[a - 1 - (kmer_i >> 4)] >> ((kmer_i & 0Xf) << 1)) & 0X3)) << (kmer_i << 1));\
    for(kmer_i = s; kmer_i < l; kmer_i++)\
        (uintf) |= (((uint64_t )((kmer_int[a - 1 - (kmer_i >> 4)] >> ((kmer_i & 0Xf) << 1)) & 0X3)) << ((kmer_i - (s)) << 1));\
	(output) = ((kmer_int[1] >> 16) & 0xff);\
	(input) = ((kmer_int[1] >> 24) & 0xff);

#define bit32_kmer(kmer_int, k, kmer_i, s)\
	for(kmer_i = 0; kmer_i < k; kmer_i++)\
    {\
        s[kmer_i] = Dna5Tochar[(((kmer_int) >> ((k - 1 - kmer_i) << 1)) & 0X3)];\
    }\
	s[kmer_i] = '\0';

#define base_edge16(input, output, edges)\
	(edges) |= (1 << (7 - charToDna5[(uint8_t)input]));\
	(edges) |= (1 << (3 - charToDna5[(uint8_t)output]));

//copy ori[m:n] to des[0]
#define seq_exact(ori,des,i,m,n)\
    for(i = m; i < n; i++)\
        des[i-(m)] = ori[i];\
    des[i-(m)] = '\0';
//copy ori to des[m:]
#define seq_exact_2(ori,des,i,m)\
    for(i = 0; ori[i] != '\0'; i++)\
        des[m + i] = ori[i];\
    des[m + i] = '\0';

//some operations on graph

#define next_kmer(kmer, kmern, k_t, outedge)\
	kmern = (((kmer & ((one_bit64 << (((k_t) << 1) - 2)) - 1)) << 2) | (outedge));

#define uint64_div(kmerw, kmerf, kmers, k_s)\
	kmerf = (kmerw >> ((k_s) << 1));\
	kmers = (kmerw & ((one_bit64 << ((k_s) << 1)) - 1));

//binary written of unipath file
#define binary64_array(array, s_b, s_b_tmp, v, a_n, a_cnt, fp)\
	if((s_b_tmp + 1) >= (a_n << 2))\
	{\
		fwrite(array, 1, a_n, fp);\
		memset(array, 0, a_n);\
		++a_cnt;\
	}\
	(s_b_tmp) = s_b - (a_n << 2) * (a_cnt);\
	array[(s_b_tmp) >> 2] |= ((uint8_t )v << (((s_b_tmp) & 0X3) << 1));\

#define binary64_array64_f(array, s_b, s_b_tmp, v, a_n, a_cnt, fp)\
	if((s_b_tmp + 1) >= (a_n << 5))\
	{\
		fwrite(array, 8, a_n, fp);\
		memset(array, 0, a_n << 3);\
		++a_cnt;\
	}\
	(s_b_tmp) = s_b - (a_n << 5) * (a_cnt);\
	array[(s_b_tmp) >> 5] |= ((uint64_t )v << ((31 - ((s_b_tmp) & 0X1f)) << 1));\

int popcount_3(uint64_t );
void strrev1(char * );
void fun_1(uint64_t );

#endif /* BIT_OPERATION_H_ */
