
#include "bit_operation.h"


char ksw_cigars[4] = {'M','I','D','S'};
uint16_t carry_ten[4] = {1, 10, 100, 1000};

uint32_t anchor_mask_boundary[11] = {0,0X3,0Xf,0X3f,0Xff,0X3ff,0Xfff,0X3fff,0Xffff,0X3ffff,0Xfffff};
uint32_t anchor_mask_boundary_re[10] = {0Xfffff,0X3ffff,0Xffff,0X3fff,0Xfff,0X3ff,0Xff,0X3f,0Xf,0X3};

uint64_t high_bit_mask = 0Xaaaaaaaaaaaaaaaa;
uint64_t low_bit_mask = 0X5555555555555555;

uint64_t high_bit_mask_lv = 0X0202020202020202;
uint64_t low_bit_mask_lv = 0X0101010101010101;

/* This table is used to transform nucleotide letters into numbers. */
const int8_t nt_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
	
//ref_load
uint8_t charToDna5_N2[] =
{
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0,
    /*    A     C           G                    N */
    /*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0,//'Z'
    /*             T */
    /*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0,
    /*    a     c           g                    n */
    /* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*             t */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

//bit operation
uint8_t charToDna5[] =
{
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
    /*    A     C           G                    N */
    /*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 127, 0, 0, 0, 0, 0,//'Z' 6 127
    /*             T */
    /*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
    /*    a     c           g                    n */
    /* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*             t */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
//modification
char charTochar[] =
{
    /*   0 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /*  16 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /*  32 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /*  48 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /*  64 */ 'A', 'A', 'G', 'C', 'G', 'A', 'A', 'G', 'A', 'A', 'A', 'G', 'A', 'A', 'N', 'A',
    /*    A     C           G                    N */
    /*  80 */ 'A', 'A', 'G', 'G', 'T', 'A', 'G', 'A', 'A', 'T', 'A', 'A', 'A', 'A', 'A', 'A',//'Z'
    /*             T */
    /*  96 */ 'A', 'a', 'G', 'c', 'G', 'A', 'A', 'g', 'A', 'A', 'A', 'G', 'A', 'A', 'n', 'A',
    /*    a     c           g                    n */
    /* 112 */ 'A', 'A', 'G', 'G', 't', 'A', 'G', 'A', 'A', 'T', 'A', 'A', 'A', 'A', 'A', 'A',
    /*             t */
    /* 128 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /* 144 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /* 160 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /* 176 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /* 192 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /* 208 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /* 224 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    /* 240 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
};

uint8_t charToDna5n[] =
{
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,//4 the last but one
    /*    A     C           G                    N */
    /*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 127, 0, 0, 0, 0, 0,//'Z'
    /*             T */
    /*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
    /*    a     c           g                    n */
    /* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*             t */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

char Dna5Tochar[5] = {'A','C','G','T','N'};

const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
const uint64_t hff = 0xffffffffffffffff; //binary: all ones
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

//This uses fewer arithmetic operations than any other known
//implementation on machines with fast multiplication.
//It uses 12 arithmetic operations, one of which is a multiply.
int popcount_3(uint64_t x) {
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
}

void strrev1(char *p)
{
  char *q = p;
  while(q && *q) ++q;
  for(--q; p < q; ++p, --q)
    *p = *p ^ *q,
    *q = *p ^ *q,
    *p = *p ^ *q;
}


uint64_t bv_64[65] = {0,1,3,7,0Xf,
					0X1f,0X3f,0X7f,0Xff,
					0X1ff,0X3ff,0X7ff,0Xfff,
					0X1fff,0X3fff,0X7fff,0Xffff,
					0X1ffff,0X3ffff,0X7ffff,0Xfffff,
					0X1fffff,0X3fffff,0X7fffff,0Xffffff,
					0X1ffffff,0X3ffffff,0X7ffffff,0Xfffffff,
					0X1fffffff,0X3fffffff,0X7fffffff,0Xffffffff,
					0X1ffffffff,0X3ffffffff,0X7ffffffff,0Xfffffffff,
					0X1fffffffff,0X3fffffffff,0X7fffffffff,0Xffffffffff,
					0X1ffffffffff,0X3ffffffffff,0X7ffffffffff,0Xfffffffffff,
					0X1fffffffffff,0X3fffffffffff,0X7fffffffffff,0Xffffffffffff,
					0X1ffffffffffff,0X3ffffffffffff,0X7ffffffffffff,0Xfffffffffffff,
					0X1fffffffffffff,0X3fffffffffffff,0X7fffffffffffff,0Xffffffffffffff,
					0X1ffffffffffffff,0X3ffffffffffffff,0X7ffffffffffffff,0Xfffffffffffffff,
					0X1fffffffffffffff,0X3fffffffffffffff,0X7fffffffffffffff,0Xffffffffffffffff
					};
uint64_t bit_tran[33] = {0,3,0Xf,
						0X3f,0Xff,
						0X3ff,0Xfff,
						0X3fff,0Xffff,
						0X3ffff,0Xfffff,
						0X3fffff,0Xffffff,
						0X3ffffff,0Xfffffff,
						0X3fffffff,0Xffffffff,
						0X3ffffffff,0Xfffffffff,
						0X3fffffffff,0Xffffffffff,
						0X3ffffffffff,0Xfffffffffff,
						0X3fffffffffff,0Xffffffffffff,
						0X3ffffffffffff,0Xfffffffffffff,
						0X3fffffffffffff,0Xffffffffffffff,
						0X3ffffffffffffff,0Xfffffffffffffff,
						0X3fffffffffffffff,0Xffffffffffffffff
						};
uint64_t bit_tran_re[33] = {0Xffffffffffffffff,0X3fffffffffffffff,
                            0Xfffffffffffffff,0X3ffffffffffffff,
                            0Xffffffffffffff,0X3fffffffffffff,
                            0Xfffffffffffff,0X3ffffffffffff,
                            0Xffffffffffff,0X3fffffffffff,
                            0Xfffffffffff,0X3ffffffffff,
                            0Xffffffffff,0X3fffffffff,
                            0Xfffffffff,0X3ffffffff,
                            0Xffffffff,0X3fffffff,
                            0Xfffffff,0X3ffffff,
                            0Xffffff,0X3fffff,
                            0Xfffff,0X3ffff,
                            0Xffff,0X3fff,
                            0Xfff,0X3ff,
                            0Xff,0X3f,
                            0Xf,3,0
                            };


uint64_t bit_assi[33] = {0,0Xc000000000000000, 0Xf000000000000000,
						   0Xfc00000000000000, 0Xff00000000000000,
						   0Xffc0000000000000, 0Xfff0000000000000,
						   0Xfffc000000000000, 0Xffff000000000000,
						   0Xffffc00000000000, 0Xfffff00000000000,
						   0Xfffffc0000000000, 0Xffffff0000000000,
						   0Xffffffc000000000, 0Xfffffff000000000,
						   0Xfffffffc00000000, 0Xffffffff00000000,
						   0Xffffffffc0000000, 0Xfffffffff0000000,
						   0Xfffffffffc000000, 0Xffffffffff000000,
						   0Xffffffffffc00000, 0Xfffffffffff00000,
						   0Xfffffffffffc0000, 0Xffffffffffff0000,
						   0Xffffffffffffc000, 0Xfffffffffffff000,
						   0Xfffffffffffffc00, 0Xffffffffffffff00,
						   0Xffffffffffffffc0, 0Xfffffffffffffff0,
						   0Xfffffffffffffffc, 0Xffffffffffffffff
						};

float mp_sub_bp[] =
{
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -3.295837, -3.295837, -3.295837, -3.295837, -3.295837, 
	/*  48 */ -3.295837, -3.295837, -3.295837, -3.295837, -3.295837, -5.693732, -5.693732, -5.693732, -5.693732, -5.693732, -5.693732, -5.693732, -5.693732, -5.693732, -5.693732, -8.005367, 
    /*  60 */ -8.005367, -8.005367, -8.005367, -8.005367, -8.005367, -8.005367, -8.005367, -8.005367, -8.005367, -10.308852, -10.308852, -10.308852, -10.308852, -10.308852, -10.308852, -10.308852, 
	/*  72 */ -10.308852, -10.308852, -10.308852, -12.611528, -12.611528, -12.611528, -12.611528, -12.611528, -12.611528, -12.611528, -12.611528, -12.611528, -12.611528, -14.914122, -14.914122, -14.914122,
    /*  84 */ -14.914122, -14.914122, -14.914122, -14.914122, -14.914122, -14.914122, -14.914122, -17.216707, -17.216707, -17.216707, -17.216707, -17.216707, -17.216707, -17.216707, -17.216707, -17.216707,
	/*  96 */ -17.216707, -19.519293, -19.519293, -19.519293, -19.519293, -19.519293, -19.519293, -19.519293, -19.519293, -19.519293, -19.519293, -21.821878, -21.821878, -21.821878, -21.821878, 0,
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
