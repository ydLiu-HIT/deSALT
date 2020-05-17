/*************************************************************************
	> File Name: test.c
	> Author: 
	> Mail: 
	> Created Time: Tue Apr 21 10:45:08 2020
 ************************************************************************/

#include<stdio.h>
#include<stdint.h>

uint32_t fix_contious_intron(uint32_t *cigar, uint32_t n_cigar)
{
    uint32_t i = 0, j, l = 0;
    int op, nxtop, op_len, conl = 0;
    while(i < n_cigar)
    {
        op = cigar[i]&0xf;
        op_len = (cigar[i] >> 4);
        if(op == 3)
        {
            j = i+1; conl = op_len;
            while(j < n_cigar)
            {
                nxtop = cigar[j]&0xf;
                if (nxtop != 3)
                {
                    if(nxtop == 2 && (j+1 < n_cigar) && ((cigar[j+1]&0xf) == 3))
                    {
                        conl += (cigar[j]>>4);
                        conl += (cigar[j+1]>>4);
                        j += 2;
                        continue;
                    }
                    cigar[l] = conl<<4|3; l++;
                    cigar[l] = cigar[j]; l++;
                    i = j + 1;
                    break;   
                }
                conl += (cigar[j]>>4);
                j += 1;
            }
        }
        else
        {
            cigar[l] = cigar[i]; 
            l++; i++;
        }
    }

    return l;
}


void main(void)
{
    uint32_t cigar[100];
    int opl[100] = {10, 1, 2, 1, 4, 100,200,300,4,1,2};
    int op[100] = {0, 1, 0, 2, 0, 3,2,3,0,1,0};
    int N = 11;
    for(int i = 0; i < N; ++i)
    {
        cigar[i] = opl[i]<<4 | op[i];
    }

    for(int i = 0; i < N; ++i)
    {
        printf("%d%c", cigar[i]>>4, "MIDNS"[cigar[i]&0xf]);
    }
    printf("\n");

    int nN = fix_contious_intron(cigar, N);
    for(int i = 0; i < nN; ++i)
    {
        printf("%d%c", cigar[i]>>4, "MIDNS"[cigar[i]&0xf]);
    }
    printf("\n");

}

