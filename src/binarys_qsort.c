#include <stdint.h>
#include "binarys_qsort.h"
#include "read_seeding.h"
#include "aln_2pass.h"

int binsearch_range(uint64_t key, uint32_t *v, int64_t n,  int64_t *range, int8_t k_off)
{
    int64_t l=0, r=n-1, m;
    uint32_t tmp = 0;
    range[0] = range[1] = -1;

    //printf("k_off = %d\n", k_off);

    if (k_off == 0)
    {
        while(l <= r)
        {
            m = (l + r)/2;
            if (key < v[m])
            {
                r = m - 1;
            }
            else if (key > v[m])
            {
                l = m + 1;
            }
            else
            {
                range[0] = range[1] = m;
                return 1;
            }
        }
    }
    else
    {
        while (l <= r)
        {
            m = (l+r)/2;
            tmp = v[m] >> k_off;
            if (tmp == key)
            {
                range[0] = range[1] = m;
                
                //run low bound
                int64_t sl=l, sr=m-1, sm;
                while (sl <= sr)
                {
                    sm = (sl+sr)/2;
                    tmp = v[sm] >> k_off;
                    if (tmp == key)
                    {
                        range[0] = sm;
                        sr = sm-1;
                    }
                    else if (tmp > key) sr = sm - 1;
                    else    sl = sm + 1;
                }

                //run upper bound
                sl = m+1; sr = r;
                while (sl <= sr)
                {
                    sm = (sl+sr)/2;
                    tmp = v[sm] >> k_off;
                    if (tmp == key)
                    {
                        range[1] = sm;
                        sl = sm+1;
                    }
                    else if (tmp > key) sr = sm - 1;
                    else    sl = sm + 1;
                }
                return 1;
            }
            else if (tmp > key) r = m - 1;
            else l = m + 1;
        }
    }

    return -1;
}

#ifdef UNPIPATH_OFF_K20
int multi_binsearch_offset64(uint32_t x, uint32_t v[], int64_t n, uint64_t offset, int64_t seed_binary[], int8_t k_r)
#else
int multi_binsearch_offset(uint32_t x, uint32_t v[], int64_t n, uint32_t offset, int64_t seed_binary[], int8_t k_r)
#endif
{
    int64_t low, high, mid;
    uint32_t temp = 0;
#ifdef UNPIPATH_OFF_K20
    uint64_t current_offset = 0;
#else
    uint32_t current_offset = 0;
#endif

    int8_t k_2r = k_r << 1;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        temp = v[mid + offset] >> k_2r;
        if(x < temp)
        {
            high = mid - 1;
        }
        else if(x > temp)
        {
            low = mid + 1;
        }
        else  /*found match*/
        {   
            //找到当前match的kmer，然后向两侧search，找到match的上下界
            if (mid == 0)
            {
                seed_binary[0] = offset;

                current_offset = mid + offset + 1;
                while((current_offset < offset + n) &&  (x == (v[current_offset] >> k_2r)))
                {
                    current_offset++;
                }
                seed_binary[1] = current_offset - 1;
            }
            else if (mid == (n - 1))
            {
                seed_binary[1] = n + offset- 1;

                current_offset = mid + offset - 1;
                while((current_offset >= offset) && (x == (v[current_offset] >> k_2r)))
                {
                    current_offset--;
                }
                seed_binary[0] = current_offset + 1;
            }
            else
            {
                //up
                current_offset = mid + offset - 1;
                while((current_offset >= offset) && (x == (v[current_offset] >> k_2r)))
                {
                    current_offset--;
                }
                seed_binary[0] = current_offset + 1;    
                
                //down
                current_offset = mid + offset + 1;
                while((current_offset < offset + n) && (x == (v[current_offset] >> k_2r)))
                {
                    current_offset++;
                }
                seed_binary[1] = current_offset - 1;
            }

            return 1;
        }
    }
    return -1;
}

int64_t binsearch_interval_unipath64(uint64_t x, uint64_t v[], uint64_t n)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        if(x < v[mid])
        {
            high = mid - 1;
        }
        else if(x > v[mid])
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid;
        }
    }

    return high;
}

int64_t binsearch_interval_unipath(uint32_t x, uint32_t v[], uint32_t n)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        if(x < v[mid])
        {
            high = mid - 1;
        }
        else if(x > v[mid])
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid;
        }
    }

    return high;
}

int compare_uniid(const void * a, const void * b)
{
    vertex_m* sm1 = (vertex_m *)a;
    vertex_m* sm2 = (vertex_m *)b;
    
    if(sm1->uid > sm2->uid)
        return 1;
    if(sm1->uid < sm2->uid)
        return -1;
    else
    {
        //according to read_pos and uni_pos_off detected if there is an inversion
        if(sm1->read_pos > sm2->read_pos)
            return 1;
        if(sm1->read_pos < sm2->read_pos)
            return -1;
        else    return 0;         
    }
}

int compare_seedid(const void* a, const void* b)
{
    vertex_m* sm1 = (vertex_m *)a;
    vertex_m* sm2 = (vertex_m *)b;

    if(sm1->seed_id > sm2->seed_id)
        return 1;
    if(sm1->seed_id < sm2->seed_id)
        return -1;
    else
    {
        if(sm1->pos_n > sm2->pos_n)
            return 1;
        if(sm1->pos_n < sm2->pos_n)
            return -1;
        else    return 0;         
    }
}

int compare_mem(const void *a , const void *b)
{
    vertex_m* vertex1 = (vertex_m *)a;
    vertex_m* vertex2 = (vertex_m *)b;

    if (vertex1->read_pos > vertex2->read_pos)
        return 1;
    else if (vertex1->read_pos < vertex2->read_pos)
        return -1;
    else
        return 0;
}

 
int compare_uniseed(const void *a , const void *b)
{
    uni_seed* seed1 = (uni_seed *)a;
    uni_seed* seed2 = (uni_seed *)b;

    if (seed1->ref_end > seed2->ref_end)
        return 1;
    else if (seed1->ref_end < seed2->ref_end)
        return -1;
    else 
    {
        if (seed1->ref_begin > seed2->ref_begin)
            return 1;
        else if (seed1->ref_begin < seed2->ref_begin)
            return -1;
        else
            return 0;
    }

    // if (seed1->read_begin > seed2->read_begin)
    //     return 1;
    // else if (seed1->read_begin < seed2->read_begin)
    //     return -1;
    // else
    // {
    //     if (seed1->ref_begin > seed2->ref_begin)
    //         return 1;
    //     else if (seed1->ref_begin < seed2->ref_begin)
    //         return -1;
    //     else
    //         return 0;
    // }
}

int compare_anchor(const void *a, const void *b)
{
    TARGET_t* anchor1 = (TARGET_t* )a;
    TARGET_t* anchor2 = (TARGET_t* )b;

    if (anchor1->ts > anchor2->ts)
        return 1;
    else if (anchor1->ts < anchor2->ts)
        return -1;
    else
    {
        if (anchor1->te > anchor2->te)
            return 1;
        else if (anchor1->te < anchor2->te)
            return -1;
        else
            return 0;
    }
}

int compare_intron(const void *a, const void *b)
{
    Ival_anno_t* intron1 = (Ival_anno_t* )a;
    Ival_anno_t* intron2 = (Ival_anno_t* )b;

    if (intron1->Is > intron2->Is)
        return 1;
    else if (intron1->Is < intron2->Is)
        return -1;
    else
    {
        if (intron1->Ie > intron2->Ie)
            return 1;
        else if (intron1->Ie < intron2->Ie)
            return -1;
        else
            return 0;
    }
}

int compare_exon(const void *a, const void *b)
{
    Anno_t* exon1 = (Anno_t* )a;
    Anno_t* exon2 = (Anno_t* )b;

    if (exon1->start > exon2->start)
        return 1;
    else if (exon1->start < exon2->start)
        return -1;
    else
    {
        if (exon1->end > exon2->end)
            return 1;
        else if (exon1->end < exon2->end)
            return -1;
        else
            return 0;
    }
   
}


