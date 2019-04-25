
#ifndef BINARY_QSORT_H_
#define BINARY_QSORT_H_

#include "seed_ali_p.h"

//binary search offset
#ifdef UNPIPATH_OFF_K20
int64_t binsearch_offset64(uint32_t x, uint32_t v[], int64_t n, uint64_t offset)
#else
int64_t binsearch_offset(uint32_t x, uint32_t v[], int64_t n, uint32_t offset)
#endif
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        if(x < v[mid + offset])
        {
            high = mid - 1;
        }
        else if(x > v[mid + offset])
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return (int64_t )(mid + offset);
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

int64_t binsearch_seed_pos_s(uint32_t x, uint32_t* bp, int64_t n)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (bp[mid] - devi))
        {
            high = mid - 1;
        }
        else if(x > (bp[mid] + devi))
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return bp[mid];
        }
    }

    return -1;
}

#ifdef UNPIPATH_OFF_K20
int64_t binsearch_seed_pos_ss64(uint64_t x, uint64_t* bp, int64_t n)
#else
int64_t binsearch_seed_pos_ss(uint32_t x, uint32_t* bp, int64_t n)
#endif
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (bp[mid] - devi))
        {
            high = mid - 1;
        }
        else if(x > (bp[mid] + devi))
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid;
        }
    }

    return -1;
}

int64_t binsearch_seed_pos_offset(uint32_t x, uint32_t* bp, int64_t n, uint32_t offset)
{
    int64_t low, high, mid;
	
	low = 0;
    high = n - 1;
	
    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (bp[mid + offset] - devi))
        {
            high = mid - 1;
        }
        else if(x > (bp[mid + offset] + devi))
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return bp[mid + offset];
        }
    }

    return -1;
}

#ifdef UNPIPATH_OFF_K20
int64_t binsearch_pair_pos_reduce64(uint64_t x, uint64_t* bp, int64_t n, int64_t offset, uint8_t tid)
#else
int64_t binsearch_pair_pos_reduce(uint32_t x, uint32_t* bp, int64_t n, uint32_t offset, uint8_t tid)
#endif
{
    int64_t low, high, mid;
	
	low = 0;
    high = n - 1;
	
    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (bp[mid + offset] - devi))
        {
            high = mid - 1;
        }
        else if(x > (bp[mid + offset] + devi))
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid + offset;
        }
    }

	r_low[tid] = low;
	
    return -1;
}

int64_t binsearch_seed_set_offset(uint32_t x, seed_sets* seedsets, int64_t n, uint32_t offset)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (seedsets[mid + offset].seed_set - uni_d))
        {
            high = mid - 1;
        }
        else if(x > (seedsets[mid + offset].seed_set + uni_d))
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid + offset;
        }
    }

    return -1;
}

int64_t binsearch_seed_set_reduce64(uint64_t x, seed_sets* seedsets, int64_t n, uint64_t offset, uint8_t tid)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (seedsets[mid + offset].seed_set - uni_d))//
        {
            high = mid - 1;
        }
        else if(x > (seedsets[mid + offset].seed_set + uni_d))//
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid + offset;
        }
    }
	
	g_low[tid] = low;
	
    return -1;
}

int64_t binsearch_seed_set_reduce(uint32_t x, seed_sets* seedsets, int64_t n, uint32_t offset, uint8_t tid)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (seedsets[mid + offset].seed_set - uni_d))//
        {
            high = mid - 1;
        }
        else if(x > (seedsets[mid + offset].seed_set + uni_d))//
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid + offset;
        }
    }
	
	g_low[tid] = low;
	
    return -1;
}

int64_t binsearch_seed_pos_reduce64(uint64_t x, uint64_t* pos, int64_t n, uint64_t offset, uint8_t tid)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (pos[mid + offset] - uni_d))
        {
            high = mid - 1;
        }
        else if(x > (pos[mid + offset] + uni_d))
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid + offset;
        }
    }
	
	//g_high = high;
	g_low[tid] = low;
	
    return -1;
}

int64_t binsearch_seed_pos_reduce(uint32_t x, uint32_t* pos, int64_t n, uint32_t offset, uint8_t tid)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (pos[mid + offset] - uni_d))
        {
            high = mid - 1;
        }
        else if(x > (pos[mid + offset] + uni_d))
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid + offset;
        }
    }
	
	//g_high = high;
	g_low[tid] = low;
	
    return -1;
}

//binary search offset
//bp[offset] -> bp[offset + n - 1]
int64_t binsearch_seed_pos(uint64_t x, uint32_t* bp, uint32_t n, uint32_t offset)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
        if(x < (bp[mid + offset] - MAX_DIFF_POS))
        {
            high = mid - 1;
        }
        else if(x > (bp[mid + offset] + MAX_DIFF_POS))
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid + offset;
        }
    }

    return -1;
}


int compare_uniid(const void * a, const void * b)
{
	seed_m* sm1 = (seed_m *)a;
    seed_m* sm2 = (seed_m *)b;
	
	if(sm1->uni_id > sm2->uni_id)
		return 1;
	if(sm1->uni_id < sm2->uni_id)
		return -1;
	else
	{
		s_uid_f[sm1->tid] = 1;
/*
		if(sm1->qu_v < sm2->qu_v)
			return 1;
		if(sm1->qu_v > sm2->qu_v)
			return -1;
		else{
			if(sm1->ref_pos_n < sm2->ref_pos_n)
				return 1;
			if(sm1->ref_pos_n > sm2->ref_pos_n)
				return -1;
			else	return 0;
		}
*/
		if(sm1->ref_pos_n < sm2->ref_pos_n)
			return 1;
		if(sm1->ref_pos_n > sm2->ref_pos_n)
			return -1;
		else	return 0;
			
	}
}

int compare_posn(const void * a, const void * b)
{
	seed_m* su1 = (seed_m *)a;
    seed_m* su2 = (seed_m *)b;
	
	if(su1->ref_pos_n > su2->ref_pos_n)
		return 1;
	if(su1->ref_pos_n < su2->ref_pos_n)
		return -1;
	else{
		if(su1->length < su2->length)
			return 1;
		if(su1->length > su2->length)
			return -1;
		else	return 0;
	}

}

int compare_plen(const void * a, const void * b)
{
	seed_pa* pa1 = (seed_pa *)a;
    seed_pa* pa2 = (seed_pa *)b;
	
	if(pa1->length < pa2->length)
		return 1;
	if(pa1->length > pa2->length)
		return -1;
	else	return 0;
	
}

int compare_plen_single(const void * a, const void * b)
{
	seed_pa_single* pa1 = (seed_pa_single *)a;
    seed_pa_single* pa2 = (seed_pa_single *)b;
	
	if(pa1->length < pa2->length)
		return 1;
	if(pa1->length > pa2->length)
		return -1;
	else	return 0;
	
}
int compare_seed_filter_posn(const void * a, const void * b)
{
	seed_pa* pa1 = (seed_pa *)a;
    seed_pa* pa2 = (seed_pa *)b;
	
	if(pa1->pos_n < pa2->pos_n)
		return -1;
	if(pa1->pos_n > pa2->pos_n)
		return 1;
	else	return 0;
	
}

int compare_sets(const void * a, const void * b)
{
	seed_sets* seedsets1 = (seed_sets* )a;
	seed_sets* seedsets2 = (seed_sets* )b;
	
	int8_t cov_r = 0;
	uint8_t cov_i = 0;
	for(cov_i = 0; cov_i < cov_a_n[seedsets1->tid][rc_cov_f[seedsets1->tid]]; cov_i++)//
	{
		if(((uint64_t )((seedsets1->cov)[cov_i])) > ((uint64_t )((seedsets2->cov)[cov_i])))
		{
			cov_r = 1;
			break;
		}
		if(((uint64_t )((seedsets1->cov)[cov_i])) < ((uint64_t )((seedsets2->cov)[cov_i])))
		{
			cov_r = -1;
			break;
		}
	}
	
	if(cov_r == 1)
		return -1;
	if(cov_r == -1)
		return 1;
	else{
		if((seedsets1->seed_set) < (seedsets2->seed_set))
			return -1;
		if((seedsets1->seed_set) > (seedsets2->seed_set))
			return 1;
		else	return 0;
	}
}
int compare_sets_s(const void * a, const void * b)
{
	seed_sets* seedsets1 = (seed_sets* )a;
	seed_sets* seedsets2 = (seed_sets* )b;
	
	int8_t cov_r = 0;
	uint8_t cov_i = 0;
	for(cov_i = 0; cov_i < cov_a_n_s[seedsets1->tid]; cov_i++)//
	{
		if(((uint64_t )((seedsets1->cov)[cov_i])) > ((uint64_t )((seedsets2->cov)[cov_i])))
		{
			cov_r = 1;
			break;
		}
		if(((uint64_t )((seedsets1->cov)[cov_i])) < ((uint64_t )((seedsets2->cov)[cov_i])))
		{
			cov_r = -1;
			break;
		}
	}
	
	if(cov_r == 1)
		return -1;
	if(cov_r == -1)
		return 1;
	else{
		if((seedsets1->seed_set) < (seedsets2->seed_set))
			return -1;
		if((seedsets1->seed_set) > (seedsets2->seed_set))
			return 1;
		else	return 0;
	}
}
int comp (const void * a, const void * b)
{
    return ( *(int*)a - *(int*)b );
}

#ifdef	ANCHOR_HASH_ALI

int compare_read_hash(const void * a, const void * b)
{
    read_h* rh1 = (read_h* )a;
    read_h* rh2 = (read_h* )b;

    if(rh1->des < rh2->des)
        return -1;
    else if(rh1->des > rh2->des)
        return 1;
    else    return 0;

}

int comepare_anchor_seed(const void * a, const void * b)
{
	anchor_seed* as1 = (anchor_seed* )a;
	anchor_seed* as2 = (anchor_seed* )b;
	
	if(as1->seed_length > as2->seed_length)
		return -1;
	else if(as1->seed_length < as2->seed_length)
		return 1;
	else	return 0;
}
#endif

#ifdef	ALTER_DEBUG
int compare_seed_length(const void * a, const void * b)
{
	seed_length_array* sa1 = (seed_length_array* )a;
	seed_length_array* sa2 = (seed_length_array* )b;
	
	if(sa1->seed_length > sa2->seed_length)
		return -1;
	else if(sa1->seed_length < sa2->seed_length)
		return 1;
	else	return 0;

}
#endif

#endif /* BINARY_QSORT_H_ */

