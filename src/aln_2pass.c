/*************************************************************************
	> File Name: aln_2pass.c
	> Author: 
	> Mail: 
 ************************************************************************/

#include <stdio.h>
#include <inttypes.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <pthread.h>

#include "binarys_qsort.h"
#include "bseq.h"
#include "aln_2pass.h"
#include "splic_junction.h"
#include "load_unipath_size.h"
#include "bit_operation.h"
#include "hash_index.h"
#include "format.h"
#include "read_seeding.h"
#include "ktime.h"

// #define DEBUG
// #define PRINT
int THREAD_READ_I;
pthread_rwlock_t RWLOCK;

Ival_Anno_t *IvalA = NULL;

int find_chr_n_by_name(char *chr_name)
{
	int i;
	for (i = 1; i < chr_file_n; ++i)
	{
		if (strcmp(chr_name, chr_names[i]) == 0)
			return i;
	}
	return 0;
}

void get_refseq(uint8_t *ref, uint32_t len, uint32_t start)
{
	uint32_t m;

    for (m = 0; m < len; ++m) 
    {
        ref[m] = (buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3;
    }
}

void get_ref_jump(uint8_t *ref, uint32_t start1, uint32_t len1, uint32_t start2, uint32_t len2)
{
	uint32_t m;
	uint32_t k = 0;

    for (m = 0; m < len1; ++m) 
    {
        ref[k++] = (buffer_ref_seq[(m + start1) >> 5] >> ((31 - ((m + start1) & 0X1f)) << 1)) & 0X3;
    }

    for (m = 0; m < len2; ++m) 
    {
        ref[k++] = (buffer_ref_seq[(m + start2) >> 5] >> ((31 - ((m + start2) & 0X1f)) << 1)) & 0X3;
    }
}

void get_ref_onebyone(uint8_t *ref, uint32_t start, uint32_t len, uint32_t pre_pos)
{
	uint32_t m;
	uint32_t k = pre_pos;
    for (m = 0; m < len; ++m) 
    {
        ref[k++] = (buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3;
    }
}

void get_junc(int chr_n, uint32_t start, uint32_t l, uint8_t *junc, uint8_t flag)
{
    if (!flag)
        return ;
    //uint32_t Sta = chr_end_n[chr_n - 1] - 1;
    int32_t end = start + l;
    int32_t i;
    int32_t low, high, mid;

    memset(junc, 0, end - start);
    
    Ival_Anno_t *IA;
    low = 0;
    IA = &IvalA[chr_n];
    high = IA->Intron_n;
    while(high > low)
    {
        mid = low + ((high - low) >> 1);
        if (IA->Ival[mid].Is >= start) high = mid;
        else low = mid + 1;
    }

    for(i = low; i < IA->Intron_n; ++i)
    {
        if (start <= IA->Ival[i].Is && end >= IA->Ival[i].Ie && IA->Ival[i].strand != 2)
        {
            if (IA->Ival[i].strand == 0)
            {
                junc[IA->Ival[i].Is - start] |= 1;
                junc[IA->Ival[i].Ie - 1 - start] |= 2;
            }
            else
            {
                junc[IA->Ival[i].Is - start] |= 8;
                junc[IA->Ival[i].Ie - 1 - start] |= 4;
            }
        }
    }
    
}


uint32_t apply_exon(param_map *opt, TARGET_t *anchor_map2ref, uint32_t map2ref_cnt, Anno_t *annotation, uint32_t total_items);

int load_anchor(TARGET_t *anchor_map2ref, uint32_t map2ref_cnt);
Ival_Anno_t *read_Annotation(param_map *opt, TARGET_t *anchor_map2ref, uint32_t map2ref_cnt, uint32_t *merge_cnt)
{
    FILE* fp_anno = fopen(opt->anno_path, "r");
	if (fp_anno == NULL)
	{
		fprintf(stderr, "[Wrong!!!] Open pre-processed annotation file %s wrong.\n", opt->anno_path);
		exit(0);
	}

	uint32_t TC_exon, TC_intron;
    fscanf(fp_anno, "%u\t%u\n", &TC_exon, &TC_intron);
    Ival_Anno_t *InterVal = (Ival_Anno_t *)calloc(chr_file_n, sizeof(Ival_Anno_t));

    int i;
    for (i = 0; i < chr_file_n; ++i)
    {
        InterVal[i].Intron_n = 0;
        InterVal[i].Ival = (Ival_anno_t *)calloc(TC_intron, sizeof(Ival_anno_t));
    }
    Anno_t *annotation = (Anno_t *)calloc(TC_exon * 2, sizeof(Anno_t));
    uint32_t idx_e = 0;
    char LINE[65535];
    char *p, *bs;
    int chr_n;
    int n_I;
    int strand;
    uint32_t Sta = 0;
    while( fgets(LINE, 65535, fp_anno) != NULL)
    {
        p = strtok(LINE, "\t");
        i = 0;
        while(p)
        {
            if(i == 0) { //seqname
                chr_n = find_chr_n_by_name(p);
                Sta = chr_end_n[chr_n - 1] - 1;
            } 
            else if (i == 2) 
            {
                if(*p == '+')
                    strand = 0;
                else if(*p == '-')
                    strand = 1;
                else
                    strand = 2;
            } 
            else if (i == 3) 
            {
                n_I = atol(p);
            } 
            else if (i == 4)
            {
                int Is, Ie;
                bs = p;
                int j;
                Is = strtol(bs, &bs, 10) + Sta; ++bs;
                Ie = strtol(bs, &bs, 10) + Sta; ++bs;
                annotation[idx_e].start = Is;
                annotation[idx_e].end = Ie;
                idx_e++;
                int st = Ie;
                for(j = 1; j < n_I; ++j)
                {
                    Is = strtol(bs, &bs, 10) + Sta; ++bs;
                    Ie = strtol(bs, &bs, 10) + Sta; ++bs;
                    annotation[idx_e].start = Is;
                    annotation[idx_e].end = Ie;
                    annotation[idx_e].strand = strand;

                    InterVal[chr_n].Ival[InterVal[chr_n].Intron_n].Is = st;
                    InterVal[chr_n].Ival[InterVal[chr_n].Intron_n].Ie = Is - 1;

                    InterVal[chr_n].Ival[InterVal[chr_n].Intron_n].strand = strand;
                    InterVal[chr_n].Intron_n += 1;

                    st = Ie;
                    idx_e++;
                }
            }
            
            p = strtok(NULL, "\t");
            i++;
        }
    }

	
    int tI1 = 0;
    int tI2 = 0;
    for(i = 0; i < chr_file_n; ++i)
    {
        tI1 += InterVal[i].Intron_n;
	    qsort(InterVal[i].Ival, InterVal[i].Intron_n, sizeof(Ival_anno_t), compare_intron);
        //remove redundant
        int j = 0;
        int nN = 1;
        if (InterVal[i].Intron_n > 2)
        {
            uint32_t ts, te;
            int strand;
            ts = InterVal[i].Ival[j].Is;
            te = InterVal[i].Ival[j].Ie;
            strand = InterVal[i].Ival[j].strand;
            j++;
            while (j < InterVal[i].Intron_n)
            {
                if ((InterVal[i].Ival[j].Is == ts) && (InterVal[i].Ival[j].Ie == te) && (InterVal[i].Ival[j].strand == strand))
                {
                    j++;
                    continue;
                }
                ts = InterVal[i].Ival[j].Is;
                te = InterVal[i].Ival[j].Ie;
                strand = InterVal[i].Ival[j].strand;
                InterVal[i].Ival[nN].Is = ts;
                InterVal[i].Ival[nN].Ie = te;
                InterVal[i].Ival[nN].strand = strand;
                
                j++;
                nN++;
            }
            InterVal[i].Intron_n = nN;
            tI2 += nN;
        }
    }
    qsort(annotation, idx_e, sizeof(Anno_t), compare_exon);
    //remove redunant
    if (idx_e > 2)
    {
        uint32_t ts, te;
        int strand;
        i = 1;
        ts = annotation[0].start;
        te = annotation[0].end;
        strand = annotation[0].strand;
        int nE = 1;
        while(i < idx_e)
        {
            if ((annotation[i].start == ts) && (annotation[i].end == te) && (annotation[i].strand == strand))
            {
                i++;
                continue;
            }
            ts = annotation[i].start;
            te = annotation[i].end;
            strand = annotation[i].strand;
            annotation[nE].start = ts;
            annotation[nE].end = te;
            annotation[nE].strand = strand;
            i++;
            nE++;
        }
        idx_e = nE;
    }

   
    *merge_cnt = apply_exon(opt, anchor_map2ref, map2ref_cnt, annotation, idx_e); //m2

    fclose(fp_anno);
    if (annotation != NULL) free(annotation);

    return InterVal;
}

static void cal_overlap(uint32_t anchor_s, uint32_t anchor_e, uint32_t anno_s, uint32_t anno_e, int *overlap_num)
{
	uint32_t small = (anchor_s < anno_s)? anno_s : anchor_s;
	uint32_t large = (anchor_e < anno_e)? anchor_e : anno_e;

	if (small >= large)
	{
		*overlap_num = 0;
	}
	else
	{
		*overlap_num = large - small + 1;
	}
}

static int have_overlap(uint32_t anchor_s, uint32_t anchor_e, uint32_t anno_s, uint32_t anno_e)
{
	uint32_t small = (anchor_s < anno_s)? anno_s : anchor_s;
	uint32_t large = (anchor_e < anno_e)? anchor_e : anno_e;

	if (small >= large)
		return 0;
	
	return 1;
}

static int find_exact_match(Anno_t *annotations, uint32_t upper, uint32_t down, uint32_t ts, uint32_t te)
{
	int j;
	int thre = 5;
	int distance = 0;
	int distance_min= 0xffff;
	int idx = -1;
	for (j = upper; j <= down; ++j)
	{
		//if ((annotations[j].start - thre > ts) || (ts > annotations[j].start + thre) || (annotations[j].end - thre > te) || (te > annotations[j].end + thre))
		if ((annotations[j].start - thre > ts) || (te > annotations[j].end + thre))			
			continue;
		else
		{
			distance = abs((int)(ts - annotations[j].start)) + abs((int)(te - annotations[j].end));
			if (distance < distance_min)
			{
				distance_min = distance;
				idx = j;
			}

		}
	}
	return idx;
}

static int load_anchor_with_gtf(TARGET_t *anchor_map2ref, Anno_t *annotations, uint32_t map2ref_cnt, uint32_t total_items, Anno_range_t *anno_range, uint32_t total_range)
{	
	uint32_t i, j, m;

	FILE* fp_temp = fopen(temp_binary_pos, "rb");

	rewind(fp_temp);
	fread(anchor_map2ref, sizeof(TARGET_t), map2ref_cnt, fp_temp);
    fclose(fp_temp);

	//sort anchor_map2ref according ts
	qsort(anchor_map2ref, map2ref_cnt, sizeof(TARGET_t), compare_anchor);

	int merge_cnt = 0;
	uint32_t ts;
	uint32_t te;
	uint32_t coverage = 1;
    int max_dis_connect = 20;
	int filter = seed_k_t + 10;
	int filter2 = seed_k_t << 1;

	ts = anchor_map2ref[0].ts;
	te = anchor_map2ref[0].te;
	i = 1;

	//merge
	while (i < map2ref_cnt)
	{
		if (anchor_map2ref[i].te <= te) //contained
		{
			coverage++;
			i++;
			continue;
		}
		else if(anchor_map2ref[i].ts <= te)
		{
			te = anchor_map2ref[i].te;
			coverage++;
			i++;
		}
		else
		{	
			if (coverage > 1 || (te - ts > filter))
			{
				anchor_map2ref[merge_cnt].ts = ts;
				anchor_map2ref[merge_cnt].te = te;
				anchor_map2ref[merge_cnt].cov = coverage;
				anchor_map2ref[merge_cnt].strand = 3; // uncertain which strand this anchor belone to
				merge_cnt++;
			}
			ts = anchor_map2ref[i].ts;
			te = anchor_map2ref[i].te;
			coverage = 1;
			i++;	
		}
	}
	//the last
	if(coverage > 1 || (te - ts > filter)) //filter anchor with low coverage
	{
		anchor_map2ref[merge_cnt].ts = ts;
		anchor_map2ref[merge_cnt].te = te;
		anchor_map2ref[merge_cnt].cov = coverage;
		anchor_map2ref[merge_cnt].strand = 3;
		merge_cnt++;
	}

	//merge, if two anchor less than 10bp, we think the two anchor belone to the same anchor ......................
	m = merge_cnt;
	merge_cnt = 0;
	ts = anchor_map2ref[0].ts;
	te = anchor_map2ref[0].te;
	coverage = anchor_map2ref[0].cov;
	for (i = 1; i < m; ++i)
	{
		if (anchor_map2ref[i].ts < te + max_dis_connect)
		{
			te = anchor_map2ref[i].te;
			coverage += anchor_map2ref[i].cov;
		}
		else
		{
			anchor_map2ref[merge_cnt].ts = ts;
			anchor_map2ref[merge_cnt].te = te;
			anchor_map2ref[merge_cnt].cov = coverage;
			anchor_map2ref[merge_cnt].strand = 3;
			merge_cnt++;

			ts = anchor_map2ref[i].ts;
			te = anchor_map2ref[i].te;
			coverage = anchor_map2ref[i].cov;
		}
	}
	anchor_map2ref[merge_cnt].ts = ts;
	anchor_map2ref[merge_cnt].te = te;
	anchor_map2ref[merge_cnt].cov = coverage;
	anchor_map2ref[merge_cnt].strand = 3;
	merge_cnt++;

	//find corresponding item in annotation by order loop
	uint32_t final_merge_cnt = 0;
    int non_exon_region_but_good = 0;
    int exon_region_but_not_complete = 0;
    int wrong_filter1 = 0;
    int wrong_filter2 = 0;
	uint32_t start_idx = 0;
    int max_overlap_num = 0;
    int max_overlap_id = -1;
    int overlap_num = 0;
	j = 0;
	/*new truck*/
	for (i = 0; i < merge_cnt; ++i)
	{
		start_idx = j;
		while (start_idx < total_range && anchor_map2ref[i].ts > anno_range[start_idx].end)
		{
			start_idx++;
		}
        cal_overlap(anchor_map2ref[i].ts, anchor_map2ref[i].te, anno_range[start_idx].start, anno_range[start_idx].end, &overlap_num);
        if (overlap_num > 0)
		{
			//have overlap
            max_overlap_id = start_idx;
            max_overlap_num = overlap_num;
			for (j = start_idx + 1; j < total_range; ++j)
			{
                cal_overlap(anchor_map2ref[i].ts, anchor_map2ref[i].te, anno_range[j].start, anno_range[j].end, &overlap_num);
                if (overlap_num == 0)
                {
                    break;
                }
                else if(max_overlap_num < overlap_num)
                {
                    max_overlap_num = overlap_num;
                    max_overlap_id = j;
                }
			}
			//query correspond items from each range trunck
			int r_max;
            r_max = find_exact_match(annotations, anno_range[max_overlap_id].upper, anno_range[max_overlap_id].down, anchor_map2ref[i].ts, anchor_map2ref[i].te);

			if (r_max >= 0)
			{
				anchor_map2ref[final_merge_cnt].ts = annotations[r_max].start;
				anchor_map2ref[final_merge_cnt].te = annotations[r_max].end;
				anchor_map2ref[final_merge_cnt].strand = annotations[r_max].strand;
				anchor_map2ref[final_merge_cnt].cov = anchor_map2ref[i].cov;
				final_merge_cnt ++;
			}
			else
			{
				if ((anchor_map2ref[i].te - anchor_map2ref[i].ts > filter2) || (anchor_map2ref[i].cov > 1)) //if satifiy, record also
				{
					anchor_map2ref[final_merge_cnt].ts = anchor_map2ref[i].ts;
					anchor_map2ref[final_merge_cnt].te = anchor_map2ref[i].te;
					anchor_map2ref[final_merge_cnt].cov = anchor_map2ref[i].cov;
					anchor_map2ref[final_merge_cnt].strand = 3; //  later will process strand == 3
					final_merge_cnt ++;
                    exon_region_but_not_complete += 1;
				}
                else
                {
                    wrong_filter1 += 1;
                }
			}	
			j = start_idx;
		}
		else
		{
			if ((anchor_map2ref[i].te - anchor_map2ref[i].ts > filter2) || (anchor_map2ref[i].cov > 1)) //if satifiy, record also
			{
				anchor_map2ref[final_merge_cnt].ts = anchor_map2ref[i].ts;
				anchor_map2ref[final_merge_cnt].te = anchor_map2ref[i].te;
				anchor_map2ref[final_merge_cnt].cov = anchor_map2ref[i].cov;
				anchor_map2ref[final_merge_cnt].strand = 3; //
				final_merge_cnt ++;
                non_exon_region_but_good += 1;
			}
            else
            {
                wrong_filter2 += 1;
            }
			j = start_idx;
		}
	}

    //re-merge the anchor_map2ref
	qsort(anchor_map2ref, final_merge_cnt, sizeof(TARGET_t), compare_anchor);
    m = final_merge_cnt;
    final_merge_cnt = 0;
    ts = anchor_map2ref[0].ts;
    te = anchor_map2ref[0].te;
    int strand = anchor_map2ref[0].strand;
    int ch = 0;
    int change = 0;
    for (i = 1; i < m; ++i)
    {
        if (anchor_map2ref[i].ts <= te + 20)
		{
			ts = (ts < anchor_map2ref[i].ts)? ts : anchor_map2ref[i].ts;
			te = (te > anchor_map2ref[i].te)? te : anchor_map2ref[i].te;
            if (anchor_map2ref[i].strand != strand)
            {
                ch = 1;
            }
		}
		else
		{
			anchor_map2ref[final_merge_cnt].ts = ts;
			anchor_map2ref[final_merge_cnt].te = te;
            anchor_map2ref[final_merge_cnt].strand = (ch == 1)? 3 : strand;
			final_merge_cnt++;
            if (ch == 1)
                change += 1;
			ts = anchor_map2ref[i].ts;
			te = anchor_map2ref[i].te;
            strand = anchor_map2ref[i].strand;
            ch = 0;
		}
    }
    anchor_map2ref[final_merge_cnt].ts = ts;
    anchor_map2ref[final_merge_cnt].te = te;
    anchor_map2ref[final_merge_cnt].strand = strand;
    final_merge_cnt++;

	return final_merge_cnt;
}

uint32_t apply_exon(param_map *opt, TARGET_t *anchor_map2ref, uint32_t map2ref_cnt, Anno_t *annotation, uint32_t total_items)
{
    //merge overlapped exon which two boundary both different as new exon
    uint32_t ts = annotation[0].start;
	uint32_t te = annotation[0].end;
    int strand = annotation[0].strand;
    uint32_t total_items_new = total_items;

    int j = 0;
    int sig = 0;
    int i = 1;
    while (i < total_items)
    {
        j = i;
        sig = 0;
        while ((ts <= annotation[j].start && te >= annotation[j].start) && (j < total_items))
        //while ((ts < annotation[j].start) && (te > annotation[j].start) && (te < annotation[j].end) && (strand == annotation[j].strand) && (i < total_items))
        {
            if ((ts < annotation[j].start) && (te < annotation[j].end) && (strand == annotation[j].strand))
            {    
                te = annotation[j].end;
                sig = 1;
            }
            j++;
        }
        if (sig)
        {
            annotation[total_items_new].strand = strand;
            annotation[total_items_new].start = ts;
            annotation[total_items_new].end = te;
            total_items_new ++;    
        }
        
        ts = annotation[j].start;
        te = annotation[j].end;
        strand = annotation[j].strand;
        i = j + 1;
    }

    total_items = total_items_new;
    qsort(annotation, total_items_new, sizeof(Anno_t), compare_exon);

	//get trunk of annotations
	Anno_range_t *anno_range = (Anno_range_t *)calloc(total_items_new, sizeof(Anno_range_t));
	uint32_t total_range = 0;
	ts = annotation[0].start;
	te = annotation[0].end;
	uint32_t upper = 0;
	for (i = 1; i < total_items_new; ++i)
	{
		if (have_overlap(ts, te, annotation[i].start, annotation[i].end))
		{
			ts = (ts < annotation[i].start)? ts : annotation[i].start;
			te = (te > annotation[i].end)? te : annotation[i].end;
		}
		else
		{
			anno_range[total_range].start = ts;
			anno_range[total_range].end = te;
			anno_range[total_range].upper = upper;
			anno_range[total_range].down = i - 1;
			total_range++;
			ts = annotation[i].start;
			te = annotation[i].end;
			upper = i;
		}
	}
	anno_range[total_range].start = ts;
	anno_range[total_range].end = te;
	anno_range[total_range].upper = upper;
	anno_range[total_range].down = i - 1;
	total_range++;

	//load anchor an compare to annotation
	uint32_t final_merge_cnt = load_anchor_with_gtf(anchor_map2ref, annotation, map2ref_cnt, total_items_new, anno_range, total_range);

	if (anno_range != NULL)	free(anno_range);
    
	return final_merge_cnt;
}

int load_anchor(TARGET_t *anchor_map2ref, uint32_t map2ref_cnt)
{
	uint32_t i, m;
	FILE* fp_temp = fopen(temp_binary_pos, "rb");

	//read anchor info from fp_temp
	rewind(fp_temp);
	fread(anchor_map2ref, sizeof(TARGET_t), map2ref_cnt, fp_temp);

	//sort anchor_map2ref according ts
	qsort(anchor_map2ref, map2ref_cnt, sizeof(TARGET_t), compare_anchor);

	int merge_cnt = 0;
	uint32_t ts;
	uint32_t te;
	uint32_t coverage = 1;
	int max_dis_connect = 30;
	int filter = seed_k_t << 1;
    //int filter = 20;

	ts = anchor_map2ref[0].ts;
	te = anchor_map2ref[0].te;
	i = 1;
	while (i < map2ref_cnt)
	{
		if (anchor_map2ref[i].te <= te) //contained
		{
			coverage++;
			i++;
			continue;
		}
        else if (anchor_map2ref[i].ts < te + max_dis_connect)
		{
			te = anchor_map2ref[i].te;
			coverage++;
			i++;
		}
		else
		{	
			if (coverage > 1 || (te - ts > filter))
			{
				anchor_map2ref[merge_cnt].ts = ts;
				anchor_map2ref[merge_cnt].te = te;
				anchor_map2ref[merge_cnt].cov = coverage;
				anchor_map2ref[merge_cnt].strand = 3; // uncertain which strand this anchor belone to
				merge_cnt++;
			}
			ts = anchor_map2ref[i].ts;
			te = anchor_map2ref[i].te;
			coverage = 1;
			i++;	
		}
	}
	//the last
	if(coverage > 1 || (te - ts > filter)) //filter anchor with low coverage
	{
		anchor_map2ref[merge_cnt].ts = ts;
		anchor_map2ref[merge_cnt].te = te;
		anchor_map2ref[merge_cnt].cov = coverage;
		anchor_map2ref[merge_cnt].strand = 3;
		merge_cnt++;
	}

	fclose(fp_temp);

	return merge_cnt;
}


int chromosome_judge(uint32_t position, uint32_t *chromosome_begin)
{
	int file_n = 0;
	
	int low = 0;
	int high = chr_file_n - 1;
	int mid;

	while ( low <= high )
	{
		mid = (low + high) >> 1;
		if(position < (chr_end_n[mid]))
		{
			high = mid - 1;
		}
		else if(position > (chr_end_n[mid]))
		{
			low = mid + 1;
		}
		else
		{
			return mid;
		}
		file_n = low;
	}
	*chromosome_begin = chr_end_n[file_n - 1]; //the reference index start from pos 0
	return file_n;
}

static inline int binarysearch_anchor(TARGET_t *anchor_map2ref, int s1, int s2, uint32_t coordinate1, uint32_t coordinate2)
{
	int low, high, mid;
	uint32_t templ, tempr;
	// int signal = -1;

	low = s1;
	high = s2;
	
	while(low <= high) //<=
	{
		mid = (low + high) >> 1;
		templ = anchor_map2ref[mid].ts;
		tempr = anchor_map2ref[mid].te;
		if(coordinate2 < templ)
		{
			high = mid - 1;
		}
		else if(coordinate1 >= tempr)
		{
				low = mid + 1;
		}
		else  //found match
		{
			return mid;
		}
	}
	return -1;
}

void fix_cigar(_aln_t *aln, const uint8_t *qseq, const uint8_t *tseq, uint32_t qlen, uint32_t tlen, int *qshift, int *tshift, uint32_t *qs, uint32_t *ts)
{
	int32_t k, toff = 0, qoff = 0, to_shrink = 0;
	*qshift = *tshift = 0;
	if (aln->n_cigar <= 1) return;
	for (k = 0; k < aln->n_cigar; ++k) { // indel left alignment
		uint32_t op = aln->cigar[k]&0xf, len = aln->cigar[k]>>4;
		if (len == 0) to_shrink = 1;
		if (op == 0) {
			toff += len, qoff += len;
		} else if (op == 1 || op == 2) { // insertion or deletion
			if (k > 0 && k < aln->n_cigar - 1 && (aln->cigar[k-1]&0xf) == 0 && (aln->cigar[k+1]&0xf) == 0) {
				int l, prev_len = aln->cigar[k-1] >> 4;
				if (op == 1) {
					for (l = 0; l < prev_len; ++l)
						if (qseq[qoff - 1 - l] != qseq[qoff + len - 1 - l])
							break;
				} else {
					for (l = 0; l < prev_len; ++l)
						if (tseq[toff - 1 - l] != tseq[toff + len - 1 - l])
							break;
				}
				if (l > 0)
					aln->cigar[k-1] -= l<<4, aln->cigar[k+1] += l<<4, qoff -= l, toff -= l;
				if (l == prev_len) to_shrink = 1;
			}
			if (op == 1) qoff += len;
			else toff += len;
		} else if (op == 3) {
			toff += len;
		}
	}

    //later give up assert function, memory consumption
	assert(qoff == qlen && toff == tlen);
	if (to_shrink) { // squeeze out zero-length operations
		int32_t l = 0;
		for (k = 0; k < aln->n_cigar; ++k) // squeeze out zero-length operations
			if (aln->cigar[k]>>4 != 0)
				aln->cigar[l++] = aln->cigar[k];
		aln->n_cigar = l;
		for (k = l = 0; k < aln->n_cigar; ++k) // merge two adjacent operations if they are the same
			if (k == aln->n_cigar - 1 || (aln->cigar[k]&0xf) != (aln->cigar[k+1]&0xf))
				aln->cigar[l++] = aln->cigar[k];
			else aln->cigar[k+1] += aln->cigar[k]>>4<<4; // add length to the next CIGAR operator
		aln->n_cigar = l;
	}

	if ((aln->cigar[0]&0xf) == 1 || (aln->cigar[0]&0xf) == 2) { // get rid of leading I or D
		int32_t l = aln->cigar[0] >> 4;
		if ((aln->cigar[0]&0xf) == 1) {
			*qs += l;
			*qshift = l;
		} else *ts += l, *tshift = l;
		--(aln->n_cigar);
		memmove(aln->cigar, aln->cigar + 1, aln->n_cigar * 4);
	}
}
void mm_update_extra(_aln_t *aln, uint8_t *qseq, uint8_t *tseq, uint32_t qlen, uint32_t tlen, uint32_t *qs, uint32_t *ts, uint8_t q, uint8_t e)
{
    uint32_t k, l, toff = 0, qoff = 0;
	int32_t s = 0, max = 0;

	if (aln->n_cigar == 0) return;

	int32_t qshift, tshift;
	fix_cigar(aln, qseq, tseq, qlen, tlen, &qshift, &tshift, qs, ts);
	qlen -= qshift; // qseq and tseq may be shifted due to the removal of leading I/D
	tlen -= tshift;
	qseq += qshift;
	tseq += tshift;

	uint32_t *cigar = aln->cigar;
	uint32_t n_cigar = aln->n_cigar;
    for (k = 0; k < n_cigar; ++k) {
        uint32_t op = cigar[k]&0xf, len = cigar[k]>>4;
        if (op == 0) 
        { // match/mismatch
            int n_ambi = 0, n_diff = 0;
            for (l = 0; l < len; ++l) 
            {
                int cq = qseq[qoff + l], ct = tseq[toff + l];
                if (ct > 3 || cq > 3) ++n_ambi; //N
                else if (ct != cq) ++n_diff; //mismatch
				s += mata_R[ct * 5 + cq];
				if (s < 0) s = 0;
				else max = max > s? max : s;
            }
            aln->blen += len - n_ambi;
            aln->mlen += len - (n_ambi + n_diff); 
            aln->n_ambi += n_ambi;

            toff += len, qoff += len;
        } else if (op == 1) { // insertion
            int n_ambi = 0;
            for (l = 0; l < len; ++l)
                if (qseq[qoff + l] > 3) ++n_ambi;
            aln->blen += len - n_ambi, aln->n_ambi += n_ambi;

			s -= q + e * len;
			if (s < 0) s = 0;
            qoff += len;
        } else if (op == 2) { // deletion
            int n_ambi = 0;
            for (l = 0; l < len; ++l)
                if (tseq[toff + l] > 3) ++n_ambi;
            aln->blen += len - n_ambi, aln->n_ambi += n_ambi;
			s -= q + e * len;
			if (s < 0) s = 0;
            toff += len;
        } else if (op == 3) { // intron
            toff += len;
        }
    }
	aln->dp_max = max;
    assert(qoff == qlen && toff == tlen);
}

void mm_append_cigar(_aln_t *aln, uint32_t n_cigar, uint32_t *cigar) // TODO: this calls the libc realloc()
{
	if (n_cigar == 0)
		return;

    if((aln->n_cigar > 0) && ((aln->cigar[aln->n_cigar - 1]&0xf) == (cigar[0]&0xf)))
    {
        aln->cigar[aln->n_cigar - 1] += cigar[0]>>4<<4;
        memcpy(aln->cigar + aln->n_cigar, cigar + 1, (n_cigar - 1)*4);
        aln->n_cigar += n_cigar - 1;
    }else
    {
        memcpy(aln->cigar + aln->n_cigar, cigar, n_cigar*4);
        aln->n_cigar += n_cigar;
    }
}

uint32_t append_intron_to_cigar(void *km, ksw_extz_t *ez, uint32_t pre_pos, uint32_t intron_pos, uint32_t intron_len)
{
	int32_t i, j;
	uint8_t op;
	uint32_t len;
	uint32_t toff = 0;
	uint32_t cigar_lt = 0;

	if ((ez->n_cigar + 3) > ez->m_cigar) 
	{
		ez->m_cigar = (ez->m_cigar)<<1;
		ez->cigar = (uint32_t*)krealloc(km, ez->cigar, (ez->m_cigar) << 2);
	}

	uint32_t *cigar = &ez->cigar[pre_pos];
	uint32_t cigar_len = ez->n_cigar - pre_pos;
	if(intron_pos == 0)
	{
		for (i = cigar_len - 1; i >= 0; --i)
		{
			cigar[i + 1] = cigar[i];
		}
		cigar[0] = intron_len<<4 | 3; //N
		ez->n_cigar += 1;
		cigar_lt = 1;
	}
	else
	{
		for (i = 0; i < cigar_len; ++i)
		{
			op = cigar[i]&0xf;
			if ((op == 0) || (op == 2))
			{
				len = cigar[i]>>4;
				toff += len;
				if (toff > intron_pos)
				{
					//find intron position
					uint32_t len2 = toff - intron_pos;
					uint32_t len1 = len - len2;
					for (j = cigar_len - 1; j > i; --j)
					{
						cigar[j + 2] = cigar[j];
					}
					cigar[i++] = len1<<4 | op;
					cigar[i++] = intron_len<<4 | 3; //N
					cigar_lt = i;
					cigar[i++] = len2<<4 | op;

					ez->n_cigar += 2;
					break;
				}else if (toff == intron_pos)
				{
					for (j = cigar_len - 1; j > i; --j)
					{
						cigar[j + 1] = cigar[j];
					}
					cigar[++i] = intron_len<<4 | 3; //N
					cigar_lt = i + 1;

					ez->n_cigar += 1;
					break;
				}
			}
		}
	}

	return cigar_lt + pre_pos;
}

static inline void mm_seq_rev(uint32_t len, uint8_t *seq)
{
    uint32_t i;
    uint8_t t;
    for (i = 0; i < len>>1; ++i)
        t = seq[i], seq[i] = seq[len - 1 - i], seq[len - 1 - i] = t;
}

static void reset_aln_t (_aln_t *aln)
{
	uint8_t i;
	for(i = 0; i < 2; i++)
	{
		aln[i].n_cigar = 0;
		aln[i].mlen = 0;
		aln[i].blen = 0;
		aln[i].n_ambi = 0;
		aln[i].dp_score = 0;
		aln[i].dp_max = 0;
		aln[i].dp_max2 = 0;
		aln[i].flag = 0;
	}
}

static int compare_idx(const void *a, const void *b)
{
	return *(int *)a - *(int *)b;
}

uint32_t local_hash_anchor(uint8_t *qseq, uint32_t qlen, uint32_t *idx_cnt_array, TARGET_t *anchor_map2ref, uint32_t key1, uint32_t key2, param_map *opt, uint8_t type, uint8_t tid, uint8_t signal)
{

	if((key1 >= key2) || (qlen < hash_kmer))
	{
		return 0;
	}

	uint32_t i, j;
	uint32_t pre_pos = 0;
	uint32_t s_s;
	uint32_t s_len;
	uint32_t n = key2 - key1;

	int* cov_score = (int* )calloc(n, sizeof(float));

	local_hash_process(qseq, qlen, cov_score, anchor_map2ref, key1, key2, tid, signal);

#ifdef PRINT
    for ( i = 0; i < key2 - key1; ++i)
    {
        fprintf(stderr, "cov_score[%d] = %d\n", i, cov_score[i]);
    }
#endif

	int real_cnt = 0;
	for (i = 0; i < n; ++i)
	{
		if (cov_score[i] == 1)
		{
			idx_cnt_array[real_cnt] = i;
			real_cnt++;
		}
	}

	
	free(cov_score);

	return real_cnt;
}


void check_more_part(uint32_t *cigar, uint32_t *n_cigar, int *score, int q2)
{
	if (*n_cigar < 4)
	{
		return;
	}
	int i,j;
	uint8_t op, op1;
	uint32_t op_len, op_len1;
	uint32_t len_ref = 0;
	uint32_t len_read = 0;

	uint32_t f_l, f_l_r;
	int left, right;
	uint32_t m_cigar = *n_cigar;
	i = 0;
	while (i < m_cigar)
	{
		op = cigar[i]&0xf;
		op_len = (cigar[i]>>4);
		if (op == 2 && op_len > 20) //if have deletion long than 5bp
		{
			len_ref = 0;
			len_read = 0;
			f_l = op_len;
			f_l_r = 0;
			left = right = i;
			//up
			j = (i == 0)? -1 : (i - 1);
			while(j >= 0)
			{
				op1 = cigar[j]&0xf;
				op_len1 = (cigar[j]>>4);
				if (op1 == 3) //the up closest intron
				{
					if (len_read <= 5)
					{
						left = j;
						f_l += len_ref + op_len1;
						f_l_r += len_read;
					}	
					break;
				}
				else if (op1 == 0)
				{
					len_read += op_len1;
					len_ref += op_len1;
				}
				else if( op1 == 1)
				{
					len_read += op_len1;
				}
				else if (op1 == 2)
				{
					len_ref += op_len1;
				}
				j--;
			}
			//down
			j = (i == m_cigar - 1)? m_cigar : (i + 1);
			len_ref = 0;
			len_read = 0;
			while (j < m_cigar)
			{
				op1 = cigar[j]&0xf;
				op_len1 = (cigar[j]>>4);
				if (op1 == 3) //the up closest intron
				{
					if (len_read <= 5)
					{
						right = j;
						f_l += len_ref + op_len1;
						f_l_r += len_read;
					}	
					break;
				}
				else if (op1 == 0)
				{
					len_read += op_len1;
					len_ref += op_len1;
				}
				else if( op1 == 1)
				{
					len_read += op_len1;
				}
				else if (op1 == 2)
				{
					len_ref += op_len1;
				}
				j++;
			}
			if (right > left)
			{
				//refine intron part
				cigar[left++] = f_l<<4 | 3; //N
				if (f_l_r > 0)
					cigar[left++] = f_l_r<<4 | 1;
				i = left;
				while (right < m_cigar - 1)
				{
					cigar[left++] = cigar[++right];
				}
				*n_cigar = left;
				m_cigar = *n_cigar;
				*score += q2;
			}
		}
		++i;
	}
}

static uint32_t refine_site(uint32_t site, uint8_t strand, uint8_t type)
{
	uint32_t m;
	int range =10;
	int len = 2*range + 10;
	uint32_t start;
	int16_t s1 = 0;
	int16_t s2 = 0;
	uint32_t r_site = site;
	uint8_t *ref = (uint8_t* )calloc(len + 2, 4);

	if (type == 0) //donor detect
	{
		start = site - range - 5;
		get_refseq(ref, len, start);
		donor_signals_detected(ref, len, start, &s1, &s2, strand);
		if (strand == 0 && s1 != -1) //+
		{
			r_site = start + s1 - 1;
		}
		else if (strand == 1 && s2 != -1) //-
		{
			r_site = start + s2 - 1;
		}
	}
	else //acceptor detec
	{
		start = site - range - 4;
		get_refseq(ref, len, start);
		acceptor_signals_detected(ref, len, start, &s1, &s2, strand);

		if (strand == 0 && s1 != -1) //+
		{
			r_site = start + s1 + 1;
		}
		else if (strand == 1 && s2 != -1) //-
		{
			r_site = start + s2 + 1;
		}
	}

	free(ref);
	return r_site;
}

int check_realign(void *km, param_map *opt, int bandwith, uint8_t *qseq, uint32_t qlen, uint32_t tlen, ksw_extz_t *ez, int pre_score, uint32_t start_pos, int intron_cnt, int splice_flag)
{
	int i, j, m;
	uint8_t re_ali = 0;
    uint8_t re_ali1 = 0;
	uint8_t op, op1;
	uint32_t op_len, op_len1;

	uint32_t m_cigar = ez->n_cigar;
	uint32_t *cigar = ez->cigar;
	
	ksw_extz_t ez_tmp;
	memset(&ez_tmp, 0, sizeof(ksw_extz_t));

    int thre = 5;

	////cal the end pos
	for (i = 0; i < m_cigar; ++i)
	{
		op = cigar[i]&0xf;
		op_len = (cigar[i]>>4);
		// if ((op == 1 || op == 2) && op_len > thre)
        if (op == 2 && op_len > thre)
		{
			re_ali = 1;
			break;
		}
	}
	if (re_ali == 0)
		return 0;

	
	uint32_t len_ref = 0;
	uint32_t len_read = 0;
	uint32_t s_s = start_pos;
	uint8_t intron_id = 0;
	uint8_t search_step = 5;
	int left, right;
	
	int startidx = 0;
	uint8_t strand = (splice_flag & MM_F_SPLICE_FOR)? 0 : 1;

	uint32_t **intron_len = (uint32_t **)calloc(intron_cnt, sizeof(uint32_t*));
	for (i = 0; i < intron_cnt; ++i)
	{
		intron_len[i] = (uint32_t* )calloc(2, 4);
	}

	
	uint8_t *ref = (uint8_t* )calloc(tlen + 50, 4);
	i = 0;
	while(i < m_cigar)
	// for (i = 0; i < m_cigar; ++i)
	{
		left = right = i;
		op = cigar[i]&0xf;
		op_len = (cigar[i]>>4);
		if (op == 3) //if intron
		{
			re_ali = 0;
			//search down by search_step to find larger deletion or insertions, if found, pass
			m = (i > m_cigar - search_step)? (m_cigar - 1) : (i + search_step);
			len_read = 0;
			for (j = i + 1; j <= m; ++j)
			{
				op1 = cigar[j]&0xf;
				op_len1 = (cigar[j]>>4);
                if (op1 == 2 && op_len1 > thre && len_read < 5) //the up closest insertion/deletion ; len_read < 5
				{
					re_ali = 1;
					break;
				}
				else if (op1 == 0 || op1 == 1)
				{
					len_read += op_len1;
				}
			}
			//record
			if (!re_ali)
			{
				for (m = startidx; m < i; ++m)
				{
					op1 = cigar[m]&0xf;
					if (op1 == 0 || op1 == 2)
						s_s += (cigar[m]>>4);
				}
				intron_len[intron_id][0] = s_s - 1;
				s_s += op_len;
				intron_len[intron_id][1] = s_s;
				intron_id++;
				startidx = i + 1;
			}
		}
		//if large D/I
        if (op == 2 && op_len > thre)
		{
			len_ref = 0;
			len_read = 0;
			//up
			j = (i == 0)? -1 : (i - 1);
			while(j >= 0)
			{
				op1 = cigar[j]&0xf;
				op_len1 = (cigar[j]>>4);
				if (op1 == 3) //the up closest intron
				{
					if (len_read <= 5)
					{
						left = j;
					}	
					break;
				}
				else if (op1 == 0 || op1 == 1)
				{
					len_read += op_len1;
				}
				j--;
				//break if condition
				if (len_read > 20)
					break;
			}
			//down
			j = (i == m_cigar - 1)? m_cigar : (i + 1);
			len_ref = 0;
			len_read = 0;
			while (j < m_cigar)
			{
				op1 = cigar[j]&0xf;
				op_len1 = (cigar[j]>>4);
				if (op1 == 3) //the up closest intron
				{
					if (len_read <= 5)
					{
						right = j;
					}	
					break;
				}
				else if (op1 == 0 || op1 == 1)
				{
					len_read += op_len1;
				}
				j++;

				if (len_read > 20)
					break;
			}
			if (left < right)
			{
				re_ali1 = 1;
				//refine intron part  >>up
				int m;
				len_ref = 0;
				uint32_t new_acceptor;
				uint32_t new_donor;
				for (m = 0; m < left; ++m)
				{
					op1 = cigar[m]&0xf;
					if (op1 == 2 || op1 == 0 || op1 == 3)
						len_ref += (cigar[m]>>4);
				}

				uint32_t tmp_start;

				if (left == i && op == 1) //the current is large 'I'
				{
					tmp_start = start_pos + len_ref + op_len;
				}
				else
				{
					tmp_start = start_pos + len_ref;
				}
				new_donor = refine_site(tmp_start, strand, 0); //donor);

				for (m = left; m <= right; ++m)
				{
					op1 = cigar[m]&0xf;
					if (op1 == 2 || op1 == 0 || op1 == 3)
						len_ref += (cigar[m]>>4);	
				}

				if (right == i && op == 1)
				{
					tmp_start = start_pos + len_ref - op_len;
				}
				else
				{
					tmp_start = start_pos + len_ref;
				}
				new_acceptor = refine_site(tmp_start, strand, 1); //donor

				if (new_acceptor <= new_donor)
				{
					if ((cigar[left]&0xf) == 3)
					{
						for (m = startidx; m < left; ++m)
						{
							op1 = cigar[m]&0xf;
							if (op1 == 0 || op1 == 2)
								s_s += (cigar[m]>>4);
						}
						intron_len[intron_id][0] = s_s - 1;
						s_s += (cigar[left]>>4);
						intron_len[intron_id][1] = s_s;
						intron_id++;
						// startidx = i + 1;
					}
					if ((cigar[right]& 0xf) == 3)
					{
						for (m = startidx; m < right; ++m)
						{
							op1 = cigar[m]&0xf;
							if (op1 == 0 || op1 == 2)
								s_s += (cigar[m]>>4);
						}
						intron_len[intron_id][0] = s_s - 1;
						s_s += (cigar[left]>>4);
						intron_len[intron_id][1] = s_s;
						intron_id++;
					}
				}
				else
				{
					intron_len[intron_id][0] = new_donor;
					intron_len[intron_id][1] = new_acceptor;
					intron_id++;
				}
				s_s = start_pos + len_ref;
				startidx = right + 1;
			}
		}
		i = right + 1;
	}
	//
	uint8_t new_intron_id = intron_id;
	if (intron_id > 1)
	{
		new_intron_id = 0;
		uint32_t s = intron_len[0][0];
		uint32_t e = intron_len[0][1];
		for (i = 1; i < intron_id; ++i)
		{
			if (intron_len[i][0] <= e)
			{
				e = intron_len[i][1];
			}
			else
			{
				if (e > s + Eindel)
				{
					intron_len[new_intron_id][0] = s;
					intron_len[new_intron_id][1] = e;
					new_intron_id ++;
					s = intron_len[i][0];
					e = intron_len[i][1];
				}
				else
				{
					s = intron_len[i][0];
					e = intron_len[i][1];
				}
			}
		}
        if (e > s + Eindel)
        {
            intron_len[new_intron_id][0] = s;
            intron_len[new_intron_id][1] = e;
            new_intron_id ++;
            s = intron_len[i][0];
            e = intron_len[i][1];
        }
	}
	else
	{
		if (intron_len[0][1] < intron_len[0][0] + Eindel)
			new_intron_id = 0;
	}

	if (re_ali1)
	{
		uint32_t pre_pos = 0;
		if (new_intron_id > 0)
		{
			s_s = start_pos;
			int l;
			for(i = 0; i < new_intron_id; ++i)
			{
				l = intron_len[i][0] - s_s + 1;
				if (l > 0)
				{
					get_ref_onebyone(ref, s_s, l, pre_pos);
					pre_pos += l;
				}
				else{
					intron_len[i][0] = s_s;
					get_ref_onebyone(ref, s_s, 1, pre_pos);
					pre_pos += 1;
				}
				s_s = intron_len[i][1];
			}
			//last
			l = tlen + start_pos - s_s ;
			if (l > 0)
			{
				get_ref_onebyone(ref, s_s, l, pre_pos);
				pre_pos += l;
			}
			else
			{
				intron_len[new_intron_id - 1][1] = tlen + start_pos;
			}
		}
		else
		{
			get_refseq(ref, tlen, start_pos);
			pre_pos = tlen;
		}
		

		align_non_splice(km, qseq, ref, qlen, pre_pos, opt, &ez_tmp, bandwith, 0, 2);
		if (ez_tmp.score < pre_score)
		{
			re_ali1 = 0;
			goto FREE;
		}

		pre_pos = 0;
        s_s = start_pos;
		for (i = 0; i < new_intron_id; ++i)
		{
			pre_pos = append_intron_to_cigar(km, &ez_tmp, pre_pos, intron_len[i][0] - start_pos + 1, intron_len[i][1] - intron_len[i][0] - 1);
            start_pos = intron_len[i][1];
		}
		//copy
		ez->score = ez_tmp.score;
		if ((ez->m_cigar) < ez_tmp.n_cigar) 
		{
			ez->m_cigar = (ez_tmp.n_cigar)<<1;
			ez->cigar = (uint32_t*)krealloc(km, ez->cigar, (ez->m_cigar) << 2);
		}
		ez->n_cigar = ez_tmp.n_cigar;
		for (i = 0; i < ez_tmp.n_cigar; ++i)
		{
			ez->cigar[i] = ez_tmp.cigar[i];
		}

	}

FREE:
    kfree(km, ez_tmp.cigar);
	free(ref);
	for (i = 0; i < intron_cnt; ++i)
		free(intron_len[i]);
	free(intron_len);
    if (re_ali1)
    	return 1;
    else
        return 0;
}

void check_cigar(uint8_t *qseq, uint8_t *tseq, uint32_t *cigar, uint32_t *n_cigar, uint32_t *qs, uint32_t *ts, int *score, uint32_t qlen, uint32_t tlen, uint32_t boundary, uint8_t type, param_map *opt)
{
	int i;
	uint8_t op;
	uint32_t op_len;
	uint32_t len_ref = 0;
	uint32_t len_read = 0;
	int score_cut = 0;
	int l;
	int cq, ct;

	uint32_t m_cigar = *n_cigar;
	int qoff = qlen - 1;
	int toff = tlen - 1;
	
	for(i = m_cigar - 1; i >= 0; --i)
	{
		op = cigar[i]&0xf;
		op_len = (cigar[i]>>4);
		if (op == 0) //M
		{
			len_ref += op_len;
			len_read += op_len;
			for (l = 0; l < op_len; ++l)
			{
				cq = qseq[qoff - l];
				ct = tseq[toff - l];
				if (cq == ct)
					score_cut -= opt->match_D;
				else
					score_cut += opt->mismatch_D;
			}
			qoff -= op_len;
			toff -= op_len;
		}
		else if (op == 1) //insertion or softclip
		{
			len_read += op_len;
			if ((i > 0) && (cigar[i-1]&0xf) == 1)
				score_cut += op_len * opt->gap_ex_D;
			else
				score_cut += opt->gap_open_D + op_len * opt->gap_ex_D;
			qoff -= op_len;
		}
		else if (op == 2) //deletion
		{
			len_ref += op_len;
			if ((i > 0) && (cigar[i-1]&0xf) == 2)
				score_cut += op_len * opt->gap_ex_D;
			else
				score_cut += opt->gap_open_D + op_len * opt->gap_ex_D;
			toff -= op_len;
		}
		else if (op == 3) // find the first intron operator
		{
			if (len_ref <= hash_kmer)
			// if (len_read <= hash_kmer)
			{
				if (type == 0) //for left ext
				{
					*qs += len_read;
					*ts += (len_ref + op_len);
				}
				else //for right ext
				{
					*qs -= len_read;
					*ts -= (len_ref + op_len);
				}
				*n_cigar = i;
				*score += score_cut;
			}
			else 
			{
				if ((type == 0) && (*ts < boundary))
				{
					*qs += len_read;
					*ts += (len_ref + op_len);
					*n_cigar = i;
					*score += score_cut;
				}
				else if ((type == 1) && (*ts > boundary))
				{
					*qs -= len_read;
					*ts -= (len_ref + op_len);
					*n_cigar = i;
					*score += score_cut;
				}
			}
			
			if (((type == 0) && (*ts >= boundary)) || ((type == 1) && (*ts <= boundary)))
				break;
			len_ref = 0;
			len_read = 0;
			score_cut = 0;
		}
	}
	
}

static int check_filter(uint32_t pos1, uint32_t pos2, uint32_t key1, uint32_t key2, uint32_t len)
{
	//filter the first and last anchor
	uint32_t begin = 0;
	if (chromosome_judge(pos1, &begin) != chromosome_judge(pos2, &begin))
		return 1;
	if ((key1 - key2 > 4) && (len < seed_k_t + 5))
		return 1;

	return 0;
}

static int call_mapping_bases(uint32_t *cigar, uint32_t n_cigar)
{
	int i;
	int op, op_len;
	int mapping_bases = 0;
	for (i = 0; i < n_cigar; ++i)
	{
		op = cigar[i]&0xf;
		op_len = (cigar[i]>>4);
		if (op == 0)
		{
			mapping_bases += op_len;
		}
	}

	return mapping_bases;
}

uint32_t remove_large_del_to_intron(uint32_t *cigar, uint32_t n_cigar)
{
    int i, j;
    int op, op_len;
    int op1, op_len1;
    int m_cigar = 0;
    int sig = 0;
    for(i = 0; i < n_cigar; ++i)
    {
        op = cigar[i]&0xf;
        op_len = (cigar[i] >> 4);

        if (op == 2 && op_len > 50)
        {
            cigar[i] = (op_len<<4) | 3; //N
            sig = 1;
        }
    }
    if (sig == 0)
        return n_cigar;

    i = 0;
    while(i < n_cigar)
    {
        op = cigar[i]&0xf;
        op_len = (cigar[i] >> 4);
        if (op == 3)
        {
            j = i + 1;
            int refl = 0;
            int readl = 0;
            while (j < n_cigar)
            {
                if (refl > 10)
                {
                    cigar[m_cigar] = cigar[i];
                    break;
                }

                op1 = cigar[j]&0xf;
                op_len1 = (cigar[j] >> 4);
                if (op1 == 0)
                {
                    refl += op_len1;
                    readl += op_len1;
                }
                else if(op1 == 1)
                {
                    readl += op_len1;
                }
                else if (op1 == 2)
                {
                    refl += op_len1;
                }
                else if (op1 == 3)
                {
                    cigar[m_cigar++] = (op_len + op_len1 + refl)<<4 | 3;
                    cigar[m_cigar] = (readl << 4) | 1;
                    i = j;
                    break;
                }
                ++j;
            }
            if (j == n_cigar)
            {
                cigar[m_cigar] = cigar[i];
            }
        }
        else
        {
            cigar[m_cigar] = cigar[i];
        }
        ++m_cigar;
        ++i;
    }
    return m_cigar;
}

void align_splic_FOR_REV(void *km, uint8_t *qseq, uint8_t *tseq, uint32_t qlen, uint32_t tlen, int splice_flag, param_map *opt, ksw_extz_t *ez, uint8_t splice_type, uint8_t *junc)
{
    if (! (opt->with_gtf))
        junc = NULL;

	if ((int64_t)tlen * qlen > opt->max_sw_mat)
	{
		ksw_reset_extz(ez);
        ez->n_cigar = 2;
        ez->cigar[0] = qlen<<4 | 1;
        ez->cigar[1] = tlen<<4 | 3;
        ez->score = 0;
	}
	else if(splice_type == 0) //left extend
	{
		ksw_exts2_sse(km, qlen, qseq, tlen, tseq, 5, mata_R, opt->gap_open_R, opt->gap_ex_R, opt->gap_open2_R, opt->noncan, opt->zdrop_R, splice_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez, junc);
	}
	else if(splice_type == 1) //right extend
	{
		ksw_exts2_sse(km, qlen, qseq, tlen, tseq, 5, mata_R, opt->gap_open_R, opt->gap_ex_R, opt->gap_open2_R, opt->noncan, opt->zdrop_R, splice_flag|KSW_EZ_EXTZ_ONLY, ez, junc);
	}
	else  //end-to-end
	{
		ksw_exts2_sse(km, qlen, qseq, tlen, tseq, 5, mata_R, opt->gap_open_R, opt->gap_ex_R, opt->gap_open2_R, opt->noncan, opt->zdrop_R, splice_flag|KSW_EZ_APPROX_MAX, ez, junc);
	}
}

void align_non_splice(void *km, uint8_t *qseq, uint8_t *tseq, uint32_t qlen, uint32_t tlen, param_map *opt, ksw_extz_t *ez, int bandwith, int flag, uint8_t type)
{
	if ((int64_t)tlen * qlen > opt->max_sw_mat)
	{
		ksw_reset_extz(ez);
        ez->n_cigar = 2;
        ez->cigar[0] = qlen<<4 | 1;
        ez->cigar[1] = tlen<<4 | 3;
        ez->score = 0;
	}
	else if (type == 0) //left extend
	{
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mata_D, opt->gap_open_D, opt->gap_ex_D, opt->gap_open2_D, opt->gap_ex2_D, bandwith, opt->zdrop_D, opt->end_bonus, flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
	}
	else if(type == 1) //right extend
	{
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mata_D, opt->gap_open_D, opt->gap_ex_D, opt->gap_open2_D, opt->gap_ex2_D, bandwith, opt->zdrop_D, opt->end_bonus, flag|KSW_EZ_EXTZ_ONLY, ez);
	}
	else if(type == 2)//end-to-end 
	{
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mata_D, opt->gap_open_D, opt->gap_ex_D, opt->gap_open2_D, opt->gap_ex2_D, bandwith, opt->zdrop_D, opt->end_bonus, flag|KSW_EZ_APPROX_MAX, ez);
	}
	else ////left extend but not reverse cigar
	{
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mata_D, opt->gap_open_D, opt->gap_ex_D, opt->gap_open2_D, opt->gap_ex2_D, bandwith, opt->zdrop_D, opt->end_bonus, flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT, ez);
	}
}

//static void merge_del_to_intron(_aln_t *aln)
//{
//    uint32_t *cigar = aln->cigar;
//    int m_cigar = aln->n_cigar;
//
//    int i;
//    int n_cigar = 0;
//    int op, op_len;
//    int thre = 30;
//    
//    //tran large 'D' to 'N' (> 30)
//    for(i = 0; i < m_cigar; ++i)
//    {
//        op = cigar[i]&0xf;
//		op_len = (cigar[i]>>4);
//        if (op == 2 && op_len > 30)
//    }
//}

static void align_core_primary(void *km, uint32_t seqlen, TARGET_t *anchor_map2ref, uint8_t *qseq0[2], uint8_t *qual0[2], _aln_t *aln, param_map *opt, ksw_extz_t *ez, ksw_extz_t *ez2, REF_t *ref_pos, REF_t *ref_temp, QUERY_t *query_pos, int *chr_n, uint32_t anchor_n, uint8_t strand, uint8_t tid, int key_total, int splice_flag, uint32_t left_bound, uint32_t right_bound)
{
	int i;
    uint32_t a_n = 0;
    int extra_flag_D = 0;
	int extra_flag_R = 0;
    int bandwith = (int)(opt->bw * 1.5 + 1.);
    uint8_t *qseq;
    uint32_t qlen, tlen;
    uint32_t key1, key2;
	uint32_t s_s;
	uint32_t s_e;
	int s_len;
	int intron_len;
	uint32_t key;
	int key_l;
	int key_r;
	uint32_t pre_pos = 0;
	uint8_t *tseq = NULL;
    uint8_t *junc = NULL;
    
	int score = KSW_NEG_INF;
	int score1 = KSW_NEG_INF;

	// if (splice_flag & MM_F_SPLICE_FOR) extra_flag_R |= strand? KSW_EZ_SPLICE_REV : KSW_EZ_SPLICE_FOR;
	// if (splice_flag & MM_F_SPLICE_REV) extra_flag_R |= strand? KSW_EZ_SPLICE_FOR : KSW_EZ_SPLICE_REV;
	extra_flag_R |= splice_flag; 
	extra_flag_R |= KSW_EZ_SPLICE_FLANK;

	uint32_t qs, qe, ts0 ,ts1, te0, te1;
	uint32_t _qs, _qe, _ts, _te;
	int a_len, b_len;
	uint32_t te_;
	uint32_t qe_;
	uint32_t te_s;
	uint32_t qe_s;
	uint16_t thre1;  //NonaSim: the avearge length of del and ins for ONT2D is 1.81bp and 1.69bp
	ksw_extz_t *ez_tmp;
	int *exon_find = (int* )calloc(e_shift,4);
	uint32_t find_cnt;

	qs = query_pos[a_n].qs;
	qe = query_pos[a_n].qe;
	ts0 = ref_pos[a_n].ts;
	te0 = ref_pos[a_n].te;
	ts1 = ref_temp[a_n].ts;
	te1 = ref_temp[a_n].te;
	a_n++;

	_qs = qs;
	_ts = ts0;


	key1 = ref_temp[0].key;
	key_l = key1;

	/*
	if two anchor in different chrosome, filter, later added
	*/
	if (qs > 0)
	{
		for(i = key1; i > 0; --i)
		{
			if (anchor_map2ref[i].ts - anchor_map2ref[i-1].te > max_intron_length)
			{
				break;
			}
		}
		key_l = (i < 0)? 0 : i;
	}

	key2 = ref_temp[anchor_n - 1].key + 1;
	key_r = key2;
	if ((seqlen - query_pos[anchor_n - 1].qe > 0) && (key2 <= key_total))
	{
		for(i = key2; i < key_total; ++i)
		{
			if (anchor_map2ref[i].ts - anchor_map2ref[i-1].te > max_intron_length)
			{
				break;
			}
		}
		// key_r = (i == key_total)? (key_total - 1) : i;
		key_r = i;
	}

	tlen = anchor_map2ref[key_r - 1].te - anchor_map2ref[key_l].ts;
	if (tlen < seqlen)
	{
		tlen = seqlen + 2*(seqlen * opt->match_R - opt->gap_open_R)/opt->gap_ex_R;
	}

    uint32_t LEN_t = tlen + 2*opt->max_extend_gap;
	tseq = (uint8_t*)kmalloc(km, LEN_t);
    if (opt->with_gtf)
        junc = (uint8_t*)kmalloc(km, LEN_t);

	//left extension
	int l;
    uint8_t *temp_ref_left = NULL;
	if (qs > 0)
	{	
		qlen = qs;
		thre1 = qlen*opt->error_overall*opt->error_ins*2 + 5; 
		a_len = (ts0 > ts1)? (ts0 - ts1) : 0;
		qseq = &qseq0[strand][0];
		// key_l = 0;
		key1 = ref_temp[0].key;
		if (((int)(qs + ts1 - ts0) > thre1) && (qs - a_len > hash_kmer) && (qs < opt->max_extend_left_right))  //local hash procedure
		{
			qlen = qs - a_len;
			find_cnt = local_hash_anchor(qseq, qlen, exon_find, anchor_map2ref, key_l, key1, opt, 1, tid, 0);

			if (find_cnt > 0)
			{
				pre_pos = 0;
				for (i = 0; i < find_cnt; ++i)
				{
					//connect the ref seq
					key = exon_find[i] + key_l;
					s_s = anchor_map2ref[key].ts;
					s_len = anchor_map2ref[key].te - anchor_map2ref[key].ts + 1;
					get_ref_onebyone(tseq, s_s, s_len, pre_pos);
					pre_pos += s_len;
				}
				if (a_len > 0)
				{
					get_ref_onebyone(tseq, ts1, a_len, pre_pos);
					pre_pos += a_len;
				}

				qlen = qs;
				qseq = &qseq0[strand][0];
				mm_seq_rev(qlen, qseq);
				mm_seq_rev(pre_pos, tseq);

                temp_ref_left = (uint8_t* )calloc(pre_pos, 1);
                for (i = 0; i < pre_pos; ++i)
                {
                    temp_ref_left[i] = tseq[i];
                }

				align_non_splice(km, qseq, tseq, qs, pre_pos, opt, ez2, bandwith, extra_flag_D, 3);
				mm_seq_rev(qlen, qseq);
				if (ez2->n_cigar > 0)
				{
					score = ez2->max;
				}
				else
				{
					score = KSW_NEG_INF;
				}
			}
			else
			{
				score = KSW_NEG_INF;
				ez2->n_cigar = 0;
			}
		}
		else
		{
			score = KSW_NEG_INF;
			ez2->n_cigar = 0;
		}

		//left extension directly
		qlen = qs;
		l = qlen;
		int l1 = (l * opt->match_R - opt->gap_open_R)/opt->gap_ex_R;
		l1 = (l1 > 0)? l1 : 0;
		l += l1;
		l = l < opt->max_extend_gap? l : opt->max_extend_gap;
		l = (l > ts0)? ts0 : l;
		if (ts0 - l < chr_end_n[*chr_n - 1])
			l = ts0 - chr_end_n[*chr_n - 1];
		get_refseq(tseq, l, ts0 - l);
        
        if (opt->with_gtf)
        {
            get_junc(*chr_n, ts0 - l, l, junc, opt->with_gtf);
            mm_seq_rev(l, junc);
        }
        
		mm_seq_rev(qlen, qseq);
		mm_seq_rev(l, tseq);
		align_splic_FOR_REV(km, qseq, tseq, qlen, l, extra_flag_R, opt, ez, 0, junc);
		mm_seq_rev(qlen, qseq);
		if(ez->n_cigar > 0)
		{
			_ts = ts0 - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1); //reference start
			_qs = qs - (ez->reach_end? 0: ez->max_q + 1); //query start

			score1 = ez->max;
		}
		else
		{
			_ts = ts0;
			_qs = qs;
		}
		
		//compare
		if (score <= score1)
		{
			if (ez->n_cigar > 0)
			{
				qlen = qs - _qs;
				tlen = ts0 - _ts;
				qseq = &qseq0[strand][_qs];

				mm_append_cigar(aln, ez->n_cigar, ez->cigar);
				
				aln->dp_score += score1;
			}
		}
		else
		{
			if (ez2->n_cigar > 0)
			{
				_ts = ts0 - (ez2->reach_end? ez2->mqe_t + 1 : ez2->max_t + 1); //reference start
				_qs = qs - (ez2->reach_end? 0: ez2->max_q + 1); //query start
				
				qlen = qs - _qs;
				tlen = ts0 - _ts;
				qseq = &qseq0[strand][_qs];
                mm_seq_rev(qlen, qseq);
				//push intron
				int ext_len = ts0 - _ts;
				pre_pos = 0;
				if (_ts < ts1) //include
				{	
					intron_len = ts1 - anchor_map2ref[exon_find[find_cnt - 1]+key_l].te - 1;
					pre_pos = append_intron_to_cigar(km, ez2, pre_pos, a_len, intron_len);
					ext_len -= a_len;
					_ts -= intron_len;

					// for(i = key1 - 1; i > key_l; --i)
					for(i = find_cnt - 1; i > 0; --i)
					{
						key = exon_find[i] + key_l;
						s_len = anchor_map2ref[key].te - anchor_map2ref[key].ts + 1;
						intron_len = anchor_map2ref[key].ts - anchor_map2ref[exon_find[i-1]+key_l].te - 1;
						if(s_len < ext_len)
						{
							pre_pos = append_intron_to_cigar(km, ez2, pre_pos, s_len, intron_len);
							_ts -= intron_len;
							ext_len -= s_len;
						}
						else
						{
							break;
						}
					}
				}
	
				check_cigar(qseq, temp_ref_left, ez2->cigar, &(ez2->n_cigar), &_qs, &_ts, &score, qlen, tlen, left_bound, 0, opt);
	
				mm_seq_rev(qlen, qseq);
				//reverse cigar
				int tmp;
				for (i = 0; i < ez2->n_cigar>>1; ++i) // reverse CIGAR
				{
					tmp = ez2->cigar[i]; 
					ez2->cigar[i] = ez2->cigar[ez2->n_cigar - 1 - i]; 
					ez2->cigar[ez2->n_cigar - 1 - i] = tmp;
				}

				//check re-align
				qseq = &qseq0[strand][_qs];
				qlen = qs - _qs;
				tlen = ts0 - _ts;
				int sig = 0;
				sig = check_realign(km, opt, bandwith, qseq, qlen, tlen, ez2, score, _ts, find_cnt + 10, splice_flag);

				mm_append_cigar(aln, ez2->n_cigar, ez2->cigar);/////
				if (sig)
					score = ez2->score;
				aln->dp_score += score;	
			}
		}
	}
	while(a_n < anchor_n)
	{
		key1 = ref_temp[a_n - 1].key;
		key2 = ref_temp[a_n].key;
		if (key2 == key1) //the same id of exon
		{
			te1 = te0;
			ref_temp[a_n].ts = ref_pos[a_n].ts;
		}
		uint8_t shift1 = 0;
		te_ = (te0 < te1)? te0 : te1;
		qe_ = qe - (te0 - te_);
		int q = qe_ - qs;
		int t = te_ - ts0;

		if (q < 5 || t < 5)
		{
			te_ = te0;
			te1 = te0;
			qe_ = qe;
		}
		else
		{
			if ((q > 2*seed_k_t) && (t > 2*seed_k_t))
			{
				shift1 = seed_k_t;
			}
			else if ((q > seed_k_t) && (t > seed_k_t))
			{
				shift1 = seed_k_t >> 1;
			}
		}
		
		uint8_t shift2 = 0;
		te_s = (ref_pos[a_n].ts < ref_temp[a_n].ts)? ref_temp[a_n].ts : ref_pos[a_n].ts;
		qe_s = query_pos[a_n].qs + (te_s - ref_pos[a_n].ts);
		uint32_t TE = (ref_pos[a_n].te < ref_temp[a_n].te)? ref_pos[a_n].te : ref_temp[a_n].te;
		int QE = query_pos[a_n].qe + (TE - ref_pos[a_n].te);
		int TT = TE - te_s;
		int QQ = QE - qe_s;

		if (TT < 5 || QQ < 5)
		{
			te_s = ref_pos[a_n].ts;
			ref_temp[a_n].ts = ref_pos[a_n].ts;
			qe_s = query_pos[a_n].qs;
		}
		else 
		{
			if ((QQ > 2*seed_k_t) && (TT > 2*seed_k_t))
			{
				shift2 = seed_k_t;
			}
			else if ((QQ > seed_k_t) && (TT > seed_k_t))
			{
				shift2 = seed_k_t >> 1;
			}
		}

		te_ -= shift1;
		qe_ -= shift1;
		te_s += shift2;
		qe_s += shift2;
			
		qseq = &qseq0[strand][qs];
		tlen = te_ - ts0 + 1;
		qlen = qe_ - qs + 1;

		get_refseq(tseq, tlen, ts0);
		align_non_splice(km, qseq, tseq, qlen, tlen, opt, ez, bandwith, extra_flag_D, 2);
		if (ez->n_cigar > 0)
		{
			mm_append_cigar(aln, ez->n_cigar, ez->cigar);
			aln->dp_score += ez->score;
		}

		a_len = te1 - te_;
		b_len = te_s - ref_temp[a_n].ts;
		qlen = qe_s - qe_ - 1;

		if(qlen == 0)
		{
			aln->cigar[aln->n_cigar++] = (te_s - te_ - 1)<<4 | 3;
			goto END;
		}
		tlen = a_len + b_len;
		if (key2 - key1 <= 1) //have no exon between
		{	
			//for non-splice
			qseq = &qseq0[strand][qe_ + 1];
			if (key2 == key1)
			{
				ez->score = KSW_NEG_INF;
			}
			else
			{
				//for non_splice aligner
				get_ref_jump(tseq, te_ + 1, a_len, ref_temp[a_n].ts, b_len);
				align_non_splice(km, qseq, tseq, qlen, tlen, opt, ez, bandwith, extra_flag_D, 2);
                
                if (ez->n_cigar > 0)
                {
    				append_intron_to_cigar(km, ez, 0, a_len, ref_temp[a_n].ts - te1 - 1);
    				check_realign(km, opt, bandwith, qseq, qlen, tlen + ref_temp[a_n].ts - te1 - 1, ez, ez->score, te_ + 1, 10, splice_flag);
                }			
			}

			//for splice aligner
			tlen = te_s - te_ - 1;			
			get_refseq(tseq, tlen, te_ + 1);
            if (opt->with_gtf)
                get_junc(*chr_n, te_ + 1, tlen, junc, opt->with_gtf);
			align_splic_FOR_REV(km, qseq, tseq, qlen, tlen, extra_flag_R, opt, ez2, 2, junc);
			
			if (ez->score > ez2->score)
			{
				ez_tmp = ez;
			}
			else if((ez->score <= ez2->score) && (ez2->score > KSW_NEG_INF))
			{
				ez_tmp = ez2;
			}
			else
			{
				ez->n_cigar = 2;
				ez->score = 0;
				//cigar
				ez->cigar[0] = qlen<<4 | 1;
				ez->cigar[1] = tlen<<4 | 2;

				ez_tmp = ez;
			}
			
			mm_append_cigar(aln, ez_tmp->n_cigar, ez_tmp->cigar);
			aln->dp_score += ez_tmp->score;

		}
		else if ((key2 - key1 == 2) && (anchor_map2ref[key1+1].te - anchor_map2ref[key1+1].ts + 1 < 51))
		{
			if (qlen > hash_kmer)
			{
				tlen = 0;
				get_ref_onebyone(tseq, te_ + 1, a_len, 0);
				tlen += a_len;
				s_s = anchor_map2ref[key1+1].ts;
				s_len = anchor_map2ref[key1+1].te - anchor_map2ref[key1+1].ts + 1;
				get_ref_onebyone(tseq, s_s, s_len, tlen);
				tlen += s_len;
				get_ref_onebyone(tseq, ref_temp[a_n].ts, b_len, tlen);
				tlen += b_len;

				qseq = &qseq0[strand][qe_ + 1];
				qlen = qe_s - qe_ - 1;
				//DNA aligner
				align_non_splice(km, qseq, tseq, qlen, tlen, opt, ez, bandwith, extra_flag_D, 2);

				if (ez->n_cigar > 0)
				{
					pre_pos = 0;
					pre_pos = append_intron_to_cigar(km, ez, pre_pos, a_len, s_s - te1 - 1);
					pre_pos = append_intron_to_cigar(km, ez, pre_pos, s_len, ref_temp[a_n].ts - anchor_map2ref[key1+1].te - 1);              
				
					uint32_t tmp = tlen + s_s - te1 - 1 + ref_temp[a_n].ts - anchor_map2ref[key1+1].te - 1;

                	check_realign(km, opt, bandwith, qseq, qlen, tmp, ez, ez->score, te_ + 1, 10, splice_flag);
                }
                else
                {
                    ez->score = KSW_NEG_INF;
                    ez->n_cigar = 0;
                }              
            }
			else
			{
				//fprintf(stderr, "have only one exon between condition 2, qlen = %d, tlen = %d\n", qlen, tlen);
				get_ref_jump(tseq, te_ + 1, a_len, ref_temp[a_n].ts, b_len);
				qseq = &qseq0[strand][qe_ + 1];
				align_non_splice(km, qseq, tseq, qlen, tlen, opt, ez, bandwith, extra_flag_D, 2);

				if (ez->n_cigar > 0)
				{
					append_intron_to_cigar(km, ez, 0, a_len, ref_temp[a_n].ts - te1 - 1);
					uint32_t tmp = tlen + ref_temp[a_n].ts - te1 - 1;
                	check_realign(km, opt, bandwith, qseq, qlen, tmp, ez, ez->score, te_ + 1, 10, splice_flag);
				}
                else
                {
                    ez->score = KSW_NEG_INF;
                    ez->n_cigar = 0;
                }
			}

			//for RNA aligner
			tlen = te_s - te_ - 1;			
			get_refseq(tseq, tlen, te_ + 1);
            if(opt->with_gtf)
                get_junc(*chr_n, te_ + l, tlen, junc, opt->with_gtf);
			align_splic_FOR_REV(km, qseq, tseq, qlen, tlen, extra_flag_R, opt, ez2, 2, junc);
			
			ez_tmp = (ez->score > ez2->score)? ez : ez2;

			mm_append_cigar(aln, ez_tmp->n_cigar, ez_tmp->cigar);
			if (ez_tmp->n_cigar > 0)
				aln->dp_score += ez_tmp->score;
		}
		else //have exons between, key2 - key1 > 1
		{
            if (qlen > Eindel && key2 - key1 < 10)
			{
				// fprintf(stderr, "key2 - key1 > 1, key1 = %u, key2 = %u, condition 1_find exons\n", key1, key2);
				key1++;
				//find exons between
				uint32_t n = key2 - key1;
				uint32_t *idx_cnt_array;
				idx_cnt_array = (int* )calloc(n, sizeof(int));

                qseq = &qseq0[strand][qe_ + 1];
                qlen = qe_s - qe_ - 1;
				
				uint32_t real_cnt = local_hash_anchor(qseq, qlen, idx_cnt_array, anchor_map2ref, key1, key2, opt, 0, tid, 1);

				if (real_cnt > 0)
				{
					pre_pos = 0;
					get_ref_onebyone(tseq, te_ + 1, a_len, 0);
					pre_pos += a_len;
					for (i = 0; i < real_cnt; ++i)
					{
						key = key1 + idx_cnt_array[i];
						s_s = anchor_map2ref[key].ts;
						s_len = anchor_map2ref[key].te - anchor_map2ref[key].ts + 1;
						get_ref_onebyone(tseq, s_s, s_len, pre_pos);
						pre_pos += s_len;
					}
					get_ref_onebyone(tseq, ref_temp[a_n].ts, b_len, pre_pos);
					pre_pos += b_len;

					align_non_splice(km, qseq, tseq, qlen, pre_pos, opt, ez, bandwith, extra_flag_D, 2);
                    uint32_t tmp = pre_pos;
					if (ez->n_cigar > 0)
					{
						//append intron
						pre_pos = 0;
						s_s = anchor_map2ref[key1 + idx_cnt_array[0]].ts;
                        tmp += s_s - te1 - 1;
						pre_pos = append_intron_to_cigar(km, ez, pre_pos, a_len, s_s - te1 - 1);

						for (i = 0; i < real_cnt - 1; ++i)
						{
							key = key1 + idx_cnt_array[i]; 
							// s_s = anchor_map2ref[key].ts;
							s_len = anchor_map2ref[key].te - anchor_map2ref[key].ts + 1;
							intron_len = anchor_map2ref[key1 + idx_cnt_array[i+1]].ts - anchor_map2ref[key].te - 1;
                            tmp += intron_len;
							pre_pos = append_intron_to_cigar(km, ez, pre_pos, s_len, intron_len);
						}
						s_s = anchor_map2ref[key1 + idx_cnt_array[real_cnt - 1]].ts;
						s_e = anchor_map2ref[key1 + idx_cnt_array[real_cnt - 1]].te;
                        tmp += ref_temp[a_n].ts - s_e - 1;
						pre_pos = append_intron_to_cigar(km, ez, pre_pos, s_e - s_s + 1, ref_temp[a_n].ts - s_e - 1);

                        check_realign(km, opt, bandwith, qseq, qlen, tmp, ez, ez->score, te_ + 1, (real_cnt<<1) + 10, splice_flag);
					}
                    else
                    {
                        ez->score = KSW_NEG_INF;
                        ez->n_cigar = 0;
                    }
				}	
                else
				{	
					get_ref_jump(tseq, te_ + 1, a_len, ref_temp[a_n].ts, b_len);
					qseq = &qseq0[strand][qe_ + 1];
					align_non_splice(km, qseq, tseq, qlen, tlen, opt, ez, bandwith, extra_flag_D, 2);

                    if (ez->cigar > 0)
                    {
					    append_intron_to_cigar(km, ez, 0, a_len, ref_temp[a_n].ts - te1 - 1);
                        uint32_t tmp = tlen + ref_temp[a_n].ts - te1 - 1;
                        check_realign(km, opt, bandwith, qseq, qlen, tmp, ez, ez->score, te_ + 1, 10, splice_flag);
                    }
                    else
                    {
                        ez->score = KSW_NEG_INF;
                        ez->n_cigar = 0;
                    }
				}
				free(idx_cnt_array);
			}
			else
            {
                get_ref_jump(tseq, te_ + 1, a_len, ref_temp[a_n].ts, b_len);
                qseq = &qseq0[strand][qe_ + 1];
                align_non_splice(km, qseq, tseq, qlen, tlen, opt, ez, bandwith, extra_flag_D, 2);

                if (ez->cigar > 0)
                {
                    append_intron_to_cigar(km, ez, 0, a_len, ref_temp[a_n].ts - te1 - 1);
                    uint32_t tmp = tlen + ref_temp[a_n].ts - te1 - 1;
                    check_realign(km, opt, bandwith, qseq, qlen, tmp, ez, ez->score, te_ + 1, 10, splice_flag);
                }
                else
                {
                    ez->score = KSW_NEG_INF;
                    ez->n_cigar = 0;
                }
            }

			//for splice aligner
            qseq = &qseq0[strand][qe_ + 1];
            qlen = qe_s - qe_ - 1;
			tlen = te_s - te_ - 1;
				
			get_refseq(tseq, tlen, te_ + 1);
            if(opt->with_gtf)
                get_junc(*chr_n, te_ + 1, tlen, junc, opt->with_gtf);
			align_splic_FOR_REV(km, qseq, tseq, qlen, tlen, extra_flag_R, opt, ez2, 2, junc);

			ez_tmp = (ez->score > ez2->score)? ez : ez2;

			mm_append_cigar(aln, ez_tmp->n_cigar, ez_tmp->cigar);
			if (ez_tmp->n_cigar > 0)
				aln->dp_score += ez_tmp->score;
		}
END:
		qs = qe_s;
		qe = query_pos[a_n].qe;
		ts0 = te_s;
		te0 = ref_pos[a_n].te;
		ts1 = ts0;
		te1 = ref_temp[a_n].te;
		a_n++;
	}

	//align the last anchor
	qlen = qe - qs + 1;
	tlen = te0 - ts0 + 1;
	qseq = &qseq0[strand][qs];
	get_refseq(tseq, tlen, ts0);

	align_non_splice(km, qseq, tseq, qlen, tlen, opt, ez, bandwith, extra_flag_D, 2);

	if (ez->n_cigar > 0)
    {
        mm_append_cigar(aln, ez->n_cigar, ez->cigar);
        aln->dp_score += ez->score;
    }

	//right extension
	_te = te0;
	_qe = qe;
    qlen = seqlen - qe - 1;
	pre_pos = 0;
	score = KSW_NEG_INF;
	score1 = KSW_NEG_INF;
    uint8_t *temp_ref_right = NULL;
	if (qlen > 0)
	{
		thre1 = qlen*opt->error_overall*opt->error_ins*2 + 5; 
		a_len = (te1 > te0)? (te1 - te0) : 0;
		key2 = ref_temp[anchor_n - 1].key + 1;
		// key_r = key_total;
		if (((int)(qlen + te0 - te1) > thre1) && (qlen - a_len > hash_kmer) && qlen < opt->max_extend_left_right)
		{

			qseq = &qseq0[strand][qe+a_len];
			qlen -= a_len;
		
			find_cnt = local_hash_anchor(qseq, qlen, exon_find, anchor_map2ref, key2, key_r, opt, 1, tid, 0);

			if (find_cnt > 0)
			{
				pre_pos = 0;
				if (a_len > 0)
				{
					get_ref_onebyone(tseq, te0+1, a_len, pre_pos);
					pre_pos += a_len;
				}
				for (i = 0; i < find_cnt; ++i)
				{
					//connect the ref seq
					key = key2 + exon_find[i];
					s_s = anchor_map2ref[key].ts;
					s_len = anchor_map2ref[key].te - anchor_map2ref[key].ts + 1;
					get_ref_onebyone(tseq, s_s, s_len, pre_pos);
					pre_pos += s_len;
				}

				qlen = seqlen - qe - 1;
				qseq = &qseq0[strand][qe+1];
                
                int i;
                temp_ref_right = (uint8_t* )calloc(pre_pos, sizeof(uint8_t));

                for (i = 0; i < pre_pos; ++i)
                {
                    temp_ref_right[i] = tseq[i];
                }

				//test ------
				align_non_splice(km, qseq, temp_ref_right, qlen, pre_pos, opt, ez2, bandwith, extra_flag_D, 1);

				if (ez2->n_cigar > 0)
				{
					score = ez2->max;
				}
				else
				{
					score = KSW_NEG_INF;
				}
			}
			else
			{
				score = KSW_NEG_INF;
				ez2->n_cigar = 0;
			}
		}
		else
		{
			score = KSW_NEG_INF;
			ez2->n_cigar = 0;
		}

		//right extension directly
		qlen = seqlen - qe - 1;
		qseq = &qseq0[strand][qe+1];
		l = qlen;
		int l2 = (l * opt->match_R - opt->gap_open_R)/opt->gap_ex_R;
		l2 = (l2 > 0)? l2 : 0;
		l += l2;
		l = l < opt->max_extend_gap? l : opt->max_extend_gap;			
		get_refseq(tseq, l, te0 + 1);
        if(opt->with_gtf)
            get_junc(*chr_n, te0 + 1, l, junc, opt->with_gtf);
		align_splic_FOR_REV(km, qseq, tseq, qlen, l, extra_flag_R, opt, ez, 1, junc);

		if(ez->n_cigar > 0)
		{
			_te = te0 + (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1); //reference end
			_qe = qe + (ez->reach_end? seqlen - qe - 1 : ez->max_q + 1); //query end;
			score1 = ez->max;
		}
		else
		{
			_te = te0;
			_qe = qe;
		}

		//compare
		if (score <= score1)
		{
			if (ez->n_cigar > 0)
			{
				qlen = _qe - qe;
				tlen = _te - te0;
				
				mm_append_cigar(aln, ez->n_cigar, ez->cigar);
				if (_qe < seqlen - 1)//have soft clip
				{
					aln->cigar[aln->n_cigar++] = (seqlen - 1 - _qe)<<4 | 4;
				}
				aln->dp_score += score1;
			}
			else
			{
				aln->cigar[aln->n_cigar++] = (seqlen - 1 - qe)<<4 | 4;
			}
		}
		else
		{
			if (ez2->n_cigar > 0)
			{
				_te = te0 + (ez2->reach_end? ez2->mqe_t + 1 : ez2->max_t + 1); //reference end
				_qe = qe + (ez2->reach_end? seqlen - qe - 1 : ez2->max_q + 1); //query end;

				qlen = _qe - qe;
				tlen = _te - te0;
				qseq = &qseq0[strand][qe+1];

				int ext_len = _te - te0;
				pre_pos = 0;
				if (_te > te1) //include
				{
					intron_len = anchor_map2ref[exon_find[0] + key2].ts - te1 - 1;
					pre_pos = append_intron_to_cigar(km, ez2, pre_pos, a_len, intron_len);
					ext_len -= a_len;
					_te += intron_len;

					for(i = 0; i < find_cnt - 1; ++i)
					{
						key = key2 + exon_find[i];
						s_len = anchor_map2ref[key].te - anchor_map2ref[key].ts + 1;
						intron_len = anchor_map2ref[key2 + exon_find[i + 1]].ts - anchor_map2ref[key].te - 1;
						if(s_len < ext_len)
						{
							pre_pos = append_intron_to_cigar(km, ez2, pre_pos, s_len, intron_len);
							_te += intron_len;
							ext_len -= s_len;
						}
						else
						{
							break;
						}
					}
				}

				check_cigar(qseq, temp_ref_right, ez2->cigar, &(ez2->n_cigar), &_qe, &_te, &score, qlen, tlen, right_bound, 1, opt);
				
				qseq = &qseq0[strand][qe+1];
				qlen = _qe - qe;
				tlen = _te - te0;
				int sig = 0;
				sig = check_realign(km, opt, bandwith, qseq, qlen, tlen, ez2, score, te0 + 1, find_cnt + 10, splice_flag);
				//append cigar to aln
				mm_append_cigar(aln, ez2->n_cigar, ez2->cigar);/////
				
				if (sig)
					score = ez2->score;

				aln->dp_score += score;
				// aln->dp_score += ez2->score;
				if (_qe < seqlen - 1)//have soft clip
				{
					aln->cigar[aln->n_cigar++] = (seqlen - 1 - _qe)<<4 | 4;
				}	
			}
			else
			{
				aln->cigar[aln->n_cigar++] = (seqlen - 1 - qe)<<4 | 4;
			}
		}
	}

	//update
	qlen = _qe - _qs + 1;
	tlen = _te - _ts + 1;
	qseq = &qseq0[strand][_qs];
    if (tlen > LEN_t)
    {
        tseq = (uint8_t *)krealloc(km, tseq, tlen);
    }
	get_refseq(tseq, tlen, _ts);
	mm_update_extra(aln, qseq, tseq, qlen, tlen, &_qs, &_ts, opt->gap_open_R, opt->gap_ex_R);

	//ending I/D 
	int op = aln->cigar[aln->n_cigar - 1]&0xf;
	int op_len;
	if ( op == 1 ) //ending I
	{
		op_len = aln->cigar[aln->n_cigar - 1] >> 4;
		aln->cigar[aln->n_cigar - 1] = op_len << 4 | 4;
	}
	else if (op == 2)// ending D
	{
		aln->n_cigar -= 1;
	}

	if (_qs > 0)
	{
		int j;
		for (j = aln->n_cigar - 1; j >= 0; --j)
		{
			aln->cigar[j + 1] = aln->cigar[j];
		}
		aln->cigar[0] = _qs<<4 | 4;
		aln->n_cigar ++;
	}
	
	//set alignment start pos
	aln->_1_based_pos = _ts;

    //remove large del 
    aln->n_cigar = remove_large_del_to_intron(aln->cigar, aln->n_cigar);

    if (temp_ref_right != NULL) free(temp_ref_right);
    if (temp_ref_left != NULL) free(temp_ref_left);
	free(exon_find);
	kfree(km, tseq);
    kfree(km, junc);
}

static int find_merge_anchor(TARGET_t *anchor_map2ref, REF_t *ref_temp, REF_t *ref_pos, QUERY_t *query_pos, uint32_t anchor_n, int primary, uint32_t *new_n, param_map *opt)
{
	uint32_t i,j;
	int key = 0;
	int key_1 = 0;
	int key_2 = merge_anchor_cnt - 1;
	// int signal = 0;
	key = binarysearch_anchor(anchor_map2ref, key_1, key_2, ref_pos[0].ts, ref_pos[0].te);
	if (key != -1)
	{
		ref_temp[0].ts = anchor_map2ref[key].ts; //remove
		ref_temp[0].te = anchor_map2ref[key].te; //remove
		ref_temp[0].key = key;
		key_1 = key;
	}else
	{
		ref_temp[0].key = -1;
	}
	key = binarysearch_anchor(anchor_map2ref, key_1, key_2, ref_pos[anchor_n - 1].ts, ref_pos[anchor_n - 1].te);
	if (key != -1)
	{
		ref_temp[anchor_n - 1].ts = anchor_map2ref[key].ts; //remove
		ref_temp[anchor_n - 1].te = anchor_map2ref[key].te; //remove
		ref_temp[anchor_n - 1].key = key;
		key_2 = key;
	}else
	{
		ref_temp[anchor_n - 1].key = -1;
	}
	

	for (i = 1; i < anchor_n/2; ++i)
	{
		key = binarysearch_anchor(anchor_map2ref, key_1, key_2, ref_pos[i].ts, ref_pos[i].te);
		if (key != -1)
		{
			ref_temp[i].ts = anchor_map2ref[key].ts; //remove
			ref_temp[i].te = anchor_map2ref[key].te; //remove
			ref_temp[i].key = key;
			key_1 = key;
		}else
		{
			ref_temp[i].key = -1;
		}
		
		key = binarysearch_anchor(anchor_map2ref, key_1, key_2, ref_pos[anchor_n - 1 - i].ts, ref_pos[anchor_n - 1 - i].te);
		if (key != -1)
		{
			ref_temp[anchor_n - 1 - i].ts = anchor_map2ref[key].ts; //remove
			ref_temp[anchor_n - 1 - i].te = anchor_map2ref[key].te; //remove
			ref_temp[anchor_n - 1 - i].key = key;
			key_2 = key;
		}else
		{
			ref_temp[anchor_n - 1 - i].key = -1;
		}
		
	}
	if (anchor_n % 2 == 1)
	{
		key = binarysearch_anchor(anchor_map2ref, key_1, key_2, ref_pos[anchor_n/2].ts, ref_pos[anchor_n/2].te);
		if (key != -1)
		{
			ref_temp[anchor_n/2].ts = anchor_map2ref[key].ts; //remove
			ref_temp[anchor_n/2].te = anchor_map2ref[key].te; //remove
			ref_temp[anchor_n/2].key = key;
		}else
		{
			ref_temp[anchor_n/2].key = -1;
		}	
	}

	//adjust anchors, if two anchor are in the sam interval of exon, then merge the two anchor
	//and remove anchor whose key = -1;
	uint32_t anchor_n_new = 0;
	uint32_t index;
    int tlen;
    int qlen;
    int thre;
	i = 0;

MULFIND:
	while(i < anchor_n)
	{
		key = ref_temp[i].key;
		if (key == -1)
		{
			i++;
			continue;
		}
		index = i;
		i++;
		while ((i < anchor_n) && (ref_temp[i].key == key))
		{
			//if the gap less than thre, then merge, else seperate,  method 2
            tlen = ref_pos[i].ts - ref_pos[i - 1].te;
            qlen = query_pos[i].qs - query_pos[i - 1].qe;
            thre = tlen*0.134*0.36*2;
            if ((abs)(tlen - qlen) > thre) //judge whether can be merge by tlen - qlen < thre
            {
                break;
            }
			i++;
		}
		query_pos[anchor_n_new].qs = query_pos[index].qs;
		query_pos[anchor_n_new].qe = query_pos[i - 1].qe;
		ref_pos[anchor_n_new].ts = ref_pos[index].ts;
		ref_pos[anchor_n_new].te = ref_pos[i - 1].te;
		ref_temp[anchor_n_new].ts = anchor_map2ref[key].ts;
		ref_temp[anchor_n_new].te = anchor_map2ref[key].te;
		ref_temp[anchor_n_new].key = key;
		anchor_n_new++;
		break;
	}
    //remove the first anchor if not satified condition
	if ((i < anchor_n) && check_filter(ref_pos[i].ts, ref_pos[anchor_n_new - 1].ts, ref_temp[i].key, ref_temp[anchor_n_new - 1].key, query_pos[anchor_n_new - 1].qe - query_pos[anchor_n_new - 1].qs))
    {
        anchor_n_new = 0;
        goto MULFIND;
    }

	while(i < anchor_n)
	{
		key = ref_temp[i].key;
		if (key == -1)
		{
			i++;
			continue;
		}
		index = i;
		i++;
		while ((i < anchor_n) && (ref_temp[i].key == key))
		{
			//if the gap less than thre, then merge, else seperate, method 2
            tlen = ref_pos[i].ts - ref_pos[i - 1].te;
            qlen = query_pos[i].qs - query_pos[i - 1].qe;
            thre = tlen*opt->error_overall*opt->error_ins*2;
            if ((abs)(tlen - qlen) > thre) //judge whether can be merge by tlen - qlen < thre
            {
                break;
            }
			i++;
		}

		int a = (int)(query_pos[index].qs - query_pos[anchor_n_new - 1].qe);
		int b = (int)(ref_pos[index].ts - ref_pos[anchor_n_new - 1].te);
		if (a <= 0 || b <= 0)
		{
            uint8_t shift;
			if (a <= 0)
				shift = abs(a) + 1;
			else
				shift = abs(b) + 1;

			query_pos[anchor_n_new - 1].qe -= shift;
			ref_pos[anchor_n_new - 1].te -= shift;

			query_pos[anchor_n_new].qs = query_pos[index].qs + shift;
			query_pos[anchor_n_new].qe = query_pos[i - 1].qe;
			ref_pos[anchor_n_new].ts = ref_pos[index].ts + shift;
			ref_pos[anchor_n_new].te = ref_pos[i - 1].te;

			ref_temp[anchor_n_new].ts = anchor_map2ref[key].ts;
			ref_temp[anchor_n_new].te = anchor_map2ref[key].te;
			ref_temp[anchor_n_new].key = key;
			anchor_n_new++;
		}
		else
		{
			query_pos[anchor_n_new].qs = query_pos[index].qs;
			query_pos[anchor_n_new].qe = query_pos[i - 1].qe;
			ref_pos[anchor_n_new].ts = ref_pos[index].ts;
			ref_pos[anchor_n_new].te = ref_pos[i - 1].te;
			ref_temp[anchor_n_new].ts = anchor_map2ref[key].ts;
			ref_temp[anchor_n_new].te = anchor_map2ref[key].te;
			ref_temp[anchor_n_new].key = key;
			anchor_n_new++;
		}	
	}
    //remove the last anchor if not satified condition
	if((anchor_n_new > 1) && check_filter(ref_pos[anchor_n_new - 1].te, ref_pos[anchor_n_new - 2].te, ref_temp[anchor_n_new - 1].key, ref_temp[anchor_n_new - 2].key, query_pos[anchor_n_new - 1].qe - query_pos[anchor_n_new - 1].qs))
        anchor_n_new--;
	
	*new_n = anchor_n_new;
	return 1;
}

static int align_core(void *km, TARGET_t *anchor_map2ref, uint32_t read_line, uint8_t tid, uint8_t strand, uint8_t *qseq0[2], uint8_t *qual0[2], _aln_t *aln, uint32_t anchor_n, int primary, param_map *opt, ksw_extz_t *ez, ksw_extz_t *ez2)
{
	REF_t *ref_temp = (REF_t* )calloc(anchor_n, sizeof(REF_t));
	uint32_t anchor_n_new = 0;
	reset_aln_t(aln);
	REF_t *ref_pos = REF_pos[tid];
	QUERY_t *query_pos = QUERY_pos[tid];
    

	int signal = find_merge_anchor(anchor_map2ref, ref_temp, ref_pos, query_pos, anchor_n, primary, &anchor_n_new, opt); 
    
	if ((anchor_n_new > 0) && (anchor_n_new < opt->max_exon_num_per_read))
	{
		uint32_t i, j;
		uint8_t which_strand;
		int key;
		int key1, key2;
		int key_l, key_r;
		seq_io *s_io = &seqio[read_line];
        char *seqname = s_io->name;
		uint32_t seqlen = s_io->read_length;
		char *read = s_io->read_seq;
		int chr_n = 0;
        

		uint32_t left_bound;
		chr_n = chromosome_judge(ref_pos[0].ts, &left_bound);
		uint32_t right_bound = chr_end_n[chr_n];
		key1 = ref_temp[0].key;
		key2 = ref_temp[anchor_n_new - 1].key;
		key_l = (key1 > e_shift)? (key1 - e_shift) : 0;
		key_r = (key2 + e_shift < merge_anchor_cnt)? (key2 + e_shift) : (merge_anchor_cnt - 1);
		int key_total = key_r - key_l + 1;
		TARGET_t *target_tmp_FOR = (TARGET_t* )malloc(key_total * sizeof(TARGET_t));
		TARGET_t *target_tmp_REV = (TARGET_t* )malloc(key_total * sizeof(TARGET_t));
        uint8_t direct_RNA = opt->transcript_strand;
        if (!direct_RNA)
        {
            if (anchor_map2ref[key1].strand == 0 && anchor_map2ref[key2].strand == 0 && anchor_map2ref[(key2+key1)>>1].strand == 0)
            {
                // fprintf(stderr, "froward\n");
                for(i = key_l; i <= key_r; ++i)
                {
                    target_tmp_FOR[i - key_l].ts = EXON_T[i].ts_f;
                    target_tmp_FOR[i - key_l].te = EXON_T[i].te_f;
                }
                for (i = 0; i < anchor_n_new; ++i)
                {
                    key = ref_temp[i].key - key_l;
                    ref_temp[i].ts = target_tmp_FOR[key].ts;
                    ref_temp[i].te = target_tmp_FOR[key].te;
                    ref_temp[i].key = key;
                }

                align_core_primary(km, seqlen, target_tmp_FOR, qseq0, qual0, &aln[0], opt, ez, ez2, ref_pos, ref_temp, query_pos, &chr_n, anchor_n_new, strand, tid, key_total, MM_F_SPLICE_FOR, left_bound, right_bound);
                
                if (aln[0].dp_score <= 0)
                {
                    free(ref_temp);
                    free(target_tmp_FOR);
                    free(target_tmp_REV);
                    return 1;
                }
                for(i = 0; i < anchor_n_new; ++i)
                {
                    key = ref_temp[i].key + key_l;
                    if (strand_arr[tid][key] != 1)
                        strand_arr[tid][key] = 0;
                    else
                        strand_arr[tid][key] = 3;
                }
                which_strand = 0;
            }
            else if (anchor_map2ref[key1].strand == 1 && anchor_map2ref[key2].strand == 1 && anchor_map2ref[(key2+key1)>>1].strand == 1)
            {
                // fprintf(stderr, "reverse\n");
                for(i = key_l; i <= key_r; ++i)
                {
                    target_tmp_REV[i - key_l].ts = EXON_T[i].ts_r;
                    target_tmp_REV[i - key_l].te = EXON_T[i].te_r;
                }
                for (i = 0; i < anchor_n_new; ++i)
                {
                    key = ref_temp[i].key - key_l;
                    ref_temp[i].ts = target_tmp_REV[key].ts;
                    ref_temp[i].te = target_tmp_REV[key].te;
                    ref_temp[i].key = key;
                }
                align_core_primary(km, seqlen, target_tmp_REV, qseq0, qual0, &aln[1], opt, ez, ez2, ref_pos, ref_temp, query_pos, &chr_n, anchor_n_new, strand, tid, key_total, MM_F_SPLICE_REV, left_bound, right_bound);
                
                if (aln[1].dp_score <= 0)
                {
                    free(ref_temp);
                    free(target_tmp_FOR);
                    free(target_tmp_REV);
                    return 1;
                }
                for(i = 0; i < anchor_n_new; ++i)
                {
                    key = ref_temp[i].key + key_l;
                    if (strand_arr[tid][key] != 0)
                        strand_arr[tid][key] = 1;
                    else
                        strand_arr[tid][key] = 3;
                }
                which_strand = 1;
            }
            else //have not been detected
            {
                for(i = key_l; i <= key_r; ++i)
                {
                    target_tmp_FOR[i - key_l].ts = EXON_T[i].ts_f;
                    target_tmp_FOR[i - key_l].te = EXON_T[i].te_f;
                }

                for (i = 0; i < anchor_n_new; ++i)
                {
                    key = ref_temp[i].key - key_l;
                    ref_temp[i].ts = target_tmp_FOR[key].ts;
                    ref_temp[i].te = target_tmp_FOR[key].te;
                    ref_temp[i].key = key;
                }

                
                align_core_primary(km, seqlen, target_tmp_FOR, qseq0, qual0, &aln[0], opt, ez, ez2, ref_pos, ref_temp, query_pos, &chr_n, anchor_n_new, strand, tid, key_total, MM_F_SPLICE_FOR, left_bound, right_bound);

                //reverse
                for(i = key_l; i <= key_r; ++i)
                {
                    target_tmp_REV[i - key_l].ts = EXON_T[i].ts_r;
                    target_tmp_REV[i - key_l].te = EXON_T[i].te_r;
                }

                for (i = 0; i < anchor_n_new; ++i)
                {
                    ref_temp[i].ts = target_tmp_REV[ref_temp[i].key].ts;
                    ref_temp[i].te = target_tmp_REV[ref_temp[i].key].te;
                }

                
                align_core_primary(km, seqlen, target_tmp_REV, qseq0, qual0, &aln[1], opt, ez, ez2, ref_pos, ref_temp, query_pos, &chr_n, anchor_n_new, strand, tid, key_total, MM_F_SPLICE_REV, left_bound, right_bound);

                which_strand = (aln[0].dp_score < aln[1].dp_score)? 1:0;
                if (aln[which_strand].dp_score <= 0)
                {
                    free(ref_temp);
                    free(target_tmp_FOR);
                    free(target_tmp_REV);
                    return 1;
                }
                //update anchor_map2ref
                int gap = opt->gap_open_D + opt->gap_ex_D;
                int thre = opt->strand_diff;
                //int thre = 30;
                //
                if ((aln[which_strand].dp_score - aln[1 - which_strand].dp_score) > thre)
                {
                    if (which_strand)
                    {
                        for (i = 0; i < anchor_n_new; ++i)
                        {
                            key = ref_temp[i].key + key_l;
                            if (strand_arr[tid][key] != 0)
                                strand_arr[tid][key] = 1;
                            else
                                strand_arr[tid][key] = 3;
                        }
                    }
                    else
                    {
                        for (i = 0; i < anchor_n_new; ++i)
                        {
                            key = ref_temp[i].key + key_l;
                            if (strand_arr[tid][key] != 1)
                                strand_arr[tid][key] = 0;
                            else
                                strand_arr[tid][key] = 3;
                        }
                    }
                }
            }
        }
        else
        {
            int flag = strand? MM_F_SPLICE_REV : MM_F_SPLICE_FOR;
            if (flag & MM_F_SPLICE_FOR)
            {
                for(i = key_l; i <= key_r; ++i)
                {
                    target_tmp_FOR[i - key_l].ts = EXON_T[i].ts_f;
                    target_tmp_FOR[i - key_l].te = EXON_T[i].te_f;
                }
                for (i = 0; i < anchor_n_new; ++i)
                {
                    key = ref_temp[i].key - key_l;
                    ref_temp[i].ts = target_tmp_FOR[key].ts;
                    ref_temp[i].te = target_tmp_FOR[key].te;
                    ref_temp[i].key = key;
                }

                align_core_primary(km, seqlen, target_tmp_FOR, qseq0, qual0, &aln[0], opt, ez, ez2, ref_pos, ref_temp, query_pos, &chr_n, anchor_n_new, strand, tid, key_total, MM_F_SPLICE_FOR, left_bound, right_bound);
                
                if (aln[0].dp_score <= 0)
                {
                    free(ref_temp);
                    free(target_tmp_FOR);
                    free(target_tmp_REV);
                    return 1;
                }
                for(i = 0; i < anchor_n_new; ++i)
                {
                    key = ref_temp[i].key + key_l;
                    if (strand_arr[tid][key] != 1)
                        strand_arr[tid][key] = 0;
                    else
                        strand_arr[tid][key] = 3;
                }
                which_strand = 0;
            }
            else
            {
                for(i = key_l; i <= key_r; ++i)
                {
                    target_tmp_REV[i - key_l].ts = EXON_T[i].ts_r;
                    target_tmp_REV[i - key_l].te = EXON_T[i].te_r;
                }
                for (i = 0; i < anchor_n_new; ++i)
                {
                    key = ref_temp[i].key - key_l;
                    ref_temp[i].ts = target_tmp_REV[key].ts;
                    ref_temp[i].te = target_tmp_REV[key].te;
                    ref_temp[i].key = key;
                }
                align_core_primary(km, seqlen, target_tmp_REV, qseq0, qual0, &aln[1], opt, ez, ez2, ref_pos, ref_temp, query_pos, &chr_n, anchor_n_new, strand, tid, key_total, MM_F_SPLICE_REV, left_bound, right_bound);
                
                if (aln[1].dp_score <= 0)
                {
                    free(ref_temp);
                    free(target_tmp_FOR);
                    free(target_tmp_REV);
                    return 1;
                }
                for(i = 0; i < anchor_n_new; ++i)
                {
                    key = ref_temp[i].key + key_l;
                    if (strand_arr[tid][key] != 0)
                        strand_arr[tid][key] = 1;
                    else
                        strand_arr[tid][key] = 3;
                }
                which_strand = 1; 
            }
        }
		
		aln[which_strand].chr_n = chr_n;
		uint32_t chr_begin = chr_end_n[chr_n - 1];
		aln[which_strand]._1_based_pos = aln[which_strand]._1_based_pos - chr_begin + 1 + 1; // change to 1_based pos

		free(target_tmp_FOR);
		free(target_tmp_REV);
	}
	else
	{
		free(ref_temp);
		return 1;
	}
	
	free(ref_temp);
	return 0;
}

static void splice_site_judge2(TARGET_t *anchor_map2ref, EXON_t *EXON_T, uint32_t merge_cnt, uint8_t offset)
{
	uint32_t i, j;
    //for(i = 0; i < merge_cnt; ++i)
    //{
    //    EXON_T[i].ts_f = anchor_map2ref[i].ts;
    //    EXON_T[i].te_f = anchor_map2ref[i].te;
    //    EXON_T[i].ts_r = anchor_map2ref[i].ts;
    //    EXON_T[i].te_r = anchor_map2ref[i].te;
    //}

    //return ;


	uint32_t donor_anchor, donor_start;
	uint32_t acceptor_anchor, acceptor_start;
	
	uint8_t detected_len = 2*offset + 10; //20bp 
	uint8_t offset_l, offset_r;
	int16_t site_FOR = 0;
	int16_t site_REV = 0;

	offset_l = 5; 
	offset_r = 2*offset - offset_l;

	uint8_t *ref = (uint8_t* )calloc(detected_len+2, 4);

	//for every site, judge the splice junction site from upstream 5bp and downstream 5bp
	//the first
	acceptor_anchor = anchor_map2ref[0].ts - 1;
	acceptor_start = (acceptor_anchor < offset_r + 4)? 0 : (acceptor_anchor - offset_r - 4);
	donor_anchor = anchor_map2ref[0].te + 1;
	donor_start = donor_anchor - offset_l - 5;

	//detect acceptor
	get_refseq(ref, detected_len, acceptor_start);
    
	acceptor_signals_detected(ref, detected_len, acceptor_start, &site_FOR, &site_REV, 2); 
	if(site_FOR != -1)
		EXON_T[0].ts_f = acceptor_start + site_FOR + 1;
	else
		EXON_T[0].ts_f = anchor_map2ref[0].ts;

	if(site_REV != -1)
		EXON_T[0].ts_r = acceptor_start + site_REV + 1;
	else
		EXON_T[0].ts_r = anchor_map2ref[0].ts;

	//detec donor
	get_refseq(ref, detected_len, donor_start);
	donor_signals_detected(ref, detected_len, donor_start, &site_FOR, &site_REV, 2);

	if (site_FOR != -1)
		EXON_T[0].te_f = donor_start + site_FOR - 1;
	else
		EXON_T[0].te_f = anchor_map2ref[0].te;

	if (site_REV != -1)
		EXON_T[0].te_r = donor_start + site_REV - 1;
	else
		EXON_T[0].te_r = anchor_map2ref[0].te;
	//other
	for (i = 1; i < merge_cnt; ++i)
	{
		acceptor_anchor = anchor_map2ref[i].ts - 1;
		acceptor_start = acceptor_anchor - offset_r - 4;
		donor_anchor = anchor_map2ref[i].te + 1;
		donor_start = donor_anchor - offset_l - 5;

		//detect acceptor
		get_refseq(ref, detected_len, acceptor_start);
		acceptor_signals_detected(ref, detected_len, acceptor_start, &site_FOR, &site_REV, 2);
		if(site_FOR != -1)
			EXON_T[i].ts_f = acceptor_start + site_FOR + 1;
		else
			EXON_T[i].ts_f = anchor_map2ref[i].ts;

		if(site_REV != -1)
			EXON_T[i].ts_r = acceptor_start + site_REV + 1;
		else
			EXON_T[i].ts_r = anchor_map2ref[i].ts;

		//detec donor
		get_refseq(ref, detected_len, donor_start);
		donor_signals_detected(ref, detected_len, donor_start, &site_FOR, &site_REV, 2);

		if (site_FOR != -1)
			EXON_T[i].te_f = donor_start + site_FOR - 1;
		else
			EXON_T[i].te_f = anchor_map2ref[i].te;

		if (site_REV != -1)
			EXON_T[i].te_r = donor_start + site_REV - 1;
		else
			EXON_T[i].te_r = anchor_map2ref[i].te;

        if (EXON_T[i].ts_f <= EXON_T[i - 1].te_f)
        {
            EXON_T[i].ts_f = anchor_map2ref[i].ts;
            EXON_T[i - 1].te_f = anchor_map2ref[i - 1].te;
        }
        if (EXON_T[i].ts_r <= EXON_T[i - 1].te_r)
        {
            EXON_T[i].ts_r = anchor_map2ref[i].ts;
            EXON_T[i - 1].te_r = anchor_map2ref[i - 1].te;
        }
	}

	free(ref);
}

static void splice_site_judge_with_gtf(TARGET_t *anchor_map2ref, EXON_t *EXON_T, uint32_t merge_cnt, uint8_t offset)
{
	uint32_t i, j;
	uint32_t donor_anchor, donor_start;
	uint32_t acceptor_anchor, acceptor_start;
	
	uint8_t detected_len = 2*offset + 10; //20bp 
	uint8_t offset_l, offset_r;
	int16_t site_FOR = 0;
	int16_t site_REV = 0;

	offset_l = 5; 
	offset_r = 2*offset - offset_l;

	uint8_t *ref = (uint8_t* )calloc(detected_len+2, 4);

	//the first
	acceptor_anchor = anchor_map2ref[0].ts - 1;
	acceptor_start = (acceptor_anchor < offset_r + 4)? 0 : (acceptor_anchor - offset_r - 4);
	donor_anchor = anchor_map2ref[0].te + 1;
	donor_start = donor_anchor - offset_l - 5;
	//detect acceptor
	get_refseq(ref, detected_len, acceptor_start);
	acceptor_signals_detected(ref, detected_len, acceptor_start, &site_FOR, &site_REV, 2); 
	if(site_FOR != -1)
		EXON_T[0].ts_f = acceptor_start + site_FOR + 1;
	else
		EXON_T[0].ts_f = anchor_map2ref[0].ts;

	if(site_REV != -1)
		EXON_T[0].ts_r = acceptor_start + site_REV + 1;
	else
		EXON_T[0].ts_r = anchor_map2ref[0].ts;

	//detec donor
	get_refseq(ref, detected_len, donor_start);
	donor_signals_detected(ref, detected_len, donor_start, &site_FOR, &site_REV, 2);

	if (site_FOR != -1)
		EXON_T[0].te_f = donor_start + site_FOR - 1;
	else
		EXON_T[0].te_f = anchor_map2ref[0].te;

	if (site_REV != -1)
		EXON_T[0].te_r = donor_start + site_REV - 1;
	else
		EXON_T[0].te_r = anchor_map2ref[0].te;

	if (anchor_map2ref[0].strand == 0)
	{
		EXON_T[0].ts_f = anchor_map2ref[0].ts;
		EXON_T[0].te_f = anchor_map2ref[0].te;
	}
	else if (anchor_map2ref[0].strand == 1)
	{
		EXON_T[0].ts_r = anchor_map2ref[0].ts;
		EXON_T[0].te_r = anchor_map2ref[0].te;
	}

	for (i = 1; i < merge_cnt; ++i)
	{	
		acceptor_anchor = anchor_map2ref[i].ts - 1;
		acceptor_start = acceptor_anchor - offset_r - 4;
		donor_anchor = anchor_map2ref[i].te + 1;
		donor_start = donor_anchor - offset_l - 5;

		//detect acceptor
		get_refseq(ref, detected_len, acceptor_start);
		acceptor_signals_detected(ref, detected_len, acceptor_start, &site_FOR, &site_REV, 2);
		if(site_FOR != -1)
			EXON_T[i].ts_f = acceptor_start + site_FOR + 1;
		else
			EXON_T[i].ts_f = anchor_map2ref[i].ts;

		if(site_REV != -1)
			EXON_T[i].ts_r = acceptor_start + site_REV + 1;
		else
			EXON_T[i].ts_r = anchor_map2ref[i].ts;

		//detec donor
		get_refseq(ref, detected_len, donor_start);
		donor_signals_detected(ref, detected_len, donor_start, &site_FOR, &site_REV, 2);

		if (site_FOR != -1)
			EXON_T[i].te_f = donor_start + site_FOR - 1;
		else
			EXON_T[i].te_f = anchor_map2ref[i].te;

		if (site_REV != -1)
			EXON_T[i].te_r = donor_start + site_REV - 1;
		else
			EXON_T[i].te_r = anchor_map2ref[i].te;

		if (anchor_map2ref[i].strand == 0)
		{
			EXON_T[i].ts_f = anchor_map2ref[i].ts;
			EXON_T[i].te_f = anchor_map2ref[i].te;
		}
		else if (anchor_map2ref[i].strand == 1)
		{
			EXON_T[i].ts_r = anchor_map2ref[i].ts;
			EXON_T[i].te_r = anchor_map2ref[i].te;
		}
        
        if (EXON_T[i].ts_f <= EXON_T[i - 1].te_f)
        {
            EXON_T[i].ts_f = anchor_map2ref[i].ts;
            EXON_T[i - 1].te_f = anchor_map2ref[i - 1].te;
        }
        if (EXON_T[i].ts_r <= EXON_T[i - 1].te_r)
        {
            EXON_T[i].ts_r = anchor_map2ref[i].ts;
            EXON_T[i - 1].te_r = anchor_map2ref[i - 1].te;
        }
	}

	free(ref);
}

void copy_aln_value(_aln_t *aln1, _aln_t *aln2, uint32_t read_line)
{
	int i;
	aln1->dp_score = aln2->dp_score;
	aln1->dp_max = aln2->dp_max;
	aln1->dp_max2 = aln2->dp_max2;
	aln1->n_ambi = aln2->n_ambi;
	aln1->mlen = aln2->mlen;
	aln1->blen = aln2->blen;
	aln1->chr_n = aln2->chr_n;
	aln1->flag = aln2->flag;
	aln1->_1_based_pos = aln2->_1_based_pos;
	// aln1->mapq = aln2->mapq;
	aln1->n_cigar = aln2->n_cigar;
	//copy cigar
	aln1->cigar = (uint32_t* )calloc(aln1->n_cigar, 4);
	for(i = 0; i < aln1->n_cigar; ++i)
	{
		aln1->cigar[i] = aln2->cigar[i];
	}
}

static void load_query_from_1pass(void *km, TARGET_t * anchor_map2ref, FILE *fp_tff, _aln_t **aln, ksw_extz_t *ez, ksw_extz_t *ez2, uint8_t tid, uint32_t seqn, param_map *opt)
{
	//read anchor info from all_anchor file, process one by one
	uint32_t i;
	uint32_t j;
	uint32_t r_i;
	uint32_t multi_n;
	uint32_t strand;
	uint32_t anchor_n;
	uint32_t primary;
	uint32_t seqi = 0;
	uint32_t seqlen;
	uint8_t splice_offset = 5;
	int max_dp;
	uint8_t max_dp_index;
	uint8_t max_dp_strand;
	uint8_t which_strand;
	uint8_t real_multi_n = 0;
	char* read;
	seq_io *s_io;
	uint8_t return_sig;

	uint8_t *qseq0[2], *qual0[2];
	qseq0[0] = (uint8_t* )kmalloc(km, readlen_max*2);
    qseq0[1] = qseq0[0] + readlen_max;

	if (opt->with_qual) 
	{
		qual0[0] = (uint8_t* )kmalloc(km, readlen_max*2);
    	qual0[1] = qual0[0] + readlen_max;
	}
	

	// rewind(fp_tff);
	while(!feof(fp_tff) && (seqi < seqn))
	{
		real_multi_n = 0;
		s_io = &seqio[seqi];
		seqlen = s_io->read_length;
		read = s_io->read_seq;

		for ( i = 0; i < seqlen; ++i)
		{
			qseq0[0][i] = nt_table[(uint8_t)read[i]];
			if (qseq0[0][i] >= 4)
				qseq0[0][i] = rand()%4;
			qseq0[1][seqlen - 1 - i] = 3 - qseq0[0][i];
		}
		if (opt->with_qual && s_io->qual) 
		{
			for (i = 0; i < seqlen; ++i)
				qual0[0][i] = qual0[1][seqlen - 1 - i] = s_io->qual[i] - 33;
		} else qual0[0] = qual0[1] = 0;
		
		max_dp = KSW_NEG_INF;
		max_dp_index = 0;
		max_dp_strand = 0;
		fscanf(fp_tff, "%u\t", &multi_n);

		if (multi_n == 0)
		{
			//unmapped read
			seqio[seqi].mapable = 0;
			seqio[seqi].aln = NULL;
            seqio[seqi].multi_n = 0;
			seqi++;
			continue;
		}

		for(j = 0; j < multi_n; ++j)
		{
			fscanf(fp_tff, "%u\t%u\t%u\t", &strand, &primary, &anchor_n);
		
			for (i = 0; i < anchor_n; ++i)
			{
				fscanf(fp_tff, "%u\t%u\t%u\t%u\t", &REF_pos[tid][i].ts, &REF_pos[tid][i].te, &QUERY_pos[tid][i].qs, &QUERY_pos[tid][i].qe);
			}

			return_sig = align_core(km, anchor_map2ref, seqi, 0, strand, qseq0, qual0, aln[j], anchor_n, primary, opt, ez, ez2);
			if (return_sig)
			{
				aln[j][0].flag = 4;
				aln[j][1].flag = 4;
			}
			else
			{
				real_multi_n++;
				which_strand = (aln[j][0].dp_score < aln[j][1].dp_score)? 1:0;
				aln[j][which_strand].flag = strand? 16:0;
				int dp = aln[j][which_strand].dp_max;
				if (dp > max_dp)
				{
					max_dp = dp;
					max_dp_index = j;
					max_dp_strand = which_strand;
				}
				else if (dp == max_dp)
				{
					int NM = aln[j][which_strand].blen - aln[j][which_strand].mlen + aln[j][which_strand].n_ambi;
					int NM2 = aln[max_dp_index][max_dp_strand].blen - aln[max_dp_index][max_dp_strand].mlen + aln[max_dp_index][max_dp_strand].n_ambi;
					if (NM < NM2)
					{
						max_dp_index = j;
						max_dp_strand = which_strand;
					}
				}
			}				
		}
		if (aln[max_dp_index][max_dp_strand].dp_max <= 0)
		{
			seqio[seqi].mapable = 0;
			seqio[seqi].aln = NULL;
            seqio[seqi].multi_n = 0;
			seqi++;
			continue;
		}
		//find the best alignment
		//write the best alignment to sam file
		if (real_multi_n > 0)
		{
			seqio[seqi].mapable = 1;
			//record alignment result for print SAM
			uint8_t m = 0;
			seqio[seqi].aln = (_aln_t* )calloc(multi_n, sizeof(_aln_t));
			copy_aln_value(&seqio[seqi].aln[m++], &aln[max_dp_index][max_dp_strand], seqi);
			int32_t tmp_score = 0;
			for(j = 0; j < multi_n; ++j)
			{
				if ((j != max_dp_index) && ((aln[j][0].flag != 4) || (aln[j][1].flag != 4)))
				{
					which_strand = (aln[j][0].dp_score < aln[j][1].dp_score)? 1:0;
					if (tmp_score < aln[j][which_strand].dp_max)
						tmp_score = aln[j][which_strand].dp_max;
					aln[j][which_strand].flag |= 0x100;
					copy_aln_value(&seqio[seqi].aln[m++], &aln[j][which_strand], seqi);
				}
			}
			assert(m==real_multi_n);
			seqio[seqi].multi_n = m;
			seqio[seqi].mapable = 1;
			seqio[seqi].mapq = 60 * (aln[max_dp_index][max_dp_strand].dp_max - tmp_score)/aln[max_dp_index][max_dp_strand].dp_max;

			//reverse read sequence according strand
			if (seqio[seqi].aln[0].flag & 16)
			{
				for(j = 0; j < seqlen; ++j)
				{
					seqio[seqi].read_seq[j] = Dna5Tochar[qseq0[1][j]];
				}
			}
		}
		else
		{
			seqio[seqi].mapable = 0;
			seqio[seqi].aln = NULL;
			seqio[seqi].multi_n = 0;
		}			
		seqi++;
	}

	kfree(km, qseq0[0]);
	if (opt->with_qual) kfree(km, qual0[0]);
}

static int aln_main(uint32_t read_line, thread_2pass_t *aux, dpSkeleton_t *dp_skeleton)
{
	uint8_t multi_n = dp_skeleton->multi_n;
	if (multi_n == 0)
	{
		seqio[read_line].mapable = 0;
		seqio[read_line].aln = NULL;
        seqio[read_line].multi_n = 0;
		return 1;
	}

	int i, j;
	uint32_t anchor_n;
	uint8_t strand;
	uint8_t primary;
	uint8_t return_sig;
	uint8_t real_multi_n = 0;
	uint8_t tid = aux->tid;
	seq_io *s_io = &seqio[read_line];
	uint32_t seqlen = s_io->read_length;
	char* read = s_io->read_seq;
	//copy
	void *km = aux->km;
	TARGET_t *anchor_map2ref = aux->anchor_map2ref;

	int max_dp = KSW_NEG_INF;
	uint8_t max_dp_index = 0;
	uint8_t max_dp_strand = 0;
	uint8_t which_strand;

	uint8_t *qseq0[2], *qual0[2];
	qseq0[0] = (uint8_t* )kmalloc(km, seqlen*2);
    qseq0[1] = qseq0[0] + seqlen;
	for ( i = 0; i < seqlen; ++i)
    {
        qseq0[0][i] = nt_table[(uint8_t)read[i]];
		if (qseq0[0][i] >= 4)
			qseq0[0][i] = rand()%4;
		
		qseq0[1][seqlen - 1 - i] = 3 - qseq0[0][i];
    }
    if (s_io->qual && aux->map->with_qual) 
    {
		qual0[0] = (uint8_t* )kmalloc(km, seqlen*2);
    	qual0[1] = qual0[0] + seqlen;
        for (i = 0; i < seqlen; ++i)
            qual0[0][i] = qual0[1][seqlen - 1 - i] = s_io->qual[i] - 33;
    } else qual0[0] = qual0[1] = 0;

	_aln_t **aln;
	aln = (_aln_t** )calloc(multi_n, sizeof(_aln_t* ));//for multiple alignment
	for(i = 0; i < multi_n; ++i)
	{
		aln[i] = (_aln_t* )calloc(2, sizeof(_aln_t));
		aln[i][0].cigar = (uint32_t* )calloc(seqlen<<1, 4);
		aln[i][1].cigar = (uint32_t* )calloc(seqlen<<1, 4);
	}
	for(i = 0; i < multi_n; ++i)
	{
		anchor_n = dp_skeleton->point[i].anchor_n;
		primary = dp_skeleton->point[i].primary;
		strand = dp_skeleton->point[i].strand;

		for(j = 0; j < anchor_n; ++j)
		{
			REF_pos[tid][j].ts = dp_skeleton->point[i].anchor_pos[j][0];
			REF_pos[tid][j].te = dp_skeleton->point[i].anchor_pos[j][1];
			QUERY_pos[tid][j].qs = dp_skeleton->point[i].anchor_pos[j][2];
			QUERY_pos[tid][j].qe = dp_skeleton->point[i].anchor_pos[j][3];
		}
		
		return_sig = align_core(km, anchor_map2ref, read_line, tid, strand, qseq0, qual0, aln[i], anchor_n, primary, aux->map, &aux->ez, &aux->ez2);
		if (return_sig)
		{
			aln[i][0].flag = 4;
			aln[i][1].flag = 4;
		}
		else
		{
			real_multi_n++;
			which_strand = (aln[i][0].dp_score < aln[i][1].dp_score)? 1:0;
			aln[i][which_strand].flag = strand? 16:0;
			int dp = aln[i][which_strand].dp_max;
			if (dp > max_dp)
			{
				max_dp = dp;
				max_dp_index = i;
				max_dp_strand = which_strand;
			}
			else if (dp == max_dp) // && max_dp >
			{
				int NM = aln[i][which_strand].blen - aln[i][which_strand].mlen + aln[i][which_strand].n_ambi;
				int NM2 = aln[max_dp_index][max_dp_strand].blen - aln[max_dp_index][max_dp_strand].mlen + aln[max_dp_index][max_dp_strand].n_ambi;
				if (NM < NM2)
				{
					max_dp_index = i;
					max_dp_strand = which_strand;
				}
			}
		}
	}
	
	if (aln[max_dp_index][max_dp_strand].dp_max <= 0)
	{
		seqio[read_line].mapable = 0;
		seqio[read_line].aln = NULL;
        seqio[read_line].multi_n = 0;
		goto FREE;
	}

	//find the best alignment
	//write the best alignment to sam file
	if (real_multi_n > 0)
	{
		//record alignment result for print SAM
		uint8_t m = 0;
		seqio[read_line].aln = (_aln_t* )calloc(multi_n, sizeof(_aln_t));
		copy_aln_value(&seqio[read_line].aln[m], &aln[max_dp_index][max_dp_strand], read_line);
		m++;
		int32_t tmp_score = 0;
		for(j = 0; j < multi_n; ++j)
		{
			if ((j != max_dp_index) && ((aln[j][0].flag != 4) || (aln[j][1].flag != 4)))
			{
				which_strand = (aln[j][0].dp_score < aln[j][1].dp_score)? 1:0;
				if (tmp_score < aln[j][which_strand].dp_max)
						tmp_score = aln[j][which_strand].dp_max;
				aln[j][which_strand].flag |= 0x100;
				copy_aln_value(&seqio[read_line].aln[m++], &aln[j][which_strand], read_line);
			}
		}
		assert(m==real_multi_n);
		seqio[read_line].multi_n = m;
		seqio[read_line].mapable = 1;
		seqio[read_line].mapq = 60 * (aln[max_dp_index][max_dp_strand].dp_max - tmp_score)/aln[max_dp_index][max_dp_strand].dp_max;
		////reverse read sequence according strand
		if (seqio[read_line].aln[0].flag & 16)
		{
			for(j = 0; j < seqlen; ++j)
			{
				seqio[read_line].read_seq[j] = Dna5Tochar[qseq0[1][j]];
			}
		}
	}
	else
	{
		seqio[read_line].mapable = 0;
		seqio[read_line].aln = NULL;
        seqio[read_line].multi_n = 0;
	}

FREE:
	for (i = 0; i < multi_n; ++i)
	{
		if (aln[i][0].cigar != NULL)	free(aln[i][0].cigar);
		if (aln[i][1].cigar != NULL)	free(aln[i][1].cigar);
		if (aln[i] != NULL)	free(aln[i]);		
	}
	if (aln != NULL)	free(aln);

	kfree(km, qseq0[0]);
	if (s_io->qual && aux->map->with_qual) kfree(km, qual0[0]);

	return 0;
}

static void *aln_main_thread(void *aux)
{
	thread_2pass_t *d = (thread_2pass_t* )aux;
	int _read_lines;

	while(1)
	{
		pthread_rwlock_wrlock(&RWLOCK);
		_read_lines = THREAD_READ_I++;
		pthread_rwlock_unlock(&RWLOCK);

		if (_read_lines < d->seqn)
		{
			aln_main(_read_lines, d, &d->dp_skeleton[_read_lines]);
		}
		else break;
	}
	return 0;
}

void load_fasta_2pass(uint32_t map2ref_cnt, param_map *opt, char *read_fastq, int *Total_mapped_reads)
{
	uint32_t read_in = opt->batch_size;
	uint32_t seqii = read_in;
	uint32_t seqi;
	uint32_t r_i, r_ii;
	bseq_file_t *bf;  //can set to extern
	void *km = 0;
	hash_kmer = opt->hash_kmer;
	e_shift = opt->e_shift;
	clock_t a = clock();

#ifdef PRINT
	thread_n = 1;
#endif
	initHashTable(hash_kmer);

#ifdef HAVE_KALLOC
    km = km_init();
#endif
	bf = bseq_open(read_fastq);
	if(bf == 0)
    {
        fprintf(stderr, "[Wrong] Wrong input file route or name: %s \n", read_fastq);
        exit(1);
    }

	FILE* fp_tff = fopen(temp_anchor_dir, "r");
	if (fp_tff == NULL)
	{
		fprintf(stderr, "[Wrong] Open the temporary file %s failed!!!\n", temp_anchor_dir);
        fprintf(stderr, "[Warring] The error may caused by two deSALT program were running in the same time, but you didn't specify the temporary file path (-f)!\n Please cheak or see more detailed from usage.\n");
		exit(0);
	}

	seqio = (seq_io* )calloc(read_in, sizeof(seq_io));

	REF_pos = (REF_t** )malloc(thread_n*sizeof(REF_t* ));
	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		REF_pos[r_i] = (REF_t* )malloc(max_exon_num_per_read*sizeof(REF_t)); //key item have no use
	}
	QUERY_pos = (QUERY_t** )malloc(thread_n*sizeof(QUERY_t* ));
	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		QUERY_pos[r_i] = (QUERY_t* )malloc(max_exon_num_per_read*sizeof(QUERY_t));
	}
	
	TARGET_t *anchor_map2ref = (TARGET_t* )malloc(map2ref_cnt*sizeof(TARGET_t));

	fprintf(stderr, "[Phase-INFO] Exons inference by skeletons of all reads\n");
	if (opt->with_gtf)
	{
		fprintf(stderr, "[Phase-INFO] Loading GTF annotations\n");
        IvalA = read_Annotation(opt, anchor_map2ref, map2ref_cnt, &merge_anchor_cnt);
	}
	else
    {
		merge_anchor_cnt = load_anchor(anchor_map2ref, map2ref_cnt);
		for (r_i = 0; r_i < merge_anchor_cnt; ++r_i)
		{
			anchor_map2ref[r_i].strand = 3;
			anchor_map2ref[r_i].cov = 0;
		}
	}

	fprintf(stderr, "[Phase-INFO] Inferring total %d isolated regions (pseudo-exons) after merging and filtering\n", merge_anchor_cnt);

	strand_arr = (uint8_t** )calloc(thread_n, sizeof(uint8_t* ));
	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		strand_arr[r_i] = (uint8_t *)malloc(merge_anchor_cnt);
		memset(strand_arr[r_i], 3, merge_anchor_cnt);
	}

	if (merge_anchor_cnt > 0)
	{
		fprintf(stderr, "[Phase-INFO] Refining pseudo-exons by scoring matrix\n");
		EXON_T = (EXON_t* )malloc(merge_anchor_cnt*sizeof(EXON_t));
		//judge the splice junction according to scoring matrix
		if (opt->read_type == 3)
			splice_offset = 10;
		else
			splice_offset = 5;

		if (opt->with_gtf)
			splice_site_judge_with_gtf(anchor_map2ref, EXON_T, merge_anchor_cnt, splice_offset);
		else
			splice_site_judge2(anchor_map2ref, EXON_T, merge_anchor_cnt, splice_offset);
	}

	int time = 0;
	pthread_rwlock_init(&RWLOCK, NULL);
	thread_2pass_t* aux;
	aux = (thread_2pass_t* )calloc(thread_n, sizeof(thread_2pass_t));
	for(r_i = 0; r_i < thread_n; ++r_i)
	{
		aux[r_i].tid = r_i;
		memset(&aux[r_i].ez, 0, sizeof(ksw_extz_t));
		memset(&aux[r_i].ez2, 0, sizeof(ksw_extz_t));
		aux[r_i].map = opt;
#ifdef HAVE_KALLOC
		aux[r_i].km = km_init();
#endif
	} 

	double t_s, t_e;
	while(seqii == read_in)
	{
		t_s = realtime();
		seqii = bseq_read_2pass(bf, read_in, seqio);
		fprintf(stderr, "[Loop-ProcessReads] The %dst loop of refined alignment procedure with %d reads......\n", time, seqii);

		THREAD_READ_I = 0;
		if (thread_n <= 1)
		{
			_aln_t **aln;
			aln = (_aln_t** )calloc(opt->top_n * 2, sizeof(_aln_t* ));//for multiple alignment
			for(r_i = 0; r_i < opt->top_n * 2; ++r_i)
			{
				aln[r_i] = (_aln_t* )calloc(2, sizeof(_aln_t));
				aln[r_i][0].cigar = (uint32_t* )calloc(readlen_max, 4);
				aln[r_i][1].cigar = (uint32_t* )calloc(readlen_max, 4);
			}

            ksw_extz_t ez;
            ksw_extz_t ez2;
            memset(&ez, 0, sizeof(ksw_extz_t));
            memset(&ez2, 0, sizeof(ksw_extz_t));
			
			load_query_from_1pass(km, anchor_map2ref, fp_tff, aln, &ez, &ez2, 0, seqii, opt);

			for (r_i = 0; r_i < opt->top_n * 2; ++r_i)
			{
				if (aln[r_i][0].cigar != NULL)	free(aln[r_i][0].cigar);
				if (aln[r_i][1].cigar != NULL)	free(aln[r_i][1].cigar);
				if (aln[r_i] != NULL)	free(aln[r_i]);
			}
			if (aln != NULL)	free(aln);
            kfree(km, ez.cigar);
            kfree(km, ez2.cigar);
		}
		else
		{
			int i, j;
			uint32_t multi_n;
			dpSkeleton_t* dp_skeleton;

			dp_skeleton = (dpSkeleton_t *)calloc(seqii, sizeof(dpSkeleton_t));
			// load all anchor pos
			seqi = 0;
			while(!feof(fp_tff) && (seqi < seqii))
			{
				fscanf(fp_tff, "%u\t", &multi_n);
				dp_skeleton[seqi].multi_n = multi_n;
				if (multi_n == 0)
				{
					// //unmapped read
					seqi++;
					continue;
				}

				anchor_t *Anchor = (anchor_t* )calloc(multi_n, sizeof(anchor_t));
				for(i = 0; i < multi_n; ++i)
				{
					fscanf(fp_tff, "%u\t%u\t%u\t", &Anchor[i].strand, &Anchor[i].primary, &Anchor[i].anchor_n);
					Anchor[i].anchor_pos = (uint32_t** )calloc(Anchor[i].anchor_n, sizeof(uint32_t* ));
					for(j = 0; j < Anchor[i].anchor_n; ++j)
					{
						Anchor[i].anchor_pos[j] = (uint32_t* )calloc(4, 4);
					}
					for (j = 0; j < Anchor[i].anchor_n; ++j)
					{
						fscanf(fp_tff, "%u\t%u\t%u\t%u\t", &Anchor[i].anchor_pos[j][0], &Anchor[i].anchor_pos[j][1], &Anchor[i].anchor_pos[j][2], &Anchor[i].anchor_pos[j][3]);
					}
				}
				dp_skeleton[seqi].point = Anchor;
				seqi++;
			}

			pthread_t* tid;
			pthread_attr_t attr;

			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			
			tid = (pthread_t* )calloc(thread_n, sizeof(pthread_t));

			for(r_i = 0; r_i < thread_n; ++r_i)
			{
				aux[r_i].seqn = seqii;
				aux[r_i].dp_skeleton = dp_skeleton;
				aux[r_i].anchor_map2ref = anchor_map2ref; //can remove

				int res = pthread_create(&tid[r_i], &attr, aln_main_thread, aux + r_i);
				if(res != 0)
				{
					fprintf(stderr, "[Wrong] Create pthread error\n");
					exit(1);
				}
			}
			for(r_i = 0; r_i < thread_n; ++r_i)	pthread_join(tid[r_i], 0);
			free(tid);
			for(r_i = 0; r_i < seqii; ++r_i)
			{
				for (r_ii = 0; r_ii < dp_skeleton[r_i].multi_n; ++r_ii)
				{
					if (dp_skeleton[r_i].point[r_ii].anchor_pos != NULL)	free(dp_skeleton[r_i].point[r_ii].anchor_pos);
				}
				if(dp_skeleton[r_i].point != NULL)	free(dp_skeleton[r_i].point);
			}
			if (dp_skeleton != NULL)	free(dp_skeleton);
		}

		//update the information of exon strands
		int forward = 0;
		int reverse = 0;
		int change = 0;
		for(r_i = 0; r_i < merge_anchor_cnt; ++r_i)
		{
			if (anchor_map2ref[r_i].strand != 3)
				continue;

			forward = 0;
			reverse = 0;
			for(r_ii = 0; r_ii < thread_n; ++r_ii)
			{
				if (strand_arr[r_ii][r_i] == 0)
					forward += 1;
				else if (strand_arr[r_ii][r_i] == 1)
					reverse += 1;
			}
			if (forward > reverse)
			{
				anchor_map2ref[r_i].strand = 0;
				change += 1;
			}
			else if (forward < reverse)
			{
				anchor_map2ref[r_i].strand = 1;
				change += 1;
			}
			
		}

		change = 0;
		for(r_i = 0; r_i < merge_anchor_cnt; ++r_i)
		{
			if (anchor_map2ref[r_i].strand != 3)
				change += 1;
		}
		t_e = realtime();
		//output SAM
		int mm = ff_print_sam (seqio, seqii, opt);
		*Total_mapped_reads += mm;
        fprintf(stderr, "[Loop-ProcessReads] Aligned %d reads to genome, and confirmed total %d exons' strand in %f seconds\n", mm, change, t_e - t_s);	


		for(r_i = 0; r_i < seqii; ++r_i)
		{
			for (r_ii = 0; r_ii < seqio[r_i].multi_n; ++r_ii)
				if (seqio[r_i].aln[r_ii].cigar != NULL) free(seqio[r_i].aln[r_ii].cigar);

			if(seqio[r_i].aln != NULL)
			{
				free(seqio[r_i].aln);
				seqio[r_i].aln = NULL;
			}
			if(seqio[r_i].name != NULL)
			{
				free(seqio[r_i].name);
				seqio[r_i].name = NULL;
			}
			if(seqio[r_i].read_seq != NULL)
			{
				free(seqio[r_i].read_seq);
				seqio[r_i].read_seq = NULL;
			}
			if(seqio[r_i].qual != NULL)
			{
				free(seqio[r_i].qual);
				seqio[r_i].qual = NULL;
			}
		}
		time++;
		for(r_i = seqii; r_i < read_in; ++r_i)
		{
			if(seqio[r_i].name != NULL)
			{
				free(seqio[r_i].name);
				seqio[r_i].name = NULL;
			}
			if(seqio[r_i].read_seq != NULL)
			{
				free(seqio[r_i].read_seq);
				seqio[r_i].read_seq = NULL;
			}
			if(seqio[r_i].qual != NULL)
			{
				free(seqio[r_i].qual);
				seqio[r_i].qual = NULL;
			}
		}
	}
	pthread_rwlock_destroy(&RWLOCK);
	for(r_i = 0; r_i < thread_n; ++r_i)
	{
        if(aux[r_i].ez.cigar)   kfree(aux[r_i].km, aux[r_i].ez.cigar);
        if(aux[r_i].ez2.cigar)  kfree(aux[r_i].km, aux[r_i].ez2.cigar);

#ifdef HAVE_KALLOC
		km_destroy(aux[r_i].km);
#endif
	}
	free(aux);	
    
	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		if (strand_arr[r_i] != NULL)	free(strand_arr[r_i]);
	}
	if (strand_arr != NULL)	free(strand_arr);

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		if (QUERY_pos[r_i] != NULL)	free(QUERY_pos[r_i]);
	}
	if(QUERY_pos != NULL)	free(QUERY_pos);

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		if (REF_pos[r_i] != NULL)	free(REF_pos[r_i]);
	}
	if(REF_pos != NULL)	free(REF_pos);

	if(anchor_map2ref != NULL)	free(anchor_map2ref);

	if (EXON_T != NULL)	free(EXON_T);
   
    if (opt->with_gtf)
    {
        for(r_i = 0; r_i < chr_file_n; ++r_i)
        {
            if (IvalA[r_i].Ival != NULL) free(IvalA[r_i].Ival);
        }
        if (IvalA != NULL) free(IvalA);
    }
    

	for (r_i = 0; r_i < read_in; ++r_i)
	{
		if (seqio[r_i].read_seq != NULL)	free(seqio[r_i].read_seq);
		if (seqio[r_i].qual != NULL)	free(seqio[r_i].qual);
		if (seqio[r_i].name != NULL)	free(seqio[r_i].name);
	}
	if(seqio != NULL)	free(seqio);

	freeHashTable();
#ifdef HAVE_KALLOC
	km_destroy(km);
#endif
	fclose(fp_tff);
	bseq_close(bf);
}
