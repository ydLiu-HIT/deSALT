/*************************************************************************
	> File Name: splic_junction.c
	> Author: 
	> Mail: 
	> Created Time: 2019年01月03日 星期四 22时05分51秒
 ************************************************************************/

#include<stdio.h>
#include <inttypes.h>
#include "splic_junction.h"

static const int16_t donor_scoring_matrix[10][4] = { 
    {-92, 168, -114, 103}, {132, 252, -109, -138},
    {-97, 273, -167, 143}, {922, -208, -132, -202},
    {-107, -288, 1104, -225}, {-344, -344, 2406, -344},
    {-344, -344, -344, 2406}, {355, -236, 543, -223},
    {939, -177, -160, -207}, {222, -206, 1145, 185}
};


static const int16_t acceptor_scoring_matrix[10][4] = {
    {-201, 431, -121, 350}, {115, 176, 199, -107},
    {61, 838, -191, -169}, {2548, -364, -364, -364},
    {-364, -364, 2548, -364}, {704, -222, 347, -272},
    {-142, -170, -117, 1030}, {-167, -92, 974, -137},
    {-33, 201, 157, 47}, {12, 485, -91, 9}
};

const int16_t donor_MIN = -2259;
const int16_t donor_MAX = 10158;
const int16_t donor_MAX_MIN = 12417;
const int16_t acceptor_MIN = -1960;
const int16_t acceptor_MAX = 9958;
const int16_t acceptor_MAX_MIN = 11918;


void acceptor_signals_detected(uint8_t* ref, uint16_t length, uint32_t detected_start, int16_t *s1, int16_t *s2, uint8_t strand)
{
    int16_t score, score1;
    int16_t score_max = 0;
	int16_t score_max1 = 0;
	int16_t site_index = 0;
	int16_t site_index1 = 0;
    uint16_t i;
    uint16_t len = length - 9;
    uint8_t j;

	if (strand == 2) //both strand
	{
		for (i = 0; i < len; ++i)
		{
			score = 0;
			score1 = 0;
			for (j = 0; j < 10; ++j)
			{
				score += acceptor_scoring_matrix[j][ref[i + j]];//+
				score1 += donor_scoring_matrix[j][3 - ref[i + 9 - j]]; //-
			}
			if (score > score_max)
			{
				score_max = score;
				site_index = i;
			}
			if(score1 > score_max1)
			{
				score_max1 = score1;
				site_index1 = i;
			}
		}
	}
	else if (strand == 0) //+
	{
		for (i = 0; i < len; ++i)
		{
			score = 0;
			for (j = 0; j < 10; ++j)
			{
				score += acceptor_scoring_matrix[j][ref[i + j]];//+
			}
			if (score > score_max)
			{
				score_max = score;
				site_index = i;
			}
		}
	}
	else if (strand == 1) //-
	{
		for (i = 0; i < len; ++i)
		{
			score1 = 0;
			for (j = 0; j < 10; ++j)
			{
				score1 += donor_scoring_matrix[j][3 - ref[i + 9 - j]]; //-
			}
			if(score1 > score_max1)
			{
				score_max1 = score1;
				site_index1 = i;
			}
		}
	}
	
	if (strand == 2)
	{
		//+
		if ((int)(100 * ((score_max - acceptor_MIN)/(float)acceptor_MAX_MIN)) >= 50)
		{
			site_index += 4; // the intron end position
		}else
		{
			site_index = -1;
		}
		*s1 = site_index;
		//-
		if ((int)(100 * ((score_max1 - donor_MIN)/(float)donor_MAX_MIN)) >= 50)
		{
			site_index1 += 4; // the intron end position
		}else
		{
			site_index1 = -1;
		}
		*s2 = site_index1;
	}
	else if (strand == 0) //+
	{
		//+
		if ((int)(100 * ((score_max - acceptor_MIN)/(float)acceptor_MAX_MIN)) >= 50)
		{
			site_index += 4; // the intron end position
		}else
		{
			site_index = -1;
		}
		*s1 = site_index;
	}
	else if (strand == 1) //-
	{
		//-
		if ((int)(100 * ((score_max1 - donor_MIN)/(float)donor_MAX_MIN)) >= 50)
		{
			site_index1 += 4; // the intron end position
		}else
		{
			site_index1 = -1;
		}
		*s2 = site_index1;
	}
	
}

void donor_signals_detected(uint8_t* ref, uint16_t length, uint32_t detected_start, int16_t *s1, int16_t *s2, uint8_t strand)
{
    int16_t score = 0;
	int16_t score1 = 0;
    int16_t score_max = 0;
	int16_t score_max1 = 0;
    int16_t site_index = 0;
	int16_t site_index1 = 0;
    uint16_t i;
    uint16_t len = length - 9;
    uint8_t j;

	if (strand == 2) //both strand
	{
		for (i = 0; i < len; ++i)
		{
			score = 0;
			score1 = 0;
			for (j = 0; j < 10; ++j)
			{
				score += donor_scoring_matrix[j][ref[i + j]];;//+
				score1 += acceptor_scoring_matrix[j][3 - ref[i + 9 - j]]; //-
			}
			if (score > score_max)
			{
				score_max = score;
				site_index = i;
			}
			if(score1 > score_max1)
			{
				score_max1 = score1;
				site_index1 = i;
			}
		}
	}
	else if (strand == 0) //+
	{
		for (i = 0; i < len; ++i)
		{
			score = 0;
			for (j = 0; j < 10; ++j)
			{
				score += donor_scoring_matrix[j][ref[i + j]];;//+
			}
			if (score > score_max)
			{
				score_max = score;
				site_index = i;
			}
		}
	}
	else if (strand == 1) //-
	{
		for (i = 0; i < len; ++i)
		{
			score1 = 0;
			for (j = 0; j < 10; ++j)
			{
				score1 += acceptor_scoring_matrix[j][3 - ref[i + 9 - j]]; //-
			}
			if(score1 > score_max1)
			{
				score_max1 = score1;
				site_index1 = i;
			}
		}
	}
	
	if (strand == 2)//both
	{
		//+
		if ((int)(100 * ((score_max - donor_MIN)/(float)donor_MAX_MIN)) >= 50)
		{
			site_index += 5; // the intron end position
		}else
		{
			site_index = -1;
		}
		*s1 = site_index;
		//-
		if ((int)(100 * ((score_max1 - acceptor_MIN)/(float)acceptor_MAX_MIN)) >= 50)
		{
			site_index1 += 5; // the intron end position
		}else
		{
			site_index1 = -1;
		}
		*s2 = site_index1;
	}
	else if (strand == 0)
	{
		//+
		if ((int)(100 * ((score_max - donor_MIN)/(float)donor_MAX_MIN)) >= 50)
		{
			site_index += 5; // the intron end position
		}else
		{
			site_index = -1;
		}
		*s1 = site_index;
	}
	else if (strand == 1)
	{
		//-
		if ((int)(100 * ((score_max1 - acceptor_MIN)/(float)acceptor_MAX_MIN)) >= 50)
		{
			site_index1 += 5; // the intron end position
		}else
		{
			site_index1 = -1;
		}
		*s2 = site_index1;
	}
}