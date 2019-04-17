/*************************************************************************
	> File Name: graph.c
	> Author: 
	> Mail: 
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "graph.h"
#include "binarys_qsort.h"
#include "load_unipath_size.h"
#include "bit_operation.h"
#include "aln_2pass.h"
#include "read_seeding.h"


#define REFLEN 0Xffffffff

double DAG_time = 0;
double findpath_time = 0;
double qsort_time = 0;

//global variant
Graph* Dp_graph = NULL;

void initGraph()
{
	uint8_t r_i;
	uint32_t r_ii;
	Dp_graph = (Graph* )calloc(thread_n,sizeof(Graph));
	if(Dp_graph == NULL) fprintf(stderr, "[Wrong] Failed to allocate graph memory!\n");

	for(r_i = 0; r_i < thread_n; ++r_i)
	{
		Dp_graph[r_i].vnode = (VNode* )calloc(new_seed_cnt*pos_n_max,sizeof(VNode));
		if(Dp_graph[r_i].vnode == NULL) printf("wrong for allocate graph memory!\n");
	}
}

void delGraph()
{	
	uint8_t r_i;
	uint32_t r_ii;
	for(r_i = 0; r_i < thread_n; ++r_i)
	{
		for (r_ii = 0; r_ii < new_seed_cnt*pos_n_max; ++r_ii)
		{
			if (Dp_graph[r_i].vnode[r_ii].preedge != NULL)	free(Dp_graph[r_i].vnode[r_ii].preedge);
		}
		if (Dp_graph[r_i].vnode != NULL)	free(Dp_graph[r_i].vnode);
	}
	if (Dp_graph != NULL)	free(Dp_graph);
}


void show_uniseed(uni_seed *uniseed, uint32_t uniseedNum)
{
	uint32_t i;
	for (i = 0; i < uniseedNum; ++i)
	{
		fprintf(stderr,"id = %u, read_begin = %u, read_end = %u, cov = %d, seed_id = %u,ref_begin  = %u, ref_end = %u\n",i, uniseed[i].read_begin, uniseed[i].read_end, uniseed[i].cov, uniseed[i].seed_id, \
			uniseed[i].ref_begin, uniseed[i].ref_end);
	}
}
/*
void show_vertexm(vertex_m *vertexArr, uint32_t vertexNum)
{
	uint32_t i;
	for (i = 0; i < vertexNum; ++i)
	{
		fprintf(stderr,"read_pos = %u, UID = %"PRId64", seed_id = %u, unipos_off = %u, length = %u\n",vertexArr[i].read_pos, vertexArr[i].uid, \
			vertexArr[i].seed_id, vertexArr[i].uni_pos_off, vertexArr[i].length);
		//output sets of ref_pos
	}
}

void show_vertexu(vertex_u *vertexArr, uint32_t vertexNum)
{
	uint32_t i;
	for (i = 0; i < vertexNum; ++i)
	{
		fprintf(stderr,"read_pos = %u, UID = %"PRId64", unipos_off = %u, length1 = %u, length2 = %u\n",vertexArr[i].read_pos, vertexArr[i].uid, \
			vertexArr[i].uni_pos_off, vertexArr[i].length1, vertexArr[i].length2);
		//output sets of ref_pos
	}
}
*/


void dynamic_programming_path(Graph *graph, uint32_t vertexNum, PATH_t *dist_path, uint8_t *out_degree);
int32_t get_longest_path(Graph *graph, uni_seed* vertexArr, uint32_t vertexNum,uint8_t rc_i, uint32_t *fillNum);

static float min(float a, float b)
{
    return a<b?a:b;
}

float creatGraph(uni_seed *vertexArr, uint32_t vertexNum, PATH_t *dist_path, uint8_t *out_degree,uint32_t *max_index, uint8_t tid, int max_read_join_gap)
{
	int32_t i,j;

	uint32_t seed_id = 0;
	int32_t weight = 0;
	uint32_t ref_pos2 = 0;
	uint32_t read_pos2 = 0;

	float max_distance = 0;
    int32_t gap;
    int32_t ove;
    int32_t ove1;
    int32_t dis1, dis2;
    int8_t param = 2;
	int8_t intron_penalty = Eindel*param/seed_k_t;
	uint32_t non_ioslated_point= 0;
	uint32_t thre_num;
	int search_step = (vertexNum < 50)? vertexNum : 50;
	
	pANode anode;
	Graph *graph = &Dp_graph[tid];

	double time1 = clock();
	qsort(vertexArr, vertexNum, sizeof(uni_seed), compare_uniseed);
	qsort_time += (clock() - time1)/CLOCKS_PER_SEC;

	if (vertexNum == 0)
	{
		return max_distance;
	}

	// re inital
	for (i = 0; i < vertexNum; ++i)
	{
        // vertexArr[i].id = i;
        weight = vertexArr[i].cov;
		graph->vnode[i].head_vertex = i;
        graph->vnode[i].weight = weight;
		graph->vnode[i].adjacent_node = 0;
		graph->vnode[i].preedge = (ANode *)calloc(search_step, sizeof(ANode));

    	dist_path[i].dist = weight;
		dist_path[i].pre_node = -1;
	}

#ifdef Annoation
	show_uniseed(vertexArr, vertexNum);	
#endif

	// add for print dist
	memset(out_degree, 0, vertexNum);

	time1 = clock();
	for ( i = 0; i < vertexNum-1; ++i)
	{
		read_pos2 = vertexArr[i].read_end;
		ref_pos2 = vertexArr[i].ref_end;
		seed_id = vertexArr[i].seed_id;

		thre_num = (vertexNum < (i + search_step))? vertexNum : (i + search_step);

		for (j = i + 1; j < thre_num; ++j)
		{
			//two mem generated from the same seed will not be connected.
			if (vertexArr[j].seed_id == seed_id)
			{
				continue;
			}

			if (vertexArr[j].ref_begin > ref_pos2 + max_intron_length)
			{
				break;
			}
        
            ove = (int32_t)(vertexArr[j].read_begin - read_pos2);
            if (ove > max_read_join_gap)
                continue;
            ove1 = (int32_t)(vertexArr[j].ref_begin - ref_pos2);

			if ((ove > 0 && ove1 > 0 && ove1 > ove/2) || (ove >= -5 && ove <= 0  && ove1 >= -5))  //5  1/error rate
            {
				//the first part is normal, if (ove - ove1) < Eindel, there is an edge
				//the last part, beaucse the begin of exon and the begin of intron have the same bases 
            	graph->vnode[j].preedge[graph->vnode[j].adjacent_node].adjvex = i;
                dis1 = vertexArr[j].read_end - read_pos2;
                dis2 = vertexArr[j].ref_end - ref_pos2;
                gap = (int32_t)(dis1 - dis2);
			
				int diff = (read_pos2 >= vertexArr[j].read_begin)? (read_pos2 + 1 - vertexArr[j].read_begin) : 0;
				weight = vertexArr[j].cov - diff;
				
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].weight = min(weight, min(dis1, dis2));
                if (gap >= 0) //in the same exon
					graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty = gap*param/(float)weight;
                else
					graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty = min(abs(gap)*param/(float)weight, intron_penalty);
                graph->vnode[j].adjacent_node++;
				out_degree[i] = 1;

        		non_ioslated_point++;
            }
            else if (ove1 == ove) 
            {
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].adjvex = i;
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].weight = vertexArr[j].cov - (read_pos2 + 1 - vertexArr[j].read_begin);
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty = 0.0;
				graph->vnode[j].adjacent_node++;
				out_degree[i] = 1;

				non_ioslated_point++;
            }
		}
	}
	DAG_time += (clock() - time1)/CLOCKS_PER_SEC;
	if (non_ioslated_point != 0)
	{
		time1 = clock();
		dynamic_programming_path(graph, vertexNum, dist_path, out_degree);
		findpath_time += (clock() - time1)/CLOCKS_PER_SEC;
	}

	for ( i = 0; i < vertexNum; ++i)
	{
		if (max_distance < dist_path[i].dist)
		{
			max_distance = dist_path[i].dist;
			*max_index = i;
		}
	}

	//free
	for (i = 0; i < vertexNum; ++i)
	{
		free(graph->vnode[i].preedge);
		graph->vnode[i].preedge = NULL;
	}
	return max_distance;
}


/*
Dplongestpath(G)
Initialize all dilg(.) values to ∞;
Let S be the set of vertices with indegree=0;
for each vertex v in S do
     dist(v)=0;
4. For each v∈V\S in Topological Sorting order do
       dilg(v)=max(u,v)∈E{dilg(u)+w(u, v)}
     let (u,v) be the edge to get the maximum     
     value;
     dad(v)=u;
5. Return the dilg(.) with maximum value.
*/

uint32_t max(uint32_t a, uint32_t b)
{
    return a>b?a:b;
}

void Find_longest_path(Graph *graph, uint32_t target, PATH_t *dist_path)
{
	float temp = 0;
    float penalty;
	float current_dist = 0;
	int32_t pre_node = 0;
    int32_t ver = 0;
    int32_t weight;
	uint32_t adjvex;

    //travel the processer of vertex target
    pANode arcnode;
	int i;
	for (i = 0; i < graph->vnode[target].adjacent_node; ++i)
	{
		arcnode = &graph->vnode[target].preedge[i];
		weight = arcnode->weight;
		penalty = arcnode->penalty;
		adjvex = arcnode->adjvex;
		temp = dist_path[adjvex].dist + weight - penalty;

		if (current_dist <= temp)
        {
            current_dist = temp;
            pre_node = adjvex;
        }
	}

	dist_path[target].dist = current_dist;   // change the distance because there are overlaps between mems
	dist_path[target].pre_node = pre_node; //the front node
}

/*
void show_path(Graph *graph, uint32_t vertexNum, PATH_t *dist_path)
{
	int32_t i,j;
	
	for ( i = 0; i < vertexNum; ++i)
	{
		//print the path
		if (graph->vnode[i].adjacent_node > 0)
		{
			fprintf(stderr, "dist = %f\n", dist_path[i].dist);
			fprintf(stderr,"%d->",i );
			j = dist_path[i].pre_node;
			fprintf(stderr,"%d->",j );
			while (graph->vnode[j].adjacent_node > 0)
			{
				fprintf(stderr, "%d->", dist_path[j].pre_node);
				j = dist_path[j].pre_node;
			}
			fprintf(stderr, "\n");
		}
	}
}
*/

void dynamic_programming_path(Graph *graph, uint32_t vertexNum, PATH_t *dist_path, uint8_t *out_degree)
{
	uint32_t i;

	// for all For each v∈V\S in Topological Sorting order do dilg(v)=max(u,v)∈E{dilg(u)+w(u, v)}
	for (i = 0; i < vertexNum; ++i)
	{
		if (graph->vnode[i].adjacent_node > 0) //calculate dist[i]
		{
			Find_longest_path(graph, i, dist_path);
		}
	}

#ifdef Annoation
	fprintf(stderr, "show the path\n");
	int j;
	for (i = 0; i < vertexNum; ++i)
	{
		j = i;
		if ((out_degree[i] == 0))
		{
			fprintf(stderr, "vertex: %d\tdist = %f\n", i, dist_path[i].dist);
			fprintf(stderr, "path: ");
			fprintf(stderr, "%d->", i);

			j = dist_path[j].pre_node;
			while(j != -1)
			{
				fprintf(stderr, "%d->", j);
				j = dist_path[j].pre_node;
			}
			fprintf(stderr, "\n");
		}
	}
#endif
}

