/*************************************************************************
	> File Name: graph.h
	> Author: 
	> Mail: 
 ************************************************************************/

#ifndef _GRAPH_H
#define _GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "read_seeding.h"

// #define UNI_S
#define SHOW_UNISEED

typedef struct ArcNode
{
	uint32_t adjvex;
    int weight; 
    float penalty;
}ANode, *pANode;

typedef struct VertexNode
{
	uint32_t head_vertex;
    int weight;
	uint32_t adjacent_node;
    ANode* preedge;
}VNode;

typedef struct dpGraph
{
	VNode* vnode;
}Graph;

extern Graph* graph;

float creatGraph(uni_seed* vertexArr, uint32_t vertexNum, PATH_t* dist_path, uint8_t *out_degree, uint32_t *max_index, uint8_t tid, int max_read_join_gap);
void initGraph();
void delGraph();
void show_vertexm(vertex_m *vertexArr, uint32_t vertexNum);
void show_vertexu(vertex_u *vertexArr, uint32_t vertexNum);
#endif
