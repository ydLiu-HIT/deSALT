/*************************************************************************
	> File Name: tgs_index.c
	> Author: 
	> Mail: 
	> Created Time: 2018年05月18日 星期五 13时51分01秒
 ************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<getopt.h>
#include<unistd.h>
#include<zlib.h>

int tgs_index_usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:	deSALT index <ref.fa> <index_route>\n");
	fprintf(stderr, "		build deBGA index file. You can get more deBGA information from https://github.com/HongzheGuo/deBGA");
	fprintf(stderr, "\n");
	return 1;
}

int get_bin_dir(char *bin, char *dir)
{
	char *end = strrchr(bin, '/');
	if (end == NULL)
		return 1;

	bin[end-bin] = '\0';
	strcpy(dir, bin);
	return 0;
}

int deBGA_index(char *dir, char *ref_fa, char *index_route)
{
	char cmd[1024];
	sprintf(cmd, "%sDeBGA index %s %s", dir, ref_fa, index_route);
	fprintf(stderr, "[tgs_index] Ececuting deBGA index ...\n");
	if (system(cmd) != 0)
	{
		fprintf(stderr, "\n[tgs_index] Indexing undoing, deBGA index exit abnormally. \n"); 
		exit(1);
	}
	fprintf(stderr, "[tgs_index] Done!\n");
}

int tgs_index(int argc, char *argv[])
{
	// clock_t t = clock();
	char dir[1024];
	char *ref_fa = 0;
	char *index_route = 0;

	if(optind + 3 > argc)	return tgs_index_usage();

	if (!get_bin_dir(argv[0], dir))
	{
		strcat(dir, "/");
	}
	ref_fa = strdup(argv[optind+1]);
	index_route = strdup(argv[optind+2]);

	deBGA_index(dir, ref_fa, index_route);

	return 0;
}
