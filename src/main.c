/*************************************************************************
	> File Name: main.c
	> Author: 
	> Mail: 
	> Created Time: 2018年05月08日 星期二 18时44分23秒
 ************************************************************************/

#include<stdio.h>
#include<string.h>
#include <time.h>

#include "desalt_index.h"
#include "read_seeding.h"
#include "ktime.h"

#define VERSION "1.5.1"
#define CONTACT "Yadong Liu <hitliuyadong1994@163.com>"

static int usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	deSALT (Third generation RNA sequence alignment)\n");
	fprintf(stderr, "Version:	%s\n", VERSION);
	fprintf(stderr, "Contact:	%s\n\n", CONTACT);

	fprintf(stderr, "Usage:		deSALT <command> [options]\n\n");
	fprintf(stderr, "Command: \n");
	fprintf(stderr, "		index	index reference sequence\n");
	fprintf(stderr, "		aln	align long RNA sequence to reference\n");

	fprintf(stderr, "\n");

	return 1;
}

int main(int argc, char *argv[])
{
	int r = 1;
	double realtime0 = realtime();
	double ts = clock();
	if (argc < 2)	return usage();
	if (strcmp(argv[1], "index") == 0)	r = desalt_index(argc, argv);
	else if (strcmp(argv[1], "aln") == 0)	r = desalt_aln(argc, argv, VERSION);
	else if (strcmp(argv[1], "--help") == 0)	return help_usage();
	else {
		fprintf(stderr, "[Waring!!!] wrong command: '%s'\n", argv[1]);
		return 1;
	}
	if (!r)
		fprintf(stderr, "[Main] Real time: %.3f sec; CPU: %.3f sec, Memory peak: %.2f GB\n", realtime() - realtime0, cputime(), peak_memory() / 1024.0 / 1024.0 / 1024.0);

	return 0;
}
