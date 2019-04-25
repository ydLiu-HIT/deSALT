/*
 * =====================================================================================
 *        Version:  0.1
 *        Created:  02/11/2014 08:20:41 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Hongzhe Guo, hzguo@hit.edu.cn
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <string.h>

#define COMMAND_VERSION "0.1"

int index_build(int argc, char *argv[]);
int load_input_map(int argc, char *argv[]);
int help_usage();

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	deBGA (De bruijn graph nucleotide alignment)\n");
	fprintf(stderr, "Version:	%s\n", COMMAND_VERSION);
	fprintf(stderr, "Contact:	Hongzhe Guo <hzguo@hit.edu.cn>\n\n");
	fprintf(stderr, "Usage:  	deBGA <command> [options]\n\n");
	fprintf(stderr, "Command:	index		index sequences in the FASTA format\n");
	fprintf(stderr, "		aln      	pair-end and single-end reads seed reduction and alignment based on De bruijn graph\n");
	
	return 1;
}

//else if (strcmp(argv[1], "alnpe") == 0) re = seed_ali(argc - 1, argv + 1);
//else if (strcmp(argv[1], "alnse") == 0) re = seed_ali_single_end(argc - 1, argv + 1);

int main(int argc, char *argv[])
{
	int re = 1;
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0) re = index_build(argc - 1, argv + 1);
	else if (strcmp(argv[1], "aln") == 0) re = load_input_map(argc - 1, argv + 1);
	else if (strcmp(argv[1], "--help") == 0) return help_usage();
	else {
		fprintf(stderr, "wrong command: '%s'\n", argv[1]);
		return 1;
	}
	
	if(re == 0)
	{
		fprintf(stderr, "Program finished\n");
	}
}




