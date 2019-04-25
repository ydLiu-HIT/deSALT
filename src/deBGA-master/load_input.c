
#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <dirent.h>
#include <inttypes.h>
#include <inttypes.h>

#include "load_input.h"

//KSEQ_INIT(gzFile, gzread)
kseq_t *seq1 = NULL;
kseq_t *seq2 = NULL;

uint64_t* buffer_ref_seq = NULL;
#ifdef UNI_SEQ64
uint64_t* buffer_seq = NULL;
#else
uint8_t* buffer_seq = NULL;
#endif

#ifdef UNPIPATH_OFF_K20
uint64_t* buffer_seqf = NULL;
uint64_t* buffer_off_g = NULL;
uint64_t* buffer_p = NULL;
uint64_t* buffer_pp = NULL;
uint64_t* buffer_hash_g = NULL;
uint64_t first_cnt_cl = 0;
uint64_t chr_end_n[MAX_CHR_NUM];
#else
uint32_t* buffer_seqf = NULL;
uint32_t* buffer_off_g = NULL;
uint32_t* buffer_p = NULL;
uint32_t* buffer_pp = NULL;
uint32_t* buffer_hash_g = NULL;
uint32_t first_cnt_cl = 0;
uint32_t chr_end_n[MAX_CHR_NUM];
#endif

//uint8_t* buffer_edge = NULL;

uint32_t* buffer_pupos = NULL;
uint32_t* buffer_puid = NULL;

uint32_t* buffer_kmer_g = NULL;
	
uint64_t result_ref_seq = 0;
uint64_t result_seq = 0;
uint64_t result_seqf = 0;
uint64_t result_edge = 0;
uint64_t result_p = 0;
uint64_t result_pp = 0;
uint64_t result_pu = 0;
uint64_t result_hash_g = 0;
uint64_t result_kmer_g = 0;
uint64_t result_off_g = 0;
uint64_t result_ref_g = 0;

const char* refseqs = "ref.seq";
const char* uniseqs = "unipath.seq";
const char* uniseq_bs = "unipath.seqb";
const char* uniseqf_bs = "unipath.seqfb";
const char* uniedges = "unipath.edge";
const char* unipus = "unipath.pu";
const char* uniposs = "unipath.pos";
const char* uniposps = "unipath.posp";
//const char* unistas = "unipath.sta";
const char* unihash_gs = "unipath_g.hash";
const char* unikmer_gs = "unipath_g.kmer";
const char* unioff_gs = "unipath_g.offset";
const char* unichrs = "unipath.chr";
const char* f_n = "N.sta";
const char* graph = "graph.sta";
const char* divs = "div/";
const char* f_size = "unipath.size";
const char* sys_c_mkdir = "mkdir ";
const char* sys_c_rm = "rm -rf ";
const char* sys_rm = "rm ";

char sam_result[ROUTE_LENGTH_MAX];
char read_fastq1[ROUTE_LENGTH_MAX];
char read_fastq2[ROUTE_LENGTH_MAX];
char index_route[ROUTE_LENGTH_MAX];
char filename_ref[ROUTE_LENGTH_MAX];
char filename_div[ROUTE_LENGTH_MAX];
char filename_sta[ROUTE_LENGTH_MAX];
char ref_seq[ROUTE_LENGTH_MAX];
char uniseq[ROUTE_LENGTH_MAX];
char uniseq_b[ROUTE_LENGTH_MAX];
char uniseqf_b[ROUTE_LENGTH_MAX];
char uniedge[ROUTE_LENGTH_MAX];
char unipu[ROUTE_LENGTH_MAX];
char unipos[ROUTE_LENGTH_MAX];
char uniposp[ROUTE_LENGTH_MAX];
//char unista[ROUTE_LENGTH_MAX];
char unihash_g[ROUTE_LENGTH_MAX];
char unikmer_g[ROUTE_LENGTH_MAX];
char unioff_g[ROUTE_LENGTH_MAX];
char unichr[ROUTE_LENGTH_MAX];
char N_route[ROUTE_LENGTH_MAX];
char unisize[ROUTE_LENGTH_MAX];
char chr_names[MAX_CHR_NUM][MAX_CHR_NAME_LENGTH];
char chr_line_content[MAX_CHR_NAME_LENGTH];

//files
FILE* fp_ref_seq = NULL;
FILE* fp_us_b = NULL;
FILE* fp_usf_b = NULL;
FILE* fp_ub = NULL;
//FILE* fp_ue = NULL;
FILE* fp_up = NULL;
FILE* fp_upp = NULL;
FILE* fp_pu = NULL;
FILE* fp_sta = NULL;
FILE* fp_chr = NULL;
FILE* fp_uh = NULL;
FILE* fp_uf = NULL;
FILE* fp_hash = NULL;
FILE* fp_kmer = NULL;
FILE* fp_off = NULL;
FILE* fp_n = NULL;
FILE* fp_us = NULL;
FILE* fp_num = NULL;
FILE* unipath_debug = NULL;

const uint8_t f = 4;
const uint8_t k = 14;
uint8_t k_t = 22;	
uint32_t chr_file_n = 1;
uint8_t thread_n = 1;
uint32_t upper_ins = 700;
uint32_t floor_ins = 300;
uint16_t readlen_max = 512;//252
uint16_t seed_l_max = 0;
uint8_t seed_l_l = 5;
uint16_t seed_step = 5; 
uint8_t cir_fix_n = 4;
uint16_t pos_n_max = 300;
uint16_t length_reduce = 22;
float va_ra = 0.02;

//important parameter 
float mis_match_r = 0.04;
float lv_rate = 0.06;//0.1
float last_circle_rate = 0;

//
float max_pair_score_r = 0.05;//0.06 0.2 
float mis_match_r_single = 0.05;//0.04 0.4
float lv_rate_anchor = 0.05;

uint8_t flag_std = 0;
uint8_t local_ksw = 0;
uint8_t mgn_flag = 1;

uint16_t cus_ali_n = 20; //150
uint16_t cus_max_output_ali = 150;//300 anchor
float max_single_score_r = 0.06;
uint16_t pr_single_outputn = 1000;
uint16_t seed_filter_pos_numn = 100;
uint16_t seed_filter_pos_num_singlen = 100;//

int8_t mat_score = 1, mis_score = 4, gapo_score = 6, gape_score = 1;

static int index_build_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: de Brijn Graph-based mapping system index building\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Hongzhe Guo <hzguo@hit.edu>\n\n");
	fprintf(stderr, "Usage:   deBGA index [options] reference.fasta <index_route> \n\n");
	fprintf(stderr, "Options: -k INT      the k-mer length of the vertices of RdBG [20-28]\n");
	
	fprintf(stderr, "\n");
	return 1;
}

static int load_input_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	de Brijn Graph-based mapping system seed reduction and alignment\n");
	fprintf(stderr, "Version:	%s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact:	Hongzhe Guo <hzguo@hit.edu>\n\n");
	fprintf(stderr, "Usage:  	deBGA aln [options] <index_route> <read pair-end1.fq> [read pair-end2.fq] <result_file.sam>\n\nOptions:\n");
	fprintf(stderr, "	-k INT	the minimum length of a valid Uni-MEM seed [21-28]\n");	
	fprintf(stderr, "	-s INT	the number of iterations of re-seeding [%u]\n", cir_fix_n);
	fprintf(stderr, "	-i INT	the minimum interval of seeding [%u]\n", seed_step);
	fprintf(stderr, "	-n INT	the maximum allowed number of hits per seed [%u]\n", pos_n_max);
	fprintf(stderr, "	-c NUM	the threshold on the edit distance for early stop [%.2f]\n", max_pair_score_r);
	fprintf(stderr, "	--cl NUM the adjusted threshold on the edit distance for early stop [%.2f]\n", last_circle_rate);
	fprintf(stderr, "	--local  the local alignment option for confident alignment\n");
	fprintf(stderr, "	--local-match NUM the score for a matched base in the local alignment [%d]\n", mat_score);
	fprintf(stderr, "	--local-mismatch NUM the penalty for a mismatched base in the local alignment [%d]\n", mis_score);
	fprintf(stderr, "	--local-gap-open NUM the penalty for a gap open in the local alignment [%d]\n", gapo_score);
	fprintf(stderr, "	--local-gap-extension NUM the penalty for gap extension in the local alignment [%d]\n", gape_score);
	fprintf(stderr, "	--stdout   (default: not set) output alignments by stdout\n");
	fprintf(stderr, "	-u INT	the upper limit of insert size (only for pair-end reads) [%u] \n", upper_ins);
	fprintf(stderr, "	-f INT	the lower limit of insert size (only for pair-end reads) [%u] \n", floor_ins);
	fprintf(stderr, "	-o INT	the maximum number of alignment output [%u]\n", cus_ali_n);
	fprintf(stderr, "	-x INT	the maximum number of alignment output for anchoring alignment [%u]\n", cus_max_output_ali);
	fprintf(stderr, "	-l INT	the maximum allowed read length [%u]\n", readlen_max);
	fprintf(stderr, "	-e INT	the budget for single-end alignment [%u]\n", seed_filter_pos_numn);
	fprintf(stderr, "	-p INT	the number of threads [%u]\n", thread_n);
	//fprintf(stderr, "	--mg 	use the mode of multi-genomes\n");

	fprintf(stderr, "	Please refer to the following link for more detailed information about the options: https://github.com/HIT-Bioinformatics/deBGA\n");
	fprintf(stderr, "\n");

	return 1;
}

int help_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	deBGA (De bruijn graph nucleotide alignment)\n");
	fprintf(stderr, "Usage:  	deBGA <command> [options]\n\n");
	fprintf(stderr, "Command:	index		index sequences in the FASTA format\n");
	fprintf(stderr, "		aln      	pair-end and single-end reads seed reduction and alignment based on De bruijn graph\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   deBGA index [options] reference.fasta <index_route> \n\n");
	fprintf(stderr, "Options: -k INT      the k-mer length of the vertices of RdBG [20-28]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:  	deBGA aln [options] <index_route> <read pair-end1.fq> [read pair-end2.fq] <result_file.sam>\n\nOptions:\n");
	fprintf(stderr, "	-k INT	the minimum length of a valid Uni-MEM seed [21-28]\n");	
	fprintf(stderr, "	-s INT	the number of iterations of re-seeding [%u]\n", cir_fix_n);
	fprintf(stderr, "	-i INT	the minimum interval of seeding [%u]\n", seed_step);
	fprintf(stderr, "	-n INT	the maximum allowed number of hits per seed [%u]\n", pos_n_max);
	fprintf(stderr, "	-c NUM	the threshold on the edit distance for early stop [%.2f]\n", max_pair_score_r);
	fprintf(stderr, "	--cl NUM the adjusted threshold on the edit distance for early stop [%.2f]\n", last_circle_rate);
	fprintf(stderr, "	--local  the local alignment option for confident alignment\n");
	fprintf(stderr, "	--local-match NUM the score for a matched base in the local alignment [%d]\n", mat_score);
	fprintf(stderr, "	--local-mismatch NUM the penalty for a mismatched base in the local alignment [%d]\n", mis_score);
	fprintf(stderr, "	--local-gap-open NUM the penalty for a gap open in the local alignment [%d]\n", gapo_score);
	fprintf(stderr, "	--local-gap-extension NUM the penalty for gap extension in the local alignment [%d]\n", gape_score);
	fprintf(stderr, "	--stdout   (default: not set) output alignments by stdout\n");
	fprintf(stderr, "	-u INT	the upper limit of insert size (only for pair-end reads) [%u] \n", upper_ins);
	fprintf(stderr, "	-f INT	the lower limit of insert size (only for pair-end reads) [%u] \n", floor_ins);
	fprintf(stderr, "	-o INT	the maximum number of alignment output [%u]\n", cus_ali_n);
	fprintf(stderr, "	-x INT	the maximum number of alignment output for anchoring alignment [%u]\n", cus_max_output_ali);
	fprintf(stderr, "	-l INT	the maximum allowed read length [%u]\n", readlen_max);
	fprintf(stderr, "	-e INT	the budget for single-end alignment [%u]\n", seed_filter_pos_numn);
	fprintf(stderr, "	-p INT	the number of threads [%u]\n", thread_n);
	//fprintf(stderr, "	--mg 	use the mode of multi-genomes\n");

	fprintf(stderr, "	Please refer to the following link for more detailed information about the options: https://github.com/HIT-Bioinformatics/deBGA\n");
	fprintf(stderr, "\n");

	
	return 1;
}

enum {
	PAR_LOCAL_KSW,
	PAR_LAST_CIRCLE,
	PAR_MULTI_GENOMES,
	PAR_MAT_SCORE,
	PAR_MIS_SCORE,
	PAR_GAPO_SCORE,
	PAR_GAPE_SCORE,
	PAR_A_ANCHOR,
	PAR_LV_ANCHOR,
	PAR_HELP,
	PAR_STD_OUT
};

static const char *short_option = "k:p:u:f:l:r:i:s:n:o:x:a:c:g:e:v:";
static struct option long_option[] = {
	{(char*)"local", no_argument,  0, PAR_LOCAL_KSW},
	{(char*)"cl", required_argument, 0, PAR_LAST_CIRCLE},
	{(char*)"mg", no_argument, 0, PAR_MULTI_GENOMES},
	{(char*)"local-match", required_argument, 0, PAR_MAT_SCORE},
	{(char*)"local-mismatch", required_argument, 0, PAR_MIS_SCORE},
	{(char*)"local-gap-open", required_argument, 0, PAR_GAPO_SCORE},
	{(char*)"local-gap-extension", required_argument, 0, PAR_GAPE_SCORE},
	{(char*)"aanchor", required_argument, 0, PAR_A_ANCHOR},
	{(char*)"vanchor", required_argument, 0, PAR_LV_ANCHOR},
	{(char*)"help", no_argument, 0, PAR_HELP},
	{(char*)"stdout", no_argument, 0, PAR_STD_OUT},
	{(char*)0, 0, 0, 0}
};

int load_input_index(int argc, char *argv[])
{
	int c = 0;

	//while ((c = getopt(argc, argv, "k:")) >= 0) {
	while((c = getopt_long(argc, argv, short_option, long_option, NULL)) != -1){
		switch (c) {
		case 'k': k_t = atoi(optarg); break;
		case PAR_HELP: return index_build_usage();
		//default: return index_build_usage();
		}
	}
	
	if (argc - optind < 2) return index_build_usage();
	
	memset(filename_ref, 0, ROUTE_LENGTH_MAX);
	memset(index_route, 0, ROUTE_LENGTH_MAX);
	
	strcpy(filename_ref, argv[optind]);
	strcpy(index_route, argv[optind + 1]);

	if(index_route[strlen(index_route) - 1] != '/')	strcat(index_route, "/");
	
	memset(ref_seq, 0, ROUTE_LENGTH_MAX);
	strcpy(ref_seq, index_route);
	strcat(ref_seq, refseqs);
	
	memset(uniseq, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseq, index_route);
	strcat(uniseq, uniseqs);
	
	memset(uniseq_b, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseq_b, index_route);
	strcat(uniseq_b, uniseq_bs);
	
	memset(uniseqf_b, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseqf_b, index_route);
	strcat(uniseqf_b, uniseqf_bs);
	
	memset(uniedge, 0, ROUTE_LENGTH_MAX);
	strcpy(uniedge, index_route);
	strcat(uniedge, uniedges);
	
	memset(unipu, 0, ROUTE_LENGTH_MAX);
	strcpy(unipu, index_route);
	strcat(unipu, unipus);
	
	memset(unipos, 0, ROUTE_LENGTH_MAX);
	strcpy(unipos, index_route);
	strcat(unipos, uniposs);
	
	memset(uniposp, 0, ROUTE_LENGTH_MAX);
	strcpy(uniposp, index_route);
	strcat(uniposp, uniposps);
	/*
	memset(unista, 0, ROUTE_LENGTH_MAX);
	strcpy(unista, index_route);
	strcat(unista, unistas);
	*/
	memset(unihash_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unihash_g, index_route);
	strcat(unihash_g, unihash_gs);
	
	memset(unikmer_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unikmer_g, index_route);
	strcat(unikmer_g, unikmer_gs);
	
	memset(unioff_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unioff_g, index_route);
	strcat(unioff_g, unioff_gs);
	
	memset(unichr, 0, ROUTE_LENGTH_MAX);
	strcpy(unichr, index_route);
	strcat(unichr, unichrs);
	
	memset(N_route, 0, ROUTE_LENGTH_MAX);
	strcpy(N_route, index_route);
	strcat(N_route, f_n);
	
	memset(filename_div, 0, ROUTE_LENGTH_MAX);
	strcpy(filename_div, index_route);
	strcat(filename_div, divs);
	
	memset(filename_sta, 0, ROUTE_LENGTH_MAX);
	strcpy(filename_sta, filename_div);
	strcat(filename_sta, graph);
	
	memset(unisize, 0, ROUTE_LENGTH_MAX);
	strcpy(unisize, index_route);
	strcat(unisize, f_size);
	
#ifdef	HANDLE_DIR
	char create_route[ROUTE_LENGTH_MAX];
	char rm_route[ROUTE_LENGTH_MAX];
	DIR* directory_pointer = NULL;

	if((directory_pointer = opendir(index_route)) != NULL)//filename_div
	{
		memset(rm_route, 0, ROUTE_LENGTH_MAX);
		strcpy(rm_route, sys_c_rm);
		strcat(rm_route, index_route);

		system(rm_route);
	}
		
	memset(create_route, 0, ROUTE_LENGTH_MAX);
	strcpy(create_route, sys_c_mkdir);
	strcat(create_route, index_route);

	if((directory_pointer = opendir(index_route))==NULL)
		system(create_route);
	
	memset(create_route, 0, ROUTE_LENGTH_MAX);
	strcpy(create_route, sys_c_mkdir);
	strcat(create_route, filename_div);
	
	if((directory_pointer = opendir(filename_div))==NULL)
		system(create_route);
#endif
		
	return 0;
}

int load_input_map(int argc, char *argv[])
{
	int c = 0;
	while((c = getopt_long(argc, argv, short_option, long_option, NULL)) != -1){
		switch (c) {
			case 'k': k_t = atoi(optarg); break;
			case 'p': thread_n = atoi(optarg); break;
			case 'u': upper_ins = atoi(optarg); break;
			case 'f': floor_ins = atoi(optarg); break;
			case 'l': readlen_max = atoi(optarg); break;
			case 'r': length_reduce = atoi(optarg); break;
			case 'i': seed_step = atoi(optarg); break;
			case 's': cir_fix_n = atoi(optarg); break;
			case 'n': pos_n_max = atoi(optarg); break;
			case 'o': cus_ali_n = atoi(optarg); break;
			case 'x': cus_max_output_ali = atoi(optarg); break;
			case 'c': max_pair_score_r = atof(optarg); lv_rate_anchor = max_pair_score_r; mis_match_r_single = max_pair_score_r; break;
			case 'e': seed_filter_pos_numn = atoi(optarg); break;
			case 'v': lv_rate = atof(optarg); break;
			case PAR_A_ANCHOR: mis_match_r_single = atof(optarg); break;
			case PAR_LV_ANCHOR: lv_rate_anchor = atof(optarg); break;
			case PAR_LAST_CIRCLE: last_circle_rate = atof(optarg); mis_match_r_single = last_circle_rate; break;
			case PAR_LOCAL_KSW: local_ksw = 1; last_circle_rate = 0; mis_match_r_single = 0.4; break;
			case PAR_MULTI_GENOMES: mgn_flag = 0;break;
			case PAR_MAT_SCORE: mat_score = atoi(optarg); break;
			case PAR_MIS_SCORE: mis_score = atoi(optarg); break;
			case PAR_GAPO_SCORE: gapo_score = atoi(optarg); break;
			case PAR_GAPE_SCORE: gape_score = atoi(optarg); break;
			case PAR_STD_OUT: flag_std = 1; break;
			case PAR_HELP: return load_input_usage();
			default: return load_input_usage();
		} 
	}   
	
	if (argc - optind < 3) return load_input_usage();
	if((k_t < 21) || (k_t > 28))
	{
		fprintf(stderr, "Input error: -k cannot be less than 21 or more than 28\n");
		exit(1);
	}
	if((thread_n < 1) || (thread_n > 32))
	{
		fprintf(stderr, "Input error: -p cannot be less than 1 or more than 32\n");
		exit(1);
	}
	if(upper_ins <= floor_ins)
	{
		fprintf(stderr, "Input error: -u should be more than floor limit\n");
		exit(1);
	}
	if(readlen_max > 2048)
	{
		fprintf(stderr, "Input error: -l cannot be more than 2048\n");
		exit(1);
	}
	if(length_reduce > 50)
	{
		fprintf(stderr, "Input error: -r cannot be more than 50\n");
		exit(1);
	}
	if((seed_step < 1) || (seed_step > 20))//5 20
	{
		fprintf(stderr, "Input error: -i cannot be more than 20 or less than 5\n");
		exit(1);
	}
	if(cir_fix_n > 10)
	{
		fprintf(stderr, "Input error: -s cannot be more than 10\n");
		exit(1);
	}
	if(cus_ali_n > 1000)//300
	{
		fprintf(stderr, "Input error: -o cannot be more than 1000\n");//300
		exit(1);
	}
	if(cus_max_output_ali > 1000)//500
	{
		fprintf(stderr, "Input error: -x cannot be more than 1000\n");//500
		exit(1);
	}
	if(max_pair_score_r > 0.35)//0.5
	{
		fprintf(stderr, "Input error: -c cannot be more than 0.5\n");
		exit(1);
	}
	if(seed_filter_pos_numn > 2000)
	{
		fprintf(stderr, "Input error: -e cannot be more than 2000\n");
		exit(1);
	}
	if(lv_rate > 0.5)
	{
		fprintf(stderr, "Input error: -v cannot be more than 0.5\n");
		exit(1);
	}
	if(last_circle_rate > 0.5)
	{
		fprintf(stderr, "Input error: --cl cannot be more than 0.5\n");
		exit(1);
	}
	
	memset(index_route, 0, ROUTE_LENGTH_MAX);
	strcpy(index_route, argv[optind]);
	
	if(index_route[strlen(index_route) - 1] != '/')	strcat(index_route, "/");
	
	memset(read_fastq1, 0, ROUTE_LENGTH_MAX);
	strcpy(read_fastq1, argv[optind + 1]);

	memset(ref_seq, 0, ROUTE_LENGTH_MAX);
	strcpy(ref_seq, index_route);
	strcat(ref_seq, refseqs);
	
	memset(uniseq_b, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseq_b, index_route);
	strcat(uniseq_b, uniseq_bs);
	
	memset(uniseqf_b, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseqf_b, index_route);
	strcat(uniseqf_b, uniseqf_bs);
	
	memset(uniedge, 0, ROUTE_LENGTH_MAX);
	strcpy(uniedge, index_route);
	strcat(uniedge, uniedges);
	
	memset(unipu, 0, ROUTE_LENGTH_MAX);
	strcpy(unipu, index_route);
	strcat(unipu, unipus);
	
	memset(unipos, 0, ROUTE_LENGTH_MAX);
	strcpy(unipos, index_route);
	strcat(unipos, uniposs);
	
	memset(uniposp, 0, ROUTE_LENGTH_MAX);
	strcpy(uniposp, index_route);
	strcat(uniposp, uniposps);
	/*
	memset(unista, 0, ROUTE_LENGTH_MAX);
	strcpy(unista, index_route);
	strcat(unista, unistas);
	*/
	memset(unihash_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unihash_g, index_route);
	strcat(unihash_g, unihash_gs);
	
	memset(unikmer_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unikmer_g, index_route);
	strcat(unikmer_g, unikmer_gs);
	
	memset(unioff_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unioff_g, index_route);
	strcat(unioff_g, unioff_gs);
	
	memset(unichr, 0, ROUTE_LENGTH_MAX);
	strcpy(unichr, index_route);
	strcat(unichr, unichrs);
	
	memset(N_route, 0, ROUTE_LENGTH_MAX);
	strcpy(N_route, index_route);
	strcat(N_route, f_n);
	
	memset(filename_div, 0, ROUTE_LENGTH_MAX);
	strcpy(filename_div, index_route);
	strcat(filename_div, divs);
	
	memset(filename_sta, 0, ROUTE_LENGTH_MAX);
	strcpy(filename_sta, filename_div);
	strcat(filename_sta, graph);
	
	memset(unisize, 0, ROUTE_LENGTH_MAX);
	strcpy(unisize, index_route);
	strcat(unisize, f_size);
	
	if ((argc - optind < 4) && (flag_std == 0))
	{
		max_single_score_r = max_pair_score_r;
		seed_filter_pos_num_singlen = seed_filter_pos_numn;

		if(strcmp(argv[optind + 2] + strlen(argv[optind + 2]) - 3, ".fq") || strcmp(argv[optind + 2] + strlen(argv[optind + 2]) - 6, ".fastq"))
		{
			fprintf(stderr, "result file cannot be .fq or .fastq\n");
			exit(1);
		}
		strcpy(sam_result, argv[optind + 2]);
		seed_ali_single_end();
	}else{
		memset(read_fastq2, 0, ROUTE_LENGTH_MAX);
		strcpy(read_fastq2, argv[optind + 2]);
		if(flag_std == 0)
			strcpy(sam_result, argv[optind + 3]);
		seed_ali();
	}	

	return 0;
}

void load_index_file()
{
	uint64_t a_size = 0;
	uint64_t us_n = 0;
	uint64_t usf_n = 0;
	uint64_t ue_n = 0;
	uint64_t up_n = 0;
	uint64_t upp_n = 0;
	uint64_t hash_n = 0;
	uint64_t kmer_n = 0;
	uint64_t off_n = 0;
	uint64_t pu_n = 0;
	uint64_t ref_seq_n = 0;
	uint64_t file_size = 0;
	
	//read index file size
	fp_num = fopen(unisize, "rb");
	if (fp_num == NULL)
    {
        fputs ("File error opening the size of index file: wrong index route or wrong index file name\n",stderr);
        exit (1);
    }
	fread(&us_n, 8, 1, fp_num);
	fread(&usf_n, 8, 1, fp_num);
	fread(&ue_n, 8, 1, fp_num);
	fread(&up_n, 8, 1, fp_num);
	fread(&upp_n, 8, 1, fp_num);
	fread(&hash_n, 8, 1, fp_num);
	fread(&kmer_n, 8, 1, fp_num);
	fread(&off_n, 8, 1, fp_num);
	fread(&pu_n, 8, 1, fp_num);
	fread(&ref_seq_n, 8, 1, fp_num);

	fclose(fp_num);
	
	//read ref seq file
	fprintf(stderr, "Load ref seq\n");
	
	fp_ref_seq = fopen(ref_seq, "rb");
	if (fp_ref_seq == NULL)
    {
        fputs ("File error opening the  seq file\n",stderr);
        exit (1);
    }
	
	fseek(fp_ref_seq, 0, SEEK_END);// non-portable
    file_size = ftell(fp_ref_seq);
	rewind(fp_ref_seq);

	ref_seq_n = file_size;
	buffer_ref_seq = (uint64_t* )calloc(ref_seq_n + 536, 1);//536 = (2048 >> 5 + 3) << 3
	
	a_size = ref_seq_n >> 3;
	result_ref_seq = fread(buffer_ref_seq, 8, a_size, fp_ref_seq);
	
	if (result_ref_seq != a_size)
    {
        fprintf(stderr, "Reading error");
        exit (3);
    }

    fclose(fp_ref_seq);

    //read input unipath seq file
    fprintf(stderr, "Load unipath seq\n");

    fp_us_b = fopen (uniseq_b, "rb" );
    if (fp_us_b == NULL)
    {
        fprintf(stderr, "File error opening the unipath seq file\n");
        exit (1);
    }

	fseek(fp_us_b, 0, SEEK_END);// non-portable
    file_size = ftell(fp_us_b);
	rewind(fp_us_b);

	us_n = file_size;
	
    // allocate memory to contain the whole file:
#ifdef UNI_SEQ64
	buffer_seq = (uint64_t* ) malloc (us_n);
#else	
    buffer_seq = (uint8_t* ) malloc (us_n);
#endif
    if (buffer_seq == NULL)
    {
        fputs ("Memory error buffer_seq\n",stderr);
        exit (2);
    }

    // copy the file into the buffer:
	
#ifdef UNI_SEQ64	
	a_size = us_n >> 3;
	result_seq = fread (buffer_seq, 8, a_size, fp_us_b);
#else	
	a_size = us_n;
    result_seq = fread (buffer_seq, 1, a_size, fp_us_b);
#endif

    if (result_seq != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }

    fclose(fp_us_b);

    //read input unipath offset file
    fprintf(stderr, "Load unipath offset\n");

    fp_usf_b = fopen(uniseqf_b, "rb");
    if(fp_usf_b == NULL)
    {
        fputs ("File error opening the unipath offset file\n",stderr);
        exit (1);
    }

	fseek(fp_usf_b, 0, SEEK_END);// non-portable
    file_size = ftell(fp_usf_b);
	rewind(fp_usf_b);

	usf_n = file_size;
	
    // allocate memory to contain the whole file:
	// copy the file into the buffer:
#ifdef UNPIPATH_OFF_K20
	buffer_seqf = (uint64_t* ) malloc (usf_n);
	if (buffer_seqf == NULL)
    {
        fprintf(stderr, "Memory error buffer_seqf");
        exit (2);
    }
	a_size = (usf_n >> 3);
    result_seqf = fread (buffer_seqf, 8, a_size, fp_usf_b);
#else	
    buffer_seqf = (uint32_t* ) malloc (usf_n);
	if (buffer_seqf == NULL)
    {
        fputs ("Memory error buffer_seqf",stderr);
        exit (2);
    }
	a_size = (usf_n >> 2);
    result_seqf = fread (buffer_seqf, 4, a_size, fp_usf_b);
#endif	
	
    if (result_seqf != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }

    fclose(fp_usf_b);

	/*
    //read input unipath edge file
    printf("Load unipath edge\n");

    fp_ue = fopen(uniedge, "rb");
    if(fp_ue == NULL)
    {
        fputs ("File error opening the unipath offset file\n",stderr);
        exit (1);
    }

    // allocate memory to contain the whole file:
    buffer_edge = (uint8_t* ) malloc (ue_n);
    if (buffer_edge == NULL)
    {
        fputs ("Memory error",stderr);
        exit (2);
    }

    // copy the file into the buffer:
    result_edge = fread (buffer_edge, 1, ue_n, fp_ue);
    if (result_edge != ue_n)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }
	
    fclose(fp_ue);
	*/
		
    //read input unipath position file
    fprintf(stderr, "Load unipath position\n");

    fp_up = fopen(unipos, "rb");
    if(fp_up == NULL)
    {
        fprintf(stderr, "File error opening the unipath position file\n");
        exit (1);
    }

	fseek(fp_up, 0, SEEK_END);// non-portable
    file_size = ftell(fp_up);
	rewind(fp_up);
	
	up_n = file_size;
	
    // allocate memory to contain the whole file:

	// copy the file into the buffer:
#ifdef	UNPIPATH_OFF_K20
	buffer_p = (uint64_t* ) malloc (up_n);
    if (buffer_p == NULL)
    {
        fputs ("Memory error buffer_p",stderr);
        exit (2);
    }
	a_size = (up_n >> 3);
	result_p = fread (buffer_p, 8, a_size, fp_up);
#else
	buffer_p = (uint32_t* ) malloc (up_n);
    if (buffer_p == NULL)
    {
        fputs ("Memory error buffer_p",stderr);
        exit (2);
    }
	a_size = (up_n >> 2);
	result_p = fread (buffer_p, 4, a_size, fp_up);
#endif
    
    if (result_p != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }
	
    fclose(fp_up);

    //read input unipath position point file
    fprintf(stderr, "Load unipath position point\n");

    fp_upp = fopen(uniposp, "rb");
    if(fp_upp == NULL)
    {
        fputs ("File error opening the unipath position point file\n",stderr);
        exit (1);
    }
	
	fseek(fp_upp, 0, SEEK_END);// non-portable
    file_size = ftell(fp_upp);
	rewind(fp_upp);

	upp_n = file_size;
	
#ifdef	UNPIPATH_OFF_K20
	// allocate memory to contain the whole file:
    a_size = (upp_n >> 3);
    buffer_pp = (uint64_t* ) malloc (upp_n);
    if (buffer_pp == NULL)
    {
        fputs ("Memory error buffer_pp",stderr);
        exit (2);
    }
    // copy the file into the buffer:
    result_pp = fread (buffer_pp, 8, a_size, fp_upp);
#else
    // allocate memory to contain the whole file:
    a_size = (upp_n >> 2);
    buffer_pp = (uint32_t* ) malloc (upp_n);
    if (buffer_pp == NULL)
    {
        fputs ("Memory error buffer_pp",stderr);
        exit (2);
    }
    // copy the file into the buffer:
    result_pp = fread (buffer_pp, 4, a_size, fp_upp);
#endif
	if (result_pp != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }
    fclose(fp_upp);
	
    //read input unipath hash file
    fprintf(stderr, "Load unipath hash\n");

    fp_hash = fopen(unihash_g, "rb");
    if(fp_hash == NULL)
    {
        fputs ("File error opening the graph hash file\n",stderr);
        exit (1);
    }
	
	fseek(fp_hash, 0, SEEK_END);// non-portable
    file_size = ftell(fp_hash);
	rewind(fp_hash);

	hash_n = file_size;
	
#ifdef	UNPIPATH_OFF_K20
	a_size = (hash_n >> 3);
    buffer_hash_g = (uint64_t* ) malloc (hash_n);
    if (buffer_hash_g == NULL)
    {
        fputs("Memory error",stderr);
        exit(2);
    }

    // copy the file into the buffer:
    result_hash_g = fread (buffer_hash_g, 8, a_size, fp_hash);
    if (result_hash_g != a_size)
    {
        fputs("Reading error",stderr);
        exit(3);
    }
#else
    a_size = (hash_n >> 2);
    buffer_hash_g = (uint32_t* ) malloc (hash_n);
    if (buffer_hash_g == NULL)
    {
        fputs("Memory error",stderr);
        exit(2);
    }

    // copy the file into the buffer:
    result_hash_g = fread (buffer_hash_g, 4, a_size, fp_hash);
    if (result_hash_g != a_size)
    {
        fputs("Reading error",stderr);
        exit(3);
    }
#endif

    fclose(fp_hash);

    //read input graph kmer file
    fprintf(stderr, "Load unipath kmer\n");

    fp_kmer = fopen(unikmer_g, "rb");
    if(fp_kmer == NULL)
    {
        fputs ("File error opening the graph hash file\n",stderr);
        exit (1);
    }
	
	fseek(fp_kmer, 0, SEEK_END);// non-portable
    file_size = ftell(fp_kmer);
	rewind(fp_kmer);

	kmer_n = file_size;
	
    a_size = (kmer_n >> 2);
    
    //a_size = kmer_num;

    buffer_kmer_g = (uint32_t* ) malloc (kmer_n);
    if (buffer_kmer_g == NULL)
    {
        fputs ("Memory error buffer_kmer_g",stderr);
        exit (2);
    }

    // copy the file into the buffer:
    result_kmer_g = fread (buffer_kmer_g, 4, a_size, fp_kmer);
    if (result_kmer_g != a_size)
    {
        fputs ("Reading error buffer_kmer_g",stderr);
        exit (3);
    }

    fclose(fp_kmer);

    //read input graph off file
    fprintf(stderr, "Load unipath off\n");

    fp_off = fopen(unioff_g, "rb");
    if(fp_off == NULL)
    {
        fprintf(stderr, "File error opening the graph hash file\n");
        exit (1);
    }
	
	fseek(fp_off, 0, SEEK_END);// non-portable
    file_size = ftell(fp_off);
	rewind(fp_off);

	off_n = file_size;
	
    // copy the file into the buffer:
#ifdef UNPIPATH_OFF_K20
	buffer_off_g = (uint64_t* ) malloc (off_n);
	if (buffer_off_g == NULL)
    {
        fputs ("Memory error",stderr);
        exit (2);
    }
	a_size = (off_n >> 3);
    result_off_g = fread (buffer_off_g, 8, a_size, fp_off);
#else	
    buffer_off_g = (uint32_t* ) malloc (off_n);
	if (buffer_off_g == NULL)
    {
        fputs ("Memory error",stderr);
        exit (2);
    }
	a_size = (off_n >> 2);
    result_off_g = fread (buffer_off_g, 4, a_size, fp_off);
#endif	

    if (result_off_g != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }

    fclose(fp_off);

	fp_chr = fopen (unichr, "r" );
    if (fp_chr == NULL)
    {
        fputs ("File error opening the chr file\n",stderr);
        exit (1);
    }
	
	uint32_t chr_line_n = 0;
    while(!feof(fp_chr))
    {
		fscanf(fp_chr, "%s", chr_line_content);
		//fgets(chr_line_content, MAX_CHR_NAME_LENGTH, fp_chr);

		if((chr_line_n & 0X1) == 0)
		{
			strcpy(chr_names[chr_file_n], chr_line_content);
			chr_names[chr_file_n][strlen(chr_names[chr_file_n])] = '\0';

		}else{
#ifdef UNPIPATH_OFF_K20
			sscanf(chr_line_content, "%"PRId64"", &chr_end_n[chr_file_n]);
#else
			sscanf(chr_line_content, "%u", &chr_end_n[chr_file_n]);
#endif
			chr_file_n++;
		}	
		
		fflush(stdout);
		
		chr_line_n++;
    }

	chr_end_n[0] = START_POS_REF + 1;

	strcpy(chr_names[chr_file_n], "*");
}


