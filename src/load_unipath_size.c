/*************************************************************************
	> File Name: load_unipath_size.c
	> Author: 
	> Mail: 
	> Created Time: 2017年03月27日 星期一 11时27分31秒
 ************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <dirent.h>
#include <inttypes.h>
#include "load_unipath_size.h"

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
uint64_t reference_len = 0;
uint64_t chr_end_n[MAX_CHR_NUM];
char chr_names[MAX_CHR_NUM][MAX_CHR_NAME_LENGTH];
char chr_line_content[MAX_CHR_NAME_LENGTH];
uint32_t chr_file_n = 1;
#else
uint32_t* buffer_seqf = NULL;
uint32_t* buffer_off_g = NULL;
uint32_t* buffer_p = NULL;
uint32_t* buffer_pp = NULL;
uint32_t* buffer_hash_g = NULL;
uint32_t reference_len = 0;
uint32_t chr_end_n[MAX_CHR_NUM];
#endif

//uint8_t* buffer_edge = NULL;

uint32_t* buffer_kmer_g = NULL;
    
uint64_t result_ref_seq = 0;
uint64_t result_seq = 0;
uint64_t result_seqf = 0;
uint64_t result_p = 0;
uint64_t result_pp = 0;
uint64_t result_pu = 0;
uint64_t result_hash_g = 0;
uint64_t result_kmer_g = 0;
uint64_t result_off_g = 0;
uint64_t result_ref_g = 0;

int load_index_file(char *index_dir){
#ifdef SIZE64
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
#else
    uint32_t a_size = 0;
    uint32_t us_n = 0;
    uint32_t usf_n = 0;
    uint32_t ue_n = 0;
    uint32_t up_n = 0;
    uint32_t upp_n = 0;
    uint32_t hash_n = 0;
    uint32_t kmer_n = 0;
    uint32_t off_n = 0;
    uint32_t pu_n = 0;
    uint32_t ref_seq_n = 0; //reference chars number
    uint32_t file_size = 0;
#endif
    //read dm file size
    char unisize[ROUTE_LENGTH_MAX] = {0};
    strcpy(unisize, index_dir);
    strcat(unisize, "unipath.size");
    // char *unisize = "/home/ydliu/deBGA/dm/unipath.size";
    FILE *fp_num = fopen(unisize,"rb");
    if(fp_num ==NULL){
        printf("File error opening file:wrong route or wrong file name\n");
        exit(1);
    }
#ifdef SIZE64
    int temp;
    temp = fread(&us_n,8,1,fp_num);
    temp = fread(&usf_n,8,1,fp_num);
    temp = fread(&ue_n, 8, 1, fp_num);
    temp = fread(&up_n, 8, 1, fp_num);
    temp = fread(&upp_n, 8, 1, fp_num);
    temp = fread(&hash_n, 8, 1, fp_num);
    temp = fread(&kmer_n, 8, 1, fp_num);
    temp = fread(&off_n, 8, 1, fp_num);
    temp = fread(&pu_n, 8, 1, fp_num);
    temp = fread(&ref_seq_n, 8, 1, fp_num);
    // printf("us_n\tusf_n\tue_n\tup_n\tupp_n\thash_n\tkmer_n\toff_n\tpu_n\tref_seq_n\n");
    // printf("%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\n",us_n,usf_n,ue_n,up_n,upp_n,hash_n,kmer_n,off_n,pu_n,ref_seq_n);
#else
    fread(&us_n,4,1,fp_num);
    fread(&usf_n,4,1,fp_num);
    fread(&ue_n, 4, 1, fp_num);
    fread(&up_n, 4, 1, fp_num);
    fread(&upp_n, 4, 1, fp_num);
    fread(&hash_n, 4, 1, fp_num);
    fread(&kmer_n, 4, 1, fp_num);
    fread(&off_n, 4, 1, fp_num);
    fread(&pu_n, 4, 1, fp_num);
    fread(&ref_seq_n, 4, 1, fp_num);
    // printf("us_n\tusf_n\tue_n\tup_n,upp_n,hash_n\tkmer_n\toff_n\tpu_n\tref_seq_n\n");
    // printf("%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",us_n,usf_n,ue_n,up_n,upp_n,hash_n,kmer_n,off_n,pu_n,ref_seq_n);
#endif

    fclose(fp_num);
    //*****************************************************
    // read ref seq file
    // printf("Load ref seq..............................\n");
    char ref_seq[ROUTE_LENGTH_MAX] = {0};
    strcpy(ref_seq, index_dir);
    strcat(ref_seq, "ref.seq");

    // char *ref_seq = "/home/ydliu/deBGA/dm/ref.seq";
    FILE *fp_ref_seq = fopen(ref_seq, "rb");
    if (fp_ref_seq == NULL)
    {
        printf("File error opening the seq file\n");
        exit (1);
    }
    fseek(fp_ref_seq, 0, SEEK_END);// 将文件指针移至文件末尾
    file_size = ftell(fp_ref_seq);//obtain文件的长度
    rewind(fp_ref_seq);//将文件指针重新指向文件的开头  =fseek(fp_ref_seq,0,SEEK_SET)
    ref_seq_n = file_size;
    // printf("the reference seq ref_seq_n = %"PRId64"\n", ref_seq_n);

    /*
        buffer_ref_seq:store reference seq;load from dm file ref.seq
    */
    buffer_ref_seq = (uint64_t* )calloc(ref_seq_n + 536, 1);//536 = (2048 >> 5 + 3) << 3 ???uint64_t  1
    a_size = ref_seq_n >> 3;
    result_ref_seq = fread(buffer_ref_seq, 8, a_size, fp_ref_seq);

    // printf("a_size = %"PRId64",result_ref_seq = %"PRId64"\n",a_size,result_ref_seq );
    
    if (result_ref_seq != a_size)
    {
        printf("Reading error\n");
        exit (3);
    }

    fclose(fp_ref_seq);
    //***********************************************
    //read input unipath seq file
    // printf("Load unipath seq.............................\n");
    char uniseq_b[ROUTE_LENGTH_MAX] = {0};
    strcpy(uniseq_b, index_dir);
    strcat(uniseq_b, "unipath.seqb");

    // char *uniseq_b = "/home/ydliu/deBGA/dm/unipath.seqb";
    FILE *fp_us_b = fopen (uniseq_b, "rb" );
    if (fp_us_b == NULL)
    {
        printf( "File error opening the unipath seq file\n");
        exit (1);
    }

    fseek(fp_us_b, 0, SEEK_END);// non-portable
    file_size = ftell(fp_us_b);
    rewind(fp_us_b);

    us_n = file_size;
    // printf("the unipath seq us_n = %"PRId64"\n", us_n);
    /*
        buffer_seq:store unipath seq(concatenate all unipaths's seq);load from dm file unipath.seqb
    */
    buffer_seq = (uint64_t* ) malloc (us_n);
    if (buffer_seq == NULL)
    {
        printf ("Memory error buffer_seq\n");
        exit (2);
    }
    a_size = us_n >> 3;
    result_seq = fread (buffer_seq, 8, a_size, fp_us_b);
    // printf("a_size = %"PRId64",result_seq = %"PRId64"\n",a_size,result_seq );
    if (result_seq != a_size)
    {
        printf ("Reading error");
        exit (3);
    }

    fclose(fp_us_b);
    // //************************************************************8
    // //read input unipath offset file
    // printf("Load unipath offset..............................\n");
    char uniseqf_b[ROUTE_LENGTH_MAX] = {0};
    strcpy(uniseqf_b, index_dir);
    strcat(uniseqf_b, "unipath.seqfb");

    // char * uniseqf_b = "/home/ydliu/deBGA/dm/unipath.seqfb";
    FILE *fp_usf_b = fopen(uniseqf_b, "rb");
    if(fp_usf_b == NULL)
    {
        printf ("File error opening the unipath offset file\n");
        exit (1);
    }

    fseek(fp_usf_b, 0, SEEK_END);// non-portable
    file_size = ftell(fp_usf_b);
    rewind(fp_usf_b);

    usf_n = file_size;
    // printf("the unipath's offset on unipath seq usf_n = %"PRId64"\n", usf_n);
    // allocate memory to contain the whole file:
    // copy the file into the buffer:
    /*
        buffer_seqf:store unipath's offset on unipath seq;load from dm file unipath.seqfb
    */
    buffer_seqf = (uint64_t* ) malloc (usf_n);
    if (buffer_seqf == NULL)
    {
        printf("Memory error buffer_seqf");
        exit (2);
    }
    a_size = (usf_n >> 3);
    result_seqf = fread (buffer_seqf, 8, a_size, fp_usf_b);
    // printf("a_size = %"PRId64",result_seqf = %"PRId64"\n",a_size,result_seqf );
    if (result_seqf != a_size)
    {
        printf ("Reading error");
        exit (3);
    }

    fclose(fp_usf_b);

    // //********************************************************
    //read input unipath position file
    // printf("Load unipath position...........................\n");
    char unipos[ROUTE_LENGTH_MAX] = {0};
    strcpy(unipos, index_dir);
    strcat(unipos, "unipath.pos");

    // char *unipos = "/home/ydliu/deBGA/dm/unipath.pos";
    FILE *fp_up = fopen(unipos, "rb");
    if(fp_up == NULL)
    {
        printf("File error opening the unipath position file\n");
        exit (1);
    }

    fseek(fp_up, 0, SEEK_END);// non-portable
    file_size = ftell(fp_up);
    rewind(fp_up);
    up_n = file_size;
    // printf("unipath's set of positions on reference up_n = %"PRId64"\n", up_n);
    /*
        buffer_p: store every unipath's set of positions on reference;load from dm file unipath.pos
    */
    // copy the file into the buffer:
    buffer_p = (uint64_t* ) malloc (up_n);
    if (buffer_p == NULL)
    {
        printf ("Memory error buffer_p");
        exit (2);
    }
    a_size = (up_n >> 3);
    result_p = fread (buffer_p, 8, a_size, fp_up);
    // printf("a_size = %"PRId64",result_p = %"PRId64"\n",a_size,result_p );
    if (result_p != a_size)
    {
        printf ("Reading error");
        exit (3);
    }
    
    fclose(fp_up);

    // //**************************************************************************
    //read input unipath position point file
    // printf( "Load unipath position point.....................\n");
    char uniposp[ROUTE_LENGTH_MAX] = {0};
    strcpy(uniposp, index_dir);
    strcat(uniposp, "unipath.posp");

    // char *uniposp = "/home/ydliu/deBGA/dm/unipath.posp";
    FILE *fp_upp = fopen(uniposp, "rb");
    if(fp_upp == NULL)
    {
        printf ("File error opening the unipath position point file\n");
        exit (1);
    }
    
    fseek(fp_upp, 0, SEEK_END);// non-portable
    file_size = ftell(fp_upp);
    rewind(fp_upp);

    upp_n = file_size;
    // printf("unipath's pointer to array buffer_p upp_n = %"PRId64"\n", upp_n);

    /*
        buffer_pp: unipath's pointer to array buffer_p;load from dm file unipath.posp
    */
    buffer_pp = (uint64_t* ) malloc (upp_n);
    if (buffer_pp == NULL)
    {
        printf ("Memory error buffer_pp");
        exit (2);
    }
    a_size = (upp_n >> 3);
    result_pp = fread (buffer_pp, 8, a_size, fp_upp);
    // printf("a_size = %"PRId64",result_pp = %"PRId64"\n",a_size,result_pp );
    if (result_pp != a_size)
    {
        printf ("Reading error");
        exit (3);
    }
    fclose(fp_upp);

    // //****************************************************************************
    // //read input unipath hash file
    // printf("Load unipath hash..............................\n");
    char unihash_g[ROUTE_LENGTH_MAX] = {0};
    strcpy(unihash_g, index_dir);
    strcat(unihash_g, "unipath_g.hash");

    // char *unihash_g = "/home/ydliu/deBGA/dm/unipath_g.hash";
    FILE *fp_hash = fopen(unihash_g, "rb");
    if(fp_hash == NULL)
    {
        printf ("File error opening the graph hash file\n");
        exit (1);
    }
    
    fseek(fp_hash, 0, SEEK_END);// non-portable
    file_size = ftell(fp_hash);
    rewind(fp_hash);

    hash_n = file_size;
    // printf("kmer's hash part hash_n = %"PRId64"\n", hash_n);
    
    /*
        buffer_hash_g: store kmer's hash part; load from dm file unipath_g.hash
    */
    buffer_hash_g = (uint64_t* ) malloc (hash_n);
    if (buffer_hash_g == NULL)
    {
        printf("Memory error");
        exit(2);
    }
    a_size = (hash_n >> 3);
    result_hash_g = fread (buffer_hash_g, 8, a_size, fp_hash);
    // printf("a_size = %"PRId64",result_hash_g = %"PRId64"\n",a_size,result_hash_g );
    if (result_hash_g != a_size)
    {
        printf("Reading error");
        exit(3);
    }
    fclose(fp_hash);

    //***************************************************************************************
    //read input graph kmer file
    // printf("Load unipath kmer.........................\n");
    char unikmer_g[ROUTE_LENGTH_MAX] = {0};
    strcpy(unikmer_g, index_dir);
    strcat(unikmer_g, "unipath_g.kmer");

    // char *unikmer_g = "/home/ydliu/deBGA/dm/unipath_g.kmer";
    FILE *fp_kmer = fopen(unikmer_g, "rb");
    if(fp_kmer == NULL)
    {
        printf ("File error opening the graph hash file\n");
        exit (1);
    }
    
    fseek(fp_kmer, 0, SEEK_END);// non-portable
    file_size = ftell(fp_kmer);
    rewind(fp_kmer);

    kmer_n = file_size;
    // printf("kmer's kmer part kmer_n = %"PRId64"\n", kmer_n);
    
    /*
        buffer_kmer_g:store kmer's kmer part;load from dm file unipath_g.kmer
    */
    buffer_kmer_g = (uint32_t* ) malloc (kmer_n);
    if (buffer_kmer_g == NULL)
    {
        printf ("Memory error buffer_kmer_g");
        exit (2);
    }

    a_size = (kmer_n >> 2);
    result_kmer_g = fread (buffer_kmer_g, 4, a_size, fp_kmer);
    // printf("a_size = %"PRId64",result_kmer_g = %"PRId64"\n",a_size,result_kmer_g );
    if (result_kmer_g != a_size)
    {
        printf ("Reading error buffer_kmer_g");
        exit (3);
    }

    fclose(fp_kmer);

    //*****************************************************************************
    //read input graph off file
    // printf("Load unipath off......................\n");
    char unioff_g[ROUTE_LENGTH_MAX] = {0};
    strcpy(unioff_g, index_dir);
    strcat(unioff_g, "unipath_g.offset");

    // char *unioff_g = "/home/ydliu/deBGA/dm/unipath_g.offset";
    FILE *fp_off = fopen(unioff_g, "rb");
    if(fp_off == NULL)
    {
        printf("File error opening the graph hash file\n");
        exit (1);
    }
    
    fseek(fp_off, 0, SEEK_END);// non-portable
    file_size = ftell(fp_off);
    rewind(fp_off);

    off_n = file_size;
    // printf("kmer's offset on unipath seq off_n = %"PRId64"\n", off_n);

    /*
        buffer_off_g: store kmer's offset on unipath seq; load from dm file unipath_g.offset
    */
    // copy the file into the buffer:
    buffer_off_g = (uint64_t* ) malloc (off_n);
    if (buffer_off_g == NULL)
    {
        printf ("Memory error");
        exit (2);
    }
    a_size = (off_n >> 3);
    result_off_g = fread (buffer_off_g, 8, a_size, fp_off);
    // printf("a_size = %"PRId64",result_off_g = %"PRId64"\n",a_size,result_off_g );

    if (result_off_g != a_size)
    {
        printf ("Reading error");
        exit (3);
    }

    fclose(fp_off);

    //********************************************************************************************
    //read chr names and length from unipath.chr
    // printf("read chr name and length.....................\n");
    char unichr[ROUTE_LENGTH_MAX] = {0};
    strcpy(unichr, index_dir);
    strcat(unichr, "unipath.chr");
    
    // char *unichr = "/home/ydliu/deBGA/dm/unipath.chr"; 
    FILE *fp_chr = fopen (unichr, "r" );
    if (fp_chr == NULL)
    {
        printf ("File error opening the chr file\n");
        exit (1);
    }
    //染色体条数
    uint32_t chr_line_n = 0;

    temp = fscanf(fp_chr,"%s",chr_line_content); //按照格式从文件获取内容，将文件中的字符串赋值给chr_line_content
    while(!feof(fp_chr))
    {
        //printf("chr_line_content = %s\n",chr_line_content);

        if ((chr_line_n & 0X1) == 0)  //交替获取染色体名和染色体长度
        {
            strcpy(chr_names[chr_file_n],chr_line_content);
            // chr_names[chr_file_n][strlen(chr_names[chr_file_n])] = '\0';
            // printf("chr_names[%u] = %s\n",chr_file_n, chr_names[chr_file_n] );
        }else{
            sscanf(chr_line_content, "%"PRId64"", &chr_end_n[chr_file_n]);

            // printf("chr_end_n[%u] = %"PRId64"\n", chr_file_n,chr_end_n[chr_file_n]);

            chr_file_n++;
        }
        fflush(stdout);
        
        chr_line_n++;
        temp = fscanf(fp_chr,"%s",chr_line_content); //按照格式从文件获取内容，将文件中的字符串赋值给chr_line_content
    }
    chr_end_n[0] = START_POS_REF + 1; //START_POS_REF = 2048 record the start position of reference

    strcpy(chr_names[chr_file_n], "*");
    reference_len = chr_end_n[chr_file_n - 1];
    // printf("total lengths of reference is %"PRId64"\n",chr_end_n[chr_file_n-1]);
    // printf("chr_line_n = %u, chr_file_n =%u \n",chr_line_n,chr_file_n );
    fclose(fp_chr);

    return 0;
}

