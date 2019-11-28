
#include <stdlib.h>
#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <dirent.h>
#include <inttypes.h>

#include "index_build.h"
#include "bit_operation.h"
#include "load_input.h"

int index_build(int argc, char *argv[])
{
	if(load_input_index(argc, argv) == 1)	return 1;
	
	printf("kmer size to build index: %u\n", k_t);
	fflush(stdout);
    //load_reffile_kmer();
	
#ifdef	LAST_DEBUG	
	load_reffile_kmer_fa();
#endif

	file_kmer_qsort();

    return 0;
}

void upper_string(char *string)
{
    while(*string)
    {
        if ( *string >= 'a' && *string <= 'z' )
        {
            *string = *string - 32;
        }
        string++;
    }
}
/*
void load_reffile_kmer()
{
    fp_n = fopen(N_route, "w");
    if(fp_n == NULL)    printf("cannot open N statistical file\n");

    uint32_t n_cnt = 0;
	uint32_t n_tr = 0;
	uint8_t n_end_f = 0;
    uint8_t n_w_f = 0;

    //read one line
    char one_line[REF_FASTA_LINE + 1] = "";
    char kmer_pre[KT_LENGTH_MAX + 1] = "";
    char kmer_com[KT_LENGTH_MAX + REF_FASTA_LINE + 1] = "";
    char kmer[KT_LENGTH_MAX + 1];
    char input = charN;
    char output;

    uint8_t line_i = 0;
    uint8_t kmer_i = 0;
    uint8_t line_l = 0;
    uint8_t l_p = 0;
    uint32_t line_t_c = 0;
    uint32_t i = 0;

    FILE* fpin_ref = NULL;

    DIR *directory_pointer;
    struct dirent *entry;

    char file_ref[ROUTE_LENGTH_MAX];

    //file creating
    uint32_t add = 0;
    char kmer_f[KF_LENGTH_MAX] = "";
    char kmer_s[KT_LENGTH_MAX] = "";

    uint32_t file_n = 0;
    char file_d[KF_LENGTH_MAX + ROUTE_LENGTH_MAX] = "";

    file_n = (1 << (f << 1));
	//modification
	//file_n = ((1 << (f << 1)) << 1);
	
    FILE** file_p = (FILE** )malloc(file_n * sizeof(FILE* ));
	if(file_p == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}
	
    for(i = 0; i < file_n; i++)
        file_p[i] = NULL;

    uint32_t line_tol = 0;
    uint32_t* file_line_cnt = (uint32_t* )calloc(file_n, 4);
	if(file_line_cnt == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}
	
    FILE* file_sta = NULL;
    file_sta = fopen(filename_sta,"w");
    if (file_sta == NULL)
    {
        fputs ("File error of creating statistical file\n",stderr);
        exit(1);
    }
	
#ifdef	UNPIPATH_OFF_K20
	uint32_t write_buff[4];
	uint64_t pos = 0;	
#else
    uint32_t write_buff[3];
	uint32_t pos = 0;
#endif
	
    uint8_t last_kmer_f = 0;
    chr_end_n[0] = 1;
    uint32_t chr_i = 0;

	//ref seq load
	uint64_t ref_seq_buffer[128];
	fp_ref_seq = fopen(ref_seq, "wb");

    if((directory_pointer = opendir(filename_ref))==NULL)
        printf( "Error opening \n ");
    else
    {
        while((entry = readdir(directory_pointer))!=NULL)
        {
            //begin with each chr file
			if((strstr(entry->d_name,".fna")!= NULL) && (strstr(entry->d_name,"NC")!= NULL))
            {
                memset(file_ref, 0, sizeof(file_ref));
                strcpy(file_ref,filename_ref);
                strcat(file_ref, entry->d_name);

                //if file is not empty
                strcpy(chr_names[chr_file_n], entry->d_name);

                //begin loading file and loading kmer into graph
                fpin_ref = fopen(file_ref, "r");
                if (fpin_ref == NULL)
                {
                    fputs ("File error of reading ref file\n",stderr);
                    exit(1);
                }

				l_p = 0;

                while ((!feof(fpin_ref)) && (fgets(one_line, REF_FASTA_LINE + 2, fpin_ref) != NULL))
                {
                    //line_c++;

                    if(strstr(one_line,identifier) != NULL)//>
                    {
                        continue;
                    }

                    line_t_c++;

                    line_l = strlen(one_line);

                    if(one_line[line_l - 1] == '\n')
                    {
                        one_line[line_l - 1] = '\0';
                        line_l--;
                    }

					//ref_load
					for(line_i = 0; line_i < line_l; line_i++)
					{
						ref_seq_buffer[(ref_seq_n & 0Xfff) >> 5] |= (((uint64_t )charToDna5_N2[(uint8_t )one_line[line_i]]) << ((31 - ((ref_seq_n & 0Xfff) & 0X1f)) << 1));

						++ref_seq_n;
						if((ref_seq_n & 0Xfff) == 0)
						{
							fwrite(ref_seq_buffer, 8, 128, fp_ref_seq);
							memset(ref_seq_buffer, 0, 1024);
						}
					}

                    memset(kmer_com, 0, KT_LENGTH_MAX + REF_FASTA_LINE);

                    seq_exact_2(kmer_pre,kmer_com,kmer_i,0)
                    seq_exact_2(one_line,kmer_com,kmer_i,l_p)

                    for(line_i = 0; line_i < line_l - k_t + l_p; line_i++,pos++)
                    {
                        seq_exact(kmer_com,kmer,kmer_i,line_i,line_i + k_t)

                        if((strstr(kmer, char_N) != NULL) || (strstr(kmer, char_n) != NULL))
                        {
							memset(kmer_pre, 0, KT_LENGTH_MAX + 1);
							kmer_pre[0] = '\0';

							input = charN;

							//N number
							++n_cnt;
							n_w_f = 1;

                            continue;
                        }

						//N position and number
						if(n_w_f)
						{
							if(n_end_f)
							{
								n_tr = n_cnt;
								n_end_f = 0;
							}
							else	n_tr = n_cnt - (k_t - 1);

							n_cnt = 0;
							n_w_f = 0;
						}

                        output = kmer_com[line_i + k_t];

                        //do something here, kmer[], pos + 1
                        strncpy(kmer_f, (const char* )kmer, f);
                        kmer_address(kmer_f, f, kmer_i, add)

                        upper_string(kmer_f);

                        memset(file_d, 0, KF_LENGTH_MAX + ROUTE_LENGTH_MAX);
                        strcpy(file_d, (const char* )filename_div);
                        strcat(file_d, (const char* )kmer_f);
                        strcat(file_d, suff);

                        if(file_p[add] == NULL)
                        {
                            file_p[add] = fopen(file_d, "wb");
                            if (file_p[add] == NULL)
                            {
                                fputs ("File error of creating divisional file\n",stderr);
                                exit(1);
                            }
                        }

                        seq_exact(kmer, kmer_s, kmer_i, f, strlen(kmer));

                        //binary file
                        //must assign in this order
                        kmer_bit32a(kmer_s, k_t - f, kmer_i, write_buff, 3)
                        write_buff[1] |= (input << 24);
                        write_buff[1] |= (output << 16);
#ifdef	UNPIPATH_OFF_K20
						write_buff[3] = (pos + 1) >> 32;
						write_buff[0] = (pos + 1) & 0Xffffffff;
						fwrite(write_buff, 4, 4, file_p[add]);
#else
                        write_buff[0] = pos + 1;
						fwrite(write_buff, 4, 3, file_p[add]);
#endif
                        
                        //

                        file_line_cnt[add]++;
                        line_tol++;
                        //

                        input = kmer_com[line_i];
                    }

                    seq_exact(kmer_com, kmer_pre, kmer_i, line_l - k_t + l_p, line_l+l_p)

                    l_p = k_t;
                };

                //write the last kmer
                if(kmer_pre[0] != '\0')
                {
                    seq_exact(kmer_pre,kmer,kmer_i,0,k_t)

                    last_kmer_f = 1;

                    if((strstr(kmer, char_N) != NULL) || (strstr(kmer, char_n) != NULL))
                    {
                        input = charN;
                        last_kmer_f = 0;

						//N number
						++n_cnt;
						n_w_f = 1;
						n_end_f = 1;
                    }

                    if(last_kmer_f)
                    {
                        output = charN;

                        //do something here, kmer[], pos + 1
                        strncpy(kmer_f, (const char* )kmer, f);
                        kmer_address(kmer_f, f, kmer_i, add)

                        upper_string(kmer_f);

                        memset(file_d, 0, KF_LENGTH_MAX + ROUTE_LENGTH_MAX);
                        strcpy(file_d, (const char* )filename_div);
                        strcat(file_d, (const char* )kmer_f);
                        strcat(file_d, suff);

                        if(file_p[add] == NULL)
                        {
                            file_p[add] = fopen(file_d, "wb");
                            if (file_p[add] == NULL)
                            {
                                fputs ("File error of creating divisional file\n",stderr);
                                exit(1);
                            }
                        }

                        seq_exact(kmer, kmer_s, kmer_i, f, strlen(kmer));

                        //binary file
                        //must assign in this order
                        kmer_bit32a(kmer_s, k_t - f, kmer_i, write_buff, 3)
                        write_buff[1] |= (input << 24);
                        write_buff[1] |= (output << 16);
                        write_buff[0] = pos + 1;

#ifdef	UNPIPATH_OFF_K20
						write_buff[3] = (pos + 1) >> 32;
						write_buff[0] = (pos + 1) & 0Xffffffff;
						fwrite(write_buff, 4, 4, file_p[add]);
#else
                        write_buff[0] = pos + 1;
						fwrite(write_buff, 4, 3, file_p[add]);
#endif
                        //

                        file_line_cnt[add]++;
                        line_tol++;
                    }
                }

				memset(kmer_pre, 0, KT_LENGTH_MAX + 1);
                kmer_pre[0] = '\0';

                input = charN;

				pos += k_t;

                //end last kmer

                fclose(fpin_ref);

                chr_end_n[chr_file_n] = pos + 1;
                chr_file_n++;

                //
                printf("Has finished loading %s\n", entry->d_name);//, pos+1
			}
        }
        closedir(directory_pointer);
    }

	fwrite(ref_seq_buffer, 8, ((ref_seq_n & 0Xfff) >> 5) + 1, fp_ref_seq);
	fclose(fp_ref_seq);

	printf("number chars of ref seq: %"PRId64"\n", (((ref_seq_n >> 5) + 1) << 3));

	FILE* fp_chr = fopen(unichr,"w");
	if (fp_chr == NULL)
    {
        fputs ("File error of creating chr info file\n",stderr);
        exit(1);
    }

#ifdef	UNPIPATH_OFF_K20
	for(chr_i = 1; chr_i < chr_file_n; chr_i++)
	{
		fprintf(fp_chr,"%s\t%"PRId64"\n",chr_names[chr_i], chr_end_n[chr_i]);
	}
	for(i = 0; i < file_n; i++)
        fprintf(file_sta,"%"PRId64"\n",file_line_cnt[i]);
    fprintf(file_sta,"\n%"PRId64"",line_tol);
    printf("total lines of ref file: %"PRId64"\n",line_tol);
#else
	for(chr_i = 1; chr_i < chr_file_n; chr_i++)
		fprintf(fp_chr,"%s\t%u\n",chr_names[chr_i], chr_end_n[chr_i]);
	for(i = 0; i < file_n; i++)
        fprintf(file_sta,"%u\n",file_line_cnt[i]);
    fprintf(file_sta,"\n%u",line_tol);
    printf("total lines of ref file: %u\n",line_tol);
#endif
	
	fclose(fp_chr);

	for(i = 0; i < file_n; i++)
    {
        if(file_p[i] != NULL)
            fclose(file_p[i]);
    }

    if(file_p)	free(file_p);
    if(file_line_cnt)	free(file_line_cnt);

    fclose(file_sta);
	fclose(fp_n);

}

*/
void load_reffile_kmer_fa()
{
    fp_n = fopen(N_route, "w");
    if(fp_n == NULL)    printf("cannot open N statistical file\n");

	ref_seq_n = START_POS_REF;

    uint32_t n_cnt = 0;
	uint32_t n_tr = 0;

	uint8_t n_end_f = 0;
    uint8_t n_w_f = 0;

    //read one line
    char one_line[REF_FASTA_LINE + 1] = "";
    char kmer_pre[KT_LENGTH_MAX + 1] = "";
    char kmer_com[KT_LENGTH_MAX + REF_FASTA_LINE + 1] = "";
    char kmer[KT_LENGTH_MAX + 1];
    char input = charN;
    char output;

    uint8_t line_i = 0;
    uint32_t kmer_i = 0;
    uint8_t line_l = 0;
    uint8_t l_p = 0;
    uint32_t line_t_c = 0;
    uint32_t i = 0;

    FILE* fpin_ref = NULL;

    //file creating
    uint32_t add = 0;
    char kmer_f[KF_LENGTH_MAX] = "";
    char kmer_s[KT_LENGTH_MAX] = "";

    uint32_t file_n = 0;
    char file_d[KF_LENGTH_MAX + ROUTE_LENGTH_MAX] = "";

    file_n = ((uint32_t )1 << (f << 1));
	//modification
	//file_n = ((1 << (f << 1)) << 1);
	printf("%u files to be written in folder div\n", file_n);
	fflush(stdout);
	
    FILE** file_p = (FILE** )malloc(file_n * sizeof(FILE* ));
	if(file_p == NULL)
	{
		printf("Fail to allocate memory for file_p\n");
		fflush(stdout);
		exit(1);
	}
    for(i = 0; i < file_n; i++)
        file_p[i] = NULL;

    FILE* file_sta = NULL;
    file_sta = fopen(filename_sta,"w");
    if (file_sta == NULL)
    {
        fputs ("File error of creating statistical file\n",stderr);
		fflush(stdout);
        exit(1);
    }

#ifdef	UNPIPATH_OFF_K20
	uint32_t write_buff[4];
	uint64_t pos = START_POS_REF;
	uint64_t pos_tr = 0;	
	uint64_t line_tol = 0;
    uint64_t* file_line_cnt = (uint64_t* )calloc(file_n, 8);
#else
    uint32_t write_buff[3];
	uint32_t pos = START_POS_REF;
	uint32_t pos_tr = 0;
	uint32_t line_tol = 0;
    uint32_t* file_line_cnt = (uint32_t* )calloc(file_n, 4);
#endif
	if(file_line_cnt == NULL)
	{
		printf("Fail to allocate memory for file_line_cnt\n");
		fflush(stdout);
		exit(1);
	}
    //read file name under path
    uint8_t last_kmer_f = 0;
    chr_end_n[0] = 1;
    uint32_t chr_i = 0;

	//ref seq load
	uint64_t ref_seq_buffer[128];

    //this is bug for previous deBGA, init
    memset(ref_seq_buffer, 0, 1024);
	fp_ref_seq = fopen(ref_seq, "wb");

    //begin loading file and loading kmer into graph
    fpin_ref = fopen(filename_ref, "r");
    if (fpin_ref == NULL)
    {
        fputs ("File error of reading ref file and wrong ref file name or route\n",stderr);
        exit(1);
    }

    l_p = 0;

	uint8_t fisrt_chr = 1;
	uint32_t chr_name_n = 0;
	uint32_t chr_posend_n = 0;
	
#ifdef	CHR_NAME_SPLIT	
	char* pch = NULL;
	char* saveptr = NULL;
	char tmp_chr_name[REF_FASTA_LINE + 1];
#endif

    while ((!feof(fpin_ref)) && (fgets(one_line, REF_FASTA_LINE + 2, fpin_ref) != NULL))
    {
        if(strstr(one_line,identifier) != NULL)
        {
			one_line[strlen(one_line) - 1] = '\0';
			
#ifdef	CHR_NAME_SPLIT
			strcpy(tmp_chr_name, one_line + 1);
			pch = strtok_r(tmp_chr_name, " ", &saveptr);

			strcpy(chr_names[chr_name_n++], pch);
#else
			strcpy(chr_names[chr_name_n++], one_line + 1);
#endif

			if(fisrt_chr)
			{
				fisrt_chr = 0;
				continue;
			}

            //write the last kmer
            if(kmer_pre[0] != '\0')
            {
                seq_exact(kmer_pre,kmer,kmer_i,0,k_t)

                last_kmer_f = 1;

                if((strstr(kmer, char_N) != NULL) || (strstr(kmer, char_n) != NULL))
                {
                    input = charN;
                    last_kmer_f = 0;

                    //N number
                    ++n_cnt;
                    n_w_f = 1;
                    n_end_f = 1;
                }

                if(last_kmer_f)
                {
                    output = charN;

                    //do something here, kmer[], pos + 1
                    strncpy(kmer_f, (const char* )kmer, f);
                    kmer_address(kmer_f, f, kmer_i, add)

                    upper_string(kmer_f);

                    memset(file_d, 0, KF_LENGTH_MAX + ROUTE_LENGTH_MAX);
                    strcpy(file_d, (const char* )filename_div);
                    strcat(file_d, (const char* )kmer_f);
                    strcat(file_d, suff);

                    if(file_p[add] == NULL)
                    {
                        file_p[add] = fopen(file_d, "wb");
                        if (file_p[add] == NULL)
                        {
                            fputs ("File error of creating divisional file\n",stderr);
                            exit(1);
                        }
                    }

                    seq_exact(kmer, kmer_s, kmer_i, f, strlen(kmer));
			
                    //binary file
                    //must assign in this order
                    kmer_bit32a(kmer_s, k_t - f, kmer_i, write_buff, 3)
                    write_buff[1] |= (input << 24);
                    write_buff[1] |= (output << 16);

#ifdef	UNPIPATH_OFF_K20
					write_buff[3] = (pos + 1) >> 32;
					write_buff[0] = (pos + 1) & 0Xffffffff;
					fwrite(write_buff, 4, 4, file_p[add]);
#else
                    write_buff[0] = pos + 1;
					fwrite(write_buff, 4, 3, file_p[add]);
					//fwrite(write_buff, sizeof(write_buff ), 1, file_p[add]);
#endif

#ifdef	DEBUG_INPUT
					if(write_buff[2] == 8258096)
						printf("\n\n2. %s:\nline: %s\n%s\n%u\n\n\n", tmp_chr_name, one_line, kmer, (uint16_t )write_buff[1]);
#endif
					
                    //
                    file_line_cnt[add]++;
                    line_tol++;
                }
            }

            memset(kmer_pre, 0, KT_LENGTH_MAX + 1);
            kmer_pre[0] = '\0';

            input = charN;

            pos += k_t;

            //end last kmer

            chr_end_n[chr_posend_n++] = pos + 1;

			printf("Has finished loading %s\n", chr_names[chr_name_n - 2]);
			fflush(stdout);
			
			l_p = 0;
        }
        else
        {
            line_t_c++;

            line_l = strlen(one_line);

            if(one_line[line_l - 1] == '\n')
            {
                one_line[line_l - 1] = '\0';
                line_l--;
            }

            //ref_load
            for(line_i = 0; line_i < line_l; line_i++)
            {
				//modification
				one_line[line_i] = charTochar[(uint8_t )one_line[line_i]];
				
				ref_seq_buffer[(ref_seq_n & 0Xfff) >> 5] |= (((uint64_t )charToDna5_N2[(uint8_t )one_line[line_i]]) << ((31 - ((ref_seq_n & 0Xfff) & 0X1f)) << 1));

                ++ref_seq_n;
                if((ref_seq_n & 0Xfff) == 0)
                {
                    fwrite(ref_seq_buffer, 8, 128, fp_ref_seq);
                    memset(ref_seq_buffer, 0, 1024);
                }
            }

            memset(kmer_com, 0, KT_LENGTH_MAX + REF_FASTA_LINE);

            seq_exact_2(kmer_pre,kmer_com,kmer_i,0)
            seq_exact_2(one_line,kmer_com,kmer_i,l_p)

            for(line_i = 0; line_i < line_l - k_t + l_p; line_i++,pos++)
            {
                seq_exact(kmer_com,kmer,kmer_i,line_i,line_i + k_t)

                if((strstr(kmer, char_N) != NULL) || (strstr(kmer, char_n) != NULL))
                {
                    memset(kmer_pre, 0, KT_LENGTH_MAX + 1);
                    kmer_pre[0] = '\0';

                    input = charN;

                    //N number
                    ++n_cnt;
                    n_w_f = 1;

                    continue;
                }

                //N position and number
                if(n_w_f)
                {
                    pos_tr = pos + 1;

                    if(n_end_f)
                    {
                        n_tr = n_cnt;
                        n_end_f = 0;
                    }
                    else	n_tr = n_cnt - (k_t - 1);


                    if(n_tr <= 0)
                    {
                    	printf("wrong n %u %u\n", n_tr, n_cnt);
                    	exit(1);
                    }
                    fwrite(&pos_tr, 4, 1, fp_n);
                    fwrite(&n_tr, 4, 1, fp_n);


                    n_cnt = 0;
                    n_w_f = 0;
                }

                output = kmer_com[line_i + k_t];

                //do something here, kmer[], pos + 1
                strncpy(kmer_f, (const char* )kmer, f);
                kmer_address(kmer_f, f, kmer_i, add)

                upper_string(kmer_f);

                memset(file_d, 0, KF_LENGTH_MAX + ROUTE_LENGTH_MAX);
                strcpy(file_d, (const char* )filename_div);
                strcat(file_d, (const char* )kmer_f);
                strcat(file_d, suff);

                if(file_p[add] == NULL)
                {
                    file_p[add] = fopen(file_d, "wb");
                    if (file_p[add] == NULL)
                    {
                        fputs ("File error of creating divisional file\n",stderr);
                        exit(1);
                    }
                }

                seq_exact(kmer, kmer_s, kmer_i, f, strlen(kmer));

                //binary file
                //must assign in this order
                kmer_bit32a(kmer_s, k_t - f, kmer_i, write_buff, 3)
                write_buff[1] |= (input << 24);
                write_buff[1] |= (output << 16);

#ifdef	UNPIPATH_OFF_K20
				write_buff[3] = (pos + 1) >> 32;
				write_buff[0] = (pos + 1) & 0Xffffffff;
				fwrite(write_buff, 4, 4, file_p[add]);
#else
                write_buff[0] = pos + 1;
				fwrite(write_buff, 4, 3, file_p[add]);
				//fwrite(write_buff, sizeof(write_buff ), 1, file_p[add]);
#endif
 
#ifdef	DEBUG_INPUT
				if(write_buff[2] == 8258096)
					printf("\n\n%2. s:\nline: %s\n%s\n%u\n\n\n", tmp_chr_name, one_line, kmer, write_buff[1]);
#endif
                //

                file_line_cnt[add]++;
                line_tol++;
                //

                input = kmer_com[line_i];
            }

            seq_exact(kmer_com, kmer_pre, kmer_i, line_l - k_t + l_p, line_l+l_p)

            l_p = k_t;
        }
    }

	//write the last kmer
    if(kmer_pre[0] != '\0')
    {
        seq_exact(kmer_pre,kmer,kmer_i,0,k_t)

        last_kmer_f = 1;

        if((strstr(kmer, char_N) != NULL) || (strstr(kmer, char_n) != NULL))
        {
            input = charN;
            last_kmer_f = 0;

            //N number
            ++n_cnt;
            n_w_f = 1;
            n_end_f = 1;
        }

        if(last_kmer_f)
        {
            output = charN;

            //do something here, kmer[], pos + 1
            strncpy(kmer_f, (const char* )kmer, f);
            kmer_address(kmer_f, f, kmer_i, add)

            upper_string(kmer_f);

            memset(file_d, 0, KF_LENGTH_MAX + ROUTE_LENGTH_MAX);
            strcpy(file_d, (const char* )filename_div);
            strcat(file_d, (const char* )kmer_f);
            strcat(file_d, suff);

            if(file_p[add] == NULL)
            {
                file_p[add] = fopen(file_d, "wb");
                if (file_p[add] == NULL)
                {
                    fputs ("File error of creating divisional file\n",stderr);
                    exit(1);
                }
            }

            seq_exact(kmer, kmer_s, kmer_i, f, strlen(kmer));


            //binary file
            //must assign in this order
            kmer_bit32a(kmer_s, k_t - f, kmer_i, write_buff, 3)
            write_buff[1] |= (input << 24);
            write_buff[1] |= (output << 16);
			
#ifdef	UNPIPATH_OFF_K20
			write_buff[3] = (pos + 1) >> 32;
			write_buff[0] = (pos + 1) & 0Xffffffff;
			fwrite(write_buff, 4, 4, file_p[add]);
#else
            write_buff[0] = pos + 1;
			fwrite(write_buff, 4, 3, file_p[add]);
			//fwrite(write_buff, sizeof(write_buff ), 1, file_p[add]);
#endif

#ifdef	DEBUG_INPUT
			if(write_buff[2] == 8258096)
				printf("\n\n2. %s:\nline: %s\n%s\n%u\n\n\n", tmp_chr_name, one_line, kmer, write_buff[1]);
#endif
            //

            file_line_cnt[add]++;
            line_tol++;
        }
    }

    pos += k_t;
    chr_end_n[chr_posend_n++] = pos + 1;
	
	printf("Has finished loading %s\n", chr_names[chr_name_n - 1]);
#ifdef	UNPIPATH_OFF_K20
	printf("total number of lines: %"PRId64"\n",line_tol);
	printf("total number of positions on ref seq: %"PRId64"\n", pos+1);
#else
	printf("total number of lines: %u\n",line_tol);
	printf("total number of positions on ref seq: %u\n", pos+1);
#endif
	fclose(fpin_ref);

	fwrite(ref_seq_buffer, 8, ((ref_seq_n & 0Xfff) >> 5) + 1, fp_ref_seq);
	fclose(fp_ref_seq);

	printf("total number of chars to allocate for ref seq: %"PRId64"\n", (((ref_seq_n >> 5) + 1) << 3));

	fflush(stdout);
	
	FILE* fp_chr = fopen(unichr,"w");
	if (fp_chr == NULL)
    {
        fputs ("File error of creating chr info file\n",stderr);
        exit(1);
    }

#ifdef	UNPIPATH_OFF_K20
	for(chr_i = 0; chr_i < chr_posend_n; chr_i++)
	{
		fprintf(fp_chr,"%s\n%"PRId64"\n",chr_names[chr_i], chr_end_n[chr_i]);
	}
	fclose(fp_chr);

    for(i = 0; i < file_n; i++)
        fprintf(file_sta,"%"PRId64"\n",file_line_cnt[i]);
    fprintf(file_sta,"\n%"PRId64"",line_tol);
	fflush(file_sta);
#else
	for(chr_i = 0; chr_i < chr_posend_n; chr_i++)
	{
		fprintf(fp_chr,"%s\n%u\n",chr_names[chr_i], chr_end_n[chr_i]);
	}
	fclose(fp_chr);

    for(i = 0; i < file_n; i++)
        fprintf(file_sta,"%u\n",file_line_cnt[i]);
    fprintf(file_sta,"\n%u",line_tol);
	fflush(file_sta);
#endif

	for(i = 0; i < file_n; i++)
    {
        if(file_p[i] != NULL)
            fclose(file_p[i]);
    }

    free(file_p);

    free(file_line_cnt);

    fclose(file_sta);

	fclose(fp_n);

}


int compare (const void * a, const void * b)
{
    k_p* kp1 = (k_p* )a;
    k_p* kp2 = (k_p* )b;

    uint16_t kmerf1 = (uint16_t )((kp1->kp)[1]);
    uint16_t kmerf2 = (uint16_t )((kp2->kp)[1]);
    uint32_t kmers1 = (uint32_t )((kp1->kp)[2]);
    uint32_t kmers2 = (uint32_t )((kp2->kp)[2]);

    if (kmerf1 > kmerf2)
        return 1;
    else if (kmerf1 < kmerf2)
        return -1;
    else
    {
        if (kmers1 > kmers2)
            return 1;
        else if (kmers1 < kmers2)
            return -1;
        else
        {
#ifdef	UNPIPATH_OFF_K20	
			/*
			uint64_t tmp_kp_pos1 = (kp1->kp)[3];
			uint64_t tmp_kp_pos2 = (kp2->kp)[3];
			tmp_kp_pos1 = ((tmp_kp_pos1 << 32) | ((kp1->kp)[0]));
			tmp_kp_pos2 = ((tmp_kp_pos2 << 32) | ((kp2->kp)[0]));

            if (tmp_kp_pos1 > tmp_kp_pos2)
                return 1;
            else if (tmp_kp_pos1 < tmp_kp_pos2)
                return -1;
            else
            {
                printf("same pos\n");
				exit(2);
                return 0;
            }
			*/
			
			if ((kp1->kp)[3] > (kp2->kp)[3])
                return 1;
            else if ((kp1->kp)[3] < (kp2->kp)[3])
                return -1;
            else
            {
				if ((kp1->kp)[0] > (kp2->kp)[0])
					return 1;
				else if ((kp1->kp)[0] < (kp2->kp)[0])
					return -1;
				else
				{
					printf("same pos\n");
					fflush(stdout);
					return 0;
				} 
            }
#else					
			if ((kp1->kp)[0] > (kp2->kp)[0])
                return 1;
            else if ((kp1->kp)[0] < (kp2->kp)[0])
                return -1;
            else
            {
                printf("same pos\n");
				fflush(stdout);
                return 0;
            }
#endif
        }
    }
}

int compare_pu (const void * a, const void * b)
{
    p_u* pu1 = (p_u *)a;
    p_u* pu2 = (p_u *)b;

    if(pu1->pos > pu2->pos)
        return 1;
    else if(pu1->pos < pu2->pos)
        return -1;
    else	return 0;
}

int compare_us (const void * a, const void * b)
{
    k_u* ku1 = (k_u *)a;
    k_u* ku2 = (k_u *)b;

    if(ku1->kmer > ku2->kmer)
        return 1;
    else if(ku1->kmer < ku2->kmer)
        return -1;
    else
    {
        return 0;
    }
}

uint32_t file_kmer_qsort()
{
#ifdef	LAST_DEBUG
    uint8_t k_i = 0;
    char char_s[32] = "";
    char char_l[64] = "";
    //uint8_t i_kmer = 0;

    uint32_t i = 0;
    uint8_t kmer_i = 0;

    uint32_t file_n = 0;
    file_n = (1 << (f << 1));
    char file_sta[ROUTE_LENGTH_MAX + KF_LENGTH_MAX];

	FILE* file_s = NULL;
	file_s = fopen(filename_sta,"r");
	if(file_s == NULL)	printf("Error of opening graph file\n");

    long int line_cnt = 0;
    uint32_t i_line = 0;
	
#ifdef	UNPIPATH_OFF_K20
	uint64_t line_tol = 0;
    uint64_t* file_line_cnt = (uint64_t* )calloc(file_n, 8);
	if(file_line_cnt == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}
#else
    uint32_t line_tol = 0;
    uint32_t* file_line_cnt = (uint32_t* )calloc(file_n, 4);
	if(file_line_cnt == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}
#endif
	
    while(!feof(file_s))
    {
        fgets(file_sta, 32, file_s);
        line_cnt = atol(file_sta);
#ifdef	UNPIPATH_OFF_K20
		if(i_line < file_n)
            file_line_cnt[i_line] = (uint64_t )line_cnt;
        if(i_line == file_n + 1)
            line_tol = (uint64_t )line_cnt;
#else
        if(i_line < file_n)
            file_line_cnt[i_line] = (uint32_t )line_cnt;
        if(i_line == file_n + 1)
            line_tol = (uint32_t )line_cnt;
#endif
        i_line++;
    }

    //begin qsorting
    char file_div[KF_LENGTH_MAX];
    char file_open[ROUTE_LENGTH_MAX + KF_LENGTH_MAX];
    FILE** file_p = (FILE** )malloc(file_n * sizeof(FILE* ));
	if(file_p == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}
	
    for(i = 0; i < file_n; i++)
        file_p[i] = NULL;

    k_p* k_p_file = NULL;

    //load into hash graph declaration
    uint32_t k_p_i = 0;
    uint64_t kmer_l = 0;
    uint32_t kmer_f = 0;
    uint32_t kmer_s = 0;
    uint64_t kmer_l_p = 0Xffffffffffffffff;
    uint32_t kmer_f_p = 0Xffffffff;
    uint64_t hash_n = ((uint64_t )1 << (k << 1));
	
	printf("hash_n: %"PRId64"\n", hash_n);//268435456
	
    //alloc the memory used to store hash graph
#ifdef	UNPIPATH_OFF_K20
	uint64_t* pos_array = (uint64_t* )calloc(line_tol, 8);
	if(pos_array == NULL)
	{
		printf("fail to allocate memory pos_array\n");
		fflush(stdout);
		exit(1);
	}
	uint64_t* kmer_point = (uint64_t* )calloc(line_tol + 1, 8);
	if(kmer_point == NULL)
	{
		printf("fail to allocate memory kmer_point\n");
		fflush(stdout);
		exit(1);
	}
	uint64_t* hash_array = (uint64_t* )calloc(hash_n + 1, 8);
	if(hash_array == NULL)
	{
		printf("fail to allocate memory\n");
		fflush(stdout);
		exit(1);
	}
#else
    uint32_t* pos_array = (uint32_t* )calloc(line_tol, 4);
	if(pos_array == NULL)
	{
		printf("fail to allocate memory pos_array\n");
		exit(1);
	}
	uint32_t* kmer_point = (uint32_t* )calloc(line_tol + 1, 4);
	if(kmer_point == NULL)
	{
		printf("fail to allocate memory kmer_point\n");
		exit(1);
	}
	
	
	uint32_t* hash_array = (uint32_t* )calloc(hash_n + 1, 4);
	if(hash_array == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}
#endif
    //actually smaller than this amount
    uint32_t* kmer_array = (uint32_t* )calloc(line_tol, 4);
	if(kmer_array == NULL)
	{
		printf("fail to allocate memory\n");
		fflush(stdout);
		exit(1);
	}

    uint8_t* edge_array = (uint8_t* )calloc(line_tol, 1);
	if(edge_array == NULL)
	{
		printf("fail to allocate memory edge_array\n");
		fflush(stdout);
		exit(1);
	}
	
	uint8_t* edge_flag = (uint8_t* )calloc(line_tol, 1);
	if(edge_flag == NULL)
	{
		printf("fail to allocate memory edge_flag\n");
		fflush(stdout);
		exit(1);
	}
	
    printf("total number of chars to allocate: %" PRIu64 "\n", (((uint64_t )line_tol) << 2));//36215262392
	fflush(stdout);
	
    uint64_t kmerf_cnt = 0;
    uint32_t hash_add = 0;

    uint8_t input = 0;
    uint8_t output = 0;
	
#ifdef	UNPIPATH_OFF_K20
	uint64_t tmp_pos64 = 0;
	uint64_t file_line_cnt_p = 0;
#else
    uint32_t file_line_cnt_p = 0;
#endif
	
    const uint32_t edge_bit = 1;
	uint32_t hash_first = 0;

    //end declaration
    for(i = 0; i < (1 << (f << 1)); i++)
    {
        bit32_kmer(i, f, kmer_i, file_div)
        strcpy(file_open, (const char* )filename_div);
        strcat(file_open, (const char* )file_div);
        strcat(file_open, suff);

        if(file_p[i] == NULL)
        {
            file_p[i] = fopen(file_open, "rb"); //all of div files are existing by default
            if(file_p[i] == NULL)
            {
				continue;
            }

            k_p_file = (k_p* )calloc(file_line_cnt[i], sizeof(k_p));
            fread(k_p_file, sizeof(k_p), file_line_cnt[i], file_p[i]);

            qsort(k_p_file, file_line_cnt[i], sizeof(k_p), compare);

            printf("has finished %u file sorting: %u\n", i, file_line_cnt[i]);
		
            //load into hash graph
            kmer_l_p = 0Xffffffffffffffff;
            kmer_f_p = 0Xffffffff;

            for(k_p_i = 0; k_p_i < file_line_cnt[i]; k_p_i++)
            {
                //add position
                //also can write to file directly
#ifdef	UNPIPATH_OFF_K20
				tmp_pos64 = (k_p_file[k_p_i].kp)[3];
				
				pos_array[k_p_i + file_line_cnt_p] = ((tmp_pos64 << 32) | ((k_p_file[k_p_i].kp)[0]));
				
#else
                pos_array[k_p_i + file_line_cnt_p] = (k_p_file[k_p_i].kp)[0];
#endif
                //get first and whole, also get input and output char
                bit32a_bit64((k_p_file[k_p_i].kp), 3, k_t - f, k_t - k, kmer_i, kmer_l, kmer_f, kmer_s, input, output)

#ifdef	DEBUG_NOT_FOUND
				//0 8258096 126 8258096
				printf("\nreadin: k_p_i %u %u %u %u %u\n", k_p_i, (uint16_t )((k_p_file[k_p_i].kp)[1]), ((k_p_file[k_p_i].kp)[2]), kmer_f, kmer_l);
#endif				
				bit32_kmer(kmer_l, k_t - f, k_i, char_l)

				if(kmer_l != kmer_l_p)
                {
                    bit32_kmer(kmer_s, k_t - k, k_i, char_s)
                    bit32_kmer(kmer_l, k_t - f, k_i, char_l)

                    //new kmer begins and add pos number to kmer
                    kmer_array[kmerf_cnt] = kmer_s;
#ifdef	DEBUG_NOT_FOUND
					printf("kmer: %u\n\n", kmer_s);
					if(kmerf_cnt > 0)
						printf("pre edge: %u\n", edge_array[kmerf_cnt - 1]);
#endif
                    kmer_point[kmerf_cnt] = k_p_i + file_line_cnt_p;

                    if(kmer_f != kmer_f_p)
                    {
                        //new hash address and add number to hash
                        hash_add = (i << ((k - f) << 1)) + kmer_f;
                        hash_array[hash_add] = kmerf_cnt;
#ifdef	DEBUG_NOT_FOUND
						printf("hash: %u %u\n\n\n", hash_add,kmerf_cnt );
#endif						
						if(kmerf_cnt == 0)	hash_first = hash_add;
                    }
                    ++kmerf_cnt;
                }
#ifdef	DEBUG_NOT_FOUND
				printf("edge: %u %u %u\n", kmerf_cnt, input, output);
#endif				
                //add kmer edge
                if(charToDna5[input] <= 3)
                    edge_array[kmerf_cnt - 1] |= (edge_bit << (7 - charToDna5[input]));
				
                if(charToDna5[output] <= 3)
                    edge_array[kmerf_cnt - 1] |= (edge_bit << (3 - charToDna5[output]));
				else	edge_flag[kmerf_cnt - 1] = 1;
				
                kmer_l_p = kmer_l;
                kmer_f_p = kmer_f;

                //
                //kmer_s_p = kmer_s;
            }

            printf("has finished %u file loading\n",i);
			fflush(stdout);
			
            //end loading

            //end output

            file_line_cnt_p += file_line_cnt[i];

            free(k_p_file);
        }
    }
    kmer_point[kmerf_cnt] = line_tol;
    hash_array[hash_n] = kmerf_cnt;

    for(i = 0; i < file_n; i++)
    {
        if(file_p[i] != NULL)
            fclose(file_p[i]);
    }
    free(file_p);

	
	//for debug from here///////////////
/*
	FILE* fp_kmer_point_w = fopen("/home/ghz/viruses_index_all/tmp_kmer_point", "wb");
    if(fp_kmer_point_w == NULL)
	{
		printf("cannot open the kmer point tmp file\n");
		//exit(1);
	}
	
	fwrite(kmer_point, 8, kmerf_cnt + 1, fp_kmer_point_w);//line_tol
	
	fclose(fp_kmer_point_w);

	fp_kmer_point_w = fopen("/home/ghz/viruses_index_all/tmp_kmer_pos", "wb");
    if(fp_kmer_point_w == NULL)
	{
		printf("cannot open the kmer pos tmp file\n");
		//exit(1);
	}
	
	fwrite(pos_array, 8, file_line_cnt_p, fp_kmer_point_w);//line_tol
	
	fclose(fp_kmer_point_w);
*/
	
	////////////////
	
    //traverse the graph and get supernode

    char ws[KT_LENGTH_MAX+2] = "";
    char wsn[KT_LENGTH_MAX+2] = "";
    char fs[KF_LENGTH_MAX+2] = "";
    char ss[KF_LENGTH_MAX+2] = "";

    uint8_t l_flag = 0;
    uint8_t n_flag = 0;
    uint8_t d_flag = 0;
    uint8_t fy_flag = 0;
    uint8_t ry_flag = 0;

    //unipath seq file
    fp_us = fopen(uniseq, "w");
    if(fp_us == NULL)
	{
		printf("cannot open the unipath seq file\n");
		fflush(stdout);
		exit(1);
	}

    //unipath offset binary file
    fp_usf_b = fopen(uniseqf_b, "wb");
    if(fp_usf_b == NULL)
	{
		printf("cannot open the unipath seq offset binary file\n");
		fflush(stdout);
		exit(1);
	}

	/*
    //unipath branch node file
    fp_ue = fopen(uniedge, "wb");
    if(fp_ue == NULL)
	{
		printf("cannot open the unipath branch node file\n");
		exit(1);
	}
	*/
	
	uint32_t ue_n = 0;
	uint32_t uni_offset_ori = 0;
    uint64_t kmer_ws = 0;
    uint64_t kmer_wsli = 0;
    uint64_t kmer_fs = 0;
    uint64_t kmer_ss = 0;
    uint64_t kmer_wsn = 0;
    uint32_t t_hash_i = 0;
    
    uint8_t node_i = 0;
    uint8_t edge_i = 0;
    uint8_t edge_out = 0;
    int64_t binary_r = 0;
    uint64_t one_bit64 = 1;
#ifdef	UNPIPATH_OFF_K20
	uint64_t hash_v_p = 0;
	uint64_t t_kmer_i = 0;
#else
    uint32_t hash_v_p = 0;
	uint32_t t_kmer_i = 0;
#endif
    uint32_t super_node_id = 0;
    uint8_t uni_edge = 0;

    //for statistic
    uint32_t fy_branch_cnt = 0;
    uint32_t ry_branch_cnt = 0;
    uint32_t x_branch_cnt = 0;
    uint32_t l_blunt_cnt = 0;
    uint32_t b_blunt_cnt = 0;
    uint32_t linear_cnt = 0;
    uint32_t e_blunt_cnt = 0;
    uint32_t es_blunt_cnt = 0;
    uint8_t li_cnt = 0;
    uint32_t kmer_intervali = 0;

    //fill the blank of hash array
    for(t_hash_i = hash_n; t_hash_i > 0; t_hash_i--)
    {
        if((hash_array[t_hash_i] == 0) && (t_hash_i > hash_first))
            hash_array[t_hash_i] = hash_v_p;
        else	hash_v_p = hash_array[t_hash_i];
    }

    printf("finish filling the hash blank\n");
	fflush(stdout);
	
	uint64_t usf_n = 1;

#ifdef UNPIPATH_OFF_K20
	uint64_t uni_max_l = 0;
    uint64_t uni_l = 0;
	uint64_t uni_tmp = 0;
	uint64_t unioff = 0;
	uint64_t uni_kmer_off = 0;
    fwrite(&unioff, 8, 1, fp_usf_b);
	uint64_t* uni_offset = (uint64_t* )calloc(kmerf_cnt, 8);
#else
	uint32_t uni_max_l = 0;
    uint32_t uni_l = 0;
	uint32_t uni_tmp = 0;
	uint32_t unioff = 0;
	uint32_t uni_kmer_off = 0;
	fwrite(&unioff, 4, 1, fp_usf_b);
	uint32_t* uni_offset = (uint32_t* )calloc(kmerf_cnt, 4);
#endif

	if(uni_offset == NULL)
	{
		printf("Failed to allocate the unipath offset array\n");
		fflush(stdout);
		exit(1);
	}

	printf("begin creating unipath: %"PRId64"\n", hash_n);
	fflush(stdout);
	
	uint32_t uni_node_cnt = 0;

    for(t_hash_i = 0; t_hash_i < hash_n; t_hash_i++)//268435456
    {
        for(t_kmer_i = hash_array[t_hash_i]; t_kmer_i < hash_array[t_hash_i + 1]; t_kmer_i++)
        {
            //get the whole kmer
            //kmer_array[t_kmer_i]: get last kmer
            kmer_ws = (((uint64_t )t_hash_i) << ((k_t - k) << 1)) + kmer_array[t_kmer_i];

            bit32_kmer(kmer_ws, k_t, kmer_i, ws)

            node_i = node_indentity(edge_array[t_kmer_i], edge_flag[t_kmer_i]);

            if(node_i == 1)	linear_cnt++;

            if(node_i != 1)
            {
                //statistic
                if(node_i == 2)	fy_branch_cnt++;
                if(node_i == 3)	ry_branch_cnt++;
                if(node_i == 4)	x_branch_cnt++;
                if(node_i == 5)	l_blunt_cnt++;
                if(node_i == 6)	b_blunt_cnt++;
                if(node_i == 7)	e_blunt_cnt++;
                if(node_i == 8)	es_blunt_cnt++;

                //unipath edge

                if((node_i != 2) && (node_i != 7))
                {
                    uni_edge |= ((edge_array[t_kmer_i] >> 4) << 4);
                }

                //trace to record super node ID
                uni_l = 0;
                fy_flag = 0;
                ry_flag = 0;
                if(node_i == 2)	fy_flag = 1;
                if((node_i == 3) || (node_i == 5))
                {
                    if(strlen(ws) != k_t)
					{
						printf("kmer length error\n");
						fflush(stdout);
						exit(1);
					}

                    uni_offset[t_kmer_i] = uni_kmer_off;

					++uni_kmer_off;

                    unioff += strlen(ws);
                    uni_l += strlen(ws);

                    ry_flag = 1;

                    fprintf(fp_us, "\n%s",ws);
                }

                d_flag = 0;
                if((node_i == 4) || (node_i == 6))
                {
                    if(strlen(ws) != k_t)
					{
						printf("kmer length error\n");
						fflush(stdout);
						exit(1);
					}

                    uni_offset[t_kmer_i] = uni_kmer_off;

					uni_kmer_off += k_t;

                    unioff += strlen(ws);
                    uni_l += strlen(ws);

                    d_flag = 1;

                    ++super_node_id;

#ifdef UNPIPATH_OFF_K20
                    fwrite(&unioff, 8, 1, fp_usf_b);
#else
					fwrite(&unioff, 4, 1, fp_usf_b);
#endif

					usf_n++;

                    if(uni_l > uni_max_l)	uni_max_l = uni_l;

                    fprintf(fp_us, "\n%s",ws);

                    uni_edge |= (edge_array[t_kmer_i] & 0Xf);
                    //fwrite(&uni_edge, 1, 1, fp_ue);
					ue_n++;

                    uni_edge = 0;
                }

                if(node_i == 8)
                {
                    if(strlen(ws) != k_t)
					{
						printf("kmer length error\n");
						fflush(stdout);
						exit(1);
					}

                    uni_offset[t_kmer_i] = uni_kmer_off;

					uni_kmer_off += k_t;

                    unioff += strlen(ws);
                    uni_l += strlen(ws);

                    ++super_node_id;
#ifdef UNPIPATH_OFF_K20
                    fwrite(&unioff, 8, 1, fp_usf_b);
#else
					fwrite(&unioff, 4, 1, fp_usf_b);
#endif
					usf_n++;

                    if(uni_l > uni_max_l)	uni_max_l = uni_l;

                    fprintf(fp_us, "\n%s",ws);

                    //fwrite(&uni_edge, 1, 1, fp_ue);
					ue_n++;

                    uni_edge = 0;

                    continue;
                }

                uni_tmp = uni_l;
                for(edge_i = 0; edge_i < 4; edge_i++) //A,C,G,T
                {
                    edge_out = (uint8_t )((edge_array[t_kmer_i] >> (3 - edge_i)) & 0x1);
                    if(edge_out == 0)
                    {
                        continue;
                    }

                    uni_l = uni_tmp;

                    next_kmer(kmer_ws, kmer_wsn, k_t, edge_i)

                    bit32_kmer(kmer_wsn, k_t, kmer_i, wsn)

                    uint64_div(kmer_wsn, kmer_fs, kmer_ss, k_t - k)

                    bit32_kmer(kmer_fs, k, kmer_i, fs)
                    bit32_kmer(kmer_ss, k_t-k, kmer_i, ss)
#ifdef	UNPIPATH_OFF_K20
					binary_r = binsearch_offset_index64(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#else
                    binary_r = binsearch_offset_index(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#endif

                    if(binary_r == -1)
                    {
						printf("kmer search error\n");
#ifdef	DEBUG_64BIT
						//504 2240 3090 3088
						printf("binsearch_offset_index1, not found %u %u %"PRId64" %"PRId64"\n", kmer_fs, kmer_ss, hash_array[kmer_fs + 1], hash_array[kmer_fs]);
						uint64_t tmp_j = 0;
						printf("\n");
						for(tmp_j = hash_array[kmer_fs]; tmp_j < hash_array[kmer_fs + 1]; tmp_j++)
							printf("%u\n", kmer_array[tmp_j]);
#endif						
						fflush(stdout);
						exit(1);
                    }

                    kmer_wsli = kmer_wsn;

                    //record seq
                    node_i = node_indentity(edge_array[binary_r], edge_flag[binary_r]);
                    n_flag = 0;	//2 4 6

                    if(((d_flag == 1) || (fy_flag == 1)) && ((node_i == 1) || (node_i == 2) || (node_i == 7)))
                    {
                        if(strlen(wsn) != k_t)
						{
							printf("kmer length error\n");
#ifdef	DEBUG_64BIT
							printf("wrong kmer length %u %u\n", strlen(wsn), k_t);
#endif							
							fflush(stdout);
							exit(1);
						}

                        uni_offset[binary_r] = uni_kmer_off;

						++uni_kmer_off;

                        unioff += strlen(wsn);
                        uni_l += strlen(wsn);

                        fprintf(fp_us, "\n%s",wsn);

                        //unipath input edge
                        uni_edge |= ((edge_array[binary_r] >> 4) << 4);

                        n_flag = 1;
                    }

                    //s_flag = 0; //3 5
                    if((ry_flag == 1) && ((node_i == 1) || (node_i == 2) || (node_i == 7)))
                    {
                        fprintf(fp_us, "%c", Dna5Tochar[edge_i]);

                        uni_offset[binary_r] = uni_kmer_off;

						++uni_kmer_off;

                        ++unioff;
                        ++uni_l;

                        //s_flag = 1;
                    }

                    l_flag = 0;
                    li_cnt = 0;
                    while(node_indentity(edge_array[binary_r], edge_flag[binary_r]) == 1)
                    {
                        ++li_cnt;
                        edge_out = edgeout_node(edge_array[binary_r]);
                        if(edge_out == 4)
                        {
                            printf("edge identification error \n");
							fflush(stdout);
                            exit(1);
                        }

                        ++unioff;
                        ++uni_l;

                        fprintf(fp_us, "%c", Dna5Tochar[edge_out]);

                        l_flag = 1;

                        next_kmer(kmer_wsli, kmer_wsn, k_t, edge_out)

                        uint64_div(kmer_wsn, kmer_fs, kmer_ss, k_t - k)

						bit32_kmer(kmer_fs, k, kmer_i, fs)
						bit32_kmer(kmer_ss, k_t-k, kmer_i, ss)
#ifdef	UNPIPATH_OFF_K20
						binary_r = binsearch_offset_index64(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#else
                        binary_r = binsearch_offset_index(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#endif						
						if(binary_r == -1)
						{
							printf("kmer search error\n");
#ifdef	DEBUG_64BIT
							printf("binsearch_offset_index2, not found %u %u %"PRId64" %"PRId64"\n", kmer_fs, kmer_ss, hash_array[kmer_fs + 1], hash_array[kmer_fs]);
#endif							
							fflush(stdout);
							exit(1);
						}
						uni_offset_ori = uni_offset[binary_r];
						uni_offset[binary_r] = uni_kmer_off;

						++uni_kmer_off;

                        kmer_wsli = kmer_wsn;
                    }

                    node_i = node_indentity(edge_array[binary_r], edge_flag[binary_r]);

                    if(((node_i == 3) || (node_i == 4) || (node_i == 8)) && (l_flag == 1))
                    {
                        fseek(fp_us, -1, SEEK_CUR);
                        --unioff;
                        --uni_l;

						uni_offset[binary_r] = uni_offset_ori;

						--uni_kmer_off;
                    }

                    if(n_flag || l_flag || ry_flag)
                    {
                        ++super_node_id;

#ifdef UNPIPATH_OFF_K20
                        fwrite(&unioff, 8, 1, fp_usf_b);
#else
						fwrite(&unioff, 4, 1, fp_usf_b);
#endif
						usf_n++;

                        if(uni_l > uni_max_l)	uni_max_l = uni_l;

                        //unipath edge
                        if((node_i == 2) || (node_i == 7))
                            uni_edge |= (edge_array[binary_r] & 0Xf);
                        else{
							if(l_flag)	uni_edge |= (1 << (3 - edge_out));
							else	uni_edge |= (edge_array[t_kmer_i] & 0Xf);
						}

                        //fwrite(&uni_edge, 1, 1, fp_ue);
						ue_n++;

                        uni_edge = 0;

						uni_kmer_off += (k_t - 1);
                    }
                }
            }

            uni_node_cnt++;
        }
		
	}
	
	if(edge_array)	free(edge_array);
	if(edge_flag)	free(edge_flag);
	
    fprintf(fp_us,"\n");

    printf("number of supernode: %u\n", super_node_id);
    printf("number of linear node: %u \n", linear_cnt);
    printf("number of forward Y branch node: %u \n", fy_branch_cnt);
    printf("number of reverse Y branch node: %u \n", ry_branch_cnt);
    printf("number of X branch node: %u \n", x_branch_cnt);
    printf("number of linear blunt node: %u \n", l_blunt_cnt);
    printf("number of branch blunt node: %u \n", b_blunt_cnt);
    printf("number of end blunt node: %u \n", e_blunt_cnt);
    printf("number of ends blunt node: %u \n", es_blunt_cnt);
	
	fflush(stdout);
	
#ifdef	UNPIPATH_OFF_K20
	printf("the max length of supernode: %"PRId64"\n", uni_max_l);
#else
    printf("the max length of supernode: %u \n", uni_max_l);
#endif
	
    //end traverse
	//fclose(fp_ue);
    fclose(fp_us);
    fclose(fp_usf_b);

	printf("begin wrtting hash and kmer array and offset on unipath into file %"PRId64"\n", kmerf_cnt);
	fflush(stdout);
	
	fp_hash = fopen(unihash_g,"wb");
	if(fp_hash == NULL)
	{
		printf("wrong file route of hash to write\n");
		fflush(stdout);
		exit(1);
	}
	
	fp_kmer = fopen(unikmer_g,"wb");
	if(fp_kmer == NULL)
	{
		printf("wrong file route of kmer to write\n");
		fflush(stdout);
		exit(1);
	}
	
	fp_off = fopen(unioff_g,"wb");
	if(fp_off == NULL)
	{
		printf("wrong file route of off_set to write\n");
		fflush(stdout);
		exit(1);
	}

	fwrite(kmer_array, 4, kmerf_cnt, fp_kmer);

#ifdef UNPIPATH_OFF_K20
	fwrite(hash_array, 8, hash_n+1, fp_hash);
	fwrite(uni_offset, 8, kmerf_cnt, fp_off);
#else
	fwrite(hash_array, 4, hash_n+1, fp_hash);
	fwrite(uni_offset, 4, kmerf_cnt, fp_off);
#endif
	if(uni_offset)	free(uni_offset);

	fclose(fp_hash);
	fclose(fp_kmer);
	fclose(fp_off);
#else
    //create the unipath position file
    printf("Begin creating unipath position and its point file\n");
	fflush(stdout);

	//for debug from here //////////////////////////////////
	uint64_t one_bit64 = 1;
	uint64_t kmer_fs = 0;
    uint64_t kmer_ss = 0;
	uint64_t kmerf_cnt = 0;
	uint64_t usf_n = 0;
	uint64_t ue_n = 0;
	int64_t binary_r = 0;
	
	char fs[KF_LENGTH_MAX+2] = "";
    char ss[KF_LENGTH_MAX+2] = "";
	
	//FILE* fp_hash = NULL;
	//FILE* fp_kmer = NULL;
	FILE* fp_kmer_point = NULL;
	
	uint64_t a_size = 0;
	int64_t result_hash_g = 0;
	int64_t result_kmer_g = 0;
	uint64_t hash_n = 268435456;
	uint64_t kmer_n = 0;
	
	hash_n = ((hash_n + 1) << 3);
	
	printf("Load unipath hash\n");
	
	fp_hash = fopen("/home/ghz/viruses_index_all/unipath_g.hash", "rb");
    if(fp_hash == NULL)
    {
        fputs ("File error opening the graph hash file\n",stderr);
        exit (1);
    }
	
	a_size = (hash_n >> 3);
    uint64_t* hash_array = (uint64_t* ) malloc (hash_n);
    if (hash_array == NULL)
    {
        fputs("Memory error hash_array",stderr);
        exit(2);
    }

    // copy the file into the buffer:
    result_hash_g = fread (hash_array, 8, a_size, fp_hash);
	
    if (result_hash_g != a_size)
    {
        fputs("Reading error",stderr);
        exit(3);
    }
	
	fclose(fp_hash);

    //read input graph kmer file
    printf("Load unipath kmer\n");

    fp_kmer = fopen("/home/ghz/viruses_index_all/unipath_g.kmer", "rb");
    if(fp_kmer == NULL)
    {
        fputs ("File error opening the graph hash file\n",stderr);
        exit (1);
    }
	
	fseek(fp_kmer, 0, SEEK_END);// non-portable
    kmer_n = ftell(fp_kmer);
	rewind(fp_kmer);
	
    a_size = (kmer_n >> 2);
    
    //a_size = kmer_num;

    uint32_t* kmer_array = (uint32_t* ) malloc (kmer_n);
    if (kmer_array == NULL)
    {
        fputs ("Memory error kmer_array", stderr);
        exit (2);
    }

    // copy the file into the buffer:
    result_kmer_g = fread (kmer_array, 4, a_size, fp_kmer);
    if (result_kmer_g != a_size)
    {
        fputs ("Reading error kmer_array", stderr);
        exit (3);
    }

    fclose(fp_kmer);
	
	
	printf("Load unipath kmer point\n");

    fp_kmer_point = fopen("/home/ghz/viruses_index_all/tmp_kmer_point", "rb");
    if(fp_kmer_point == NULL)
    {
        fputs ("File error opening the graph hash file\n",stderr);
        exit (1);
    }
	
	fseek(fp_kmer_point, 0, SEEK_END);// non-portable
    kmer_n = ftell(fp_kmer_point);
	rewind(fp_kmer_point);
	
    a_size = (kmer_n >> 3);
    
    //a_size = kmer_num;

    uint64_t* kmer_point = (uint64_t* ) malloc (kmer_n);
    if (kmer_point == NULL)
    {
        fputs ("Memory error kmer_point",stderr);
        exit (2);
    }

    // copy the file into the buffer:
    result_kmer_g = fread (kmer_point, 8, a_size, fp_kmer_point);
    if (result_kmer_g != a_size)
    {
        fputs ("Reading error kmer_point",stderr);
        exit (3);
    }

    fclose(fp_kmer_point);

	printf("Load unipath kmer pos\n");

    fp_kmer_point = fopen("/home/ghz/viruses_index_all/tmp_kmer_pos", "rb");
    if(fp_kmer_point == NULL)
    {
        fputs ("File error opening the graph hash file\n",stderr);
        exit (1);
    }
	
	fseek(fp_kmer_point, 0, SEEK_END);// non-portable
    kmer_n = ftell(fp_kmer_point);
	rewind(fp_kmer_point);
	
    a_size = (kmer_n >> 3);
    
    //a_size = kmer_num;

    uint64_t* pos_array = (uint64_t* ) malloc (kmer_n);
    if (pos_array == NULL)
    {
        fputs ("Memory error pos_array",stderr);
        exit (2);
    }

    // copy the file into the buffer:
    result_kmer_g = fread (pos_array, 8, a_size, fp_kmer_point);
    if (result_kmer_g != a_size)
    {
        fputs ("Reading error pos_array",stderr);
        exit (3);
    }

    fclose(fp_kmer_point);

	int uni_max_l = 600000;
	uint32_t kmer_i = 0;
	////////////////////////////////
#endif	

    fp_us = fopen(uniseq, "r");
    if(fp_us == NULL)
	{
		printf("cannot read the unipath seq file\n");
		fflush(stdout);
		exit(1);
	}

    //unipath statistics
	/*
    fp_sta = fopen(unista, "wb");
    if(fp_sta == NULL)
	{
		printf("cannot open the unipath statistical file\n");
		fflush(stdout);
		exit(1);
	}
	*/
    //unipath position file
    fp_up = fopen(unipos, "wb");
    if(fp_up == NULL)
	{
		printf("cannot open the unipath position file\n");
		fflush(stdout);
		exit(1);
	}

    //unipath position point file
    fp_upp = fopen(uniposp, "wb");
    if(fp_upp == NULL)
	{
		printf("cannot open the unipath position point file\n");
		fflush(stdout);
		exit(1);
	}

	uint64_t up_n = 0;
	uint64_t upp_n = 0;

	uint8_t uni_de = 0;
    uint32_t pos_l_i = 0;

    char* uni_in_seq = (char* )malloc(uni_max_l + 2);
	if(uni_in_seq == NULL)
	{
		printf("fail to allocate memory for unipath seq input array\n");
		fflush(stdout);
		exit(1);
	}
	
    char uni_fkmer[KT_LENGTH_MAX];

    uint32_t uni_in_i = 0;
    uint32_t line_l = 0;
    uint32_t first_i = 0;
	uint32_t line_cntu = 0;
    uint64_t cnt_cl_cnt = 0;
    //the number of kmers in unipath seq hash kmer array
    uint64_t us_h_n = 0;
	
#ifdef	UNPIPATH_OFF_K20
	fwrite(&first_cnt_cl, 8, 1, fp_upp);
	uint64_t* first = NULL;
    uint64_t* first_pos = NULL;
	uint64_t first_cnt = 0;
	uint64_t r_cnt = 0;
	uint64_t uni_pos_i = 0;
	uint64_t pos_first = 0;
	uint64_t* r = NULL;
#else
    fwrite(&first_cnt_cl, 4, 1, fp_upp);
	uint32_t* first = NULL;
    uint32_t* first_pos = NULL;
	uint32_t first_cnt = 0;
	uint32_t r_cnt = 0;
	uint32_t uni_pos_i = 0;
	uint32_t pos_first = 0;
	uint32_t* r = NULL;
#endif
	int64_t second = 0;
	uint64_t uni_kmer = 0;
	
	upp_n++;

    //transfer to binary seq file
#ifdef UNI_SEQ64
	uint64_t uni_arr[UNI_SEQ_WRI_ARR];
#else
    uint8_t uni_arr[UNI_SEQ_WRI_ARR];
#endif
    uint32_t s_b_tmp = 0;
    uint64_t a_cnt = 0;
    uint32_t seq_i = 0;
    uint32_t w_offset = 0;
	

    fp_us_b = fopen(uniseq_b, "wb");
    if(fp_us_b == NULL)
	{
		printf("cannot open the unipath seq file\n");
		fflush(stdout);
		exit(1);
	}		
#ifdef UNI_SEQ64
	memset(uni_arr, 0, (UNI_SEQ_WRI_ARR << 3));
#else
    memset(uni_arr, 0, (UNI_SEQ_WRI_ARR));
#endif
    //repetitive read the last line

    fgets(uni_in_seq, uni_max_l + 2, fp_us);
    while((!feof(fp_us)) && (fgets(uni_in_seq, uni_max_l + 2, fp_us) > 0))
    {
        line_cntu++;
        line_l = strlen(uni_in_seq) - 1;

        if(!line_l)	continue;

        if(line_l < k_t)
        {
            printf("Error of unipath seq file: %u %u %s\n", line_cntu, line_l, uni_in_seq);
			fflush(stdout);
            exit(1);
        }
        else
        {
            us_h_n += (line_l - k_t + 1);

            //transfer to binary seq file
            for(seq_i = 0; seq_i < line_l; seq_i++)
            {
#ifdef UNI_SEQ64
				binary64_array64_f(uni_arr, w_offset, s_b_tmp, charToDna5[(uint8_t )uni_in_seq[seq_i]], UNI_SEQ_WRI_ARR, a_cnt, fp_us_b)
#else
                binary64_array(uni_arr, w_offset, s_b_tmp, charToDna5[(uint8_t )uni_in_seq[seq_i]], UNI_SEQ_WRI_ARR, a_cnt, fp_us_b)
#endif
				w_offset++;
            }

            strncpy(uni_fkmer, uni_in_seq, k_t);

            kmer_bit64(uni_kmer, k_t, kmer_i, uni_fkmer)

            uint64_div(uni_kmer, kmer_fs, kmer_ss, k_t - k)

            //debug
            bit32_kmer(kmer_fs, k, kmer_i, fs)
            bit32_kmer(kmer_ss, k_t-k, kmer_i, ss)
#ifdef	UNPIPATH_OFF_K20
			binary_r = binsearch_offset_index64(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#else
            binary_r = binsearch_offset_index(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#endif
			if(binary_r == -1)
			{
				printf("kmer search error\n");
#ifdef	DEBUG_64BIT
				printf("binsearch_offset_index, not found %u %u %"PRId64" %"PRId64"\n", kmer_fs, kmer_ss, hash_array[kmer_fs + 1], hash_array[kmer_fs]);
#endif				
				fflush(stdout);
				exit(1);
			}
            first_cnt = kmer_point[binary_r + 1] - kmer_point[binary_r];
			
#ifdef	UNPIPATH_OFF_K20
			first = (uint64_t* )calloc(first_cnt, 8);
			first_pos = (uint64_t* )calloc(first_cnt, 8);
#else
			first = (uint32_t* )calloc(first_cnt, 4);
			first_pos = (uint32_t* )calloc(first_cnt, 4);
#endif
          

            first_i = 0;
            for(uni_pos_i = kmer_point[binary_r]; uni_pos_i < kmer_point[binary_r + 1]; uni_pos_i++)
            {
                first[first_i] = pos_array[uni_pos_i];
				
		
                first_pos[first_i] = pos_array[uni_pos_i];

                first_i++;
            }
	
            uni_de = 0;

            for(uni_in_i = 1; uni_in_i < line_l - k_t + 1; uni_in_i++)
            {
                strncpy(uni_fkmer, uni_in_seq + uni_in_i, k_t);
                kmer_bit64(uni_kmer, k_t, kmer_i, uni_fkmer)

                uint64_div(uni_kmer, kmer_fs, kmer_ss, k_t - k)
				
#ifdef	UNPIPATH_OFF_K20
				second = binsearch_offset_index64(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#else
                second = binsearch_offset_index(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#endif
				
                if(second != -1)
                {
#ifdef	UNPIPATH_OFF_K20
					r_cnt = uni_pos_combine64(kmer_point, pos_array, first, first_cnt, second, &r);
#else
                    r_cnt = uni_pos_combine(kmer_point, pos_array, first, first_cnt, second, &r);
#endif
                }
                else
				{
					printf("kmer search error\n");
#ifdef	DEBUG_64BIT
					printf("error -1 binsearch_offset_index, not found %u %u %"PRId64" %"PRId64"\n", kmer_fs, kmer_ss, hash_array[kmer_fs + 1], hash_array[kmer_fs]);
#endif					
					fflush(stdout);
					exit(1);
				}

                if(!r_cnt)
                {		
                    if(!uni_de)
                    {
                        for(uni_pos_i = 0; uni_pos_i < first_cnt; uni_pos_i++)
                        {
                            for(pos_l_i = 0; pos_l_i < chr_file_n; pos_l_i++)
                            {
                                if(first_pos[uni_pos_i] < chr_end_n[pos_l_i])
                                    break;
                            }
						}

                        uni_de = 1;
                    }
                }

                first_cnt = r_cnt;
                first = r;

				if(!r_cnt)	break;
					
            }

            for(uni_pos_i = 0; uni_pos_i < first_cnt; uni_pos_i++)
            {
				pos_first = first[uni_pos_i] + 1 - uni_in_i;

				if(pos_first <= START_POS_REF)
				{
					printf("combine pos error\n");
					fflush(stdout);
					exit(1);
				}
#ifdef	UNPIPATH_OFF_K20
				fwrite(&pos_first, 8, 1, fp_up);
#else
				fwrite(&pos_first, 4, 1, fp_up);
#endif
				up_n++;
            }

            first_cnt_cl += first_cnt;
#ifdef	UNPIPATH_OFF_K20
			fwrite(&first_cnt_cl, 8, 1, fp_upp);
#else
            fwrite(&first_cnt_cl, 4, 1, fp_upp);
#endif
			upp_n++;

            cnt_cl_cnt++;

            free(first);

            free(first_pos);
        }
    }

#ifdef UNI_SEQ64
    fwrite(uni_arr, 8, UNI_SEQ_WRI_ARR, fp_us_b);
#else
	fwrite(uni_arr, 1, UNI_SEQ_WRI_ARR, fp_us_b);
#endif

	//free hash graph memory
	if(pos_array)	free(pos_array);
    if(kmer_array)	free(kmer_array);
    if(kmer_point)	free(kmer_point);
    if(hash_array)	free(hash_array);
	if(uni_in_seq)	free(uni_in_seq);
    //for index file close
    fclose(fp_us_b);
    fclose(fp_up);
    fclose(fp_upp);
    fclose(fp_us);
	
	printf("begin writing info about size of index files\n");
	fflush(stdout);
	
	fp_num = fopen(unisize, "wb");
	if(fp_num == NULL)	printf("error of creating index size file\n");

#ifdef UNI_SEQ64
	sta_num_write = ((a_cnt + 1) * (UNI_SEQ_WRI_ARR << 3));
#else
	sta_num_write = ((a_cnt + 1) * (UNI_SEQ_WRI_ARR));
#endif

	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef UNPIPATH_OFF_K20
	sta_num_write = (usf_n << 3);
#else
	sta_num_write = (usf_n << 2);
#endif

	fwrite(&sta_num_write, 8, 1, fp_num);

	sta_num_write = ue_n;
	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef	UNPIPATH_OFF_K20
	sta_num_write = (up_n << 3);
#else
	sta_num_write = (up_n << 2);
#endif
	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef	UNPIPATH_OFF_K20
	sta_num_write = (upp_n << 3);
#else	
	sta_num_write = (upp_n << 2);
#endif
	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef	UNPIPATH_OFF_K20
	sta_num_write = ((hash_n + 1) << 3);
#else
	sta_num_write = ((hash_n + 1) << 2);
#endif
	fwrite(&sta_num_write, 8, 1, fp_num);

	sta_num_write = (kmerf_cnt << 2);
	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef UNPIPATH_OFF_K20
	sta_num_write = (kmerf_cnt << 3);
#else
	sta_num_write = (kmerf_cnt << 2);
#endif
	fwrite(&sta_num_write, 8, 1, fp_num);
	
	sta_num_write = 0;
	fwrite(&sta_num_write, 8, 1, fp_num);
	
	sta_num_write = (((ref_seq_n >> 5) + 1) << 3);
	fwrite(&sta_num_write, 8, 1, fp_num);
	fflush(fp_num);
	
	fflush(fp_num);
	
	fclose(fp_num);

		
#ifdef	HANDLE_DIR
	//delete the div tmp folder
	printf("deleting div temporary index file folder\n");
	
	char rm_route[ROUTE_LENGTH_MAX];
	DIR* directory_pointer = NULL;

	if((directory_pointer = opendir(filename_div)) != NULL)
	{
		memset(rm_route, 0, ROUTE_LENGTH_MAX);
		strcpy(rm_route, sys_c_rm);
		strcat(rm_route, filename_div);
        
		system(rm_route);
	}
	
	memset(rm_route, 0, ROUTE_LENGTH_MAX);
	strcpy(rm_route, sys_rm);
	strcat(rm_route, uniseq);
	system(rm_route);
#endif
	
    printf("finish index building\n");

    return w_offset;
}

//create pos-unipath table
void build_pos_unipath()
{
	/*
    printf("begin creating position unipath file\n");
#ifdef	UNPIPATH_OFF_K20
	printf("number of unipath positions: %"PRId64"\n",first_cnt_cl);
#else
    printf("number of unipath positions: %u\n",first_cnt_cl);
#endif
    //unipath position file
    fp_up = fopen(unipos, "rb");
    if(fp_up == NULL)	printf("cannot read the unipath position file\n");

    //unipath position point file
    fp_upp = fopen(uniposp, "rb");
    if(fp_upp == NULL)	printf("cannot read the unipath position point file\n");

    //position unipath file
    fp_pu = fopen(unipu, "wb");
    if(fp_pu == NULL)	printf("cannot open the  position unipath file\n");

    uint32_t p_cnt = 0;
    uint32_t p_cnt_i = 0;
    uint32_t p_cnt_pre = 0;
    uint32_t p_i = 0;
    
    p_u* p_u_file = NULL;

    p_u_file = (p_u* )calloc(first_cnt_cl, sizeof(p_u));

    while((!feof(fp_upp)) && (fread(&p_cnt, 4, 1, fp_upp) > 0))
    {
        for(p_i = 0; p_i < p_cnt - p_cnt_pre; p_i++)
        {
            fread(&(p_u_file[p_u_c].pos), 4, 1, fp_up);
            p_u_file[p_u_c].uniid = p_cnt_i;

            p_u_c++;
        }
        p_cnt_pre = p_cnt;
        p_cnt_i++;
    }

    qsort(p_u_file, p_u_c, sizeof(p_u), compare_pu);

    for(p_i = 0; p_i < p_u_c; p_i++)
    {
        fwrite(&(p_u_file[p_i].pos),4, 1, fp_pu);
        fwrite(&(p_u_file[p_i].uniid),4, 1, fp_pu);
    }
	
	free(p_u_file);
    fclose(fp_up);
    fclose(fp_upp);
    fclose(fp_pu);
	
	
	uint64_t p_u_c = 0;
	sta_num_write = (p_u_c << 3);
	fwrite(&sta_num_write, 8, 1, fp_num);
	*/
}

#ifdef	UNPIPATH_OFF_K20
uint64_t uni_pos_combine64(uint64_t kmer_point[], uint64_t pos_array[], uint64_t first[], uint64_t first_cnt, int64_t second, uint64_t** r)
{
    uint64_t pos_i = 0;
    uint64_t r_cnt = 0;
    uint64_t alloc_cnt = 0;
    uint64_t r_i = 0;
    int64_t b_r = 0;

    alloc_cnt = kmer_point[second + 1] - kmer_point[second];

    r_cnt = (first_cnt > alloc_cnt ? first_cnt : alloc_cnt);
    (*r) = (uint64_t* )calloc(r_cnt, 8);

    //can change to less time usage
    for(pos_i = 0; pos_i < first_cnt; pos_i++)
    {
        if((b_r = binsearch_offset_index_com64(first[pos_i] + 1, pos_array, alloc_cnt, kmer_point[second])) != -1)
        {
            (*r)[r_i] = first[pos_i] + 1;
            r_i++;
        }
    }
	
#ifdef	DEBUG_COMBINE
	if(r_i == 0)
	{
		printf("new\nposes:\n");
		for(pos_i = 0; pos_i < first_cnt; pos_i++)
			printf("%"PRId64"\n", first[pos_i]);
		printf("\n");
		
		printf("poses:\n");
		uint64_t tmp_se = 0;
		for(tmp_se = kmer_point[second]; tmp_se < kmer_point[second + 1]; tmp_se++)
			printf("%"PRId64"\n", pos_array[tmp_se]);
	}
#endif
		
    free(first);
	
	
    return r_i;
}
#else
uint32_t uni_pos_combine(uint32_t kmer_point[], uint32_t pos_array[], uint32_t first[], uint32_t first_cnt, uint32_t second, uint32_t** r)
{
    uint32_t pos_i = 0;
    uint32_t r_cnt = 0;
    uint32_t alloc_cnt = 0;
    uint32_t r_i = 0;
    int64_t b_r = 0;

    alloc_cnt = kmer_point[second + 1] - kmer_point[second];
    r_cnt = (first_cnt > alloc_cnt ? first_cnt : alloc_cnt);
    (*r) = (uint32_t* )calloc(r_cnt, 4);

    //can change to less time usage
    for(pos_i = 0; pos_i < first_cnt; pos_i++)
    {
        if((b_r = binsearch_offset_index(first[pos_i] + 1, pos_array, alloc_cnt, kmer_point[second])) != -1)
        {
            (*r)[r_i] = first[pos_i] + 1;
            r_i++;
        }
    }

    free(first);

    return r_i;
}
#endif
//below is some operations on graph

//get the number of 1 in one value
uint8_t one_number_uint(uint8_t n)
{
    n=(n & 0x55555555)+((n>>1) & 0x55555555);
    n=(n & 0x33333333)+((n>>2) & 0x33333333);
    n=(n & 0x0F0F0F0F)+((n>>4) & 0x0F0F0F0F);
    n=(n & 0x00FF00FF)+((n>>8) & 0x00FF00FF);
    n=(n & 0x0000FFFF)+((n>>16) & 0x0000FFFF);

    return n;
}

//identify branching node
uint8_t node_indentity(uint8_t edge, uint8_t edge_flag)
{
	uint8_t node_i = 0;
	
	if(edge_flag)	node_i =  2;
	else{
		uint8_t edge_out = 0;
		uint8_t edge_in = 0;
		edge_out = (edge & 0xf);
		edge_in = (edge >> 4);
    
		if((one_number_uint(edge_in) == 1) && (one_number_uint(edge_out) == 1))	node_i =  1; //linear
		if((one_number_uint(edge_in) == 1) && (one_number_uint(edge_out) > 1))	node_i =  2; //forward Y branch
		if((one_number_uint(edge_in) > 1) && (one_number_uint(edge_out) == 1))	node_i =  3; //reverse Y branch
		if((one_number_uint(edge_in) > 1) && (one_number_uint(edge_out) > 1))	node_i =  4; //X branch
		if((one_number_uint(edge_in) == 0) && (one_number_uint(edge_out) == 1))	node_i =  5; //linear blunt
		if((one_number_uint(edge_in) == 0) && (one_number_uint(edge_out) > 1))	node_i =  6; //branch blunt
		if((one_number_uint(edge_in) == 1) && (one_number_uint(edge_out) == 0))	node_i =  7; //end blunt
		if((one_number_uint(edge_in) > 1) && (one_number_uint(edge_out) == 0))	node_i =  8; //ends blunt
	}
    
    return node_i;
}

//get the unique char of out edge
uint8_t edgeout_node(uint8_t edge)
{
    uint8_t outedge = 4;
    uint8_t edge_out = 0;
    edge_out = (edge & 0xf);

    if(edge_out == 1)	outedge = 3;
    if(edge_out == 2)	outedge = 2;
    if(edge_out == 4)	outedge = 1;
    if(edge_out == 8)	outedge = 0;

    return outedge;
}

//binary search offset
int64_t binsearch_offset_index_com64(uint64_t x, uint64_t v[], uint64_t n, uint64_t offset)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
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
            return mid + offset;
        }
    }

    return -1;
}
int64_t binsearch_offset_index64(uint32_t x, uint32_t v[], uint64_t n, uint64_t offset)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
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
            return mid + offset;
        }
    }

    return -1;
}
int64_t binsearch_offset_index(uint32_t x, uint32_t v[], uint32_t n, uint32_t offset)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = (low + high) >> 1;
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
            return mid + offset;
        }
    }

    return -1;
}


