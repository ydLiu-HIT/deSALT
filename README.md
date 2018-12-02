# deSALT:
deSALT - De Bruijn graph-based Spliced Aligner for Long Transcriptome reads

---
### Getting started
    git clone https://github.com/ydLiu-HIT/deSALT/
    cd deSALT
    make
    ./deSALT index ...
    ./deSALT aln ...

---
### Introduction
deSALT(de Bruijn graph-based Spliced Aligner for Long Transcriptome reads) is a novel alignment approach with faster speed and sensitive exon identification. Taking the advantages of its novel two pass alignment strategy based on de Bruijn graph-based index, efficient alignment skeleton generation, sensitive exon identification and specifically designed local alignment, deSALT is a fast and accurate RNA-seq long read alignment approach. It has ability to produce high quality full-length read alignment, which is effective to recover the exons and splicing junctions along the entire reads.

We benchmarked deSALT with simulated and real datasets having various read length and sequencing error rates. For simulated  benchmark, we developed experiments from four different aspects.

    1. simulate reads with different abundance of gene expressions
    2. simulate reads with

---
### Synopsis
Reference genome indexing
```
deSALT index ref.fa <index_route>
```
	
Read alignment
```
deSALT aln <index_route> read.fa/fq
```

---
### Commands and options
```
Algorithm options:

	-t --thread           [INT]	Number of threads. [1]
	-K --index-kmer       [INT]	K-mer length of deBGA index. [22]
	-k --seeding-kmer     [INT]	K-mer length of seeding process. [15]
	-a --local-hash-kmer  [INT]	K-mer length of local hash process. [8]
	-s --seed-step        [INT]	Interval of seeding. [%u]\n", SEED_STEP);
    	-B --batch-size       [INT]	The number of reads to be processed in one loop. [65535]
	-n --max-uni-pos      [INT]	Maximum allowed number of hits per seed. [50]
	-l --max-readlen      [INT]	Maximum allowed read length. [1000000]
	-i --min-frag-dis     [INT]	Maximum allowed distance of two fragment can be connect. [%u]\n", MIN_FRAG_DIS);
	-I --max-intron-len   [INT]	Maximum allowed intron length. [%u]\n", SPLICDISTANCE);
	-c --min-chain-score  [INT]	Minimal skeleton score(match bases minus gap penalty). [%u]\n", MIN_CHAIN_SCORE);
	-g --max-read-gap     [INT]	Maximum allowed gap in read when chaining. [%u]\n", MAX_READ_JOIN_GAP);
	-p --secondary-ratio  [FLOAT]	Min secondary-to-primary score ratio. [%.2f]\n", SECONDARY_TO_PRIMARY);
	-p --e-shift          [INT]	The shift of downstream and upstream when alignment. [%u]\n", E_SHIFT);
    	-G --gtf              [STR]	Provided an annotation file for precise intron donor and acceptor sites.\n");
    	                           	The release of annotation file and reference genome must the same!\n\n");
	-x --read-type        [STR]	Specifiy the type of reads and set multiple paramters unless overriden.\n");
	                           	[null] default parameters.\n");
	                           	ccs (PacBio SMRT CCS reads): error rate 1%%\n");
	                              	clr (PacBio SMRT CLR reads): error rate 15%%\n");
	                              	ont1d (Oxford Nanopore 1D reads): error rate > 20%%\n");
	                               	ont2d (Oxford Nanopore 2D reads): error rate > 12%%\n");

Scoring options:

	-O --open-pen         [INT(,INT)]	
					Gap open penealty. [%u,%u]\n", GAP_OPEN, GAP_OPEN2);
	-E --ext-pen          [INT(,INT)]	
					Gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2}. [%u,%u]\n", GAP_EXT, GAP_EXT2);
	-m --match-score      [INT]    	Match score for SW-alginment. [%u]\n", MATCH_SCORE);
	-M --mis-score        [INT]    	Mismatch score for SW-alignment. [%u]\n", MISMATCH_SCORE);
	-z --zdrop            [INT(,INT)]
					Z-drop score for splice/non-splice alignment. [%u]\n", ZDROP_SCORE);
	-w --band-width       [INT]    	Bandwidth used in chaining and DP-based alignment. [%u]\n\n", BANDWIDTH);

Output options:

	-N --top-num-aln      [INT]    	Max allowed number of secondary alignment. [%u]\n", TOP_NUM_ALN);
	-Q --without-qual              	Don't output base quality in SAM\n");
	-f --temp-file-perfix [STR]    	Perfix of temp file during the program. [%s]\n", TEMP_FILE_PERFIRX);
		                        If you run more than one tgs program in the same time, \n");
		                        you must point at a different perfix of temp file for each program!\n");
	-o --output           [STR]     Output file (SAM format). [%s]\n", OUTPUT);
```


