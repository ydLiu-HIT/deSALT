# deSALT:
deSALT - De Bruijn graph-based Spliced Aligner for Long Transcriptome reads

---
### Getting started
    git clone https://github.com/ydLiu-HIT/deSALT/
    cd deSALT
    make
    ./deSALT index ref.fa index_route
    ./deSALT aln index_route read.fq

---
### Introduction
deSALT(de Bruijn graph-based Spliced Aligner for Long Transcriptome reads) is a novel alignment approach with faster speed and sensitive exon identification. Taking the advantages of its novel two pass alignment strategy based on de Bruijn graph-based index, efficient alignment skeleton generation, sensitive exon identification and specifically designed local alignment, deSALT is a fast and accurate RNA-seq long read alignment approach. It has ability to produce high quality full-length read alignment, which is effective to recover the exons and splicing junctions along the entire reads.

We benchmarked deSALT with 36 simulated datasets having various read length, sequencing error rates and sequencing depth (simulation workflow links). deSALT also assessed the ability of aligners by two real RNA-seq datasets produced by Oxford Nanopore and PacBio platform. One from the well-studied human sample NA12878 by ONT technology and another from real mouse sample (SRR6238555) by SMRT technology.

deSALT is open source and free for non-commerical use which is mainly designed by Yadong Liu & Bo Liu and developed by Yadong Liu in Center for Bioinformatics, Harbin Institute of Technology, China.

---
### Memory usage
The memory usage of deSALT can fit the configuraions of most modern servers and workstations. Its peak memory footprint depends on the size of reference genome mainly due to the generation of RdBG-index. 35 Gigabytes, 31 Gigabytes and 3.5 Gigabytes are required for Homo Sapiens(GRCh38), Mus Musculus(GRCm38) and Drosophila melanogaster(DM6) geomes, on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04. For instance, the peak memory is about 37.65 Gigabytes for human NA12878 dataset.

---
### Installation
Current version of deSALT needs to be run on Linux operating system. The source code is written in C, and can be directly download from: https://github.com/ydLiu-HIT/deSALT. The makefile is attached. Use the make command for generating the executable file.

Moreover, in current version of deSALT, we employed the deBGA mapper (https://github.com/HongzheGuo/deBGA) for generation RdBG-index. To be more user-firendly, we have built the source code of deBGA (version 0.1) into that of deSALT. And set START_POS_REF = 0 instead of START_POS_REF = 2048 in load_input.h source file of deBGA, i.e. we set the reference genome starting from 0 rather than 2048.

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
	-s --seed-step        [INT]	Interval of seeding. [5]
    	-B --batch-size       [INT]	The number of reads to be processed in one loop. [65535]
	-n --max-uni-pos      [INT]	Maximum allowed number of hits per seed. [50]
	-l --max-readlen      [INT]	Maximum allowed read length. [1000000]
	-i --min-frag-dis     [INT]	Maximum allowed distance of two fragment can be merge. [20]
	-I --max-intron-len   [INT]	Maximum allowed intron length. [200000]
	-c --min-chain-score  [INT]	Minimal skeleton score(match bases minus gap penalty). [30]
	-g --max-read-gap     [INT]	Maximum allowed gap in read when chaining. [2000]
	-p --secondary-ratio  [FLOAT]	Min secondary-to-primary score ratio. [0.9]
	-p --e-shift          [INT]	The shift of downstream and upstream when alignment. [5]
    	-G --gtf              [STR]	Provided an annotation file for precise intron donor and acceptor sites.
    	                           	The release of annotation file and reference genome must the same!
	-x --read-type        [STR]	Specifiy the type of reads and set multiple paramters unless overriden.
	                           	[null] default parameters.
	                           	ccs (PacBio SMRT CCS reads): error rate 1%
	                              	clr (PacBio SMRT CLR reads): error rate 15%
	                              	ont1d (Oxford Nanopore 1D reads): error rate > 20%
	                               	ont2d (Oxford Nanopore 2D reads): error rate > 12%

Scoring options:

	-O --open-pen         [INT(,INT)]	
					Gap open penealty. [2, 32]
	-E --ext-pen          [INT(,INT)]	
					Gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2}. [1, 0]
	-m --match-score      [INT]    	Match score for SW-alginment. [1]
	-M --mis-score        [INT]    	Mismatch score for SW-alignment. [2]
	-z --zdrop            [INT(,INT)]
					Z-drop score for splice/non-splice alignment. [400]
	-w --band-width       [INT]    	Bandwidth used in chaining and DP-based alignment. [500]

Output options:

	-N --top-num-aln      [INT]    	Max allowed number of secondary alignment. [5]
	-Q --without-qual              	Don't output base quality in SAM.
	-f --temp-file-perfix [STR]    	Perfix of temp file during the program. [./1pass_anchor]
		                        If you run more than one tgs program in the same time,
		                        you must point at a different perfix of temp file for each program!
	-o --output           [STR]     Output file (SAM format). [./aln.sam]
```


