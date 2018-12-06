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


