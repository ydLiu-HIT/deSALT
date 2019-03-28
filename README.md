# deSALT:
deSALT - De Bruijn graph-based Spliced Aligner for Long Transcriptome reads

## Getting started
    git clone https://github.com/ydLiu-HIT/deSALT/
    cd deSALT/src
    make
    ./deSALT index ref.fa index_route
    ./deSALT aln index_route read.fq

## Introduction
deSALT(de Bruijn graph-based Spliced Aligner for Long Transcriptome reads) is a novel alignment approach with faster speed and sensitive exon identification. Taking the advantages of its novel two pass alignment strategy based on de Bruijn graph-based index, efficient alignment skeleton generation, sensitive exon identification and specifically designed local alignment, deSALT is a fast and accurate RNA-seq long read alignment approach. It has ability to produce high quality full-length read alignment, which is effective to recover the exons and splicing junctions along the entire reads.

We benchmarked deSALT with 36 simulated datasets having various read length, sequencing error rates and sequencing depth (simulation workflow links). deSALT also assessed the ability of aligners by two real RNA-seq datasets produced by Oxford Nanopore and PacBio platform. One from the well-studied human sample NA12878 by ONT technology and another from real mouse sample (SRR6238555) by SMRT technology.

deSALT is open source and free for non-commerical use which is mainly designed by Yadong Liu & Bo Liu and developed by Yadong Liu in Center for Bioinformatics, Harbin Institute of Technology, China.

## Memory usage
The memory usage of deSALT can fit the configuraions of most modern servers and workstations. Its peak memory footprint depends on the size of reference genome mainly due to the generation of RdBG-index. 35 Gigabytes, 31 Gigabytes and 3.5 Gigabytes are required for Homo Sapiens(GRCh38), Mus Musculus(GRCm38) and Drosophila melanogaster(DM6) geomes, on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04. For instance, the peak memory is about 37.65 Gigabytes for human NA12878 dataset.

## Installation
Current version of deSALT needs to be run on Linux operating system. The source code is written in C, and can be directly download from: https://github.com/ydLiu-HIT/deSALT. The makefile is attached. Use the make command for generating the executable file.

Moreover, in current version of deSALT, we employed the deBGA mapper (https://github.com/HongzheGuo/deBGA) for generation RdBG-index. To be more user-firendly, we have built the source code of deBGA (version 0.1) into that of deSALT. I correct some bugs of deBGA and set `START_POS_REF = 0` instead of `START_POS_REF = 2048` in `load_input.h` source file of deBGA, i.e. we set the reference genome starting from 0 rather than 2048.

## Synopsis
Reference genome indexing
```
deSALT index ref.fa <index_route>
```
	
Read alignment
```
deSALT aln <index_route> read.fa/fq
```

## Commands and options
```
Algorithm options:

	-t --thread           [INT]	Number of threads. [4]
	-K --index-kmer       [INT]	K-mer length of RdBG-index, the default index kmer-size of deBGA.[22]
	-k --seeding-kmer     [INT]	K-mer length of seeding process (no long than RdBG-index). [15]
	-a --local-hash-kmer  [INT]	K-mer length of local hash process. In order to detect spanning exons in 2-pass
					alignment, a local hash query procedure is needed. The hash kmer is recommend no 
					more than 10 bp. [8]
	-s --seed-step        [INT]	The interval of seeding. deSALT extracts kmer at every s bp. [5]
    	-B --batch-size       [INT]	The counts of reads to be processed in one loop. For occupuying less memory, 
					deSALT take only 655350 reads into the memory every time. [655350]
	-n --max-uni-pos      [INT]	Maximum allowed number of hits per seed. If one seed in unipath has more than 50
					copies in reference genome, we will ingore the seed. [50]
	-l --max-readlen      [INT]	Maximum allowed read length. [1000000]
	-i --min-frag-dis     [INT]	Maximum allowed distance of two fragment can be merge. [20]
	-I --max-intron-len   [INT]	Maximum allowed intron length. [200000]
	-c --min-chain-score  [INT]	Minimal skeleton score(match bases minus gap penalty). [30]
	-d --strand-diff      [INT]     The minimal difference of dp score by two strand to make sure the transcript strand. [10]
	-g --max-read-gap     [INT]	Maximum allowed gap in read when generating skeleton. [2000]
	-p --secondary-ratio  [FLOAT]	Min secondary-to-primary score ratio. An alignment can be regard as a secondary
					alignment if secondary_score / primary_score > 0.9. [0.9]
	-p --e-shift          [INT]	The number of downstream (upstream) exons will be processed when left (right) extension. [5]
    	-G --gtf              [STR]	Provided an annotation file for precise intron donor and acceptor sites.
    	                           	The version information of annotation file and reference genome must the same!
	-x --read-type        [STR]	Specifiy the type of reads and set multiple paramters unless overriden.
	                           	[null] default parameters. error rate 13% 
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
	-f --temp-file-perfix [STR]    	Route of temporary files after the first-pass alignment. [./skeletons]
		                        If you run more than one deSALT program in the same time,
		                        you must point out different routes of temporary files for each program!
					If no, every deSALT program will write temporary data to the same file which
					will cause crash of program in 2-pass alignment due to inconsistent temporary data.
	-o --output           [STR]     Output file (SAM format). [./aln.sam]
```
## Important options
#### 1. Three different kmer length in deSALT.

`-K index-kmer`: the kmer length to construct the reference de Bruijn graph index(RdBG-index), which organize the reference by unitigs. The default length is 22bp with length range from 20-28bp.

`-k seeding-kmer`: a smaller kmer length than index-kmer for the seeding process. Due the high error rate of long reads (except PacBio ROI reads), large kmers are hard to locate reads in reference genomes. Based on experience, a 15-18bp seeding kmer is best which take cares the search space and enough hits for skeletons generation.

`-a local-hash-kmer`: if one read has a spanning exon due to there are no seed matches between read and spanning exon in the 2-pass alignment, a extreme small kmer is needed to find matches.

In general, `index-kmer > seeding-kmer > local_hash_kmer`. Considering that `seeding-kmer` is smaller than `index-kmer`, when we do a binary search for seeding process, the base of seeding-kmer is the perfix of index-kmer. Thus, a seeding-kmer will have at most 4(|index-kmer| - |seeding-kmer|). If we use a large index-kmer and a small seeding-kmer, the search space for seeding will be increased fast. Take the consider of time consumption, we give a suggestion of the corresponding index-kmer length and seeding-kmer length by the following table.

|seeding-kmer | index-kmer|
|:------:|:------:|
|18|22|
|17|22|
|16|20 / 22|
|15|20 / 22|
|14|20|

#### 2. Different specified temporary file path.
`-f temp-file-perfix:` route of temporary files after the first-pass alignment,if users run more than one deSALT program in the same time in the same folder,users must point out different routes of temporary files for each single program! If no, every deSALT program will write temporary data to the same file which will cause crash of program in 2-pass alignment due to inconsistent temporary data. If uses run two deSALT program at the same time within the same folder, different temporary should be specified like follows:
```
deSAL aln -f tmp_path1 -o out1.sam index_route read1.fq   #the first deSALT program
deSAL aln -f tmp_path2 -o out2.sam index_route read2.fq   #the second deSALT program
```


## Simulation benchmarking
In the simulation study, we simulated 36 RNA-seq long read datasets with various sequencing error rates and read lengths (refers to supplementary) to mimic the datasets from various platforms, i.e., ONT 1D reads (error rate: 25%, mean read length: 7800 bp), ONT 2D reads (error rate: 12%, mean read length: 7800 bp), PacBio subreads (error rate: 15%, mean read length: 8000 bp) and PacBio ROI reads (error rate: 2%, mean read length: 2000 bp). For each of the platforms, there are respectively 9 datasets from 3 species (human, mouse and fruitfly) and in 3 sequencing depths (4X, 10X, and 30X). All the datasets were produced by PBSim based on Ensembl gene annotations (human: GRCh38, version 94, mouse: GRCm38, version 94 and fruitfly: BDGP6, version 94).

Due to there is no well-studied simulator for noisy long RNA-seq reads, we used PBSIM to generate synthetic datasets by a set of transcripts generated from a particular reference genome and corresponding annotations inspired by RNAseqEval project (https://github.com/kkrizanovic/RNAseqEval). Detailed description of synthetic dataset preparation can be found at https://github.com/ydLiu-HIT/deSALT/blob/master/simulation/Sim_data_generation.md.

The simulated datasets and description used for benchmarking are available at https://drive.google.com/drive/folders/16RpDYkdTCwOHmvoWNehnUp7nxrNt_Q7T?usp=sharing

## Contact
For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn or ydliu@hit.edu.cn
