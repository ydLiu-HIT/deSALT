# deSALT:
[![DOI](https://zenodo.org/badge/159485852.svg)](https://zenodo.org/badge/latestdoi/159485852)

deSALT - De Bruijn graph-based Spliced Aligner for Long Transcriptome reads

![deSALT](https://github.com/ydLiu-HIT/deSALT/blob/master/img/deSALT_fig.png)

**The latest version deSALT-v1.5.3 have been published. New version adress bugs pushed by users(e.g. [#19](https://github.com/ydLiu-HIT/deSALT/issues/19), [#17](https://github.com/ydLiu-HIT/deSALT/issues/17), [#16](https://github.com/ydLiu-HIT/deSALT/issues/16), [#20](https://github.com/ydLiu-HIT/deSALT/issues/20), [#25](https://github.com/ydLiu-HIT/deSALT/issues/25)), and improve the speed of deSALT by 15-20%**

**The version deSALT-v1.5.2 some bugs in Annotation_Load.py, and make deSALT usage more readability.**

## Getting started
    git clone --recursive https://github.com/ydLiu-HIT/deSALT.git
    cd deSALT/src/deBGA-master/
    make   ## built deBGA for RdBG-index
    cd ..
    make   ## built deSALT for alignment
    
    ./deSALT index ref.fa index_route
    ./deSALT aln index_route read.fq
    
    or 
    run deSALT directly in the same folder (Executable programs have been built in advance.)
    
    ## install by conda
    conda install -c bioconda desalt

### **[Special emphasis]**

1. In current version of deSALT, [deBGA](https://github.com/HongzheGuo/deBGA) is employed for generation of RdBG-index. Some bugs have been corrected and some parameters have been reset (i.e. `START_POS_REF = 0` replaced `START_POS_REF = 2048` in `load_input.h`). 
**Strongly recommend using the executable program in deSALT.**
2. The input reference genome for indexing required the sequence **cutted with a fixed length each line, the length should be no longer than 500bp**. So before indexing, user should check the linewidth of your genome, if large than 500bp, you can change the linewidth with the script I support in the folder (changelinewidth.py), e.g.

	`python changelinewidth.py your_genome.fa changed.fa linewidth(e.g. 80)`

## Introduction
deSALT(de Bruijn graph-based Spliced Aligner for Long Transcriptome reads) is a novel alignment approach with faster speed and sensitive exon identification. Taking the advantages of its novel two pass alignment strategy based on de Bruijn graph-based index, efficient alignment skeleton generation, sensitive exon identification and specifically designed local alignment, deSALT is a fast and accurate RNA-seq long read alignment approach. It has ability to produce high quality full-length read alignment, which is effective to recover the exons and splicing junctions along the entire reads. The workflow of deSALT can be found in `img` folder.

We benchmarked deSALT with 60 simulated datasets having various read length, sequencing error rates and sequencing depth (https://github.com/ydLiu-HIT/deSALT/blob/master/simulation/Sim_data_generation.md). We assessed the aligners with three real sequencing datasets. The first two datasets are from a well-studied CEPH sample (NA12878), and respectively produced by ONT cDNA sequencing (containing 15,152,101 reads and 14,134,831,170 bases in total) and ONT direct RNA sequencing (containing 10,302,647 reads and 10,614,186,428 bases in total). The two datasets are available at https://github.com/nanopore-wgs-consortium/NA12878. The third dataset is from a mouse sample produced by the PacBio platform (SRA Accession Number: SRR6238555; containing 2,269,795 reads and 3,213,849,871 bases in total).

deSALT is open source and free for non-commerical use which is mainly designed by Yadong Liu & Bo Liu and developed by Yadong Liu in Center for Bioinformatics, Harbin Institute of Technology, China.

## Memory usage
deSALT fits most modern servers and workstations and the peak memory footprint depends on the size of reference genome assembly. In practice, 35 Gigabytes, 31 Gigabytes and 3.5 Gigabytes are required for Homo Sapiens(GRCh38), Mus Musculus(GRCm38) and Drosophila melanogaster(DM6) geomes, on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04.

It is worthnoting that, the construction of RdBG-index for large genomes could cost a couple of hours (129, 112 and 4 minutes for human, mouse and fruit fly, respectively) and several tens of GB RAM space (73GB, 63GB and 5.5GB for human, mouse and fruit fly, respectively), depending on the number of distinct k-mers in the genome. This is mainly due to that it needs to extract and sort all the k-mers to construct RdBG at first. However, the index needs only to be built once before use, and we also provide pre-built RdBG-indexes of human, mouse and fruit fly in google drive, users can download the RdBG-index directly.

 - https://drive.google.com/file/d/11E2j1X5jGqKNVtyNHsfPjqnu-fVBgyoi/view?usp=sharing, human GRCh38
 - https://drive.google.com/file/d/1tipOySE-_tmLI4jiy3GZdTj08RhJMXcl/view?usp=sharing, mouse, GRCm38
 - https://drive.google.com/file/d/1ZUg-Yc7oRQQjjJdh4_IJcDGomIVNfJBs/view?usp=sharing, fruit fly, DM6


## Installation
Current version of deSALT has been tested on 64-bit Linux. The source code is written in C, and can be directly download from: https://github.com/ydLiu-HIT/deSALT. The makefile is attached. Use the make command for generating the executable file.

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
Usage: deSALT aln [options] -f <temporary file> <index_route> <read.fa/fq>

	-f <temporary file>           The temporary file for storing alignment skeletons in first pass.
				      If users run two deSALT program in the same time, -f option is necessary.
	<index_route>                 The path of RdBG index.
	<read.fq/fa>                  The input reads in fasta or fastq format.

Algorithm options:

	-t --thread           [INT]	Number of threads. [4]
	-k --index-kmer       [INT]	K-mer length of RdBG-index, the default index kmer-size of deBGA.[21,22]
	-l --seeding-lmer     [INT]	K-mer length of seeding process (no long than RdBG-index). [15]
	-a --local-hash-kmer  [INT]	K-mer length of local hash process. In order to detect spanning exons in 2-pass
					alignment, a local hash query procedure is needed. The hash kmer is recommend no 
					more than 10 bp. [8]
	-s --seed-step        [INT]	The interval of seeding. deSALT extracts kmer at every s bp. [5]
    	-B --batch-size       [INT]	The counts of reads to be processed in one loop. For occupuying less memory, 
					deSALT take only 655350 reads into the memory every time. [655350]
	-n --max-uni-pos      [INT]	Maximum allowed number of hits per seed. If one seed in unipath has more than 50
					copies in reference genome, we will ingore the seed. [50]
	-L --max-readlen      [INT]	Maximum allowed read length. [1000000]
	-i --min-frag-dis     [INT]	Maximum allowed distance of two fragment can be merge. [20]
	-I --max-intron-len   [INT]	Maximum allowed intron length. [200000]
	-c --min-chain-score  [INT]	Minimal skeleton score(match bases minus gap penalty). [30]
	-d --strand-diff      [INT]     The minimal difference of dp score by two strand to make sure the transcript strand. [10]
	-g --max-read-gap     [INT]	Maximum allowed gap in read when generating skeleton. [2000]
	-p --secondary-ratio  [FLOAT]	Min secondary-to-primary score ratio. An alignment can be regard as a secondary
					alignment if secondary_score / primary_score > 0.9. [0.9]
	-e --e-shift          [INT]	The number of downstream (upstream) exons will be processed when left (right) extension. [5]
	-T --trans-strand               Find splicing sites in transcript strand.
    	-G --gtf              [STR]	Provided an annotation file for precise intron donor and acceptor sites.
    	                           	Convert GTF file(now support GTF format only) to fixed format of deSALT by Annotation_Load.py
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
    	-R --noncan           [INT]     Penalty score for non-canonical splice junction sites. [9]

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
### 1. Three different kmer length in deSALT.

`-k index-kmer`: the kmer length to construct the reference de Bruijn graph index(RdBG-index), which organize the reference by unitigs. The default length is 22bp with length range from 20-28bp.

`-l seeding-lmer`: a smaller kmer length than index-kmer for the seeding process. Due the high error rate of long reads (except PacBio ROI reads), large kmers are hard to locate reads in reference genomes. Based on experience, a 15-18bp seeding kmer is best which take cares the search space and enough hits for skeletons generation.

`-a local-hash-kmer`: if one read has a spanning exon due to there are no seed matches between read and spanning exon in the 2-pass alignment, a extreme small kmer is needed to find matches.

In general, `index-kmer > seeding-lmer > local_hash_kmer`. Considering that `seeding-lmer` is smaller than `index-kmer`, when we do a binary search for seeding process, the base of seeding-lmer is the perfix of index-kmer. Thus, a seeding-lmer will have at most **4<sup>(|index-kmer| - |seeding-lmer|)</sup>**. If we use a large index-kmer and a small seeding-lmer, the search space for seeding will be increased fast. Take the consider of time consumption, we give a suggestion of the corresponding index-kmer length and seeding-lmer length by the following table. A smaller index-kmer can be faster and the accuracy remain similar. 

|seeding-lmer | index-kmer|
|:------:|:------:|
|18|22|
|17|22|
|16|21 / 22|
|15|21 / 22|
|14|21|

**What's more, with the limitation of RdBG-index kmer can not be less than 21, two binary search steps are needed for seeding process, so a smaller seeding-lmer(e.g. -l 14) will take more time for alignment than a larger seeding-lmer(e.g. -l 15)**

**Additional, a smaller seed step(`-s`) and a smaller chain score(`-c`) will get a better result, but with the cost of more time.**

### 2. For Iso-seq, Direct RNA-seq, the parameter `-T` can be applied to detect splicing junction site in forward transcript strand only.

### 3. Align different kinds of reads with various sequencing error rates.
`-x read-type:` deSALT can process reads from four main stream platforms,  i.e., ONT 1D reads (error rate: 25%), ONT 2D (1D2) reads (error rate: 12%), PacBio subreads (error rate: 15%) and PacBio ROI reads (error rate: 1%). The total sequencing error rates and the ratios of the sequencing errors (represented as mismatches: insertions: deletions) are configured by referring to previous studies[1-2].

**For error-prone (ONT1D) reads, options `-l 14 -s 2 -x ont1d` are highly recommend to improve the accuracy of exons recovery and full length of transcripts generation.** Of course, it will cost more time than default parameters, but not too much.

**For low error rate (CCS) reads, options `-x ccs -O6,24 -M4` are recommend to give alignments with fewer mismatches/gaps and to open introns more freely.**

### 4. Different specified temporary file path.
`-f temp-file-perfix:` route of temporary files after the first-pass alignment,if users run more than one deSALT program in the same time in the same folder,users must point out different routes of temporary files for each single program! If no, every deSALT program will write temporary data to the same file which will cause crash of program in 2-pass alignment due to inconsistent temporary data. If uses run two deSALT program at the same time within the same folder, different temporary should be specified like follows:
```
deSAL aln -f tmp_path1 -o out1.sam index_route read1.fq   #the first deSALT program
deSAL aln -f tmp_path2 -o out2.sam index_route read2.fq   #the second deSALT program
```

### 5. Alignment with annotations
```
python Annotation_Load.py genome.gtf genome.info   #the annotation file should be in GTF format
deSALT aln -G genome.info index_route read.fa
```

## Simulation benchmarking
All the benchmarks were implemented on a server with Intel Xeon E4280 CPU at 2.0GHZ and 1 Terabytes RAM, running Linux Ubuntu 16.04. The simulated datasets were generated from the reference of three organisms: Homo sapiens GRCh38 (human), Mus musculus GRCm38 (mouse), and Drosophila melanogaster r6 (fruit fly), with corresponding Ensembl gene annotations. There are in total 60 datasets used for the benchmark, and each of them was generated by a specific combination of sequencing model, simulated transcriptome and coverage. 
6 sequencing models were built according to previous studies, to comprehensively benchmark the aligners on the datasets produced by various long read sequencing platforms. For PacBio platforms, there were two models built with fixed parameters: “PacBio ROI reads” (error rate = 2%, mean read length = 2000 bp) and “PacBio subreads” (error rate = 15%, mean read length = 8000 bp). For ONT platforms, four models were built by which two of them were also with fix parameters: “ONT 2D reads” (error rate = 13%, mean read length = 7800 bp) and “ONT 1D reads” (error rate = 25%, mean read length = 7800 bp). And the other two models, “PS-ONT reads” and “NS-ONT reads” were automatically built by PBSim and NanoSim based on a real ONT sequencing dataset (SRA Accession Number: SRR2848544), respectively. we used the two parameter-based models, “ONT 2D reads” and “ONT 1D reads”, as a complement, where the 25% and 12% error rates coincide with typical error rates of ONT 2D (1D2) and 1D reads. The parameters and command lines of PBSim and NanoSim are in Supplementary Notes.

Detailed description of synthetic dataset preparation can be found at https://github.com/ydLiu-HIT/deSALT/blob/master/simulation/Sim_data_generation.md.

The simulated datasets and description used for benchmarking are available at https://drive.google.com/drive/folders/1jk1ddv_QGozumnO_S_f-1lJlOehL3SyW?usp=sharing

## Evaluation on simulated and real datasets
For synthetic dataset, deSALT compares the alignment files (SAM or BAM) to the simulation data generation by PBSIM or NanoSim which have ground truth. In order to reveal the performance of aligners, we take the potential structure of simulation data into consideration and evaluate the results from four aspects. Detailed description of synthetic dataset evaluation can be found at https://github.com/ydLiu-HIT/deSALT/blob/master/evaluation/data_evaluation.md

As for the evaluation of real datasets, we compare the alignment files to corresponding annotations. For each alignment, we find the most overlapped transcript in annotations as aligned transcript, then detect the overlapped exons and calculate the boundaries to make a decision whether the alignment is a good alignment. The details also refer to https://github.com/ydLiu-HIT/deSALT/blob/master/evaluation/data_evaluation.md.

## Citation
Liu, B., Liu, Y., Li, J. et al. deSALT: fast and accurate long transcriptomic read alignment with de Bruijn graph-based index. Genome Biol 20, 274 (2019) doi:10.1186/s13059-019-1895-9

[![DOI](https://zenodo.org/badge/159485852.svg)](https://zenodo.org/badge/latestdoi/159485852)

## Contact
For advising, bug reporting and requiring help, please post on GitHub Issue or contact ydliu@hit.edu.cn.

## Thanks
deSALT relies on the hard work of other projects:
 - The reference de bruijn graph index(RdBG-index):https://github.com/HongzheGuo/deBGA
 -  Dynamic programming in the second phase:https://github.com/lh3/ksw2

## Reference
[1] Weirather JL et al. Comprehensive comparison of Pacific Biosciences and Oxford Nanopore Technologies and their applications to transcriptome analysis. F1000Res (2017), 6: 100. 

[2] Carneiro MO et al. Pacific biosciences sequencing technology for genotyping and variation discovery in human data. BMC Genomics (2012), 13:375. 
