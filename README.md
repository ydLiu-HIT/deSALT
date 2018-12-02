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
