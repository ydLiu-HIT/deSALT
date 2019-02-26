#!/bin/bash

python ../dm/Annotation_Load.py /data/ydliu/human_ONT_2D_simulate/generator_largeset/Homo_sapiens.GRCh38.92.chr_gfread.gtf
python /home/ydliu/RNAseqEval_old/generate_transcriptome.py AS.gtf /data/ydliu/human_ONT_2D_simulate/hg38_filter.fa AS_transcriptome.fa
python /home/ydliu/RNAseqEval_old/generate_transcriptome.py SS.gtf /data/ydliu/human_ONT_2D_simulate/hg38_filter.fa SS_transcriptome.fa
python /home/ydliu/RNAseqEval_old/generate_transcriptome.py small.gtf /data/ydliu/human_ONT_2D_simulate/hg38_filter.fa small_transcriptome.fa
cat AS_transcriptome.fa SS_transcriptome.fa small_transcriptome.fa >simulate_human.fa
python /home/ydliu/RNAseqEval_old/samscripts/src/fastqfilter.py minlen 200 simulate_human.fa for_simulate_human.fa

cd ONT2D
bash pipeline.sh
cd ..

cd ONT1D
bash pipeline.sh
cd ..

cd CCS
bash pipeline.sh
cd ..

cd CLR
bash pipeline.sh
cd ..

