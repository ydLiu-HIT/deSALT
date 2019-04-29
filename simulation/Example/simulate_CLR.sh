mkdir group1 group2 group3

cd group1
pbsim   ../../transcriptome_for_simulation.fa \
        --prefix SimG1_S \
        --data-type CLR \
        --model_qc ../../../pbsim-1.0.3-Linux-amd64/data/model_qc_clr \
        --length-mean 3000 \
        --length-sd 2300 \
        --length-min 100 \
        --length-max 25000 \
        --difference-ratio 1:12:2 \
        --accuracy-mean 0.85 \
        --accuracy-min 0.75 \
        --depth 4

cd ../group2/
pbsim   ../../transcriptome_for_simulation.fa \
        --prefix SimG2_S \
        --data-type CLR \
        --model_qc ../../../pbsim-1.0.3-Linux-amd64/data/model_qc_clr \
        --length-mean 3000 \
        --length-sd 2300 \
        --length-min 100 \
        --length-max 25000 \
        --difference-ratio 1:12:2 \
        --accuracy-mean 0.85 \
        --accuracy-min 0.75 \
        --depth 10

cd ../group3/
pbsim   ../../transcriptome_for_simulation.fa \
        --prefix SimG3_S \
        --data-type CLR \
        --model_qc ../../../pbsim-1.0.3-Linux-amd64/data/model_qc_clr \
        --length-mean 3000 \
        --length-sd 2300 \
        --length-min 100 \
        --length-max 25000 \
        --difference-ratio 1:12:2 \
        --accuracy-mean 0.85 \
        --accuracy-min 0.75 \
        --depth 30

cd ..
cat group1/*.fastq > dataset_sim_dm_CLR_g1.fastq
cat group2/*.fastq > dataset_sim_dm_CLR_g2.fastq
cat group3/*.fastq > dataset_sim_dm_CLR_g3.fastq

python ../tran_qname.py dataset_sim_dm_CLR_g1.fastq SimG1_S g1.fastq
mv g1.fastq dataset_sim_dm_CLR_g1.fastq
python ../tran_qname.py dataset_sim_dm_CLR_g2.fastq SimG2_S g2.fastq
mv g2.fastq dataset_sim_dm_CLR_g2.fastq
python ../tran_qname.py dataset_sim_dm_CLR_g3.fastq SimG3_S g3.fastq
mv g3.fastq dataset_sim_dm_CLR_g3.fastq

rm group1/*.fastq
rm group2/*.fastq
rm group3/*.fastq

cat dataset_sim_dm_CLR_g1.fastq dataset_sim_dm_CLR_g2.fastq dataset_sim_dm_CLR_g3.fastq > dataset_sim_dm_CLR.fastq
