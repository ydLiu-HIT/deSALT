
python ../Annotation_grouping.py Drosophila_melanogaster.BDGP6.94.gtf
python ../generate_transcriptome.py Drosophila_melanogaster.BDGP6.94_AS.gtf Drosophila_melanogaster.BDGP6.fa AS_transcriptome.fa
python ../generate_transcriptome.py Drosophila_melanogaster.BDGP6.94_SS.gtf Drosophila_melanogaster.BDGP6.fa SS_transcriptome.fa
python ../generate_transcriptome.py Drosophila_melanogaster.BDGP6.94_short.gtf Drosophila_melanogaster.BDGP6.fa short_transcriptome.fa
cat AS_transcriptome.fa SS_transcriptome.fa short_transcriptome.fa >merge_transcriptome.fa

python ../../samscripts/src/fastqfilter.py minlen 200 merge_transcriptome.fa transcriptome_for_simulation.fa

rm Drosophila_melanogaster.BDGP6.94_AS.gtf Drosophila_melanogaster.BDGP6.94_SS.gtf Drosophila_melanogaster.BDGP6.94_short.gtf
rm AS_transcriptome.fa SS_transcriptome.fa short_transcriptome.fa merge_transcriptome.fa

mkdir CCS CLR ONT2D ONT1D

cd CCS
bash ../simulate_CCS.sh
cd ..
#
#cd CLR 
#bash ../simulate_CLR.sh
#cd ..
#
#cd ONT2D
#bash ../simulate_ONT2D.sh
#cd ..
#
#cd ONT1D
#bash ../simulate_ONT1D.sh
#cd ..
