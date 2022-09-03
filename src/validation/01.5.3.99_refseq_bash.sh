#!bin/bash/

#wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz
#gunzip gene2refseq.gz 
cd /binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/01_crossMapping/01.5_analysis

for i in {1..7}
do
time head -n $((100000*$i)) refesqs.txt|tail -n 100000 |grep -f - gene2refseq|awk '{print $1"\t"$2"\t"$4"\t"$14"\t"$19"\t"$13"\t"$14}'>>mappings.tsv
done

sort -u mappings.tsv > mappings_uniq.tsv

#wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz
