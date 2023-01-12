#!/bin/bash

# $1 puede ser Mm o Hs actualmente

species=$1
cd ./resources/$species
echo 'creating ensembl file...'
wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz -O ensembl.gtf.gz
gunzip ensembl.gtf.gz
cat ensembl.gtf| grep -v "#" | awk 'BEGIN{FS="\t"} $3=="gene"{print $9}' | awk 'BEGIN{FS="; ";OFS="\t"}{print $1,$3,$6}' | sed -E 's/(gene_id )|(gene_name )|(havana_gene )|(")//g' >ensembl.txt
rm ensembl.gtf
echo 'done ensembl'
echo 'creating gencode file...'
# Comprehensive gene annotation ALL in gff3 format
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.chr_patch_hapl_scaff.annotation.gff3.gz -O gencode.gff.gz
gunzip gencode.gff.gz
cat gencode.gff| grep -v "#" | awk 'BEGIN{FS="\t"} $3 == "gene"{print $9}' | awk 'BEGIN{FS=";";OFS="\t";ORS="\t"}{for(i=1;i<=NF;i++){if($i~/ID|gene_name|havana_gene/){print $i} }{print "\n"} }' | sed -E 's/(ID=)|(gene_name=)|(havana_gene=)//g' | sed 's/^\t//g' | gawk 'BEGIN{FS="\t";OFS="\t"}{b=gensub(/(.*)\.[0-9]/, "\\1", "g", $1);c=gensub(/(.*)\.[0-9]/, "\\1", "g", $3);print b,$2,c}' >gencode.txt 
rm gencode.gff
echo 'done gencode'
echo 'creating ncbi file...'
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz -O ncbi.gff.gz
gunzip ncbi.gff.gz
cat ncbi.gff | grep -v "#" | awk 'BEGIN{FS="\t"}{if($3=="gene"||$3=="pseudogene"){print $9}}' |gawk 'BEGIN{FS=";";ORS="\t"}{for(i=1;i<NF;i++){if($i~/ID=|Dbxref=|Name=/){print $i}}{print "\n"} }'| sed 's/^\t//g' |  sed -E 's/(ID=gene-)|(Dbxref=GeneID:)|(Name=)//g'  | sed 's/\t$//g' | sed -E 's/(.*\t[0-9]+).*(\t.*)/\1 \2/g' | gawk 'BEGIN{FS="\t";OFS="\t"}{b=gensub(/(.*)-.*/,"\\1","g",$1);print b,$2,$3}' | uniq > ncbi.txt
rm ncbi.gff
echo 'done ncbi'
cd ../..

Rscript --vanilla geneIdCreator.R $species
cd ./resources/$species
rm ensembl.txt gencode.txt ncbi.txt

echo 'Done. Bye'
