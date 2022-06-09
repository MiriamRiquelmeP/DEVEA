#!/bin/bash

gff=$1 # 'ftp://ftp.ensembl.org/pub/release-100/gff3/mus_musculus/Mus_musculus.GRCm38.100.chr.gff3.gz'
specie=$2 # 'Mm'
curl -L ${gff} >${specie}.gtf.gz
zcat ${specie}.gtf.gz | awk '$3=="gene"{print $1,$4,$5,$9}' | awk 'BEGIN{OFS="\t"}{split($4,a,";");print a[1],$1,$2,$3}' | sed 's/ID=gene://g' >./resources/${specie}/cytoband/${specie}_annot.txt
rm ${specie}.gtf.gz
