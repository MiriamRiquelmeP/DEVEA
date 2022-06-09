# Create/update species

for example: *Homo sapiens*

```bash
~$ createTree.sh  Hs
~$ createAnnotation.sh http://ftp.ensembl.org/pub/release-....gff3.gz Hs
```

Run R script RscriptUpdateResources.R to create/update all data resources: 

```bash
~$ Rscript --vanilla  RscriptUpdateResources.R 
```



