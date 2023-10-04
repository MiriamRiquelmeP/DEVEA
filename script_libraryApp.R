install.packages("dplyr")
install.packages("rvest")
library(dplyr)
library(rvest)
if(!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
install.packages("httpuv")
install.packages("MASS")
library (MASS)
require(dplyr)

InstalledPackages <- function(packages = NULL){
  CRANpackages <- available.packages() %>%
    as.data.frame() %>%
    dplyr::select(Package) %>%
    mutate(source = 'CRAN')
  
  url <- 'https://www.bioconductor.org/packages/release/bioc/'
  biocPackages <- url %>%
    read_html() %>%
    html_table() %>%
    .[[1]] %>%
    dplyr::select(Package) %>%
    mutate(source = 'BioConductor')
  
  all_packages <- bind_rows(CRANpackages, biocPackages)
  rownames(all_packages) <- NULL
  if(is.null(packages)){
    packages <- installed.packages() %>% as.data.frame()
  }else{
    packages <- data.frame(Package=packages, Priority = NA)
  }
  installed <- all_packages %>% filter( Package %in% packages$Package )
  NotIdentified <- packages[!(packages$Package %in% installed$Package), ] %>% filter(is.na(Priority))
  CRAN <- installed %>% filter(source=="CRAN") %>% dplyr::select(Package)
  Bioconductor <- installed %>% filter(source=="BioConductor") %>% dplyr::select(Package)
  pkg <- list(cran = CRAN$Package, bioc = Bioconductor$Package, nid = NotIdentified$Package )
  return(pkg)
}

packages <- c("AnnotationDbi","chorddiag","DESeq2","DT","ensembldb","EnsDb.Mmusculus.v79","EnsDb.Hsapiens.v86","EnsDb.Rnorvegicus.v79","fgsea","ggpubr","ggrepel",
              "grid","gridExtra","heatmaply","limma","mychordplot","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","org.At.tair.db","BSgenome.Athaliana.TAIR.TAIR9", "pheatmap","plotly","purrr",
              "randomcoloR","RColorBrewer","rgl","rglwidget","scales","shinyalert","shinyBS","shinybusy","shinydashboard","shinydashboardPlus",
              "shinyjs","shinythemes","shinyWidgets","shinymanager","stringr","tidyverse","tidytext","visNetwork","wordcloud","flexdashboard",
              "ggraph","igraph","edgeR","PoiClaClu","ggplot2","knitr","apeglm","calibrate","AnnotationHub",'topGO',"Rgraphviz","BiocParallel",
              "clusterProfiler","pathview","gage","gageData","EnhancedVolcano","gplots","htmltools","dplyr","tidyr","GO.db","GOplot","karyoploteR")

pkg <- InstalledPackages(packages = packages)
pkg$cran
pkg$bioc
pkg$nid

for(i in pkg$cran){ if(!require(i, character.only = TRUE)){install.packages(i)} }
for(i in pkg$bioc){ if(!require(i, character.only = TRUE)){BiocManager::install(i)} } 

install.packages("devtools")
devtools::install_github("mattflor/chorddiag")
devtools::install_github("fpsanz/mychordplot")
devtools::install_github("jokergoo/circlize")




