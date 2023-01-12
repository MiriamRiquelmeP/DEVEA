alibrary(dplyr)
library(rvest)

InstalledPackages <- function(packages = NULL){
  CRANpackages <- available.packages() %>%
    as.data.frame() %>%
    select(Package) %>%
    mutate(source = 'CRAN')
  
  url <- 'https://www.bioconductor.org/packages/release/bioc/'
  biocPackages <- url %>%
    read_html() %>%
    html_table() %>%
    .[[1]] %>%
    select(Package) %>%
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
  CRAN <- installed %>% filter(source=="CRAN") %>% select(Package)
  Bioconductor <- installed %>% filter(source=="BioConductor") %>% select(Package)
  pkg <- list(cran = CRAN$Package, bioc = Bioconductor$Package, nid = NotIdentified$Package )
  return(pkg)
}

# /usr/local/lib/R/site-library
# /home/kirk/R/x86_64-pc-linux-gnu-library/4.0
#remotes::install_version("RSQLite", version = "2.2.5")
#packages <- system("ls /home/miriam/R/x86_64-pc-linux-gnu-library/3.6", intern = TRUE)

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

devtools::install_github("mattflor/chorddiag")
devtools::install_github("fpsanz/mychordplot")




