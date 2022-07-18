library(dplyr)
library(rvest)

paquetesInstalados <- function(paquetes = NULL){
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
  if(is.null(paquetes)){
    paquetes <- installed.packages() %>% as.data.frame()
  }else{
    paquetes <- data.frame(Package=paquetes, Priority = NA)
  }
  installed <- all_packages %>% filter( Package %in% paquetes$Package )
  noIdentificados <- paquetes[!(paquetes$Package %in% installed$Package), ] %>% filter(is.na(Priority))
  CRAN <- installed %>% filter(source=="CRAN") %>% select(Package)
  Bioconductor <- installed %>% filter(source=="BioConductor") %>% select(Package)
  pkg <- list(cran = CRAN$Package, bioc = Bioconductor$Package, nid = noIdentificados$Package )
  return(pkg)
}

# /usr/local/lib/R/site-library
# /home/kirk/R/x86_64-pc-linux-gnu-library/4.0

#remotes::install_version("RSQLite", version = "2.2.5")

#paquetes <- system("ls /home/miriam/R/x86_64-pc-linux-gnu-library/3.6", intern = TRUE)

paquetes <- c("AnnotationDbi","chorddiag","DESeq2","DT","EnsDb.Mmusculus.v79","EnsDb.Hsapiens.v86","EnsDb.Rnorvegicus.v79","fgsea","ggpubr","ggrepel",
              "grid","gridExtra","heatmaply","limma","mychordplot","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","pheatmap","plotly","purrr",
              "randomcoloR","RColorBrewer","rgl","rglwidget","scales","shinyalert","shinyBS","shinybusy","shinydashboard","shinydashboardPlus",
              "shinyjs","shinythemes","shinyWidgets","shinymanager","stringr","tidyverse","tidytext","visNetwork","wordcloud","flexdashboard",
              "ggraph","igraph","edgeR","PoiClaClu","ggplot2","knitr","apeglm","calibrate","AnnotationHub",'topGO',"Rgraphviz","BiocParallel",
              "clusterProfiler","pathview","gage","gageData","EnhancedVolcano","gplots","htmltools","dplyr","tidyr","GO.db","GOplot","karyoploteR")

pkg <- paquetesInstalados(paquetes = paquetes)
pkg$cran
pkg$bioc
pkg$nid

for(i in pkg$cran){ if(!require(i, character.only = TRUE)){install.packages(i)} }
for(i in pkg$bioc){ if(!require(i, character.only = TRUE)){BiocManager::install(i)} } 

devtools::install_github("mattflor/chorddiag")
devtools::install_github("fpsanz/mychordplot")




