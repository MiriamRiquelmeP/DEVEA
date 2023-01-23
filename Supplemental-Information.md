DEVEA: an interactive shiny application for **D**ifferential
**E**xpression, data **V**isualisation and **E**nrichment **A**nalysis
of **transcriptomics data**. 

- Supplemental Materials info author:

        Riquelme-Perez M, Perez-Sanz F, Deleuze JF, Escartin C, Bonnet E, Brohard S.
        Supplementary material


# Required R packages:

The following R packages, together with their dependencies are required for running DEVEA either locally or as an 
interactive R/Shiny application. All the packages are at the latest version which is compatible with the R version installed on the machine.

- General functionalities:

htmltools, knitr, BiocParallel, BiocManager, shinyBS, shinyWidgets, shinyalert, shinybusy, shinydashboard, shinydashboardPlus, 
shinyjs, shinymanager, shinythemes, DT, rgl, rglwidget, calibrate, stringr, tidytext, tidyverse, dplyr, rvest, httpuv, MASS.

- Differential expression and enrichment analysis:

DESeq2, apeglm, limma, edgeR, clusterProfiler, fgsea, gage, gageData, pathview, topGO, EnsDb.Hsapiens.v86, 
EnsDb.Mmusculus.v79, EnsDb.Rnorvegicus.v79, org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, org.At.tair.db, BSgenome.Athaliana.TAIR.TAIR9,
AnnotationDbi, AnnotationHub.

- General tools for visualization:

RColorBrewer, Rgraphviz, randomcoloR, ggplot2, ggpubr, igraph, ggraph, ggrepel, gplots, grid, 
gridExtra, plotly, purrr, scales, heatmaply, pheatmap, PoiClaClu, EnhancedVolcano, mychordplot, wordcloud, visNetwork.



# Running DEVEA locally on a machine:

## From the terminal

To run DEVEA locally on a Linux machine from the terminal: 

   1. Download the source code from GitHub:
        
> git clone https://github.com/MiriamRiquelmeP/DEVEA.git
 
> cd ./DEVEA

   2. Open **"script_libraryApp.R"** on Rstudio to install all R libraries that DEVEA needs to work. 
        
   3. Back to the terminal:

   - If you have a **Counting Matrix + Sample Information (CM + SI)** or **DESeqDataSet Object (DO)**, you can access **DESeqDEVEA**:
                
> R -e "shiny::runApp('~/DESeqDEVEA')"

   - If you have a **Gene List + Statistical Values (GL + SV)** or **Gene List (GL)**, you can access **simpleDEVEA**:
                
> R -e "shiny::runApp('~/simpleDEVEA')"


## From Rstudio

To run DEVEA from RStudio on a local computer:

   1. **Download** the compressed folder from the DEVEA github repository (https://github.com/MiriamRiquelmeP/DEVEA/).
        
   2. **Unzip** the package and **access** the main folder.
        
   3. From Rstudio, open **script_libraryApp.R** inside DEVEA to install all R libraries that the tool needs to work. 
        
   4. Depending on your input data:
                
        If you have a **Counting Matrix + Sample Information (CM + SI)** or **DESeqDataSet Object (DO)**, you can open app.R from **DESeqDEVEA** folder.
        
        If you have a **Gene List + Statistical Values (GL + SV)** or **Gene List (GL)**, you can open app.R from **simpleDEVEA** folder.
        
   5. Click **Run** on RStudio.
        

