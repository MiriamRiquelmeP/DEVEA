## NEW FUNCTIONS ###############################################
## validate geneList #########
validateGeneList <- function(list, specie, annotation){
    aa <- gsub(",$", "", 
               gsub(",,", ",", 
                              gsub("\n", ",", list, perl = TRUE),
                    perl =TRUE),
               perl = TRUE)
    aa <- gsub(" ", "", aa)
    geneList <-data.frame(gene = unlist(strsplit(aa, ",")), stringsAsFactors = FALSE)
    geneList <- geneList %>% drop_na(gene)  %>% filter(gene != "") %>% as.data.frame()
    geneList <- data.frame(gene = unique(geneList$gene))
    return(geneList)
}
## Formatear datos y prepararlos para enrichment #########################
formatData <- function(genelist, specie, annotation){
  require("EnsDb.Mmusculus.v79")
  require("org.Mm.eg.db")
  require("EnsDb.Hsapiens.v86")
  require("org.Hs.eg.db")
  require("EnsDb.Rnorvegicus.v79")
  require("org.Rn.eg.db")
  genelist <- as.data.frame(genelist)
  if(annotation=="ensg"){ann="ENSEMBL"}
    else{ann="SYMBOL"}
  ## listado simple
  if( dim(genelist)[2]==1 ){
    names(genelist)<-ann
  }
  ## listado triple (gen, logfc, padj)
  if( dim(genelist)[2]==3){
    names(genelist) <- c(ann, "log2FC","pval")
  }

  if(specie=="Mm"){
      ensdb <- EnsDb.Mmusculus.v79
      orgdb <- org.Mm.eg.db
  }
    else{
        ensdb <- EnsDb.Hsapiens.v86
        orgdb <- org.Hs.eg.db
    }
  # genelist[[ann]] <- toupper(genelist[[ann]]) # todo a mayusculas
  # if( ann=="SYMBOL" & (specie == "Mm" | specie == "Rn") ){ # si es symbol y mouse o rattus a capital
  #   genelist[[ann]] <- stringr::str_to_title(genelist[[ann]])
  #   }
  annot <- NULL
  genes <- as.character(genelist[[ann]])
  if(ann=="ENSEMBL"){
    annot$ENSEMBL <- data.frame(ENSEMBL = genes, stringsAsFactors = FALSE)
    annot$SYMBOL <-  mapIds(ensdb, keys=genes, column="SYMBOL",keytype="GENEID")
    annot$SYMBOL1 <- mapIds(orgdb, keys = genes, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
    annot$description <- mapIds(orgdb, keys = genes, column = 'GENENAME', keytype = 'ENSEMBL', multiVals = 'first')
    annot <- as.data.frame(annot)
    consensus <- data.frame('Symbol'= ifelse(!is.na(annot$SYMBOL), as.vector(annot$SYMBOL),
                                             ifelse(!is.na(annot$SYMBOL1),as.vector(annot$SYMBOL1),
                                                    as.vector(annot$ENSEMBL))), stringsAsFactors = F)
    annot$consensus <- consensus$Symbol
    entrez1 <- mapIds(orgdb, keys = annot$consensus, column = "ENTREZID", keytype = "SYMBOL")
    entrez2 <- mapIds(orgdb, keys = as.character(annot$ENSEMBL),
                      column = "ENTREZID", keytype = "ENSEMBL")
    annot$entrez1 <- entrez1
    annot$entrez2 <- entrez2
    ENTREZID <- ifelse(!is.na(annot$entrez1), annot$entrez1, annot$entrez2)
    annot$ENTREZID <- ENTREZID
    annot <- annot[c("ENSEMBL","consensus","ENTREZID")]
    names(annot) <- c("ENSEMBL", "SYMBOL", "ENTREZID")
  }
  if(ann=="SYMBOL"){
    annot <- data.frame(SYMBOL = genes, stringsAsFactors = FALSE)
    annot$ENSEMBL <-  mapIds(orgdb, keys=genes, column="ENSEMBL",keytype="SYMBOL")
    entrez1 <- mapIds(orgdb, keys = annot$SYMBOL, column = "ENTREZID", keytype = "SYMBOL")
    entrez2 <- mapIds(orgdb, keys = as.character(annot$ENSEMBL),
                      column = "ENTREZID", keytype = "ENSEMBL")
    annot$entrez1 <- entrez1
    annot$entrez2 <- unlist(unname(lapply(entrez2, function(x){ifelse(is.null(x), NA, x)})))
    ENTREZID <- ifelse(!is.na(annot$entrez1), annot$entrez1, annot$entrez2)
    annot$ENTREZID <- ENTREZID
    annot <- annot[c("ENSEMBL","SYMBOL","ENTREZID")]
  }
  if(dim(genelist)[2]==3){
      annot$logFC <- as.numeric(genelist[ ,2])
      annot$pval <- as.numeric(genelist[ ,3])
  }
  return(as.data.frame(annot) )
}