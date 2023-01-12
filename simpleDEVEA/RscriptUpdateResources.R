getGOlevel <- function(specie){
  require(GO.db)
    bp <- "GO:0008150"
    mf <- "GO:0003674"
    cc <- "GO:0005575"
    getAllBPChildren <- function(goids){
    ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
    ans <- ans[!is.na(ans)]
    return(ans)
    }
    getAllCCChildren <- function(goids){
    ans <- unique(unlist(mget(goids, GOCCCHILDREN), use.names=FALSE))
    ans <- ans[!is.na(ans)]
    return(ans)
    }
    getAllMFChildren <- function(goids){
    ans <- unique(unlist(mget(goids, GOMFCHILDREN), use.names=FALSE))
    ans <- ans[!is.na(ans)]
    return(ans)
    }
    reskk <- list()
    dfBP <- data.frame()
    dfCC <- data.frame()
    dfMF <- data.frame()
    kk <- bp
    for (i in seq_len(100)) {
        reskk[[i]] <- getAllBPChildren(kk)
        if (is.logical(reskk[[i]])) {
            reskk <- reskk[-i]
            break()
        }
        dfBP <- rbind(dfBP, data.frame(id = reskk[[i]], level = i))
        kk <- reskk[[i]]
    }
    kk <- cc
    reskk <- list()
    for (i in seq_len(100)) {
        reskk[[i]] <- getAllCCChildren(kk)
        if (is.logical(reskk[[i]])) {
            reskk <- reskk[-i]
            break()
        }
        dfCC <- rbind(dfCC, data.frame(id = reskk[[i]], level = i))
        kk <- reskk[[i]]
    }
    reskk <- list()
    kk <- mf
    for (i in seq_len(100)) {
        reskk[[i]] <- getAllMFChildren(kk)
        if (is.logical(reskk[[i]])) {
            reskk <- reskk[-i]
            break()
        }
        dfMF <- rbind(dfMF, data.frame(id = reskk[[i]], level = i))
        kk <- reskk[[i]]
    }
    GOlevel <- data.frame(rbind(dfBP,dfCC, dfMF))
    GOlevel <- rbind(
        GOlevel,
        data.frame(
            id = c("GO:0008150","GO:0003674","GO:0005575"),
            level= 0))
    GOlevel <- aggregate(level~id, GOlevel, function(x)x[which.min(abs(x))])
    filePath <- paste0("./resources/",specie,"/GO/GOlevels.Rds")
    saveRDS(GOlevel, filePath)
}


cytoBandCreate <- function(specie = "Mm"){
  library(dplyr)
  #para ideobandas
  specieUCSC <- switch(specie, Hs = "hg38", Mm = "mm10", Rn = "rn6",
                       Ss = "susScr11", Bt = "bosTau9", Mmu = "rheMac10",
                       Dr = "danRer11")
  ruta <- paste0("https://api.genome.ucsc.edu/getData/track?genome=",
                specieUCSC,
                ";track=cytoBandIdeo")
  dat <- curl::curl(url = ruta)
  open(dat)
  out <- jsonlite::read_json(ruta, simplifyVector = TRUE) #readLines(dat)
  datos<- bind_rows(out$cytoBandIdeo)
  #datos <- jsonlite::prettify(out)
  #dfIdeoBand <- jsonlite::fromJSON(datos)
  #dfIdeoBanda <- as.data.frame.list(dfIdeoBand$cytoBandIdeo[[1]])
  dfIdeoBanda <- datos
  dfIdeoBandaLimpio <- dfIdeoBanda %>% dplyr::select(-c(4, 5)) %>%
    group_by(chrom) %>%
    summarise(min = min(chromStart), max = max(chromEnd)) %>%
    dplyr::filter(!grepl("_", chrom)) %>% mutate(chrom = sub("chr", "", chrom))
  dfIdeoSort <- dfIdeoBandaLimpio[gtools::mixedorder(dfIdeoBandaLimpio$chrom), ] %>%
    mutate(chrom = paste0("chr", chrom)) %>% as.data.frame()
  dfRanges <- regioneR::toGRanges(dfIdeoSort)
  saveRDS(dfRanges, paste0("./resources/",specie,"/cytoband/",specie,"_genomicRanges.Rds"))
  saveRDS(dfIdeoBanda, paste0("./resources/",specie,"/cytoband/",specie,"_cytoBand.Rds"))
  # la anotacion aún no está arreglada, de momento por bash
  # curl -L ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.chr.gff3.gz >mm10.gtf.gz
  # zcat mm10.gtf.gz | awk '$3=="gene"{print $1,$4,$5,$9}' | awk 'BEGIN{OFS="\t"}{split($4,a,";");print a[1],$1,$2,$3}' | sed 's/ID=gene://g' >Mm_annot.txt
}


updateDatabases <- function(species){
  require("limma")
    species.KEGG <- NULL
        if (is.null(species.KEGG)) {
        species <- match.arg(species, c("Ag", "At", "Bt", "Ce",
            "Dm", "Dr", "EcK12", "EcSakai", "Gg", "Hs", "Mm",
            "Mmu", "Pf", "Pt", "Rn", "Ss", "Xl"))
        species.KEGG <- switch(species, Ag = "aga", At = "ath",
            Bt = "bta", Ce = "cel", Cf = "cfa", Dm = "dme", Dr = "dre",
            EcK12 = "eco", EcSakai = "ecs", Gg = "gga", Hs = "hsa",
            Mm = "mmu", Mmu = "mcc", Pf = "pfa", Pt = "ptr",
            Rn = "rno", Ss = "ssc", Xl = "xla")
        }
    # Kegg
    GeneID.PathID <- getGeneKEGGLinks(species.KEGG, convert = FALSE)
    filename <- paste0("./resources/",species,"/KEGG/KeggLinks.Rds")
    saveRDS(GeneID.PathID, filename)
    
    PathID.PathName <- getKEGGPathwayNames(species.KEGG,
           remove.qualifier = TRUE)
    filename <- paste0("./resources/",species,"/KEGG/KeggPathwayNames.Rds")
    saveRDS(PathID.PathName, filename)
    # GO
    orgPkg <- paste0("org.", species, ".eg.db")
    obj <- paste0("org.", species, ".egGO2ALLEGS")
    egGO2ALLEGS <- tryCatch(getFromNamespace(obj, orgPkg), error = function(e) FALSE)
    GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[, c("gene_id",
            "go_id", "Ontology")]
    filename <- paste0("./resources/",species,"/GO/GOlinks.Rds")
    saveRDS(GeneID.PathID, filename)
}

# Función para crear dataset para hacer GSEA pathway ##################
buildKeggDataset <- function(specie="Mm"){
    require("dplyr")
    #GeneID.PathID <- getGeneKEGGLinks(specie, convert = FALSE)
    #PathName <- getKEGGPathwayNames(specie, remove.qualifier = TRUE)
    GeneID.PathID <- readRDS(paste0("./resources/",specie,"/KEGG/KeggLinks.Rds"))
    PathName <- readRDS(paste0("./resources/",specie,"/KEGG/KeggPathwayNames.Rds"))
    PathName$Id <- paste(PathName$PathwayID,PathName$Description,sep="_")
    dataSet <- left_join(GeneID.PathID, PathName, by = c("PathwayID"="PathwayID"))
    dataSet$Id <- gsub("path:","",dataSet$Id)
    dataSet <- dataSet[,c(4,1)]
    saveRDS(dataSet,paste0("./resources/",specie,"/GSEA/keggDataGSEA.Rds"))
}


# getGOlevel("Mm")
# updateDatabases("Mm")
# cytoBandCreate("Mm")
# buildKeggDataset("Mm")
# 
# getGOlevel("Hs")
# updateDatabases("Hs")
# cytoBandCreate("Hs")
# buildKeggDataset("Hs")
# 
# getGOlevel("Rn")
# updateDatabases("Rn")
# cytoBandCreate("Rn")
# buildKeggDataset("Rn")

getGOlevel("At")
updateDatabases("At")
cytoBandCreate("At")
buildKeggDataset("At")
