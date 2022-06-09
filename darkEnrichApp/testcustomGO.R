library(tidyverse)
library(limma)
library(DESeq2)
source("utils.R")


datos <- readRDS("./data/deseq.Rds")

colData(datos)@listData <- colData(datos)@listData %>%
    as.data.frame() %>% mutate_at(vars(-sizeFactor, contains('replaceable')), as.character) %>% ##aqui##
    mutate_at(vars(-sizeFactor, contains('replaceable')), as.factor) %>% as.list()
#ressh <- as.data.frame(lfcShrink(datos, coef=(as.numeric(design())+1), type="apeglm", parallel = TRUE))
ressh <- as.data.frame(results(datos, contrast = list(resultsNames(datos)[2] ))) #09/02/2020
ressh <- ressh %>% dplyr::select(-stat) #09/02/2020
conversionids <- geneIdConverter2(rownames(ressh), "Mm" )
padjNAtrue <- which(is.na(ressh$padj)) 
if(length(padjNAtrue)!=0 ){
    conversionRes <- conversion$ids[-padjNAtrue,]
    ressh <- ressh[-padjNAtrue,]
}else{
    conversionRes <- conversionids
}
ressh$baseMean <- round(ressh$baseMean,4)
ressh$lfcSE <- round(ressh$lfcSE,4)
ressh$log2FoldChange <- round(ressh$log2FoldChange,4)
ressh <- cbind(`Description`=conversionRes$description, ressh)
ressh <- cbind(`ENTREZ`=conversionRes$ENTREZID, ressh)
ressh <- cbind(`ENSEMBL` = conversionRes$ENSEMBL, ressh)
#ressh <- cbind(`GeneName_Symbol`=conversion$consensus, ressh)
ressh <- cbind(`GeneName_Symbol`=conversionRes$SYMBOL, ressh) #
ressh <-  ressh %>% dplyr::select(-c(pvalue))
spc = "Mus_musculus"
links = paste0("<a href='http://www.ensembl.org/",spc,"/Gene/Summary?db=core;g=",
               rownames(ressh),"' target='_blank'>",rownames(ressh),"</a>")
ressh <- cbind(`User_GeneId`= links, ressh)

datagenesUp <- getSigUpregulated(ressh, 0.05, 0.5, "Mm" ) 
datagenesDown <- getSigDownregulated(ressh, 0.05, -0.5, "Mm" ) 
datagenesall <- rbind(datagenesUp, datagenesDown)

data <- datagenesall
species = "Mm"
kkcustomGO <- function(data, universe = NULL, species = "Hs", prior.prob = NULL,
                       covariate = NULL, plot = FALSE, coef = 1, FDR = 0.05, golevelFile) {
    if (!is.data.frame(data)) {
        stop("de should be a data.frame with firt column as symbol
        Id and second column as entrez Id.")
    }
    if( dim(data)[2] != 2){
        stop("de should be a data.frame with firt column as symbol
        Id and second column as entrez Id.")
    }
    golevelFile <- paste0("./resources/",species,"/GO/GOlevels.Rds")
    names(data) <- c("SYMBOL","ENTREZID")
    de <- data
    de <- as.character(de$ENTREZID)
    suppressPackageStartupMessages(OK <- requireNamespace("GO.db",quietly = TRUE))
    if (!OK) stop("GO.db package required but not installed (or can't be loaded)")
    suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly = TRUE))
    if (!OK) stop("AnnotationDbi package required but not installed (or can't be loaded)")
    orgPkg <- paste0("org.", species, ".eg.db")
    suppressPackageStartupMessages(OK <- requireNamespace(orgPkg, quietly = TRUE))
    if (!OK)
        stop(orgPkg, " package required but not not installed (or can't be loaded)")
    obj <- paste0("org.", species, ".egGO2ALLEGS")
    egGO2ALLEGS <- tryCatch(getFromNamespace(obj, orgPkg), error = function(e) FALSE)
    if (is.logical(egGO2ALLEGS))
        stop("Can't find gene ontology mappings in package ", orgPkg)
    if (is.null(universe)) {
        #GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[, c("gene_id",
        #  "go_id", "Ontology")]
        file1 <- paste0("./resources/",species,"/GO/GOlinks.Rds")
        GeneID.PathID <- readRDS(file1)
        i <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
        GeneID.PathID <- GeneID.PathID[i, ]
        universe <- unique(GeneID.PathID[, 1])
        prior.prob <- covariate <- NULL
    } else {
        universe <- as.character(universe)
        lu <- length(universe)
        if (!is.null(prior.prob) && length(prior.prob) != lu)
            stop("universe and prior.prob must have same length")
        if (!is.null(covariate) && length(covariate) != lu)
            stop("universe and covariate must have same length")
        if (anyDuplicated(universe)) {
            i <- !duplicated(universe)
            if (!is.null(covariate))
                covariate <- covariate[i]
            if (!is.null(prior.prob))
                prior.prob <- prior.prob[i]
            universe <- universe[i]
        }
        i <- (universe %in% AnnotationDbi::Lkeys(egGO2ALLEGS))
        universe <- universe[i]
        if (!is.null(covariate))
            covariate <- covariate[i]
        if (!is.null(prior.prob))
            prior.prob <- prior.prob[i]
        AnnotationDbi::Lkeys(egGO2ALLEGS) <- universe
        GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[, c("gene_id",
                                                                 "go_id", "Ontology")]
        d <- duplicated(GeneID.PathID[, c("gene_id", "go_id")])
        GeneID.PathID <- GeneID.PathID[!d, ]
    }
    if (is.list(de)) {
        if (is.data.frame(de))
            stop("de should be a list of character vectors. It should not be a data.frame.")
    } else {
        de <- list(DE = de)
    }
    nsets <- length(de)
    if (!all(vapply(de, is.vector, TRUE)))
        stop("components of de should be vectors")
    for (s in 1:nsets) de[[s]] <- unique(as.character(de[[s]]))
    names(de) <- trimWhiteSpace(names(de))
    NAME <- names(de)
    i <- which(NAME == "" | is.na(NAME))
    NAME[i] <- paste0("DE", i)
    names(de) <- makeUnique(NAME)
    NGenes <- length(universe)
    if (NGenes < 1L)
        stop("No annotated genes found in universe")
    for (s in 1:nsets) de[[s]] <- de[[s]][de[[s]] %in% universe]
    i <- GeneID.PathID[, 1] %in% universe
    if (sum(i) == 0L)
        stop("Pathways do not overlap with universe")
    GeneID.PathID <- GeneID.PathID[i, ]
    if (!is.null(covariate)) {
        if (!is.null(prior.prob))
            message("prior.prob being recomputed from covariate")
        covariate <- as.numeric(covariate)
        isDE <- (universe %in% unlist(de))
        o <- order(covariate)
        prior.prob <- covariate
        span <- approx(x = c(20, 200), y = c(1, 0.5), xout = sum(isDE),
                       rule = 2)$y
        prior.prob[o] <- tricubeMovingAverage(isDE[o], span = span)
        if (plot)
            barcodeplot(covariate, index = isDE, worm = TRUE, span.worm = span,
                        main = "DE status vs covariate")
    }
    if (is.null(prior.prob)) {
        X <- matrix(1, nrow(GeneID.PathID), nsets + 1)
        colnames(X) <- c("N", names(de))
    } else {
        names(prior.prob) <- universe
        X <- matrix(1, nrow(GeneID.PathID), nsets + 2)
        X[, nsets + 2] <- prior.prob[GeneID.PathID[, 1]]
        colnames(X) <- c("N", names(de), "PP")
    }
    for (s in 1:nsets) X[, s + 1] <- (GeneID.PathID[, 1] %in% de[[s]])
    S <- rowsum(X, group = GeneID.PathID[, 2], reorder = FALSE)
    PValue <- matrix(0, nrow = nrow(S), ncol = nsets)
    colnames(PValue) <- paste("P", names(de), sep = ".")
    nde <- lengths(de, use.names = FALSE)
    if (!is.null(prior.prob)) {
        SumPP <- sum(prior.prob)
        M2 <- NGenes - S[, "N"]
        Odds <- S[, "PP"]/(SumPP - S[, "PP"]) * M2/S[, "N"]
        if (!requireNamespace("BiasedUrn", quietly = TRUE))
            stop("BiasedUrn package required but is not installed (or can't be loaded)")
        for (j in seq_len(nsets)) for (i in seq_len(nrow(S))) PValue[i,
                                                                     j] <- BiasedUrn::pWNCHypergeo(S[i, 1L + j], S[i, "N"],
                                                                                                   M2[i], nde[j], Odds[i], lower.tail = FALSE) + BiasedUrn::dWNCHypergeo(S[i,
                                                                                                                                                                           1L + j], S[i, "N"], M2[i], nde[j], Odds[i])
        S <- S[, -ncol(S)]
    } else {
        for (j in seq_len(nsets)) PValue[, j] <- phyper(S[, 1L + j] -
                                                            0.5, nde[j], NGenes - nde[j], S[, "N"], lower.tail = FALSE)
    }
    GOID <- rownames(S)
    TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = GOID,
                                                   columns = "TERM"))
    m <- match(GOID, GeneID.PathID[, 2])
    Ont <- GeneID.PathID[m, 3]
    ## aportacion
    kk <- GeneID.PathID  # [ which(GeneID.PathID$gene_id %in% de[[1]]), ]
    kkb <- left_join(data, kk, by = c("ENTREZID" = "gene_id"))
    kkb <- kkb[which(!is.na(kkb$ENTREZID)), ]
    kkk <- kkb %>% group_by(go_id) %>%
        summarize(genes = paste(ENTREZID, collapse = ", "))
    ##
    Results <- data.frame(Term = TERM[, 2], Ont = Ont, S, PValue, stringsAsFactors = FALSE)
    Results$go_id <- rownames(Results)
    resultado <- left_join(Results, kkk, by = c(go_id = "go_id"))
    resultado <- resultado[which(resultado$P.DE < 0.05), ]
    GOlevel <- readRDS(golevelFile)
    GOlevel$id <- as.character(GOlevel$id)
    resultado <- left_join(resultado, GOlevel, by = c("go_id"="id"))
    resultado <- resultado %>% arrange(P.DE)
    return(resultado)
}


which(is.na(goall$Term))
