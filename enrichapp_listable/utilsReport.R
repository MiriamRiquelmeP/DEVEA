# Chorplot ##########################################
legendChorplotReport <- function(enrichdf){
    labels <- enrichdf$Pathway
    colours = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(labels))
    par(bg="#edf0f2", mar=c(0.5,0.5,0.5,0.5))
    # plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    # legend("topleft", legend = labels, pch=16, pt.cex = 1.4,
    #        cex = 0.8, bty="n", col = colours, bg="#edf0f2",
    #        text.col="black", border = "#edf0f2")
    lng <- length(labels)
    plot(x=rep(0, lng), y = seq_len(lng)*(1/lng),
     pch=16, xlim = 0:1, ylim=0:1, yaxt="n", xaxt="n", bty="n", col = colours )
text(x=rep(0.01,lng), y =  seq_len(lng)*(1/lng), labels = labels, col = "black", adj=0 )
}

chordPlotReport <- function(enrichdf, nRows = 10, ont=NULL,  orderby=NULL) {
  if(! "dplyr" %in% .packages()) require("dplyr")
  if(! "tidyr" %in% .packages()) require("tidyr")
  if(! "chordiag" %in% .packages()) require("chorddiag")

  name <- match.arg(names(enrichdf)[1], c("Pathway", "Term"))

  if(dim(enrichdf)[2]==7){
    if(!is.null(ont)){
      enrichdf = enrichdf[enrichdf$Ont==ont, ]
    } else{
      stop("Ontology must be provided if enrichdf is customGO object")
    }
  }

  if(!is.null(orderby)){
    orderby = match.arg(orderby, c("DE", "P.DE", "N", name))
    if(orderby=="P.DE" | orderby == name){
      enrichdf <- enrichdf %>% arrange(get(orderby))
    } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
  }

  kgg <- enrichdf
  n <- nRows

  tmp <- kgg[1:n, c(name, "genes") ]
  tmp2 <- tmp %>% separate_rows(genes, convert = TRUE)
  kk <- tidyr::pivot_wider(tmp2, names_from = name, values_from = genes)

  ns <- NULL
  lenVect <- NULL
  comunTotal <- NULL
  for (i in seq(1, dim(kk)[2])) {
    pt <- unlist(kk[1, i][[1]])
    nd <- NULL
    lenVect <- append(lenVect, length(pt))
    for (j in seq(1, dim(kk)[2])) {
      ns <-
        append(ns, length(intersect(unlist(kk[1, i][[1]]), unlist(kk[1, j][[1]]))))
      nd <- append(ns, intersect(unlist(kk[1, i][[1]]), unlist(kk[1, j][[1]])))
      # if (i != j) {
      #   pt <- pt[!(pt %in% unlist(kk[1, j][[1]]))]
      # }
    }
    comunTotal <- append(comunTotal, length(pt[!pt %in% nd]))
    #nd <- append(nd, length(pt))
  }
  m <- matrix(ns, nrow = dim(kk)[2])
  
  dimnames(m) <- list(have = names(kk), prefer = names(kk))

  diag(m) <- (as.vector(lenVect) - comunTotal )
  
  p <- chorddiag(
    m,
    type = "directional",
    groupnamePadding = 20,
    margin = 100,
    groupColors = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n),
    showGroupnames = FALSE,
    tooltipFontsize = 12,
    showTicks = FALSE
  )
  return(p)
}

# plot para exportar cnet a cytoscape ###############
customCnet2Cytoscape <- function(kgg, category=NULL, nPath=NULL, byDE=FALSE){
    if(! "ggraph" %in% .packages()) require("ggraph")
    if(! "igraph" %in% .packages()) require("igraph")
    if(! "dplyr" %in% .packages()) require("dplyr")
    if(! "tidyr" %in% .packages()) require("tidyr")
    ## check if cytoscape is running
    a <- tryCatch( expr= { cytoscapePing() },
                   error = function(e){print("Error!! Cytoscape must be running")}
                   )
    ## define palette
    color_palette <- function(colors) colorRampPalette(colors)(n = 299)
    if(byDE){
        kgg <- kgg %>% arrange(-DE)
    }
    ## check parameters
    if(!is.null(category)){
        if(is.numeric(category)){
            tmp <- kgg[category, c(1,6)]
        } else{
            tmp <- kgg[kgg$Pathway %in% category, c(1,6)]
        }
    } else{
        tmp <- kgg[,c(1,6)]
    }
    if(!is.null(nPath)){
        if(is.numeric(nPath)){
            tmp <- tmp[1:nPath, ]
        } else{ stop("nPath must be numeric or NULL")}
    }
    # prepare data
    tmp2 <- tmp %>%  separate_rows( genes, convert=TRUE)
    g <- igraph::graph.data.frame(tmp2, directed=FALSE)
    size <- tmp2 %>% dplyr::select(Pathway) %>% group_by(Pathway) %>%  summarise(n=n())
    size <- left_join(tmp, size, by=c("Pathway"))
    size <- size$n
    V(g)$size <- min(size)/2
    n <- dim(tmp)[1]
    V(g)$size[1:n] <- size
    V(g)$pval <- NA
    V(g)$pval[1:n] <- kgg$p-value[1:n]
    edge_layer <- geom_edge_link(alpha=.8, colour='darkgrey')
    fc <- V(g)$pval
    V(g)$color <- fc
    createNetworkFromIgraph(g, "customIgraph")
}

# Plot para plotear cnet para kegg ###############
customCnetKeggReport <- function(kgg, category=NULL, nPath=NULL, byDE=FALSE, nr, genesUp, genesDown){
    if(! "ggraph" %in% .packages()) require("ggraph")
    if(! "igraph" %in% .packages()) require("igraph")
    if(! "dplyr" %in% .packages()) require("dplyr")
    if(! "tidyr" %in% .packages()) require("tidyr")
    kgg <- kgg[nr,]
    color_palette <- function(colors) colorRampPalette(colors)(n = 299)
    if(byDE){
        kgg <- kgg %>% arrange(-DE)
    }
    if(!is.null(category)){
        if(is.numeric(category)){
            tmp <- kgg[category,]
        } else{
            tmp <- kgg[kgg$Pathway %in% category,]
        }
    } else{
        tmp <- kgg
    }
    if(!is.null(nPath)){
        if(is.numeric(nPath)){
            tmp <- tmp[1:nPath, ]
        } else{ stop("nPath must be numeric or NULL")}
    }
    genesUp$dir <- "#ffa200"
    genesDown$dir <- "#91ebff"
    genesAll <- rbind(genesUp,genesDown)
    pval <- tmp[,c(1,4)]
    tmp <- tmp[,c(1,6)]
    tmp2 <- tmp %>%  separate_rows( genes, convert=TRUE)
    g <- igraph::graph.data.frame(tmp2, directed=FALSE)
    size <- tmp2 %>% dplyr::select(Pathway) %>% group_by(Pathway) %>%  summarise(n=n())
    size <- left_join(tmp, size, by=c("Pathway"))
    pval <- left_join(tmp, pval, by=c("Pathway"))
    size <- size$n
    pval <- pval$P.DE
    V(g)$size <- min(size)/2
    n <- dim(tmp)[1]
    V(g)$size[1:n] <- size
    V(g)$pval <- NA
    V(g)$pval[1:n] <- pval
    edge_layer <- geom_edge_link(alpha=.8, colour='darkgrey')
    fc <- V(g)$pval
    V(g)$color <- fc
    palette <- color_palette(c("red", "blue"))
    p <- ggraph(g, layout="stress", circular=FALSE) +
        edge_layer +
        geom_node_point( aes_(color=~pval, size=~size) ) +
        scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
        theme_void()
    ##
    ptmp <- p$data
    ptmp2 <- left_join(ptmp, genesAll, by = c("name"="ENTREZID"))
    genesColor <- ptmp2$dir[!is.na(ptmp2$dir)]
    ##
    p <- p + geom_node_text(aes_(label=~name), size=3, data = p$data[1:n,]) +
        scale_color_gradientn(name = "pval", colors=palette, na.value = genesColor )
    return(p)
}


# Plot para plotear cnet para GO ###############
customCnetGoReport <- function(gos, category=NULL, nTerm=NULL, byDE=FALSE, ont="BP"){
    if(! "ggraph" %in% .packages()) require("ggraph")
    if(! "igraph" %in% .packages()) require("igraph")
    if(! "dplyr" %in% .packages()) require("dplyr")
    if(! "tidyr" %in% .packages()) require("tidyr")
    gos <- gos[gos$Ont == ont, ]
    color_palette <- function(colors) colorRampPalette(colors)(n = 299)
    if(byDE){
        gos <- gos %>% arrange(-DE)
    }
    if(!is.null(category)){
        if(is.numeric(category)){
            tmp <- gos[category,]
        } else{
            tmp <- gos[gos$Pathway %in% category,]
        }
    } else{
        tmp <- gos
    }
    if(!is.null(nTerm)){
        if(is.numeric(nTerm)){
            tmp <- tmp[1:nTerm, ]
        } else{ stop("nTerm must be numeric or NULL")}
    }
    pval <- tmp[,c(1,5)]
    tmp <- tmp[,c(1,7)]
    tmp2 <- tmp %>%  separate_rows( genes, convert=TRUE)
    g <- igraph::graph.data.frame(tmp2, directed=FALSE)
    size <- tmp2 %>% dplyr::select(Term) %>% group_by(Term) %>%  summarise(n=n())
    size <- left_join(tmp, size, by=c("Term"))
    pval <- left_join(tmp, pval, by=c("Term"))
    size <- size$n
    pval <- pval$P.DE
    V(g)$size <- min(size)/2
    n <- dim(tmp)[1]
    V(g)$size[1:n] <- size
    V(g)$pval <- NA
    V(g)$pval[1:n] <- pval
    edge_layer <- geom_edge_link(alpha=.8, colour='darkgrey')
    fc <- V(g)$pval
    V(g)$color <- fc
    palette <- color_palette(c("red", "blue"))
    p <- ggraph(g, layout="stress", circular=FALSE) +
        edge_layer +
        geom_node_point( aes_(color=~pval, size=~size) ) +
        scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
        theme_void()
    p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,]) +
        scale_color_gradientn(name = "pval", colors=palette, na.value = "#E5C494")
    return(p)
}

# Función para hacer enrich GO ################
customGOReport <- function(data, universe = NULL, species = "Hs", prior.prob = NULL,
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

# Función para hacer enrich kegg ################
customKeggReport <- function(data, universe = NULL, restrict.universe = FALSE,
    species = "Hs", species.KEGG = NULL, convert = FALSE, gene.pathway = NULL,
    pathway.names = NULL, prior.prob = NULL, covariate = NULL,
    plot = FALSE, ...) {
    if( !is.data.frame(data)){
        stop("de should be a data.frame with firt column as symbol Id and second
             column as entrez Id.")
    }
    if( dim(data)[2] != 2){
        stop("de should be a data.frame with firt column as symbol
        Id and second column as entrez Id.")
    }
    names(data) <- c("SYMBOL","ENTREZID")
    de <- data
    de <- as.vector(de[,2])
    de <- list(DE = de)
    nsets <- length(de)
    if (!all(vapply(de, is.vector, TRUE))) {
        stop("components of de should be vectors")
    }
    for (s in 1:nsets) de[[s]] <- unique(as.character(de[[s]]))
    names(de) <- trimWhiteSpace(names(de))
    NAME <- names(de)
    i <- which(NAME == "" | is.na(NAME))
    NAME[i] <- paste0("DE", i)
    names(de) <- makeUnique(NAME)
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
    if (is.null(gene.pathway)) {
        #GeneID.PathID <- getGeneKEGGLinks(species.KEGG, convert = convert)
         file1 <- paste0("./resources/",species,"/KEGG/KeggLinks.Rds")
         GeneID.PathID <- readRDS(file1)
    } else {
        GeneID.PathID <- gene.pathway
        d <- dim(GeneID.PathID)
        if (is.null(d)) {
            stop("gene.pathway must be data.frame or matrix")
        }
        if (d[2] < 2) {
            stop("gene.pathway must have at least 2 columns")
        }
        isna <- rowSums(is.na(GeneID.PathID[, 1:2])) > 0.5
        GeneID.PathID <- GeneID.PathID[!isna, ]
        ID.ID <- paste(GeneID.PathID[, 1], GeneID.PathID[, 2],
            sep = ".")
        if (anyDuplicated(ID.ID)) {
            d <- duplicated(ID.ID)
            GeneID.PathID <- GeneID.PathID[!d, ]
        }
    }
    if (is.null(pathway.names)) {
        #PathID.PathName <- getKEGGPathwayNames(species.KEGG,
         #   remove.qualifier = TRUE)
        file2 <- paste0("./resources/",species,"/KEGG/KeggPathwayNames.Rds")
        PathID.PathName <- readRDS(file2)
    } else {
        PathID.PathName <- pathway.names
        d <- dim(PathID.PathName)
        if (is.null(d)) {
            stop("pathway.names must be data.frame or matrix")
        }
        if (d[2] < 2) {
            stop("pathway.names must have at least 2 columns")
        }
        isna <- rowSums(is.na(PathID.PathName[, 1:2])) > 0.5
        PathID.PathName <- PathID.PathName[!isna, ]
    }
    if (is.null(universe)) {
        universe <- unique(GeneID.PathID[, 1])
        prior.prob <- covariate <- NULL
    } else {
        universe <- as.character(universe)
        lu <- length(universe)
        if (!lu) {
            stop("No genes in universe")
        }
        if (!is.null(prior.prob) && length(prior.prob) != lu) {
            stop("universe and prior.prob must have same length")
        }
        if (!is.null(covariate) && length(covariate) != lu) {
            stop("universe and covariate must have same length")
        }
        if (restrict.universe) {
            i <- universe %in% GeneID.PathID[, 1]
            universe <- universe[i]
            if (!is.null(prior.prob)) {
                prior.prob <- prior.prob[i]
            }
            if (!is.null(covariate)) {
                covariate <- covariate[i]
            }
        }
    }
    if (anyDuplicated(universe)) {
        d <- duplicated(universe)
        if (!is.null(covariate)) {
            covariate <- rowsum(covariate, group = universe,
                reorder = FALSE)
            n <- rowsum(rep_len(1L, length(universe)), group = universe,
                reorder = FALSE)
            covariate <- covariate/n
        }
        if (!is.null(prior.prob)) {
            prior.prob <- rowsum(prior.prob, group = universe,
                reorder = FALSE)
            n <- rowsum(rep_len(1L, length(universe)), group = universe,
                reorder = FALSE)
            prior.prob <- prior.prob/n
        }
        universe <- universe[!d]
    }
    NGenes <- length(universe)
    if (NGenes < 1L) {
        stop("No annotated genes found in universe")
    }
    for (s in 1:nsets) de[[s]] <- de[[s]][de[[s]] %in% universe]
    i <- GeneID.PathID[, 1] %in% universe
    if (sum(i) == 0L) {
        stop("Pathways do not overlap with universe")
    }
    GeneID.PathID <- GeneID.PathID[i, ]
    if (!is.null(covariate)) {
        if (!is.null(prior.prob)) {
            message("prior.prob being recomputed from covariate")
        }
        covariate <- as.numeric(covariate)
        isDE <- (universe %in% unlist(de))
        o <- order(covariate)
        prior.prob <- covariate
        span <- approx(x = c(20, 200), y = c(1, 0.5), xout = sum(isDE),
            rule = 2)$y
        prior.prob[o] <- tricubeMovingAverage(isDE[o], span = span)
        if (plot) {
            barcodeplot(covariate, index = isDE, worm = TRUE,
                span.worm = span, main = "DE status vs covariate")
        }
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
    for (s in 1:nsets) {
        X[, s + 1] <- (GeneID.PathID[, 1] %in% de[[s]])
    }
    S <- rowsum(X, group = GeneID.PathID[, 2], reorder = FALSE)
    PValue <- matrix(0, nrow = nrow(S), ncol = nsets)
    colnames(PValue) <- paste("P", names(de), sep = ".")
    nde <- lengths(de, use.names = FALSE)
    if (!is.null(prior.prob)) {
        SumPP <- sum(prior.prob)
        M2 <- NGenes - S[, "N"]
        Odds <- S[, "PP"]/(SumPP - S[, "PP"]) * M2/S[, "N"]
        if (!requireNamespace("BiasedUrn", quietly = TRUE)) {
            stop("BiasedUrn package required but is not installed (or can't be loaded)")
        }
        for (j in seq_len(nsets)) {
            for (i in seq_len(nrow(S))) {
                PValue[i, j] <- BiasedUrn::pWNCHypergeo(S[i,
                  1L + j], S[i, "N"], M2[i], nde[j], Odds[i],
                  lower.tail = FALSE) + BiasedUrn::dWNCHypergeo(S[i,
                  1L + j], S[i, "N"], M2[i], nde[j], Odds[i])
            }
        }
        S <- S[, -ncol(S)]
    } else {
        for (j in seq_len(nsets)) {
            PValue[, j] <- phyper(S[, 1L + j] - 0.5, nde[j],
                NGenes - nde[j], S[, "N"], lower.tail = FALSE)
        }
    }

    g <- rownames(S)
    m <- match(g, PathID.PathName[, 1])
    path <- PathID.PathName[m, 2]
    ##
    kk <- GeneID.PathID
    data <- data.frame(ENTREZID = data[,2])
    data$ENTREZID <- as.character(data$ENTREZID)
    kkb <- inner_join(kk, data, by = c(GeneID = "ENTREZID"))
    kkb <- kkb[!duplicated(kkb), ]
    # kkb <- kkb[ which( !is.na( kkb$SYMBOL ) ), ]
    kkk <- kkb %>% group_by(PathwayID) %>% summarize(genes = paste(GeneID,
        collapse = ","))
    ##

    Results <- data.frame(Pathway = path, S, PValue, stringsAsFactors = FALSE)
    rownames(Results) <- g
    Results$pathID <- rownames(Results)
    resultado <- left_join(Results, kkk, by = c(pathID = "PathwayID"))
    resultado <- resultado[which(resultado$P.DE < 0.05), ]
    resultado <- resultado %>% arrange(P.DE)
    return(resultado)
}

# Función para crear tablas  con desplegable de genes ############
datatable2Report <- function(x, vars = NULL, opts = NULL, ...) {
  names_x <- names(x)
  if (is.null(vars))
    stop("'vars' must be specified!")
  pos <- match(vars, names_x)
  if (any(map_chr(x[, pos], typeof) == "list"))
    stop("list columns are not supported in datatable2()")
  pos <- pos[pos <= ncol(x)] + 1
  rownames(x) <- NULL
  if (nrow(x) > 0)
    x <- cbind(` ` = "&oplus;", x)
  if(dim(x)[2]==7){
      js <- c("function(row, data) {",
                "if (data[5]>1000 | data[5]<1){",
                "$('td:eq('+(4)+')', row).html(data[5].toExponential(3));",
                "}",
                "}")
  }else if(dim(x)[2]==9){
      js <- c("function(row, data) {",
                "if (data[6]>1000 | data[6]<1){",
                "$('td:eq('+(5)+')', row).html(data[6].toExponential(3));",
                "}",
                "}")
  }
  opts <- c(opts,
            list(
              columnDefs = list(list(visible = FALSE, targets = c(0,pos)),
                                list(orderable = FALSE,
                                     className = "details-control", targets = 1),
                                list(className = "dt-left", targets = 1:2),
                                list(className = "dt-right",targets = 3:ncol(x) 
                                     )
                                ),
              lengthMenu = list(c(10,25,50,100,-1), c(10,25,50,100,"All")),
              rowCallback = JS(js),
              dom = "Bfrtipl"
              ) )
  datatable(x, ..., extensions = "Buttons", options = opts, 
            callback = JS(.callback2(x = x, pos = c(0, pos) ) ) )
}

# funcion auxiliar de datatable2 ##############################
.callback2 <- function(x, pos = NULL) {
    part1 <- "table.column(1).nodes().to$().css({cursor: 'pointer'});"
    part2 <- .child_row_table2(x, pos = pos)
    part3 <- "
   table.on('click', 'td.details-control', function() {
    var td = $(this), row = table.row(td.closest('tr'));
    if (row.child.isShown()) {
      row.child.hide();
      td.html('&oplus;');
    } else {
      row.child(format(row.data())).show();
      td.html('&ominus;');
    }
  });"
    paste(part1, part2, part3)
}

# funcion auxiliar de datatable2 ##############################
.child_row_table2 <- function(x, pos = NULL) {
    names_x <- paste0(names(x), ":")
    text <- "var format = function(d) {
    text = '<div><table >' +"
    for (i in seq_along(pos)) {
        text <- paste(text, glue::glue("'<tr>' +
          '<td>' + '{names_x[pos[i]]}' + '</td>' +
          '<td>' + d[{pos[i]}] + '</td>' +
        '</tr>' + "))
    }
    paste0(text, "'</table></div>'
      return text;};")
}

# funcion que preparar los datos de enrich go para pasárlos a datatable2 ###############
go2DTReport <- function(enrichdf, data, orderby = NULL, nrows = NULL) {
    if(!is.data.frame(enrichdf) | !is.data.frame(data)){
        stop("enrichdf and data should be data.frame")
    }
    names(data) <- c("SYMBOL", "ENTREZID")
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DE", "P.DE", "N", "Term"))
        if(orderby=="P.DE" | orderby =="Term"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    enrichdf2 <- enrichdf %>%
        mutate(url = paste0("<a href='https://www.ebi.ac.uk/QuickGO/term/",
                            go_id,"' target='_blank'>",go_id,"</a>"))
    CAup <- enrichdf2[, c(1, 2, 3, 4, 5, 7, 8, 9)]
    CAup$genes <- gsub(",", ", ", CAup$genes)
    for (i in seq(1:length(CAup$genes))) {
       mg <- as.numeric(unlist(strsplit(CAup$genes[i], ", ")))
       mg2 <- match(mg, data$ENTREZID)
       CAup$genes[i] <- paste0(data[mg2, 1], collapse = ", ")
    }
    splitGenes <- strsplit(CAup$genes, ", ")
    CAup$genes <- lapply(
        splitGenes, function(x){
            paste(sort(x), collapse = ", ")
            }
        )
    if(!is.null(nrows) & is.numeric(nrows)){
        CAup <- CAup[1:nrows, ]
    }
    #CAup <- CAup %>% mutate(p-value = format(p-value, scientific = T, digits = 4))
    return(CAup)
}

# Recupera todos los ids de GO y el nivel al que pertenecen
# Ejemplo de uso:
# GOlevel = getGOlevel()
# Guardar lo que genera en resources/GOlevel.Rds #############
getGOlevelReport <- function(specie){
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
    filePath <- paste0("./resources/",specie,"/GO/GOlevel.Rds")
    saveRDS(GOlevel, filePath <- paste0("./resources/",specie,"/GO/GOlevels.Rds"))
}

# funcion que preparar los datos de enrich kegg para pasárlos a datatable2 ###############
kegg2DTReport <- function(enrichdf, data, orderby = NULL, nrows = NULL) {
    if(!is.data.frame(enrichdf) | !is.data.frame(data)){
        stop("enrichdf and data should be data.frame")
    }
    names(data) <- c("SYMBOL", "ENTREZID")
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DE", "P.DE", "N", "Pathway"))
        if(orderby=="P.DE" | orderby =="Pathway"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    enrichdf2 <- enrichdf %>%
        mutate(pathway = limma::strsplit2(enrichdf$pathID, ":")[, 2],
               genesURL = gsub(",", "+", genes)) %>%
        mutate(url = paste0("<a href='https://www.kegg.jp/kegg-bin/show_pathway?",
        pathway, "/", genesURL,"' target='_blank'>", pathway, "</a>"))
    CAup <- enrichdf2[, c(1, 2, 3, 4, 6, 9)]
    CAup$genes <- gsub(",", ", ", CAup$genes)
    for (i in seq(1:length(CAup$genes))) {
        mg <- as.numeric(unlist(strsplit(CAup$genes[i], ", ")))
        mg2 <- match(mg, data$ENTREZID)
        CAup$genes[i] <- paste0(data[mg2, 1], collapse = ", ")
    }
    splitGenes <- strsplit(CAup$genes, ", ")
    CAup$genes <- lapply(
        splitGenes, function(x){
            paste(sort(x), collapse = ", ")
            }
        )
    if(!is.null(nrows) & is.numeric(nrows)){
        CAup <- CAup[1:nrows, ]
    }
    #CAup <- CAup %>% mutate(p-value = format(p-value, scientific = T, digits = 4))
    return(CAup)
}

# Plot barras de GO ####################
plotGOReport <- function(enrichdf, nrows = 30, orderby="p-val", ont, colors=NULL){
    require(plotly)
    if(!is.data.frame(enrichdf)){
        stop("enrichdf should be data.frame")
    }
    if(!exists("ont") | !(ont %in% c("BP","MF","CC")) ){
        stop("A valid value should be provided for 'ont'")
    }
    names(enrichdf) <- gsub("P.DE", "p-val", names(enrichdf) )
    names(enrichdf) <- gsub("DE", "DEG", names(enrichdf) )
    dataTitle <- list(BP=c("Biological Process", colors),
                      MF=c("Molecular Function", colors),
                      CC=c("Cellular Component", colors))
    enrichdf <- enrichdf[enrichdf$Ont == ont, ]
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DEG", "p-val", "N", "Term"))
        if(orderby=="DEG" | orderby =="Term"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    # p <- enrichdf[1:nrows,] %>%
    #     plot_ly(x=~DEG, y=~go_id, text=~Term, type = "bar",
    #             marker = list(color=dataTitle[[ont]][2]),
    #             orientation = "v",
    #             hovertext = paste0(enrichdf$Pathway,"\np-val: ",format(enrichdf$`p-val`, scientific = T, digits = 4))) %>%
    #     layout(margin = list(l=100), yaxis = list(title=""),
    #            title=dataTitle[[ont]][1], xaxis = list(tickvals = c(1:max(enrichdf$DEG) ) ))
    p <- enrichdf[1:nrows,] %>% 
      ggplot( aes( y = DEG, x = go_id, 
                   text = paste0("p-val: ",format(`p-val`, scientific = T, digits = 4)) )) +
              geom_bar(position = "stack", stat = "identity", fill = colors) + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = "red") +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    p <- p %>% ggplotly(tooltip = "all")
    return(p)
}
# Plot barras de GOAll ####################
plotGOAllReport <- function(enrichdf, nrows = 30, orderby="p-val",
                      ont, genesUp = NULL, genesDown = NULL, colors = NULL){
    require(plotly)
    require(ggplot2)
    if(!is.data.frame(enrichdf)){
        stop("enrichdf should be data.frame")
    }
    if(!exists("ont") | !(ont %in% c("BP","MF","CC")) ){
        stop("A valid value should be provided for 'ont'")
    }
    names(enrichdf) <- gsub("P.DE", "p-val", names(enrichdf) )
    names(enrichdf) <- gsub("DE", "DEG", names(enrichdf) )
    dataTitle <- list(BP=c("Biological Process", colors ),
                      MF=c("Molecular Function",colors ),
                      CC=c("Cellular Component",colors  ))
    enrichdf <- enrichdf[enrichdf$Ont == ont, ]
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DEG", "p-val", "N", "Term"))
        if(orderby=="p-val" | orderby =="Term"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    enrichAll <- enrichdf
    enrichAll <- enrichAll %>% rowwise() %>%
        mutate(numUp = length(which(strsplit(genes,", ")[[1]] %in% genesUp$ENTREZID ))) %>% 
        mutate(numDown = length(which(strsplit(genes,", ")[[1]] %in% genesDown$ENTREZID ))) %>% 
        mutate(numDownNeg = -length(which(strsplit(genes,", ")[[1]] %in% genesDown$ENTREZID )))
    enrichAll <- enrichAll %>% dplyr::select(go_id, numUp, numDown, numDownNeg, `p-val`)
    df <- data.frame(
        Regulation = rep(c("Up", "Down"), each = nrows),
        goId = enrichAll$go_id,
        x = rep(1:nrows, 2),
        DEG = c(enrichAll$numUp, enrichAll$numDown),
        DENeg = c(enrichAll$numUp, enrichAll$numDownNeg),
        `p_val` = rep(enrichAll$`p-val`, 2)
    )
    colorfill <- c(dataTitle[[ont]][2:3])
    r <- ggplot(df, aes(x = goId, y = DENeg, fill = Regulation,
                        text =paste0("p-val: ",format(p_val, scientific = T, digits = 4)) )) +
        geom_bar(stat = "identity", position = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    r <- r %>% plotly::ggplotly(tooltip = "all")
    p <- ggplot(df, aes(fill = Regulation, y = DEG, x = goId,
                        text =paste0("p-val: ",format(p_val, scientific = T, digits = 4)) )) +
        geom_bar(position = "dodge", stat = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    p <- p %>% ggplotly(tooltip = "all")
    q <- ggplot(df, aes(fill = Regulation, y = DEG, x = goId,
                        text =paste0("p-val: ",format(p_val, scientific = T, digits = 4)) )) +
        geom_bar(position = "stack", stat = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    q <- q %>% ggplotly(tooltip = "all")
    return(list(p,q,r) ) 
}

# Plot barras de Kegg ###########################
plotKeggReport <- function(enrichdf, nrows = 30, orderby="p-val", colors = NULL){
    require(plotly)
    if(!is.data.frame(enrichdf)){
        stop("enrichdf should be data.frame")
    }
    names(enrichdf) <- gsub("P.DE", "p-val", names(enrichdf) )
    names(enrichdf) <- gsub("DE", "DEG", names(enrichdf) )
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DEG", "p-val", "N", "Pathway"))
        if(orderby=="p-val" | orderby =="Pathway"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    enrichdf <- enrichdf[1:nrows,]
    p <- enrichdf %>%
        plot_ly(x=~DEG, y=~pathID, text=~Pathway, type = "bar",
                marker = list(color=colors),
                orientation = "v",
                hovertext = paste0(enrichdf$Pathway,"\np-val: ",format(enrichdf$`p-val`, scientific = T, digits = 4))) %>%
        layout(margin = list(l=100), yaxis = list(title=""),
               title="Kegg pathways", xaxis = list(tickvals = c(1:max(enrichdf$DEG) ) ) )
        
    return(p)
}

# Plot barras de KeggALL ###################
plotKeggAllReport <- function(enrichdf, nrows = 10, orderby = "p-val", 
                        genesUp = NULL, genesDown = NULL, colors = NULL){
    require(plotly)
    require(ggplot2)
        if(!is.data.frame(enrichdf)){
        stop("enrichdf should be data.frame")
        }
    names(enrichdf) <- gsub("P.DE", "p-val", names(enrichdf) )
    names(enrichdf) <- gsub("DE", "DEG", names(enrichdf) )
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DEG", "p-val", "N", "Pathway"))
        if(orderby=="p-val" | orderby =="Pathway"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    enrichAll <- enrichdf[1:nrows, ]
    enrichAll <- enrichAll %>% rowwise() %>%
        mutate(numUp = length(which(strsplit(genes,",")[[1]] %in% genesUp$ENTREZID ))) %>% 
        mutate(numDown = length(which(strsplit(genes,",")[[1]] %in% genesDown$ENTREZID ))) %>% 
        mutate(numDownNeg = -length(which(strsplit(genes,",")[[1]] %in% genesDown$ENTREZID )))
    enrichAll <- enrichAll %>% dplyr::select(pathID, numUp, numDown, numDownNeg, `p-val`)
    enrichAll$pathID <- sub(  "path:", "", enrichAll$pathID )
    df <- data.frame(
        Regulation = rep(c("Up", "Down"), each = nrows),
        pathId = enrichAll$pathID,
        x = rep(1:nrows, 2),
        DEG = c(enrichAll$numUp, enrichAll$numDown),
        DENeg = c(enrichAll$numUp, enrichAll$numDownNeg),
        `p_val` = rep(enrichAll$`p-val`, 2)
    )
    colorfill <- colors #c("#eb3f3f","#3e90bd")
    r <- ggplot(df, aes(x = pathId, y = DENeg, fill = Regulation,
                        text =paste0("p-val: ",format(p_val, scientific = T, digits = 4)) )) +
        geom_bar(stat = "identity", position = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    r <- r %>% plotly::ggplotly(tooltip = "all" )
    p <- ggplot(df, aes(fill = Regulation, y = DEG, x = pathId,
                        text =paste0("p-val: ",format(p_val, scientific = T, digits = 4)) )) +
        geom_bar(position = "dodge", stat = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    p <- p %>% ggplotly(tooltip = "all")
    q <- ggplot(df, aes(fill = Regulation, y = DEG, x = pathId, 
                        text =paste0("p-val: ",format(p_val, scientific = T, digits = 4)) )) +
        geom_bar(position = "stack", stat = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    q <- q %>% ggplotly(tooltip = "all")

    return(list(p, q, r))
}

# Función sin uso actualmente -creo- #############
loadGenes <- function(filegenes){
  load(filegenes)
  auxgenes <- genes
}

# PCA de un objeto DESeq #####################

    plotPCAReport <- function(object, intgroup = "condition", ntop = 500,
                   returnData = TRUE, labels = NULL, customColor = NULL){
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select, ]))
  # the contribution to   the total variance for each component
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <-
    as.data.frame(colData(object)[, intgroup, drop = FALSE])
  # add the intgroup factors together to create a new grouping factor
  # group <- if (length(intgroup) > 1) {
  #     factor(apply(intgroup.df, 1, paste, collapse = ":"))
  # } else {
  #     colData(object)[[intgroup]]
  # }
  if(length(intgroup)>1){
    colgroup <- factor(intgroup.df[ ,intgroup[1] ] )
    shapegroup <- factor(intgroup.df[ ,intgroup[2] ] )
    d <-
      data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = colgroup,
        shape = shapegroup,
        intgroup.df,
        name = colnames(object),
        labels = colData(object)[[labels]]
      )
  } else{
    colgroup <- factor(intgroup.df[ ,intgroup[1] ] )
    d <-
      data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = colgroup,
        intgroup.df,
        name = colnames(object),
        labels = colData(object)[[labels]]
      )
  }
  d$group <- as.factor(d$group)
  # assembly the data for the plot
  if(is.null(customColor)){
    getPalette <- colorRampPalette(c("#f7837b","#1cc3c8"))
    colours <- getPalette(length(levels(d$group)))
  }else{
        colours <- customColor
    }
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    #return(d)
  }
 if(length(intgroup)>1){
    p <- ggplot(data = d,
                aes_string(x = "PC1", y = "PC2", color = "group", shape = "shape")) +
      geom_point(size = 3) +
      ggtitle("PCA for top 500 genes on normalized data") +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      scale_color_manual(values = colours, name = intgroup[1]) +
      scale_shape_manual(values = seq_len(length(d$shape)), name=intgroup[2] )+
      #coord_fixed() +
      ggrepel::geom_text_repel(aes(label = labels, #paste("",d$name, sep = ""),
                                   size = "tam"),
                               show.legend = FALSE, size=3, nudge_y = 0.1) +
      # scale_size_manual("tam", c(1)) +
      theme(text = element_text(size=20))
    }
  else{
    p <- ggplot(data = d,
                aes_string(x = "PC1", y = "PC2", color = "group")) +
      geom_point(size = 3) +
      ggtitle("PCA for top 500 genes on normalized data") +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      scale_color_manual(values = colours, name = intgroup[1]) +
      #coord_fixed() +
      ggrepel::geom_text_repel(aes(label = labels,# paste("",d$name, sep = ""),
                                   size = "tam"),
                               show.legend = FALSE, nudge_y = 0.1) +
      # scale_size_manual(labels = c("tam"), values = c(1)) +
      theme(text = element_text(size=20))
  }
    return(p)
}

pca3dplotReport <- function(object, intgroup = "condition", ntop = 500,
                   returnData = TRUE){
    rv <- rowVars(assay(object))
   select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  
    colgroup <- factor(intgroup.df[ ,intgroup[1] ] )
    d <-data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        PC3 = pca$x[, 3],
        group = colgroup,
        intgroup.df,
        name = colnames(object),
        labels = colData(object)[[intgroup[1]]]
      )

  # assembly the data for the plot
  
  getPalette <- colorRampPalette(c("#f7837b","#1cc3c8"))
  colours <- getPalette(length(levels(d$group)))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    #return(d)
  }
    return(d)
  }


# Función para recuperar los genes up de un objeto DEseq #############
# actualmente para p-val <= 0.05 fijo.
getSigUpregulated <- function(dds, pval=0.05, logfc=0, specie="Mm"){
  #res.sh <- lfcShrink(dds, coef=2, type="apeglm", res = results(dds))
  rk <- as.data.frame(dds)
  rk <- rk[rk$log2FoldChange >logfc & rk$padj<=pval,]
  rk <- rk[ order(rk$padj, decreasing = TRUE), ]
  annot <- geneIdConverter(rownames(rk), specie)
  return(data.frame(SYMBOL = annot$consensus, ENTREZID = annot$ENTREZID, stringsAsFactors = F) )
}

# Función para recuperar los genes down de un objeto DEseq #############
# actualmente para p-val <= 0.05 fijo.
getSigDownregulated <- function(dds, pval=0.05, logfc=0, specie="Mm"){
  #res.sh <- lfcShrink(dds, coef=2, type="apeglm", res = results(dds))
  rk <- as.data.frame(dds)
  rk <- rk[rk$log2FoldChange <logfc & rk$padj<=pval,]
  rk <- rk[ order(rk$padj, decreasing = TRUE), ]
  annot <- geneIdConverter(rownames(rk), specie)
  return(data.frame(SYMBOL = annot$consensus, ENTREZID = annot$ENTREZID, stringsAsFactors = F) )
}

# Convertidor de nombres de genes ###################
# Se le pasa un vector de ensembl y devuelve un df con varios nombres
geneIdConverter <- function(genes, specie="Mm"){ # genes = vector of ensembl gene ids (sólo para Mm por ahora)
  require("EnsDb.Mmusculus.v79")
  require("org.Mm.eg.db")
  require("EnsDb.Hsapiens.v86")
  require("org.Hs.eg.db")
  require("EnsDb.Rnorvegicus.v79")
  require("org.Rn.eg.db") 
  
  if(specie=="Mm"){
      ensdb <- EnsDb.Mmusculus.v79
      orgdb <- org.Mm.eg.db
  }
    else{
        ensdb <- EnsDb.Hsapiens.v86
        orgdb <- org.Hs.eg.db
    }
  annot <- NULL
  annot$ENSEMBL <- genes
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
  return(annot)
}

# Dotplot de objeto enrich kegg ##########################
dotPlotkeggReport <- function(data, n = 20){
  names(data) <- gsub("P.DE", "p-val", names(data) )
  names(data) <- gsub("DE", "DEG", names(data) )
  data$ratio <- data$DEG/data$N
  data <- data[order(data$ratio, decreasing = F), ]
  data <- data[seq_len(n),]
  data$Pathway <- factor(data$Pathway, levels = data$Pathway)
  p <- ggplot(data, aes(y=Pathway, x=ratio, color=`p-val`))+
    geom_point(aes(size=DEG)  )+
    scale_radius()+
    theme_bw()+
    labs(x = "ratio (DEG/N)") +
    scale_color_continuous(low = "red", high = "blue",
                           guide = guide_colorbar(reverse = TRUE))+
    theme(text = element_text(size=12))
  return(p)
}

# Dotplot de objeto enrich GO ##########################
dotPlotGOReport <- function(data, n = 20){
    names(data) <- gsub("P.DE", "p-val", names(data) )
  names(data) <- gsub("DE", "DEG", names(data) )
  data$ratio <- data$DEG/data$N
  data <- data[order(data$ratio, decreasing = F), ]
  data <- data[seq_len(n),]
  data$Term <- factor(data$Term, levels = data$Term)
  p <- ggplot(data, aes(y=Term, x=ratio, color=`p-val`))+
    geom_point(aes(size=DEG), stat="identity")+
    scale_radius()+
    theme_bw()+
    labs(x = "ratio (DEG/N)") +
    scale_color_continuous(low = "red", high = "blue", 
                           guide = guide_colorbar(reverse = TRUE))+
    theme(text = element_text(size=12))+
    scale_y_discrete(labels = function(Term){stringr::str_wrap(Term,50)})
  return(p)
}

# Heatmap de objeto enrich kegg ##########################
heatmapKeggReport <- function(kdt, nr){
  kdt <- kdt[nr, ]
  colourCount <- length(unique(kdt$DEG)) # number of levels
  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))
  kdt <- kdt %>% dplyr::select(Pathway, genes, DEG) %>% 
    separate_rows(genes) %>%
    mutate(Pathway = fct_inorder(Pathway)) %>% 
    mutate(Pathway = fct_rev(Pathway)) %>% 
    mutate(genes = fct_infreq(genes)) %>% 
    mutate(DEG = factor(DEG))
    yNum <- length(unique(kdt$Pathway))
    if(yNum <=35){ySize=12}else if(yNum>35 | yNum <=50){ySize=10}else{ySize=0}
    xNum <- length(unique(kdt$genes))
    if(xNum <=60){xSize=8}else if(xNum>60 | yNum <=80){xSize=7}else{xSize=0}
    kdt %>% ggplot(aes_(~genes, ~Pathway)) + 
    geom_tile(aes_(fill = ~DEG), color = 'black', size =0.2) +
    xlab(NULL) + ylab(NULL) +
    theme_minimal() +
    theme(panel.grid.major = element_line(colour = "gray88", size = 0.8),
          axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=xSize))+
    # scale_fill_continuous(low="blue", high="red", name = "N")
    # scale_fill_brewer(palette = "YlOrRd")
    scale_fill_manual(values = getPalette(colourCount))+
    theme(text = element_text(size=ySize, angle=0), plot.margin = unit(c(15,25,15,15), "pt"))
}

####
# HeatMapKegg con logFC ######
####
heatmapKeggLogFCReport <- function(kdt, res, nr){
    kdt <- kdt[nr, ]
    kk <- kdt %>% dplyr::select(Pathway, genes) %>% separate_rows(Pathway, genes, sep=",")
    kk$genes <- gsub(" ", "", kk$genes)
    resSig <- res[ which(res$SYMBOL %in% kk$genes), ]
    kk2 <- resSig %>% dplyr::select(SYMBOL, logFC, pval)
    kk3 <- left_join(kk, kk2, by = c("genes"="SYMBOL"))
    kk3$pval <- format(kk3$pval, scientific = TRUE, digits = 3)
    yNum <- length(unique(kdt$Pathway))
    if(yNum <=35){ySize=10}else if(yNum>35 | yNum <=50){ySize=8}else{ySize=0}
    xNum <- length(unique(kdt$genes))
    if(xNum <=60){xSize=7}else if(xNum>60 | yNum <=80){xSize=6}else{xSize=0}
    
    kk3 %>% ggplot(aes_(~genes, ~Pathway)) + 
    geom_tile(aes_(fill = ~logFC, label= ~pval), color = 'black', size =0.2) +
    xlab(NULL) + ylab(NULL) +
    theme_minimal() +
    theme(panel.grid.major = element_line(colour = "gray88", size = 0.8),
          axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=xSize))+
    scale_fill_gradient2(low="blue", mid = "gray88", high="red", name = "Log2FC")+
    # scale_fill_brewer(palette = "YlOrRd")
    # scale_fill_manual(values = getPalette(colourCount))+
    theme(text = element_text(size=ySize, angle=0), plot.margin = unit(c(15,25,15,15), "pt"))
}
heatmapKeggReport <- function(kdt, nr){
  kdt <- kdt[nr, ]
  colourCount <- length(unique(kdt$DEG)) # number of levels
  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))
  kdt <- kdt %>% dplyr::select(Pathway, genes, DEG) %>% 
    separate_rows(genes) %>%
    mutate(Pathway = fct_inorder(Pathway)) %>% 
    mutate(Pathway = fct_rev(Pathway)) %>% 
    mutate(genes = fct_infreq(genes)) %>% 
    mutate(DEG = factor(DEG))
    yNum <- length(unique(kdt$Pathway))
    if(yNum <=35){ySize=12}else if(yNum>35 | yNum <=50){ySize=10}else{ySize=0}
    xNum <- length(unique(kdt$genes))
    if(xNum <=60){xSize=8}else if(xNum>60 | yNum <=80){xSize=7}else{xSize=0}
    kdt %>% ggplot(aes_(~genes, ~Pathway)) + 
    geom_tile(aes_(fill = ~DEG), color = 'black', size =0.2) +
    xlab(NULL) + ylab(NULL) +
    theme_minimal() +
    theme(panel.grid.major = element_line(colour = "gray88", size = 0.8),
          axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=xSize))+
    # scale_fill_continuous(low="blue", high="red", name = "N")
    # scale_fill_brewer(palette = "YlOrRd")
    scale_fill_manual(values = getPalette(colourCount))+
    theme(text = element_text(size=ySize, angle=0), plot.margin = unit(c(15,25,15,15), "pt"))
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

# Función para hacer GSEA pathway #################################
gseaKeggReport <- function(res, specie){
  pathwayDataSet <- readRDS(paste0("./resources/",specie,"/GSEA/keggDataGSEA.Rds"))
  res.sh <- res
  #res.sh <- as.data.frame(lfcShrink(dds, coef=2, type="apeglm", res = results(dds)))
  res.sh <- res.sh[order(res.sh$log2FoldChange, decreasing = TRUE), ]
  res.sh$ENSEMBL <- rownames(res.sh)
  geneRank <- geneIdConverter( res.sh$ENSEMBL)
  resRank <- left_join(res.sh, geneRank, by=c("ENSEMBL"="ENSEMBL"))
  resRank <- resRank[!is.na(resRank$ENTREZID), c("ENTREZID","log2FoldChange") ]
  vectRank <- resRank$log2FoldChange
  attr(vectRank, "names") <- as.character(resRank$ENTREZID)
  mygsea <- clusterProfiler::GSEA(vectRank, 
                                  TERM2GENE = pathwayDataSet, 
                                  by="fgsea", pvalueCutoff = 0.1)
  mygsea <- DOSE::setReadable(mygsea, "org.Mm.eg.db", "ENTREZID")
  return(mygsea)
}

# Función para actualizar las bases de datos de kegg y GO #############
# esto mejora la velocidad de los enrich en unos 10 segs
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
## VOLCANY volcanoplot interactivo ###########################
volcanyReport <- function(res, padj = NULL, fcup = NULL,
                    fcdown = NULL,
                    col = NULL, genes= NULL ){

geneNames <- data.frame(ens=as.character(rownames(res)), symbol=as.character(res$GeneName_Symbol), stringsAsFactors = F)
res2 <- res %>% dplyr::select(GeneName_Symbol, log2FoldChange, padj) %>% as.data.frame()
rownames(res2) <- rownames(res)
res2$group <- "NS"
res2[which( res2$padj<=padj & ( res2$log2FoldChange>fcdown & res2$log2FoldChange<fcup) ), "group"] <- "Only Padj"
res2[which( res2$padj<=padj & res2$log2FoldChange>fcup), "group"] <- "FC Up & Padj"
res2[which( res2$padj<=padj & res2$log2FoldChange<fcdown), "group"] <- "FC Down & Padj"
res2[which( res2$padj>padj & ( res2$log2FoldChange<fcdown | res2$log2FoldChange>fcup ) ), "group"] <- "Only FC"

if(length(which(res2$group=="NS"))>1000 ){
  ns <- sample(which(res2$group=="NS"), length(which(res2$group=="NS"))-1000  )
  res2 <- res2[-c(ns), ]
  }
if(length(which(res2$group=="Only Padj"))>1000){
  pj <- sample(which(res2$group=="Only Padj"), length(which(res2$group=="Only Padj"))-1000)
  res2 <- res2[-c(ns), ]
}
if(length(which(res2$group=="Only FC"))>1000){
  fc <- sample(which(res2$group=="Only FC"), length(which(res2$group=="Only FC"))-1000)
  res2 <- res2[-c(fc), ]
}

if(!is.null(genes)){
  topTotal <- res2[which(res2$GeneName_Symbol %in% genes), ]
  topTotal$GeneName_Symbol <- as.character(topTotal$GeneName_Symbol)
  a <- list()
  for (i in seq_len(nrow(topTotal))) {
    m <- topTotal[i,]
    a[[i]] <- list(
      x = m[["log2FoldChange"]],
      y = -log10(m[["padj"]]),
      text = m[["GeneName_Symbol"]],
      xref = "x",
      yref = "y",
      xanchor = "left",
      yanchor = "bottom",
      showarrow = FALSE,
      arrowhead = 4,
      arrowsize = 0.5,
      ax = 20,
      ay = -40
    )
  }
}
line <- list(
  type = "line",
  line = list(color = "black", dash = "dash", width=0.7),
  xref = "paper"
)
lines <- list()
for (i in c(0.05)) {
  line[["x0"]] <- 0
  line[["x1"]] <- 1
  line[c("y0", "y1")] <- -log10(i)
  lines <- c(lines, list(line))
}
linev <- list(
  type = "line",
  line = list(color = "black", dash = "dash",width=0.7),
  yref = "paper"
)
linesv <- list()
for(j in c(-0.5,0.5)){
  linev[["y0"]] <- 0
  linev[["y1"]] <- 1
  linev[c("x0","x1")] <- j
  linesv <- c(linesv, list(linev))
}
linesT <- c(lines, linesv)


pal <- c(NS="gray", "Only Padj"="#7cccc3", "Only FC"="#d99c01", "FC Up & Padj"=col[1],
         "FC Down & Padj"=col[2] )
p <- res2 %>% plot_ly(x = ~log2FoldChange, y = ~(-log10(padj)), text = ~GeneName_Symbol,
                      mode = "markers", color = ~group, type = "scatter", size=I(5), colors=pal )%>% 
  layout(shapes = linesT ) %>% 
  layout(xaxis = list(zeroline=FALSE, title = "Log2 fold change"), yaxis=list(zeroline=FALSE, title = "-Log10 Padj"))
if(!is.null(genes)){ p <- p %>% layout(annotations = a) }
return(p)
}

# Volcano plot Miriam

VolcanoMiriReport <- function(res, padj, fcdown, fcup){
  res$log10FDR <- -log10(res$padj)
  res$sig <- as.factor((res$log2FoldChange > fcup & res$padj < padj) | (res$log2FoldChange < fcdown & res$padj < padj))
  ggplot(res,aes(x=log2FoldChange, y=log10FDR, color=sig)) +
  geom_point() +
  coord_cartesian() +
  ylab("-log10 FDR") +
  xlab("log2 fold change")
}

# Customized Volcano Plot ###############
CustomVolcanoReport <- function (toptable, lab, x, y, selectLab = NULL, xlim = c(min(toptable[[x]], 
                           na.rm = TRUE), max(toptable[[x]], na.rm = TRUE)), 
                           ylim = c(0, max(-log10(toptable[[y]])+1, na.rm = TRUE) ), xlab = bquote(~Log[2] ~ "fold change"), 
                           ylab = bquote(~-Log[10] ~ italic(P)), axisLabSize = 18, 
                           title = "Volcano plot highlighting the different groups of signification", subtitle = "", caption = paste0("Total = ", 
                           nrow(toptable), " variables"), titleLabSize = 18, subtitleLabSize = 14, 
                           captionLabSize = 14, pCutoff = 1e-05, pLabellingCutoff = pCutoff, 
                           FCcutoffDOWN = -1, FCcutoffUP = 1 , cutoffLineType = "longdash", cutoffLineCol = "black", 
                           cutoffLineWidth = 0.4, transcriptPointSize = 0.8, transcriptLabSize = 3, 
                           transcriptLabCol = "black", transcriptLabFace = "plain", 
                           transcriptLabhjust = 0, transcriptLabvjust = 1.5, pointSize = 2, 
                           labSize = 3, labCol = "black", labFace = "plain", labhjust = 0, 
                           labvjust = 1.5, boxedlabels = FALSE, boxedLabels = FALSE, 
                           shape = 19, shapeCustom = NULL, col = c("grey30", "forestgreen", 
                           "royalblue", "red2"), colCustom = NULL, colAlpha = 1/2, 
                           legend = c("NS", "Log2 FC", "P", "P & Log2 FC"), legendLabels = c("NS", 
                           expression(Only ~ log[2]~FC), "Only p-adjusted", expression(p-adjusted ~ and ~ log[2]~FC)), 
                           legendPosition = "top", legendLabSize = 14, 
                           legendIconSize = 4, legendVisible = TRUE, shade = NULL, 
                           shadeLabel = NULL, shadeAlpha = 1/2, shadeFill = "grey", 
                           shadeSize = 0.01, shadeBins = 2, drawconnectors = FALSE, 
                           drawConnectors = FALSE, widthConnectors = 0.5, typeConnectors = "closed", 
                           endsConnectors = "first", lengthConnectors = unit(0.01, "npc"), colConnectors = "grey10", 
                           hline = NULL, hlineType = "longdash", 
                           hlineCol = "black", hlineWidth = 0.4, vline = NULL, vlineType = "longdash", 
                           vlineCol = "black", vlineWidth = 0.4, gridlines.major = TRUE, 
                           gridlines.minor = TRUE, border = "partial", borderWidth = 0.8, 
                           borderColour = "black") 
{
  
  if (!is.numeric(toptable[[x]])) {
    stop(paste(x, " is not numeric!", sep = ""))
  }
  if (!is.numeric(toptable[[y]])) {
    stop(paste(y, " is not numeric!", sep = ""))
  }
  i <- xvals <- yvals <- Sig <- NULL
  if (!missing("transcriptPointSize")) {
    warning(paste0("transcriptPointSize argument deprecated in v1.4", 
                   " - please use pointSize"))
    pointSize <- transcriptPointSize
  }
  if (!missing("transcriptLabSize")) {
    warning(paste0("transcriptLabSize argument deprecated in v1.4", 
                   " - please use labSize"))
    labSize <- transcriptLabSize
  }
  if (!missing("transcriptLabCol")) {
    warning(paste0("transcriptLabCol argument deprecated in v1.4", 
                   " - please use labCol"))
    labCol <- transcriptLabCol
  }
  if (!missing("transcriptLabFace")) {
    warning(paste0("transcriptLabFace argument deprecated in v1.4", 
                   " - please use labFace"))
    labFace <- transcriptLabFace
  }
  if (!missing("transcriptLabhjust")) {
    warning(paste0("transcriptLabhjust argument deprecated in v1.4", 
                   " - please use labhjust"))
    labhjust <- transcriptLabhjust
  }
  if (!missing("transcriptLabvjust")) {
    warning(paste0("transcriptLabvjust argument deprecated in v1.4", 
                   " - please use labvjust"))
    labvjust <- transcriptLabvjust
  }
  if (!missing("boxedlabels")) {
    warning(paste0("boxedlabels argument deprecated in v1.4", 
                   " - please use boxedLabels"))
    boxedLabels <- boxedlabels
  }
  if (!missing("drawconnectors")) {
    warning(paste0("drawconnectors argument deprecated since v1.2", 
                   " - please use drawConnectors"))
    drawConnectors <- drawconnectors
  }
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[toptable[[x]] >= FCcutoffUP] <- "FC"
  toptable$Sig[toptable[[x]] <= FCcutoffDOWN] <- "FC"
  toptable$Sig[(toptable[[y]] < pCutoff)] <- "P"
  toptable$Sig[(toptable[[y]] < pCutoff) & (toptable[[x]] >= FCcutoffUP)] <- "FC_Pup"
  toptable$Sig[(toptable[[y]] < pCutoff) & (toptable[[x]] <= FCcutoffDOWN)] <- "FC_Pdown"
  toptable$Sig <- factor(toptable$Sig, levels = c("NS", "FC", "P", "FC_Pup","FC_Pdown"))
  
  if (min(toptable[[y]], na.rm = TRUE) == 0) {
    warning(paste("One or more p-values is 0.", "Converting to 10^-1 * current", 
                  "lowest non-zero p-value..."), call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(toptable[which(toptable[[y]] != 0), y], na.rm = TRUE) * 10^-1
  }
  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(), 
          plot.title = element_text(angle = 0, size = titleLabSize, 
          face = "bold", vjust = 1), plot.subtitle = element_text(angle = 0, 
          size = subtitleLabSize, face = "plain", vjust = 1), 
          plot.caption = element_text(angle = 0, size = captionLabSize, 
          face = "plain", vjust = 1), axis.text.x = element_text(angle = 0, size = axisLabSize, vjust = 1), 
          axis.text.y = element_text(angle = 0, 
          size = axisLabSize, vjust = 1), axis.title = element_text(size = axisLabSize), 
          legend.position = legendPosition, legend.key = element_blank(), 
          legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = legendLabSize), 
          title = element_text(size = legendLabSize), legend.title = element_blank())
  if (!is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(names(shapeCustom))), alpha = colAlpha, 
                 size = pointSize, na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom))), 
                 alpha = colAlpha, shape = shape, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(guide = TRUE)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           4) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(Sig)), alpha = colAlpha, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = c(NS = shape[1], FC = shape[2], 
                                    P = shape[3], FC_Pup = shape[4], FC_Pdown = shape[4] ),
                         labels = c(NS = legendLabels[1], 
                                    FC = legendLabels[2], P = legendLabels[3], FC_Pup = legendLabels[4],
                                    FC_Pdown = legendLabels[4]), 
                                    guide = TRUE)
  }
  else if (is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig), shape = factor(names(shapeCustom))), 
                 alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
      scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                    P = col[3], FC_Pup = col[4], FC_Pdown = col[5] ),
                         labels = c(NS = legendLabels[1],
                                    FC = legendLabels[2], P = legendLabels[3],
                                    FC_Pup = legendLabels[4],
                                    FC_Pdown = legendLabels[4])) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(shape = shape, 
               size = legendIconSize))) + geom_point(aes(color = factor(Sig)), 
               alpha = colAlpha, shape = shape, size = pointSize, 
               na.rm = TRUE, show.legend = legendVisible) + scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                    P = col[3], FC_Pup = col[4], FC_Pdown = col[5] ),
                         labels = c(NS = legendLabels[1],
                                    FC = legendLabels[2], P = legendLabels[3],
                                    FC_Pup = legendLabels[4],
                                    FC_Pdown = legendLabels[4]))
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(shape = c(NS = shape[1], 
      FC = shape[2], P = shape[3], FC_Pup = shape[4],FC_Pdown = shape[4]), size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig), shape = factor(Sig)), 
                 alpha = colAlpha, size = pointSize, na.rm = TRUE, 
                 show.legend = legendVisible) + scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                    P = col[3], FC_Pup = col[4], FC_Pdown = col[5] ),
                         labels = c(NS = legendLabels[1],
                                    FC = legendLabels[2], P = legendLabels[3],
                                    FC_Pup = legendLabels[4],
                                    FC_Pdown = legendLabels[4])) + 
      scale_shape_manual(values = c(NS = shape[1], FC = shape[2], 
                                    P = shape[3], FC_Pup = shape[4], FC_Pdown = shape[4]), guide = FALSE)
  }
  plot <- plot + xlab(xlab) + ylab(ylab) + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) + geom_vline(xintercept = c(FCcutoffDOWN,FCcutoffUP), linetype = cutoffLineType, colour = cutoffLineCol, 
                                        size = cutoffLineWidth) + geom_hline(yintercept = -log10(as.numeric(pCutoff)), 
                                        linetype = cutoffLineType, colour = cutoffLineCol, size = cutoffLineWidth)
  plot <- plot + labs(title = title, subtitle = subtitle, caption = caption)
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline, linetype = vlineType, 
                              colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline), 
                              linetype = hlineType, colour = hlineCol, size = hlineWidth)
  }
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = borderColour, 
                                                     fill = NA, size = borderWidth))
  }
  else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(size = borderWidth, 
                                                  colour = borderColour), panel.border = element_blank(), 
                         panel.background = element_blank())
  }
  else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  if (boxedLabels == FALSE) {
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_text_repel(data = subset(toptable, 
                  toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN))),
                  aes(label = subset(toptable, toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
                  size = labSize, segment.color = colConnectors, 
                  segment.size = widthConnectors, arrow = arrow(length = lengthConnectors, 
                                     type = typeConnectors, ends = endsConnectors), 
                                     hjust = labhjust, vjust = labvjust, colour = labCol, 
                                     fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_text_repel(data = subset(toptable, 
                                     !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                     !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                     segment.color = colConnectors, segment.size = widthConnectors, 
                                     arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                   ends = endsConnectors), hjust = labhjust, 
                                     vjust = labvjust, colour = labCol, fontface = labFace, 
                                     na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                               !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                               !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                               check_overlap = TRUE, hjust = labhjust, vjust = labvjust, 
                               colour = labCol, fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                                             toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN))), 
                                             aes(label = subset(toptable, toptable[[y]] < 
                                             pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
                               size = labSize, check_overlap = TRUE, hjust = labhjust, 
                               vjust = labvjust, colour = labCol, fontface = labFace, 
                               na.rm = TRUE)
    }
  }
  else {
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    toptable[[y]] < pLabellingCutoff & (toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)), 
                                                    aes(label = subset(toptable, 
                                                    toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
                                      size = labSize, segment.color = colConnectors, 
                                      segment.size = widthConnectors, arrow = arrow(length = lengthConnectors, 
                                                                                    type = typeConnectors, ends = endsConnectors), 
                                      hjust = labhjust, vjust = labvjust, colour = labCol, 
                                      fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                      !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                      !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                      segment.color = colConnectors, segment.size = widthConnectors, 
                                      arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                    ends = endsConnectors), hjust = labhjust, 
                                      vjust = labvjust, colour = labCol, fontface = labFace, 
                                      na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                hjust = labhjust, vjust = labvjust, colour = labCol, 
                                fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN))),
                                aes(label = subset(toptable, toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
                                size = labSize, hjust = labhjust, vjust = labvjust, 
                                colour = labCol, fontface = labFace, na.rm = TRUE)
    }
  }
  if (!is.null(shade)) {
    plot <- plot + stat_density2d(data = subset(toptable, 
                                                rownames(toptable) %in% shade), fill = shadeFill, 
                                  alpha = shadeAlpha, geom = "polygon", contour = TRUE, 
                                  size = shadeSize, bins = shadeBins, show.legend = FALSE, 
                                  na.rm = TRUE) + scale_fill_identity(name = shadeLabel, 
                                                                      labels = shadeLabel, guide = "legend")
  }
  return(plot)
}


# MA plot ##################
.levels <- function (x) 
    {
    if (!is.factor(x)) 
        x <- as.factor(x)
    levels(x)
    }

MAReport <- function (data, fdr = 0.05, fcDOWN = -1, fcUP = 1, genenames = NULL, detection_call = NULL, 
          size = NULL, font.label = c(12, "plain", "black"), label.rectangle = FALSE, 
          palette = c("#f7837b", "#1cc3c8", "darkgray"), top = 15, 
          select.top.method = c("padj", "fc"), main = NULL, xlab = "Log2 mean expression", 
          ylab = "Log2 fold change", ggtheme = theme_classic(), ...) 
{
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame", 
                              "DE_Results", "DESeqResults"))) 
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  if (!is.null(detection_call)) {
    if (nrow(data) != length(detection_call)) 
      stop("detection_call must be a numeric vector of length = nrow(data)")
  }
  else if ("detection_call" %in% colnames(data)) {
    detection_call <- as.vector(data$detection_call)
  }
  else detection_call = rep(1, nrow(data))
  if (is.null(list(...)$legend)) 
    legend <- c(0.12, 0.9)
  ss <- base::setdiff(c("baseMean", "log2FoldChange", "padj"), 
                      colnames(data))
  if (length(ss) > 0) 
    stop("The colnames of data must contain: ", paste(ss, collapse = ", "))
  if (is.null(genenames)) 
    genenames <- rownames(data)
  else if (length(genenames) != nrow(data)) 
    stop("genenames should be of length nrow(data).")
  sig <- rep(3, nrow(data))
  sig[which(data$padj <= fdr & data$log2FoldChange < 0 & data$log2FoldChange <= (as.numeric(fcDOWN)) & detection_call == 1)] = 2
  sig[which(data$padj <= fdr & data$log2FoldChange > 0 & data$log2FoldChange >= (as.numeric(fcUP)) & detection_call == 1)] = 1
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange, padj = data$padj, sig = sig)
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- .levels(data$sig) %>% as.numeric()
  palette <- palette[.lev]
  new.levels <- c(paste0("Up: ", sum(sig == 1)), paste0("Down: ", sum(sig == 2)), "NS") %>% .[.lev]
  data$sig <- factor(data$sig, labels = new.levels)
  select.top.method <- match.arg(select.top.method)
  if (select.top.method == "padj") 
    data <- data[order(data$padj), ]
  else if (select.top.method == "fc") 
    data <- data[order(abs(data$lfc), decreasing = TRUE), 
                 ]
  labs_data <- stats::na.omit(data)
  labs_data <- subset(labs_data, padj <= fdr & name != "" & 
                        (lfc >= fcUP | lfc <=fcDOWN) )
  labs_data <- utils::head(labs_data, top)
  font.label <- ggpubr:::.parse_font(font.label)
  font.label$size <- ifelse(is.null(font.label$size), 12, 
                            font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", 
                             font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", 
                            font.label$face)
  set.seed(42)
  mean <- lfc <- sig <- name <- padj <- NULL
  p <- ggplot(data, aes(x = log2(mean + 1), y = lfc)) + geom_point(aes(color = sig), size = size)
  if (label.rectangle) {
    p <- p + ggrepel::geom_label_repel(data = labs_data, 
                                       mapping = aes(label = name), box.padding = unit(0.35,"lines"), point.padding = unit(0.3, "lines"), 
                                       force = 1, fontface = font.label$face, size = font.label$size/3, 
                                       color = font.label$color)
  }
  else {
    p <- p + ggrepel::geom_text_repel(data = labs_data, 
                                      mapping = aes(label = name), box.padding = unit(0.35,  "lines"), point.padding = unit(0.3, "lines"), 
                                      force = 1, fontface = font.label$face, size = font.label$size/3, 
                                      color = font.label$color)
  }
  p <- p + scale_x_continuous(breaks = seq(0, max(log2(data$mean + 
                                                         1)), 2)) + labs(x = xlab, y = ylab, title = main, color = "") + 
    geom_hline(yintercept = c(0, (fcDOWN), (fcUP)), linetype = c(1, 2, 2), color = c("black", "black", "black"))
  p <- ggpar(p, palette = palette, ggtheme = ggtheme, ...)
  return(p)
}


# VST ###############  SIN USO CREO !!!!!
VSTReport <- function (object, blind = TRUE, nsub = 1000, fitType = "parametric") 
{
  if (nrow(object) < nsub) {
    stop("less than 'nsub' rows,\n  it is recommended to use varianceStabilizingTransformation directly")
  }
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.matrix(object)) {
    matrixIn <- TRUE
    object <- DESeqDataSetFromMatrix(object, DataFrame(row.names = colnames(object)), 
                                     ~1)
  }
  else {
    if (blind) {
      design(object) <- ~1
    }
    matrixIn <- FALSE
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  baseMean <- rowMeans(counts(object, normalized = TRUE))
  if (sum(baseMean > 5) < nsub) {
    stop("less than 'nsub' rows with mean normalized count > 5, \n  it is recommended to use varianceStabilizingTransformation directly")
  }
  object.sub <- object[baseMean > 5, ]
  baseMean <- baseMean[baseMean > 5]
  o <- order(baseMean)
  idx <- o[round(seq(from = 1, to = length(o), length = nsub))]
  object.sub <- object.sub[idx, ]
  object.sub <- estimateDispersionsGeneEst(object.sub, quiet = TRUE)
  object.sub <- estimateDispersionsFit(object.sub, fitType = fitType, 
                                       quiet = TRUE)
  suppressMessages({
    dispersionFunction(object) <- dispersionFunction(object.sub)
  })
  vsd <- varianceStabilizingTransformation(object, blind = FALSE)
  if (matrixIn) {
    return(assay(vsd))
  }
  else {
    return(vsd)
  }
}

# Heatmap #############
heatReport <- function (vsd, n = 40, intgroup = "AAV", sampleName = "condition",
                      specie="Mm", customColor = c("red","blue")) 
    {
      require("EnsDb.Mmusculus.v79")
      require("org.Mm.eg.db")
      require("EnsDb.Hsapiens.v86")
      require("org.Hs.eg.db")
      require("EnsDb.Rnorvegicus.v79")
      require("org.Rn.eg.db") 
      
      if(specie=="Mm"){
        ensdb <- EnsDb.Mmusculus.v79
        orgdb <- org.Mm.eg.db
      }
      else{
        ensdb <- EnsDb.Hsapiens.v86
        orgdb <- org.Hs.eg.db
      }
    #vsd <- vst(data)
    topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), n)
    mat  <- assay(vsd)[ topVarGenes, ]
    mat  <- mat - rowMeans(mat)
    if (!all(intgroup %in% names(colData(vsd)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }

    if(length(intgroup)>1){
      df <- as.data.frame(colData(vsd)[, intgroup[1:2], drop = FALSE])
    } else{
      df <- as.data.frame(colData(vsd)[, intgroup, drop = FALSE])
    }
    
    annot <- NULL
    annot$ENSEMBL <- rownames(mat)
    annot$SYMBOL <-  mapIds(ensdb, keys=rownames(mat),
                            column="SYMBOL",keytype="GENEID")
    annot$SYMBOL1 <- mapIds(orgdb, keys = rownames(mat),
                            column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first') 
    annot$description <- mapIds(orgdb, keys = rownames(mat),
                                column = 'GENENAME', keytype = 'ENSEMBL', multiVals = 'first')
    annot <- as.data.frame(annot)
    consensus <- data.frame('Symbol'= ifelse(!is.na(annot$SYMBOL), as.vector(annot$SYMBOL),
                                             ifelse(!is.na(annot$SYMBOL1),as.vector(annot$SYMBOL1),
                                                    as.vector(annot$ENSEMBL))), stringsAsFactors = F)
    ann_colors<-list()
    ann_colors[[intgroup[1]]] <- customColor
    names(ann_colors[[intgroup[1]]]) <- c(levels(df[[intgroup[1]]]))
    sizesDf <- data.frame( ch = c(rep(14,20), rep(12,20),rep(10,10), 
                                  rep(8,10), rep(7,10), rep(6,10), rep(5,20), rep(4,20)), 
                           fsr = c(rep(10,50), rep(8,10), rep(7,10), rep(6,10), rep(0.1, 40) ))
    ch <- sizesDf$ch[ nrow(mat) ]
    fsr <- sizesDf$fsr[ nrow(mat) ]
    pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE,
             show_colnames=TRUE, show_rownames = TRUE, annotation_col = df,
             labels_col = as.character(vsd[[sampleName]]),
             labels_row = as.character(consensus$Symbol),
             cellwidth = 14, cellheight = ch,
             fontsize_row = fsr,
             annotation_colors = ann_colors,
             main = "Heatmap top variant genes on normalized data")
}
## New heatmap plotly
heat2Report <- function (vsd, n = 40, intgroup = NULL, sampleName = NULL,
                      specie="Mm", customColor = NULL ) 
    {
      require("EnsDb.Mmusculus.v79")
      require("org.Mm.eg.db")
      require("EnsDb.Hsapiens.v86")
      require("org.Hs.eg.db")
      require("EnsDb.Rnorvegicus.v79")
      require("org.Rn.eg.db") 
      
      if(specie=="Mm"){
        ensdb <- EnsDb.Mmusculus.v79
        orgdb <- org.Mm.eg.db
      } else{
        ensdb <- EnsDb.Hsapiens.v86
        orgdb <- org.Hs.eg.db
      }
    #vsd <- vst(data)
    topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), n)
    mat  <- assay(vsd)[ topVarGenes, ]
    mat  <- mat - rowMeans(mat)
    if (!all(intgroup %in% names(colData(vsd)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }

    if(length(intgroup)>1){
      df <- as.data.frame(colData(vsd)[, intgroup[1:2], drop = FALSE])
    } else{
      df <- as.data.frame(colData(vsd)[, intgroup, drop = FALSE])
    }
    
    annot <- NULL
    annot$ENSEMBL <- rownames(mat)
    annot$SYMBOL <-  mapIds(ensdb, keys=rownames(mat),
                            column="SYMBOL",keytype="GENEID")
    annot$SYMBOL1 <- mapIds(orgdb, keys = rownames(mat),
                            column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first') 
    annot$description <- mapIds(orgdb, keys = rownames(mat),
                                column = 'GENENAME', keytype = 'ENSEMBL', multiVals = 'first')
    annot <- as.data.frame(annot)
    consensus <- data.frame('Symbol'= ifelse(!is.na(annot$SYMBOL), as.vector(annot$SYMBOL),
                                             ifelse(!is.na(annot$SYMBOL1),as.vector(annot$SYMBOL1),
                                                    as.vector(annot$ENSEMBL))), stringsAsFactors = F)
    #colores definidos por el usuario para primera variable
    ann_colors<-list()
    ann_colors[[intgroup[1]]] <- customColor
    names(ann_colors[[intgroup[1]]]) <- c(levels(df[[intgroup[1]]]))
    
    # colores calculados para segunda variable
    if(length(intgroup)>1){
    numCol <- length( unique( df[,intgroup[[2]]]))
    ann_colors[[intgroup[2]]] <- colorRampPalette( brewer.pal(11, "Spectral") )(numCol)
    names(ann_colors[[intgroup[2]]]) <- c(levels(as.factor(df[[intgroup[2]]]) ))}

    ann <- c(ann_colors[[intgroup[1]]], ann_colors[[intgroup[2]]])
    
    sizesDf <- data.frame( ch = c(rep(14,20), rep(12,20),rep(10,10), 
                                  rep(8,10), rep(7,10), rep(6,10), rep(5,20), rep(4,20)), 
                           fsr = c(rep(10,50), rep(8,10), rep(7,10), rep(6,10), rep(1, 40) ))
    ch <- sizesDf$ch[ nrow(mat) ]
    fsr <- sizesDf$fsr[ nrow(mat) ]
    if(nrow(mat)>80){labrow = rep(NA, nrow(mat))}else{labrow = as.character(consensus$Symbol)}
    heatmaply(mat, labRow = labrow, col_side_colors = df,
              col_side_palette = ann, labCol = as.character(vsd[[sampleName]] ), fontsize_row = fsr,
              margins = c(50,50,20,0)  )
}

# cluster #############

clusterReport <- function(vsd, intgroup = "condition")
  {
  #vsd <- vst(data)
  sampleDists_vsd <- dist(t(assay(vsd)))
  sampleDistMatrix_vsd <- as.matrix( sampleDists_vsd )
  rownames(sampleDistMatrix_vsd) <- vsd[[intgroup]]
  colnames(sampleDistMatrix_vsd) <- vsd[[intgroup]]
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  heatmaply(sampleDistMatrix_vsd, colors = colors, margins = c(50,50,50,0)  )
  # pheatmap(sampleDistMatrix_vsd,
  #          clustering_distance_rows = sampleDists_vsd,
  #          clustering_distance_cols = sampleDists_vsd,
  #          col = colors, main = 'Heatmap clustering of samples on normalized data')
}






#############  TOP6 genes #########################


plotCountsSymbolReport <- function (dds, gene, intgroup = "condition", normalized = TRUE,
                              transform = TRUE, main, xlab = "group", returnData = FALSE,
                              replaced = FALSE, pc, specie, ...){
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) &
                                                         (gene >= 1 & gene <= nrow(dds)))))
  if (!all(intgroup %in% names(colData(dds))))
    stop("all variables in 'intgroup' must be columns of colData")
  if (!returnData) {
    if (!all(sapply(intgroup, function(v) is(colData(dds)[[v]],
                                             "factor")))) {
      stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
    }
  }
  if (missing(pc)) {
    pc <- if(transform) 0.5
    else 0
  }
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  if(specie=="Mm"){
      ensdb <- EnsDb.Mmusculus.v79
      orgdb <- org.Mm.eg.db
  }
    else{
        ensdb <- EnsDb.Hsapiens.v86
        orgdb <- org.Hs.eg.db
    }
  geneE <- mapIds(orgdb, keys = gene, column = "ENSEMBL", keytype = "SYMBOL")
  cnts <- counts(dds, normalized = normalized, replaced = replaced)[geneE,]
  intgroup <- intgroup[1]
  group <- if (length(intgroup) == 1) {
    colData(dds)[[intgroup]]
  }
  else if (length(intgroup) == 2) {
    lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]),
                              levels(colData(dds)[[intgroup[2]]]), function(x,
                                                                            y) paste(x, y, sep = ":"))))
    droplevels(factor(apply(as.data.frame(colData(dds)[,
                                                       intgroup, drop = FALSE]), 1, paste, collapse = ":"),
                      levels = lvls))
  }
  else {
    factor(apply(as.data.frame(colData(dds)[, intgroup,
                                            drop = FALSE]), 1, paste, collapse = ":"))
  }
  #rownames(cnts) <- res$GeneName_Symbol
  data <- data.frame(count = cnts + pc, group = as.integer(group))
  logxy <- if (transform) "y" else ""
  ylab <- ifelse(normalized, "normalized counts")
  #colors = c("#008000","#800080")
  if (returnData)
    return(data.frame(count = data$count, colData(dds)[intgroup]))
  plot(data$group + runif(ncol(dds), -0.05, 0.05), data$count, #col=colors[(dds)[intergroup]],
       xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n",
       xlab = xlab, ylab = ylab, main = paste0("Expression of ", gene), ...)
  axis(1, at = seq_along(levels(group)), levels(group))
  #text(data$group + runif(ncol(dds), -0.05, 0.05), data$count, labels=colnames(dds))
}


# Boxplot Violin plot ###########################
boxViolinReport <- function(datos=NULL, vsd=NULL, names=NULL, boxplotswitch=NULL,
                      intgroup=NULL, customColor = NULL){
    data <- vsd
    df <- assay(data)
    condition <- colData(data)[[intgroup[1]]] # variables()
    samples <- colData(data)[[names[1]]]
    df <- as.data.frame(t(df))
    df$condition <- as.character(condition)
    df$samples <- as.character(samples)
    dflong <- df %>% pivot_longer(c(-condition,-samples),
                            names_to = "genes",
                            values_to = "count")
    p <- dflong %>% ggplot(aes(x = samples, y = count, fill = condition))
    if(!isTRUE(boxplotswitch)){
        p <- p + geom_boxplot() + 
            ggtitle('Box plot for the normalized counts per gene and sample') +
            theme(axis.text.x = element_text(angle = 90), text = element_text(size = 10)) +
            scale_fill_manual( values = customColor )
        return(p)
    } else {
        p <- p + geom_violin() + 
            ggtitle('Violin plot for the normalized counts per gene and sample') +
            theme(axis.text.x = element_text(angle = 90), text = element_text(size = 10))+
          scale_fill_manual( values = customColor )
        return(p)
        }
          }


choices_brewer <- list(
  "Blues" = brewer_pal(palette = "Blues")(9),
  "Greens" = brewer_pal(palette = "Greens")(9),
  "Reds" = brewer_pal(palette = "Reds")(9),
  "Oranges" = brewer_pal(palette = "Oranges")(9),
  "Purples" = brewer_pal(palette = "Purples")(9),
  "Greys" = brewer_pal(palette = "Greys")(9)
)

choices_brewer2 <- list(
  as.list(rev(brewer_pal(palette = "Blues")(9))),
  as.list(rev(brewer_pal(palette = "Greens")(9))),
  as.list(rev(brewer_pal(palette = "Reds")(9))),
  as.list(rev(brewer_pal(palette = "Oranges")(9))),
  as.list(rev(brewer_pal(palette = "Purples")(9))),
  as.list(rev(brewer_pal(palette = "Greys")(9)))
)


## karyoplotter ########################
customkaryploterReport <- function(genome = "mm10", plot.type = 1, ideogram.plotter = kpAddCytobands, labels.plotter = kpAddChromosomeNames, chromosomes = "auto", zoom = NULL, cytobands = NULL, plot.params = NULL, use.cache = TRUE, main = NULL, bg="white" ){
  require("karyoploteR")
    if (is.null(genome)) 
        stop("genome cannot be NULL.")
    if (!is.null(cytobands)) {
        cytobands <- tryCatch(toGRanges(cytobands), error = function(e) {
        })
        if (!methods::is(cytobands, "GRanges")) 
            stop("'cytobands' must be NULL, a GRanges or something accepted by regioneR::toGRanges")
    }
    if (!is.null(zoom)) {
        zoom <- regioneR::toGRanges(zoom)
        if (!methods::is(zoom, "GRanges")) 
            stop("'zoom' must be NULL or a GRanges object")
        if (length(zoom) > 1) {
            warning("The zoom parameter has more than one region. Only the first one will be used.")
            zoom <- zoom[1]
        }
    }
    if (is.null(plot.params)) {
        plot.params <- getDefaultPlotParams(plot.type)
    }
    gr.genome <- NULL
    genome.name <- NULL
    if (methods::is(genome, "GRanges")) {
        gr.genome <- genome
    } else {
        if (methods::is(genome, "BSgenome")) {
            genome <- seqinfo(genome)
        }
        if (methods::is(genome, "Seqinfo")) {
            gr.genome <- as(genome, "GRanges")
            genome.name <- genome(genome)[1]
        }
        else {
            if (methods::is(genome, "character") & use.cache == 
                TRUE) {
                if (genome %in% names(karyoploteR:::data.cache[["genomes"]])) {
                    gr.genome <- karyoploteR:::data.cache[["genomes"]][[genome]]
                    genome.name <- genome
                }
            }
        }
    }
    if (is.null(gr.genome)) {
        gr.genome <- tryCatch(regioneR::getGenomeAndMask(genome = genome, 
            mask = NA)$genome, error = function(e) {
            stop("It was not possible to identify or load the requested genome. ", 
                e)
        })
    }
    chr.names <- as.character(GenomeInfoDb::seqnames(gr.genome))
    if (any(duplicated(chr.names))) {
        stop(paste0("There are duplicate chromosome names in the genome. Chromosome names must be unique. Chromosome names are: ", 
            paste0(chr.names, collapse = ", ")))
    }
    names(gr.genome) <- chr.names
    if (!is.null(zoom)) {
        if (!IRanges::overlapsAny(zoom, gr.genome)) {
            stop("You are trying to set the zoom to a region not part of the current genome.")
        }
        else {
            chromosomes <- as.character(GenomeInfoDb::seqnames(zoom))
        }
    }
    if (!is.null(chromosomes) && any(chromosomes != "all")) {
        if (length(chromosomes) == 1 && (chromosomes %in% c("canonical", 
            "autosomal", "auto"))) {
            if (!is.null(genome.name) && is.character(genome.name)) {
                if (chromosomes == "auto") 
                    chromosomes <- "canonical"
                tryCatch(expr = {
                    gr.genome <- filterChromosomes(gr.genome, 
                        organism = genome.name, chr.type = chromosomes)
                }, error = function(e) {
                    message("WARNING: There was an error when filtering the chromosomes and selecting only ", 
                        chromosomes, " chromosomes.  Falling back to using the unfiltered genome. \n", 
                        e)
                })
            }
            else {
                if (chromosomes != "auto") {
                    message("NOTE: It is only possible to filter chromosomes using named ", 
                        "chromosome classes (i.e. 'canonical', 'autosomal') when the genome ", 
                        "is specified by name (i.e. 'hg19', 'mm10'). Please, either ", 
                        "use a genome specified by name or explicitly select the ", 
                        "chromosomes to plot (i.e. chromosomes=c('chr1', 'chr2') ). ", 
                        " Falling back to using the unfiltered genome.")
                }
            }
        }
        else {
            if (!all(chromosomes %in% as.character(seqnames(gr.genome)))) {
                message("NOTE: Not all requested chromosomes are part of the genome. Trying to filter as requested. ", 
                    "   * Requested Chromosomes: ", paste0(chromosomes, 
                        collapse = ", "), "   * Chromosomes in the genome: ", 
                    paste0(as.character(seqnames(gr.genome)), 
                        collapse = ", "))
            }
            tryCatch(expr = {
                gr.genome <- filterChromosomes(gr.genome, keep.chr = chromosomes)[chromosomes]
            }, error = function(e) {
                message("WARNING: There was an error when filtering the chromosomes. Falling back to using the unfiltered genome. \n", 
                    e)
            })
        }
    }
    if (length(gr.genome) == 0) {
        stop("The genome has no chromosomes left after filtering. Cannot plot with no chromosomes.")
    }
    if (is.null(cytobands)) {
        if (!is.null(genome.name) && all(is.character(genome.name))) {
            cytobands <- getCytobands(genome.name)
        }
        else {
            if (all(is.character(genome))) {
                cytobands <- getCytobands(genome)
            }
        }
        if (!is.null(cytobands) && length(cytobands) > 0) {
            cytobands <- GenomeInfoDb::keepSeqlevels(cytobands, 
                value = GenomeInfoDb::seqlevels(gr.genome), 
                pruning.mode = "coarse")
        }
    }
    kp <- list()
    class(kp) <- "KaryoPlot"
    kp$plot.params <- plot.params
    if (!is.null(genome.name)) {
        kp$genome.name <- genome.name
    } else {
        kp$genome.name <- "custom"
    }
    kp$chromosomes <- as.character(GenomeInfoDb::seqlevels(gr.genome))
    kp$chromosome.lengths <- stats::setNames(end(gr.genome), 
        seqnames(gr.genome))
    kp$genome <- gr.genome
    kp$cytobands <- cytobands
    kp$plot.type <- plot.type
    if (is.null(zoom)) {
        kp$plot.region <- kp$genome
        kp$zoom <- FALSE
    } else {
        kp$plot.region <- zoom
        names(kp$plot.region) <- as.character(seqnames(kp$plot.region))
        kp$zoom <- TRUE
    }
    coordChangeFunctions <- karyoploteR:::getCoordChangeFunctions(karyoplot = kp)
    kp$coord.change.function <- coordChangeFunctions$coordChangeFunction
    kp$ideogram.mid <- coordChangeFunctions$ideogramMid
    kp$chromosome.height <- coordChangeFunctions$chr.height
    kp$graphical.par <- list()
    kp$graphical.par$old.par <- graphics::par(no.readonly = TRUE)
    graphics::par(mar = c(0, 0, 0, 0) + 0.1)
    graphics::par(bg=bg)
    graphics::par(col = "white")
    kp$beginKpPlot <- function() {
        graphics::par(kp$graphical.par$new.par)
    }
    kp$endKpPlot <- function() {
        graphics::par(kp$graphical.par$old.par)
    }
    on.exit(kp$endKpPlot())
    pp <- plot.params
    if (plot.type %in% c(1, 2, 6)) {
        xlim <- c(0, 1)
        ylim <- c(0, pp$bottommargin + length(gr.genome) * kp$chromosome.height + 
            pp$topmargin)
    } else {
        if (plot.type %in% c(3, 4, 5, 7)) {
            xlim <- c(0, 1)
            ylim <- c(0, pp$bottommargin + kp$chromosome.height + 
                pp$topmargin)
        }
    }
    graphics::plot(0, type = "n", xlim = xlim, ylim = ylim, 
        axes = FALSE, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
    kp$plot <- list()
    p <- graphics::par("usr")
    kp$plot$xmin <- p[1]
    kp$plot$xmax <- p[2]
    kp$plot$ymin <- p[3]
    kp$plot$ymax <- p[4]
    kp$graphical.par$new.par <- graphics::par(no.readonly = TRUE)
    if (!is.null(ideogram.plotter)) {
        kpAddCytobands(kp, color.table = )
    }
    if (!is.null(labels.plotter)) {
        kpAddChromosomeNames(kp, col="black")
    }
    if (!is.null(main)) {
        kpAddMainTitle(kp, main)
    }
    return(kp)
}
  
  
krtpReport <- function(res, specie="Mm", pval, fcdown, 
                 fcup, bg="white", coldown="#1f31ff", colup="#DC143C", annotation){
  require(karyoploteR)
  fileAnnot <- paste0("./resources/",specie,"/cytoband/",specie,"_annot.txt")
  annot <- read.table(fileAnnot, header = F, sep = "\t")
  res2 <- res[ res$pval <pval & (res$logFC<(fcdown) | res$logFC>fcup),]
  res3 <- as.data.frame(res2)
  if(annotation=="ensg"){
    res3$genes <- res3$ENSEMBL
    genes <- left_join(annot, res3, by = c("V1"="genes"))
    }
  if(annotation=="symbol"){
    res3$genes <- res3$ENSEMBL
    genes <- left_join(annot, res3, by = c("V1"="genes"))
  }
  sig <- which( !is.na(genes$pval) )
  genes <- genes[sig,]
  if(annotation =="ensg"){genesv1 <- genes$V1}
  if(annotation =="symbol"){genesv1 <- genes$V1}
  A <- data.frame(chr = paste0("chr",genes$V2), start = genes$V3,
                  end=genes$V4, x = genesv1, y = genes$logFC)
  genesSig <- toGRanges(A)
  one <- getDefaultPlotParams(2)
  one$ideogramheight <- 300
  one$data2height <- 300
  dfRanges <- readRDS(paste0("./resources/",specie,"/cytoband/",specie,"_genomicRanges.Rds") )
  dfIdeoBanda <- readRDS(paste0("./resources/",specie,"/cytoband/",specie,"_cytoBand.Rds"))
  kp <- customkaryploterReport(genome = dfRanges, cytobands = dfIdeoBanda,
                         plot.params = one, plot.type = 2, bg=bg, use.cache = FALSE) %>%
  kpPlotRegions(genesSig[genesSig$y>0,], col = colup, data.panel = 1)  %>%
      kpPlotRegions(genesSig[genesSig$y<0,], col = coldown, data.panel = 2)
  #return(kp)
}

## funcion para generar datos para ideograma para una especie concreta #############

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
  out <- readLines(dat)
  datos <- jsonlite::prettify(out)
  dfIdeoBand <- jsonlite::fromJSON(datos)
  dfIdeoBanda <- as.data.frame.list(dfIdeoBand$cytoBandIdeo[[1]])
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
# customVisNet #######################################
customVisNetReport <- function( enrich, kggDT, nTerm = NULL, up = NULL, down = NULL ){
    require(visNetwork)
    require(scales)
    enrich$genes <- kggDT$genes
    enrich <- enrich %>% arrange(P.DE)
    enrich <- enrich[nTerm, ]
    enrich$genes <- gsub(",", ";", enrich$genes)
    enrich$genes <- gsub(" ", "", enrich$genes)
    if( dim(enrich)[2]==8 ){
        names(enrich) <- c("Term","Ont","N","DE","P.DE","id","genes","level")
    } else if(dim(enrich)[2]==6){
        names(enrich) <- c("Term","N","DE","P.DE","id","genes")
    }
    edges <- separate_rows(enrich, Term, genes, sep=";")
    if( dim(enrich)[2]==8 ){
        edgesf <- edges[, c(1, 7)]
        edges <- edges %>% dplyr::select(-Ont)    
    } else if(dim(enrich)[2]==6){
        edgesf <- edges[, c(1, 6)]
    }
    names(edgesf) <- c("from", "to")
    edgesf$to <- gsub(" ", "", edgesf$to)
    nd1 <- edges[, 1:5]
    nd1$P.DE <- nd1$P.DE
    nd2 <- as.data.frame( 
    cbind( Term = edges$genes, N = NA, DE = NA, P.DE = NA, id = NA))
    nd2$Term <- gsub(" ", "", nd2$Term)
    nd2$P.DE <- NA
    nd2$id <- nd2$Term
    nd3 <- rbind(nd1, nd2)
    nd3 <- dplyr::distinct(nd3)
        
    nd3$title <- paste0( nd3$id, "<br>", "pval: ",
                         formatC(nd3$P.DE, format = "e", digits = 3),"<br>",
                         "DE: ", nd3$DE )
    nd3$DE <- ifelse( is.na(nd3$DE), 
                      min( as.numeric( nd3$DE[ !is.na( nd3$DE ) ] ) ),
                      nd3$DE )
    nd3$title <- gsub("<br>pval:   NA<br>DE: NA", "", nd3$title)
        
    pvalCol <- rescale(nd3$P.DE, to = c(0, 1))
    colores <- scales::seq_gradient_pal("red", "blue")(pvalCol)
    colores <- ifelse(is.na(colores), "#bbceed", colores)
        
    nodesf <- data.frame( id = nd3$Term, label = nd3$Term, group = 1,
            value = nd3$DE, color = colores, shadow = F, title = nd3$title,
            stringsAsFactors = F )
    if(!is.null(up) ){
      genesUp <- which(nodesf$id %in% up)
      nodesf$color[genesUp] <- "#ffa200"
      }
    if(!is.null(down)){
      genesDown <- which(nodesf$id %in% down)
      nodesf$color[genesDown] <- "#91ebff"
      }
    return(list(nodes = nodesf, edges = edgesf))
} 
# GoBarplot ########################
goBarplotReport <- function(enrichGO=NULL, resGO=NULL, genes=NULL,
                      category=NULL, nrows=NULL ){
    require(GOplot)
    go <- enrichGO
    go <- go[ go$Ont==category, ]
    if(is.null(nrows) | length(nrows)<2 ){
        totalRows <- min(90, dim(go)[1] )
        go <- go[ seq_len(totalRows), ]
    } else{
        go <- go[nrows, ]
        }
    res <- resGO
    go2 <- go %>% group_by(Ont) %>% as.data.frame()
    goDT <- go2DT(go2, genes)
    # preparar tabla GO
    go2$genes <- goDT$genes
    go2 <- go2 %>% dplyr::select(Ont,go_id,Term,genes, P.DE)
    names(go2) <- c("Category","ID", "Term", "Genes", "adj_pval")
    #preparar tabla genelist
    names(res)
    res2 <- res %>% dplyr::select(SYMBOL, logFC, pval)
    names(res2) <- c("ID","logFC","adj.P.Val")
    # crear objeto circ
    circ <- circle_dat(go2, res2)
    GOBar(circ)
}
# data2circle ##############################
data2circleReport <- function(go=NULL, res=NULL, genes=NULL){
  go2 <- go %>% group_by(Ont) %>% as.data.frame()
  goDT <- go2DT(go2, genes)
  # preparar tabla GO
  go2$genes <- goDT$genes
  go2 <- go2 %>% dplyr::select(Ont,go_id,Term,genes, P.DE)
  names(go2) <- c("Category","ID", "Term", "Genes", "adj_pval")
  #preparar tabla genelist
  names(res)
  res2 <- res %>% dplyr::select(SYMBOL, logFC, pval)
  names(res2) <- c("ID","logFC","adj.P.Val")
  # crear objeto circ
  require(GOplot)
  circ <- circle_dat(go2, res2)
  return(circ)
}

# GO circle ################################
circleReport <- function (data, title, nsub, rad1, rad2, table.legend = F, zsc.col, 
          lfc.col, label.size, label.fontface) {
    require(GOplot)
    xmax <- y1 <- zscore <- y2 <- ID <- logx <- logy2 <- logy <- logFC <- NULL
    if (missing(title))
        title <- ""
    if (missing(nsub))
        if (dim(data)[1] > 10)
            nsub <- 10
        else
            nsub <- dim(data)[1]
        if (missing(rad1))
            rad1 <- 2
        if (missing(rad2))
            rad2 <- 3
        if (missing(zsc.col))
            zsc.col <- c("red", "white", "blue")
        if (missing(lfc.col))
            lfc.col <- c("cornflowerblue", "firebrick1")
        else
            lfc.col <- rev(lfc.col)
        if (missing(label.size))
            label.size = 5
        if (missing(label.fontface))
            label.fontface = "bold"
        data$adj_pval <- -log(data$adj_pval, 10)
        suby <- data[!duplicated(data$term),]
        if (is.numeric(nsub) == T) {
            suby <- suby[1:nsub,]
        }
        else {
            if (strsplit(nsub[1], ":")[[1]][1] == "GO") {
                suby <- suby[suby$ID %in% nsub,]
            }
            else {
                suby <- suby[suby$term %in% nsub,]
            }
            nsub <- length(nsub)
        }
        N <- dim(suby)[1]
        r_pval <- round(range(suby$adj_pval), 0) + c(-2, 2)
        ymax <- c()
        for (i in 1:length(suby$adj_pval)) {
            val <- (suby$adj_pval[i] - r_pval[1]) / (r_pval[2] - r_pval[1])
            ymax <- c(ymax, val)
        }
        df <- data.frame(x = seq(0, 10 - (10 / N), length = N), 
                         xmax = rep(10 / N -0.2, N),
                         y1 = rep(rad1, N), y2 = rep(rad2, N),
                         ymax = ymax, zscore = suby$zscore, ID = suby$ID)
        scount <- data[!duplicated(data$term), which(colnames(data) == "count")][1:nsub]
        idx_term <- which(!duplicated(data$term) == T)
        xm <- c()
        logs <- c()
        for (sc in 1:length(scount)) {
            idx <- c(idx_term[sc], idx_term[sc] + scount[sc] - 1)
            val <- stats::runif(scount[sc], df$x[sc] + 0.06, (df$x[sc] + df$xmax[sc] - 0.06))
            xm <- c(xm, val)
            r_logFC <- round(range(data$logFC[idx[1]:idx[2]]), 0) + c(-1, 1)
            for (lfc in idx[1]:idx[2]) {
                val <- (data$logFC[lfc] - r_logFC[1]) / (r_logFC[2] - r_logFC[1])
                logs <- c(logs, val)
            }
        }
        cols <- c()
        for (ys in 1:length(logs))
            cols <- c(cols, ifelse(data$logFC[ys] > 0, "upregulated", "downregulated"))
        dfp <- data.frame( logx = xm, 
                           logy = logs,
                           logFC = factor(cols),
                           logy2 = rep(rad2, length(logs)))
        c <- ggplot() + 
          geom_rect( data = df, 
                     aes( xmin = x, xmax = x + xmax, ymin = y1,
                          ymax = y1 + ymax, fill = zscore), colour = "black") +
            geom_rect(data = df, 
                      aes( xmin = x, xmax = x + xmax, ymin = y2, ymax = y2 + 1),
                      fill = "gray70") + 
          geom_rect( data = df, 
                     aes( xmin = x, xmax = x + xmax, ymin = y2 + 0.5, ymax = y2 + 0.5 ),
                          colour = "white" ) + 
          geom_rect( data = df,
                          aes(xmin = x, xmax = x + xmax, ymin = y2 + 0.25, ymax = y2 + 0.25),
                          colour = "white") + 
          geom_rect(data = df,
                    aes( xmin = x, xmax = x + xmax, ymin = y2 + 0.75, ymax = y2 + 0.75),
                          colour = "white" ) +
          geom_text(data = df, 
                    aes(x = x + (xmax/2), y = y2 + 1.3, label = ID, 
                        angle = 360 - (x = x + (xmax/2))/(10/360)),
                    size = label.size, fontface = label.fontface) +
          coord_polar() + labs(title = title) + 
          ylim(1, rad2 + 1.6) + xlim(0, 10) + 
          GOplot:::theme_blank + 
          scale_fill_gradient2("z-score", space = "Lab", low = zsc.col[3],
                               mid = zsc.col[2], high = zsc.col[1], 
                               guide = guide_colourbar(title.position = "top", title.hjust = 0.5),
                               breaks = c(min(df$zscore), max(df$zscore)),
                               labels = c("decreasing", "increasing")) +
          theme( legend.position = "bottom", 
                 legend.background = element_rect(fill = "transparent"),
                 legend.box = "horizontal", legend.direction = "horizontal") +
            geom_point( data = dfp,
                        aes(x = logx, y = logy2 + logy),
                        pch = 21, fill = "transparent", colour = "black", size = 2)+
          geom_point(data = dfp, 
                     aes(x = logx, y = logy2 + logy, colour = logFC), size = 2.0) +
            scale_colour_manual(values = lfc.col, 
                                guide = guide_legend(title.position = "top", title.hjust = 0.5))
        if (table.legend) {
            table <- GOplot:::draw_table(suby)
            graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
            grid.arrange(c, table, ncol = 2)
        }
        else {
            c + theme(
                plot.background = element_rect(fill = "white"),
                panel.background = element_rect(fill = "white")
            )
        }
}
# Leyenda para los cnet interactivos #########
visnetLegendReport <- function(kggDT = NULL, rows = NULL){
      mydf <- data.frame(id = rep(1, 100), sales = 1:100)
    minVal <- format( min( kggDT$`p-value`[rows] ), scientific = T, digits = 2)
    maxVal <- format( max( kggDT$`p-value`[rows] ), scientific = T, digits = 2)
    p <- ggplot(mydf) +
      geom_tile(aes(x = 1, y=sales, fill = sales), show.legend=FALSE) +
      scale_x_continuous(limits=c(0.5,1.5),breaks=1)+
      scale_y_continuous(breaks = c(1,100), 
                         labels = c(minVal,maxVal), position = "right")+
      scale_fill_gradient(high = "blue", low = "red") +
      theme_void() +
      theme(axis.text.y.right = element_text(hjust = 0, size = 12, colour = "#cdcdcd"),
            plot.background = element_rect(fill= "#2d3741", color = NA ) )
    return(p)
}
# popUpModal para report ######################
    popupModal <- function() {
      modalDialog(
          title = "Report configuration",
          size = "l",
          fluidRow(column(width=11,
          tabsetPanel(
              tabPanel("Preview",
                       checkboxGroupButtons(
                           size="sm",
                           individual = TRUE,
                           inputId = "modalPreview",
                           label = "Select preview elements to report",
                           choices = c("PCA", "BoxPlot", "Heatmap", "Cluster","Top6",
                                       "Top1", "Karyoplot","Volcano","MA"),
                           status = "primary",
                           checkIcon = list(
                               yes = icon("ok",
                                          lib = "glyphicon"),
                               no = icon("remove",
                                         lib = "glyphicon")
                           )
                       )
                    ),
              tabPanel("Kegg",
                       checkboxGroupButtons(
                           size = "sm",
                           individual = TRUE,
                           inputId = "modalkeggAll",
                           label = "Select elements to report Kegg All",
                           choices = c("Table", "Barplot", "Chorplot", "Dotplot", "Heatmap", "Netplot"),
                           status = "primary",
                           checkIcon = list(
                               yes = icon("ok",
                                          lib = "glyphicon"),
                               no = icon("remove",
                                         lib = "glyphicon")
                           )
                       ),
                       checkboxGroupButtons(
                           size = "sm",
                           individual = TRUE,
                           inputId = "modalkeggUp",
                           label = "Select elements to report Kegg Up",
                           choices = c("Table", "Barplot", "Chorplot", "Dotplot", "Heatmap", "Netplot"),
                           status = "primary",
                           checkIcon = list(
                               yes = icon("ok",
                                          lib = "glyphicon"),
                               no = icon("remove",
                                         lib = "glyphicon")
                           )
                       ),
                       checkboxGroupButtons(
                           size = "sm",
                           individual = TRUE,
                           inputId = "modalkeggDown",
                           label = "Select elements to report Kegg Down",
                           choices = c("Table", "Barplot", "Chorplot", "Dotplot", "Heatmap", "Netplot"),
                           status = "primary",
                           checkIcon = list(
                               yes = icon("ok",
                                          lib = "glyphicon"),
                               no = icon("remove",
                                         lib = "glyphicon")
                           )
                       )
                       ), # fin tabpanel KEGG
              tabPanel("GO",
                    checkboxGroupButtons(
                           size = "sm",
                           individual = TRUE,
                           inputId = "modalGOAll",
                           label = "Select elements to report GO All",
                           choices = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
                           status = "primary",
                           checkIcon = list(
                               yes = icon("ok",
                                          lib = "glyphicon"),
                               no = icon("remove",
                                         lib = "glyphicon")
                           )
                       ),
                    checkboxGroupButtons(
                           size = "sm",
                           individual = TRUE,
                           inputId = "modalGOUp",
                           label = "Select elements to report GO Up",
                           choices = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
                           status = "primary",
                           checkIcon = list(
                               yes = icon("ok",
                                          lib = "glyphicon"),
                               no = icon("remove",
                                         lib = "glyphicon")
                           )
                       ),
                    checkboxGroupButtons(
                           size = "sm",
                           individual = TRUE,
                           inputId = "modalGODown",
                           label = "Select elements to report GO Down",
                           choices = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
                           status = "primary",
                           checkIcon = list(
                               yes = icon("ok",
                                          lib = "glyphicon"),
                               no = icon("remove",
                                         lib = "glyphicon")
                           )
                       )
              ), #fin tabpanel GO
              tabPanel("GSEA",
                  checkboxGroupButtons(
                           size = "sm",
                           individual = TRUE,
                           inputId = "modalGSEA",
                           label = "Select elements to report GSEA",
                           choices = c("Table", "GSEA plot"),
                           status = "primary",
                           checkIcon = list(
                               yes = icon("ok",
                                          lib = "glyphicon"),
                               no = icon("remove",
                                         lib = "glyphicon")
                           )
                       )
              ) #fin de tabpanel GSEA
                       ) # fin tabsetpanel
          )
      ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("ok", "Apply"),
          uiOutput("downloadhtml")
        )
      )
    }
## myggwordcloud ############################
myggwordcloudReport <- function(data){
  df <- data
  text_df <- tibble( text = paste0( df$Term, collapse = " "), line = 1)
  unigrama <- text_df %>% 
    unnest_tokens(input = text, output = bigram, token = "ngrams", n = 1 )
  counter <- unigrama %>% 
    dplyr::count(bigram, sort = TRUE)
  letras <- c(letters, LETTERS)
  bigram_filter <- counter %>% 
    filter(!bigram %in% stop_words$word) %>% 
    filter(!bigram %in% letras)
  bigram_filter <- bigram_filter[ - which(!is.na(extract_numeric(bigram_filter$bigram))  ) ,]
  wordcloud::wordcloud(bigram_filter$bigram, bigram_filter$n, random.order = F, random.color = T,
                       min.freq = 2, max.words = 200, scale = c(6,1),
                       colors = distinctColorPalette(length(unique(bigram_filter$n))) )
}