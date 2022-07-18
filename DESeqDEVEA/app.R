library(AnnotationDbi)
library(chorddiag) #devtools::install_github("mattflor/chorddiag")
library(DESeq2)
library(DT)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Rnorvegicus.v79)
library(fgsea)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)
library(heatmaply)
library(limma)
library(`mychordplot`) #install_github("fpsanz/mychordplot")
#remotes::install_version("RSQLite", version = "2.2.5")
library(org.Hs.eg.db) #Homo sapiens
library(org.Mm.eg.db) #Mus musculus
library(org.Rn.eg.db) #Rattus norvegicus
library(pheatmap)
library(plotly)
library(purrr)
library(randomcoloR)
library(RColorBrewer)
library(rgl)
library(rglwidget)
library(scales)
library(shinyalert)
library(shinyBS)
library(shinybusy)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(shinymanager)
library(stringr)
library(tidyverse)
library(tidytext)
library(visNetwork)
library(wordcloud)
library(scales)
source("utils.R")
source("updatepopModals.R")
options(shiny.maxRequestSize = 3000*1024^2)
  
### HEADER ############ 
header <- dashboardHeader(title = "Go DESeq DEVEA", 
                          titleWidth = 300, 
                          dropdownMenuOutput("messageMenu"),
                          tags$li(class="dropdown", 
                                  actionButton("notesButton","Report notes"),
                                  style="margin-top:8px; margin-right: 5px"),
                          tags$li(class = "dropdown",
                                  actionButton("report2", "HTML report"),
                                  style="margin-top:8px; margin-right: 10px"),
                          tags$li(class="dropdown", 
                                  actionButton("moreinfo","Tutorial",
                                  style = "background-color: #8ec4d9"),
                                  style="margin-top:8px; margin-right: 5px"),
                          tags$li(class = "dropdown",
                                  actionButton("aboutButton", "About",
                                  style = "background-color: #8ec4d9"),
                                  style="margin-top:8px; margin-right: 15px"),
                          tags$li(class = "dropdown", actionBttn(inputId = "resetbutton",
                                                                 label = "Reset App", style="simple", size="sm",
                                                                 color ="danger"),
                                  style="margin-top:8px; margin-right: 10px"
                          )
)
### SIDEBAR ##########
sidebar <- dashboardSidebar(useShinyalert(),
                            useShinyjs(),
                            tags$br(),
                            sidebarMenu(id="menupreview",
                              menuItem("Input data",
                                       tabName = "info",
                                       icon = icon("info"))),
                            # sidebarMenu(
                            #     menuItem(
                            #         pickerInput(
                            #             inputId = "specie",
                            #             label = "1. Select species",
                            #             choices = list( "Human" = "Hs",
                            #                             "Mouse" = "Mm"),
                            #             options = list(title = "species"),
                            #             selected = NULL
                            #         ) 
                            #     )
                            # ),
                            # sidebarMenu(menuItem(uiOutput("matrixDeseq"))),
                            # sidebarMenu(fluidRow(
                            #             column(width=1,menuItem(uiOutput("circleinfoDO"))),
                            #             column(width=11,menuItem(uiOutput("tooltipDO")),
                            #             menuItem(uiOutput("deseqFile"))))),
                            # sidebarMenu(fluidRow(
                            #             column(width=1,menuItem(uiOutput("circleinfoCM"))),
                            #             column(width=11,menuItem(uiOutput("tooltipCM")),
                            #             menuItem(uiOutput("countFile"))))),
                            # sidebarMenu(fluidRow(
                            #             column(width = 11,
                            #               menuItem(uiOutput("sampleFile"))))),
                            #sidebarMenu(id="menuimport",sidebarMenuOutput("importvw")),
                            sidebarMenu(id = "previewMenu", sidebarMenuOutput("prevw")),
                            sidebarMenu("", sidebarMenuOutput("menu")),
                            # tags$br(),
                            # sidebarMenu(menuItem(uiOutput("design"))),
                            tags$div(
                            tags$a(href='https://jacob.cea.fr/drf/ifrancoisjacob/Pages/Departements/MIRCen/themes/astrocytes-reactifs-biomarqueurs-imagerie-cibles-therapeutiques.aspx', target="_blank",
                                   tags$img(src='mircen.png',width='49%',
                                            style="padding: 5px; position: absolute; bottom:0px; left:0") ),
                            tags$a(href='http://www.bioinformatica.imib.es', target="_blank",
                                   tags$img(src='IMIB_color_gris.svg',width='51%',
                                            style="padding: 5px; float: right; bottom:5px;") ),
                            tags$a(href='https://jacob.cea.fr/drf/ifrancoisjacob/Pages/Departements/CNRGH/LABORATOIRES/Bio-analyse.aspx', target="_blank",
                                   tags$img(src='cnrgh.png',width='49%',
                                            style="padding: 5px; float: left; bottom:15px;") ),
                            
                            style = "position: absolute; bottom:0;width:100%;")
                            
)
                      

### BODY ###############
body <- dashboardBody(
    tags$script(HTML("$('body').addClass('fixed');")),
    add_busy_gif(src="dna-mini.gif", position = "full-page", width = 10, height = 10 ),
     # tags$head(
     #   tags$link(rel = "stylesheet", type = "text/css", href = "customDark.css")
     # ),
    includeCSS("./www/customDark.css"),
    setShadow(class = "shiny-plot-output"),
    setShadow( class = "box"),
    setShadow( class = "svg-container"),
    bsAlert("alert"),
    tabItems(
    # Initial INFO
      tabItem( tabName = "info",
               source(file = "ui-info-tab.R",
                      local = TRUE,
                      encoding = "UTF-8"
                )$value),
      # import tab ######
      # tabItem(tabName = "tabimport",
      #         source(file = "ui-import-tab.R",
      #                local=TRUE,
      #                encoding = "UTF-8"
      #         )$value),
      # preview tab ######
      tabItem(tabName = "preview",
              source(file = "ui-preview-tab.R",
                     local=TRUE,
                     encoding = "UTF-8"
              )$value),
      # kegg tab content #####
      tabItem(tabName = "kegg",
              source(file = "ui-kegg-tab.R",
                     local = TRUE,
                     encoding = "UTF-8"
              )$value),
      # GO tab GO tab ######
      tabItem( tabName = "go",
               source(file = "ui-go-tab.R",
                      local = TRUE,
                      encoding = "UTF-8",
               )$value),
      # GSEA tab ######
      tabItem( tabName = "gsea",
               source(file = "ui-gsea-tab.R",
                      local = TRUE,
                      encoding = "UTF-8",
               )$value)
  ), # fin tab items
    bsModal("modalNotes", "Notes", "notesButton",size="large",
          textAreaInput("textNotes", "Comments", width = "850px", height = "575px", placeholder = "
          Input comments about your analysis here will be incorporated into the report. \n
          Markdown notation is accepted. 
          To avoid creating first level sections, please use double ## or triple ### to define subsections instead of a single #.\n
          Example:\n
          ## Enrichment
          Some notes for enrichment\n
          ### Subsection GO
          Some notes for GO\n
          ### Subsection KEGG
          Some notes for KEGG" ))
) # fin dashboardbody

########################################## UI #################################################

ui <- dashboardPage(title="Go DESeq DEVEA",
                    header,
                    sidebar,
                    body
) # fin del UI

# ui <- secure_app(ui, enable_admin = TRUE, theme = shinythemes::shinytheme("darkly"),
#                  head_auth = HTML("<style>
#                  .panel-auth{
#                                   background-color: #343e48 !important;
#                                   }
#                                   </style>"
#                                   ),
#                  tags_bottom = tagList(tags$div(style = "text-align: center;",
#                    tags$image(
#                      height = 40,
#                      src = "mircen.png",
#                      style = "padding-right: 10px; padding-top: 10px;"
#                    ),
#                    tags$image(
#                      height = 50,
#                      src = "imib.png"#,
#                    ))
#                   )
#                  )
########################################## SERVER #################################################
server <- function(input, output, session) {
  
  # res_auth <- secure_server(
  #   timeout = 0,
  #   check_credentials = check_credentials(
  #       "~/.users/users.sqlite",
  #       passphrase = readRDS("~/.users/dbpass.Rds")
  #   )
  # )

  observeEvent(input$aboutButton, {
    shinyalert("DEVEA app 2022", HTML("Authors:<br>
    Miriam Riquelme Pérez 
    <a href='https://www.linkedin.com/in/miriam-riquelme-perez/' target='_blank'> 
    <img src='linkedin_little.svg'> </a> <a href='mailto:miriam.riquelmep@gmail.com'>
    <img src='email.svg'></a><br>
    Fernando Pérez Sanz 
    <a href='https://www.linkedin.com/in/fernandoperez72/' target='_blank'> 
    <img src='linkedin_little.svg'> 
    </a> <a href='mailto:fernando.perez@ffis.es'> <img src='email.svg'></a><br>
    For any suggestion or bug, please contact us"),
               imageUrl = "dna-svg-small-13.gif", 
               imageWidth = 200, imageHeight = 100, html=TRUE)})

  observeEvent(input$resetbutton,{
    session$reload()
  })
  
  # observeEvent(input$moreinfo,{
  #   showModal(
  #     modalDialog(
  #       size="l",
  #       tags$iframe(src="tutorial.html",  width="850px", height="700px")
  #     )
  #   )
  # })
  shinyjs::onclick("moreinfo", runjs("window.open('tutorial.html','_blank')") )
  
# Definir reactiveVariables globales ##############
  coloresPCA <- reactiveValues(niveles=NULL, numNiveles=NULL)
  countdata <- reactiveValues() # para convertir count matrix y sample data
  countfile <- reactiveValues() # para leer count matrix y sample data
  data <- reactiveValues() # genes
  datos <- reactiveValues(dds=NULL) #objetos dds post DESeq()
  fcRange <- reactiveValues() # min y max fc
  go <- reactiveValues() # enrich GO
  goDT <- reactiveValues() #pretabla GO
  gsea <- reactiveValues() # objeto GSEA
  kgg <- reactiveValues() # enrich kegg
  kggDT <- reactiveValues() # pretabla kegg
  logfcRange <- reactiveValues() # min y max logfc
  numgenesDE <- reactiveValues(up=NULL, down=NULL)
  res <- reactiveValues()
  rlog <- reactiveValues()
  validateCountData <- reactiveValues(ok=FALSE) #para validar count y sample ok
  vals <- reactiveValues()
  vsd <- reactiveValues()
  svg <- reactiveValues()
  padjNA <- reactiveValues()
  conversion <- reactiveValues()
  
  observeEvent(input$deseqFile, {
    datos$dds <- readRDS(input$deseqFile$datapath)
  })
  
  observeEvent(datos$dds,{
        validate(need(datos$dds, ""))
          if(!is(datos$dds, "DESeqDataSet") | !("results" %in% mcols(mcols(datos$dds))$type) ){
          createAlert(session, "alert", "fileAlert",title = "Oops!!", 
          content = "File must be a DESeqDataSet class object
           and you should have run DESeq()", append=FALSE, style = "error")}
  })
  coloresPCA$colores <- reactive({
        tmp = NULL
      for(i in seq_len(coloresPCA$numNiveles)){
          tmp <- c(tmp, input[[ coloresPCA$niveles[i] ]] )
      }
      return(tmp)
  })

  # Acciones al pulsar boton generar DESeq ###################
  # observeEvent(input$applyTest, {
  #   countdata$sample <- countdata$sample %>% mutate_at(.vars = vars(-1), ~as.factor(.) )
  #   deseqObj <- DESeqDataSetFromMatrix(countData = countdata$count, 
  #                                       colData = DataFrame(countdata$sample),
  #                                       design = ~1)
  #   design(deseqObj) <- as.formula(paste("~", input$testVariablePicker ))
  #   if(input$testAlgorithmPicker == "wald"){
  #     datos$dds <- DESeq(deseqObj, test = "Wald")
  #   }
  #   if(input$testAlgorithmPicker == "lrt"){
  #     datos$dds <- DESeq(deseqObj, test = "LRT", reduced = ~1, parallel = TRUE)
  #     }
  # })
  observeEvent(input$testVariablePicker, {
    if(input$testVariablePicker!=""){
        countdata$sample <- countdata$sample %>% mutate_at(.vars = vars(-1), ~as.factor(.) )
        deseqObj <- DESeqDataSetFromMatrix(countData = countdata$count, 
                                           colData = DataFrame(countdata$sample),
                                           design = ~1)
        design(deseqObj) <- as.formula(paste("~", input$testVariablePicker ))
          #    if(input$testAlgorithmPicker == "wald"){
          datos$dds <- DESeq(deseqObj, test = "Wald")
      }
    #    }
    # if(input$testAlgorithmPicker == "lrt"){
    #   datos$dds <- DESeq(deseqObj, test = "LRT", reduced = ~1, parallel = TRUE)
    # }
  })
  # Acciones al cargar fichero deseq ##########################
  observeEvent(design(), {
    validate(need(design(), ""))
        closeAlert(session, "alert")
          colData(datos$dds)@listData <- colData(datos$dds)@listData %>%
              as.data.frame() %>% mutate_at(vars(-sizeFactor, contains('replaceable')), as.character) %>% ##aqui##
              mutate_at(vars(-sizeFactor, contains('replaceable')), as.factor) %>% as.list()
        #res$sh <- as.data.frame(lfcShrink(datos$dds, coef=(as.numeric(design())+1), type="apeglm", parallel = TRUE))
        res$sh <- as.data.frame(results(datos$dds, contrast = list(resultsNames(datos$dds)[as.numeric(design())+1] ))) #09/02/2020
        res$sh <- res$sh %>% dplyr::select(-stat) #09/02/2020
        conversion$ids <- geneIdConverter2(rownames(res$sh), specie() )
        padjNA$true <- which(is.na(res$sh$padj)) 
        if(length(padjNA$true)!=0 ){
          conversionRes <- conversion$ids[-padjNA$true,]
          res$sh <- res$sh[-padjNA$true,]
        }else{
          conversionRes <- conversion$ids
        }
        res$sh$baseMean <- round(res$sh$baseMean,4)
        res$sh$lfcSE <- round(res$sh$lfcSE,4)
        res$sh$log2FoldChange <- round(res$sh$log2FoldChange,4)
        res$sh <- cbind(`Description`=conversionRes$description, res$sh)
        res$sh <- cbind(`ENTREZ`=conversionRes$ENTREZID, res$sh)
        res$sh <- cbind(`ENSEMBL` = conversionRes$ENSEMBL, res$sh)
        #res$sh <- cbind(`GeneName_Symbol`=conversion$consensus, res$sh)
        res$sh <- cbind(`GeneName_Symbol`=conversionRes$SYMBOL, res$sh) #
        res$sh <-  res$sh %>% dplyr::select(-c(pvalue))
        if(specie() == "Mm" ){spc = "Mus_musculus"}
        else if(specie() == "Hs") {spc = "Homo_sapiens"}
        else{spc = "Rattus_norvegicus"}
        links = paste0("<a href='http://www.ensembl.org/",spc,"/Gene/Summary?db=core;g=",
                       rownames(res$sh),"' target='_blank'>",rownames(res$sh),"</a>")
        res$sh <- cbind(`User_GeneId`= links, res$sh)
        res$lostgene <- length(which(is.na(res$sh$ENTREZ)))
        #vsd$data <- vst(datos$dds)
        rlog$datos <- rlog(datos$dds)
        logfcRange$min <- min(res$sh$log2FoldChange)
        logfcRange$max <- max(res$sh$log2FoldChange)
        fcRange$min <- ifelse(logfcRange$min<0, -(2^abs(logfcRange$min)), 2^abs(logfcRange$min))
        fcRange$max <- ifelse(logfcRange$max<0, -(2^abs(logfcRange$max)), 2^abs(logfcRange$max))
        closeAlert(session, "fileAlert")
        updateTabItems(session, "previewMenu", "preview")
  })
  
  # Acciones al pulsar el boton enrich #####################
  observeEvent(input$runEnrich, {
    if( is.null( fileuniverse() ) ){
        bckgnd <- NULL
    }else{ 
        universe <- read.table(fileuniverse()$datapath, header = F)
        bckgnd <- geneIdConverter2( universe[,1], specie = specie() )
    }
    data$genesUp <- getSigUpregulated(res$sh, padj(), logfc()[2], specie() ) 
    data$genesDown <- getSigDownregulated(res$sh, padj(), logfc()[1], specie() ) 
    data$genesall <- rbind(data$genesUp, data$genesDown)
    
    kgg$all <- customKegg(data$genesall, species = specie(), universe = bckgnd$ENTREZID ) 
    kggDT$all <- kegg2DT(kgg$all, data$genesall)
    
    kgg$up <- customKegg(data$genesUp, species = specie(), universe = bckgnd$ENTREZID) 
    kggDT$up <- kegg2DT(kgg$up, data$genesUp)
    
    kgg$down <- customKegg(data$genesDown, species = specie(), universe = bckgnd$ENTREZID) 
    kggDT$down <- kegg2DT(kgg$down, data$genesDown)
    
    go$all <- customGO(data$genesall, species = specie(), universe = bckgnd$ENTREZID )
    goDT$all <- go2DT(enrichdf = go$all, data = data$genesall )
    
    go$up <- customGO(data$genesUp, species = specie(), universe = bckgnd$ENTREZID )
    goDT$up <- go2DT(enrichdf = go$up, data = data$genesUp )
    
    go$down <- customGO(data$genesDown, species = specie(), universe = bckgnd$ENTREZID )
    goDT$down <- go2DT(enrichdf = go$down, data = data$genesDown )
    updateTabItems(session, "previewMenu", "preview")
  })
  
  # Acciones al pulsar applyButton ################
  padjVal <- reactiveValues(val=0.05)
  fcVal <- reactiveValues( val=c(-1.5, 1.5) )
  logfcVal <- reactiveValues(val=c(-0.5,0.5))
  padj <- reactive({padjVal$val})
  logfc <- reactive({logfcVal$val})
  fc <- reactive({fcVal$val})
  
  observeEvent(input$applyParam,{
      padjVal$val <- input$padj
      if( isTRUE( fc_switch()) ){
        logfcTmp <- vector()
        fcVal$val <-input$fc
        logfcTmp[1] <- ifelse(fc()[1]<0, -log2(abs(fc()[1])), log2(abs(fc()[1])) )
        logfcTmp[2] <- ifelse(fc()[2]<0, -log2(abs(fc()[2])), log2(abs(fc()[2])) )
      } else {
        logfcTmp <- input$logfc
      }
      if(logfcTmp[1]==logfcTmp[2]){
          logfcVal$val <- c(0,0)
          }else{
            logfcVal$val <- logfcTmp
          }
    })
  
  # Acciones al seleccionar variables ################
  observeEvent(input$variables, {
        if(!is.factor(colData(datos$dds)[[ variables()[1] ]] ) ){
          coloresPCA$niveles <- as.character(unique( colData(datos$dds)[[ variables()[1] ]] ))
      } else {
          coloresPCA$niveles <- as.character(levels(colData(datos$dds)[[ variables()[1] ]]))
      }
      coloresPCA$numNiveles <- length(coloresPCA$niveles)
  })
  
  # generate reactive variable ###################
  rowsAll <- reactive({input$tableAll_rows_selected})
  rowsUp <- reactive({input$table_rows_selected})
  rowsdown <- reactive({input$tableDown_rows_selected})
  
  bprowsall <- reactive({input$tableBPall_rows_selected}) 
  mfrowsall <- reactive({input$tableMFall_rows_selected})
  ccrowsall <- reactive({input$tableCCall_rows_selected})
  
  bprowsup <- reactive({input$tableBP_rows_selected})
  mfrowsup <- reactive({input$tableMF_rows_selected})
  ccrowsup <- reactive({input$tableCC_rows_selected})
  
  bprowsdown <- reactive({input$tableBPdown_rows_selected})
  mfrowsdown <- reactive({input$tableMFdown_rows_selected})
  ccrowsdown <- reactive({input$tableCCdown_rows_selected})
  
  variables <- reactive({input$variables})
  genesVolcano <- reactive({input$genesVolcano})
  samplename <- reactive({input$samplename})
  gsearow <- reactive({input$gseaTable_rows_selected})
  specie <- reactive({input$specie})
  
  biologicalText <- reactive({input$biologicalText})
  explainPreview <- reactive({input$explainPreview})
  keggAllText <- reactive({input$keggAllText})
  GSEAText <- reactive({input$GSEAText})
  numheatmap <- reactive({input$numheatmap})
  gene <- reactive({
      if(!is.null(input$genetop1)){input$genetop1}else{
          as.character(res$sh$GeneName_Symbol[ which(!( res$sh$padj>padj() &
                                                            res$sh$log2FoldChange>logfc()[1] &
                                                            res$sh$log2FoldChange<logfc()[2] )) ])[1]
      }
  })
  
  typeBarKeggAll <- reactive({input$selectkeggall})
  typeBarBpAll <- reactive({input$selectbpall})
  typeBarMfAll <- reactive({input$selectmfall})
  typeBarCcAll <- reactive({input$selectccall})
  pca3d <- reactive({input$pca3d})
  boxplotswitch <- reactive({input$boxplotswitch})
  design <- reactive({input$designPicker})
  fc_switch <- reactive({input$fc_switch})
  
  fileuniverse <- reactive({input$fileuniverse})
  universe <- reactive({input$universe})
  
  # Matrix or DESeq ##############
  output$matrixDeseq <- renderUI({
      validate(need(specie(), ""))
      pickerInput("matrixDeseq",
          label = "2. Select input mode",
          choices = list("Count matrix"="cm", "DESeq object"="do"),
          selected = NULL,
          options = list(title = "mode")
          )
  })
  # ...... Count/sample File .......############
  
  # Render infoCM
  output$infoCM <- renderUI({
    validate(need(input$deseqFile, "")) 
    })
  
  # Render infoDO
  output$infoCM <- renderUI({
    validate(need(input$countFile,"")) 
  })
  
  
  # Tooltip de CountMAtrix ##################
  output$circleinfoCM <- renderUI({
    validate(need(input$matrixDeseq =="cm", ""))
    circleButton(
      inputId = "infoCM",
      icon = icon("info"),
      size = "xs",
      status = "primary"
    ) 
  })
  output$tooltipCM <- renderUI({
    validate(need(input$matrixDeseq =="cm", ""))
    bsTooltip(
      "infoCM",
      paste0("The accepted formats are .txt, .tsv, .xlsx"),
      trigger = "hover",
      placement = "right"
    )
  })
   # Render fileInput ########
      output$countFile <- renderUI({
      validate(need(input$matrixDeseq =="cm", ""))
      fileInput("countFile",
          "3.1. Choose file with raw counts",
          placeholder = "Counts",
          accept = c(".txt", ".tsv", ".xlsx") )
  })
     # Leer fichero count data #######
        observeEvent(input$countFile, {
        if(grepl("spreadsheet",input$countFile$type)){
            countfile$count <- readxl::read_xlsx(input$countFile$datapath)
            countdata$count <- countfile$count %>% data.frame(., row.names = 1)
        }else{
            countfile$count <- data.table::fread(input$countFile$datapath, nThread = 4,
                              header = TRUE,
                              )
            countdata$count <- countfile$count %>% data.frame(., row.names = 1)
        }
  })
    # Render sampleInput ########
      output$sampleFile <- renderUI({
      validate(need(input$matrixDeseq =="cm", ""))  
      fileInput("sampleFile",
          "3.2. Choose file with sample information",
          placeholder = "Sample",
          accept = c(".txt", ".tsv", ".xlsx") )
  })
     # Leer fichero sample data ########
        observeEvent(input$sampleFile, {
        if(grepl("spreadsheet",input$sampleFile$type)){
            countfile$sample <- readxl::read_xlsx(input$sampleFile$datapath)
            countdata$sample <- countfile$sample %>% data.frame()
            rownames(data$sample) <- data$sample[,1]
        }else{
            file$sample <- data.table::fread(input$sampleFile$datapath, nThread = 4,
                              header = TRUE,
                              )
            countdata$sample <- countfile$sample %>% data.frame()
            rownames(countdata$sample) <- countdata$sample[,1]
        }
  })
    observe({
      validate(need(input$sampleFile, ""),
               need(input$countFile, ""))
      samplesCount <- sort(colnames(countdata$count))
      samplesSample <- sort(countdata$sample[,1])
      if( is_empty(which(samplesSample != samplesCount ) )){
        countdata$count <- countdata$count %>% dplyr::select(countdata$sample[,1])
        validateCountData$ok = TRUE
      } else {
      shinyalert("Sorry!!", "At least one sample name is inconsistent between the two tables", type = "error")
        validateCountData$ok =FALSE
      }
    })
  ######........###############################
  # tooltip DEseqFile ############
    output$circleinfoDO <- renderUI({
      validate(need(input$matrixDeseq =="do", ""))
       circleButton(inputId = "infoDO",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                    )
    })
    output$tooltipDO <- renderUI({
      validate(need(input$matrixDeseq =="do", ""))
      bsTooltip(
      "infoDO",
      paste0("The file must be compressed as .RDS"),
      trigger = "hover",
      placement = "right"
      )
    })
  # InputFile #################
  output$deseqFile <- renderUI({
      validate(need(input$matrixDeseq =="do", ""))
      fileInput("deseqFile",
          "3. Choose DESeq object",
          placeholder = "RDS DESeq",
          accept = ".Rds")
  })
  
  # InputDesign de objeto deseq ###########################
  output$design <- renderUI({
        validate(need(datos$dds,""), 
                 need(input$matrixDeseq=="do", "") )
          opciones <- as.list(seq_len(length(resultsNames(datos$dds)[-1] )))
          names(opciones) <- resultsNames(datos$dds)[-1]
          pickerInput(
          inputId = "designPicker",
          label = "4. Select design",
          choices = opciones,
          options = list(title = "Design"),
          selected = NULL
          ) 
          })
  
  # InputDesign de countMatrix y countData
    output$designMatrix <- renderUI({
      validate(need(datos$dds,""),
               need(countdata$sample,""),
               need(countdata$count,""))
      opciones <- as.list(seq_len(length(resultsNames(datos$dds)[-1] )))
      names(opciones) <- resultsNames(datos$dds)[-1]
      pickerInput(
        inputId = "designPicker",
        label = "5. Select design",
        choices = opciones,
        options = list(title = "Design"),
        selected = NULL
      )
    })
    
    
  # side bar menu ####################
  output$menu <- renderMenu({
      validate(need(kgg$all, ""))
      sidebarMenu(
          menuItem(
              "Kegg Enrichment",
              tabName = "kegg",
              icon = icon("chart-bar")
          ),
          menuItem(
              "GO Enrichment",
              tabName = "go",
              icon = icon("chart-bar")
          ),
          menuItem("GSEA",
                   tabName = "gsea",
                   icon = icon("chart-line"))
      )
      })
  
  output$prevw <- renderMenu({
    validate(need(res$sh, ""))
    sidebarMenu(
      menuItem("Preview dataset",
               tabName = "preview",
               icon = icon("eye"))
    )
  })
  
  # output$importvw <- renderMenu({
  #   validate(need(countdata$sample,""),
  #            need(countdata$count,""))
  #   sidebarMenu(
  #   menuItem("5. Import view",
  #            tabName = "tabimport",
  #            icon = icon("file-import")))
  # })
  # ui selector sample groups ###################
  output$sampleGroup <- renderUI({
    validate(need(datos$dds, ""))
    nvars <- colData(datos$dds) %>% 
      as.data.frame() %>% 
      dplyr::select(-any_of(c("sizeFactor", "replaceable"))) %>% 
      names()
    def_nvar <- nvars[which(nvars %in% as.character(datos$dds@design)[2] ) ]
    selectInput("variables", label="Select condition[s] of interest to highlight",
                choices = nvars,
                selected = def_nvar,
                multiple = TRUE)
  })
  # ........................####
  # Variable Selector ###########
    output$testVariable <- renderUI({
        validate(need(countdata$sample,""))
        opciones <- as.list(names(countdata$sample))
        pickerInput(
          inputId = "testVariablePicker",
          label = "4. Select column for contrast",
          choices = opciones,
          options = list(title = "Variable"),
          selected = NULL
        ) 
          })
  # Tabla colData ################
  output$coldataTable <- renderDT({
    validate(need(countdata$sample,""))
    countdata$sample %>% datatable(options = list(scrollX=TRUE, scrollY="400px"))
  })
    # Tabla countData ################
  output$expressionTable <- DT::renderDataTable({
    validate(need(countdata$sample,""))
    countdata$count %>% head(10) %>% DT::datatable(options = list(scrollX=TRUE, scrollY="400px"))
  })
  
  # ........................####
  # ui selector de genes para volcano plot #######################
  output$geneSelector <- renderUI({
    validate(need(res$sh, ""),
             need(padj(),""))
    # genes <- as.character(res$sh$GeneName_Symbol[ which(!( res$sh$padj>padj() &
    #                                                          res$sh$log2FoldChange>logfc()[1] &
    #                                                          res$sh$log2FoldChange<logfc()[2] )) ])
    genes <- as.character(res$sh$GeneName_Symbol)
    selectInput("genesVolcano", label="Select gene[s] to include label",
                choices = genes,
                multiple = TRUE)
  })
  # ui selector sample name ###################
  output$samplesName <- renderUI({
    validate(need(datos$dds, ""))
    nvars <- colData(datos$dds) %>% 
      as.data.frame() %>% 
      dplyr::select(-any_of(c("sizeFactor", "replaceable"))) %>% 
      names()
    selectInput("samplename", label="Select column for sample name or gathering",
                choices = nvars,
                multiple = FALSE)
  })
  
  # Deslizador fc/logfc según switch #################
  output$fc_control <- renderUI({
    if(isTRUE(fc_switch())){
      validate(need(datos$dds, ""))
      valmin <- ifelse(input$logfc[1]<0, -2^(abs(input$logfc[1] )), 2^(abs(input$logfc[1])) )
      valmax <- ifelse(input$logfc[2]<0, -2^(abs(input$logfc[2] )), 2^(abs(input$logfc[2])) )
      sliderInput("fc", label = "Select FC range to remove (keeps the tails)",
                  min=round(fcRange$min,3), max=round(fcRange$max, 3),
                  value = c(valmin, valmax), step = 0.1 )
    } else {
      validate(need(datos$dds, ""),
               need(fc(), ""))
        if(is.null(input$fc[1]) ){
          valmin = -0.5
          valmax = 0.5
        } else{
          valmin <- ifelse(input$fc[1]<0, -log2(abs(input$fc[1] )), log2(abs(input$fc[1])) )
          valmax <- ifelse(input$fc[2]<0, -log2(abs(input$fc[2] )), log2(abs(input$fc[2])) )
      }
      sliderInput("logfc", label = "Select logFC range to remove (keeps the tails)",
                min=round(logfcRange$min,3), max=round(logfcRange$max, 3),
                value = c(valmin, valmax), 
                step = 0.1 )
    }
  })
  
  # ui selector padj #################################
  output$padj <- renderUI({
    validate(need(datos$dds,""))
    sliderInput("padj", label = "Select p-adjusted threshold", min = 0, max=0.5,
                value=0.05, step = 0.005 )
  })
  # ui selector Colores para PCA y demás #######################
  output$colorPalettes <- renderUI({
      validate(need(datos$dds, ""),
               need(variables(), ""),
               need(coloresPCA$numNiveles, ""))
      l1 <- rep(1:6, times = coloresPCA$numNiveles / 6 , length.out = coloresPCA$numNiveles)
      l2 <- rep(1:9, each = 6, length.out = coloresPCA$numNiveles)
      selectores <- lapply(seq_len(coloresPCA$numNiveles), function(x){
          spectrumInput(
            inputId = paste0(coloresPCA$niveles[x]),
            label = paste0("Select ",coloresPCA$niveles[x]," color"),
            selected = choices_brewer2[[l1[x]]][l2[x]],
            choices = choices_brewer2,
            width = "50%",
            options = list(`toggle-palette-more-text` = "Show more")
            )
      })
  })
 
  # infoboxes ###############################
  output$allbox <- renderInfoBox({
      validate(need(res$sh, ""),
               need(padj(), ""),
               need(logfc(), ""))
      numall <- nrow( res$sh[ ((res$sh$log2FoldChange >= logfc()[2] |
                                    res$sh$log2FoldChange< logfc()[1]) &
                                   res$sh$padj <= padj() ),] ) 
      infoBox("All DE genes", numall, icon = icon("arrows-alt-v"), color = "light-blue", fill = TRUE)
  })
  output$upbox <- renderInfoBox({
      validate(need(res$sh, ""),
               need(padj(), ""),
               need(logfc(), ""))
      numup <- nrow( res$sh[(res$sh$log2FoldChange >= logfc()[2]) & (res$sh$padj <= padj()), ]) 
      numgenesDE$up <- numup
      infoBox("Upregulated genes", numup, icon = icon("thumbs-up", lib = "glyphicon"), color = "light-blue", fill=TRUE)
  })
  output$downbox <- renderInfoBox({
      validate(need(res$sh, ""),
               need(padj(), ""),
               need(logfc(), ""))
      numdown <- nrow( res$sh[(res$sh$log2FoldChange <= logfc()[1]) & (res$sh$padj <= padj()), ])
      numgenesDE$down <- numdown
      infoBox("Downregulated genes", numdown, icon = icon("thumbs-down", lib = "glyphicon"), color = "light-blue", fill=TRUE)
  })
  
  output$fcdown <- renderUI({
        validate(need(logfcRange$min, ""),
                 need(logfc(),""))
        initMin <- round( logfcRange$min, 2)
        initMax <- round( logfcRange$max, 2)
        if(logfc()[1]>=0){
            fgColor="#6baed6"
            inputColor="white"
            bgColor ="#46505a"
            rotation="clockwise"
            min=0
            max=initMax
            angleOffset = 0
        }else{
            fgColor="#46505a"
            inputColor="white"
            bgColor ="#e6550d"
            rotation="clockwise"
            min=initMin
            max=0
            angleOffset=180
        }
      knobInput(
          inputId = "myKnobdown",
          label = "Lower logFC cutoff",
          readOnly = TRUE,
          value = round(logfc()[1],2),
          min = min,
          max=max,
          rotation=rotation,
          displayPrevious = TRUE,
          fgColor = fgColor,
          inputColor = inputColor,
          bgColor = bgColor,
          width = "80%",
          height = "80%"
      )
  })
  output$fcup <- renderUI({
        validate(need(logfcRange$min, ""),
                 need(logfc(),""))
        initMin <- round( logfcRange$min, 2)
        initMax <- round( logfcRange$max, 2)
        if(logfc()[2]>=0){
            fgColor="#6baed6"
            inputColor="white"
            bgColor ="#46505a" 
            rotation="clockwise"
            min=0
            max=initMax
            angleOffset = 0
        }else{
            fgColor="#46505a"
            inputColor="white"
            bgColor ="#e6550d"
            rotation="clockwise"
            min=initMin
            max=0
            angleOffset=180
        }
      knobInput(
          inputId = "myKnobup",
          label = "Upper LogFC cutoff",
          readOnly = TRUE,
          value = round(logfc()[2], 2),
          min = min,
          max=max,
          rotation=rotation,
          displayPrevious = TRUE,
          fgColor = fgColor,
          inputColor = inputColor,
          bgColor = bgColor,
          width = "80%",
          height = "80%"
      )
  })
  output$pval <- renderUI({
      validate(need(padj(), ""))
      fgColor = "#74c476"
      inputColor = "white"
      bgColor = "#46505a"
      knobInput(
          inputId = "myKnobpval",
          label = "P.adj cutoff",
          readOnly = TRUE,
          value = round(padj(), 2),
          min = 0,
          max = 0.2,
          rotation = "clockwise",
          displayPrevious = TRUE,
          fgColor = fgColor,
          inputColor = inputColor,
          bgColor = bgColor,
          width = "80%",
          height = "80%"
      )
  })

  # preview table ###################
  output$preview <- DT::renderDataTable(server = FALSE,{
    validate(need(datos$dds, "Load file to render table"),
             need(res$sh, "Load file to render table"))
    res.sh <- res$sh
    res.sh <- res.sh[ ((res.sh$log2FoldChange >= logfc()[2] |
                          res.sh$log2FoldChange < logfc()[1]) &
                         res.sh$padj <= padj() ),]
    tituloTabla <- paste0("Table: Expression values | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list(list(extend="csv",filename="expressionValues"), list(extend="excel",filename="expressionValues")),
           text="Download", title=tituloTabla ) )
    res.sh <- res.sh %>% dplyr::select(-lfcSE) 
    res.sh %>% mutate_if(is.character, as.factor) %>% 
      datatable( extensions = "Buttons", escape = FALSE,
               rownames = FALSE,
               colnames = c("User GeneID","SYMBOL","ENSEMBL","ENTREZ","Description","baseMean","log2FoldChange","pAdj"),
               filter = list(position="top", clear=FALSE),
               options = list(order = list(list(7, 'asc')),
                 lengthMenu = list(c(10,25,50,100,-1), c(10,25,50,100,"All")),
                 columnDefs = list(list(orderable = TRUE,
                                        className = "details-control",
                                        targets = 1),
                                   list(className = "dt-right", targets = 1:(ncol(res.sh)-1))
                 ),
                 scrollY = "400px",
                 rowCallback = JS(
                   "function(row, data) {",
                   "for (i = 6; i < 9; i++) {",
                   "if (data[i]>1000 | data[i]<1){",
                   "$('td:eq('+i+')', row).html(data[i].toExponential(3));",
                   "}",
                   "}",
                   "}"),
                 dom = "Bfrtipl",
                 buttons = customButtons,
                 list(pageLength = 10, white_space = "normal")
               )
    )   %>% 
      formatStyle('log2FoldChange',
                      backgroundColor = styleInterval(c(logfc()[1],logfc()[2]) , 
                                                      c(input$downColor,"gray",input$upColor) ) ) %>% 
      formatStyle('baseMean', background = styleColorBar(res.sh$baseMean, "#357E43") )
  })
  
  output$lostgenes <- renderText({
    print( paste0(res$lostgene," out of ", dim(res$sh)[1], " have no ENTREZ ID. These genes will be missing in enrichment analysis"))
  })
  
# preview samples ###################
  output$samples <- DT::renderDataTable(server = FALSE,{
    validate(need(datos$dds, "Load file to render table"))
    metadata <- as.data.frame(colData(datos$dds)) %>% dplyr::select(-any_of(c("sizeFactor", "replaceable")))
    tituloTabla <- paste0("Table: ColData | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list(list(extend="csv",filename="coldata"), list(extend="excel",filename="coldata")),
           text="Download", title=tituloTabla ) )
    
    datatable( metadata, extensions = "Buttons",
               rownames=FALSE,
               filter = list(position="top", clear=FALSE),
               options = list(
                 columnDefs = list(list(orderable = TRUE,
                                        className = "details-control",
                                        targets = 1),
                                   list(className = "dt-right", targets = 1:(ncol(metadata)-1))
                 ),
                 dom = "Bfrtipl",
                 buttons = customButtons,
                 list(pageLength = 10, white_space = "normal")
               )
    )
  })  
# ............ #############
  # view pca plot data ###################
output$pca3 <- renderUI({
      if (!isTRUE(pca3d())) {
          plotlyOutput("pca", width = "100%", height = "800px")
      } else{
          rglwidgetOutput("pca3d", width = "500px", height = "500px")
      }
  })

output$dimensions <- renderUI({
    #validate(need(ndmax(),""))
    if( ncol(assay(rlog$datos) ) < 5 ){
      ndmax=ncol(assay(rlog$datos))
    }else{
      ndmax <- 5
    }
    selectizeInput("ndmax", "Select components to plot", choices = seq_len(ndmax),
                   multiple = TRUE, selected = c(1,2), options = list(maxItems=2))
  })
  
output$pca <- renderPlotly({
    validate(need(!isTRUE(pca3d()), ""),
             need(datos$dds, ""),
             need(variables(), "Select condition to render PCA"),
             need(samplename(), ""),
             need(coloresPCA$colores(), ""),
             need(length(input$ndmax)==2, ""))
    p <- plotPCA(
        rlog$datos,
        intgroup = variables(),
        labels = samplename(),
        customColor = coloresPCA$colores(),
        axes = as.numeric(input$ndmax)
    ) +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
        scale_size_manual(values = 3) +
        theme(text = element_text(size = 16))
    svg$pca <- p
    print(p)
})

output$downPCA <- downloadHandler(
  filename = "pca.svg",
  content = function(file){
    ggsave(file, svg$pca, "svg", units = "in", width = 10, height = 10)}
)

# output$pca3d <- renderRglwidget({
#     validate(need(datos$dds, ""),
#              need(variables(),"Select condition to render PCA" ),
#              need(coloresPCA$colores(), "" ))
#     d <- pca3dplot(rlog$datos, intgroup = variables(), ntop = 500,
#                    returnData = TRUE )
#     levels(d$labels) <- coloresPCA$colores()
#     try(rgl.close(), silent = TRUE)
#     rgl.open(useNULL = TRUE) 
#     x = d$PC1; y = d$PC2; z = d$PC3
#     plot3d(x,y,x, size = 2, type="s", col = (d$labels),
#            box=FALSE, axes=FALSE, xlab = names(d)[1],
#            ylab=names(d)[2], names(d)[3])
#     bg3d(sphere = FALSE, fogtype = "none", color = "#46505a")
#     rgl.lines(c(min(x), max(x)), c(0, 0), c(0, 0), color = "white")
#     rgl.lines(c(0, 0), c(min(y),max(y)), c(0, 0), color = "white")
#     rgl.lines(c(0, 0), c(0, 0), c(min(z),max(z)), color = "white")
#     rglwidget()
#   })
#   
  # view Volcano plot data ###################
   output$volcano <- renderPlot( {
    validate(need(res$sh, "Load file to render plot"))
    res <-  res$sh
    svg$volcano <- CustomVolcano(res, lab = as.character(res$GeneName_Symbol),
                  selectLab = genesVolcano(),
                    x = 'log2FoldChange',
                    y = 'padj',
                    pCutoff = padj(),
                    FCcutoffUP = logfc()[2],
                    FCcutoffDOWN = logfc()[1],
                    drawConnectors=TRUE,
                    #xlim = c(-8, 8),
                    col = c("gray", "#7cccc3", "#d99c01", input$upColor, input$downColor))
    svg$volcano
    })
output$downVolcano <- downloadHandler(
  filename = "volcano.svg",
  content = function(file){
    ggsave(file, svg$volcano, "svg", units = "in", width = 12, height = 7)}
)

xy <- reactive({
  res <- res$sh
  res$`-log10padj` <- (-log10(res$padj)) 
  nearPoints(res, input$plot_click1, xvar = "log2FoldChange", yvar = "-log10padj")
})
output$texto1 <- renderTable( digits = -2, {
        xy <- xy()
        xy[,c(2,5,7,8)]
    })

# view MA plot data ###################
  output$MA <- renderPlot( {
    validate(need(res$sh, "Load file to render plot"),
             need(logfc(), ""))
    svg$maplot <- MA(res$sh, 
                 main = 'MA plot applying the DESeq2 Shrinkage normalization for Foldchange',
       fdr = padj(), fcDOWN = logfc()[1], fcUP = logfc()[2] , size = 1.5,
       palette = c(input$upColor, input$downColor, "gray"),
       genenames = res$sh$GeneName_Symbol,
       legend = "center", top = 15, select.top.method = c('padj','fc'),
       font.label = c("plain", 12),
       font.legend = c("plain", 15),
       font.main = c("plain",22),
       cex.axis = 1.1, cex.lab = 1.3,
       ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5)),
       usergenes = genesVolcano()
    )
    svg$maplot
  })

output$MAdownload <- downloadHandler(
  filename = "ma.svg",
  content = function(file){
    ggsave(file, svg$maplot, "svg", units = "in", width = 12, height = 7)}
)
clicked <- reactive({
  res <- res$sh
  res$`log2(mean - 1)` <- log2(res$baseMean + 1)
  nearPoints(res, input$plot_click2, xvar = "log2(mean - 1)", yvar = "log2FoldChange")
})

output$texto2 <- renderTable( digits = -2, {
  clicked <- clicked()
  clicked[,c(2,5,7,8)]
} )

  # view HEATMAP data ###################
  output$heat <- renderPlotly( {
    validate(need(datos$dds, ""),
             need(rlog$datos, "Load file to render plot"),
             need(variables(),"Load condition to render plot" ),
             need(samplename(),"Load condition to render plot" ) )
    p <- heat2(rlog$datos, n=numheatmap(), intgroup = variables(), sampleName = samplename(),
         specie=specie(), customColor = coloresPCA$colores(), annot=conversion$ids )
    q <- heat2(rlog$datos, n=numheatmap(), intgroup = variables(), sampleName = samplename(),
               specie=specie(), customColor = coloresPCA$colores(), ggplt = TRUE, annot=conversion$ids )
    svg$heat <- q
    print(p)
  })

output$downHeat <- downloadHandler(
  filename = "heat.svg",
  content = function(file){
    ggsave(file, svg$heat, "svg", units = "in", width = 10, height = 10)}
)
  # view CLUSTER data ###################
  output$cluster <- renderPlotly( {
    validate(
      need(datos$dds, ""),
      need(rlog$datos, "Load file to render plot"),
      need(variables(), "Load condition to render plot"),
      need(samplename(), "Load condition to render plot")
    )
    p <- cluster(rlog$datos, intgroup = samplename())
    q <- cluster(rlog$datos, intgroup = samplename(), ggplt=TRUE)
    svg$cluster <- q
    print(p)
  })

output$downCluster <- downloadHandler(
  filename = "cluster.svg",
  content = function(file){
    ggsave(file, svg$cluster, "svg", units = "in", width = 10, height = 10)}
)

# view TOP6 data ###################
  output$top6 <- renderPlotly( {
    validate(need(datos$dds, ""),
             need(res$sh, "Load file to render plot"),
             need(variables(),"Load condition to render plot" ),
             need(coloresPCA$colores(), ""))
    topGenes <- rownames(res$sh)[order(res$sh$padj)][1:6]
    topSymbol <- as.character(res$sh$GeneName_Symbol)[order(res$sh$padj)][1:6]
    z <- lapply(topGenes, function(x) plotCounts(dds=datos$dds, gene=x,
                                                 res=res$sh, intgroup = variables(),
                                                 returnData = TRUE))
    z <- lapply(z, function(x){x %>% group_by(!!as.name(variables()[1])) %>%
        mutate(mean = round(mean(count),2), sem = round(sd(count)/sqrt(n()),3 ), n = n() ) %>% 
        mutate(text = paste0("Mean: ",mean,"\n","SEM: ",sem)) } )
    z <- do.call(rbind, z)
    z$symbol <- rep(topSymbol, each =(nrow(z)/6) ) 
    z[[variables()[1]]] <- as.factor(z[[variables()[1]]])
    p <- ggplot(z, aes_(as.name(variables()), ~count, colour = as.name(variables()[1] ), text = ~text ) ) + 
      scale_y_log10() +
      geom_point(position = position_jitter(width = 0.1, height = 0), size = 2) +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 1))+
      facet_wrap(~symbol) + scale_color_manual( values = coloresPCA$colores() ) +
      ggtitle("Expression of top 6 most significant genes") +
      theme(plot.title = element_text(hjust = 0.5))
    svg$topsix <- p
    p %>% ggplotly(tooltip = c("x","y","text"))
    })

output$downTopsix <- downloadHandler(
  filename = "topsix.svg",
  content = function(file){
    ggsave(file, svg$topsix, "svg", units = "in", width = 10, height = 10)}
)

# ui selector de genes para top1 #######################
output$gene <- renderUI({
  validate(need(res$sh, ""),
           need(padj(),""))
  genes <- as.character(res$sh$GeneName_Symbol[ which(!( res$sh$padj>padj() &
                                                           res$sh$log2FoldChange>logfc()[1] &
                                                           res$sh$log2FoldChange<logfc()[2] )) ])
  selectInput("genetop1", label="Select gene to label",
              choices = genes,
              multiple = FALSE)
})
# view TOP1 data ###################  
  output$top1 <- renderPlotly( {
    validate(need(datos$dds, ""),
             need(res$sh, "Load file to render plot"),
             need(variables(),"Load condition to render plot" ),
             need(coloresPCA$colores(), ""),
             need(gene(), "Enter a gene of interest in Ensembl or symbol name")
             )
    gene <- gene()
    generow <- which(conversion$ids == gene, arr.ind = TRUE)[1,1] #07/02/2021
    gene <- conversion$ids[generow,1] #07/02/2021
    z <- plotCounts(dds = datos$dds,gene = gene, returnData = TRUE,
                    intgroup = variables()[1]) #07/02/2021
    symbol <- conversion$ids$SYMBOL[generow] #07/02/2021
    # if (grepl("^ENS", gene, ignore.case = TRUE)) {
    #     gene <- toupper(gene)
    #     z <- plotCounts(dds = datos$dds,gene = gene, returnData = TRUE,
    #                     intgroup = variables()[1])
    #     symbol = as.character(res$sh$User_GeneId[rownames(res$sh) == gene])
    # } else{
    #     if (specie() == "Mm") {
    #         gene <- stringr::str_to_title(gene)
    #     }
    #     else{
    #         gene = toupper(gene)
    #     }
    #     z <- plotCountsSymbol(dds = datos$dds, gene = gene, returnData = TRUE,
    #                           intgroup = variables()[1], specie=specie())
    #     symbol <- gene
    # }
    z <- z %>% group_by(!!as.name(variables()[1])) %>%
        mutate(mean = round(mean(count),2), sem = round(sd(count)/sqrt(n()),3 ), n = n() ) %>% 
        mutate(text = paste0("Mean: ",mean,"\n","SEM: ",sem))
    p <- z %>% ggplot(aes_(as.name(variables()[1]), ~count, colour = as.name(variables()[1] ),
                           text = ~text) ) + scale_y_log10() +
        geom_point(position = position_jitter(width = 0.1, height = 0), size = 2)+
        scale_color_manual( values = coloresPCA$colores() )+
      theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 1))+
        ggtitle(paste0(symbol) )
    # texto del plot1 #############
    output$top1text <- renderUI({
        validate(need(gene(), ""))
        texto <- as.data.frame(unique(z[ ,c(variables()[1],"text") ] ))
        txt <- paste0(apply(texto, 1, function(x){x} ), collapse = "<br/><br/>")
        txt <- gsub("\n","<br/>",txt)
        gene_fc_padj <- res$sh[which(res$sh==gene, arr.ind = TRUE)[1], c("log2FoldChange","padj")]
        txt <- paste0("Log2FC: ",round(gene_fc_padj[1],2),
                      "<br/>","Padj: ", 
                      format(gene_fc_padj[2], digits=3, scientific=TRUE),
                      "<br/><br/>",txt)
        tags$h5(HTML(txt))
    })
    svg$topone <- p
    p %>% ggplotly(tooltip = c("x","y","text"))
  })

output$downTopone <- downloadHandler(
  filename = "topone.svg",
  content = function(file){
    ggsave(file, svg$topone, "svg", units = "in", width = 10, height = 10)}
)
## karyoplot ######################################
output$karyoPlot <- renderPlot(bg = "#46505a", {
    validate(need(res$sh, "Load file to render plot"))
    krtp(res$sh, specie = specie(), pval = padj(), fcdown = logfc()[1],
         fcup = logfc()[2], bg="#46505a", coldown=input$downColor , colup=input$upColor )
})
output$downKrpt <- downloadHandler(
  filename = "karyoplot.png",
  content = function(file){
    png(file)
    krtp(res$sh, specie = specie(), pval = padj(), fcdown = logfc()[1],
         fcup = logfc()[2], bg="#46505a", coldown= input$downColor , colup=input$upColor )
    dev.off()
  }
)
 # Boxviolin plot #################################
  output$boxviolin <- renderPlotly({
          validate(need(datos$dds, "Load file and condition to render Volcano"),
                   need(variables(),"Load condition to render plot" ),
                   need(rlog$datos, ""),
                   need(variables(), ""),
                   need(samplename(),"" ),
                   need(coloresPCA$colores(), ""))
          p <- boxViolin( names = samplename() , vsd=rlog$datos, boxplotswitch=boxplotswitch(),
                    intgroup=variables(), customColor = coloresPCA$colores() ) 
          p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 1))
          svg$boxviolin <- p
          print(p)
  })
output$downViolin <- downloadHandler(
  filename = "box_violin.svg",
  content = function(file){
    ggsave(file, svg$boxviolin, "svg", units = "in", width = 10, height = 10)}
)
# ............ ###############################
# Funcion tamaño plot ########################
myHeightfunction <- function(filas) {
  if (length(filas) <= 10) {
    return(400)
  } else if (length(filas) > 10 & length(filas) <= 35) {
    return(600)
  } else{
    return(800)
  }
  }
# KEGG table All #####################################
  output$tableAll <- DT::renderDT(server=FALSE,{
    validate(need(kgg$all, "Load file to render table"))
    names(kggDT$all)[names(kggDT$all) == "DE"] <- "DEG"
    names(kggDT$all)[names(kggDT$all) == "P.DE"] <- "p-value"
    tituloTabla <- paste0("Table: Kegg all genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="keggAll"),list(extend="excel",filename="keggAll") ),
           text="Download", title=tituloTabla ) )
      
    datatable2(
      kggDT$all, 
      vars = c("genes"),
      filter = list(position="top", clear=FALSE),
      escape = FALSE,
      opts = list(order = list(list(5, 'asc')),
        pageLength = 10, white_space = "normal",
        scrollY = "400px",
        buttons = customButtons
        ) )
  }) 

  tableallProxy <- dataTableProxy("tableAll")
  observeEvent(input$resettableall, {
    tableallProxy %>% selectRows(NULL)
  })
# KEGG barplot All ################
  output$keggPlotAll <- renderPlotly({
    validate(need(kgg$all, "Load file to render BarPlot"),
             need(rowsAll(), "Select the paths of interest to render BarPlot") )
    rowsAll <- rowsAll()
    if(is.null(rowsAll)){
        if( dim(kgg$all)[1]<10 ){rowsAll <-  seq_len(nrow(kgg$all)) }
        else{ rowsAll <-  seq_len(10)  }
        }
    p <- plotKeggAll(enrichdf = kgg$all[rowsAll,], nrows = length(rowsAll),
                genesUp = data$genesUp, genesDown = data$genesDown,
                colors = c(input$downColor, input$upColor))
    if (typeBarKeggAll() == "Dodge") {
      plt <- p[[1]]; svg$keggall <- p[[1]]
    } else if (typeBarKeggAll() == "Stack") {
      plt <- p[[2]]; svg$keggall <- p[[2]]
    } else {
      plt <- p[[3]]; svg$keggall <- p[[3]]
    }
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( rowsAll() )
    plt$x$layout$height <- myHeightfunction(rowsAll() )
    plt
  })

output$barKeggAll <- downloadHandler(
  filename = "barkeggall.svg",
  content = function(file){
    ggsave(file, svg$keggall, "svg", width = 10, units = "in") }
)
  # KEGG chordiag plot All ###############
  output$keggChordAll <- renderMychordplot({
    validate(need(kgg$all, "Load file to render ChordPlot"),
             need(length(rowsAll())>3, "Select at least 4 paths of interest to render ChordPlot"))
    rowsAll <- rowsAll()
    if(is.null(rowsAll)){
        if( dim(kgg$all)[1]<10 ){rowsAll <-  seq_len(nrow(kgg$all)) }
        else{ rowsAll <-  seq_len(10)  }
    }
    mychordplot(kgg$all[rowsAll, c("Pathway","genes") ], div="keggChordAll" )
    #chordPlot(kgg$all[rowsAll, ], nRows = length(rowsAll), orderby = "P.DE")
  })

  output$legendChorAll <- renderPlot(bg = "#37414b",{
    validate(need(kgg$all, "Load file to render ChordPlot"),
             need(length(rowsAll())>3, "Select at least 4 paths of interest to render ChordPlot") )
    rowsAll <- rowsAll()
    if(is.null(rowsAll)){
        if( dim(kgg$all)[1]<10 ){rowsAll <-  seq_len(nrow(kgg$all)) }
        else{ rowsAll <-  seq_len(10)  }
        }
    legendChorplot(kgg$all[rowsAll, ] )
  })
  

  # KEGG dotplot All ###################
  output$keggDotAll <- renderPlotly({
    validate(need(kgg$all, "Load file and select to render dotPlot"),
             need(rowsAll(), "Select the paths of interest to render DotPlot"))
    rowsAll <- rowsAll()
    if(is.null(rowsAll)){
      if( dim(kgg$all)[1]<10 ){rowsAll <-  seq_len(nrow(kgg$all)) }
      else{ rowsAll <-  seq_len(10)  }
    }
    plt <- dotPlotkegg(kgg$all[rowsAll,], n = length(rowsAll))
    svg$dotKeggAll <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( rowsAll() )
    plt$x$layout$height <- myHeightfunction(rowsAll() )
    plt
  } )
  
  output$dotkeggAll <- downloadHandler(
    filename = "dotKeggAll.svg",
    content = function(file){
    ggsave(file, svg$dotKeggAll, device = "svg", width = 10, units = "in") }
  )
  
  # KEGG heatmap All #################
  output$heatmapKeggAll <- renderPlotly({
    validate(need(kgg$all, "Load file and select to render Heatmap"),
             need(rowsAll(), "Select the paths of interest to render HeatMap"),
             need(kggDT$all, ""))
    plt <- heatmapKeggLogFC(kggDT$all, res$sh, rowsAll() ) 
    svg$heatKeggAll <- list(kggDT$all, res$sh, rowsAll())
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( rowsAll() )
    plt$x$layout$height <- myHeightfunction(rowsAll() )
    plt
  })
  
    output$heatKeggAll <- downloadHandler(
    filename = "heatKeggAll.svg",
    content = function(file){
    #p <- heatmapKegg( svg$heatKeggAll[[1]], svg$heatKeggAll[[2]] )sustituido 06/09/21
    p <- heatmapKeggLogFC(kggDT$all, res$sh, rowsAll() ) #añadido 06/09/21
    ggsave(filename= file, plot = p, device = "svg", width = 10, units = "in") }
  )
  # KEGG cnet All #################
  output$legend <- renderPlot({
    validate(need(kgg$all, "Load file and select to render Net Plot"),
             need(rowsAll(), "Select the paths of interest to render NetPlot"),
             need(kggDT$all, ""))
    visnetLegend(kggDT = kggDT$all , rows = rowsAll() )
  })

  output$keggAllNet <- renderUI({
    if(!isTRUE( input$keggAllNet_switch ) ){
      plotOutput("cnetKeggAll", height = "600px")
    } else{
      visNetworkOutput("visnetKeggAll", height = "600px")
    }
  })

  output$cnetKeggAll <- renderPlot({
    validate(need(kgg$all, "Load file and select to render Net Plot"),
             need(rowsAll(), "Select the paths of interest to render NetPlot"))
    p <- customCnetKegg(kgg$all, rowsAll(), genesUp = data$genesUp, genesDown = data$genesDown)
    svg$cnetKeggAll <- p
    print(p)
  })

  output$visnetKeggAll <- renderVisNetwork({
    validate(need(kgg$all, "Load file and select to render Net Plot"),
             need(rowsAll(), "Select the paths of interest to render NetPlot"),
             need(kggDT$all, ""))
    visData <- customVisNet(kgg$all, nTerm=rowsAll(), kggDT$all,
                             up = data$genesUp$SYMBOL, down = data$genesDown$SYMBOL )
    visNetwork(visData$nodes, visData$edges, background = "#ffffff") %>%
    visOptions(highlightNearest = list(enabled=TRUE, hover=TRUE),
                nodesIdSelection = TRUE)
  })
  
  output$cnetKggAll <- downloadHandler(
    filename = "cnetKeggAll.svg",
    content = function(file){
    ggsave(filename = file, plot = svg$cnetKeggAll, device = "svg", width = 12, height = 10, units = "in") }
  )
  # ............ ###############################
  # KEGG table up#####################################
  output$table <- DT::renderDataTable(server=FALSE,{
    validate(need(kgg$up, "Load file to render table"))
    names(kggDT$up)[names(kggDT$up) == "DE"] <- "DEG"
    names(kggDT$up)[names(kggDT$up) == "P.DE"] <- "p-value"
    tituloTabla <- paste0("Table: Kegg up-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="keggUp"),list(extend="excel",filename="keggUp") ),
           text="Download", title=tituloTabla ) )

    datatable2(
      kggDT$up,
      vars = c("genes"),
      filter = list(position="top", clear=FALSE),
      escape = FALSE,
      opts = list(order = list(list(5, 'asc')),
        pageLength = 10, white_space = "normal",
        scrollY = "400px",
        buttons = customButtons))
  }) 
  
  tableProxy <- dataTableProxy("table")
  observeEvent(input$resettable, {
    tableProxy %>% selectRows(NULL)
  })
  # KEGG barplot up################
  output$keggPlot <- renderPlotly ({
    validate(need(kgg$up, "Load file to render BarPlot"), 
             need(rowsUp(), "Select the paths of interest to render BarPlot"))
    rowsUp <- rowsUp()
    if(is.null(rowsUp)){
        if( dim(kgg$up)[1]<10 ){rowsUp <-  seq_len(nrow(kgg$up)) }
        else{ rowsUp <-  seq_len(10)  }
        }
      p <- plotKegg(enrichdf = kgg$up[rowsUp,], nrows = length(rowsUp), colors = c(input$upColor))
      svg$barkeggup <- p[[2]] 
      
      plt <- p[[1]]
      plt$height <- myHeightfunction( rowsUp() )
      plt$x$layout$height <- myHeightfunction(rowsUp() )
      plt
  })
  
    output$barKeggUp <- downloadHandler(
    filename = "barkeggup.svg",
    content = function(file){
    ggsave(file, svg$barkeggup, "svg", width = 10, height = 10, units = "in") }
  )
  # KEGG chordiag plot up ###############
  output$keggChord <- renderMychordplot({
    validate(need(kgg$up, "Load file to render ChordPlot"),
             need(length(rowsUp())>3, "Select at least 3 paths of interest to render ChordPlot"))
    rowsUp<- rowsUp()
    if(is.null(rowsUp)){
        if( dim(kgg$up)[1]<10 ){rowsUp <-  seq_len(nrow(kgg$up)) }
        else{ rowsUp <-  seq_len(10)  }
    }
    mychordplot(kgg$up[rowsUp, c("Pathway","genes") ], div="keggChord" )
  })
 output$legendChorUp <- renderPlot(bg = "#37414b",{
    validate(need(kgg$up, "Load file to render ChordPlot"),
             need(length(rowsUp())>3, ""))
    rowsUp <- rowsUp()
    if(is.null(rowsUp)){
        if (dim(kgg$up)[1] < 10) {rowsUp <-  seq_len(nrow(kgg$up))}
        else{rowsUp <-  seq_len(10)}
        }
    legendChorplot(kgg$up[rowsUp, ] )
  })
  # KEGG dotplot UP ################### 
  output$keggDotUp <- renderPlotly({
    validate(need(kgg$up, "Load file and select to render dotPlot"),
             need(rowsUp(), "Select the paths of interest to render DotPlot"))
    rowsUp <- rowsUp()
    if(is.null(rowsUp)){
      if (dim(kgg$up)[1] < 10) {rowsUp <-  seq_len(nrow(kgg$up))}
      else{rowsUp <-  seq_len(10)}
    }
    plt <- dotPlotkegg(kgg$up[rowsUp,], n = length(rowsUp))
    svg$dotKeggUp <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( rowsUp() )
    plt$x$layout$height <- myHeightfunction(rowsUp() )
    plt
  })
  
  output$dotKeggUp <- downloadHandler(
    filename = "dotKeggUp.svg",
    content = function(file){
    ggsave(file, svg$dotKeggUp, device = "svg", width = 10, units = "in") }
  )
  # KEGG heatmap Up #################
  output$heatmapKeggUp <- renderPlotly({
    validate(need(kgg$up, "Load file and select to render Heatmap"),
             need(rowsUp(), "Select the paths of interest to render HeatMap"),
             need(kggDT$up, ""))
    p <- heatmapKeggLogFC(kggDT$up, res$sh, rowsUp())
    svg$heatKeggUp <- list(kggDT$up, res$sh, rowsUp())
    plt <- ggplotly(p)
    plt$height <- myHeightfunction( rowsUp() )
    plt$x$layout$height <- myHeightfunction(rowsUp() )
    plt
  })

  output$heatKeggUp <- downloadHandler(
    filename = "heatKeggUp.svg",
    content = function(file){
    #p <- heatmapKegg(svg$heatKeggUp[[1]],svg$heatKeggUp[[2]] )
    p <- heatmapKeggLogFC(kggDT$up, res$sh, rowsUp() )
    ggsave(filename= file, plot = p, device = "svg", width = 10, units = "in") }
  )
  
  
  
  # KEGG cnet Up #################
   output$keggUpNet <- renderUI({
    if(!isTRUE( input$keggUpNet_switch ) ){
      plotOutput("cnetKeggUp", height = "600px")
    } else{
      visNetworkOutput("visnetKeggUp", height = "600px")
    }
  })
  output$cnetKeggUp <- renderPlot({
    validate(need(kgg$up, "Load file and select to render Net Plot"),
             need(rowsUp(), "Select the paths of interest to render NetPlot"))
    p <- customCnetKegg(kgg$up, rowsUp(), genesUp = data$genesUp, genesDown = data$genesDown)
    svg$cnetKeggUp <- p
    print(p)
  })
  output$visnetKeggUp <- renderVisNetwork({
    validate(need(kgg$up, "Load file and select to render Net Plot"),
             need(rowsUp(), "Select the paths of interest to render NetPlot"),
             need(kggDT$up, ""))
    visData <- customVisNet(kgg$up, nTerm=rowsUp(), kggDT$up,
                             up = data$genesUp$SYMBOL, down = data$genesDown$SYMBOL )
    visNetwork(visData$nodes, visData$edges, background = "#ffffff") %>%
    visOptions(highlightNearest = list(enabled=TRUE, hover=TRUE),
                nodesIdSelection = TRUE)
  })
  output$cnetkeggUp <- downloadHandler(
    filename = "cnetKeggUp.svg",
    content = function(file){
    ggsave(filename = file, plot = svg$cnetKeggUp, device = "svg", width = 12, height = 10, units = "in") }
  )
  # ............ ###############################
  # KEGG table down #####################################
  output$tableDown <- DT::renderDataTable(server=FALSE,{
    validate(need(kgg$down, "Load file to render table"))
    names(kggDT$down)[names(kggDT$down) == "DE"] <- "DEG"
    names(kggDT$down)[names(kggDT$down) == "P.DE"] <- "p-value"
    tituloTabla <- paste0("Table: Kegg down-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="keggDown"),list(extend="excel",filename="keggDown") ),
           text="Download",  title=tituloTabla ) )

    datatable2(
      kggDT$down,
      vars = c("genes"),
      filter = list(position="top", clear=FALSE),
      escape = FALSE,
      opts = list(order = list(list(5, 'asc')),
        pageLength = 10, white_space = "normal",
        scrollY = "400px",
        buttons = customButtons))
  }) 
  
  tableDownProxy <- dataTableProxy("tableDown")
  observeEvent(input$resettableDown, {
    tableDownProxy %>% selectRows(NULL)
  })
  
  # KEGG barplot down ################
  output$keggPlotDown <- renderPlotly ({
    validate(need(kgg$down, "Load file to render BarPlot"),
             need(rowsdown(), "Select the paths of interest to render BarPlot"))
      rowsdown <- rowsdown()
    if(is.null(rowsdown)){
        if( dim(kgg$down)[1]<10 ){rowsdown <-  seq_len(nrow(kgg$down)) }
        else{ rowsdown <-  seq_len(10)  }
        }
    p <- plotKegg(enrichdf = kgg$down[rowsdown,], nrows = length(rowsdown),
                  colors = c(input$downColor))
    svg$barkeggdown <- p[[2]] 
    plt <- p[[1]]
    plt$height <- myHeightfunction( rowsdown() )
    plt$x$layout$height <- myHeightfunction(rowsdown() )
    plt

  })
  
  output$barKeggDown <- downloadHandler(
    filename = "barkeggdown.svg",
    content = function(file){
    ggsave(file, svg$barkeggdown, "svg", width = 10, height = 10, units = "in") }
  )
  # KEGG chordiag plot down ###############
  output$keggChordDown <- renderMychordplot({
    validate(need(kgg$down, "Load file to render ChordPlot"),
             need(length(rowsdown())>3, "Select at least 3 paths of interest to render ChordPlot"))
    rowsdown <- rowsdown()
    if(is.null(rowsdown)){
        if( dim(kgg$down)[1]<10 ){rowsdown <-  seq_len(nrow(kgg$down)) }
        else{ rowsdown <-  seq_len(10)  }
    }
    mychordplot(kgg$down[rowsdown, c("Pathway","genes") ], div="keggChordDown" )
    #chordPlot(kgg$down[rowsdown, ], nRows = length(rowsdown), orderby = "P.DE")
  })
  output$legendChorDown <- renderPlot(bg = "#37414b",{
    validate(need(kgg$down, "Load file to render ChordPlot"),
             need(length(rowsdown())>3, ""))
    rowsdown <- rowsdown()
    if(is.null(rowsdown)){
        if( dim(kgg$down)[1]<10 ){rowsdown <-  seq_len(nrow(kgg$down)) }
        else{ rowsdown <-  seq_len(10)  }
        }
    legendChorplot(kgg$down[rowsdown, ] )
  })
  # KEGG dotplot Down ################### 
  output$keggDotDown <- renderPlotly({
    validate(need(kgg$down, "Load file to render DotPlot"),
             need(rowsdown(), "Select the paths of interest to render dotPlot"))
    rowsdown <- rowsdown()
    if(is.null(rowsdown)){
      if( dim(kgg$down)[1]<10 ){rowsdown <-  seq_len(nrow(kgg$down)) }
      else{ rowsdown <-  seq_len(10)  }
    }
    plt <- dotPlotkegg(kgg$down[rowsdown,], n = length(rowsdown))
    svg$dotKeggDown <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( rowsdown() )
    plt$x$layout$height <- myHeightfunction(rowsdown() )
    plt
  })
 
  output$dotKeggDown <- downloadHandler(
    filename = "dotKeggDown.svg",
    content = function(file){
    ggsave(file, svg$dotKeggDown, device = "svg", width = 10, units = "in") }
  )

  # KEGG heatmap Down #################
  output$heatmapKeggDown <- renderPlotly({
    validate(need(kgg$down, "Load file to render Heatmap"),
             need(rowsdown(), "Select the paths of interest to render Heatmap"))
    p <- heatmapKeggLogFC(kggDT$down, res$sh, rowsdown())
    svg$heatKeggDown <- list(kggDT$down, res$sh, rowsdown())
    plt <- ggplotly(p)
    plt$height <- myHeightfunction( rowsdown() )
    plt$x$layout$height <- myHeightfunction(rowsdown() )
    plt
  })

  output$heatKeggDown <- downloadHandler(
    filename = "heatKeggDown.svg",
    content = function(file){
    #p <- heatmapKegg(svg$heatKeggDown[[1]],svg$heatKeggDown[[2]] )
    p <- heatmapKeggLogFC(kggDT$down, res$sh, rowsdown() )
    ggsave(filename= file, plot = p, device = "svg", width = 10, units = "in") }
  )
    # KEGG cnet Down #################
   output$keggDownNet <- renderUI({
    if(!isTRUE( input$keggDownNet_switch ) ){
      plotOutput("cnetKeggDown", height = "600px")
    } else{
      visNetworkOutput("visnetKeggDown", height = "600px")
    }
  })
  output$cnetKeggDown <- renderPlot({
    validate(need(kgg$down, "Load file and select to render Net Plot"),
             need(rowsdown(), "Select the paths of interest to render NetPlot"))
    p <- customCnetKegg(kgg$down, rowsdown(), genesUp = data$genesUp,
                        genesDown = data$genesDown)
    svg$cnetKeggDown <- p
    print(p)
  })
  output$visnetKeggDown <- renderVisNetwork({
    validate(need(kgg$down, "Load file and select to render Net Plot"),
             need(rowsdown(), "Select the paths of interest to render NetPlot"),
             need(kggDT$down, ""))
    visData <- customVisNet(kgg$down, nTerm=rowsdown(), kggDT$down,
                             up = data$genesUp$SYMBOL, down = data$genesDown$SYMBOL )
    visNetwork(visData$nodes, visData$edges, background = "#ffffff") %>%
    visOptions(highlightNearest = list(enabled=TRUE, hover=TRUE),
                nodesIdSelection = TRUE)
  })
  output$cnetkeggDown <- downloadHandler(
    filename = "cnetKeggDown.svg",
    content = function(file){
    ggsave(filename = file, plot = svg$cnetKeggDown, device = "svg", width = 10, height = 10, units = "in") }
  )
  # ............ ###############################
  # GO table BP ALL #####################
  output$tableBPall <- DT::renderDataTable(server=FALSE,{
    validate(need(goDT$all, "Load file to render table"))
    goDT <- goDT$all 
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level) 
    tituloTabla <- paste0("Table: GO-BP all genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="BPall"),list(extend="excel",filename="BPall") ),
           text="Download", title=tituloTabla ) )
      
    datatable2(goDT[goDT$Ont=="BP",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(order = list(list(6, 'asc')),
                 pageLength = 10, white_space = "normal",
                 scrollY = "400px",
                 buttons = customButtons))
  })
  tableBPallProxy <- dataTableProxy("tableBPall")
  observeEvent(input$resettableBPall, {
    tableBPallProxy %>% selectRows(NULL)
  })
  # GO plots BP all #####################
  output$plotBPall <- renderPlotly({
    validate(need(go$all, "Load file to render plot"),
             need(bprowsall(), "Select at least one row to plot") )
    bprowsall <- bprowsall()
    gosBP <- go$all[go$all$Ont=="BP",]
    if(is.null(bprowsall)){
        if( dim(gosBP)[1]<10 ){bprowsall <-  seq_len(nrow(gosBP)) }
        else{ bprowsall <-  seq_len(10)  }
    }
    p <- plotGOAll(enrichdf = gosBP[bprowsall, ], nrows = length(bprowsall), ont="BP", 
              genesUp = data$genesUp, genesDown = data$genesDown,
              colors = c(input$downColor, input$upColor))
    if (typeBarBpAll() == "Dodge") {
      plt <- p[[1]]; svg$barbpall <- p[[1]]
    }
    else if (typeBarBpAll() == "Stack") {
      plt <- p[[2]]; svg$barbpall <- p[[2]]
    }
    else {
      plt <- p[[3]] ;svg$barbpall <- p[[3]] 
    }
    plt$height <- myHeightfunction( bprowsall() )
    plt$x$layout$height <- myHeightfunction(bprowsall() )
    plt
  })
    output$barBpAll <- downloadHandler(
    filename = "barbpall.svg",
    content = function(file){
      ggsave(file, svg$barbpall, "svg", width = 10, units = "in") }
  )
  # GO BP dotplot all ################### 
    
  output$BPDotall <- renderPlotly({
    validate(need(go$all, "Load file to render dotPlot"),
             need(bprowsall(), "Select the terms of interest to render DotPlot"))
   bprowsall <- bprowsall()
    gosBP <- go$all[go$all$Ont=="BP",]
    if(is.null(bprowsall)){
      if( dim(gosBP)[1]<10 ){bprowsall <-  seq_len(nrow(gosBP)) }
      else{ bprowsall <-  seq_len(10)  }
    }
    plt <- dotPlotGO(gosBP[bprowsall,], n = length(bprowsall))
    svg$dotbpall <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( bprowsall() )
    plt$x$layout$height <- myHeightfunction(bprowsall() )
    plt
  })
  

  output$dotBpAll <- downloadHandler(
    filename = "dotbpall.svg",
    content = function(file){
      ggsave(file, svg$dotbpall, device = "svg", width = 10, units = "in") }
  )
  # GO gobarplot BP all #######################
  output$gobarplotAllBP <- renderPlot({
    validate(need(go$all, "Load file to render dotPlot"), 
             need(length(bprowsall())>=2, "Select at least 2 row" ) )
    bprowsall <- bprowsall()
    p <- goBarplot(enrichGO = go$all, resGO = res$sh, genes= data$genesall,
              category = "BP", nrows = bprowsall)
    svg$gobarbpall <- p
    print(p)
  })
  
  output$gobarBpAll <- downloadHandler(
    filename = "gobarbpall.svg",
    content = function(file){
      ggsave(file, svg$gobarbpall, device = "svg", width = 10, units = "in") }
  )
  # GO circle BP all #####################
  output$goCircleAllBP <- renderPlot({
    validate(need(go$all, "Load file to render dotPlot"),
             need(res$sh,""),
             need( length(bprowsall())>=4 , "Select at least 4 rows"))
    bprowsall <- bprowsall()
    if(length(bprowsall)>=4){
      go <- go$all[go$all$Ont=="BP",]
      circ <- data2circle(go=go[bprowsall, ], res=res$sh, genes=data$genesall)
      p <- circle(circ, label.size = 3, nsub = length(bprowsall), table.legend = FALSE)
      svg$cirbpall <- p
      print(p)
    }
  })
  output$cirBpAll <- downloadHandler(
    filename = "cirbpall.svg",
    content = function(file){
      ggsave(file, svg$cirbpall, device = "svg", width = 10, units = "in") }
  )
  # GO cloud BP all #######################
  output$cloudBPAll <- renderPlot({
    validate(need(go$all, "Load file to render dotPlot"))
    goall <- go$all[go$all$Ont=="BP" & go$all$level>=input$bpallLevel, ]
    myggwordcloud(goall, bg = "#343e48")
  })
  
  output$cloudbpall <- downloadHandler(
    filename = "cloudbpall.svg",
    content = function(file){
      svg(file, width = 8, height = 6)
      myggwordcloud(go$all[go$all$Ont=="BP" & go$all$level>=input$bpallLevel, ])
      dev.off()
    }
  )
  # ............ ###############################
  # GO table MF all #####################
  output$tableMFall <- DT::renderDataTable(server = FALSE,{
    validate(need(goDT$all, "Load file to render table"))
    goDT <- goDT$all
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-MF all genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="MFall"),list(extend="excel",filename="MFall") ),
           text="Download", title=tituloTabla ) )

    datatable2(goDT[goDT$Ont=="MF",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(order = list(list(6, 'asc')),
                           pageLength = 10, white_space = "normal",
                           buttons = customButtons,
                           scrollY = "400px"
                           #ajax = list(serverSide = TRUE, processing = TRUE)
                           )
    )
  })
  tableMFallProxy <- dataTableProxy("tableMFall")
  observeEvent(input$resettableMFall, {
    tableMFallProxy %>% selectRows(NULL)
  })
  # GO plots MF all  #####################
  output$plotMFall <- renderPlotly({
    validate(need(go$all, "Load file to render plot"),
             need(mfrowsall(), "Select at least one row to plot") )
    mfrowsall <- mfrowsall()
    gosMF <- go$all[go$all$Ont=="MF",]
    if(is.null(mfrowsall)){
        if( dim(gosMF)[1]<10 ){mfrowsall <-  seq_len(nrow(gosMF)) }
        else{ mfrowsall <-  seq_len(10)  }
    }
    p <- plotGOAll(enrichdf = gosMF[mfrowsall, ], nrows = length(mfrowsall), ont="MF", 
                   genesUp = data$genesUp, genesDown = data$genesDown,
                   colors = c(input$downColor, input$upColor))
    if( typeBarMfAll() == "Dodge") {
          plt <- p[[1]]; svg$barmfall <- p[[1]]
          }
        else if ( typeBarMfAll() == "Stack") {
          plt <- p[[2]]; svg$barmfall <- p[[2]]
        }
        else { 
          plt <- p[[3]]; svg$barmfall <- p[[3]]
          }
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( mfrowsall() )
    plt$x$layout$height <- myHeightfunction(mfrowsall() )
    plt
  })
  output$barMfAll <- downloadHandler(
    filename = "barmfall.svg",
    content = function(file){
      ggsave(file, svg$barmfall, "svg", width = 10, units = "in") }
  )
  # GO MF dotplot all ################### 
  output$MFDotall <- renderPlotly({
    validate(need(go$all, "Load file to render dotPlot"),
             need(mfrowsall(), "Select the terms of interest to render DotPlot"))
    mfrowsall <- mfrowsall()
    gosMF <- go$all[go$all$Ont=="MF",]
    if(is.null(mfrowsall)){
      if( dim(gosMF)[1]<10 ){mfrowsall <-  seq_len(nrow(gosMF)) }
      else{ mfrowsall <-  seq_len(10)  }
    }
    plt <- dotPlotGO(gosMF[mfrowsall,], n = length(mfrowsall))
    svg$dotmfall <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( mfrowsall() )
    plt$x$layout$height <- myHeightfunction(mfrowsall() )
    plt
  })
  
  output$dotMfAll <- downloadHandler(
    filename = "dotmfall.svg",
    content = function(file){
      ggsave(file, svg$dotmfall, device = "svg", width = 10, units = "in") }
  )
  # GO gobarplot MF all ####################
  output$gobarplotAllMF <- renderPlot({
    validate(need(go$all, "Load file to render dotPlot"),
             need(length(mfrowsall())>=2, "Select at least 2 row"))
    mfrowsall <- mfrowsall()
    p <- goBarplot(enrichGO = go$all, resGO = res$sh, genes= data$genesall,
              category = "MF", nrows = mfrowsall)
    svg$gobarmfall <- p
    print(p)
  })
  
  output$gobarMfAll <- downloadHandler(
    filename = "gobarmfall.svg",
    content = function(file){
      ggsave(file, svg$gobarmfall, device = "svg", width = 10, units = "in") }
  )
  # GO circle MF all #####################
  output$goCircleAllMF <- renderPlot({
    validate(need(go$all, "Load file to render dotPlot"),
             need(res$sh,""),
             need( length(mfrowsall())>=4 , "Select at least 4 rows"))
    mfrowsall <- mfrowsall()
    if(length(mfrowsall)>=4){
      go <- go$all[go$all$Ont=="MF",]
      circ <- data2circle(go=go[mfrowsall, ], res=res$sh, genes=data$genesall)
      p <- circle(circ, label.size = 3, nsub = length(mfrowsall), table.legend = FALSE)
      svg$cirmfall <- p
      print(p)
    }
  })
  output$cirMfAll <- downloadHandler(
    filename = "cirmfall.svg",
    content = function(file){
      ggsave(file, svg$cirmfall, device = "svg", width = 10, units = "in") }
  )
  # GO cloud MF all #######################
  output$cloudMFAll <- renderPlot({
    validate(need(go$all, "Load file to render dotPlot"))
    goall <- go$all[go$all$Ont=="MF" & go$all$level>=input$mfallLevel, ]
    myggwordcloud(goall, bg = "#343e48")
  })
  
  output$cloudmfall <- downloadHandler(
    filename = "cloudmfall.svg",
    content = function(file){
      svg(file, width = 8, height = 6)
      myggwordcloud(go$all[go$all$Ont=="MF" & go$all$level>=input$mfallLevel, ])
      dev.off()
    }
  )
  # ............ ###############################
  # GO table CC all #####################
  output$tableCCall <- DT::renderDataTable(server=FALSE,{
    validate(need(goDT$all, "Load file to render table"))
    goDT <- goDT$all
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-CC all genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="CCall"),list(extend="excel",filename="CCall") ),
           text="Download", title=tituloTabla ) )
      
    datatable2(goDT[goDT$Ont=="CC",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(order = list(list(6, 'asc')),
                           pageLength = 10, white_space = "normal",
                           buttons = customButtons,
                           scrollY = "400px"
                           #ajax = list(serverSide = TRUE, processing = TRUE)
                           )
    )
  })
  tableCCallProxy <- dataTableProxy("tableCCall")
  observeEvent(input$resettableCCall, {
    tableCCallProxy %>% selectRows(NULL)
  })
  # GO plots CC all #####################
  output$plotCCall <- renderPlotly({
    validate(need(go$all, "Load file to render plot"),
             need(ccrowsall(), "Select at least one row to plot") )
    ccrowsall <- ccrowsall()
    gosCC <- go$all[go$all$Ont=="CC",]
    if(is.null(ccrowsall)){
        if( dim(gosCC)[1]<10 ){ccrowsall <-  seq_len(nrow(gosCC)) }
        else{ ccrowsall <-  seq_len(10)  }
    }
    p <- plotGOAll(enrichdf = gosCC[ccrowsall, ], nrows = length(ccrowsall), ont="CC", 
                   genesUp = data$genesUp, genesDown = data$genesDown,
                   colors = c(input$downColor, input$upColor))
    if( typeBarCcAll() == "Dodge") {
      plt <- p[[1]]; svg$barccall <- p[[1]]
      }
    else if ( typeBarCcAll() == "Stack") {
      plt <- p[[2]] ; svg$barccall <- p[[2]]
      }
    else {
      plt <- p[[3]] ; svg$barccall <- p[[3]]
      }
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( ccrowsall() )
    plt$x$layout$height <- myHeightfunction(ccrowsall() )
    plt
  })
    output$barCcAll <- downloadHandler(
    filename = "barccall.svg",
    content = function(file){
      ggsave(file, svg$barccall, "svg", width = 10, units = "in") }
  )
  # GO CC dotplot all ################### 
  output$CCDotall <- renderPlotly({
    validate(need(go$all, "Load file to render dotPlot"),
             need(ccrowsall(), "Select the terms of interest to render DotPlot"))
    ccrowsall <- ccrowsall()
    gosCC <- go$all[go$all$Ont=="CC",]
    if(is.null(ccrowsall)){
      if( dim(gosCC)[1]<10 ){ccrowsall <-  seq_len(nrow(gosCC)) }
      else{ ccrowsall <-  seq_len(10)  }
    }
    plt <- dotPlotGO(gosCC[ccrowsall,], n = length(ccrowsall))
    svg$dotccall <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( ccrowsall() )
    plt$x$layout$height <- myHeightfunction(ccrowsall() )
    plt
  })
  
  output$dotCcAll <- downloadHandler(
    filename = "dotccall.svg",
    content = function(file){
      ggsave(file, svg$dotccall, device = "svg", width = 10, units = "in") }
  )
  # GO gobarplot CC all #######################
  output$gobarplotAllCC <- renderPlot({
    validate(need(go$all, "Load file to render dotPlot"),
             need(length(ccrowsall())>=2, "Select at least 2 row") )
    ccrowsall <- ccrowsall()
    p <- goBarplot(enrichGO = go$all, resGO = res$sh, genes= data$genesall,
              category = "CC", nrows = ccrowsall)
    svg$gobarmfall <- p
    print(p)
  })
  
  output$gobarMfAll <- downloadHandler(
    filename = "gobarmfall.svg",
    content = function(file){
      ggsave(file, svg$gobarmfall, device = "svg", width = 10, units = "in") }
  )
  # GO circle CC all #####################
  output$goCircleAllCC <- renderPlot({
    validate(need(go$all, "Load file to render dotPlot"),
             need(res$sh,""),
             need( length(ccrowsall())>=4, "Select at least 4 rows"))
    ccrowsall <- ccrowsall()
    if(length(ccrowsall)>=4){
        go <- go$all[go$all$Ont=="CC",]
      circ <- data2circle(go=go[ccrowsall, ], res=res$sh, genes=data$genesall)
      p <- circle(circ, label.size = 3, nsub = length(ccrowsall), table.legend = FALSE)
      svg$circcall <- p
      print(p)
    }
  })
  output$cirCcAll <- downloadHandler(
    filename = "circcall.svg",
    content = function(file){
      ggsave(file, svg$circcall, device = "svg", width = 10, units = "in") }
  )
    # GO cloud CC all #######################
  output$cloudCCAll <- renderPlot({
    validate(need(go$all, "Load file to render dotPlot"))
    goall <- go$all[go$all$Ont=="CC" & go$all$level>=input$ccallLevel, ]
    myggwordcloud(goall, bg = "#343e48")
  })
  
  output$cloudccall <- downloadHandler(
    filename = "cloudccall.svg",
    content = function(file){
      svg(file, width = 8, height = 6)
      myggwordcloud(go$all[go$all$Ont=="CC" & go$all$level>=input$ccallLevel, ])
      dev.off()
    }
  )
  # ............ ###############################
  # GO table BP UP#####################
  output$tableBP <- DT::renderDataTable(server=FALSE,{
    validate(need(goDT$up, "Load file to render table"))
    goDT <- goDT$up
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-BP up-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="BPup"),list(extend="excel",filename="BPup") ),
           text="Download", title=tituloTabla ) )
    
    datatable2(goDT[goDT$Ont=="BP",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(order = list(list(6, 'asc')),
                           pageLength = 10, white_space = "normal",
                           scrollY = "400px",
                           buttons = customButtons))
  })
  tableBPProxy <- dataTableProxy("tableBP")
  observeEvent(input$resettableBP, {
    tableBPProxy %>% selectRows(NULL)
  })
  # GO plots BP UP #####################
  output$plotBP <- renderPlotly({
    validate(need(go$up, "Load file to render plot"),
             need(bprowsup(), "Select at least one row to plot"))
      bprowsup <- bprowsup()
    gosBP <- go$up[go$up$Ont=="BP",]
    if(is.null(bprowsup)){
      if( dim(gosBP)[1]<10 ){bprowsup <-  seq_len(nrow(gosBP)) }
      else{ bprowsup <-  seq_len(10)  }
    }
    p <- plotGO(enrichdf = gosBP[bprowsup, ], nrows = length(bprowsup), ont="BP",
           colors = c(input$upColor) )
    svg$barbpup <- p
    plt <- p
    
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( ccrowsall() )
    plt$x$layout$height <- myHeightfunction(ccrowsall() )
    plt
  })
  
  output$barBpUp <- downloadHandler(
    filename = "barbpup.svg",
    content = function(file){
      ggsave(file, svg$barbpup, "svg", width = 10, units = "in") }
  )
  # GO BP dotplot up ################### 
  output$BPDotUp <- renderPlotly({
    validate(need(go$up, "Load file to render dotPlot"),
             need(bprowsup(), "Select the terms of interest to render DotPlot"))
    bprowsup <- bprowsup()
    gosBP <- go$up[go$up$Ont=="BP",]
    if(is.null(bprowsup)){
      if( dim(gosBP)[1]<10 ){bprowsup <-  seq_len(nrow(gosBP)) }
      else{ bprowsup <-  seq_len(10)  }
    }
    plt <- dotPlotGO(gosBP[bprowsup,], n = length(bprowsup))
    svg$dotbpup <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( bprowsup() )
    plt$x$layout$height <- myHeightfunction(bprowsup() )
    plt
  })
  
  output$dotBpUp <- downloadHandler(
    filename = "dotbpup.svg",
    content = function(file){
      ggsave(file, svg$dotbpup, device = "svg", width = 10, units = "in") }
  )
  
  # GO gobarplot BP Up #######################
  output$gobarplotUpBP <- renderPlot({
    validate(need(go$up, "Load file to render dotPlot"),
             need( length(bprowsup())>=2, "Select at least 2 row") )
    bprowsup <- bprowsup()
    p <- goBarplot(enrichGO = go$up, resGO = res$sh, genes= data$genesUp,
              category = "BP", nrows = bprowsup)
    svg$gobarbpup <- p
    print(p)
  })
  
  output$gobarBpUp <- downloadHandler(
    filename = "gobarbpup.svg",
    content = function(file){
      ggsave(file, svg$gobarbpup, device = "svg", width = 10, units = "in") }
  )
  
    # GO circle BP Up #####################
  output$goCircleUpBP <- renderPlot({
    validate(need(go$up, "Load file to render dotPlot"),
             need(res$sh,""),
             need( length(bprowsup())>=4 , "Select at least 4 rows"))
    bprowsup <- bprowsup()
    if(length(bprowsup)>=4){
        go <- go$up[go$up$Ont=="BP",]
      circ <- data2circle(go=go[bprowsup, ], res=res$sh, genes=data$genesUp)
      p <- circle(circ, label.size = 3, nsub = length(bprowsup), table.legend = FALSE)
      svg$cirbpup <- p
      print(p)
    }
  })
  output$cirBpUp <- downloadHandler(
    filename = "cirbpup.svg",
    content = function(file){
      ggsave(file, svg$cirbpup, device = "svg", width = 10, units = "in") }
  )
  
    # GO cloud BP UP  #######################
  output$cloudBPUp <- renderPlot({
    validate(need(go$up, "Load file to render dotPlot"))
    goup <- go$up[go$up$Ont=="BP" & go$up$level>=input$bpupLevel, ]
    myggwordcloud(goup, bg = "#343e48")
  })
  
  output$cloudbpup <- downloadHandler(
    filename = "cloudbpup.svg",
    content = function(file){
      svg(file, width = 8, height = 6)
      myggwordcloud(go$up[go$up$Ont=="BP" & go$up$level>=input$bpupLevel, ])
      dev.off()
    }
  )
  
  # ............ ###############################
  # GO table MF UP #####################
  output$tableMF <- DT::renderDataTable(server=FALSE, {
    validate(need(goDT$up, "Load file to render table"))
    goDT <- goDT$up
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-MF up-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="MFup"),list(extend="excel",filename="MFup") ),
           text="Download", title=tituloTabla ) )
    
    datatable2(goDT[goDT$Ont=="MF",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(order = list(list(6, 'asc')),
                           pageLength = 10, white_space = "normal",
                           buttons = customButtons,
                           scrollY = "400px"
                           #ajax = list(serverSide = TRUE, processing = TRUE)
                           )
    )
  })
  tableMFProxy <- dataTableProxy("tableMF")
  observeEvent(input$resettableMF, {
    tableMFProxy %>% selectRows(NULL)
  })
  # GO plots MF UP #####################
  output$plotMF <- renderPlotly({
    validate(need(go$up, "Load file to render plot"),
             need(mfrowsup(), "Select at least one row to plot"))
    mfrowsup <- mfrowsup()
    gosMF <- go$up[go$up$Ont=="MF",]
    if(is.null(mfrowsup)){
        if( dim(gosMF)[1]<10 ){mfrowsup <-  seq_len(nrow(gosMF)) }
        else{ mfrowsup <-  seq_len(10)  }
    }
    p <- plotGO(enrichdf = gosMF[mfrowsup, ], nrows = length(mfrowsup), ont = "MF",
           colors = c(input$upColor) )
    svg$barmfup <- p
    plt <- p
    
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( mfrowsup() )
    plt$x$layout$height <- myHeightfunction(mfrowsup() )
    plt
    
  })
  
  output$barMfUp <- downloadHandler(
    filename = "barmfup.svg",
    content = function(file){
      ggsave(file, svg$barmfup, "svg", width = 10, units = "in") }
  )
  # GO MF dotplot up ################### 
  output$MFDotUp <- renderPlotly({
    validate(need(go$up, "Load file to render dotPlot"),
             need(mfrowsup(), "Select the terms of interest to render DotPlot"))
    mfrowsup <- mfrowsup()
    gosMF <- go$up[go$up$Ont=="MF",]
    if(is.null(mfrowsup)){
      if( dim(gosMF)[1]<10 ){mfrowsup <-  seq_len(nrow(gosMF)) }
      else{ mfrowsup <-  seq_len(10)  }
    }
    plt <- dotPlotGO(gosMF[mfrowsup,], n = length(mfrowsup))
    svg$dotmfup <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( mfrowsup() )
    plt$x$layout$height <- myHeightfunction(mfrowsup() )
    plt
  })
  
  output$dotMfUp <- downloadHandler(
    filename = "dotmfup.svg",
    content = function(file){
      ggsave(file, svg$dotmfup, device = "svg", width = 10, units = "in") }
  )
  
  # GO gobarplot MF Up #######################
  output$gobarplotUpMF <- renderPlot({
    validate(need(go$up, "Load file to render dotPlot"),
             need(length(mfrowsup())>=2, "Select at least 2 row") )
      mfrowsup <- mfrowsup()
    p <- goBarplot(enrichGO = go$up, resGO = res$sh, genes= data$genesUp,
              category = "MF", nrows = mfrowsup)
    svg$gobarmfup <- p
    print(p)
  })
  
  output$gobarMfUp <- downloadHandler(
    filename = "gobarmfup.svg",
    content = function(file){
      ggsave(file, svg$gobarmfup, device = "svg", width = 10, units = "in") }
  )
  
  # GO circle MF Up #####################
  output$goCircleUpMF <- renderPlot({
    validate(need(go$up, "Load file to render dotPlot"),
             need(res$sh,""),
             need( length(mfrowsup())>=4 , "Select at least 4 rows"))
    mfrowsup <- mfrowsup()
    if(length(mfrowsup)>=4){
        go <- go$up[go$up$Ont=="MF",]
      circ <- data2circle(go=go[mfrowsup, ], res=res$sh, genes=data$genesUp)
      p <- circle(circ, label.size = 3, nsub = length(mfrowsup), table.legend = FALSE)
      svg$cirmfup <- p
      print(p)
    }
  })
  output$cirMfUp <- downloadHandler(
    filename = "cirmfup.svg",
    content = function(file){
      ggsave(file, svg$cirmfup, device = "svg", width = 10, units = "in") }
  )
  
  # GO cloud MF UP  #######################
  output$cloudMFUp <- renderPlot({
    validate(need(go$up, "Load file to render dotPlot"))
    goup <- go$up[go$up$Ont=="MF" & go$up$level>=input$mfupLevel, ]
    myggwordcloud(goup, bg = "#343e48")
  })
  
  output$cloudmfup <- downloadHandler(
    filename = "cloudmfup.svg",
    content = function(file){
      svg(file, width = 8, height = 6)
      myggwordcloud(go$up[go$up$Ont=="MF" & go$up$level>=input$mfupLevel, ])
      dev.off()
    }
  )
  # ............ ###############################
  # GO table CC UP #####################
  output$tableCC <- DT::renderDataTable(server=FALSE,{
    validate(need(goDT$up, "Load file to render table"))
    goDT <- goDT$up
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-CC up-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="CCup"),list(extend="excel",filename="CCup") ),
           text="Download",  title=tituloTabla ) )
    
    datatable2(goDT[goDT$Ont=="CC",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(order = list(list(6, 'asc')),
                           pageLength = 10, white_space = "normal",
                           buttons = customButtons,
                           scrollY = "400px"
                           #ajax = list(serverSide = TRUE, processing = TRUE)
                           )
    )
  })
  tableCCProxy <- dataTableProxy("tableCC")
  observeEvent(input$resettableCC, {
    tableCCProxy %>% selectRows(NULL)
  })
  # GO plots CC UP #####################
  output$plotCC <- renderPlotly({
    validate(need(go$up, "Load file to render plot"),
             need(ccrowsup(), "Select at least one row to plot"))
    ccrowsup <- ccrowsup()
    gosCC <- go$up[go$up$Ont=="CC",]
    if(is.null(ccrowsup)){
        if( dim(gosCC)[1]<10 ){ccrowsup <-  seq_len(nrow(gosCC)) }
        else{ ccrowsup <-  seq_len(10)  }
    }
    p <- plotGO(enrichdf = gosCC[ccrowsup,], nrows = length(ccrowsup), ont="CC",
           colors = c(input$upColor))
    svg$barccup <- p
    plt <- p
    
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( ccrowsup() )
    plt$x$layout$height <- myHeightfunction(ccrowsup() )
    plt
  })
  
  output$barCcUp <- downloadHandler(
    filename = "barccup.svg",
    content = function(file){
      ggsave(file, svg$barccup, "svg", width = 10, units = "in") }
  )
  
  # GO CC dotplot up ################### 
  output$CCDotUp <- renderPlotly({
    validate(need(go$up, "Load file to render dotPlot"),
             need(ccrowsup(), "Select the terms of interest to render DotPlot"))
    ccrowsup <- ccrowsup()
    gosCC <- go$up[go$up$Ont=="CC",]
    if(is.null(ccrowsup)){
      if( dim(gosCC)[1]<10 ){ccrowsup <-  seq_len(nrow(gosCC)) }
      else{ ccrowsup <-  seq_len(10)  }
    }
    plt <- dotPlotGO(gosCC[ccrowsup,], n = length(ccrowsup))
    svg$dotccup <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( ccrowsup() )
    plt$x$layout$height <- myHeightfunction(ccrowsup() )
    plt
  })
  
  output$dotCcUp <- downloadHandler(
    filename = "dotccup.svg",
    content = function(file){
      ggsave(file, svg$dotccup, device = "svg", width = 10, units = "in") }
  )
  
  # GO gobarplot CC Up #######################
  output$gobarplotUpCC <- renderPlot({
    validate(need(go$up, "Load file to render dotPlot"),
             need(length(ccrowsup())>=2, "Select at least 2 row"))
    ccrowsup <- ccrowsup()
    p <- goBarplot(enrichGO = go$up, resGO = res$sh, genes= data$genesUp,
              category = "CC", nrows = ccrowsup)
    svg$gobarccup <- p
    print(p)
  })
  
  output$gobarCcUp <- downloadHandler(
    filename = "gobarccup.svg",
    content = function(file){
      ggsave(file, svg$gobarccup, device = "svg", width = 10, units = "in") }
  )
  # GO circle CC Up #####################
  output$goCircleUpCC <- renderPlot({
    validate(need(go$up, "Load file to render dotPlot"),
             need(res$sh,""),
             need( length(ccrowsup())>=4 , "Select at least 4 rows"))
      ccrowsup <- ccrowsup()
    if(length(ccrowsup)>=4){
        go <- go$up[go$up$Ont=="CC",]
      circ <- data2circle(go=go[ccrowsup, ], res=res$sh, genes=data$genesUp)
      p <- circle(circ, label.size = 3, nsub = length(ccrowsup), table.legend = FALSE)
      svg$circcup <- p
      print(p)
    }
  })
  
  output$cirCcUp <- downloadHandler(
    filename = "circcup.svg",
    content = function(file){
      ggsave(file, svg$circcup, device = "svg", width = 10, units = "in") }
  )
  
    # GO cloud CC UP  #######################
  output$cloudCCUp <- renderPlot({
    validate(need(go$up, "Load file to render dotPlot"))
    goup <- go$up[go$up$Ont=="CC" & go$up$level>=input$ccupLevel, ]
    myggwordcloud(goup, bg = "#343e48")
  })
  
  output$cloudccup <- downloadHandler(
    filename = "cloudccup.svg",
    content = function(file){
      svg(file, width = 8, height = 6)
      myggwordcloud(go$up[go$up$Ont=="CC" & go$up$level>=input$ccupLevel, ])
      dev.off()
    }
  )
  # ............ ###############################
  # GO table BP DOWN #####################
  output$tableBPdown <- DT::renderDataTable(server=FALSE,{
    validate(need(goDT$down, "Load file to render table"))
    goDT <- goDT$down
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-BP down-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="BPdown"),list(extend="excel",filename="BPdown") ),
           text="Download", title=tituloTabla ) )
    
    datatable2(goDT[goDT$Ont=="BP",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(order = list(list(6, 'asc')),
                           buttons = customButtons,
                           scrollY = "400px",
                           pageLength = 10, white_space = "normal"
                           #ajax = list(serverSide = TRUE, processing = TRUE)
                           )
    )
  })
  tableBPdownProxy <- dataTableProxy("tableBPdown")
  observeEvent(input$resettableBPdown, {
    tableBPdownProxy %>% selectRows(NULL)
  })
  # GO plots BP DOWN #####################
  output$plotBPdown <- renderPlotly({
    validate(need(go$down, "Load file to render plot"),
             need(bprowsdown(), "Select at least one row to plot"))
    bprowsdown <- bprowsdown()
    gosBP <- go$down[go$down$Ont=="BP",]
    if(is.null(bprowsdown)){
        if( dim(gosBP)[1]<10 ){bprowsdown <-  seq_len(nrow(gosBP)) }
        else{ bprowsdown <-  seq_len(10)  }
    }
    p <- plotGO(enrichdf = gosBP[bprowsdown, ], nrows = length(bprowsdown), ont="BP",
           colors = c(input$downColor))
    svg$barbpdown <- p
    plt <- p
    
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( bprowsdown() )
    plt$x$layout$height <- myHeightfunction(bprowsdown() )
    plt
  })
  
  output$barBpDown <- downloadHandler(
    filename = "barbpdown.svg",
    content = function(file){
      ggsave(file, svg$barbpdown, "svg", width = 10, units = "in") }
  )
  
  # GO BP dotplot down ################### 
  output$BPDotDown <- renderPlotly({
    validate(need(go$down, "Load file to render dotPlot"),
             need(bprowsdown(), "Select the terms of interest to render DotPlot"))
    bprowsdown <- bprowsdown()
    gosBP <- go$down[go$down$Ont=="BP",]
    if(is.null(bprowsdown)){
      if( dim(gosBP)[1]<10 ){bprowsdown <-  seq_len(nrow(gosBP)) }
      else{ bprowsdown <-  seq_len(10)  }
    }
    plt <- dotPlotGO(gosBP[bprowsdown,], n = length(bprowsdown))
    svg$dotbpdown <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( bprowsdown() )
    plt$x$layout$height <- myHeightfunction(bprowsdown() )
    plt
  })
  
  
  output$dotBpDown <- downloadHandler(
    filename = "dotbpdown.svg",
    content = function(file){
      ggsave(file, svg$dotbpdown, device = "svg", width = 10, units = "in") }
  )
  # GO gobarplot BP down #######################
  output$gobarplotDownBP <- renderPlot({
    validate(need(go$down, "Load file to render dotPlot"),
             need(length(bprowsdown())>=2, "Select at least 2 row"))
    bprowsdown <- bprowsdown()
    p <- goBarplot(enrichGO = go$down, resGO = res$sh, genes= data$genesDown,
              category = "BP", nrows = bprowsdown)
    svg$gobarbpdown <- p
    print(p)
  })
  
  output$gobarBpDown <- downloadHandler(
    filename = "gobarbpdown.svg",
    content = function(file){
      ggsave(file, svg$gobarbpdown, device = "svg", width = 10, units = "in") }
  )
  # GO circle BP Down #####################
  output$goCircleDownBP <- renderPlot({
    validate(need(go$down, "Load file to render dotPlot"),
             need(res$sh,""),
             need( length(bprowsdown())>=4 , "Select at least 4 rows"))
    bprowsdown <- bprowsdown()
    if(length(bprowsdown)>=4){
        go <- go$down[go$down$Ont=="BP",]
      circ <- data2circle(go=go[bprowsdown, ], res=res$sh, genes=data$genesDown)
      p <- circle(circ, label.size = 3, nsub = length(bprowsdown), table.legend = FALSE)
      svg$cirbpdown <- p
      print(p)
    }
  })
  output$cirBpDown <- downloadHandler(
    filename = "cirbpdown.svg",
    content = function(file){
      ggsave(file, svg$cirbpdown, device = "svg", width = 10, units = "in") }
  )
  # GO cloud BP Down #######################
  output$cloudBPDown <- renderPlot({
    validate(need(go$down, "Load file to render dotPlot"))
    godown <- go$down[go$down$Ont=="BP" & go$down$level>=input$bpdownLevel, ]
    myggwordcloud(godown, bg = "#343e48")
  })
  
  output$cloudbpdown <- downloadHandler(
    filename = "cloudbpdown.svg",
    content = function(file){
      svg(file, width = 8, height = 6)
      myggwordcloud(go$down[go$down$Ont=="BP" & go$down$level>=input$bpdownLevel, ])
      dev.off()
    }
  )
  # ............ ###############################
  # GO table MF DOWN #####################
  output$tableMFdown <- DT::renderDataTable(server=FALSE, {
    validate(need(goDT$down, "Load file to render table"))
    goDT <- goDT$down
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-MF down-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="MFdown"),list(extend="excel",filename="MFdown") ),
           text="Download", title=tituloTabla ) )
    
    datatable2(goDT[goDT$Ont=="MF",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(order = list(list(6, 'asc')),
                           pageLength = 10, white_space = "normal",
                           buttons = customButtons,
                           scrollY = "400px"
                           #ajax = list(serverSide = TRUE, processing = TRUE)
                           )
    )
  })
  tableMFdownProxy <- dataTableProxy("tableMFdown")
  observeEvent(input$resettableMFdown, {
    tableMFdownProxy %>% selectRows(NULL)
  })
  # GO plots MF DOWN #####################
  output$plotMFdown <- renderPlotly({
    validate(need(go$down, "Load file to render plot"),
             need(mfrowsdown(), "Select at least one row to plot"))
    mfrowsdown <- mfrowsdown()
    gosMF <- go$down[go$down$Ont=="MF",]
    if(is.null(mfrowsdown)){
        if( dim(gosMF)[1]<10 ){mfrowsdown <-  seq_len(nrow(gosMF)) }
        else{ mfrowsdown <-  seq_len(10)  }
    }
    p <- plotGO(enrichdf = gosMF[mfrowsdown, ], nrows = length(mfrowsdown), ont = "MF",
           colors = c(input$downColor) )
    svg$barmfdown <- p
    plt <- p
    
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( mfrowsdown() )
    plt$x$layout$height <- myHeightfunction(mfrowsdown() )
    plt
  })
  
  output$barMfDown <- downloadHandler(
    filename = "barmfdown.svg",
    content = function(file){
      ggsave(file, svg$barmfdown, "svg", width = 10, units = "in") }
  )
  # GO MF dotplot down ################### 
  output$MFDotDown <- renderPlotly({
    validate(need(go$down, "Load file to render dotPlot"),
             need(mfrowsdown(), "Select the terms of interest to render DotPlot"))
    mfrowsdown <- mfrowsdown()
    gosMF <- go$down[go$down$Ont=="MF",]
    if(is.null(mfrowsdown)){
      if( dim(gosMF)[1]<10 ){mfrowsdown <-  seq_len(nrow(gosMF)) }
      else{ mfrowsdown <-  seq_len(10)  }
    }
    plt <- dotPlotGO(gosMF[mfrowsdown,], n = length(mfrowsdown))
    svg$dotmfdown <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( mfrowsdown() )
    plt$x$layout$height <- myHeightfunction(mfrowsdown() )
    plt
  })
  
  
  output$dotMfDown <- downloadHandler(
    filename = "dotmfdown.svg",
    content = function(file){
      ggsave(file, svg$dotmfdown, device = "svg", width = 10, units = "in") }
  )
  # GO gobarplot MF down #######################
  output$gobarplotDownMF <- renderPlot({
    validate(need(go$down, "Load file to render dotPlot"),
             need(length(mfrowsdown())>=2, "Select at least 2 row"))
    mfrowsdown <- mfrowsdown()
    p <- goBarplot(enrichGO = go$down, resGO = res$sh, genes= data$genesDown,
              category = "MF", nrows = mfrowsdown)
    svg$gobarmfdown <- p
    print(p)
  })
  
  output$gobarMfDown <- downloadHandler(
    filename = "gobarmfdown.svg",
    content = function(file){
      ggsave(file, svg$gobarmfdown, device = "svg", width = 10, units = "in") }
  )
  # GO circle MF Down #####################
  output$goCircleDownMF <- renderPlot({
    validate(need(go$down, "Load file to render dotPlot"),
             need(res$sh,""),
             need( length(mfrowsdown())>=4 , "Select at least 4 rows"))
    mfrowsdown <- mfrowsdown()
    if(length(mfrowsdown)>=4){
        go <- go$down[go$down$Ont=="MF",]
      circ <- data2circle(go=go[mfrowsdown, ], res=res$sh, genes=data$genesDown)
      p <- circle(circ, label.size = 3, nsub = length(mfrowsdown), table.legend = FALSE)
      svg$cirmfdown <- p
      print(p)
    }
  })
  output$cirMfDown <- downloadHandler(
    filename = "cirmfdown.svg",
    content = function(file){
      ggsave(file, svg$cirmfdown, device = "svg", width = 10, units = "in") }
  )
    # GO cloud MF Down #######################
  output$cloudMFDown <- renderPlot({
    validate(need(go$down, "Load file to render dotPlot"))
    godown <- go$down[go$down$Ont=="MF" & go$down$level>=input$mfdownLevel, ]
    myggwordcloud(godown, bg = "#343e48")
  })
  
  output$cloudmfdown <- downloadHandler(
    filename = "cloudmfdown.svg",
    content = function(file){
      svg(file, width = 8, height = 6)
      myggwordcloud(go$down[go$down$Ont=="MF" & go$down$level>=input$mfdownLevel, ])
      dev.off()
    }
  )
  # ............ ###############################
  # GO table CC DOWN #####################
  output$tableCCdown <- DT::renderDataTable(server=FALSE,{
    validate(need(goDT$down, "Load file to render table"))
    goDT <- goDT$down
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-CC down-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="CCdown"),list(extend="excel",filename="CCdown") ),
           text="Download", title=tituloTabla ) )
    
    datatable2(goDT[goDT$Ont=="CC",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(order = list(list(6, 'asc')),
                           pageLength = 10, white_space = "normal",
                           buttons = customButtons,
                           scrollY = "400px"
                           #ajax = list(serverSide = TRUE, processing = TRUE)
                           )
    )
  })
  tableCCdownProxy <- dataTableProxy("tableCCdown")
  observeEvent(input$resettableCCdown, {
    tableCCdownProxy %>% selectRows(NULL)
  })
  # GO plots CC DOWN #####################
  output$plotCCdown <- renderPlotly({
    validate(need(go$down, "Load file to render plot"),
    need(ccrowsdown(), "Select at least one row to plot"))
    ccrowsdown <- ccrowsdown()
    gosCC <- go$down[go$down$Ont=="CC",]
    if(is.null(ccrowsdown)){
      if( dim(gosCC)[1]<10 ){ccrowsdown <-  seq_len(nrow(gosCC)) }
      else{ ccrowsdown <-  seq_len(10)  }
    }
    p <- plotGO(enrichdf = gosCC[ccrowsdown,], nrows = length(ccrowsdown), ont="CC",
           colors = c(input$downColor) )
    svg$barccdown <- p
    plt <- p
    
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( ccrowsdown() )
    plt$x$layout$height <- myHeightfunction(ccrowsdown() )
    plt
  })
  
  output$barCcDown <- downloadHandler(
    filename = "barccdown.svg",
    content = function(file){
      ggsave(file, svg$barccdown, "svg", width = 10, units = "in") }
  )
  # GO CC dotplot down ################### 
  output$CCDotDown <- renderPlotly({
    validate(need(go$down, "Load file to render dotPlot"),
             need(ccrowsdown(), "Select the terms of interest to render DotPlot"))
    ccrowsdown <- ccrowsdown()
    gosCC <- go$down[go$down$Ont=="CC",]
    if(is.null(ccrowsdown)){
        if( dim(gosCC)[1]<10 ){ccrowsdown <-  seq_len(nrow(gosCC)) }
        else{ ccrowsdown <-  seq_len(10)  }
    }
    plt <- dotPlotGO(gosCC[ccrowsdown,], n = length(ccrowsdown))
    svg$dotccdown <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( ccrowsdown() )
    plt$x$layout$height <- myHeightfunction(ccrowsdown() )
    plt
  })
  
  output$dotCcDown <- downloadHandler(
    filename = "dotccdown.svg",
    content = function(file){
      ggsave(file, svg$dotccdown, device = "svg", width = 10, units = "in") }
  )
  # GO gobarplot CC down #######################
  output$gobarplotDownCC <- renderPlot({
    validate(need(go$down, "Load file to render dotPlot"),
             need(length(ccrowsdown())>=2, "Select at least 2 row"))
    ccrowsdown <- ccrowsdown()
    p <- goBarplot(enrichGO = go$down, resGO = res$sh, genes= data$genesDown,
              category = "CC", nrows = ccrowsdown)
    svg$gobarccdown <- p
    print(p)
  })
  
  output$gobarCcDown <- downloadHandler(
    filename = "gobarccdown.svg",
    content = function(file){
      ggsave(file, svg$gobarccdown, device = "svg", width = 10, units = "in") }
  )
  # GO circle CC Down #####################
  output$goCircleDownCC <- renderPlot({
    validate(need(go$down, "Load file to render dotPlot"),
             need(res$sh,""),
             need( length(ccrowsdown() )>=4 , "Select at least 4 rows"))
    ccrowsdown <- ccrowsdown()
    if(length(ccrowsdown)>=4){
        go <- go$down[go$down$Ont=="CC",]
      circ <- data2circle(go=go[ccrowsdown, ], res=res$sh, genes=data$genesDown)
      p <- circle(circ, label.size = 3, nsub = length(ccrowsdown), table.legend = FALSE)
      svg$circcdown <- p
      print(p)
    }
  })
  output$cirCcDown <- downloadHandler(
    filename = "circcdown.svg",
    content = function(file){
      ggsave(file, svg$circcdown, device = "svg", width = 10, units = "in") }
  )
    # GO cloud CC Down #######################
  output$cloudCCDown <- renderPlot({
    validate(need(go$down, "Load file to render dotPlot"))
    godown <- go$down[go$down$Ont=="CC" & go$down$level>=input$ccdownLevel, ]
    myggwordcloud(godown, bg = "#343e48")
  })
  
  output$cloudccdown <- downloadHandler(
    filename = "cloudccdown.svg",
    content = function(file){
      svg(file, width = 8, height = 6)
      myggwordcloud(go$down[go$down$Ont=="CC" & go$down$level>=input$ccdownLevel, ])
      dev.off()
    }
  )
  # GSEA......... ###############################
  output$gseaSelectize <- renderUI({
    if(specie()=="Mm"){
      datasets <- list.files("./resources/Mm/GSEA/")
    }else if(specie()=="Hs"){
      datasets <- list.files("./resources/Hs/GSEA/")
    }else{datasets <- list.files("./resources/Rn/GSEA")}
    pickerInput(inputId = "gseadb", label = "Select GSEA dataset",
                   choices = datasets, 
                options = list(title = "dataset"),
                selected = NULL )
  })
  # GSEA table ##########################
  output$gseaTable <- DT::renderDataTable(server=FALSE, {
    validate(need(res$sh, "Load file to render table"))
    validate(need(input$gseadb!="","Select dataset"))
    gsea$gsea <- gseaKegg(res$sh, specie(), gseadb = input$gseadb )
    mygsea <- gsea$gsea
    if( length(which(mygsea@result$p.adjust<=0.05)) == 0 ){
        createAlert(session, anchorId = "gsea", title = "Oops!!", 
          content = "Sorry, I didn't get any significant results for this analysis",
          append=FALSE, style = "info")
    } else{
    table <- mygsea@result[mygsea@result$p.adjust<=0.05 ,2:9] %>% 
      mutate_at(vars(3:7), ~round(., 4))

    tituloTabla <- paste0("Table: GSEA pathway | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="GSEAkegg"),list(extend="excel",filename="GSEAkegg") ),
           text="Download", title=tituloTabla ) )
    
    DT::datatable( table,
                   rownames=FALSE,
                   filter = list(position="top", clear=FALSE),
                   options = list(order = list(list(4, 'asc')),
                     lengthMenu = list(c(10,25,50,100,-1), c(10,25,50,100,"All")),
                     columnDefs = list(list(orderable = FALSE,
                                            className = "details-control",
                                            targets = 1)
                     ),
                     dom = "Bfrtipl",
                     buttons = customButtons,
                     scrollY = "400px",
                     list(pageLength = 10, white_space = "normal")
                   )
    )
    }
  })
  
  gseaTableProxy <- dataTableProxy("gseaTable")
  observeEvent(input$resetgseaTable, {
    gseaTableProxy %>% selectRows(NULL)
  })
  # GSEA plot ##########################
  output$gseaPlot <- renderPlot({
    validate(need(gsea$gsea, "Load file to render table"))
    gseanr <- gsearow()
    if(is.null(gseanr)){gseanr <- c(1)}
    mygsea <- gsea$gsea
    if( length(which(mygsea@result$p.adjust<=0.05)) == 0 ){
        createAlert(session, anchorId = "gseaPlot", title = "Oops!!", 
          content = "Sorry, I didn't get any significant results for this analysis",
          append=FALSE, style = "info")
    } else{
        p <- enrichplot::gseaplot2(gsea$gsea, geneSetID = gseanr, pvalue_table = TRUE, ES_geom = "line")
        svg$gseaplot <- p
        print(p)
        }
  })
  output$gseaButton <- downloadHandler(
    filename = "gseaplot.svg",
    content = function(file){
      ggsave(file, svg$gseaplot, device = "svg", width = 10, units = "in") }
  )
  # ............ ###############################
  # generate report #############################

  observeEvent(input$report2, {
     showModal( modalDialog(
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
                                    selected = c("PCA", "BoxPlot", "Heatmap", "Cluster","Top6",
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
                                    choices = c("Table", "Barplot", "Chorplot", "Dotplot",
                                                "Heatmap", "Netplot"),
                                    selected = c("Table", "Barplot", "Chorplot", "Dotplot",
                                                 "Heatmap", "Netplot"),
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
                                    choices = c("Table", "Barplot", "Chorplot", "Dotplot",
                                                "Heatmap", "Netplot"),
                                    selected = c("Table", "Barplot", "Chorplot", "Dotplot",
                                                 "Heatmap", "Netplot"),
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
                                    choices = c("Table", "Barplot", "Chorplot", "Dotplot",
                                                "Heatmap", "Netplot"),
                                    selected = c("Table", "Barplot", "Chorplot", "Dotplot",
                                                 "Heatmap", "Netplot"),
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
                                    choices = c("Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
                                    selected = c("Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
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
                                    choices = c("Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
                                    selected = c("Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
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
                                    choices = c("Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
                                    selected = c("Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
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
                                    selected = c("Table", "GSEA plot"),
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
         actionButton("unselect","Select/Unselect all"),
         modalButton("Cancel"),
         actionButton("ok", "Apply"),
         uiOutput("downloadhtml")
       )
     ) )
   })

  
  observeEvent(input$unselect, {
    if (input$unselect > 0) {
      if (input$unselect %% 2 == 0) {
        selectPopUpModal(session = session)
      } else{
        unselectPopUpModal(session = session)
      }
    }
  })
  
  
  applyPress <- reactiveValues(ok = FALSE)
  observeEvent(input$ok, {
    applyPress$ok <- TRUE
    vals$preview <- input$modalPreview
    vals$keggAll <- input$modalkeggAll
    vals$keggUp <- input$modalkeggUp
    vals$keggDown <- input$modalkeggDown
    vals$GOAll <- input$modalGOAll
    vals$GOUp <- input$modalGOUp
    vals$GODown <- input$modalGODown
    vals$GSEA <- input$modalGSEA
  })
  
  output$downloadhtml <- renderUI({
    validate(need(isTRUE(applyPress$ok), ""))
    downloadButton("download", "Download report")
  })
  
    output$download <- downloadHandler(
    filename = "report.html",
    content = function(file) {
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      file.copy("report.css", file.path(tempdir(), "report.css"), overwrite = TRUE)
      file.copy("utilsReport.R", file.path(tempdir(),"utils.R"), overwrite = TRUE)
      file.copy("resources/", tempdir(), overwrite = TRUE, recursive = TRUE)
      file.copy("resources/dna-svg-small-13.gif",
      file.path(tempdir(), "resources/dna-svg-small-13.gif"), overwrite = TRUE)
      ## inicializar variables preview
      pcaObj <- boxObj <- heatObj <- clusterObj <- top6Obj <- top1Obj <- FALSE
      karyObj <- FALSE
      volcObj <- maObj <- FALSE
      ## inicializar variables kegg
      tablekgaObj <- barkgaObj <- chorkgaObj <- dotkgaObj <- heatkgaObj <- netkgaObj <- FALSE
      tablekguObj <- barkguObj <- chorkguObj <- dotkguObj <- heatkguObj <- netkguObj <- FALSE
      tablekgdObj <- barkgdObj <- chorkgdObj <- dotkgdObj <- heatkgdObj <- netkgdObj <- FALSE
      ## inicializar variables Go
      tablegoaObj <- bargoaObj <- dotgoaObj <- gobargoaObj <- gocirclegoaObj <- cloudgoaObj <- FALSE
      tablegouObj <- bargouObj <- dotgouObj <- gobargouObj <- gocirclegouObj <- cloudgouObj <- FALSE
      tablegodObj <- bargodObj <- dotgodObj <- gobargodObj <- gocirclegodObj <- cloudgodObj <- FALSE
      ## inicializar variables GSEA
      tablegseaObj <- plotgseaObj <- FALSE 
      ## Asigna variables
          rlogdatos <- rlog$datos; colorespca <- coloresPCA$colores();
          variables <- variables(); samplename <- samplename() 
          vsddata <- rlog$datos; boxplotswitch <- boxplotswitch()
          specie <- specie(); numheatmap <- numheatmap()
          ressh <- res$sh; datosdds <- datos$dds; gene <- gene(); 
          padj <- padj(); logfc <- logfc(); genesvolcano <- genesVolcano();
          upcolor <- input$upColor; downcolor <- input$downColor
          kggall <- kgg$all; genesdeup <- numgenesDE$up; genesdedown <- numgenesDE$down
          kggdtall <- kggDT$all; datagenesup <- data$genesUp; datagenesdown <- data$genesDown
          typebarkeggall <- typeBarKeggAll()
          typebarbpall <- typeBarBpAll(); typebarmfall <- typeBarMfAll();
          typebarccall <- typeBarCcAll()
          kggup <- kgg$up; kggdown <- kgg$down; kggdtup <- kggDT$up; kggdtdown <- kggDT$down; 
          goall <- go$all; godtall <- goDT$all; 
          goup <- go$up; godtup <- goDT$up; 
          godown <- go$down; godtdown <- goDT$down; 
          bprowsall <- bprowsall(); mfrowsall <- mfrowsall(); ccrowsall <- ccrowsall()
          bprowsup <- bprowsup(); mfrowsup <- mfrowsup(); ccrowsup <- ccrowsup()
          bprowsdown <- bprowsdown(); mfrowsdown <- mfrowsdown(); ccrowsdown <- ccrowsdown()
          gsearow <- gsearow(); gseagsea <- gsea$gsea
          textnotes <- input$textNotes
      #nrows
          nrowsall <- rowsAll()
          if(!is.null(kggDT$all)){
            if(is.null(nrowsall)){ 
              nrowsall <-  ( if( dim(kggDT$all)[1]<10) seq_len(nrow(kggDT$all)) else seq_len(10) ) }
          }
          nrowsup <- rowsUp()
          if(!is.null(kggDT$up)){
            if(is.null(nrowsup)){ 
              nrowsup <-  ( if( dim(kggDT$up)[1]<10) seq_len(nrow(kggDT$up)) else seq_len(10) ) }
          }
          nrowsdown <- rowsdown()
          if(!is.null(kggDT$down)){
            if(is.null(nrowsdown)){ 
              nrowsdown <-  ( if( dim(kggDT$down)[1]<10) seq_len(nrow(kggDT$down)) else seq_len(10) ) }
          }
      if(!is.null(vals$preview)){      #para preview
        if( ("PCA" %in% vals$preview) ){ pcaObj <- TRUE}
        if("BoxPlot" %in% vals$preview){boxObj <- TRUE}
        if("Heatmap" %in% vals$preview){heatObj <- TRUE}
        if("Cluster" %in% vals$preview){clusterObj <- TRUE}
        if("Top6" %in% vals$preview){top6Obj <- TRUE}
        if("Top1" %in% vals$preview){top1Obj <- TRUE}
        if("Karyoplot" %in% vals$preview){karyObj <- TRUE}
        if("Volcano" %in% vals$preview){volcObj <- TRUE}
        if("MA" %in% vals$preview){ maObj <- TRUE}
      }
      if(!is.null(vals$keggAll)){ #para keggAll
        if("Table" %in% vals$keggAll){ tablekgaObj <- TRUE }
        if("Barplot" %in% vals$keggAll){ barkgaObj <- TRUE }
        if("Chorplot" %in% vals$keggAll){ chorkgaObj <- TRUE }
        if("Dotplot" %in% vals$keggAll){ dotkgaObj <- TRUE }
        if("Heatmap" %in% vals$keggAll){ heatkgaObj <- TRUE }
        if("Netplot" %in% vals$keggAll){ netkgaObj <- TRUE }
      }
      if(!is.null(vals$keggUp)){ #para keggUp
        if("Table" %in% vals$keggUp){ tablekguObj <- TRUE }
        if("Barplot" %in% vals$keggUp){ barkguObj <- TRUE }
        if("Chorplot" %in% vals$keggUp){ chorkguObj <- TRUE }
        if("Dotplot" %in% vals$keggUp){ dotkguObj <- TRUE }
        if("Heatmap" %in% vals$keggUp){ heatkguObj <- TRUE }
        if("Netplot" %in% vals$keggUp){ netkguObj <- TRUE }
      }
      if(!is.null(vals$keggDown)){ #para keggDown
        if("Table" %in% vals$keggDown){ tablekgdObj <- TRUE }
        if("Barplot" %in% vals$keggDown){ barkgdObj <- TRUE }
        if("Chorplot" %in% vals$keggDown){ chorkgdObj <- TRUE }
        if("Dotplot" %in% vals$keggDown){ dotkgdObj <- TRUE }
        if("Heatmap" %in% vals$keggDown){ heatkgdObj <- TRUE }
        if("Netplot" %in% vals$keggDown){ netkgdObj <- TRUE }
      }
      if(!is.null(vals$GOAll)){#para GoAll
        if("Table" %in% vals$GOAll){ tablegoaObj <- TRUE }
        if("Barplot" %in% vals$GOAll){ bargoaObj <- TRUE }
        if("Dotplot" %in% vals$GOAll){ dotgoaObj <- TRUE }
        if("GObarplot" %in% vals$GOAll){ gobargoaObj <- TRUE }
        if("WordCloud" %in% vals$GOAll){ cloudgoaObj <- TRUE}
        if("GOcircleplot" %in% vals$GOAll){ gocirclegoaObj <- TRUE }
      }
      if(!is.null(vals$GOUp)){#para GoUp
        if("Table" %in% vals$GOUp){ tablegouObj <- TRUE }
        if("Barplot" %in% vals$GOUp){ bargouObj <- TRUE }
        if("Dotplot" %in% vals$GOUp){ dotgouObj <- TRUE }
        if("GObarplot" %in% vals$GOUp){ gobargouObj <- TRUE }
        if("WordCloud" %in% vals$GOUp){ cloudgouObj <- TRUE}
        if("GOcircleplot" %in% vals$GOUp){ gocirclegouObj <- TRUE }
      }
      if(!is.null(vals$GODown)){#para GoDown
        if("Table" %in% vals$GODown){ tablegodObj <- TRUE }
        if("Barplot" %in% vals$GODown){ bargodObj <- TRUE }
        if("Dotplot" %in% vals$GODown){ dotgodObj <- TRUE }
        if("GObarplot" %in% vals$GODown){ gobargodObj <- TRUE }
        if("WordCloud" %in% vals$GODown){ cloudgodObj <- TRUE}
        if("GOcircleplot" %in% vals$GODown){ gocirclegodObj <- TRUE }
      }
      if(!is.null(vals$GSEA)){#para GSEA
        if("Table" %in% vals$GSEA){ tablegseaObj <- TRUE}
        if("GSEA plot" %in% vals$GSEA){ plotgseaObj <- TRUE}
        }

      params <- list( values = vals, 
                      pcaObj = pcaObj, rlog = rlogdatos, colorespca = colorespca,
                     variables = variables, samplename = samplename,
                     vsd = vsddata, boxplotswitch = boxplotswitch, boxObj = boxObj,
                     specie = specie, numheatmap = numheatmap, heatObj = heatObj,
                     clusterObj = clusterObj, 
                     ressh = ressh, datosdds = datosdds, top6Obj = top6Obj, 
                     gene = gene, top1Obj = top1Obj,
                     karyObj = karyObj, padj =padj, logfc = logfc,
                     volcObj = volcObj, genesvolcano = genesvolcano, 
                     upcolor = upcolor, downcolor = downcolor, 
                     maObj = maObj, 
                     tablekgaObj = tablekgaObj, kggall = kggall, genesdedown = genesdedown,
                     genesdeup = genesdeup, kggdtall = kggdtall,
                     barkgaObj = barkgaObj, nrowsall = nrowsall, datagenesdown = datagenesdown, 
                     datagenesup = datagenesup, typebarkeggall = typebarkeggall,
                     chorkgaObj = chorkgaObj, dotkgaObj = dotkgaObj, heatkgaObj = heatkgaObj,
                     netkgaObj = netkgaObj, tablekguObj = tablekguObj, barkguObj = barkguObj,
                     chorkguObj = chorkguObj, dotkguObj =dotkguObj, heatkguObj = heatkguObj,
                     netkguObj = netkguObj, tablekgdObj = tablekgdObj, barkgdObj = barkgdObj,
                     chorkgdObj = chorkgdObj, dotkgdObj = dotkgdObj, heatkgdObj = heatkgdObj,
                     netkgdObj = netkgdObj, kggup = kggup, kggdown = kggdown, kggdtup = kggdtup, 
                     kggdtdown = kggdtdown, nrowsup = nrowsup, nrowsdown = nrowsdown, 
                     typebarbpall=typebarbpall, typebarmfall=typebarmfall, typebarccall=typebarccall,
                     tablegoaObj = tablegoaObj, bargoaObj=bargoaObj, dotgoaObj=dotgoaObj,
                     gobargoaObj=gobargoaObj,gocirclegoaObj=gocirclegoaObj, tablegouObj = tablegouObj,
                     bargouObj=bargouObj, dotgouObj=dotgouObj, gobargouObj=gobargouObj,
                     gocirclegouObj=gocirclegouObj, tablegodObj = tablegodObj, bargodObj=bargodObj,
                     dotgodObj=dotgodObj, gobargodObj=gobargodObj, gocirclegodObj=gocirclegodObj,
                     cloudgoaObj=cloudgoaObj, cloudgouObj=cloudgouObj, cloudgodObj=cloudgodObj,
                     goall = goall, godtall=godtall, goup = goup, godtup=godtup,
                     godown = godown, godtdown=godtdown,
                     bprowsall=bprowsall, mfrowsall=mfrowsall, ccrowsall=ccrowsall,
                     bprowsup=bprowsup, mfrowsup=mfrowsup, ccrowsup=ccrowsup,
                     bprowsdown=bprowsdown, mfrowsdown=mfrowsdown, ccrowsdown=ccrowsdown,
                     gsearow = gsearow, gseagsea = gseagsea, tablegseaObj = tablegseaObj,
                     plotgseaObj = plotgseaObj, textnotes = textnotes)
      
      params <- c(params, list(tempdir=tempdir() ))
      removeModal()
      applyPress$ok <- FALSE
      rmarkdown::render(
        tempReport,
        output_file = file,
        params = params,
        envir = new.env(parent = globalenv( ))
      )
    } )
}


shinyApp(ui, server)
