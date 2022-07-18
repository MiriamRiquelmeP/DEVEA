library(shinydashboard)
library(AnnotationDbi)
library(org.Mm.eg.db) #Mus musculus
library(org.Hs.eg.db) #Homo sapiens
#library(org.Dr.eg.db) #Danio rerio (zebra fish)
library(org.Rn.eg.db) #Rattus norvegicus
#library(org.Mmu.eg.db) #Macaca mulata
library(chorddiag)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Rnorvegicus.v79)
library(limma)
library(tidyverse)
library(DT)
library(RColorBrewer)
library(purrr)
library(plotly)
library(ggpubr)
library(DESeq2)
library(fgsea)
library(shinyalert)
library(shinyBS)
library(shinyWidgets)
library(shinydashboardPlus)
library(pheatmap)
library(heatmaply)
library(shinyjs)
library(shinythemes)
library(shinymanager)
library(rgl)
library(rglwidget)
library(scales)
library(stringr)
library(shinybusy)
library(visNetwork)
library(ggrepel)
library(circlize)
library(mychordplot)
#library(ggwordcloud)
#library(wordcloud2)
library(randomcoloR)
library(tidytext)
source("global.R")
source("UpdatepopModals.R")
source("utils.R")
options(shiny.maxRequestSize = 3000*1024^2)


### HEADER ############ 
header <- dashboardHeader(title = "Go Simple DEVEA", 
                          titleWidth = 300, 
                          dropdownMenuOutput("messageMenu"),
                          tags$li(class="dropdown", 
                                  actionButton("notesButton","Report notes"),
                                  style="margin-top:8px; margin-right: 5px"),
                          tags$li(class = "dropdown",
                                  #actionButton("report", "HTML report"),
                                  uiOutput("report"),
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
                            sidebarMenu(
                              menuItem(
                              "Input data",
                              tabName = "info",
                              icon = icon("info")
                            )),
                            #sidebarMenu(
                            #    menuItem(
                                #     pickerInput(
                                #         inputId = "specie",
                                #         label = "Select species",
                                #         choices = list( "Human" = "Hs", "Mouse" = "Mm", "Rat" = "Rn"),
                                #         options = list(title = "species"),
                                #         selected = "Mm"
                                #     ) 
                                # ),
                                # menuItem(
                                #   pickerInput(
                                #     inputId = "annotation",
                                #     label = "Select annotation gene",
                                #     choices = list("Ensembl" = "ensg", "Symbol"="symbol"),
                                #     options = list(title="annotation"),
                                #     selected = "ensg"
                                #   )
                                # ),
                                sidebarMenu("", sidebarMenuOutput("prevw")),
                                sidebarMenu("", sidebarMenuOutput("menuKegg")),
                                sidebarMenu("", sidebarMenuOutput("menuGO")),
                                sidebarMenu("", sidebarMenuOutput("menuGSEA")),
                            #),
                            tags$div(
                            #   box(width = 12,
                            #     h5(strong("Generate report"), align = 'center'),
                            #     sidebarMenu( 
                            #       menuItem(
                            #         fluidRow(column(12, align = "center", offset=0,
                            #                         uiOutput("report")))))
                            # ),
                            
                            tags$a(href='https://jacob.cea.fr/drf/ifrancoisjacob/Pages/Departements/MIRCen/themes/astrocytes-reactifs-biomarqueurs-imagerie-cibles-therapeutiques.aspx', target="_blank",
                                   tags$img(src='mircen.png',width='49%',
                                            style="padding: 5px; position: absolute; bottom:0px; left:0") ),
                            tags$a(href='http://www.bioinformatica.imib.es', target="_blank",
                                   tags$img(src='IMIB_color_gris.svg',width='51%',
                                            style="padding: 5px; float: right; bottom:5px;") ),
                            tags$a(href='https://jacob.cea.fr/drf/ifrancoisjacob/Pages/Departements/CNRGH/LABORATOIRES/Bio-analyse.aspx', target="_blank",
                                   tags$img(src='cnrgh.png',width='49%',
                                            style="padding: 5px; float: left; bottom:15px;") ),
                            
                            style = "position: absolute; bottom:0;width:100%;"
                            ) #fin div
                            ) # fin sideDashBoard

### BODY ###############
body <- dashboardBody(
      tags$script(HTML("$('body').addClass('fixed');")),
      add_busy_gif(src="dna-mini.gif", position = "full-page", width = 10, height = 10 ),
     #htmltools::includeCSS("./www/customDark.css"),
      #tags$head(
        #HTML("<link rel = 'stylesheet' type = 'text/css' href = 'customDark.css'>"),
      #),
    setShadow(class = "shiny-plot-output"),
    setShadow( class = "box"),
    setShadow( class = "svg-container"),
    tags$head(
      HTML("<link rel = 'stylesheet' type = 'text/css' href = 'customDark.css'>"),
      tags$style(HTML(".irs-min, .irs-max {
              color: rgb(215,215,215) !important;
              background-color: rgb(45,55,65) !important;
          }"))
      ),
  bsAlert("alert"),
  tabItems(
    # Initial INFO
    tabItem(tabName = "info",
            fluidRow(column(width = 4,
              box(width=12, 
                    pickerInput(
                          inputId = "specie",
                          label = "Select species",
                          choices = list( "Human" = "Hs", "Mouse" = "Mm", "Rat" = "Rn"),
                          options = list(title = "species"),
                          selected = "Mm"
                      ),
                    pickerInput(
                      inputId = "annotation",
                      label = "Select gene annotation",
                      choices = list("Ensembl" = "ensg", "Symbol"="symbol"),
                      options = list(title="annotation"),
                      selected = "ensg"
                      ),
                  uiOutput(outputId = "geneList"),
                  column(
                    width = 1,
                    circleButton(
                      inputId = "informationGL",
                      icon = icon("info"),
                      size = "xs",
                      status = "primary"
                    ),
                    bsTooltip(
                      "informationGL",
                      paste0("The accepted formats are .txt, .tsv, .xlsx"),
                      
                      trigger = "hover",
                      placement = "right"
                    )),
                  #menuItem(uiOutput("circleinfoGL")),
                  #menuItem(uiOutput("tooltipGL")),
                  uiOutput(outputId = "geneFile"),
                  uiOutput(outputId = "geneButton")
                  ),
              box(width=12,
                  dataTableOutput(outputId = "table1col"))
              ),
            column(width = 8,
            box(width=12,
                status = "info",
                title = h1(strong("Welcome to Go Simple DEVEA app 2022!") ),
                h3("Data requirements:"),
                p(HTML(paste0("Option 1. ", '<B>Gene list and statistical values (GL + SV) mode: </B>')), HTML(paste0("consists of a ", '<B>gene list (GL)</B>')),
                  HTML(paste0("associated with ", '<B>statistical values (SV).</B>')),  
                  "The first column should contain gene names, the second column the fold-change and the 
                  third column the statistical adjusted p-value, in this precise order. 
                  The values will be used without further modifications by the tool. 
                  This table can be uploaded as a file or directly integrated in the 
                  application by copy/paste in the dedicated field, following the same order as the file mode.", 
                  br(),
                  br(),
                  HTML(paste0("Option 2. ", '<B>Gene list (GL) mode: </B>')), HTML(paste0("based on a unique ", '<B>gene list (GL)</B>')), 
                  " containing the favourite gene names for the analysis. The gene list can be 
                  uploaded directly to the application by copy/paste in the dedicated field. "),
                br(),
                p("Note that demo data is also available here to test the performance.
                For further details on how to use DEVEA,
                please take advantage of every symbol of info" , icon("info-circle"), "that you might find,
                see our ",
                  a("paper", href = "https://f1000research.com/"), 
                  "and a detailed walkthrough in the 'Tutorial' section. Source code at ",
                  a("GitHub.", href = "https://github.com/MiriamRiquelmeP/Full-EnrichApp"),
                  "For feedbacks, please contact us (information available in the 'About' section).") ),
            
            menuItem(uiOutput("circleinfoBack")),
            menuItem(uiOutput("tooltipBack")),
            uiOutput("universefile"),
            uiOutput("enrichbutton")
            )
            )
    ),
    # preview tab
    tabItem(tabName = "preview",
            source(file = "ui-preview-tab.R",
            local=TRUE,
            encoding = "UTF-8"
            )$value),
    # kegg tab content
    tabItem(tabName = "kegg",
            source(file = "ui-kegg-tab.R",
                   local = TRUE,
                   encoding = "UTF-8"
                   )$value),
    # GO tab GO tab
    tabItem( tabName = "go",
             source(file = "ui-go-tab.R",
                    local = TRUE,
                    encoding = "UTF-8",
                    )$value),
    # GSEA tab
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

ui <- dashboardPage(title="Go Simple DEVEA",
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
# 
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
# )

########################################## SERVER #################################################
server <- function(input, output, session) {
  
   #    res_auth <- secure_server(
   #   check_credentials = check_credentials(
   #       "/datos/repos/darkEnrichApp/users.sqlite",
   #       passphrase = readRDS("/datos/repos/enrichapp_listable/dbpass.Rds")
   #   )
   # )
  #res_auth <- secure_server(
  #   check_credentials = check_credentials(
  #       "/datos/repos/darkEnrichApp/users.sqlite",
  #       passphrase = readRDS("/datos/repos/enrichapp_listable/dbpass.Rds")
  #   )
  # )
    
    
    enrichflag=NULL
    
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
  #       #tags$iframe(src="https://155.54.120.105/shiny/enrich_listable/pres1.html",  width="850px", height="700px")
  #       tags$iframe(src="pres2.html",  width="850px", height="700px")
  #     )
  #   )
  # })
  
  shinyjs::onclick("moreinfo", runjs("window.open('tutorial.html','_blank')") )
  
  
  # variables reactivas ######
  annotation <- reactive({input$annotation})
  data <- reactiveValues(df=NULL, dfilt=NULL)
  df3cols <- reactiveValues(TF=FALSE)
  enrichflag <- reactiveValues(one=NULL, three=NULL)
  fc_switch <- reactive({input$fc_switch})
  fcRange <- reactiveValues() # min y max fc
  go <- reactiveValues(all=NULL)
  goDT <- reactiveValues(all=NULL)
  gene <- reactiveValues(lost=NULL)
  genes <- reactiveValues()
  genesVolcano <- reactive({input$genesVolcano})
  gsea <- reactiveValues(gsea=NULL)
  kgg <- reactiveValues(all=NULL)
  kggDT <- reactiveValues(all=NULL)
  logfcRange <- reactiveValues() # min y max logfc
  numgenesDE <- reactiveValues(up=NULL, down=NULL)
  specie <- reactive({input$specie})
  typeBarKeggAll <- reactive({input$selectkeggall})
  validatedGene <- reactiveValues(list=NULL)
  vals <- reactiveValues()
  svg <- reactiveValues()
  
  gsearow <- reactive({input$gseaTable_rows_selected}) 
  
  fileuniverse <- reactive({input$fileuniverse})
  funiverse <- reactive({input$universef})
  
  ## Leer data ##########################
  observeEvent(input$geneButton,{
    # comprobaciones listado manual
    if( !is.null(input$geneList)){
      if(input$geneList!="" ){
        if (annotation() == "ensg") {
          if (length(which(grepl("^ENS", input$geneList, ignore.case = TRUE))) <
              length(input$geneList)  ) {
            shinyalert("Oops!!", "One or more genes are not 
            in ENSEMBL format",
                       type = "error")
          } else{
            validatedGene$list <- validateGeneList(input$geneList, specie(), annotation() )
            data$df <- formatData( validatedGene$list, specie(), annotation() )
            # lost <- which(is.na(data$df$ENTREZID))
            # gene$lost <- data$df$ENSEMBL[lost]
            # if(length(lost)!=0){ data$df <- data$df[-lost, ] }
            data$df$SYMBOL <- ifelse(is.na(data$df$SYMBOL), data$df$ENSEMBL, data$df$SYMBOL )
            }
        }
        if (annotation() == "symbol") {
          if (length(which(grepl("^ENS", input$geneList, ignore.case = TRUE))) ==  length(input$geneList) ) {
            shinyalert("Oops!!", "Looks like this entry is 
                       ENSEMBL please check your selection", 
                       type = "error")
          } else{
            validatedGene$list <- validateGeneList(input$geneList, specie(), annotation() )
            data$df <- formatData( validatedGene$list, specie(), annotation() )
            # lost <- which(is.na(data$df$ENTREZID))
            # gene$lost <- data$df$SYMBOL[lost]
            # if(length(lost)!=0){ data$df <- data$df[-lost, ] }
          }
          }
        }
      }
    ## comprobaciones cargar fichero
    if(!is.null(input$geneFile$datapath)){
      if(input$geneFile$datapath != "" ){
          if(input$geneFile$type == "text/plain"){
                validatedGene$list <- as.data.frame(read.table(input$geneFile$datapath,
                                                 header = F, sep = "\t"))
          }else
          if(input$geneFile$type == "text/csv"){
              validatedGene$list <- as.data.frame( read.table(input$geneFile$datapath,
                                                              header = F, sep = ";") )
          } else{
            validatedGene$list <- as.data.frame( readxl::read_xlsx(input$geneFile$datapath, sheet = 1)) } 
        if (annotation() == "ensg") {
          if(length(which(grepl("^ENS", validatedGene$list[,1], ignore.case = TRUE))) < nrow(validatedGene$list) ) {
            shinyalert("Oops!!", "One or more genes are not 
            in ENSEMBL format",
                       type = "error" )
          } else{
                data$df <- formatData( validatedGene$list, specie(), annotation() )
                # lost <- which(is.na(data$df$ENTREZID))
                # gene$lost <- data$df$ENSEMBL[lost]
                # if(length(lost)!=0){ data$df <- data$df[-lost, ] }
                data$df$SYMBOL <- ifelse(is.na(data$df$SYMBOL), data$df$ENSEMBL, data$df$SYMBOL )
          }
        }
        if (annotation() == "symbol") {
          if (length(which(grepl("^ENS", validatedGene$list[,1], ignore.case = TRUE))) ==  nrow(validatedGene$list)) {
            shinyalert("Oops!!", "Looks like this entry is 
                       ENSEMBL please check your selection", 
                       type = "error" )
          } else{
                data$df <- formatData( validatedGene$list, specie(), annotation() )
                # lost <- which(is.na(data$df$ENTREZID))
                # gene$lost <- data$df$SYMBOL[lost]
                # if(length(lost)!=0){ data$df <- data$df[-lost, ] }
          }
        }
      }
    }
    if( dim(data$df)[2]==5 ){
      if( length( which(data$df$logFC > 0)) == 0 | length( which(data$df$logFC < 0)) ==0){
        shinyalert("Warning", "It seems that the data do not have both up and down regulated genes.
                   It will be considered as a simple gene list.", 
                       type = "warning")
        df3cols$TF <- FALSE
        data$df <- data$df %>% select(-logFC, -pval)
      }else{
        df3cols$TF <- TRUE
        logfcRange$min <- min(data$df$logFC)
        logfcRange$max <- max(data$df$logFC)
        fcRange$min <- ifelse(logfcRange$min<0, -(2^abs(logfcRange$min)), 2^abs(logfcRange$min))
        fcRange$max <- ifelse(logfcRange$max<0, -(2^abs(logfcRange$max)), 2^abs(logfcRange$max))
      }}
    
    updateTabItems(session, "previewMenu", "preview")
  })
  
  ## Pulsar Enrich Button para listado simple ##################################
  observeEvent(input$enrichButtons,{
    if( dim(data$df)[2]==3 ){
        if( is.null( fileuniverse() ) ){
            bckgnd <- NULL
        }else{ 
            universe <- read.table(fileuniverse()$datapath, header = F)
            bckgnd <- geneIdConverter2( universe[,1], specie = specie() )
        }
      lost <- which(is.na(data$df$ENTREZID))
      gene$lost <- data$df$SYMBOL[lost]
      if(length(lost)!=0){ data$dfilt <- data$df[-lost, ] }else{data$dfilt <- data$df}
      kgg$all <- customKegg(data$dfilt[,c("SYMBOL","ENTREZID") ], species = specie(), universe = bckgnd$ENTREZID  )
      kggDT$all <- kegg2DT(kgg$all, data$dfilt[,c("SYMBOL","ENTREZID") ] )
      go$all <- customGO(data$dfilt[,c("SYMBOL","ENTREZID") ], species = specie(), universe = bckgnd$ENTREZID )
      goDT$all <- go2DT(enrichdf = go$all, data = data$dfilt[,c("SYMBOL","ENTREZID") ] )
      enrichflag$one <- TRUE
      updateTabItems(session, "KeggTab", "kegg")
      hideTab(inputId = "keggTabSetPanel", target = "keggDownTab")
      hideTab(inputId = "keggTabSetPanel", target = "keggUpTab")
      hideTab(inputId = "goTabSetPanel", target = "goUpTab")
      hideTab(inputId = "goTabSetPanel", target = "goDownTab")
      hideTab(inputId = "boxPanelBP", target = "gobarplotallbp")
      hideTab(inputId = "boxPanelMF", target = "gobarplotallmf")
      hideTab(inputId = "boxPanelCC", target = "gobarplotallcc")
      hideTab(inputId = "boxPanelBP", target = "gocirplotallbp")
      hideTab(inputId = "boxPanelMF", target = "gocirplotallmf")
      hideTab(inputId = "boxPanelCC", target = "gocirplotallcc")
    }  
    })
  ## Pulsar Enrich Button para listado 3 columnas ##################################
  observeEvent(input$enrichButton,{
    if( dim(data$df)[2]==5 ){
        if( is.null( funiverse() ) ){
            bckgnd <- NULL
        }else{ 
            universe <- read.table(funiverse()$datapath, header = F)
            bckgnd <- geneIdConverter2( universe[,1], specie = specie() )
        }
      lost <- which(is.na(data$df$ENTREZID))
      gene$lost <- data$df$SYMBOL[lost]
      if(length(lost)!=0){ data$dfilt <- data$df[-lost, ] }else{data$dfilt <- data$df}
      genes$Up <- data$dfilt[data$dfilt$logFC >= logfc()[2] & data$dfilt$pval <= padj(),
                          c("SYMBOL","ENTREZID")]
      genes$Down <- data$dfilt[data$dfilt$logFC <= logfc()[1] & data$dfilt$pval <= padj(),
                          c("SYMBOL","ENTREZID")]
      genes$all <- rbind(genes$Up, genes$Down)
      kgg$all <- customKegg(genes$all, species = specie(), universe = bckgnd$ENTREZID ) 
      kggDT$all <- kegg2DT(kgg$all, genes$all)
      kgg$up <- customKegg(genes$Up, species = specie(), universe = bckgnd$ENTREZID ) 
      kggDT$up <- kegg2DT(kgg$up, genes$Up)
      kgg$down <- customKegg(genes$Down, species = specie(), universe = bckgnd$ENTREZID ) 
      kggDT$down <- kegg2DT(kgg$down, genes$Down)
      go$all <- customGO(genes$all, species = specie(), universe = bckgnd$ENTREZID )
      goDT$all <- go2DT(enrichdf = go$all, data = genes$all )
      go$up <- customGO(genes$Up, species = specie(), universe = bckgnd$ENTREZID )
      goDT$up <- go2DT(enrichdf = go$up, data = genes$Up )
      go$down <- customGO(genes$Down, species = specie(), universe = bckgnd$ENTREZID )
      goDT$down <- go2DT(enrichdf = go$down, data = genes$Down )
      enrichflag$three <- TRUE
      updateTabItems(session, "KeggTab", "kegg")
    }
  })
## ........................ #####################
##datatable preview 1 columna #####################
  output$table1col <- DT::renderDataTable(server = FALSE,{
  shiny::validate(need(data$df, ""))
    shiny::validate(need( dim(data$df)[2]==3,"")) 
  customButtons <- list(
        list(extend = "copy", title="Preview table"),
        list(extend="collection", 
             buttons = list( list(extend="csv",filename="coldata"),list(extend="excel",filename="coldata") ),
             text="Download", filename="coldata", title="Preview table" ) )
    
    datatable( data$df, extensions = "Buttons", caption ="Preview table for input gene list",
               rownames=FALSE,
               filter = list(position="top", clear=FALSE),
               options = list(
                 dom = "Bfrtipl",
                 lengthMenu = list(c(10,25,50,100,-1), c(10,25,50,100,"All")),
                  columnDefs = list(list(orderable = TRUE,
                                        className = "details-control",
                                        targets = 1)),
                 buttons = customButtons,
                 scrollY = "400px",
                 list(pageLength = 10, white_space = "normal")
               )
    )
  })
## Cosas a renderizar en preview si dflist 3 columns ##################
  ## sidebar menu preview ###################
  output$prevw <- renderMenu({
      shiny::validate(need(isTRUE(df3cols$TF), ""))
      sidebarMenu(id = "previewMenu",
          menuItem(
              "Preview",
              tabName = "preview",
              icon = icon("chart-bar")
          )
          )
          })
  
  ## Gene List ####################
  output$geneList <- renderUI({
    shiny::validate(need(specie(),""))
    shiny::validate(need(annotation(),""))
    textAreaInput(inputId = "geneList", label = "Input simple gene list ...", resize = "vertical")
  })
  
  ## Gene File #####################
  
  output$infoCM <- renderUI({
    shiny::validate(need(specie(),""))
    shiny::validate(need(annotation(),""))
  })
  
  output$circleinfoGL <- renderUI({
    shiny::validate(need(specie(),""))
    shiny::validate(need(annotation(),""))
    circleButton(inputId = "infoGL",
                 icon = icon("info"),
                 size = "xs",
                 status = "primary"
    )
  })
  output$tooltipGL <- renderUI({
    shiny::validate(need(specie(),""))
    shiny::validate(need(annotation(),""))
    bsTooltip(
      "infoGL",
      paste0("The accepted formats are .txt, .tsv, .xlsx"),
      trigger = "hover",
      placement = "right"
    )
  })
  
  
  output$geneFile <- renderUI({
    shiny::validate(need(specie(),""))
    shiny::validate(need(annotation(),""))
    fileInput(inputId = "geneFile", label="...or upload file with gene list")
  })
  
  ## Boton validar ###############
  output$geneButton <- renderUI({
    shiny::validate(need(specie(),""))
    shiny::validate(need(annotation(),""))
    actionButton("geneButton", label = "Click to validate data")
  })
  ## upload universe file #################
  
  output$circleinfoBack <- renderUI({
    shiny::validate(need(data$df, ""))
    shiny::validate(need(!isTRUE(df3cols$TF), ""))
    circleButton(
      inputId = "informaBack",
      icon = icon("info"),
      size = "xs",
      status = "primary"
    ) 
  })
  output$tooltipBack <- renderUI({
    shiny::validate(need(data$df, ""))
    shiny::validate(need(!isTRUE(df3cols$TF), ""))
    bsTooltip(
      "informaBack",
      paste0("By default, genes used as enrichment ",
             "background consist of the whole stated species genome. ", 
             "Add your own universe list as a file with your favourite genes. ",
             "Click Run enrichment to unlock the enrichment analysis tabs."),
      trigger = "hover",
      placement = "right"
    )
  })
      
  output$universefile <- renderUI({
      shiny::validate(need(data$df, ""))
      shiny::validate(need(!isTRUE(df3cols$TF), ""))
      fileInput("fileuniverse", "Background dataset", placeholder = "Leave empty to use entire database")
  })
  ## boton enrich #########################
  output$enrichbutton <- renderUI({
    shiny::validate(need(data$df, ""))
    shiny::validate(need(!isTRUE(df3cols$TF), ""))
    actionBttn("enrichButtons", label = "Click to run enrichment", size="lg", color="default", icon = icon("images"))
    })
  

## sidebar menu kegg ###################
  output$menuKegg <- renderMenu({
      shiny::validate(need(kgg$all, ""))
      sidebarMenu(id = 'KeggTab',
          menuItem(
              "Kegg Enrichment",
              tabName = "kegg",
              icon = icon("chart-bar")
          )
          )
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
  
      
  ## side menubar GO #########################
      output$menuGO <- renderMenu({
      shiny::validate(need(go$all,""))
      sidebarMenu(  
        menuItem(
              "GO Enrichment",
              tabName = "go",
              icon = icon("chart-bar")
          )
        )
      })
  ## sidebar menu GSEA #####################333
  output$menuGSEA <- renderMenu({
    shiny::validate(need(go$all,""))
      sidebarMenu(  
          menuItem("GSEA",
                   tabName = "gsea",
                   icon = icon("chart-line"))
        )
      })
    # ui selector de genes para volcano plot #######################
  output$geneSelector <- renderUI({
    shiny::validate(need(data$df, ""))
    genes <- as.character(data$df$SYMBOL[ which(!( data$df$pval>padj() &
                                                             data$df$logFC>logfc()[1] &
                                                             data$df$logFC<logfc()[2] )) ])
    selectInput("genesVolcano", label="Select gene[s] to label",
                choices = genes,
                multiple = TRUE)
  })
  
  # Deslizador fc/logfc según switch #################
  output$fc_control <- renderUI({
    if(isTRUE(fc_switch())){
      shiny::validate(need(data$df, ""))
      valmin <- ifelse(input$logfc[1]<0, -2^(abs(input$logfc[1] )), 2^(abs(input$logfc[1])) )
      valmax <- ifelse(input$logfc[2]<0, -2^(abs(input$logfc[2] )), 2^(abs(input$logfc[2])) )
      sliderInput("fc", label = "Select FC range to remove (keeps the tails)",
                  min=round(fcRange$min,3), max=round(fcRange$max, 3),
                  value = c(valmin, valmax), step = 0.1 )
    } else {
      shiny::validate(need(data$df, ""))
      shiny::validate(need(fc(), ""))
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
    shiny::validate(need(data$df,""))
    sliderInput("padj", label = "Select p-adjusted threshold", min = 0, max=0.2,
                value=0.05, step = 0.005 )
  })
  
    # infoboxes ###############################
  output$allbox <- renderInfoBox({
      shiny::validate(need(data$df, ""))
      shiny::validate(need(padj(), ""))
      shiny::validate(need(logfc(), ""))
      numall <- nrow( data$df[ ((data$df$logFC >= logfc()[2] |
                                    data$df$logFC< logfc()[1]) &
                                   data$df$pval <= padj() ),] ) 
      infoBox("All DE genes", numall, icon = icon("arrows-alt-v"), color = "light-blue", fill = TRUE)
  })
  output$upbox <- renderInfoBox({
      shiny::validate(need(data$df, ""))
      shiny::validate(need(padj(), ""))
      shiny::validate(need(logfc(), ""))
      numup <- nrow( data$df[(data$df$logFC >= logfc()[2]) & (data$df$pval <= padj()), ]) 
      numgenesDE$up <- numup
      infoBox("Upregulated genes", numup, icon = icon("thumbs-up", lib = "glyphicon"), color = "light-blue", fill=TRUE)
  })
  output$downbox <- renderInfoBox({
      shiny::validate(need(data$df, ""))
      shiny::validate(need(padj(), ""))
      shiny::validate(need(logfc(), ""))
      numdown <- nrow( data$df[(data$df$logFC <= logfc()[1]) & (data$df$pval <= padj()), ])
      numgenesDE$down <- numdown
      infoBox("Downregulated genes", numdown, icon = icon("thumbs-down", lib = "glyphicon"), color = "light-blue", fill=TRUE)
  })
  
  output$fcdown <- renderUI({
        shiny::validate(need(logfcRange$min, ""))
        shiny::validate(need(logfc(),""))
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
          #angleArc = 180,
          #angleOffset = angleOffset,
          width = "80%",
          height = "80%"
      )
  })
  output$fcup <- renderUI({
        shiny::validate(need(logfcRange$min, ""))
        shiny::validate(need(logfc(),""))
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
          value = round(logfc()[2], 2),
          readOnly = TRUE,
          min = min,
          max=max,
          rotation=rotation,
          displayPrevious = TRUE,
          fgColor = fgColor,
          inputColor = inputColor,
          bgColor = bgColor,
          #angleArc = 180,
          #angleOffset = angleOffset,
          width = "80%",
          height = "80%"
      )
  })
  output$pval <- renderUI({
      shiny::validate(need(padj(), ""))
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
          #angleArc = 180,
          #angleOffset = 90,
          width = "80%",
          height = "80%"
      )
  })
#....................... ####
## table preview ######################################
output$tablepreview <- DT::renderDataTable(server=FALSE,{
  shiny::validate(need(data$df, ""))
  customButtons <- list(
        list(extend = "copy", title="Preview table"),
        list(extend="collection", 
             buttons = list( list(extend="csv",filename="preview"),list(extend="excel",filename="preview") ),
             text="Download", filename="preview", title="Preview table" ) )
    
    datatable( data$df, extensions = "Buttons", caption ="Preview table for input gene list",
               rownames=FALSE,
               filter = list(position="top", clear=FALSE),
               options = list(
                 dom = "Bfrtipl",
                 lengthMenu = list(c(10,25,50,100,-1), c(10,25,50,100,"All")),
                  columnDefs = list(list(orderable = TRUE,
                                        className = "details-control",
                                        targets = 1)),
                 buttons = customButtons,
                 scrollY = "400px",
                 list(pageLength = 10, white_space = "normal")
               )
    )
})

output$lostgene <- renderText({
  shiny::validate(need(data$df, ""))
  lost <- length(which(is.na(data$df$ENTREZID)))
  print(paste0(lost," out of ", dim(data$df)[1]," genes have no ENTREZ Id. Theses genes will be missing in enrichment analysis" ) )
})
# ....................... ####
  # volcano plot #########
  output$volcano <- renderPlot( {
    shiny::validate(need(data$df, "Load file to render plot"))
    res <-  data$df
    svg$volcano <- CustomVolcano(res, lab = as.character(res$SYMBOL),
                  selectLab = genesVolcano(),
                    x = 'logFC',
                    y = 'pval',
                    pCutoff = padj(),
                    FCcutoffUP = logfc()[2],
                    FCcutoffDOWN = logfc()[1],
                  drawconnectors = TRUE,
                    #xlim = c(-8, 8),
                    col = c("gray", "#7cccc3", "#d99c01", input$upColor, input$downColor))
    svg$volcano
    })

output$downVolcano <- downloadHandler(
  filename = "volcano.svg",
  content = function(file){
    ggsave(file, svg$volcano, "svg",width = 10, units = "in")}
)

xy <- reactive({
  res <- data$df
  res$`-log10padj` <- (-log10(res$pval)) 
  nearPoints(res, input$plot_click1, xvar = "logFC", yvar = "-log10padj")
})
output$texto1 <- renderTable( digitgeoms = -2, {
        xy <- xy()
        xy[,c(2,4,5,6)]
    })
#....................... ####
## karyoplot ######################################
output$karyoPlot <- renderPlot({
    shiny::validate(need(data$df, "Load file to render plot"))
    krtp(data$df, specie = specie(), pval = padj(), fcdown = logfc()[1],
         fcup = logfc()[2], bg="#46505a", coldown="#4ADBFF" , colup="#f7665c", annotation=annotation() )
})

output$downKrpt <- downloadHandler(
  filename = "karyoplot.png",
  content = function(file){
    png(file)
    krtp(data$df, specie = specie(), pval = padj(), fcdown = logfc()[1],
         fcup = logfc()[2], bg="#46505a", coldown="#4ADBFF" , colup="#f7665c",
         annotation=annotation() )
    dev.off()
    }
)
# .......................####
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
  # variables KEGG ALL ##########################
  rowsAll <- reactive({input$tableAll_rows_selected})
 # KEGG table all #####################################
  output$tableAll <- DT::renderDataTable(server=FALSE,{
    shiny::validate(need(kgg$all, "Load file to render table"))
    names(kggDT$all)[names(kggDT$all) == "DE"] <- "DEG"
    names(kggDT$all)[names(kggDT$all) == "P.DE"] <- "p-value"
    tituloTabla <- paste0("Table: Kegg all genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="keggAll"),list(extend="excel",filename="keggAll") ),
           text="Download", filename="keggall", title=tituloTabla ) )

    datatable2(
      kggDT$all,
      vars = c("genes"),
      filter = list(position="top", clear=FALSE),
      escape = FALSE,
      opts = list(order = list(list(5, 'asc')),
        pageLength = 10, white_space = "normal",
        scrollY = "400px",
        buttons = customButtons))
  }) 

tableallProxy <- dataTableProxy("tableAll")
observeEvent(input$resettableall, {
  tableallProxy %>% selectRows(NULL)
})
  # KEGG barplot all################
  output$keggPlotAll <- renderPlotly ({
    shiny::validate(need(kgg$all, "Load file to render BarPlot"), 
                    need(rowsAll(), "Select the paths of interest to render BarPlot"))
    rowsAll <- rowsAll()
    if(is.null(rowsAll)){
        if( dim(kgg$all)[1]<10 ){rowsAll <-  seq_len(nrow(kgg$all)) }
        else{ rowsAll <-  seq_len(10)  }
    }
    if (isTRUE(df3cols$TF)) {
      p <- plotKeggAll( enrichdf = kgg$all[rowsAll,], nrows = length(rowsAll),
                        genesUp = genes$Up, genesDown = genes$Down, 
                        colors = c(input$downColor, input$upColor) )
      if (typeBarKeggAll() == "Dodge") {
        plt <- p[[1]]; svg$keggall <- p[[1]]
      }
      else if (typeBarKeggAll() == "Stack") {
        plt <- p[[2]]; svg$keggall <- p[[2]]
      }
      else {
        plt <- p[[3]]; svg$keggall <- p[[3]]
      }
    } else{
      # caso de que sea sólo una lista simple
      p <- plotKegg( enrichdf = kgg$all[rowsAll,], nrows = length(rowsAll), 
                     colors = "#045a8d" )
      plt <- p[[1]]
      svg$keggall <- p[[2]]
    }
    
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
    plt$height <- myHeightfunction( rowsAll() )
    plt$x$layout$height <- myHeightfunction(rowsAll() )
    plt
  })

output$barKeggAll <- downloadHandler(
  filename = "barkeggall.svg",
  content = function(file){
  ggsave(file, svg$keggAll, "svg", width = 10, units = "in") }
)

# KEGG chordiag plot all ###############
  output$keggChordAll <- renderMychordplot({
      shiny::validate(need(kgg$all, "Load file to render ChordPlot"),
               need(length(rowsAll())>3, "Select at least 4 paths of interest to render ChordPlot"))
      rowsAll <- rowsAll()
      if(is.null(rowsAll)){
          if( dim(kgg$all)[1]<10 ){rowsAll <-  seq_len(nrow(kgg$all)) }
          else{ rowsAll <-  seq_len(10)  }
      }
      mychordplot(kgg$all[rowsAll, c("Pathway","genes") ], div="keggChordAll" )
    # p <- chordPlot(kgg$all[rowsAll, ], nRows = length(rowsAll), orderby = "P.DE")
    # svg$chordAll <- list(p$x$matrix, rowsAll)
    # p
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
 # download chorplot All
#  output$chordKeggAll <- downloadHandler(
#   filename = "chordKeggAll.svg",
#   content = function(file){
#     svg(file)
#     chordDiagram(svg$chordAll[[1]], transparency = 0.3, big.gap = 1,
#                  annotationTrack = c("grid"),
#                  grid.col = colorRampPalette( 
#                    RColorBrewer::brewer.pal(11, "Spectral"))(length(svg$chordAll[[2]])) )
#     circos.track(track.index = 1, panel.fun = function(x, y) {
#     circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
#         facing = "clockwise", niceFacing = TRUE, adj = c(0, 3))},
#     bg.border = NA)
#     dev.off()
#      }
# )
 
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
    shiny::validate(need(kgg$all, "Load file and select to render Heatmap"))
    shiny::validate(need(rowsAll(), "Select the paths of interest to render HeatMap"))
    shiny::validate(need(kggDT$all, ""))
    p <- heatmapKegg(kggDT$all, rowsAll())
    svg$heatKeggAll <- list(kggDT$all, rowsAll())
    plt <- ggplotly(p)
    plt$height <- myHeightfunction( rowsAll() )
    plt$x$layout$height <- myHeightfunction(rowsAll() )
    plt
  })
 
  output$heatKeggAll <- downloadHandler(
    filename = "heatKeggAll.svg",
    content = function(file){
    p <- heatmapKegg(svg$heatKeggAll[[1]],svg$heatKeggAll[[2]] )
    ggsave(filename= file, plot = p, device = "svg", width = 10, units = "in") }
  )
  
# KEGG cnet All #################
  output$legendAll <- renderPlot({
    shiny::validate(need(kgg$all, "Load file and select to render Net Plot"))
    shiny::validate(need(rowsAll(), "Select the paths of interest to render NetPlot"))
    shiny::validate(need(kggDT$all, ""))
    visnetLegend(kggDT = kggDT$all , rows = rowsAll() )
  })
   output$keggAllNet <- renderUI({
    if(!isTRUE( input$keggAllNet_switch ) ){
      plotOutput("cnetAllKegg", height = "600px")
    } else{
      visNetworkOutput("visnetKeggAll", height = "600px")
    }
  })
  output$cnetAllKegg <- renderPlot({
    shiny::validate(need(kgg$all, "Load file and select to render Net Plot"))
    shiny::validate(need(rowsAll(), "Select the paths of interest to render NetPlot"))
    p <- customCnetKegg(kgg$all, rowsAll(), genesUp = genes$Up, genesDown = genes$Down)
    svg$cnetKeggAll <- p
    print(p)
  })
  output$visnetKeggAll <- renderVisNetwork({
    shiny::validate(need(kgg$all, "Load file and select to render Net Plot"))
    shiny::validate(need(rowsAll(), "Select the paths of interest to render NetPlot"))
    shiny::validate(need(kggDT$all, ""))
    visData <- customVisNet(kgg$all, nTerm=rowsAll(), kggDT$all,
                             up = genes$Up$SYMBOL, down = genes$Down$SYMBOL )
    visNetwork(visData$nodes, visData$edges, background = "#ffffff") %>%
    visOptions(highlightNearest = list(enabled=TRUE, hover=TRUE),
                nodesIdSelection = TRUE)
  })

  output$cnetKeggAll <- downloadHandler(
    filename = "cnetKeggAll.svg",
    content = function(file){
    ggsave(filename = file, plot = svg$cnetKeggAll, device = "svg", width = 10, height = 10, units = "in") }
  )
  
# ....................... ####
  # variables KEGG UP ##########################
  rowsUp <- reactive({input$table_rows_selected})
 # KEGG table up #####################################
  output$table <- DT::renderDataTable(server=FALSE,{
    shiny::validate(need(kgg$up, "Load file to render table"))
    names(kggDT$up)[names(kggDT$up) == "DE"] <- "DEG"
    names(kggDT$up)[names(kggDT$up) == "P.DE"] <- "p-value"
    tituloTabla <- paste0("Table: Kegg up-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="keggUp"),list(extend="excel",filename="keggUp") ),
           text="Download", filename="keggup", title=tituloTabla ) )
    
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
    shiny::validate(need(kgg$up, "Load file to render BarPlot"), 
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
      shiny::validate(need(kgg$up, "Load file to render ChordPlot"),
                      need(length(rowsUp())>3, "Select at least 3 paths of interest to render ChordPlot"))
      rowsUp <- rowsUp()
    if (is.null(rowsUp)) {
      if (dim(kgg$up)[1] < 10) {
        rowsUp <-  seq_len(nrow(kgg$up))
      }
      else{
        rowsUp <-  seq_len(10)
      }
    }
    mychordplot(kgg$up[rowsUp, c("Pathway", "genes")], div = "keggChord")
  })
 output$legendChorUp <- renderPlot(bg = "#37414b",{
    shiny::validate(need(kgg$up, "Load file to render ChordPlot"),
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
    shiny::validate(need(kgg$up, "Load file and select to render Heatmap"))
    shiny::validate(need(rowsUp(), "Select the paths of interest to render HeatMap"))
    shiny::validate(need(kggDT$up, ""))
    p <- heatmapKegg(kggDT$up, rowsUp())
    svg$heatKeggUp <- list(kggDT$up, rowsUp())
    plt <- ggplotly(p)
    plt$height <- myHeightfunction( rowsUp() )
    plt$x$layout$height <- myHeightfunction(rowsUp() )
    plt
  })

  output$heatKeggUp <- downloadHandler(
    filename = "heatKeggUp.svg",
    content = function(file){
    p <- heatmapKegg(svg$heatKeggUp[[1]],svg$heatKeggUp[[2]] )
    ggsave(filename= file, plot = p, device = "svg", width = 10, units = "in") }
  )
# KEGG cnet Up #################
  output$legendUp <- renderPlot({
    shiny::validate(need(kgg$up, "Load file and select to render Net Plot"))
    shiny::validate(need(rowsUp(), "Select the paths of interest to render NetPlot"))
    shiny::validate(need(kggDT$up, ""))
    visnetLegend(kggDT = kggDT$up , rows = rowsUp() )
  })
   output$keggUpNet <- renderUI({
    if(!isTRUE( input$keggUpNet_switch ) ){
      plotOutput("cnetKeggUp", height = "600px")
    } else{
      visNetworkOutput("visnetKeggUp", height = "600px")
    }
  })
  output$cnetKeggUp <- renderPlot({
    shiny::validate(need(kgg$up, "Load file and select to render Net Plot"))
    shiny::validate(need(rowsUp(), "Select the paths of interest to render NetPlot"))
    p <- customCnetKegg(kgg$up, rowsUp(), genesUp = data$dfilt, genesDown = NULL)
    svg$cnetKeggUp <- p
    print(p)
  })
  output$visnetKeggUp <- renderVisNetwork({
    shiny::validate(need(kgg$up, "Load file and select to render Net Plot"))
    shiny::validate(need(rowsUp(), "Select the paths of interest to render NetPlot"))
    shiny::validate(need(kggDT$up, ""))
    visData <- customVisNet(kgg$up, nTerm=rowsUp(), kggDT$up,
                             up = genes$Up$SYMBOL, down = NULL )
    visNetwork(visData$nodes, visData$edges, background = "#ffffff") %>%
    visOptions(highlightNearest = list(enabled=TRUE, hover=TRUE),
                nodesIdSelection = TRUE)
  })
  
  output$cnetkeggUp <- downloadHandler(
    filename = "cnetKeggUp.svg",
    content = function(file){
    ggsave(filename = file, plot = svg$cnetKeggUp, device = "svg", width = 10, height = 10, units = "in") }
  )
  # ....................... ####
  # variables KEGG Down ##########################
  rowsDown <- reactive({input$tableDown_rows_selected})
 # KEGG table down #####################################
  output$tableDown <- DT::renderDataTable(server=FALSE,{
    shiny::validate(need(kgg$down, "Load file to render table"))
    names(kggDT$down)[names(kggDT$down) == "DE"] <- "DEG"
    names(kggDT$down)[names(kggDT$down) == "P.DE"] <- "p-value"
    tituloTabla <- paste0("Table: Kegg down-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="keggDown"),list(extend="excel",filename="keggDown") ),
           text="Download", filename="keggdown", title=tituloTabla ) )
    
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
  # KEGG barplot down################
  output$keggPlotDown <- renderPlotly ({
    shiny::validate(need(kgg$down, "Load file to render BarPlot"),
                    need(rowsDown(), "Select the paths of interest to render BarPlot"))
    rowsDown <- rowsDown()
    if(is.null(rowsDown)){
        if( dim(kgg$down)[1]<10 ){rowsDown <-  seq_len(nrow(kgg$down)) }
        else{ rowsDown <-  seq_len(10)  }
        }
    p <- plotKegg(enrichdf = kgg$down[rowsDown,], nrows = length(rowsDown), colors = c(input$downColor))
    svg$barkeggdown <- p[[2]] 
    plt <- p[[1]]
    plt$height <- myHeightfunction( rowsDown() )
    plt$x$layout$height <- myHeightfunction(rowsDown() )
    plt

  })
  
  output$barKeggDown <- downloadHandler(
    filename = "barkeggdown.svg",
    content = function(file){
    ggsave(file, svg$barkeggdown, "svg", width = 10, height = 10, units = "in") }
  )
  
  # KEGG chordiag plot down ###############
  output$keggChordDown <- renderMychordplot({
      shiny::validate(need(kgg$down, "Load file to render ChordPlot"),
                      need(length(rowsDown())>3, "Select at least 3 paths of interest to render ChordPlot"))
    rowsDown<- rowsDown()
    if(is.null(rowsDown)){
        if( dim(kgg$down)[1]<10 ){rowsDown <-  seq_len(nrow(kgg$down)) }
        else{ rowsDown <-  seq_len(10)  }
    }
    mychordplot(kgg$down[rowsDown, c("Pathway", "genes")], div = "keggChordDown")
  })
 output$legendChorDown <- renderPlot(bg = "#37414b",{
    shiny::validate(need(kgg$down, "Load file to render ChordPlot"),
                    need(length(rowsDown())>3, ""))
    rowsDown <- rowsDown()
    if(is.null(rowsDown)){
        if (dim(kgg$down)[1] < 10) {rowsDown <-  seq_len(nrow(kgg$down))}
        else{rowsDown <-  seq_len(10)}
        
        }
    legendChorplot(kgg$down[rowsDown, ] )
  })
  
 
  # KEGG dotplot Down ################### 
  output$keggDotDown <- renderPlotly({
    validate(need(kgg$down, "Load file to render DotPlot"),
             need(rowsDown(), "Select the paths of interest to render dotPlot"))
    rowsDown <- rowsDown()
    if(is.null(rowsDown)){
      if( dim(kgg$down)[1]<10 ){rowsDown <-  seq_len(nrow(kgg$down)) }
      else{ rowsDown <-  seq_len(10)  }
    }
    plt <- dotPlotkegg(kgg$down[rowsDown,], n = length(rowsDown))
    svg$dotKeggDown <- plt
    plt <- ggplotly(plt)
    plt$height <- myHeightfunction( rowsDown() )
    plt$x$layout$height <- myHeightfunction(rowsDown() )
    plt
  })
 
  output$dotKeggDown <- downloadHandler(
    filename = "dotKeggDown.svg",
    content = function(file){
    ggsave(file, svg$dotKeggDown, device = "svg", width = 10, units = "in") }
  )
  # KEGG heatmap Down #################
  output$heatmapKeggDown <- renderPlotly({
    shiny::validate(need(kgg$down, "Load file and select to render Heatmap"))
    shiny::validate(need(rowsDown(), "Select the paths of interest to render HeatMap"))
    shiny::validate(need(kggDT$down, ""))
    p <- heatmapKegg(kggDT$down, rowsDown())
    svg$heatKeggDown <- list(kggDT$down, rowsDown())
    plt <- ggplotly(p)
    plt$height <- myHeightfunction( rowsDown() )
    plt$x$layout$height <- myHeightfunction(rowsDown() )
    plt
  })

  output$heatKeggDown <- downloadHandler(
    filename = "heatKeggDown.svg",
    content = function(file){
    p <- heatmapKegg(svg$heatKeggDown[[1]],svg$heatKeggDown[[2]] )
    ggsave(filename= file, plot = p, device = "svg", width = 10, units = "in") }
  )
 
# KEGG cnet Down #################
  output$legendDown <- renderPlot({
    shiny::validate(need(kgg$down, "Load file and select to render Net Plot"))
    shiny::validate(need(rowsDown(), "Select the paths of interest to render NetPlot"))
    shiny::validate(need(kggDT$down, ""))
    visnetLegend(kggDT = kggDT$down , rows = rowsDown() )
  })
   output$keggDownNet <- renderUI({
    if(!isTRUE( input$keggDownNet_switch ) ){
      plotOutput("cnetKeggDown", height = "600px")
    } else{
      visNetworkOutput("visnetKeggDown", height = "600px")
    }
  })
  output$cnetKeggDown <- renderPlot({
    shiny::validate(need(kgg$down, "Load file and select to render Net Plot"))
    shiny::validate(need(rowsDown(), "Select the paths of interest to render NetPlot"))
    p <- customCnetKegg(kgg$down, rowsDown(), genesDown = data$dfilt, genesUp = NULL)
    svg$cnetKeggDown <- p
    print(p)
    })
  output$visnetKeggDown <- renderVisNetwork({
    shiny::validate(need(kgg$down, "Load file and select to render Net Plot"))
    shiny::validate(need(rowsDown(), "Select the paths of interest to render NetPlot"))
    shiny::validate(need(kggDT$down, ""))
    visData <- customVisNet(kgg$down, nTerm=rowsDown(), kggDT$down,
                             down = genes$Down$SYMBOL, up = NULL )
    visNetwork(visData$nodes, visData$edges, background = "#ffffff") %>%
    visOptions(highlightNearest = list(enabled=TRUE, hover=TRUE),
                nodesIdSelection = TRUE)
  })
  output$cnetkeggDown <- downloadHandler(
    filename = "cnetKeggDown.svg",
    content = function(file){
    ggsave(filename = file, plot = svg$cnetKeggDown, device = "svg", width = 10, height = 10, units = "in") }
  )
 # ....................... ####
 # variables GO ###################################
  bprowsall <- reactive({input$tableBPall_rows_selected}) 
  mfrowsall <- reactive({input$tableMFall_rows_selected})
  ccrowsall <- reactive({input$tableCCall_rows_selected})
  
  bprowsup <- reactive({input$tableBP_rows_selected})
  mfrowsup <- reactive({input$tableMF_rows_selected})
  ccrowsup <- reactive({input$tableCC_rows_selected})
  
  bprowsdown <- reactive({input$tableBPdown_rows_selected})
  mfrowsdown <- reactive({input$tableMFdown_rows_selected})
  ccrowsdown <- reactive({input$tableCCdown_rows_selected})
  
  typeBarBpAll <- reactive({input$selectbpall})
  typeBarMfAll <- reactive({input$selectmfall})
  typeBarCcAll <- reactive({input$selectccall})
  
 # ....................... ####
   # GO table BP ALL #####################
  output$tableBPall <- DT::renderDataTable(server=FALSE,{
    shiny::validate(need(goDT$all, "Load file to render table"))
    goDT <- goDT$all 
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level) 
    tituloTabla <- paste0("Table: GO-BP all genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="BPall"),list(extend="excel",filename="BPall") ),
           text="Download", filename="BPall", title=tituloTabla ) )
    
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
    shiny::validate(need(go$all, "Load file to render plot"),
                    need(bprowsall(), "Select at least one row to plot") )
    bprowsall <- bprowsall()
    gosBP <- go$all[go$all$Ont=="BP",]
    if(is.null(bprowsall)){
        if( dim(gosBP)[1]<10 ){bprowsall <-  seq_len(nrow(gosBP)) }
        else{ bprowsall <-  seq_len(10)  }
    }
    if(isTRUE(df3cols$TF)){
          p <- plotGOAll(enrichdf = gosBP[bprowsall, ], nrows = length(bprowsall), ont="BP", 
                    genesUp = genes$Up, genesDown = genes$Down,
                    colors = c(input$downColor, input$upColor))
          if( typeBarBpAll() == "Dodge") {
            plt <- p[[1]]; svg$barbpall <- p[[1]]
            }
          else if ( typeBarBpAll() == "Stack") {
            plt <- p[[2]]; svg$barbpall <- p[[2]]
            }
          else {
            plt <- p[[3]] ;svg$barbpall <- p[[3]] 
          }
    } else{
      p <- plotGO(enrichdf = gosBP[bprowsall, ], nrows = length(bprowsall), ont="BP",
               colors = "#045a8d" )
        svg$barbpall <- p
        plt <- p
    }
    
    plt <- plt %>% plotly::ggplotly(tooltip = "all" )
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
    shiny::validate(need(go$all, "Load file to render dotPlot"), 
                    need(length(bprowsall())>=2, "Select at least 2 row" ) )
    bprowsall <- bprowsall()
    p <- goBarplot(enrichGO = go$all, resGO = data$dfilt, genes= genes$all,
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
    shiny::validate(need(go$all, "Load file to render dotPlot"))
    shiny::validate(need(data$dfilt,""))
    shiny::validate(need( length(bprowsall())>=4 , "Select at least 4 rows"))
    bprowsall <- bprowsall()
    if(length(bprowsall)>=4){
      go <- go$all[go$all$Ont=="BP",]
      circ <- data2circle(go=go[bprowsall, ], res=data$dfilt, genes=genes$all)
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
  output$cloudBPAll <- renderPlotly({
    validate(need(go$all, "Load file to render dotPlot"))
    goall <- go$all[go$all$Ont=="BP"& go$all$level>=input$bpallLevel, ]
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
  # ...................... #############
  # GO table MF all #####################
  output$tableMFall <- DT::renderDataTable(server=FALSE,{
    shiny::validate(need(goDT$all, "Load file to render table"))
    goDT <- goDT$all
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-MF all genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="MFall"),list(extend="excel",filename="MFall") ),
           text="Download", filename="MFall", title=tituloTabla ) )
        
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
    shiny::validate(need(go$all, "Load file to render plot"),
                    need(mfrowsall(), "Select at least one row to plot") )
    mfrowsall <- mfrowsall()
    gosMF <- go$all[go$all$Ont=="MF",]
    if(is.null(mfrowsall)){
        if( dim(gosMF)[1]<10 ){mfrowsall <-  seq_len(nrow(gosMF)) }
        else{ mfrowsall <-  seq_len(10)  }
    }
    if(isTRUE(df3cols$TF)){
        p <- plotGOAll(enrichdf = gosMF[mfrowsall, ], nrows = length(mfrowsall), ont="MF", 
                       genesUp = genes$Up, genesDown = genes$Down,
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
    }
    else{
      p <- plotGO(enrichdf = gosMF[mfrowsall, ], nrows = length(mfrowsall), ont="MF",
               colors = "#045a8d" )
      svg$barmfall <- p
      plt <- p
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
    shiny::validate(need(go$all, "Load file to render dotPlot"),
                    need(length(mfrowsall())>=2, "Select at least 2 row"))
    mfrowsall <- mfrowsall()
    p <- goBarplot(enrichGO = go$all, resGO = data$dfilt, genes= genes$all,
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
    shiny::validate(need(go$all, "Load file to render dotPlot"))
    shiny::validate(need(data$dfilt,""))
    shiny::validate(need( length(mfrowsall())>=4 , "Select at least 4 rows"))
    mfrowsall <- mfrowsall()
    if(length(mfrowsall)>=4){
        go <- go$all[go$all$Ont=="MF",]
      circ <- data2circle(go=go[mfrowsall, ], res=data$dfilt, genes=genes$all)
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
    shiny::validate(need(goDT$all, "Load file to render table"))
    goDT <- goDT$all
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-CC all genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="CCall"),list(extend="excel",filename="CCall") ),
           text="Download", filename="CCall", title=tituloTabla ) )
    
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
    shiny::validate(need(go$all, "Load file to render plot"),
                    need(ccrowsall(), "Select at least one row to plot") )
    ccrowsall <- ccrowsall()
    gosCC <- go$all[go$all$Ont=="CC",]
    if(is.null(ccrowsall)){
        if( dim(gosCC)[1]<10 ){ccrowsall <-  seq_len(nrow(gosCC)) }
        else{ ccrowsall <-  seq_len(10)  }
    }
    if(isTRUE(df3cols$TF)){
    p <- plotGOAll(enrichdf = gosCC[ccrowsall, ], nrows = length(ccrowsall), ont="CC", 
                   genesUp = genes$Up, genesDown = genes$Down,
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
    }
    else{
        p <- plotGO(enrichdf = gosCC[ccrowsall, ], nrows = length(ccrowsall), ont="CC",
               colors = "#045a8d" )
        svg$barccall <- p
        plt <- p
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
    shiny::validate(need(go$all, "Load file to render dotPlot"))
    ccrowsall <- ccrowsall()
    p <- goBarplot(enrichGO = go$all, resGO = data$dfilt, genes= genes$all,
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
    shiny::validate(need(go$all, "Load file to render dotPlot"))
    shiny::validate(need(data$dfilt,""))
    shiny::validate(need( length(ccrowsall())>=4 , "Select at least 4 rows"))
    ccrowsall <- ccrowsall()
    if(length(ccrowsall)>=4){
        go <- go$all[go$all$Ont=="CC",]
      circ <- data2circle(go=go[ccrowsall, ], res=data$dfilt, genes=genes$all)
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
    shiny::validate(need(goDT$up, "Load file to render table"))
    goDT <- goDT$up
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-BP up-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="BPup"),list(extend="excel",filename="BPup") ),
           text="Download", filename="BPup", title=tituloTabla ) )
    
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
    shiny::validate(need(go$up, "Load file to render plot"),
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
    shiny::validate(need(go$up, "Load file to render dotPlot"),
                    need( length(bprowsup())>=2, "Select at least 2 row") )
    bprowsup <- bprowsup()
    p <- goBarplot(enrichGO = go$up, resGO = data$dfilt, genes= genes$Up,
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
    shiny::validate(need(go$up, "Load file to render dotPlot"))
    shiny::validate(need(data$dfilt,""))
    shiny::validate(need( length(bprowsup())>=4 , "Select at least 4 rows"))
    bprowsup <- bprowsup()
    if(length(bprowsup)>=4){
        go <- go$up[go$up$Ont=="BP",]
      circ <- data2circle(go=go[bprowsup, ], res=data$dfilt, genes=genes$Up)
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
  output$tableMF <- DT::renderDataTable(server=FALSE,{
    shiny::validate(need(goDT$up, "Load file to render table"))
    goDT <- goDT$up
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-MF up-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="MFup"),list(extend="excel",filename="MFup") ),
           text="Download", filename="MFup", title=tituloTabla ) )
  
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
    shiny::validate(need(go$up, "Load file to render plot"),
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
    shiny::validate(need(go$up, "Load file to render dotPlot"),
                    need(length(mfrowsup())>=2, "Select at least 2 row") )
    mfrowsup <- mfrowsup()
    p <- goBarplot(enrichGO = go$up, resGO = data$dfilt, genes= genes$Up,
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
    shiny::validate(need(go$up, "Load file to render dotPlot"))
    shiny::validate(need(data$dfilt,""))
    shiny::validate(need(  length(mfrowsup())>=4  , "Select at least 4 rows"))
    mfrowsup <- mfrowsup()
    if(length(mfrowsup)>=4){
        go <- go$up[go$up$Ont=="MF",]
      circ <- data2circle(go=go[mfrowsup, ], res=data$dfilt, genes=genes$Up)
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
    shiny::validate(need(goDT$up, "Load file to render table"))
    goDT <- goDT$up
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-CC up-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="CCup"),list(extend="excel",filename="CCup") ),
           text="Download", filename="CCup", title=tituloTabla ) )

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
    shiny::validate(need(go$up, "Load file to render plot"),
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
    shiny::validate(need(go$up, "Load file to render dotPlot"),
                    need(length(ccrowsup())>=2, "Select at least 2 row"))
    ccrowsup <- ccrowsup()
    p <- goBarplot(enrichGO = go$up, resGO = data$dfilt, genes= genes$Up,
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
    shiny::validate(need(go$up, "Load file to render dotPlot"))
    shiny::validate(need(data$dfilt,""))
    shiny::validate(need( length(ccrowsup())>=4 , "Select at least 4 rows"))
    ccrowsup <- ccrowsup()
    if(length(ccrowsup)>=4){
        go <- go$up[go$up$Ont=="CC",]
      circ <- data2circle(go=go[ccrowsup, ], res=data$dfilt, genes=genes$Up)
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
    shiny::validate(need(goDT$down, "Load file to render table"))
    goDT <- goDT$down
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-BP down-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="BPdown"),list(extend="excel",filename="BPdown") ),
           text="Download", filename="BPdown", title=tituloTabla ) )
   
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
    shiny::validate(need(go$down, "Load file to render plot"),
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
    shiny::validate(need(go$down, "Load file to render dotPlot"),
                    need(length(bprowsdown())>=2, "Select at least 2 row"))
    bprowsdown <- bprowsdown()
    p <- goBarplot(enrichGO = go$down, resGO = data$dfilt, genes= genes$Down,
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
    shiny::validate(need(go$down, "Load file to render dotPlot"))
    shiny::validate(need(data$dfilt,""))
    shiny::validate(need(  length(bprowsdown())>=4 , "Select at least 4 rows"))
    bprowsdown <- bprowsdown()
    if(length(bprowsdown)>=4){
        go <- go$down[go$down$Ont=="BP",]
      circ <- data2circle(go=go[bprowsdown, ], res=data$dfilt, genes=genes$Down)
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
  output$tableMFdown <- DT::renderDataTable(server=FALSE,{
    shiny::validate(need(goDT$down, "Load file to render table"))
    goDT <- goDT$down
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-MF down-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection",
           buttons = list( list(extend="csv",filename="MFdown"),list(extend="excel",filename="MFdown") ),
           text="Download", filename="MFdown", title=tituloTabla ) )

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
    shiny::validate(need(go$down, "Load file to render plot"),
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
    shiny::validate(need(go$down, "Load file to render dotPlot"),
                    need(length(mfrowsdown())>=2, "Select at least 2 row"))
    mfrowsdown <- mfrowsdown()
    p <- goBarplot(enrichGO = go$down, resGO = data$dfilt, genes= genes$Down,
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
    shiny::validate(need(go$down, "Load file to render dotPlot"))
    shiny::validate(need(data$dfilt,""))
    shiny::validate(need( length(mfrowsdown())>=4 , "Select at least 4 rows"))
    mfrowsdown <- mfrowsdown()
    if(length(mfrowsdown)>=4){
        go <- go$down[go$down$Ont=="MF",]
      circ <- data2circle(go=go[mfrowsdown, ], res=data$dfilt, genes=genes$Down)
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
    shiny::validate(need(goDT$down, "Load file to render table"))
    goDT <- goDT$down
    names(goDT)[names(goDT) == "DE"] <- "DEG"
    names(goDT)[names(goDT) == "P.DE"] <- "p-value"
    names(goDT)[names(goDT) == "level"] <- "Ont.level"
    goDT$Ont.level = as.integer(goDT$Ont.level)
    tituloTabla <- paste0("Table: GO-CC down-regulated genes | ","log2FC: ",logfc()[1],"_",logfc()[2]," | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <-  list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="CCdown"),list(extend="excel",filename="CCdown") ),
           text="Download", filename="CCdown", title=tituloTabla ) )

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
    shiny::validate(need(go$down, "Load file to render plot"),
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
    shiny::validate(need(go$down, "Load file to render dotPlot"),
                    need(length(ccrowsdown())>=2, "Select at least 2 row"))
    ccrowsdown <- ccrowsdown()
    p <- goBarplot(enrichGO = go$down, resGO = data$dfilt, genes= genes$Down,
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
    shiny::validate(need(go$down, "Load file to render dotPlot"))
    shiny::validate(need(data$dfilt,""))
    shiny::validate(need( length(ccrowsdown() )>=4 , "Select at least 4 rows"))
    ccrowsdown <- ccrowsdown()
    if(length(ccrowsdown)>=4){
        go <- go$down[go$down$Ont=="CC",]
      circ <- data2circle(go=go[ccrowsdown, ], res=data$dfilt, genes=genes$Down)
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
      godown <- go$down[go$down$Ont=="CC" & go$down$level>=input$ccdownLevel, ]
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
  output$gseaTable <- DT::renderDataTable(server=FALSE,{
    shiny::validate(need(data$dfilt, "Load file to render table"))
    shiny::validate(need(input$gseadb!="","Select dataset"))
    gsea$gsea <- gseaKegg(data$dfilt[, c("ENTREZID","logFC")], specie(), gseadb = input$gseadb )
    mygsea <- gsea$gsea
    if( length(which(mygsea@result$p.adjust<=0.05)) == 0 ){
        createAlert(session, anchorId = "gsea", title = "Oops!!", 
          content = "Sorry, I didn't get any significant results for this analysis",
          append=FALSE, style = "info")
    } else{
    table <- mygsea@result[mygsea@result$p.adjust<=0.05 ,2:9] %>% 
      mutate_at(vars(3:7), ~round(., 4))

    tituloTabla <- paste0("Table: GSEA pathway | ",
                          "log2FC: ",logfc()[1],"_",logfc()[2],
                          " | ","padj: ",padj()," | ",
                          "Num genes Up/down: ",numgenesDE$up,"/",numgenesDE$down)
    customButtons <- list(
      list(extend = "copy", title=tituloTabla),
      list(extend="collection", 
           buttons = list( list(extend="csv",filename="GSEAkegg"),list(extend="excel",filename="GSEAkegg") ),
           text="Download", filename="GSEAkegg", title=tituloTabla ) )
    
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
  # GSEA plot ##########################
  output$gseaPlot <- renderPlot({
    shiny::validate(need(gsea$gsea, "Load file to render table"))
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

  output$report <- renderUI({
      actionButton("report1", "html report")
  })
  
  observeEvent(input$report1, {
     if(isTRUE(enrichflag$one)){
         enrichflag$three <- FALSE
         showModal(popupModal1())
     }else{
        if(isTRUE(enrichflag$three)){
            enrichflag$one <- FALSE
            showModal(popupModal3())
        }
     }
    })
  
  observeEvent(input$unselect3,{
    if(input$unselect3 >0){
      if(input$unselect3 %% 2 == 0 ){
        selectPopUpModal3(session = session)
      }else{
        unselectPopUpModal3(session = session)
      }
    }
  })
  
  # observeEvent(input$report1, {
  #     showModal(popupModal1())
  #   })

  observeEvent(input$unselect1,{
    if(input$unselect1 >0){
      if(input$unselect1 %% 2 == 0 ){
        selectPopUpModal1(session = session)
      }else{
        unselectPopUpModal1(session = session)
      }
    }
  })

  applyPress <- reactiveValues(ok=FALSE)
  observeEvent(input$ok,{
        applyPress$ok <- TRUE
        if(isTRUE(enrichflag$three)){
          vals$preview <- input$modalPreview
          vals$keggUp <- input$modalkeggUp
          vals$keggDown <- input$modalkeggDown
          vals$GOUp <- input$modalGOUp
          vals$GODown <- input$modalGODown
          vals$GSEA <- input$modalGSEA
        }
        vals$keggAll <- input$modalkeggAll
        vals$GOAll <- input$modalGOAll
        #removeModal()
  })
 
  output$downloadhtml <- renderUI({
    shiny::validate(need(isTRUE(applyPress$ok), ""))
    downloadButton("download", "Download report")
    })
 ## report ##################
    output$download <- downloadHandler(
    filename = "report.html",
    content = function(file) {
      removeModal()
      applyPress$ok <- FALSE
      tempReport <- file.path(tempdir(), "report.Rmd")
      if( dim(data$df)[2]==5 ){
      file.copy("report.Rmd", tempReport, overwrite = TRUE) }else{
        file.copy("report1col.Rmd", tempReport, overwrite = TRUE)
      }
      file.copy("report.css", file.path(tempdir(), "report.css"), overwrite = TRUE)
      file.copy("utilsReport.R", file.path(tempdir(),"utils.R"), overwrite = TRUE)
      file.copy("resources/", tempdir(), overwrite = TRUE, recursive = TRUE)
      file.copy("resources/dna-svg-small-13.gif",
      file.path(tempdir(), "resources/dna-svg-small-13.gif"), overwrite = TRUE)
      ## inicializar variables preview
      volcObj <- karyObj <- FALSE
      ## inicializar variables kegg
      tablekgaObj <- barkgaObj <- chorkgaObj <- dotkgaObj <- heatkgaObj <- netkgaObj <- FALSE
      tablekguObj <- barkguObj <- chorkguObj <- dotkguObj <- heatkguObj <- netkguObj <- FALSE
      tablekgdObj <- barkgdObj <- chorkgdObj <- dotkgdObj <- heatkgdObj <- netkgdObj <- FALSE
      ## inicializar variables Go
      cloudgoaObj <- tablegoaObj <- bargoaObj <- dotgoaObj <- gobargoaObj <- gocirclegoaObj <- FALSE
      cloudgouObj <- tablegouObj <- bargouObj <- dotgouObj <- gobargouObj <- gocirclegouObj <- FALSE
      cloudgodObj <- tablegodObj <- bargodObj <- dotgodObj <- gobargodObj <- gocirclegodObj <- FALSE
      ## inicializar variables GSEA
      tablegseaObj <- plotgseaObj <- FALSE 
      bprowsall <- bprowsall(); mfrowsall <- mfrowsall(); ccrowsall <- ccrowsall()
      bprowsup <- bprowsup(); mfrowsup <- mfrowsup(); ccrowsup <- ccrowsup()
      bprowsdown <- bprowsdown(); mfrowsdown <- mfrowsdown(); ccrowsdown <- ccrowsdown()
      gsearow <- gsearow()

          nrowsall <- rowsAll()
          if(!is.null(kggDT$all)){
            if(is.null(nrowsall)){ 
              nrowsall <-  ( if( dim(kggDT$all)[1]<10) seq_len(nrow(kggDT$all))
                             else seq_len(10) ) }
          }
          nrowsup <- rowsUp()
          if(!is.null(kggDT$up)){
            if(is.null(nrowsup)){ 
              nrowsup <-  ( if( dim(kggDT$up)[1]<10) seq_len(nrow(kggDT$up))
                            else seq_len(10) ) }
          }
          nrowsdown <- rowsDown()
          if(!is.null(kggDT$down)){
            if(is.null(nrowsdown)){ 
              nrowsdown <-  ( if( dim(kggDT$down)[1]<10) seq_len(nrow(kggDT$down))
                              else seq_len(10) ) }
          }
      if(!is.null(vals$preview)){      #para preview
        if("Volcano" %in% vals$preview){volcObj <- TRUE}
        if("Karyoplot" %in% vals$preview){ karyObj <- TRUE}
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
        if("WordCloud" %in% vals$GOAll){cloudgoaObj <- TRUE}
        if("Table" %in% vals$GOAll){ tablegoaObj <- TRUE }
        if("Barplot" %in% vals$GOAll){ bargoaObj <- TRUE }
        if("Dotplot" %in% vals$GOAll){ dotgoaObj <- TRUE }
        if("GObarplot" %in% vals$GOAll){ gobargoaObj <- TRUE }
        if("GOcircleplot" %in% vals$GOAll){ gocirclegoaObj <- TRUE }
      }
      if(!is.null(vals$GOUp)){#para GoUp
        if("WordCloud" %in% vals$GOUp){cloudgouObj <- TRUE}
        if("Table" %in% vals$GOUp){ tablegouObj <- TRUE }
        if("Barplot" %in% vals$GOUp){ bargouObj <- TRUE }
        if("Dotplot" %in% vals$GOUp){ dotgouObj <- TRUE }
        if("GObarplot" %in% vals$GOUp){ gobargouObj <- TRUE }
        if("GOcircleplot" %in% vals$GOUp){ gocirclegouObj <- TRUE }
      }
      if(!is.null(vals$GODown)){#para GoDown
        if("WordCloud" %in% vals$GODown){cloudgodObj <- TRUE}
        if("Table" %in% vals$GODown){ tablegodObj <- TRUE }
        if("Barplot" %in% vals$GODown){ bargodObj <- TRUE }
        if("Dotplot" %in% vals$GODown){ dotgodObj <- TRUE }
        if("GObarplot" %in% vals$GODown){ gobargodObj <- TRUE }
        if("GOcircleplot" %in% vals$GODown){ gocirclegodObj <- TRUE }
      }
      if(!is.null(vals$GSEA)){#para GSEA
        if(is.null(gsea$gsea)){
            tablegseaObj <- FALSE; plotgseaObj <- FALSE  
        }else{
        if( length(which(gsea$gsea@result$p.adjust<=0.05)) == 0 ){
            tablegseaObj <- FALSE; plotgseaObj <- FALSE  
        }else{
          if("Table" %in% vals$GSEA){ tablegseaObj <- TRUE}
          if("GSEA plot" %in% vals$GSEA){ plotgseaObj <- TRUE}
        }
        }
          }

      params <- list( values = vals, datadf = data$df, datadfilt = data$dfilt, annotation=annotation(),
                     specie = specie(), padj =padj(), logfc = logfc(),
                     volcObj = volcObj, genesvolcano = genesVolcano(), 
                     upcolor = input$upColor, downcolor = input$downColor, 
                     karyObj = karyObj,
                     datagenesup = genes$Up, datagenesdown = genes$Down,
                     tablekgaObj = tablekgaObj, kggall = kgg$all, genesdedown = numgenesDE$down,
                     genesdeup = numgenesDE$up, kggdtall = kggDT$all,
                     barkgaObj = barkgaObj, nrowsall = nrowsall, typebarkeggall = typeBarKeggAll(),
                     chorkgaObj = chorkgaObj, dotkgaObj = dotkgaObj, heatkgaObj = heatkgaObj,
                     netkgaObj = netkgaObj, tablekguObj = tablekguObj, barkguObj = barkguObj,
                     chorkguObj = chorkguObj, dotkguObj =dotkguObj, heatkguObj = heatkguObj,
                     netkguObj = netkguObj, tablekgdObj = tablekgdObj, barkgdObj = barkgdObj,
                     chorkgdObj = chorkgdObj, dotkgdObj = dotkgdObj, heatkgdObj = heatkgdObj,
                     cloudgoaObj = cloudgoaObj, cloudgouObj = cloudgouObj, cloudgodObj = cloudgodObj,
                     netkgdObj = netkgdObj, kggup = kgg$up, kggdown = kgg$down, kggdtup = kggDT$up, 
                     kggdtdown = kggDT$down, nrowsup = nrowsup, nrowsdown = nrowsdown, 
                     typebarbpall=typeBarBpAll(), typebarmfall=typeBarMfAll(),
                     typebarccall=typeBarCcAll(),
                     tablegoaObj = tablegoaObj, bargoaObj=bargoaObj, dotgoaObj=dotgoaObj,
                     gobargoaObj=gobargoaObj,gocirclegoaObj=gocirclegoaObj, tablegouObj = tablegouObj,
                     bargouObj=bargouObj, dotgouObj=dotgouObj, gobargouObj=gobargouObj,
                     gocirclegouObj=gocirclegouObj, tablegodObj = tablegodObj, bargodObj=bargodObj,
                     dotgodObj=dotgodObj, gobargodObj=gobargodObj, gocirclegodObj=gocirclegodObj,
                     goall = go$all, godtall=goDT$all, goup = go$up, godtup=goDT$up,
                     godown = go$down, godtdown=goDT$down,
                     bprowsall=bprowsall, mfrowsall=mfrowsall, ccrowsall=ccrowsall,
                     bprowsup=bprowsup, mfrowsup=mfrowsup, ccrowsup=ccrowsup,
                     bprowsdown=bprowsdown, mfrowsdown=mfrowsdown, ccrowsdown=ccrowsdown,
                     gsearow = gsearow, gseagsea = gsea$gsea, tablegseaObj = tablegseaObj,
                     plotgseaObj = plotgseaObj, textnotes = input$textNotes)
      
      params <- c(params, list(tempdir=tempdir() ))
      rmarkdown::render(
        tempReport,
        output_file = file,
        params = params,
        envir = new.env(parent = globalenv( ))
      )
    } )
  

  
     
}
shinyApp(ui, server)
