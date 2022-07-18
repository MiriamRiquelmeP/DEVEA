library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(shinyjs)
library(networkD3)
library(dplyr)
library(tmap)
library(jsonlite)
library(datasets)
source("global.R")
data("World")
dataWorld <- mapData()
ui <- dashboardPage(
  dashboardHeader(title = "DEVEA shiny app lobby", titleWidth = "300px",
                  tags$li(class = "dropdown", actionButton("statButton", "Stats",
                                   style="margin-top:8px; margin-right: 15px"))
                  ),
  dashboardSidebar(disable = FALSE),
  dashboardBody(
      setShadow(class = "box"),
      tags$head(tags$style(HTML('
                                .bg-limon {
                                background-color: #295e7d;
                                }
                                .box-body {
                                color: #ffffff;
                                }'))),
            useShinyjs(),
      fluidRow(
          column(width = 6,
                 flipBox(
                     width = 12,
                     id = 1,
                     front = div(
                       style = "width: 95%; margin: auto;",
                       class = "text-center",
                       h2("Transcriptomics interactive tools set"),
                       p(tags$br(),
                         tags$br(),
                         tags$b("DEVEA"), "is an R shiny application tool developed to ", 
                         HTML(paste0("perform ", '<B>D</B>',"ifferential ", '<B>E</B>',"xpression, ")),
                         HTML(paste0('<B>V</B>',"isualisation and ", '<B>E</B>',"nrichment ", '<B>A</B>',"nalysis."))),
                       tags$br(),
                       p("Its intuitive and easy-to-manipulate interface facilitates gene expression 
                         visualization, statistical comparisons and further meta-analysis such as 
                         enrichment analysis, without bioinformatics expertise."),
                       tags$br(),
                       p("DEVEA performs an extended analysis from different input formats at distinct 
                         stages on transcriptomics data, producing a wide variety of dynamic exploratory graphs  
                         and statistical results from different comparisons of interest. Moreover, it 
                         generates an extensive pathway analysis from the selected set of significant 
                         features with interactive tables and plots. 
                         Finally, a thorough and customizable HTML report can be extracted for 
                         further exploration outside the software."),
                       p("The application can be run from the two links in the boxes below."),
                       p("For details, please see our ", a("paper", href="https://f1000research.com/articles/11-711"),
                         "and a detailed ", a("demo.", href="https://shiny.imib.es/DESeqDevea/tutorial.html"),
                       "Source code available at our ", a("GitHub repository.", href="https://github.com/MiriamRiquelmeP/DEVEA"),
                       tags$br(),
                         p("Check all the possible outcomes extracted according to the type of input data 
                         and the overall DEVEA workflow in the graphs on the right."), 
                         p("Click", "here", style = "color:steelblue", "to discover more about us!"),
                         tags$br(),),
                     ),
                     back = div(
                       class = "text-center",
                       h1("About us"),
                       tagList(
                         column(
                           width = 12,
                           align = "center",
                           HTML("<img src='dna-svg-small-13.gif', width='100px'><br>
                  <h4>Main authors:<br><br>
    Miriam Riquelme-Perez 
    <a href='https://www.linkedin.com/in/miriam-riquelme-perez/' target='_blank'> 
    <img src='linkedin_little.svg'> </a> <a href='mailto:miriam.riquelmep@gmail.com'>
    <img src='email.svg'></a><br>
    Fernando Pérez Sanz 
    <a href='https://www.linkedin.com/in/fernandoperez72/' target='_blank'> 
    <img src='linkedin_little.svg'> 
    </a> <a href='mailto:fernando.perez@ffis.es'> <img src='email.svg'></a></h4><br>
    For any suggestion or bug, please contact us.")
                         )
                       )
                     )
                 )
          ),
          column(width = 6, align = "center",
                 
                 flipBox(width = 12, id = 2, 
                         front = HTML("<img src='schema_final.png', height = '525px'>"),
                         back = HTML("<img src='workflow.png', height = '525px'")
                        )
       )
      ),
      fluidRow(width =12, column(width = 6, offset = 0, style='padding:10px;')),
      fluidRow(
        column(width =12, offset = 2,
          miBoxPlus(
              title = "Go DESeq DEVEA",
              closable = FALSE,
              width = 4,
              background = "limon",
              height = "250px",
              tags$p("Enter a counting matrix (CM) associated with sample information (SI) or a DESeqDataSet object (DO) 
                     with the designs or contrasts of interest to run the analysis. 
                     Performs a complete differential expression analysis, data visualization and enrichment analysis 
                     of transcriptomics data based on ", a("DESeq2", href="https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html"),"
                     calculation, Kyoto Encyclopedia of Genes and Genomes (", a("Kegg", href="https://www.genome.jp/kegg/pathway.html"), 
                     ") pathways database, Gene Ontology (", a("GO", href="http://geneontology.org"), ") 
                     terms resource and Gene Set Enrichment Analysis (", a("GSEA", href="https://www.gsea-msigdb.org/gsea/index.jsp"),")."),
              tags$br(),
              fluidRow(
                  column(
                      align = "center",
                      width = 2,
                      dropdown(
                          style = "unite",
                          icon = icon("info"),
                          status = "primary",
                          width = "600px",
                          animate = animateOptions(
                              enter = animations$fading_entrances$fadeInRightBig,
                              exit = animations$fading_exits$fadeOutLeftBig
                          ),
                          box(
                            style = "width: 90%; margin: auto;",
                              width = 12,
                              solidHeader = TRUE,
                              status = "primary",
                              height = "150px",
                              title = tags$strong("DESeq DEVEA 2022"),
                              background = "light-blue",
                              tags$h4(
                                  "File formats allowed: .txt, .tsv or .xlsx files containing the counting matrix plus 
                                  the sample data information or a compressed .RDS file containing the DESeq object.
                                  Annotation gene names: Ensembl and Gene Symbol."
                              )
                          )
                      ),
                      # fin dropdown
                  ),
                  #fin column
                  column(
                      width = 8,
                      align = "center",
                      actionBttn(
                          inputId = "enrich",
                          label = "Go DESeq DEVEA",
                          style = "jelly",
                          color = "primary",
                      )
                  ) #fin column
              ) # fin fluid row
          ),# fin miboxplus 3
          miBoxPlus(
            title = "Go Simple DEVEA",
            closable = FALSE,
            width = 4,
            background = "limon",
            height = "250px",
            tags$p("Enter a gene list (GL) without or with associated statistical values (GL + SV) 
            to run this analysis. Manage your threshold if GL + SV are provided to select your significant features
            or use your whole GL and get a full functional enrichment analysis 
            composed of tables and graphics based on Kyoto Encyclopedia of Genes and Genomes 
                   (", a("Kegg", href="https://www.genome.jp/kegg/pathway.html"), 
                   ") pathways database, Gene Ontology (", a("GO", href="http://geneontology.org"), ") 
                     terms resource and Gene Set Enrichment Analysis (", a("GSEA", href="https://www.gsea-msigdb.org/gsea/index.jsp"),")."),
            tags$br(),
            fluidRow(
              column(
                align = "center",
                width = 2,
                dropdown(
                  style = "unite",
                  icon = icon("info"),
                  status = "primary",
                  width = "600px",
                  animate = animateOptions(
                    enter = animations$fading_entrances$fadeInLeftBig,
                    exit = animations$fading_exits$fadeOutRightBig
                  ),
                  box(
                    width = 12,
                    status = "primary",
                    height = "150px",
                    solidHeader = TRUE,
                    title = tags$strong("Simple DEVEA 2022"),
                    background = "light-blue",
                    tags$h4(
                      "File formats allowed: .txt, .tsv or .xlsx files containing the gene list with statistical values associated or the simple gene list.
                      Annotation gene names: Ensembl and Gene Symbol."
                    )
                  )
                ),# fin dropdown
              ),#fin column
              column(
                width = 8,
                align = "center",
                actionBttn(
                  inputId = "simple",
                  label = "Go Simple DEVEA",
                  style = "jelly",
                  color = "primary",
                )
              ) #fin column
            ) # fin fluid row
          )) 
      ), # fin fluidrow 
      fluidRow(
      #tags$div(
        
        tags$a(href='https://jacob.cea.fr/drf/ifrancoisjacob/Pages/Departements/MIRCen/themes/astrocytes-reactifs-biomarqueurs-imagerie-cibles-therapeutiques.aspx', target="_blank",
               tags$img(src='mircen.png',width='60%',
                        style="padding: 10px; position: absolute; bottom:-95px; left:250px") ),
        tags$a(href='http://www.bioinformatica.imib.es', target="_blank",
               tags$img(src='IMIB_color_gris.svg',width='40%',
                        style="padding: 10px; position: absolute; bottom:-105px; right:300px;") ),
        tags$a(href='https://jacob.cea.fr/drf/ifrancoisjacob/Pages/Departements/CNRGH/LABORATOIRES/Bio-analyse.aspx', target="_blank",
               tags$img(src='cnrgh.png',width='50%',
                        style="padding: 10px; position: absolute; bottom:-110px; right:90px;") ),
        
        style = "position: relative; bottom:0; width:20%; margin-left: auto; margin-right: auto;"
      ) #fin div
  ) # fin dashboardbody
) # fin dashboarpage

# Plug router into Shiny server.
server <- function(input, output) {
  shinyjs::onclick("enrich", runjs("window.open('https://shiny.imib.es/DESeqDEVEA/','_blank')") )
  shinyjs::onclick("simple", runjs("window.open('https://shiny.imib.es/simpleDEVEA/','_blank')") )
  
  
  observeEvent(input$statButton, {
    showModal(popupModal())
  })
  
  #########################
  output$distPlot <- renderTmap({
    tm_shape(dataWorld) + tm_polygons("count") +
      tm_layout(legend.format = list(format="d") ) +
      tm_view(set.view = 0.5, view.legend.position = c("left","bottom") ) +
      tmap_options(basemaps = "OpenStreetMap")
  
  })
  #########################
  output$visits <- renderValueBox({
    valueBox(sum(dataWorld$count, na.rm = T), "Visits", width = 1, icon = icon("eye"), color="orange")
  })
  

 output$sankey <- renderSankeyNetwork({
 links <- data.frame(
      source=c("User DEseq Object","User DEseq Object","User DEseq Object","User DEseq Object",
               "ColData",
               "Expression Matrix",
               "User Gene List", "Gene List",
               "User Expression Matrix","User Expression Matrix",
               "Expression Matrix",
               "Pval + logFC",
               "User Gene List + Pval + logFC","User Gene List + Pval + logFC",
               "Gene List","Gene List","Gene List",
               "Expression Plots","Sample Info. Plots","Statistics Plots","Statistics Plots","Expression Plots","Expression Plots","Sample Info. Plots",
               "Expression Plots","Sample Info. Plots","Expression Plots","Sample Info. Plots","Expression Plots","Sample Info. Plots","Statistics Plots",
               "Enrichment Plots","Enrichment Plots","Enrichment Plots","Enrichment Plots","Enrichment Plots",
               "Enrichment Plots","Enrichment Plots","Enrichment Plots","Statistics Plots","Enrichment Plots","Statistics Plots",
               "Enrichment Plots","Statistics Plots"),
      target=c("Gene List","Pval + logFC","Expression Matrix","ColData",
               "Sample Info. Plots",
               "Expression Plots",
               "Gene List", "Enrichment Plots",
               "Gene List", "Expression Matrix",
               "Sample Info. Plots",
               "Statistics Plots",
               "Gene List","Pval + logFC",
               "Statistics Plots","Expression Plots","Sample Info. Plots",
               "PCA","PCA","Volcano","MA","MA","Boxplot/Violin","Boxplot/Violin",
               "Heatmap","Heatmap","Cluster","Cluster","TopGene","TopGene","Karyo plot",
               "Kegg Barplot","Kegg Chorplot","Kegg Dotplot","Kegg Heatmap","Kegg Netplot",
               "GO Bar","GO Dot","GO plotBar","GO plotBar","GO Circle","GO Circle",
               "GSEA","GSEA"),
      value=rep(3,44)
      )
    # From these flows we need to create a node data frame: it lists every entities involved in the flow
    nodes <- data.frame(
      name=c(as.character(links$source),
      as.character(links$target)) %>% unique()
    )
    # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
    links$IDsource <- match(links$source, nodes$name)-1
    links$IDtarget <- match(links$target, nodes$name)-1
    links$group <- "blue"
    #my_color <- 'd3.scaleOrdinal() .domain(["blue"]) .range(["blue"])'

    # Make the Network
    sankeyNetwork(Links = links, Nodes = nodes,
                  Source = "IDsource", Target = "IDtarget",
                  Value = "value", NodeID = "name",
                  sinksRight=TRUE, fontSize = 14, margin=0) #, colourScale = brewer.pal(9,'PuBu'))

  }) #Fin sankey
}

# Run server in a standard way.
shinyApp(ui, server)
