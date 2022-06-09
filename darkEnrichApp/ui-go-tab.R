fluidPage(
      tabsetPanel(
## All genes ##########################################
## BP all #######################
        tabPanel( "All DE genes",          # pestaña All
            #tags$br(),
            tabsetPanel(  # tabsetpanel all
                tabPanel("Biological Proccess",
                         fluidRow(# primera fila all
                             column(
                                 width = 9, offset = 2,
                                 box(
                                     title = "Biological proccess terms",
                                     solidHeader = FALSE,
                                     status = "primary",
                                     width = NULL,
                                     DTOutput("tableBPall"),
                                     actionButton("resettableBPall","Clear selection")
                                 ) # caja para la tabla
                             )), 
            fluidRow( # 2 fila
                column(
                  circleButton(
                    inputId = "informationGO1",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                  ),
                  bsTooltip(
                    "informationGO1",
                    paste0("Please, select rows with the paths of interest ",
                           "in the table of pathways above to activate all the plots. ",
                           "Some plots may need up to 4 paths to work properly. ",
                           "For further interpretation on the plots check ",
                           a("GOplot R package", href="https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html"),
                           " and ", a("clusterProfiler R package vignettes", 
                                      href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                    trigger = "hover",
                    placement = "left"
                  ),
                  width = 9, offset = 2,
                    tabBox( width = 12,  height = "650px",# caja con pestañas para los plots
                            tabPanel(title = "GO term cloud", height = "600px", width = "100%",
                                     fluidRow(
                                       column(width=12, downloadButton("cloudbpall","Download SVG"))),
                                     fluidRow(
                                       column(
                                         circleButton(
                                           inputId = "informationCloud1",
                                           icon = icon("info"),
                                           size = "xs",
                                           status = "primary"
                                         ),
                                         bsTooltip(
                                           "informationCloud1",
                                           paste0("The cloud will represent the observed terms with a min frequence of 2 from the selected level ",
                                                  "to the largest number from the Ont.level column in this ontology. ",
                                                  "The higher the level, the more semantically specific the term is. ",
                                                  "In an event of no possible cloud, the default one from level 0 is displayed."),
                                           trigger = "hover",
                                           placement = "left"
                                         ),
                                         width=2, 
                                              selectInput( "bpallLevel", label = "Select starting GO level", 
                                                           choices = list(0,1,2,3,4,5,6,7,8,9,10),
                                                           multiple = FALSE) ),
                                       column(width=10, plotOutput("cloudBPAll", height = "600px" )))
                                     ),
                            tabPanel(title = "Barplot",
                                     fluidRow(column(
                                         width = 3,
                                         radioGroupButtons(
                                             inputId = "selectbpall",
                                             label = "Select bar plot type",
                                             choices = c("Dodge", "Stack", "Opposite"),
                                             selected = "Dodge",
                                             size = "sm",
                                             status = "primary",
                                             checkIcon = list(
                                                 yes = icon("ok",
                                                            lib = "glyphicon"),
                                                 no = icon("remove",
                                                           lib = "glyphicon")
                                             )
                                         )
                                     ),
                                     column(width=2, offset = 4,
                                                     downloadButton("barBpAll","Download SVG"))
                                     ), # fin fluidRow, column & radioGroupButtons
                                     fluidRow(class = "text-center",
                                              column(
                                                  align = "center",
                                                  plotlyOutput("plotBPall"),
                                                  width = 9
                                              ))),  #barplot
                            tabPanel(title = "Dotplot",
                                     fluidRow(column(width=2,
                                            downloadButton("dotBpAll","Download SVG"))),
                                     plotlyOutput("BPDotall")
                                     ), # dotplot
                            tabPanel(title = "GoBarplot",
                                     fluidRow(column(width=2,
                                            downloadButton("gobarBpAll","Download SVG"))),
                                     plotOutput("gobarplotAllBP")
                                     ),
                            tabPanel(title = "GoCirclePlot",
                                     fluidRow(column(width=2,
                                            downloadButton("cirBpAll","Download SVG"))),
                                     plotOutput("goCircleAllBP")
                                     )
                        )
                    )
            )),
                ### MF all #############################
            #tags$br(),
            tabPanel("Molecular function",  # MF all
             fluidRow(  # primera fila mf all 
                column( width = 9, offset = 2,
                    box(title = "Molecular function terms",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("tableMFall"),
                        actionButton("resettableMFall","Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column(
                  circleButton(
                    inputId = "informationGO2",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                  ),
                  bsTooltip(
                    "informationGO2",
                    paste0("Please, select rows with the paths of interest ",
                           "in the table of pathways above to activate all the plots. ",
                           "Some plots may need up to 4 paths to work properly. ",
                           "For further interpretation on the plots check ",
                           a("GOplot R package", href="https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html"),
                           " and ", a("clusterProfiler R package vignettes", 
                                      href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                    trigger = "hover",
                    placement = "left"
                  ),
                    width = 9, offset = 2,
                    tabBox( width = 12,  height = "650px",# caja con pestañas para los plots
                            tabPanel(title = "GO term cloud",width = "100%", height = "600px",
                                     fluidRow(
                                       column(width = 12,downloadButton("cloudmfall","Download SVG"))),
                                     fluidRow(
                                       column(
                                         circleButton(
                                         inputId = "informationCloud2",
                                         icon = icon("info"),
                                         size = "xs",
                                         status = "primary"
                                       ),
                                       bsTooltip(
                                         "informationCloud2",
                                         paste0("The cloud will represent the observed terms with a min frequence of 2 from the selected level ",
                                                "to the largest number from the Ont.level column in this ontology. ",
                                                "The higher the level, the more semantically specific the term is. ",
                                                "In an event of no possible cloud, the default one from level 0 is displayed."),
                                         trigger = "hover",
                                         placement = "left"
                                       ),
                                       width=2, 
                                              selectInput( "mfallLevel", label = "Select starting GO level", 
                                                           choices = list(0,1,2,3,4,5,6,7,8,9,10),
                                                           multiple = FALSE) ),
                                       column(width = 10,
                                         plotOutput("cloudMFAll", height = "600px")))), 
                            tabPanel(title = "Barplot",
                                     fluidRow(column(
                                         width = 3,
                                         radioGroupButtons(
                                             inputId = "selectmfall",
                                             label = "Select bar plot type",
                                             choices = c("Dodge", "Stack", "Opposite"),
                                             selected = "Dodge",
                                             size = "sm",
                                             status = "primary",
                                             checkIcon = list(
                                                 yes = icon("ok",
                                                            lib = "glyphicon"),
                                                 no = icon("remove",
                                                           lib = "glyphicon")
                                             )
                                         )
                                     ),
                                     column(width=2, offset = 3,
                                                     downloadButton("barMfAll","Download SVG"))
                                     ), # fin fluidRow, column & radioGroupButtons
                                     fluidRow(class = "text-center",
                                              column(
                                                  align = "center",
                                                  plotlyOutput("plotMFall"),
                                                  width = 9
                                              ))),  #barplot
                            tabPanel(title = "Dotplot",
                                     fluidRow(column(width=2,
                                            downloadButton("dotMfAll","Download SVG"))),
                                     plotlyOutput("MFDotall")
                                     ), # dotplot
                            tabPanel(title = "GoBarplot",
                                     fluidRow(column(width=2,
                                            downloadButton("gobarMfAll","Download SVG"))),
                                     plotOutput("gobarplotAllMF")
                                     ),
                            tabPanel(title = "GoCirclePlot",
                                     fluidRow(column(width=2,
                                            downloadButton("cirMfAll","Download SVG"))),
                                     plotOutput("goCircleAllMF"))
                        )
                    )
            )),
                ### CC all #############################
                tabPanel("Cellular component",   ##tab CC all
                fluidRow(  # primera fila all CC
                column( width = 9, offset = 2,
                    box(title = "Cellular component terms",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("tableCCall"),
                        actionButton("resettableCCall","Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column(
                  circleButton(
                    inputId = "informationGO3",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                  ),
                  bsTooltip(
                    "informationGO3",
                    paste0("Please, select rows with the paths of interest ",
                           "in the table of pathways above to activate all the plots. ",
                           "Some plots may need up to 4 paths to work properly. ",
                           "For further interpretation on the plots check ",
                           a("GOplot R package", href="https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html"),
                           " and ", a("clusterProfiler R package vignettes", 
                                      href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                    trigger = "hover",
                    placement = "left"
                  ),
                    width = 9, offset = 2,
                    tabBox( width = 12,  height = "650px",# caja con pestañas para los plots
                            tabPanel(title = "GO term cloud", width = "100%", height = "600px",
                                    fluidRow(downloadButton("cloudccall","Download SVG")),
                                    fluidRow(
                                      column(
                                        circleButton(
                                          inputId = "informationCloud3",
                                          icon = icon("info"),
                                          size = "xs",
                                          status = "primary"
                                        ),
                                        bsTooltip(
                                          "informationCloud3",
                                          paste0("The cloud will represent the observed terms with a min frequence of 2 from the selected level ",
                                                 "to the largest number from the Ont.level column in this ontology. ",
                                                 "The higher the level, the more semantically specific the term is. ",
                                                 "In an event of no possible cloud, the default one from level 0 is displayed."),
                                          trigger = "hover",
                                          placement = "left"
                                        ),
                                        width=2, 
                                             selectInput( "ccallLevel", label = "Select starting GO level", 
                                                          choices = list(0,1,2,3,4,5,6,7,8,9,10),
                                                          multiple = FALSE) ),
                                      column(width=10,plotOutput("cloudCCAll", height = "600px") ))
                                    ),
                            tabPanel(title = "Barplot",
                                     fluidRow(column(
                                         width = 3,
                                         radioGroupButtons(
                                             inputId = "selectccall",
                                             label = "Select bar plot type",
                                             choices = c("Dodge", "Stack", "Opposite"),
                                             selected = "Dodge",
                                             size = "sm",
                                             status = "primary",
                                             checkIcon = list(
                                                 yes = icon("ok",
                                                            lib = "glyphicon"),
                                                 no = icon("remove",
                                                           lib = "glyphicon")
                                             )
                                         )
                                     ),
                                     column(width=2, offset = 3,
                                                     downloadButton("barCcAll","Download SVG"))
                                     ), # fin fluidRow, column & radioGroupButtons
                                     fluidRow(class = "text-center",
                                              column(
                                                  align = "center",
                                                  plotlyOutput("plotCCall"),
                                                  width = 9
                                              ))),  #barplot
                            tabPanel(title = "Dotplot",
                                     fluidRow(column(width=2,
                                            downloadButton("dotCcAll","Download SVG"))),
                                     plotlyOutput("CCDotall")
                                     ), # dotplot
                            tabPanel(title = "GoBarplot",
                                     fluidRow(column(width=2,
                                            downloadButton("gobarCcAll","Download SVG"))),
                                     plotOutput("gobarplotAllCC")
                                     ),
                            tabPanel(title = "GoCirclePlot",
                                     fluidRow(column(width=2,
                                            downloadButton("cirCcAll","Download SVG"))),
                                     plotOutput("goCircleAllCC")
                            )
                           
                        ) #tabbox
                    ) #column 
            ) #fluidrow
            ) # tabpanel CC 
            ) # tabsetpanel all
        ), #fin tab all genes
###### Upregulated genes #####################################
## BP up #############################
        tabPanel( "Upregulated genes",          # pestaña upregulates
            tabsetPanel(  # tabsetpanel up
             tabPanel("Biological Proccess",
            fluidRow(  # primera fila
              column( width = 9, offset = 2,
                    box(title = "Biological proccess terms",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("tableBP"),
                        actionButton("resettableBP","Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column(
                  circleButton(
                    inputId = "informationGO4",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                  ),
                  bsTooltip(
                    "informationGO4",
                    paste0("Please, select rows with the paths of interest ",
                           "in the table of pathways above to activate all the plots. ",
                           "Some plots may need up to 4 paths to work properly. ",
                           "For further interpretation on the plots check ",
                           a("GOplot R package", href="https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html"),
                           " and ", a("clusterProfiler R package vignettes", 
                                      href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                    trigger = "hover",
                    placement = "left"
                  ),
                    width = 9, offset = 2,
                    tabBox( width = 12, height = "650px", # caja con pestañas para los plots
                        tabPanel(title = "GO term cloud",  height = "600px",
                                 fluidRow(downloadButton("cloudbpup","Download SVG")),
                                 fluidRow(
                                   column(
                                     circleButton(
                                       inputId = "informationCloud4",
                                       icon = icon("info"),
                                       size = "xs",
                                       status = "primary"
                                     ),
                                     bsTooltip(
                                       "informationCloud4",
                                       paste0("The cloud will represent the observed terms with a min frequence of 2 from the selected level ",
                                              "to the largest number from the Ont.level column in this ontology. ",
                                              "The higher the level, the more semantically specific the term is. ",
                                              "In an event of no possible cloud, the default one from level 0 is displayed."),
                                       trigger = "hover",
                                       placement = "left"
                                     ),
                                     width=2, 
                                          selectInput( "bpupLevel", label = "Select starting GO level", 
                                                       choices = list(0,1,2,3,4,5,6,7,8,9,10),
                                                       multiple = FALSE) ),
                                   column(width=10,plotOutput("cloudBPUp", height = "600px") ))
                                 ),
                        tabPanel(title = "Barplot",
                                 fluidRow(column(width=2,
                                            downloadButton("barBpUp","Download SVG"))),
                                  plotlyOutput("plotBP")
                                              ),  #barplot
                        tabPanel(title = "Dotplot",
                                 fluidRow(column(width=2,
                                            downloadButton("dotBpUp","Download SVG"))),
                                 plotlyOutput("BPDotUp")
                                 ), # dotplot
                        tabPanel(title = "GoBarplot",
                                 fluidRow(column(width=2,
                                            downloadButton("gobarBpUp","Download SVG"))),
                                 plotOutput("gobarplotUpBP")
                                     ),
                        tabPanel(title = "GoCirclePlot",
                                 fluidRow(column(width=2,
                                            downloadButton("cirBpUp","Download SVG"))),
                                 plotOutput("goCircleUpBP")
                                     )
                        )
                    )
            )),
## MF up #############################
tabPanel("Molecular function",
fluidRow(  # primera fila
                column( width = 9, offset = 2,
                    box(title = "Molecular function terms",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("tableMF"),
                        actionButton("resettableMF","Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column( 
                  circleButton(
                    inputId = "informationGO5",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                  ),
                  bsTooltip(
                    "informationGO5",
                    paste0("Please, select rows with the paths of interest ",
                           "in the table of pathways above to activate all the plots. ",
                           "Some plots may need up to 4 paths to work properly. ",
                           "For further interpretation on the plots check ",
                           a("GOplot R package", href="https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html"),
                           " and ", a("clusterProfiler R package vignettes", 
                                      href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                    trigger = "hover",
                    placement = "left"
                  ),
                    width = 9, offset = 2,
                    tabBox( width = 12,  height = "650px",# caja con pestañas para los plots
                        tabPanel(title = "GO term cloud",  height = "600px", width = "100%",
                                 fluidRow(downloadButton("cloudmfup","Download SVG")),
                                 fluidRow( 
                                   column(
                                     circleButton(
                                       inputId = "informationCloud5",
                                       icon = icon("info"),
                                       size = "xs",
                                       status = "primary"
                                     ),
                                     bsTooltip(
                                       "informationCloud5",
                                       paste0("The cloud will represent the observed terms with a min frequence of 2 from the selected level ",
                                              "to the largest number from the Ont.level column in this ontology. ",
                                              "The higher the level, the more semantically specific the term is. ",
                                              "In an event of no possible cloud, the default one from level 0 is displayed."),
                                       trigger = "hover",
                                       placement = "left"
                                     ),
                                     width=2, 
                                          selectInput( "mfupLevel", label = "Select starting GO level", 
                                                       choices = list(0,1,2,3,4,5,6,7,8,9,10),
                                                       multiple = FALSE) ),
                                   column(width=10, plotOutput("cloudMFUp", height = "600px") ))
                                 ),
                        tabPanel(title = "Barplot",
                                 fluidRow(column(width=2,
                                            downloadButton("barMfUp","Download SVG"))),
                                  plotlyOutput("plotMF")
                                              ),  #barplot
                        tabPanel(title = "Dotplot",
                                 fluidRow(column(width=2,
                                            downloadButton("dotMfUp","Download SVG"))),
                                 plotlyOutput("MFDotUp")
                                 ), # dotplot
                        tabPanel(title = "GoBarplot",
                                 fluidRow(column(width=2,
                                            downloadButton("gobarMfUp","Download SVG"))),
                                 plotOutput("gobarplotUpMF")
                                     ),
                        tabPanel(title = "GoCirclePlot",
                                 fluidRow(column(width=2,
                                            downloadButton("cirMfUp","Download SVG"))),
                                 plotOutput("goCircleUpMF")
                                     )
                        )
                    )
            )),
## CC up ###################################
        tabPanel("Cellular component",
            fluidRow(  # primera fila
                        column( width = 9, offset = 2,
                            box(title = "Cellular component terms",
                                solidHeader = FALSE,
                                status = "primary",
                                width = NULL,
                                DTOutput("tableCC"),
                                actionButton("resettableCC","Clear selection")
                            ) # caja para la tabla
                            )
                        ),
            fluidRow( # 2 fila
                column(
                  circleButton(
                    inputId = "informationGO6",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                  ),
                  bsTooltip(
                    "informationGO6",
                    paste0("Please, select rows with the paths of interest ",
                           "in the table of pathways above to activate all the plots. ",
                           "Some plots may need up to 4 paths to work properly. ",
                           "For further interpretation on the plots check ",
                           a("GOplot R package", href="https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html"),
                           " and ", a("clusterProfiler R package vignettes", 
                                      href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                    trigger = "hover",
                    placement = "left"
                  ),
                    width = 9, offset = 2,
                    tabBox( width = 12, height = "650px", # caja con pestañas para los plots
                        tabPanel(title = "GO term cloud",  height = "600px", width = "100%" ,
                                 fluidRow(downloadButton("cloudccup","Download SVG")),
                                 fluidRow(
                                   column(
                                     circleButton(
                                       inputId = "informationCloud6",
                                       icon = icon("info"),
                                       size = "xs",
                                       status = "primary"
                                     ),
                                     bsTooltip(
                                       "informationCloud6",
                                       paste0("The cloud will represent the observed terms with a min frequence of 2 from the selected level ",
                                              "to the largest number from the Ont.level column in this ontology. ",
                                              "The higher the level, the more semantically specific the term is. ",
                                              "In an event of no possible cloud, the default one from level 0 is displayed."),
                                       trigger = "hover",
                                       placement = "left"
                                     ),
                                     width=2, 
                                          selectInput( "ccupLevel", label = "Select starting GO level", 
                                                       choices = list(0,1,2,3,4,5,6,7,8,9,10),
                                                       multiple = FALSE) ),
                                   column(width=10, plotOutput("cloudCCUp", height = "600px") ))
                                 ),
                        tabPanel(title = "Barplot",
                                 fluidRow(column(width=2,
                                            downloadButton("barCcUp","Download SVG"))),
                                  plotlyOutput("plotCC")
                                              ),  #barplot
                        tabPanel(title = "Dotplot",
                                 fluidRow(column(width=2,
                                            downloadButton("dotCcUp","Download SVG"))),
                                 plotlyOutput("CCDotUp")
                                 ), # dotplot
                        tabPanel(title = "GoBarplot",
                                 fluidRow(column(width=2,
                                            downloadButton("gobarCcUp","Download SVG"))),
                                 plotOutput("gobarplotUpCC")
                                     ),
                        tabPanel(title = "GoCirclePlot",
                                 fluidRow(column(width=2,
                                            downloadButton("cirCcUp","Download SVG"))),
                                 plotOutput("goCircleUpCC")
                                     )
                            )
                        )
                    )
        ))), #fin tab upregulated genes
## Down regulated genes ############################################
## BP down ##############################
        tabPanel( "Downregulated genes",          # pestaña downregulates
            tabsetPanel(
                tabPanel( "Biological proccess",
            fluidRow(  # primera fila
                column( width = 9, offset = 2,
                    box(title = "Biological proccess terms",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("tableBPdown"),
                        actionButton("resettableBPdown","Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column( 
                  circleButton(
                    inputId = "informationGO7",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                  ),
                  bsTooltip(
                    "informationGO7",
                    paste0("Please, select rows with the paths of interest ",
                           "in the table of pathways above to activate all the plots. ",
                           "Some plots may need up to 4 paths to work properly. ",
                           "For further interpretation on the plots check ",
                           a("GOplot R package", href="https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html"),
                           " and ", a("clusterProfiler R package vignettes", 
                                      href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                    trigger = "hover",
                    placement = "left"
                  ),
                    width = 9, offset = 2,
                    tabBox( width = 12,  height = "650px",# caja con pestañas para los plots
                        tabPanel(title = "GO term cloud",  height = "600px", width = "100%",
                                 fluidRow(downloadButton("cloudbpdown","Download SVG")),
                                 fluidRow(
                                   column(
                                     circleButton(
                                       inputId = "informationCloud7",
                                       icon = icon("info"),
                                       size = "xs",
                                       status = "primary"
                                     ),
                                     bsTooltip(
                                       "informationCloud7",
                                       paste0("The cloud will represent the observed terms with a min frequence of 2 from the selected level ",
                                              "to the largest number from the Ont.level column in this ontology. ",
                                              "The higher the level, the more semantically specific the term is. ",
                                              "In an event of no possible cloud, the default one from level 0 is displayed."),
                                       trigger = "hover",
                                       placement = "left"
                                     ),
                                     width=2, 
                                          selectInput( "bpdownLevel", label = "Select starting GO level", 
                                                       choices = list(0,1,2,3,4,5,6,7,8,9,10),
                                                       multiple = FALSE) ),
                                   column(width=10, plotOutput("cloudBPDown", height = "600px") ))
                                 ),
                        tabPanel(title = "Barplot",
                                 fluidRow(column(width=2,
                                            downloadButton("barBpDown","Download SVG"))),
                                  plotlyOutput("plotBPdown")
                                              ),  #barplot
                        tabPanel(title = "Dotplot",
                                 fluidRow(column(width=2,
                                            downloadButton("dotBpDown","Download SVG"))),
                                 plotlyOutput("BPDotDown")
                                 ), # dotplot
                        tabPanel(title = "GoBarplot",
                                 fluidRow(column(width=2,
                                            downloadButton("gobarBpDown","Download SVG"))),
                                 plotOutput("gobarplotDownBP")
                                     ),
                        tabPanel(title = "GoCirclePlot",
                                 fluidRow(column(width=2,
                                            downloadButton("cirBpDown","Download SVG"))),
                                 plotOutput("goCircleDownBP")
                                     )
                        )
                    )
            )),
## MF down ###################################
        tabPanel( "Molecular function",
            fluidRow(  # primera fila
                column( width = 9, offset = 2,
                    box(title = "Molecular functions terms",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("tableMFdown"),
                        actionButton("resettableMFdown","Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column(
                  circleButton(
                    inputId = "informationGO8",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                  ),
                  bsTooltip(
                    "informationGO8",
                    paste0("Please, select rows with the paths of interest ",
                           "in the table of pathways above to activate all the plots. ",
                           "Some plots may need up to 4 paths to work properly. ",
                           "For further interpretation on the plots check ",
                           a("GOplot R package", href="https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html"),
                           " and ", a("clusterProfiler R package vignettes", 
                                      href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                    trigger = "hover",
                    placement = "left"
                  ),
                  width = 9, offset = 2,
                    tabBox( width = 12, height = "650px", # caja con pestañas para los plots
                        tabPanel(title = "GO term cloud",  height = "600px",width = "100%",
                                 fluidRow(downloadButton("cloudmfdown","Download SVG")),
                                 fluidRow( 
                                   column(
                                     circleButton(
                                     inputId = "informationCloud8",
                                     icon = icon("info"),
                                     size = "xs",
                                     status = "primary"
                                   ),
                                   bsTooltip(
                                     "informationCloud8",
                                     paste0("The cloud will represent the observed terms with a min frequence of 2 from the selected level ",
                                            "to the largest number from the Ont.level column in this ontology. ",
                                            "The higher the level, the more semantically specific the term is. ",
                                            "In an event of no possible cloud, the default one from level 0 is displayed."),
                                     trigger = "hover",
                                     placement = "left"
                                   ),
                                   width=2, 
                                          selectInput( "mfdownLevel", label = "Select starting GO level", 
                                                       choices = list(0,1,2,3,4,5,6,7,8,9,10),
                                                       multiple = FALSE) ),
                                   column(width=10, plotOutput("cloudMFDown", height = "600px") ))
                                 ),
                        tabPanel(title = "Barplot",
                                 fluidRow(column(width=2,
                                            downloadButton("barMfDown","Download SVG"))),
                                  plotlyOutput("plotMFdown")
                                              ),  #barplot
                        tabPanel(title = "Dotplot",
                                  fluidRow(column(width=2,
                                            downloadButton("dotMfDown","Download SVG"))),
                                 plotlyOutput("MFDotDown")
                                 ), # dotplot
                        tabPanel(title = "GoBarplot",
                                 fluidRow(column(width=2,
                                            downloadButton("gobarMfDown","Download SVG"))),
                                 plotOutput("gobarplotDownMF")
                                     ),
                        tabPanel(title = "GoCirclePlot",
                                 fluidRow(column(width=2,
                                            downloadButton("cirMfDown","Download SVG"))),
                                 plotOutput("goCircleDownMF")
                                     )
                        )
                    )
            )),
## CC down ###################################
                tabPanel("Cellular component",
                fluidRow(  # primera fila
                  column( width = 9, offset = 2,
                    box(title = "Cellular component terms",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("tableCCdown"),
                        actionButton("resettableCCdown","Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column( 
                  circleButton(
                    inputId = "informationGO9",
                    icon = icon("info"),
                    size = "xs",
                    status = "primary"
                  ),
                  bsTooltip(
                    "informationGO9",
                    paste0("Please, select rows with the paths of interest ",
                           "in the table of pathways above to activate all the plots. ",
                           "Some plots may need up to 4 paths to work properly. ",
                           "For further interpretation on the plots check ",
                           a("GOplot R package", href="https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html"),
                           " and ", a("clusterProfiler R package vignettes", 
                                      href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                    trigger = "hover",
                    placement = "left"
                  ),
                  width = 9, offset = 2,
                    tabBox( width = 12, height = "650px",# caja con pestañas para los plots
                        tabPanel(title = "GO term cloud",  height = "600px",width = "100%",
                                 fluidRow( downloadButton("cloudccdown","Download SVG")),
                                 fluidRow( 
                                   column(
                                     circleButton(
                                       inputId = "informationCloud9",
                                       icon = icon("info"),
                                       size = "xs",
                                       status = "primary"
                                     ),
                                     bsTooltip(
                                       "informationCloud9",
                                       paste0("The cloud will represent the observed terms with a min frequence of 2 from the selected level ",
                                              "to the largest number from the Ont.level column in this ontology. ",
                                              "The higher the level, the more semantically specific the term is. ",
                                              "In an event of no possible cloud, the default one from level 0 is displayed."),
                                       trigger = "hover",
                                       placement = "left"
                                     ),
                                     width=2, 
                                          selectInput( "ccdownLevel", label = "Select starting GO level", 
                                                       choices = list(0,1,2,3,4,5,6,7,8,9,10),
                                                       multiple = FALSE) ),
                                   column(width=10, plotOutput("cloudCCDown", height = "600px") ))
                                 ),
                        tabPanel(title = "Barplot",
                                 fluidRow(column(width=2,
                                            downloadButton("barCcDown","Download SVG"))),
                                  plotlyOutput("plotCCdown")
                                              ),  #barplot
                        tabPanel(title = "Dotplot",
                                 fluidRow(column(width=2,
                                            downloadButton("dotCcDown","Download SVG"))),
                                 plotlyOutput("CCDotDown")
                                 ), # dotplot
                        tabPanel(title = "GoBarplot",
                                  fluidRow(column(width=2,
                                            downloadButton("gobarCcDown","Download SVG"))),
                                 plotOutput("gobarplotDownCC")
                                     ),
                        tabPanel(title = "GoCirclePlot",
                                 fluidRow(column(width=2,
                                            downloadButton("cirCcDown","Download SVG"))),
                                 plotOutput("goCircleDownCC")
                                     )
                        )
                    )
            ))
) # fin tabsetpanel downregulated
        ) #fin tab downrergulated genes
    ) # fin tabsetpanel
)#fin fluidpage