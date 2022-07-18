fluidPage(
    tabsetPanel(
        tabPanel( "All DE genes",          # pestaña All #################
            tags$br(),
            fluidRow(  # primera fila
                column( width = 9, offset = 2,
                    box(title = "Table of pathways",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("tableAll"),
                        actionButton("resettableall", "Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column( 
                    circleButton(
                        inputId = "information5",
                        icon = icon("info"),
                        size = "xs",
                        status = "primary"
                    ),
                    bsTooltip(
                        "information5",
                        paste0("Please, select rows with the paths of interest ",
                               "in the table of pathways above to activate all the plots. ",
                               "Some plots may need up to 4 paths to work properly. ",
                               "For further interpretation on the plots check ",
                               a("clusterProfiler R package vignette", 
                               href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                        trigger = "hover",
                        placement = "left"
                    ),
                    width = 9, offset = 2,
                    tabBox( width = 12, # caja con pestañas para los plots
                            tabPanel(title = "Barplot",
                                    fluidRow(column(
                                         width = 3,
                                         radioGroupButtons(
                                             inputId = "selectkeggall",
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
                                     column(width=2, offset = 3, downloadButton("barKeggAll","Download SVG")),
                                     ), # fin fluidRow, column & radioGroupButtons
                                     fluidRow(class = "text-center",
                                              column(
                                                  align = "center",
                                                  plotlyOutput("keggPlotAll"),
                                                  width = 9
                                              ))),  #barplot
                        tabPanel(title = "Chordplot",
                                 tagList(fluidRow(column(width = 8,
                                                         mychordplotOutput("keggChordAll",
                                                                         width = "100%",
                                                                         height = "600px") 
                                                         ),
                                                  column(width = 4,
                                                         plotOutput("legendChorAll", width="100%")
                                                         )) )
                                 ), #cordplot
                        tabPanel(title = "Dotplot",
                                 fluidRow(column(width=2,
                                                 downloadButton("dotkeggAll","Download SVG"))),
                                 plotlyOutput("keggDotAll")
                                 ), # dotplot
                        tabPanel(title = "Heatmap",
                                 fluidRow(column(
                                     width=2,
                                     downloadButton(outputId = "heatKeggAll","Download SVG")
                                 )),
                                 plotlyOutput("heatmapKeggAll", height = "600px")
                                 ), # heatmap
                        tabPanel(title = "Netplot",
                                 column(width = 1,
                                        switchInput(
                                            size = "mini",
                                            inputId = "keggAllNet_switch",
                                            offLabel = "Static",
                                            onLabel = "Interactive"),
                                        plotOutput("legend")
                                        ),
                                 column(width = 11,
                                        uiOutput("keggAllNet")
                                        ),
                                 fluidRow(column(
                                     width=2,
                                     downloadButton(outputId = "cnetKggAll","Download SVG")
                                 ))
                                 #plotOutput("cnetKeggAll")
                        #          ), # cnetplot
                        # tabPanel(title = "VisNetPlot",
                        #          visNetworkOutput("visnetKeggAll", height = "600px")
                                  ) #visnetall
                        )
                    )
            )
        ), #fin tab all genes ..................##############
        tabPanel( "Upregulated genes",          # pestaña Up ##############
            tags$br(),
            fluidRow(  # primera fila
                column( width = 9, offset = 2,
                    box(title = "Table of pathways",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("table"),
                        actionButton("resettable", "Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column( 
                    circleButton(
                        inputId = "information6",
                        icon = icon("info"),
                        size = "xs",
                        status = "primary"
                    ),
                    bsTooltip(
                        "information6",
                        paste0("Please, select rows with the paths of interest ",
                               "in the table of pathways above to activate all the plots. ",
                               "Some plots may need up to 4 paths to work properly. ",
                               "For further interpretation on the plots check ",
                               a("clusterProfiler R package vignette", 
                                 href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                        trigger = "hover",
                        placement = "left"
                    ),
                    width = 9, offset = 2,
                    tabBox( width = 12, # caja con pestañas para los plots
                        tabPanel(title = "Barplot",
                                 fluidRow(column(width=3,
                                                 downloadButton("barKeggUp","Download SVG"))),
                                 plotlyOutput("keggPlot")
                                              ),  #barplot
                        tabPanel(title = "Chordplot",
                                 tagList(fluidRow(
                                     column(
                                         width = 8,
                                         mychordplotOutput("keggChord",
                                                         width = "100%",
                                                         height = "600px")
                                     ),
                                     column(
                                         width = 4,
                                         plotOutput("legendChorUp", 
                                                    width ="100%",
                                                    height ="600px")
                                     )
                                 ))
                                 ), #cordplot
                        tabPanel(title = "Dotplot",
                                 fluidRow(column(width=2,
                                         downloadButton("dotKeggUp","Download SVG"))),
                                 plotlyOutput("keggDotUp")
                                 ), # dotplot
                        tabPanel(title = "Heatmap",
                                 fluidRow(column(width=2,
                                         downloadButton("heatKeggUp","Download SVG"))),
                                 plotlyOutput("heatmapKeggUp", height = "600px")
                                 ), # heatmap
                        tabPanel(title = "Netplot",
                                 column(width = 1,
                                        switchInput(
                                            size = "mini",
                                            inputId = "keggUpNet_switch",
                                            offLabel = "Static",
                                            onLabel = "Interactive")
                                        ),
                                 column(width = 11,
                                        uiOutput("keggUpNet"),
                                        fluidRow(column(width=2,
                                            downloadButton("cnetkeggUp","Download SVG"))),
                                        )
                                  ) #visnetup
                        )
                    )
            )
        ), #fin tab Up genes ................... #####
        tabPanel( "Downregulated genes",   # pestaña Down ##########
            tags$br(),
            fluidRow(  # primera fila
                column( width = 9, offset = 2,
                    box(title = "Table of pathways",
                        solidHeader = FALSE,
                        status = "primary",
                        width = NULL,
                        DTOutput("tableDown"),
                        actionButton("resettableDown", "Clear selection")
                    ) # caja para la tabla
                    )
                ),
            fluidRow( # 2 fila
                column( 
                    circleButton(
                        inputId = "information7",
                        icon = icon("info"),
                        size = "xs",
                        status = "primary"
                    ),
                    bsTooltip(
                        "information7",
                        paste0("Please, select rows with the paths of interest ",
                                "in the table of pathways above to activate all the plots. ",
                               "Some plots may need up to 4 paths to work properly. ",
                               "For further interpretation on the plots check ",
                               a("clusterProfiler R package vignette", 
                                 href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")),
                        trigger = "hover",
                        placement = "left"
                    ),
                    width = 9, offset = 2,
                    tabBox( width = 12, # caja con pestañas para los plots
                        tabPanel(title = "Barplot",
                                 fluidRow(column(width=3,
                                    downloadButton("barKeggDown","Download SVG"))),
                                 plotlyOutput("keggPlotDown")
                                              ),  #barplot
                        tabPanel(title = "Chordplot",
                                 tagList(fluidRow(
                                     column(
                                         width = 8,
                                         mychordplotOutput("keggChordDown",
                                                         width = "100%",
                                                         height = "600px")
                                     ),
                                     column(
                                         width = 4,
                                         plotOutput("legendChorDown",
                                                    width ="100%",
                                                    height ="600px")
                                     )
                                     ))
                                 ), #cordplot
                        tabPanel(title = "Dotplot",
                                 fluidRow(column(width=2,
                                            downloadButton("dotKeggDown","Download SVG"))),
                                 plotlyOutput("keggDotDown")
                                 ), # dotplot
                        tabPanel(title = "Heatmap",
                                 fluidRow(column(width=2,
                                            downloadButton("heatKeggDown","Download SVG"))),
                                 plotlyOutput("heatmapKeggDown", height = "600px")
                                 ), # heatmap
                        tabPanel(title = "Netplot",
                                 column(width = 1,
                                        switchInput(
                                            size = "mini",
                                            inputId = "keggDownNet_switch",
                                            offLabel = "Static",
                                            onLabel = "Interactive")
                                        ),
                                 column(width = 11,
                                        uiOutput("keggDownNet"),
                                        fluidRow(column(width=2,
                                            downloadButton("cnetkeggDown","Download SVG"))),
                                        )
                                  ) #cnetplot
                        )
                    )
            )
        ) #fin tab Down genes ............. ####
    ) # fin tabsetpanel
) #fin fluidpage    


