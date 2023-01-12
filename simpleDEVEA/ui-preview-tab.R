fluidPage(
    fluidRow(infoBoxOutput("allbox", width = 4),
             infoBoxOutput("downbox", width = 4),
             infoBoxOutput("upbox", width = 4)),
    fluidRow(column(width = 2, offset = 3,
                 uiOutput("fcdown")),
             column(width=2,
                 uiOutput("fcup")),
             column(width = 2,
                 uiOutput("pval"))
             ),
    # fin fluidrow boxinfos
    tags$br(),
    fluidRow(
        column(
            width = 5,
            offset = 0,
          circleButton(
                inputId = "information1",
                icon = icon("info"),
                size = "xs",
                status = "primary"
            ),
            bsTooltip(
                    "information1",
                    paste0("Customize here the statistical cutoff to consider a ",
                           "gene as differentially expressed. The current threshold can be checked above. ",
                           "Click to apply values after every modification to update."),
                    trigger = "hover",
                    placement = "right"
                ),
        box(
            align = 'center',
            title = "Cutoff values",
            solidHeader = FALSE,
            status = "primary",
            width = 12,
            switchInput(inputId = "fc_switch",
                        offLabel = "log2FC",
                        onLabel = "FC",
                        size = "mini"
                        ),
            tags$br(),
            uiOutput("fc_control"),
            tags$br(),
            uiOutput("padj"),
            tags$br(),
            actionButton("applyParam", label = "Click to apply values")
        ),
        fluidRow(
          
          column(
            width = 10, offset = 3,
            circleButton(
              inputId = "informaBack",
              icon = icon("info"),
              size = "xs",
              status = "primary"
            ),
            bsTooltip(
              "informaBack",
              paste0("By default, genes used as enrichment ",
                     "background consist of the whole stated species genome. ", 
                     "Add your own universe list as a file with your favourite genes. ",
                     "Click Run enrichment to unlock the enrichment analysis tabs."),
              
              trigger = "hover",
              placement = "right"
            ),
            fileInput("universef", label="Backgroung dataset", placeholder = "Leave empty to use entire database"),  
            
            fluidRow(column(width = 12, offset = 0,
                            strong("Click to compute enrichment"),
                            tags$br(),
                            actionBttn("enrichButton", label = "Run enrichment", 
                                       size = "lg",color = "default",icon = icon("images") )
            ))
            
          )
          
          ) # end fluidrow
        ), #end column cutoff p-values controls

    column(
      width = 6,
      box(width=12,
        DT::dataTableOutput("tablepreview")
      )
    )
    ), 
  fluidRow(
      column( width = 3,
              circleButton(
                inputId = "infoColor",
                icon = icon("info"),
                size = "xs",
                status = "primary"
              ),
              bsTooltip(
                "infoColor",
                paste0("Customize color applied for ",
                       "up and/or down-regulated genes."),
                
                trigger = "hover",
                placement = "right"
              ),
              box( title = "Customize plots",
                  status = "info",
                  width = NULL,
                tagList(
                    tags$p("Volcano plots color scheme"),
             spectrumInput(
                 inputId = "upColor",
                 label = "Pick upregulated color:",
                 selected = "#b30000",
                 width = "60%",
                 choices = list(
                     list("#b30000","#e34a33","#fc8d59","#fdbb84","#fdd49e","#fef0d9"),
                     list('#045a8d','#2b8cbe','#74a9cf','#a6bddb','#d0d1e6','#f1eef6'),
                     list('#006d2c', '#2ca25f', '#66c2a4', '#99d8c9', '#ccece6','#edf8fb'),
                     list('#252525', '#636363', '#969696', '#bdbdbd', '#d9d9d9', '#f7f7f7')
                 ),
                 options = list(`toggle-palette-more-text` = "Show more")
             ),
             spectrumInput(
                 inputId = "downColor",
                 label = "Pick downregulated color:",
                 selected = '#045a8d',
                 width = "60%",
                 choices = list(
                     list("#b30000","#e34a33","#fc8d59","#fdbb84","#fdd49e","#fef0d9"),
                     list('#045a8d','#2b8cbe','#74a9cf','#a6bddb','#d0d1e6','#f1eef6'),
                     list('#006d2c', '#2ca25f', '#66c2a4', '#99d8c9', '#ccece6','#edf8fb'),
                     list('#252525', '#636363', '#969696', '#bdbdbd', '#d9d9d9', '#f7f7f7')
                 ),
                 options = list(`toggle-palette-more-text` = "Show more")
             ),
             circleButton(
               inputId = "informationVol1",
               icon = icon("info"),
               size = "xs",
               status = "primary"
             ),
             bsTooltip(
               "informationVol1",
               paste0("Add specific label for individual genes by including the name on the bar below."),
               trigger = "hover",
               placement = "right"
             ),
             uiOutput("geneSelector")
             )
             )
  ),
      column(
     width = 9,
     tabBox(
         width = 12,
         title = "",
         tabPanel(
                title = "Volcano plot",
                
                circleButton(
                  inputId = "informationVol",
                  icon = icon("info"),
                  size = "xs",
                  status = "primary"
                ),
                bsTooltip(
                  "informationVol",
                  paste0("Click over the dots to explore the gene values below the graph. ",
                         "The color scale for up/down genes and the statistics threshold ",
                         "are consistent with the previous selection above."),
                  trigger = "hover",
                  placement = "right"
                ),
                
                plotOutput("volcano", click = "plot_click1" , width = "100%", height = "600px"),
                column(width=8,
                  tableOutput("texto1")
                ),
                column(width = 4,
                       downloadButton("downVolcano","Download SVG")
                       )
            ),
         tabPanel(title = "KaryoPlot",
                                          circleButton(
                                            inputId = "karyoInfo",
                                            icon = icon("info"),
                                            size = "xs",
                                            status = "primary"
                                          ),
                                          bsTooltip(
                                            "karyoInfo",
                                            paste0("Genes up regulated will be shown above and down regolated below chromosome"),
                                            trigger = "hover",
                                            placement = "right"
                                          ),
                                          tagList(
                                            fluidRow(
                                              column(width=10, offset=1,
                                                     plotOutput("karyoPlot", height = "800px") )),
                                            downloadButton("downKrpt","Download PNG")
                                            )
                                          )
 )))

) # fin page



