fluidPage(
  h3("Gene Set Enrichment Analysis"),
  fluidRow(
    column(
      circleButton(
        inputId = "informationGSEA1",
        icon = icon("info"),
        size = "xs",
        status = "primary"
      ),
      bsTooltip(
        "informationGSEA1",
        paste0("Select the database to apply Gene Set Enrichment Analysis algorithm on it. ",
               "For further interpretation on the plot check ",
               a("GSEA user guide", 
                 href="https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html")),
        trigger = "hover",
        placement = "right"
      ),
      width = 3,
      box(
        title = "",
        width = 12,
        status = "info",
        uiOutput("gseaSelectize")
      )
    ),
    column(width = 9,
           box(title = "Table of GSEA pathways",
               solidHeader = FALSE,
               status = "primary",
               width = NULL,
               bsAlert("gsea"),
               DTOutput("gseaTable"),
               actionButton("resetgseaTable","Clear selection")
           )
    )
  ),
  fluidRow(column( 
    circleButton(
      inputId = "information21",
      icon = icon("info"),
      size = "xs",
      status = "primary"
    ),
    bsTooltip(
      "information21",
      "Select up to 3 pathways at the same time from the table above to visualize overlapping results.",
      trigger = "hover",
      placement = "left"
    ),
    width = 9, offset = 3,
    box(title = "GSEA plot",
        
        solidHeader = FALSE,
        status = "primary",
        width = NULL,
        bsAlert("gseaPlot"),
        plotOutput("gseaPlot"),
        fluidRow(column(width = 2,
                        downloadButton("gseaButton","Download SVG") ))
    )
  )
  )
)

