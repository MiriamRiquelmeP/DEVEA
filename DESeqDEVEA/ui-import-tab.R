fluidPage(
    fluidRow(
        column(width = 4,
               circleButton(
                   inputId = "infoTEST",
                   icon = icon("info"),
                   size = "xs",
                   status = "primary"
               ),
               bsTooltip(
                   "infoTEST",
                   paste0("Select Wald's test or Likelihood Ratio Test.",
                          "Wald's test performs pairwise test using first category as reference.",
                          "LTR performs comparison between full and reduced model"),
                   trigger = "hover",
                   placement = "left"
               ),
               box(width = 12,
                #uiOutput("testVariable"),
                htmlOutput("textTestAlgorithm"),
                
                uiOutput("testAlgorithm"),
                uiOutput("testButton")
               )
               )) #,
    #     column(width = 8,
    #            box(width = 12, 
    #                title = "Sample table",
    #                DTOutput("coldataTable"))
    #            ) ),
    # fluidRow(
    #     column(width = 12,
    #            box(width = NULL,
    #                title = "Expression matrix (first 10 rows)",
    #                div(style = 'overflow-x: scroll', DT::dataTableOutput("expressionTable")) )
    #            )
    # )
)