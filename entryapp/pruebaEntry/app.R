#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 



        # Show a plot of the generated distribution
        mainPanel(
            actionBttn(
                inputId = "simple",
                label = "Go Simple EnrichApp",
                style = "jelly",
                color = "primary",
            ),
            actionBttn(
                inputId = "enrich",
                label = "Go enrich EnrichApp",
                style = "jelly",
                color = "primary",
            )
        )
    )

# Define server logic required to draw a histogram
server <- function(input, output) {
    shinyjs::onclick("enrich", shinyjs::runjs("window.open('http://155.54.120.105/shiny/enrichappDark/','_blank')") )
    shinyjs::onclick("simple", shinyjs::runjs("window.open('http://155.54.120.105/shiny/enrich_listable/','_blank')") )
    
    

}

# Run the application 
shinyApp(ui = ui, server = server)
