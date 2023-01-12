miBoxPlus <- function (..., title = NULL, footer = NULL, status = NULL, solidHeader = FALSE, 
    background = NULL, width = 6, height = NULL, collapsible = FALSE, 
    collapsed = FALSE, closable = TRUE, enable_label = FALSE, 
    label_text = NULL, label_status = "primary", enable_dropdown = FALSE, 
    dropdown_icon = "wrench", dropdown_menu = NULL, enable_sidebar = FALSE, 
    sidebar_content = NULL, sidebar_width = 25, sidebar_background = "#222d32", 
    sidebar_start_open = FALSE, footer_padding = TRUE) 
{
    if (sidebar_width < 0 || sidebar_width > 100) 
        stop("The box sidebar should be between 0 and 100")
    boxClass <- "box"
    if (solidHeader || !is.null(background)) {
        boxClass <- paste(boxClass, "box-solid")
    }
    if (!is.null(status)) {
        validateStatusPlus(status)
        boxClass <- paste0(boxClass, " box-", status)
    }
    if (collapsible && collapsed) {
        boxClass <- paste(boxClass, "collapsed-box")
    }
    if (!is.null(background)) {
        #validateColor(background)
        boxClass <- paste0(boxClass, " bg-", background)
    }
    if (enable_sidebar) {
        if (sidebar_start_open) {
            boxClass <- paste0(boxClass, " direct-chat direct-chat-contacts-open")
        }
        else {
            boxClass <- paste0(boxClass, " direct-chat")
        }
    }
    style <- NULL
    if (!is.null(height)) {
        style <- paste0("height: ", shiny::validateCssUnit(height))
    }
    titleTag <- NULL
    if (!is.null(title)) {
        titleTag <- shiny::tags$h3(class = "box-title", title)
    }
    boxToolTag <- NULL
    if (collapsible || closable) {
        boxToolTag <- shiny::tags$div(class = "box-tools pull-right")
    }
    collapseTag <- NULL
    if (collapsible) {
        buttonStatus <- status %OR% "default"
        collapseIcon <- if (collapsed) 
            "plus"
        else "minus"
        collapseTag <- shiny::tags$button(class = paste0("btn btn-box-tool"), 
            `data-widget` = "collapse", shiny::icon(collapseIcon))
    }
    closableTag <- NULL
    if (closable) {
        closableTag <- shiny::tags$button(class = "btn btn-box-tool", 
            `data-widget` = "remove", type = "button", shiny::tags$i(shiny::icon("times")))
    }
    labelTag <- NULL
    if (enable_label) {
        labelTag <- dashboardLabel(label_text, status = label_status)
    }
    dropdownTag <- NULL
    if (enable_dropdown) {
        dropdownTag <- shiny::tags$div(class = "btn-group", 
            shiny::tags$button(type = "button", class = "btn btn-box-tool dropdown-toggle", 
                `data-toggle` = "dropdown", shiny::icon(dropdown_icon)), 
            shiny::tagList(dropdown_menu))
    }
    sidebarTag <- NULL
    if (enable_sidebar) {
        sidebarTag <- shiny::tags$button(class = "btn btn-box-tool", 
            `data-widget` = "chat-pane-toggle", `data-toggle` = "tooltip", 
            `data-original-title` = "More", title = NA, type = "button", 
            shiny::icon("info"))
    }
    boxToolTag <- shiny::tagAppendChildren(boxToolTag, labelTag, 
        dropdownTag, sidebarTag, collapseTag, closableTag)
    headerTag <- NULL
    if (!is.null(titleTag) || !is.null(collapseTag)) {
        headerTag <- shiny::tags$div(class = "box-header", titleTag, 
            boxToolTag)
    }
    boxPlusTag <- shiny::tags$div(class = if (!is.null(width)) 
        paste0("col-sm-", width), shiny::tags$div(class = boxClass, 
        style = if (!is.null(style)) 
            style, headerTag, shiny::tags$div(class = "box-body", 
            ..., if (enable_sidebar) {
                shiny::tags$div(style = "z-index: 10000;", class = "direct-chat-contacts", 
                    shiny::tags$ul(class = "contacts-list", shiny::tags$li(style = paste0("width: ", 
                        sidebar_width, "%;"), sidebar_content)))
            }), if (!is.null(footer)) 
            shiny::tags$div(class = if (isTRUE(footer_padding)) 
                "box-footer"
            else "box-footer no-padding", footer)))
    translation_rate <- paste0(100 - sidebar_width, "%")
    shiny::tagList(shiny::singleton(shiny::tags$head(shiny::tags$style(shiny::HTML(paste0(".direct-chat-contacts {\n                 -webkit-transform: translate(100%, 0);\n                 -ms-transform: translate(100%, 0);\n                 -o-transform: translate(100%, 0);\n                 transform: translate(100%, 0);\n                 position: absolute;\n                 top: 0;\n                 bottom: 0;\n                 height: 100%;\n                 width: 100%;\n                 background: ", 
        sidebar_background, ";\n                 color: #fff;\n                 overflow: auto;\n              }\n              .direct-chat-contacts-open .direct-chat-contacts {\n                -webkit-transform: translate(", 
        translation_rate, ", 0);\n                -ms-transform: translate(", 
        translation_rate, ", 0);\n                -o-transform: translate(", 
        translation_rate, ", 0);\n                transform: translate(", 
        translation_rate, ", 0);\n              }\n              "))))), 
        boxPlusTag)
}

popupModal <- function(){
    modalDialog(
        title = "Who is visiting us",
        size = "m",
        fluidRow(
            column(width = 12,
                tmapOutput("distPlot", width = "550px", height = "400px"),
                )),
            tags$br(),
        fluidRow(
            column(width = 12,
                valueBoxOutput("visits"),
                actionButton("statButton", "Details", onclick = "window.open('report_stats.html','_blank')"))
        )
    )
}



mapData <- function(){
    require(tmap)
    require(tidyverse)
    data("World")
    js <- jsonlite::fromJSON("www/report_stats.json")
    kk <- js$geolocation$data
    df <- rbind_pages(kk$items) %>% as.matrix() %>% as.data.frame()
    df <- df %>% select(visitors.count, data)
    names(df) <- c("count","country")
    df$country <- sub(" ",";", df$country)
    df <- df %>% tidyr::separate(country, c("countryId", "country"), sep=";")
    world2 <- left_join(World, df, by = c("name" = "country") )
    world2 <- world2 %>%  mutate(count = as.numeric(as.character(count)) )
    return(world2)
}


