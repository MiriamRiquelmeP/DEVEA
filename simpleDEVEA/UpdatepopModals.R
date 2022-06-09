unselectPopUpModal3 <- function(session){
    updateCheckboxGroupButtons(
        session = session,
        size = "sm",
        inputId = "modalPreview",
        label = "Select preview elements to report",
        choices = c("Karyoplot", "Volcano"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        session = session,
        size = "sm",
        inputId = "modalkeggAll",
        label = "Select elements to report Kegg All",
        choices = c("Table", "Barplot", "Chorplot", "Dotplot",
                    "Heatmap", "Netplot"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalkeggUp",
        label = "Select elements to report Kegg Up",
        choices = c("Table", "Barplot", "Chorplot", "Dotplot",
                    "Heatmap", "Netplot"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalkeggDown",
        label = "Select elements to report Kegg Down",
        choices = c("Table", "Barplot", "Chorplot", "Dotplot",
                    "Heatmap", "Netplot"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalGOAll",
        label = "Select elements to report GO All",
        choices = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalGOUp",
        label = "Select elements to report GO Up",
        choices = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalGODown",
        label = "Select elements to report GO Down",
        choices = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalGSEA",
        label = "Select elements to report GSEA",
        choices = c("Table", "GSEA plot"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
}

#######################################################
selectPopUpModal3 <- function(session){
    updateCheckboxGroupButtons(
        session = session,
        size = "sm",
        inputId = "modalPreview",
        label = "Select preview elements to report",
        choices = c("Karyoplot", "Volcano"),
        selected = c("Karyoplot", "Volcano"),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        session = session,
        size = "sm",
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
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
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
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
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
    
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalGOAll",
        label = "Select elements to report GO All",
        choices = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        selected = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalGOUp",
        label = "Select elements to report GO Up",
        choices = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        selected = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalGODown",
        label = "Select elements to report GO Down",
        choices = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        selected = c("WordCloud","Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
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
}

###################################
selectPopUpModal1 <- function(session){
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
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
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalGOAll",
        label = "Select elements to report GO All",
        choices = c("WordCloud","Table", "Barplot", "Dotplot"),
        selected = c("WordCloud","Table", "Barplot", "Dotplot"),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
}

###################################
unselectPopUpModal1  <- function(session){
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalkeggAll",
        label = "Select elements to report Kegg All",
        choices = c("Table", "Barplot", "Chorplot", "Dotplot",
                    "Heatmap", "Netplot"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalGOAll",
        label = "Select elements to report GO All",
        choices = c("WordCloud","Table", "Barplot", "Dotplot"),
        selected = c(),
        status = "primary",
        checkIcon = list(
            yes = icon("ok",
                       lib = "glyphicon"),
            no = icon("remove",
                      lib = "glyphicon")
        )
    )
}


