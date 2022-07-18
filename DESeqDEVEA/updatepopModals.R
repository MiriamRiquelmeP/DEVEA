selectPopUpModal <- function(session){
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalPreview",
        label = "Select preview elements to report",
        choices = c(
            "PCA",
            "BoxPlot",
            "Heatmap",
            "Cluster",
            "Top6",
            "Top1",
            "Karyoplot",
            "Volcano",
            "MA"
        ),
        selected = c(
            "PCA",
            "BoxPlot",
            "Heatmap",
            "Cluster",
            "Top6",
            "Top1",
            "Karyoplot",
            "Volcano",
            "MA"
        ),
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
        choices = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        selected = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
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
        choices = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        selected = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
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
        choices = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
        selected = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
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
##############################################
unselectPopUpModal <- function(session){
    updateCheckboxGroupButtons(
        size = "sm",
        session = session,
        inputId = "modalPreview",
        label = "Select preview elements to report",
        choices = c(
            "PCA",
            "BoxPlot",
            "Heatmap",
            "Cluster",
            "Top6",
            "Top1",
            "Karyoplot",
            "Volcano",
            "MA"
        ),
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
        choices = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
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
        choices = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
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
        choices = c("WordCloud", "Table", "Barplot", "Dotplot", "GObarplot", "GOcircleplot"),
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