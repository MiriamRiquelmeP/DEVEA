#!/bin/bash

cp -r www/ /srv/shiny-server/enrichappDark
cp -r resources/ /srv/shiny-server/enrichappDark
cp ./{app.R,report.css,README.md,report.Rmd,ui-go-tab.R,ui-gsea-tab.R,ui-import-tab.R,ui-kegg-tab.R,ui-preview-tab.R,utils.R,utilsReport.R,updatepopModals.R} /srv/shiny-server/enrichappDark/
