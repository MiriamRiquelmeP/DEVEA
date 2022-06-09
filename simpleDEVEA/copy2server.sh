#!/bin/bash

cp -r www/ /srv/shiny-server/enrich_listable
cp -r resources/ /srv/shiny-server/enrich_listable
cp ./{app.R,global.R,report.css,README.md,report1col.Rmd,report.Rmd,ui-go-tab.R,ui-gsea-tab.R,ui-kegg-tab.R,ui-preview-tab.R,utils.R,utilsReport.R,UpdatepopModals.R} /srv/shiny-server/enrich_listable/
