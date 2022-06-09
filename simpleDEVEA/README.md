# Variables enrichApp

## Pestaña preview

| Elemento         | Variables   |
| ---------------- | ----------- |
| tabla samples    | **samples** |
| tabla resultados | **preview** |
| gráfico PCA      | **pca**     |

## Pestaña kegg

| Elemento  | All_together       | upRegulated       | downregulated       |
| --------- | ------------------ | ----------------- | ------------------- |
| Tabla     | **tableAll**       | **table**         | **tableDown**       |
| Barplot   | **keggPlotAll**    | **keggPlot**      | **keggPlotDown**    |
| Chordplot | **keggChordAll**   | **keggChord**     | **keggChordDown**   |
| Dotplot   | **keggDotAll**     | **keggDotUp**     | **keggDotDown**     |
| HeatMap   | **heatmapKeggAll** | **heatmapKeggUp** | **heatmapKeggDown** |
| NetPlot   | **cnetKeggAll**    | **cnetKeggUp**    | **cnetKeggDown**    |
|           |                    |                   |                     |

## Pestaña GO BP

| ELEMENTO | all_together   | upregulated | downregulated   |
| -------- | -------------- | ----------- | --------------- |
| Tabla    | **tableBPall** | **tableBP** | **tableBPdown** |
| Barplot  | **plotBPall**  | **plotBP**  | **plotBPdown**  |
| Dotplot  | **BPDotall**   | **BPDotUp** | **BPDotdown**   |

## Pestaña GO MF

| ELEMENTO | all_together   | upregulated | downregulated   |
| -------- | -------------- | ----------- | --------------- |
| Tabla    | **tableMFall** | **tableMF** | **tableMFdown** |
| Barplot  | **plotMFall**  | **plotMF**  | **plotMFdown**  |
| Dotplot  | **MFDotall**   | **MFDotUp** | **MFDotdown**   |

## Pestaña GO CC

| ELEMENTO | all_together   | upregulated | downregulated   |
| -------- | -------------- | ----------- | --------------- |
| Tabla    | **tableCCall** | **tableCC** | **tableCCdown** |
| Barplot  | **plotCCall**  | **plotCC**  | **plotCCdown**  |
| Dotplot  | **CCDotall**   | **CCDotUp** | **CCDotdown**   |



## Pestaña GSEA

| ELEMENTO | All genes deseq |
| -------- | --------------- |
| Tabla    | **gseaTable**   |
| Barplot  | **gseaPlot**    |



## Variables observeEvent

| variable | valores                         | contenido              |
| -------- | ------------------------------- | ---------------------- |
| data     | $genesUp, $genesDown, $genesall | listados de genes      |
| goDT     | $up, $down                      | pretabla GO            |
| go       | $up, $down                      | tabla enrich GO        |
| kgg      | $up, $down, $all                | tabla enrich Kegg      |
| kggDT    | $up, $down, $all                | pretabla enrich Kegg   |
| datos    | $dds                            | objeto DESeq importado |
| gsea     | $gsea                           | objeto GSEA            |

## Variables reactive

| Variable   | contenido                          |
| ---------- | ---------------------------------- |
| rowsAll    | filas seleccionadas de tableAll    |
| rows       | filas seleccionadas de table       |
| rowdown    | filas seleccionadas de tableDown   |
| bprows     | filas seleccionadas de tableBP     |
| mfrows     | filas seleccionadas de tableMF     |
| ccrows     | filas seleccionadas de tableCC     |
| bprowndown | filas seleccionadas de tableBPdown |
| mfrowsdown | filas seleccionadas de tableMFdown |
| ccrowsdown | filas seleccionadas de tableCCdown |
| variables  | Variables seleccionadas para PCA   |
| gsearow    | filas seleccionadas para GSEA      |
| bprowsall  | filas seleccionadas de tableBPall  |
| mfrowsall  | filas seleccionadas de tableMFall  |
| ccrowsall  | filas seleccionadas de tableCCall  |
|            |                                    |



## Diario de acciones

01/03/2020: FPS: 

* Arreglado tabla samples para que no salgan las dos últimas columnas
* Arreglado tabla results para que salga con notación científica p-val ordenable
* Añadida y funcional, pestaña kegg: All DE genes
* Centrado boton download report

03/03/2020: FPS:

* Movido seleccionable de variables al body
* Primeras pruebas plantilla pdf

## TODOs

* Hacer all together en GO
* Poner en marcha report
* Sugerencias??¿??¿?