# Resources

* Carpeta specie (**Hs**, **Mm**, ...)

  º Subcarpetas **GO**, **KEGG**, **GSEA**, **cytoband**

Ejecutar el script createTree.sh con parámetro nombre de la especie para crear la estructura

```bash
$ createTree.sh Mm
```

## Crear recursos en carpeta GO y KEGG

* Ejecutar la función de *utils.R*  con la especie que sea

  ```R
  >getGOlevels(species="Mm")
  ```

  esto crea un Rds con los niveles de los GO

* Ejecutar:

  ```R
  > updateDatabases(species="Mm")
  ```

  Esto crea *GOlinks.Rds*,  *KeggLinks.Rds*, *KeggPathwayNames.Rds* en carpetas GO y KEGG de la especie en cuestion.

## Crear recursos GSEA

> TODO: gseaKegg ahora mismo sólo funciona para ratón

* Ejecutar siempre después de *updateDatabases()* :

  ```R
  > buildKeggDataset(specie="Mm")
  ```

## Crear recursos cytoband

* Ejecutar:

  ```R
  > cytoBandCreate(specie="Mm")
  ```

  Esto crear dos Rds con las coordenadas de los cromosomas y de las citobandas para el ideograma

* Ejecutar el script:

  ```bash
  $ createAnnotation.sh full_ftp_ensembl_path_gff3 Mm
  ```

  Los parámetros son url del ftp de ensembl al fichero gff3 (no gtf) el fichero debe ser:

  > ftp://ftp.ensembl.org/pub/release-100/gff3/mus_musculus/Mus_musculus.GRCm38.100.chr.gff3.gz

  ... de ese tipo (ojo hay que fijarse en que sea ese **specie.Assembly.ensVersion.chr.gff3.gz** )

  Esto descarga el fichero se queda con las columnas necesarias y lo mueve a su sitio.