Annotation of ORFs
---
* EggNog

eggNOG release 5.0 was used as the reference for annotation.

 - Source: http://eggnog5.embl.de/download/eggnog_5.0/

eggNOG-mapper v2.1.7 was used to perform annotation. Command:

```
emapper.py --cpu # -i input.faa -o output
```

---
* KEGG

KOfam release 2022-03-01 was used as the reference for annotation.

 - Source: https://www.genome.jp/ftp/db/kofam/archives/2022-03-01/

KofamScan v1.3.0 was used to perform annotation. Command:

```
exec_annotation -f detail-tsv --no-report-unannotated -o output.tsv input.faa
```

KEGG release 102.0+ (2022-05-03) was used to construct high-level hierarchies
of functional catalogs.

 - Source: From the KEGG website using the REST API.


