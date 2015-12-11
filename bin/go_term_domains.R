#!/usr/bin/env Rscript

library(GO.db)

keySet = keys(GO.db, keytype="GOID")
table = select(GO.db, keyType="GOID", keys=keySet, columns="ONTOLOGY")
write.table(table, "go_term_domains.txt", row.names=FALSE, quote=FALSE,
		sep="\t")
