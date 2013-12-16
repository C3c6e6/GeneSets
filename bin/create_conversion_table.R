#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE);
annotationPackage = args[1]
from = args[2]
to = args[3]

require(annotationPackage, character.only=TRUE)
keySet = keys(eval(parse(text=annotationPackage)), from)
conversionTable = select(eval(parse(text=annotationPackage)), keys=keySet,
	keytype=from, columns=to)

write.table(conversionTable, stdout(), row.names=FALSE, quote=FALSE, sep="\t");
