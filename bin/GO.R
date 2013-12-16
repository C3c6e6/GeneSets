#!/usr/bin/env Rscript
library(GO.db)

"%d%" <- setdiff
"%i%" <- intersect
"%u%" <- union

getTermDefinitions <- function() {
    keySet = keys(GO.db, "GOID")
    conversionTable = select(GO.db, keys=keySet, keytype="GOID",
            columns="DEFINITION")
    definitions = conversionTable$DEFINITION
    names(definitions) <- conversionTable$GOID
    definitions
}

organismDBI2AnnotationTable <- function(annotationPackageName) {
    require(annotationPackageName, character.only=TRUE)
    message("Querying organismDBI...")
    organismDBIData = select(eval(parse(text=annotationPackageName)), 
            keys=keys(eval(parse(text=annotationPackageName)), "GOID"), 
            columns=c("ENTREZID", "TERM"), keytype="GOID")
    message("Constructing preliminary table...")
    organismDBIData = organismDBIData[!is.na(organismDBIData$ENTREZID),]
    preliminaryTable = data.frame(geneID = organismDBIData$ENTREZID, 
            termID = organismDBIData$GOID, termName = organismDBIData$TERM,
            dbName = organismDBIData$ONTOLOGY, stringsAsFactors=FALSE)
    message("Querying GO.db...")
    offspringList = c(as.list(GOBPOFFSPRING), as.list(GOCCOFFSPRING),
            as.list(GOMFOFFSPRING))
    goTerms = as.list(GOTERM)
    message("Constructing uncovered term table... ", appendLF=FALSE)
    omittedTermIDs = names(offspringList) %d% unique(preliminaryTable$termID)
    message("adding ", length(omittedTermIDs), " omitted terms")
    omittedTermTable = do.call(rbind, lapply(omittedTermIDs, function(x) { 
            t = goTerms[[x]]; 
            data.frame(geneID=NA, termID=GOID(t), termName=Term(t), 
                    dbName=Ontology(t), stringsAsFactors=FALSE)}))
    message("merging...")
    preliminaryTable = rbind(preliminaryTable, omittedTermTable)
    message("splitting...")
    tableSplit = split(preliminaryTable, preliminaryTable$termID)
    message("extending...")
    do.call(rbind, lapply(tableSplit, expandWithTermOffspring, tableSplit, 
            offspringList))
}

expandWithTermOffspring <- function(subTable, tableSplit, offspringList) {
    termID = unique(as.character(subTable$termID))
    termName = unique(as.character(subTable$termName))
    dbName = unique(as.character(subTable$dbName))
    offspring = offspringList[[termID]] %i% names(tableSplit)
    extension = do.call( rbind, lapply(offspring, function(x)
            data.frame(geneID = tableSplit[[x]]$geneID, termID = termID,
                    termName = termName, dbName = dbName,
                    stringsAsFactors = FALSE)))
    expandedTable = rbind(subTable, extension)
    expandedTable[!is.na(expandedTable$geneID),]
}

#organismFile = commandArgs(trailing=TRUE)[1]
annotationPackage = commandArgs(trailing=TRUE)[1]
packageDir = commandArgs(trailing=TRUE)[2]
#organisms = sub(" ", ".", readLines(organismFile))
termDefinitions = getTermDefinitions()
outputDir = sprintf("%s/GeneSets.%s/data/", packageDir, annotationPackage)
annotation = organismDBI2AnnotationTable(annotationPackage)
annotation$description = termDefinitions[annotation$termID]
for (db in unique(annotation$dbName)) {
    write.table(annotation[annotation$dbName == db,],
            sprintf("%s/GO%s.txt", outputDir, db), sep="\t", quote=FALSE,
            row.names=FALSE)
}

