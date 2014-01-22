#' Parse BioPax level 3 data files: following databases are supported:
#' INOH, BioModel, NCI-PID (only NCI-Nature curated data) and PANTHER
#'
#' @param inputFolder Path to folder which contains the .owl files (e.g.: "Bioinformatics/testData/PANTHER").
#' In all NCI-Nature curated data files the '&' needs to be removed, as '&' are not supported in XML.
#' @param dbName Name of the database (e.g.: INOH, BioModel, NCI-PID, PANTHER)
#' @param organismFile File with organism names (e.g.: "Bioinformatics/GeneSets/organisms")
#' @param organismTranslationFile File with organism translation (e.g.: "Bioinformatics/GeneSets/names.dmp"
#' @param outFolder Path to folder containing organsim folders (e.g.: "Bioinformatics/testData/Output/").
#' Output is written to folder/geneSets.species/data/dbName.txt' files
#' @author Heidi Lischer

library(rBiopaxParser)
library(XML)

parse_BioPax3 <- function(inputFolder, dbName, organismFile, organismTranslationFile, outFolder){
    owlFiles <- list.files(inputFolder, "*.owl", full.names=T)

    organismList <- read.delim(organismFile, header=F, stringsAsFactors=F)
    orgTransl <- read.delim(organismTranslationFile, stringsAsFactors=F)
    organismList$taxId <- sapply(organismList[,1], function(x){orgTransl[orgTransl[,3] == x,][1,1]} )
    organismList$outFolder <- sapply(organismList[,1], function(x){paste(outFolder, "geneSets.", gsub(" ", ".", x), "\\data", sep="")})
    rownames(organismList) <- organismList$taxId

    #make sure to write into a new file if parsing the fisrt .owl file
    newFile <- TRUE

    for(owlFile in owlFiles){
        message(paste("Parse ", owlFile, "...", sep=""))
        #get pathways and their genes
        biopax <- readBiopax(owlFile, verbose=F)
        pathways <- listPathways(biopax)
        if(!is.null(pathways)){
            onePathway <- T
            if(dbName == "PANTHER" & nrow(pathways) > 1){
                onePathway <- F
            }
            geneList <- lapply(pathways$id, getUniProtIds, biopax, dbName, onePathway)
            names(geneList) <- pathways$id

            #get gene translations and organisms
            translation <- uniProtIdTranslation(unique(stack(geneList)$values))
            translation$taxId <- getOrganism(translation$geneID)[rownames(translation),"taxId"]

            geneList <- lapply(geneList, function(x, translation){translation[translation$uniProtId %in% x,]},translation)

            #store data into list
            for(i in 1:nrow(organismList)){
                organism <- organismList$taxId[i]
                organismDataList <- data.frame()
                for(j in 1:nrow(pathways)){
                    pathwayId <- pathways[j,"id"]
                    pathwayGenes <- geneList[[pathwayId]]
                    pathwayGenes <- pathwayGenes[pathwayGenes$taxId == organism,]

                    if(nrow(pathwayGenes) > 0){
                        pathwayGenes$termID <- pathwayId
                        pathwayGenes$termName <- pathways[j, "name"]
                        pathwayGenes$dbName <- dbName
                        pathwayGenes$description <- ""
                        organismDataList <- rbind(organismDataList, pathwayGenes[, c("geneID", "termID", "termName", "dbName", "description")])
                    }
                }

                #write to the output file
                if(newFile){
                    newFile <- FALSE
                    if(!file.exists(organismList$outFolder[i])) stop(paste("folder",organismList$outFolder[i], "does not exist!"))
                    write.table(organismDataList, paste(organismList$outFolder[i], "\\", dbName ,".txt", sep=""),
                                row.names=FALSE, quote=FALSE, sep="\t")
                } else {
                    write.table(organismDataList, paste(organismList$outFolder[i], "\\", dbName ,".txt", sep=""),
                                row.names=FALSE, quote=FALSE, sep="\t", append=TRUE, col.names=FALSE)
                }
            }
        }
    }
}


getUniProtIds <- function(pathwayId, biopax, dbName, onePathway){
    geneSets <- pathway2Geneset(biopax, pwid=pathwayId)
    if(dbName == "PANTHER"){
        if(onePathway & (is.null(geneSets) || nrow(geneSets[geneSets$class == "Protein",]) == 0)){
            UniProtIds <- parsePANTHER(listInstances(biopax, class="Protein"), biopax)
        } else {
            UniProtIds <- parsePANTHER2(geneSets, biopax)
        }
    } else {
        if(dbName == "BioModel"){
            UniProtIds <- parseBioModel(geneSets, biopax)
        } else if (dbName == "NCI-PID"){
            UniProtIds <- parseNCIPID(geneSets, biopax)
        } else {
             UniProtIds <- parseINOH(geneSets, biopax)
        }
    }
    return(UniProtIds)
}


parseINOH <- function(geneSets, biopax){
    UniProtIds <- vector()
    if(!is.null(geneSets)){
        for(i in 1:nrow(geneSets)){
            refs <- getReferencedIDs(biopax, id=geneSets$id[i],  recursive = F, onlyFollowProperties="xref")
            UniProtIds <- c(UniProtIds, gsub("UniProt_", "", refs[grepl("UniProt_", refs)], "_"))
        }
        UniProtIds <- unique(UniProtIds)
    }
    return(UniProtIds)
}


parseBioModel <- function(geneSets, biopax){
    UniProtIds <- vector()
    if(!is.null(geneSets)){
        for(i in 1:nrow(geneSets)){
            refs <- getReferencedIDs(biopax, id=geneSets$name[i],  recursive = F, onlyFollowProperties="xref")
            UniProtIds <- c(UniProtIds, gsub("http://identifiers.org/uniprot/", "", refs[grepl("uniprot", refs)], "_"))
        }
        UniProtIds <- unique(UniProtIds)
    }
    return(UniProtIds)
}


parseNCIPID <- function(geneSets, biopax){
    UniProtIds <- vector()
    if(!is.null(geneSets)){
        for(i in 1:nrow(geneSets)){
            refs <- getXrefAnnotations(biopax, id=geneSets$id[i])
            newIds <- gsub("UniProt:", "", refs[grepl("UniProt:", refs)], "_")
            if(length(newIds) == 0){
                newIds <- refs[grepl("EntrezGene:", refs)]
            }
            UniProtIds <- c(UniProtIds, newIds)
        }
        UniProtIds <- unique(UniProtIds)
    }
    return(UniProtIds)
}


parsePANTHER <- function(proteinList, biopax){
    UniProtIds <- vector()
    if(!is.null(proteinList)){
        for(i in 1:nrow(proteinList)){
            refs <- getXrefAnnotations(biopax, id=proteinList$id[i])
            UniProtIds <- c(UniProtIds, gsub("uniprot:", "", refs[grepl("uniprot:", refs)], "_"))
        }
        UniProtIds <- unique(UniProtIds)
    }
    return(UniProtIds)
}


parsePANTHER2 <- function(geneSets, biopax){
    UniProtIds <- vector()
    if(!is.null(geneSets)){
        for(i in 1:nrow(geneSets)){
            refs <- getXrefAnnotations(biopax, id=geneSets$id[i])
            UniProtIds <- c(UniProtIds, gsub("uniprot:", "", refs[grepl("uniprot:", refs)], "_"))
        }
        UniProtIds <- unique(UniProtIds)
    }
    return(UniProtIds)
}


#translate UniProtId to ENTREZ gene id
uniProtIdTranslation <- function(UniProtIds){
    #only get 80 entries at once (else there is a possibility for an error)
    urlPath <- "http://www.uniprot.org/mapping/?query="
    urlOptions <- "&from=ACC+ID&to=P_ENTREZGENEID&format=tab"
    translation <- data.frame()
    index <- 1
    while(index+80 < length(UniProtIds)){
        translation <- rbind(translation,read.delim(paste(urlPath,paste(UniProtIds[index:(index+79)], collapse="+"),urlOptions, sep=""),stringsAsFactors=F))
        index <- index + 80
    }
    translation <- rbind(translation,read.delim(paste(urlPath,paste(UniProtIds[index:length(UniProtIds)], collapse="+"),urlOptions, sep=""),stringsAsFactors=F))
    #add genes already in entrez
    entrezGenes <- UniProtIds[grepl("EntrezGene:", UniProtIds)]
    entrezIds <- gsub("EntrezGene:", "", entrezGenes, "_")
    translation <- rbind(translation, data.frame(From=entrezGenes, To=entrezIds))

    colnames(translation) <- c("uniProtId", "geneID")
    #remove duplicated entrez ids
    translation <- translation[!duplicated(translation$geneID),]
    rownames(translation) <- translation$geneID
    return(translation)
}



#get organism information from ENTREZ ids
getOrganism <- function(entrezIds){
    urlPath <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id="
    geneString <- paste(entrezIds, collapse=",")
    xmlDoc <- xmlTreeParse(paste(urlPath, geneString, sep=""), isURL=TRUE, useInternalNodes = TRUE)

    #parse XML
    geneNodeSet <- getNodeSet(xmlDoc, "//DocSum")
    translation <- data.frame()
    if(length(geneNodeSet)!=0){
        translation <- data.frame(
            geneID=unlist(xmlSApply(geneNodeSet, function(x) xpathApply(x, c("Id"), xmlValue))),
            taxId=unlist(xmlSApply(geneNodeSet, function(x) xpathApply(x, c("Item[@Name='TaxID']"), xmlValue))),
            stringsAsFactors=FALSE)
    }

    #add genes with no query results
    if(nrow(translation) != length(entrezIds)){
        genesMissing <- geneSet[!entrezIds %in% translation$geneID]
        translation <-rbind(translation, data.frame(geneID=genesMissing, taxId=NA))
    }
    rownames(translation) <- translation$geneID
    return(translation)
}
