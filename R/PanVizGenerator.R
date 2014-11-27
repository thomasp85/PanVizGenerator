#' Parses compliant csv into a list with pangenome matrix and gene information
#' 
#' @param file A data.frame as provided by shinys fileInput widget
#' 
#' @return A list with two elements. pangenome is a data.frame with count data. 
#' It contains a column for each strain and a row for each gene group. genes is 
#' a data.frame with a row for each gene group and columns name, go and 
#' optionally ec
#' 
#' @noRd
#' 
parsePangenome <- function(file) {
    name <- file$name
    path <- file$datapath
    knownFileTypes <- c('csv')
    fileType <- tolower(tools::file_ext(name))
    if(!fileType %in% knownFileTypes) {
        stop('Unknown filetype. Use', knownFileTypes)
    }
    data <- tryCatch(
        expr = {
            switch(
                fileType,
                
                csv = read.csv(path, header = TRUE, check.names=FALSE, comment.char='', stringsAsFactors=FALSE)
            )
        },
        error = function(e) {
            stop('Unable to parse file')
        }
    )
    if(nrow(data) == 0) {
        stop('Missing data in file')
    }
    if(!all(c('name', 'go') %in% tolower(names(data)))) {
        stop('Missing gene information. Data must contain \'Name\' and \'GO\' columns')
    }
    descriptorCols <- tolower(names(data)) %in% c('name', 'go', 'ec')
    descriptors <- data[, descriptorCols]
    data <- data[, !descriptorCols, drop=FALSE]
    if(ncol(data) < 2) {
        stop('Data must contain at least two organisms')
    }
    if(!all(sapply(data, class) %in% c('numeric', 'integer'))) {
        stop('Pangenome data must be a numeric matrix')
    }
    stripEmpty(list(pangenome=data, genes=descriptors))
}

#' Create a javascript file defining the necessary pangenomic data
#' 
#' @param data A list as returned by parsePangenome
#' 
#' @param location The location to put the file
#' 
#' @param dist The distance measure to use
#' 
#' @param clust The clustering function to use
#' 
#' @param center logical. Should data be centered in the pca
#' 
#' @param scale logical. Should data be scaled in the pca
#' 
#' @return Used for the side effect
#' 
#' @noRd
#' 
createPanData <- function(data, location, dist, clust, center, scale) {
    varNames <- c('obo', 'dimReduc', 'root', 'pan', 'geneInfo')
    obo <- system.file('extdata', 'go-basic.obo', package = 'PanVizGenerator')
    obo <- scan(obo, what='character', quiet=TRUE, sep='\n')
    obo <- obo[!grepl('^synonym|^xref', obo)]
    varData <- c(
        toJSON(obo),
        createScatter(data, dist, center, scale),
        createCluster(data, dist, clust),
        createPangenome(data),
        createGeneInfo(data)
    )
    data <- paste0('var ', varNames, '=', varData, ';')
    
    write(data, file=file.path(location, 'data.js'))
}

#' Create json string with pangenome matrix information
#' 
#' @param data A list as returned by parsePangenome
#' 
#' @return A string with a json representation of data
#' 
#' @importFrom jsonlite toJSON
#' 
#' @noRd
#' 
createPangenome <- function(data) {
    toJSON(as.list(data$pangenome), pretty=FALSE)
}

#' Create json string with gene group information
#' 
#' @param data A list as returned by parsePangenome
#' 
#' @return A string with a json representation of data
#' 
#' @importFrom jsonlite toJSON
#' 
#' @noRd
#' 
createGeneInfo <- function(data) {
    genes <- data$genes
    names(genes) <- tolower(names(genes))
    genes$domain <- setDomain(data$pangenome)
    genes$name[is.na(genes$name) | genes$name == ''] <- 'Unknown'
    if(is.null(genes$ec)) genes$ec <- ''
    genes$go[is.na(genes$go)] <- ''
    genes$ec[is.na(genes$ec)] <- ''
    
    genes$go <- strsplit(genes$go, ';\\s*')
    genes$ec <- strsplit(genes$ec, ';\\s*')
    
    genes <- lapply(1:nrow(genes), function(i) {
        list(
            name=unbox(genes$name[i]),
            go=genes$go[[i]],
            ec=genes$ec[[i]],
            domain=unbox(genes$domain[i])
        )
    })
    names(genes) <- NULL
    
    toJSON(genes, pretty=FALSE)
}

#' Create json string with pca and mds information
#' 
#' @param data A list as returned by parsePangenome
#' 
#' @param dist The distance measure to use
#' 
#' @param center logical. Should data be centered in the pca
#' 
#' @param scale logical. Should data be scaled in the pca
#' 
#' @return A string with a json representation of data
#' 
#' @importFrom jsonlite toJSON
#' @importFrom pcaMethods pca
#' 
#' @noRd
#' 
createScatter <- function(data, dist, center, scale) {
    ### Multidimensional scaling
    mds <- cmdscale(d=dist(t(data$pangenome), method=dist), k=2)
    mds <- split(data.frame(name=rownames(mds), x=mds[,1], y=mds[,2], stringsAsFactors=F), 1:nrow(mds))
    mds <- lapply(mds, unbox)
    names(mds) <- NULL
    
    ### PCA
    pca <- pca(t(data$pangenome), method='svd', scale=ifelse(scale, 'uv', 'none'), center=center)
    pca <- split(data.frame(name=rownames(pca@scores), x=pca@scores[,1], y=pca@scores[,2], stringsAsFactors=F), 1:nrow(pca@scores))
    pca <- lapply(pca, unbox)
    names(pca) <- NULL
    
    toJSON(list(MDS=mds, PCA=pca), pretty=FALSE)
}

#' Create json string with pangenome matrix information
#' 
#' @param data A list as returned by parsePangenome
#' 
#' @param dist The distance measure to use
#' 
#' @param clust The clustering function to use
#' 
#' @return A string with a json representation of data
#' 
#' @importFrom jsonlite toJSON
#' 
#' @noRd
#' 
createCluster <- function(data, dist, clust) {
    den <- hclust(d=dist(t(data$pangenome), method=dist), method=clust)
    den <- denToList(den)
    
    toJSON(den, pretty=FALSE)
}

#' Recursively converts a dendrogram to a list
#' 
#' @param tree Clustering data matching the structure of the dendrogram class
#' but without a class attribute
#' 
#' @return A list that describes the hierarchical relationship of the data. Each
#' non-leaf node contains the elements height, leaf and children. Each leaf node
#' contains the elements name, height and leaf
#' 
#' @noRd
#' 
formatNode <- function(tree){
    ans <- list()
    if(length(tree) == 1){
        ans$name <- unbox(attr(tree, 'label'))
        ans$height <- unbox(attr(tree, 'height'))
        ans$leaf <- unbox(TRUE)
    } else {
        ans$height <- unbox(attr(tree, 'height'))
        ans$leaf <- unbox(FALSE)
        ans$children <- lapply(tree, formatNode)
    }
    ans
}

#' Convert a clustering to a list for JSON serialization
#' 
#' @param den A dendrogram of class dendrogram or which can be coerced to such
#' 
#' @return A list that describes the hierarchical relationship of the data. Each
#' non-leaf node contains the elements height, leaf and children. Each leaf node
#' contains the elements name, height and leaf
#' 
#' @noRd
#' 
denToList <- function(den){
    if(class(den) != 'dendrogram'){
        den <- as.dendrogram(den)
    }
    class(den) <- NULL
    formatNode(den)
}

#' Calculate domain membership of a pangenome
#' 
#' @param pangenome A pangenome matrix with strains as columns and gene groups 
#' as rows
#' 
#' @return A vector containing 'Core', 'Accessory' or 'Singleton' for each row
#' in the input
#' 
#' @noRd
#' 
setDomain <- function(pangenome) {
    nStrains <- ncol(pangenome)
    nOccur <- apply(pangenome, 1, function(x) sum(x != 0))
    ifelse(nOccur != nStrains, ifelse(nOccur == 1, 'Singleton', 'Accessory'), 'Core')
}

#' Removes empty gene groups and strains from pangenome data
#' 
#' @param pangenome A list as returned by parsePangenome
#' 
#' @return A list as the input but with empty genes and strains removed
#' 
#' @noRd
#' 
stripEmpty <- function(pangenome) {
    emptyGenes <- apply(pangenome$pangenome, 1, sum) == 0
    emptyStrains <- apply(pangenome$pangenome, 2, sum) == 0
    if(sum(emptyGenes) == nrow(pangenome$pangenome)) {
        stop('Pangenome is empty')
    }
    list(pangenome=pangenome$pangenome[!emptyGenes, !emptyStrains, drop=FALSE], genes=pangenome$genes[!emptyGenes, , drop=FALSE])
}