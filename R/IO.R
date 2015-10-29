#' @include aaa.R
NULL

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
createPanData <- function(data, location, dist = 'canberra', clust = 'ward.D2', 
                          center = TRUE, scale = TRUE) {
    varNames <- c('go', 'dimReduc', 'root', 'pan', 'geneInfo')
    data <- checkTerms(data)
    varData <- c(
        createGO(),
        createScatter(data, dist, center, scale),
        createCluster(data, dist, clust),
        createPangenome(data),
        createGeneInfo(data)
    )
    data <- paste0('var ', varNames, '=', varData, ';')
    
    write(data, file = file.path(location, 'data.js'))
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
    toJSON(as.list(data$pangenome), pretty = FALSE)
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
    
    genes <- lapply(1:nrow(genes), function(i) {
        list(
            name = unbox(genes$name[i]),
            go = genes$go[[i]],
            ec = genes$ec[[i]],
            domain = unbox(genes$domain[i])
        )
    })
    names(genes) <- NULL
    
    toJSON(genes, pretty = FALSE)
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
#' @importFrom jsonlite toJSON unbox
#' @importFrom pcaMethods pca
#' 
#' @noRd
#' 
createScatter <- function(data, dist = 'canberra', center = TRUE, 
                          scale = TRUE) {
    ### Multidimensional scaling
    mds <- cmdscale(d = dist(t(data$pangenome), method = dist), k = 2)
    mds <- split(data.frame(name = rownames(mds), x = mds[,1], y = mds[,2], 
                            stringsAsFactors = FALSE), 
                 1:nrow(mds))
    mds <- lapply(mds, unbox)
    names(mds) <- NULL
    
    ### PCA
    pca <- pca(t(data$pangenome), method = 'svd', 
               scale = ifelse(scale, 'uv', 'none'), 
               center = center)
    pca <- split(data.frame(name = rownames(pca@scores), x = pca@scores[,1], 
                            y = pca@scores[,2], stringsAsFactors = FALSE), 
                 1:nrow(pca@scores))
    pca <- lapply(pca, unbox)
    names(pca) <- NULL
    
    toJSON(list(MDS = mds, PCA = pca), pretty = FALSE)
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
createCluster <- function(data, dist = 'canberra', clust = 'ward.D2') {
    den <- hclust(d = dist(t(data$pangenome), method = dist), method = clust)
    den <- denToList(den)
    
    toJSON(den, pretty = FALSE)
}

#' Create GO data
#' 
#' @importFrom jsonlite toJSON
#' 
#' @noRd
#' 
createGO <- function() {
    go <- loadGO()
    go <- pruneGO(go, relations = c('is_a', 'replaced_by'), 
                  metadata = c('id', 'name', 'namespace', 'def', 'is_obsolete', 
                               'alt_id', 'subset'))
    toJSON(go, dataframe = 'columns', rownames = FALSE)
}

# HELPERS

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
    if (length(tree) == 1) {
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
    if (class(den) != 'dendrogram') {
        den <- as.dendrogram(den)
    }
    class(den) <- NULL
    formatNode(den)
}
availableTerms <- function(go) {
    unique(c(go$vertices$id, unlist(go$vertices$alt_id)))
}
checkTerms <- function(data) {
    go <- loadGO()
    terms <- unique(unlist(data$genes$go))
    if (!all(terms %in% availableTerms(go))) {
        if (!checkVersion()) {
            message('Missing terms in current GO. Getting latest version...')
            getGO()
            go <- loadGO()
        }
        missingTerms <- terms[!terms %in% availableTerms(go)]
        if (length(missingTerms) > 0) {
            message('Removing unknowm terms:\n', 
                    paste(missingTerms, collapse = '\n'))
            goInd <- rep(seq_along(data$genes$go), lengths(data$genes$go))
            goTerms <- unlist(data$genes$go)
            keep <- goTerms %in% go$vertices$id
            goInd <- goInd[keep]
            goTerms <- goTerms[keep]
            data$genes$go <- lapply(seq_len(nrow(data$genes)), function(x) {
                character()
            })
            data$genes$go[unique(goInd)] <- split(goTerms, goInd)
        }
    }
    data
}