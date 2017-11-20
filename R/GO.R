#' Download and store the Gene Ontology
#' 
#' This function downloads a copy of the gene ontology and formats it for quick 
#' access. It can optionally check the current version and only download a new
#' if it supersede the current one.
#' 
#' @param mode Either 'force' (default) or 'auto'. If 'force' the gene ontology
#' will get fetched no matter what. If 'auto' it will only get fetched if the
#' current local version is older than the one available on the net.
#' 
#' @return This function is called for its side effects
#' 
#' @references \url{http://geneontology.org}
#' 
#' @examples 
#' if (interactive()) {
#'     getGO('auto')
#' }
#' 
#' @importFrom utils download.file
#' 
#' @export
#' 
getGO <- function(mode = 'force') {
    get <- !hasGO()
    if (!get) {
        if (mode == 'force') {
            get <- TRUE
        } else if (mode == 'auto') {
            get <- !checkVersion()
        }
    }
    if (get) {
        oboPath <- file.path(
            system.file('extdata', package = 'PanVizGenerator'),
            'go-basic.obo'
        )
        download.file('http://geneontology.org/ontology/go-basic.obo', oboPath)
        message('Parsing obo file... ', appendLF = FALSE)
        go <- parseOBO(oboPath)
        cacheGO(go)
        message('Done')
    }
}

# INTERNALS

oboVersion <- function(obo) {
    if (grepl('https?://', obo)) {
        con <- url(obo)
    } else {
        con <- file(obo)
    }
    open(con)
    on.exit({
        close(con)
    })
    chunkSize <- 5
    start <- 0
    while (TRUE) {
        rec <- scan(con, character(), sep = '\n', skip = start, 
                    nlines = chunkSize, quiet = TRUE)
        if (any(grepl('data-version', rec))) {
            break
        } else {
            start <- start + chunkSize
        }
    }
    sub('data-version:\\s+', '', rec[grepl('data-version', rec)])
}
hasGO <- function() {
    oboFile <- system.file('extdata', 'go-basic.obo', 
                           package = 'PanVizGenerator')
    oboFile != ''
}
checkVersion <- function() {
    if (!hasGO()) return(FALSE)
    oboFile <- system.file('extdata', 'go-basic.obo', 
                           package = 'PanVizGenerator')
    currentVersion <- oboVersion(oboFile)
    publicVersion <- oboVersion('http://geneontology.org/ontology/go-basic.obo')
    as.Date(currentVersion, 'releases/%Y-%m-%d') >= 
        as.Date(publicVersion, 'releases/%Y-%m-%d')
}
parseOBO <- function(file) {
    obo <- readLines(file, warn = FALSE)
    obo <- obo[seq(which(obo == '[Term]')[1], length(obo))]
    elemStart <- grep('^\\[.+\\]$', obo)
    obo <- split(obo, rep.int(seq_along(elemStart), 
                              diff(c(elemStart, length(obo) + 1))))
    obo <- obo[sapply(obo, function(elem) elem[1] == '[Term]')]
    termInd <- rep.int(seq_along(obo), lengths(obo))
    obo <- unlist(obo, use.names = FALSE)
    type <- sub('^(\\w+):.*$', '\\1', obo)
    type[!grepl('^\\w+?:', obo)] <- NA
    typeSplit <- split(seq_along(type), type)
    obo <- sub('^\\w+?:\\s+', '', obo)
    obo <- sub('^\"(.+)\".*$', '\\1', obo)
    obo <- sub(' ! .+$', '', obo)
    
    vertices <- data.frame(
        id = obo[typeSplit$id],
        name = obo[typeSplit$name],
        namespace = obo[typeSplit$namespace],
        def = obo[typeSplit$def],
        alt_id = NA,
        comment = NA,
        is_obsolete = FALSE,
        subset = NA,
        synonym = NA,
        xref = NA,
        stringsAsFactors = FALSE
    )
    vertices$alt_id[unique(termInd[typeSplit$alt_id])] <- 
        split(obo[typeSplit$alt_id], termInd[typeSplit$alt_id])
    vertices$comment[termInd[typeSplit$comment]] <- 
        obo[typeSplit$comment]
    vertices$is_obsolete[termInd[typeSplit$is_obsolete]] <- 
        TRUE
    vertices$subset[unique(termInd[typeSplit$subset])] <- 
        split(obo[typeSplit$subset], termInd[typeSplit$subset])
    vertices$synonym[unique(termInd[typeSplit$synonym])] <- 
        split(obo[typeSplit$synonym], termInd[typeSplit$synonym])
    vertices$xref[unique(termInd[typeSplit$xref])] <- 
        split(obo[typeSplit$xref], termInd[typeSplit$xref])
    
    isAEdges <- data.frame(
        from = termInd[typeSplit$is_a],
        to = match(obo[typeSplit$is_a], vertices$id),
        type = 'is_a',
        stringsAsFactors = FALSE
    )
    considerEdges <- data.frame(
        from = termInd[typeSplit$consider],
        to = match(obo[typeSplit$consider], vertices$id),
        type = 'consider',
        stringsAsFactors = FALSE
    )
    replacedByEdges <- data.frame(
        from = termInd[typeSplit$replaced_by],
        to = match(obo[typeSplit$replaced_by], vertices$id),
        type = 'replaced_by',
        stringsAsFactors = FALSE
    )
    relationship <- unlist(strsplit(obo[typeSplit$relationship], ' '))
    relationshipEdges <- data.frame(
        from = termInd[typeSplit$relationship],
        to = match(relationship[
            seq(2, by = 2, length.out = length(relationship)/2)
        ], vertices$id),
        type = relationship[seq(1, by = 2, 
                                length.out = length(relationship)/2)],
        stringsAsFactors = FALSE
    )
    edges <- rbind(
        isAEdges,
        considerEdges,
        replacedByEdges,
        relationshipEdges
    )
    edges <- edges[!is.na(edges$from) & !is.na(edges$to), ]
    row.names(edges) <- NULL
    list(
        vertices = vertices,
        edges = edges
    )
}
cacheGO <- function(go) {
    file <- file.path(system.file('extdata', package = 'PanVizGenerator'), 
                      'go-cache.RDS')
    saveRDS(go, file)
}
loadGO <- function() {
    if (!hasGO()) {
        getGO()
    }
    file <- file.path(system.file('extdata', package = 'PanVizGenerator'), 
                      'go-cache.RDS')
    if (!file.exists(file)) {
        go <- parseOBO(system.file('extdata', 'go-basic.obo', 
                                   package = 'PanVizGenerator'))
        cacheGO(go)
    } else {
        go <- readRDS(file)
    }
    go
}
#' @importFrom igraph graph_from_data_frame V V<-
goToGraph <- function(go) {
    names(go$vertices) <- sub('^name$', '.name', names(go$vertices))
    go$vertices <- data.frame(name = seq_len(nrow(go$ver)), go$vertices, 
                              stringsAsFactors = FALSE)
    gr <- graph_from_data_frame(go$edges, vertices = go$vertices)
    V(gr)$name <- V(gr)$id
    gr
}
#' @importFrom igraph as_data_frame
graphToGo <- function(graph) {
    go <- as_data_frame(graph, 'both')
    go$edges$from <- match(go$edges$from, go$vertices$name)
    go$edges$to <- match(go$edges$to, go$vertices$name)
    go$vertices$name <- go$vertices$.name
    go$vertices$.name <- NULL
    edgeColOrder <- order(match(names(go$edges), c('from', 'to', 'type')))
    go$edges <- go$edges[, edgeColOrder]
    vertexColOrder <- order(match(names(go$vertices),
                                  c("id", "name", "namespace", "def", "alt_id",
                                    "comment", "is_obsolete", "subset",
                                    "synonym", "xref")))
    go$vertices <- go$vertices[, vertexColOrder]
    rownames(go$vertices) <- NULL
    rownames(go$edges) <- NULL
    go
}
#' @importFrom igraph induced_subgraph V E vertex_attr vertex_attr<- vertex_attr_names gorder
pruneGO <- function(go, namespace = NULL, subset = NULL, relations = NULL, 
                    metadata = NULL) {
    if (!inherits(go, 'igraph')) {
        go <- goToGraph(go)
        asList <- TRUE
    } else {
        asList <- FALSE
    }
    
    if (!is.null(namespace)) {
        possibles <- unique(V(go)$namespace)
        namespace <- unlist(lapply(namespace, function(n) {
            possibles[grep(n, possibles)]
        }))
        go <- induced_subgraph(go, which(V(go)$namespace %in% namespace))
    }
    if (!is.null(subset)) {
        possibles <- unique(unlist(V(go)$subset))
        subset <- unlist(lapply(subset, function(n) {
            possibles[grep(n, possibles)]
        }))
        subsetInd <- rep(seq_len(gorder(go)), lengths(V(go)$subset))
        subsetMatch <- unlist(V(go)$subset) %in% subset
        subsetMatch <- sapply(split(subsetMatch, subsetInd), any)
        go <- induced_subgraph(go, which(subsetMatch))
    }
    if (!is.null(relations)) {
        possibles <- unique(E(go)$type)
        relations <- unlist(lapply(relations, function(n) {
            possibles[grep(n, possibles)]
        }))
        go <- go - E(go)[!E(go)$type %in% relations]
    }
    if (!is.null(metadata)) {
        possibles <- vertex_attr_names(go)
        metadata <- unlist(lapply(metadata, function(n) {
            possibles[grep(n, possibles)]
        }))
        metadata <- c('name', unique(sub('^name$', '.name', metadata)))
        vertex_attr(go) <- vertex_attr(go)[metadata]
    }
    
    if (asList) {
        go <- graphToGo(go)
    }
    go
}