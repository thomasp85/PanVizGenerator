#' @include aaa.R
#' @include IO.R
NULL

#' Convert pangenome into PanViz
#' 
#' This method contains the main functionality of PanVizGenerator. It takes the
#' pangenome data and properly formats it, combines it with the PanViz code and 
#' creates the PanViz files needed for the visualization. Per default everything 
#' is consolidated into the PanViz.html file but data and code can also be 
#' destributed into separate files. Currently pangenome data from 
#' \code{\link[FindMyFriends:FindMyFriends-package]{FindMyFriends}} is supported
#' natively (And FindMyFriends support grouping based on other algorithms), 
#' while more general use is supported by supplying a pangenome matrix along 
#' with functional annotation of each row. A last option is to provide the file
#' path to a csv file containing the pangenome matrix along with the functional
#' annotation.
#' 
#' @param object The object containing the pangenome data. If supplied as 
#' pangenome matrix it is expected that rows are gene groups and columns are 
#' genomes.
#' 
#' @param location The path to write the resulting PanViz files to.
#' 
#' @param name Depending on object either the name of the column with the gene 
#' group names or a vector of gene group names. See details.
#' 
#' @param go Depending on object either the name of the column with the gene 
#' group gene ontology annotation or a vector/list of gene group ontologies. 
#' See details.
#' 
#' @param ec Depending on object either the name of the column with the gene 
#' group E.C. annotation or a vector/list of gene group enzyme numbers. See 
#' details.
#' 
#' @param ignore Columns in the csv file to ignore, either given as column names
#' or indexes
#' 
#' @param consolidate Logical. Should all data and code be consolidated into the
#' PanViz.html file or spread out to multiple files.
#' 
#' @param showcase Logical. Should the resulting PanViz.html be opened in the 
#' default browser upon completion.
#' 
#' @param useDescription Logical. Should the description column in orgInfo 
#' be used in favor of group names (if description is NA it falls back to group 
#' name)
#' 
#' @param ... Parameters passed along to the clustering and dimensionality 
#' reduction functions. See details.
#' 
#' @details
#' The calculation of mds/pca as well as the hierarchical clustering can be 
#' controlled with the use of \code{dist} (default: 'canberra') for setting the
#' method used in the \code{\link[stats]{dist}} calls, \code{clust} (default: 
#' 'ward.D2') for setting the method used in the \code{\link[stats]{hclust}} 
#' calls, \code{center} (default: TRUE) to control whether variables should be 
#' centered prior to doing PCA and \code{scale} (default: TRUE) to control 
#' whether scaling should be performed prior to PCA.
#' 
#' When using \code{panviz} with a csv file the \code{name}, \code{go} and 
#' \code{ec} parameters should point to the columns in the csv file containing 
#' the respective information, either by name or index. If E.C. annotation is 
#' not given in the csv file it can be set to NA. For column with multiple 
#' possible values (go and ec) any delimiter can be used but ',', '.', and 
#' numbers (don't know why you would use numbers as delimiter anyway).
#' 
#' For panviz with a matrix or data.frame the \code{name}, \code{go} and 
#' \code{ec} parameters should contain the actual annotation as character 
#' vectors or lists of strings. For character vectors the same delimiting 
#' restrictions exists as for csv files described above. \code{ec} and 
#' \code{name} can be omitted. If \code{name} is missing the rownames of the 
#' matrix or data.frame will be used instead - if these are not present an error
#' will be thrown. \code{name}, \code{go} and \code{ec} must match the number of 
#' rows in the pangenome matrix if given.
#' 
#' @return NULL. This function is called for its side effects. If 
#' \code{showcase=TRUE} the resulting PanViz.html file will be opened in the 
#' default browser.
#' 
#' @seealso \code{\link{PanVizGenerator}} for a shiny interface to converting 
#' csv files.
#' 
#' @examples 
#' if(interactive()) {
#'     exampleFile <- system.file('extdata', 'exampleData.csv')
#'     panviz(exampleFile, location = tempdir(), showcase = TRUE)
#' }
#' 
#' @export
#' 
setGeneric('panviz', def = function(object, ...) {
        standardGeneric('panviz')
    }
)
#' @describeIn panviz Method for file paths
#' 
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' 
setMethod(
    'panviz', 'character',
    function(object, name = 'name', go = 'go', ec = 'ec', ignore, location, 
             consolidate = TRUE, showcase = FALSE, ...) {
        knownFiletypes <- c('csv')
        filetype <- file_ext(object)
        if (!filetype %in% knownFiletypes) {
            stop('Unknown filetype. Use: ', 
                 paste(knownFiletypes, collapse = ', '))
        }
        data <- tryCatch(
            switch(
                filetype,
                
                csv = read.csv(object, header = TRUE, check.names = FALSE, 
                               comment.char = '', stringsAsFactors = FALSE)
            ),
            error = function(e) {
                stop('Unable to parse file')
            }
        )
        if (nrow(data) == 0) {
            stop('Missing data in file')
        }
        if (!missing(ignore)) {
            if (inherits(ignore, 'character')) {
                ignore <- match(tolower(ignore), tolower(names(data)))
            }
            data <- data[, -ignore]
        }
        if (inherits(name, 'character')) {
            name <- which(tolower(names(data)) == tolower(name))
            if (length(name) != 1) {
                stop('name must uniquely match one column in object')
            }
        }
        if (inherits(go, 'character')) {
            go <- which(tolower(names(data)) == tolower(go))
            if (length(go) != 1) {
                stop('go must uniquely match one column in object')
            }
        }
        if (inherits(ec, 'character')) {
            ec <- which(tolower(names(data)) == tolower(ec))
            if (length(ec) > 1) {
                stop('ec must uniquely match one column in object')
            } else if (length(ec) == 0) {
                ec <- NA
            }
        }
        annotColumns <- c(name, go, ec)
        annotColumns <- annotColumns[!is.na(annotColumns)]
        panviz(
            data[, -annotColumns], 
            name = data[[name]], 
            go = data[[go]], 
            ec = data[[ec]], 
            location, 
            consolidate, 
            showcase, 
            ...
        )
    }
)
#' @describeIn panviz Method for pgVirtual subclasses from FindMyFriends
#' 
#' @importClassesFrom FindMyFriends pgVirtual
#' @importFrom FindMyFriends groupInfo groupNames
#' 
setMethod(
    'panviz', 'pgVirtual',
    function(object, location, useDescription = TRUE, consolidate = TRUE, 
             showcase = FALSE, ...) {
        gInfo <- groupInfo(object)
        if (all(is.na(gInfo$GO))) {
            stop('object not annotated with gene ontologies')
        }
        name <- ifelse(is.na(gInfo$description), groupNames(object), 
                       gInfo$description)
        mat <- as(object, 'matrix')
        panviz(mat, name = name, go = gInfo$GO, ec = gInfo$EC, 
               location = location, consolidate = consolidate, 
               showcase = showcase, ...)
    }
)
#' @describeIn panviz Method for pangenome matrix as numeric/integer matrix
#' @importFrom utils browseURL
#' 
setMethod(
    'panviz', 'matrix',
    function(object, name, go, ec, location, consolidate = TRUE, 
             showcase = FALSE, ...) {
        if (!dir.exists(location)) {
            if (!dir.create(location)) {
                stop('Unable to create folder: ', location)
            }
        }
        data <- formatData(object, name, go, ec)
        createPanData(data, location, ...)
        copyPanVizFiles(location)
        if (consolidate) {
            consolidateFiles(location)
        }
        if (showcase) {
            browseURL(file.path(location, 'PanViz.html'))
        }
        invisible(NULL)
    }
)
#' @describeIn panviz Method for pangenome matrix as data.frame (will be coerced 
#' to matrix)
#' 
setMethod(
    'panviz', 'data.frame',
    function(object, ...) {
        panviz(as.matrix(object), ...)
    }
)

# INTERNAL
formatGO <- function(go, empty = character()) {
    formatAnnot(go, '[^0-9]', 'GO:', empty = empty)
}
formatEC <- function(ec, empty = character()) {
    formatAnnot(ec, '[^0-9.]', 'EC:', empty = empty)
}
formatAnnot <- function(annot, nonChars, prefix, empty = character()) {
    nElem <- length(annot)
    if (!inherits(annot, 'list')) {
        annot <- strsplit(annot, paste0(nonChars, '+'))
    }
    annotInd <- rep(seq_along(annot), lengths(annot))
    annot <- unlist(annot)
    annotStrings <- !is.na(annot) & annot != ''
    annotInd <- annotInd[annotStrings]
    annot <- annot[annotStrings]
    annot <- sub(paste0('^', nonChars, '*'), prefix, annot)
    annotRes <- lapply(1:nElem, function(x) empty)
    annot <- split(annot, annotInd)
    annotRes[as.integer(names(annot))] <- annot
    annotRes
}
copyPanVizFiles <- function(location) {
    files <- list.files(system.file('PanViz', package = 'PanVizGenerator'), 
                        full.names = TRUE)
    file.copy(files, location, overwrite = TRUE)
}
consolidateFiles <- function(dir) {
    cssFile <- file.path(dir, 'style.css')
    jsFile <- file.path(dir, 'script.js')
    dataFile <- file.path(dir, 'data.js')
    htmlFile <- file.path(dir, 'PanViz.html')
    if (any(!file.exists(c(cssFile, jsFile, dataFile, htmlFile)))) {
        stop('Missing files. Dir must contain style.css, script.js, data.js and PanViz.html')
    }
    html <- readLines(htmlFile, warn = FALSE)
    cssLoc <- which(grepl('<link href="style.css" rel="stylesheet">', html))
    jsLoc <- which(
        grepl('<script type="text/javascript" src="script.js"></script>', html)
    )
    dataLoc <- which(
        grepl('<script type="text/javascript" src="data.js"></script>', html)
    )
    html[cssLoc] <- paste0(
        '<style>\n', 
        paste(readLines(cssFile, warn = FALSE), collapse = "\n"), 
        '\n</style>'
    )
    html[jsLoc] <- paste0(
        '<script>\n', 
        paste(readLines(jsFile, warn = FALSE), collapse = "\n"), 
        '\n</script>'
    )
    html[dataLoc] <- paste0(
        '<script>\n', 
        paste(readLines(dataFile, warn = FALSE), collapse = "\n"), 
        '\n</script>'
    )
    unlink(cssFile)
    unlink(jsFile)
    unlink(dataFile)
    write(html, file = htmlFile)
}
formatData <- function(object, name, go, ec) {
    if (!inherits(object, 'matrix')) stop('object must be a matrix')
    mode(object) <- 'integer'
    if (any(is.na(object))) {
        stop('object must be a matrix of integers or convertibel to one')
    }
    if (is.null(colnames(object))) {
        stop('columns in object must be named')
    }
    if (ncol(object) < 3) {
        stop('Pangenome must contain at least 3 genomes')
    }
    if (missing(name) || is.null(name)) {
        if (is.null(rownames(object))) {
            stop('Gene groups must be named')
        } else {
            name <- rownames(object)
        }
    }
    if (missing(ec) || is.null(ec)) {
        ec <- as.list(rep(NA, nrow(object)))
    }
    elementLengths <- c(nrow(object), length(name), length(go), length(ec))
    if (max(elementLengths) - min(elementLengths) != 0) {
        stop('Number of gene groups must correspond between object, name, go and ec')
    }
    go <- formatGO(go)
    ec <- formatEC(ec)
    data <- list(
        pangenome = as.data.frame(object), 
        genes = data.frame(
            name = unlist(name), 
            go = I(go), 
            ec = I(ec), 
            stringsAsFactors = FALSE
        )
    )
    stripEmpty(data)
}