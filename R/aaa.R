#' @import methods
NULL

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
    ifelse(nOccur != nStrains, 
           ifelse(nOccur == 1, 'Singleton', 'Accessory'), 
           'Core')
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
    if (sum(emptyGenes) == nrow(pangenome$pangenome)) {
        stop('Pangenome is empty')
    }
    list(
        pangenome = pangenome$pangenome[!emptyGenes, !emptyStrains, 
                                        drop = FALSE], 
        genes = pangenome$genes[!emptyGenes, , drop = FALSE]
    )
}