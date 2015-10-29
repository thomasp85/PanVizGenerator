context("Gene ontology manipulation")

goFile <- file.path(system.file('extdata', package = 'PanVizGenerator'), 
                    'go-basic.obo')
unlink(goFile)
cacheFile <- file.path(system.file('extdata', package = 'PanVizGenerator'), 
                    'go-cache.obo')
unlink(cacheFile)

test_that("GO existance query works", {
    expect_false(hasGO())
    expect_false(checkVersion())
    getGO()
    expect_true(hasGO())
    expect_true(checkVersion())
})

test_that("GO loading works", {
    GO <- loadGO()
    unlink(cacheFile)
    GO1 <- loadGO()
    expect_equal(GO, GO1)
    expect_is(GO, 'list')
    expect_equal(length(GO), 2)
    expect_named(GO, c('edges', 'vertices'), ignore.order = TRUE)
    expect_named(GO$edges, c('to', 'from', 'type'), ignore.order = TRUE)
    expect_named(GO$vertices, c("id", "name", "namespace", "def", "alt_id", "comment", "is_obsolete", 
                                "subset", "synonym", "xref"), ignore.order = TRUE)
    expect_is(GO$edges$from, 'integer')
    expect_is(GO$edges$to, 'integer')
    expect_is(GO$edges$type, 'character')
    expect_is(GO$vertices$id, 'character')
    expect_is(GO$vertices$name, 'character')
    expect_is(GO$vertices$namespace, 'character')
    expect_is(GO$vertices$def, 'character')
    expect_is(GO$vertices$alt_id, 'list')
    expect_is(GO$vertices$comment, 'character')
    expect_is(GO$vertices$is_obsolete, 'logical')
    expect_is(GO$vertices$subset, 'list')
    expect_is(GO$vertices$synonym, 'list')
    expect_is(GO$vertices$xref, 'list')
})

test_that("Graph conversion works", {
    library(igraph)
    GO <- loadGO()
    GOgraph <- goToGraph(GO)
    expect_is(GOgraph, 'igraph')
    expect_equal(GO, graphToGo(GOgraph))
    expect_equal(gorder(GOgraph), nrow(GO$vertices))
    expect_equal(gsize(GOgraph), nrow(GO$edges))
    expect_true(is.directed(GOgraph))
    expect_true(is.named(GOgraph))
    expect_true(is.dag(GOgraph))
    expect_equal(V(GOgraph)$name, GO$vertices$id)
    expect_is(GO, class(pruneGO(GO)))
    expect_is(GOgraph, class(pruneGO(GOgraph)))
    expect_false('def' %in% names(pruneGO(GO, metadata = c('id', 'name'))$vertices))
    expect_equal('biological_process', unique(pruneGO(GO, namespace = 'bio')$vertices$namespace))
    expect_true(all(sapply(pruneGO(GO, subset = 'prok')$vertices$subset, function(x) {
        "gosubset_prok" %in% x
    })))
    expect_equal('is_a', unique(pruneGO(GO, relations = 'is_a')$edges$type))
})
