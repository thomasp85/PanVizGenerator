context("Create JSON data")

library('digest')

pg <- system.file('extdata', 'exampleData.csv', package = 'PanVizGenerator')
pg <- read.csv(pg, header = TRUE, check.names = FALSE, 
               comment.char = '', stringsAsFactors = FALSE)
name <- pg$name
go <- pg$go
ec <- pg$ec
pg <- as.matrix(pg[, -(1:3)])
data <- formatData(pg, name, go, ec)

test_that("createPangenome works", {
  expect_equal(digest(createPangenome(data)), "499671f6a853915c06e7b08df3288d9d")
})

test_that("createGeneInfo works", {
    expect_equal(digest(createGeneInfo(data)), "f0c0df28f9a76fa0b0633d6923f467fc")
})

test_that("createScatter works", {
    expect_equal(digest(suppressWarnings(createScatter(data))), "43c144d9f16242c28ee878c4b442fe73")
    expect_equal(digest(suppressWarnings(createScatter(data, dist = 'binary'))), "3eeae5f3b9ce4e6c63faabe990faa7fe")
    expect_error(createScatter(data, dist = 'test'))
    expect_equal(digest(createScatter(data, scale = FALSE)), "0bfd5ea88c1b760d27573d9aa20ccb54")
    expect_equal(suppressWarnings(digest(createScatter(data, center = FALSE))), "182889b54c676ef1c7504c9c43631392")
})

test_that("createCluster works", {
    expect_equal(digest(createCluster(data)), "48ed102e83adc31b2e7ebb49e9be1c91")
    expect_equal(digest(createCluster(data, dist = 'binary')), "bd562222abf288d6481577376735c73a")
    expect_error(createCluster(data, dist = 'test'))
    expect_equal(digest(createCluster(data, clust = 'complete')), "abe93133a59712ea9c9cedb27d996148")
    expect_error(createCluster(data, clust = 'test'))
})

test_that("createGO works", {
    GO <- createGO()
    expect_is(GO, "json")
    expect_equal(length(GO), 1)
    GOrev <- jsonlite::fromJSON(GO)
    expect_named(GOrev, c('edges', 'vertices'), ignore.order = TRUE)
    expect_named(GOrev$edges, c('to', 'from', 'type'), ignore.order = TRUE)
    expect_named(GOrev$vertices, c('id', 'name', 'namespace', 'def', 'is_obsolete', 'alt_id', 'subset'), ignore.order = TRUE)
    expect_true(max(c(GOrev$edges$from, GOrev$edges$to)) <= length(GOrev$vertices$id))
})
