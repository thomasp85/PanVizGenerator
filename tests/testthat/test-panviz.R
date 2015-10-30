context("panviz conversion")

library(FindMyFriends)
outputDir <- tempfile(pattern = 'PanViz')
dir.create(outputDir)

csvFile <- system.file('extdata', 'exampleData.csv', 
                       package = 'PanVizGenerator')
pangenome <- read.csv(csvFile, quote='', stringsAsFactors = FALSE)
name <- pangenome$name
go <- pangenome$go
ec <- pangenome$ec
pangenome <- pangenome[, !names(pangenome) %in% c('name', 'go', 'ec')]
pangenomeFMF <- .loadPgExample(withNeighborhoodSplit = TRUE)
annotation <- readAnnot(system.file('extdata', 
                                    'examplePG', 
                                    'example.annot', 
                                    package = 'FindMyFriends'))
head(annotation)
pangenomeFMFannot <- addGroupInfo(pangenomeFMF, annotation, key = 'name')

test_that("Data formatting works", {
    expect_equal(formatGO(go), formatGO(strsplit(go, '; ')))
    expect_equal(formatGO(go), formatGO(gsub('GO:', '', go)))
    expect_equal(formatEC(ec), formatEC(strsplit(ec, '; ')))
    expect_equal(formatEC(ec), formatEC(gsub('EC:', '', ec)))
    expect_is(formatGO(go), 'list')
    expect_equal(length(go), length(formatGO(go)))
    expect_true(all(sapply(formatGO(go), class) == 'character'))
    expect_is(formatEC(ec), 'list')
    expect_equal(length(ec), length(formatEC(ec)))
    expect_true(all(sapply(formatEC(ec), class) == 'character'))
    expect_error(formatData(pangenome, name, go, ec))
    data <- formatData(as.matrix(pangenome), name, go, ec)
    expect_named(data, c('pangenome', 'genes'))
    expect_is(data$pangenome, 'data.frame')
    expect_equal(dim(data$pangenome), dim(pangenome))
    expect_is(data$genes, 'data.frame')
    expect_named(data$genes, c('name', 'go', 'ec'))
    expect_equal(nrow(data$pangenome), nrow(data$genes))
    expect_equal(sapply(data$genes, class), c(name='character', go='AsIs', ec='AsIs'))
})

unlink(list.files(outputDir, full.names = TRUE))
test_that("data.frame input works", {
    panviz(pangenome, name=name, go=go, ec=ec, location = outputDir, scale = FALSE)
    expect_equal(sort(list.files(outputDir)), c("PanViz.html", "README"))
})

unlink(list.files(outputDir, full.names = TRUE))
test_that("matrix input works", {
    panviz(as.matrix(pangenome), name=name, go=go, ec=ec, location = outputDir, scale = FALSE)
    expect_equal(sort(list.files(outputDir)), c("PanViz.html", "README"))
})

unlink(list.files(outputDir, full.names = TRUE))
test_that("csv input works", {
    panviz(csvFile, location = outputDir, scale = FALSE)
    expect_equal(sort(list.files(outputDir)), c("PanViz.html", "README"))
})

unlink(list.files(outputDir, full.names = TRUE))
test_that("FMF input works", {
    expect_error(panviz(pangenomeFMF, location = outputDir))
    panviz(pangenomeFMFannot, location = outputDir, scale = FALSE)
    expect_equal(sort(list.files(outputDir)), c("PanViz.html", "README"))
})

unlink(list.files(outputDir, full.names = TRUE))
test_that("File consolidation works", {
    panviz(as.matrix(pangenome), name=name, go=go, ec=ec, location = outputDir, consolidate = FALSE, scale = FALSE)
    expect_true(all(list.files(outputDir) %in% c("data.js", "PanViz.html", "README", "script.js", "style.css")))
    consolidateFiles(outputDir)
    expect_equal(sort(list.files(outputDir)), c("PanViz.html", "README"))
    unlink(list.files(outputDir, full.names = TRUE))
    expect_error(consolidateFiles(outputDir))
})
