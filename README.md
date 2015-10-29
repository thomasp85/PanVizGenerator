# PanVizGenerator

[![Build Status](https://travis-ci.org/thomasp85/PanVizGenerator.svg?branch=master)](https://travis-ci.org/thomasp85/PanVizGenerator)
[![codecov.io](https://codecov.io/github/thomasp85/PanVizGenerator/coverage.svg?branch=master)](https://codecov.io/github/thomasp85/PanVizGenerator?branch=master)

This R package is a companion to the 
[PanViz](https://github.com/thomasp85/PanViz) javascript visualization for 
functionally annotated pangenomes. While PanViz is fully self-contained once it 
is created, this package takes care of converting your pangenome data into a 
PanViz file.

PanVizGenerator supports generic pangenome notation in the form of 
presence/absence matrices (also known as pangenome matrices) as well as a direct
link to the [FindMyFriends](http://bioconductor.org/packages/FindMyFriends/) 
class system to allow for direct creation of PanViz from your FindMyFriends 
analysis.

## Usage
There are two general approaches to generating PanViz visualizations with
PanVizGenerator. Either using the `panviz()` method or by summoning a shiny app
using the `PanVizGenerator()` function.

Simple example with csv input:

```r
csvFile <- system.file('extdata', 'exampleData.csv', 
                       package = 'PanVizGenerator')
outputDir <- tempdir()

# Generate the visualization
panviz(csvFile, location = outputDir)
```

Or using FindMyFriends objects:

```r
library(FindMyFriends)
# Load an example pangenome
pangenome <- .loadPgExample(withNeighborhoodSplit = TRUE)

# Add functional annotation from a Blast2GO analysis
annotation <- readAnnot(system.file('extdata', 
                                    'examplePG', 
                                    'example.annot', 
                                    package = 'FindMyFriends'))
pangenome <- addGroupInfo(pangenome, annotation, key = 'name')

# Generate the visualization
panviz(pangenome, location = outputDir)
```

## Installation
This package is intended for future inclusion into the Bioconductor project. 
Until then it can be installed using devtools:

```r
if (!require(devtools)) {
    install.packages('devtools')
    library(devtools)
}
install_github('thomasp85/PanVizGenerator')
```
