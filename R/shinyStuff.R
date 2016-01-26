#' @include aaa.R
#' @include panviz.R
NULL

#' Launch the PanVizGenerator shiny app
#' 
#' This function launches a shiny based GUI that lets the user create panviz 
#' visualizations from csv files. The same formatting restriction as in 
#' \code{\link{panviz}} applies to the csv file. Furthermore the GUI also 
#' contains descriptions of the ideas behind the visualization as well as a 
#' demonstration video.
#' 
#' @param ... Parameters passed on to \code{\link[shiny]{runApp}}
#' 
#' @return This function is called for its side effects.
#' 
#' @seealso \code{\link{panviz}} for an R API for the same functionality.
#' 
#' @examples 
#' if (interactive()) {
#'     PanVizGenerator()
#' }
#' 
#' @importFrom shiny runApp
#' 
#' @export
#' 
PanVizGenerator <- function(...) {
    runApp(system.file(package = 'PanVizGenerator'), ...)
}

# INTERNALS
#' Server logic for the PanVizGenerator shiny app
#' 
#' @importFrom shiny isolate observe downloadHandler outputOptions
#' @importFrom utils zip
#' 
#' @noRd
#' 
server <- function(input, output, session) {
    options(shiny.maxRequestSize = 100*1024 ^ 2)
    globalVars <- new.env(parent = emptyenv())
    observe({
        if (input$generateBtn == 0) return(NULL)
        file <- isolate({input$panSelect})
        if (is.null(file)) return(NULL)
        
        tempFolder <- tempfile('pan')
        panFolder <- file.path(tempFolder, 'PanViz')
        workFolder <- file.path(tempFolder, 'wd')
        fileName <- file.path(workFolder, file$name)
        success <- tryCatch(
            expr = {
                dir.create(tempFolder)
                dir.create(panFolder)
                dir.create(workFolder)
                file.copy(file$datapath, fileName)
                TRUE
            },
            error = function(e) {
                session$sendCustomMessage(
                    'toggleError', 
                    list(state = TRUE, message = 'Unable to create work folder')
                )
                return(FALSE)
            }
        )
        if (!success) return(NULL)
        
        success <- tryCatch(
            expr = {
                panviz(fileName, 
                       location = panFolder, 
                       dist = isolate(input$distance), 
                       clust = isolate(input$clustering), 
                       center = isolate(input$center), 
                       scale = isolate(input$scale))
                TRUE
            },
            error = function(e) {
                session$sendCustomMessage(
                    'toggleError', 
                    list(state = TRUE, message = e$message)
                )
                return(FALSE)
            }
        )
        if (!success) return(NULL)
        
        assign('path', tempFolder, envir = globalVars)
        session$sendCustomMessage('toggleSuccess', TRUE)
    })
    
    output$exampleDownload <- downloadHandler(
        filename = 'Pangenome.csv',
        content = function(file) {
            exampleFile <- system.file(
                'extdata', 
                'exampleData.csv', 
                package = 'PanVizGenerator'
            )
            print(exampleFile)
            print(file)
            file.copy(from = exampleFile, to = file)
        }
    )
    outputOptions(output, 'exampleDownload', suspendWhenHidden = FALSE)
    output$download <- downloadHandler(
        filename = 'PanViz.zip',
        content = function(file) {
            oldWd <- setwd(globalVars$path)
            on.exit(setwd(oldWd))
            zip(file, 'PanViz/')
            if (file.exists(paste0(file, ".zip"))) {
                file.rename(paste0(file, ".zip"), file)
            }
        },
        contentType = "application/zip"
    )
    outputOptions(output, 'download', suspendWhenHidden = FALSE)
}