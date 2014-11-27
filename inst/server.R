library(PanVizGenerator)
library(shiny)

shinyServer(function(input, output, session) {
    options(shiny.maxRequestSize=100*1024^2)
    panviz <- reactiveValues(path=NA)
    observe({
        #browser()
        if(input$generateBtn == 0) return(NULL)
        
        file <- isolate({input$panSelect})
        
        if(is.null(file)) return(NULL)
        
        data <- tryCatch(
            expr = PanVizGenerator:::parsePangenome(file),
            error = function(e) {
                session$sendCustomMessage('toggleError', list(state=TRUE, message=e$message))
                return(NA)
            }
        )
        if(class(data) != 'list') return(NULL)
        
        tempFolder <- tempfile('pan')
        panFolder <- file.path(tempFolder, 'PanViz')
        success <- tryCatch(
            expr = {
                dir.create(tempFolder)
                dir.create(panFolder)
                TRUE
            },
            error = function(e) {
                session$sendCustomMessage('toggleError', list(state=TRUE, message='Unable to create working folder'))
                return(FALSE)
            }
        )
        if(!success) return(NULL)
        
        success <- tryCatch(
            expr = {
                PanVizGenerator:::createPanData(
                    data, 
                    panFolder, 
                    dist=isolate(input$distance), 
                    clust=isolate(input$clustering), 
                    center=isolate(input$center), 
                    scale=isolate(input$scale)
                )
                TRUE
            },
            error = function(e) {
                session$sendCustomMessage('toggleError', list(state=TRUE, message=e$message))
                return(FALSE)
            }
        )
        if(!success) return(NULL)
        
        oldWarn <- options(warn=2)$warn
        success <- tryCatch(
            expr = {
                files <- list.files('PanViz/', full.names=TRUE)
                file.copy(from=files, to=panFolder, overwrite=TRUE)
                TRUE
            },
            error = function(e) {
                session$sendCustomMessage('toggleError', list(state=TRUE, message='Could not copy PanViz files'))
                return(FALSE)
            }
        )
        options(warn=oldWarn)
        if(!success) return(NULL)
        
        panviz$path <- tempFolder
        session$sendCustomMessage('toggleSuccess', TRUE)
    })
    
    output$download <- downloadHandler(
        filename = 'PanViz.zip',
        content = function(file) {
            oldWd <- setwd(isolate(panviz$path))
            on.exit(setwd(oldWd))
            zip(file, 'PanViz/')
            if (file.exists(paste0(file, ".zip"))) {
                file.rename(paste0(file, ".zip"), file)
            }
        },
        contentType = "application/zip"
    )
    outputOptions(output, 'download', suspendWhenHidden = FALSE)
})