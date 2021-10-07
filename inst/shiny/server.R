library(shiny)


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  
  
  
  
  outplotfn <- eventReactive(eventExpr=input$goButton, ignoreNULL = FALSE, valueExpr={
    
    file <- isolate(input$FileOpen)
    center <- isolate(input$center)
    rotation <- isolate(as.numeric(input$rotation))
    ijtiff <- isolate(input$ijtiff)
    angles <- isolate(as.numeric(input$angles))
    distances <- isolate(as.numeric(input$distances))
    twosteps <- isolate(input$twosteps)
    
    #' file <- list(datapath="/Users/marcgirondot/Desktop/ProblemesBP/astrochelys_femur.tif")
    #' file <- list(datapath="/Users/marcgirondot/Downloads/femur45l-astrochelys.tif")
    #' center <- "ontogenic"
    #' rotation <- 0
    #' ijtiff <- TRUE
    #' angles <- 60
    #' distances <- 100
    #' twosteps <- TRUE
    
    if (is.null(file)) {
      
      oldpar <- par(no.readonly = TRUE)    # code line i
      on.exit(par(oldpar))            # code line i + 1 
      
      
      par(mar=c(0, 0, 0, 0))
      plot(x=c(0, 1), y=c(0, 1), axes=FALSE,
           xaxt="n", yaxt="n", main="",
           xlab = "", ylab = "",
           xaxs="i", yaxs="i", type="n")
      text(x = 0.5, y=0.6, labels = "Load an image to be analyzed",
           col="red", cex = 1.6)
    } else {
      # specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
      
      
      bone <- BP_OpenImage(file=file$datapath, ijtiff = ijtiff)
      name <- attributes(bone)$name
      bone <- BP_DetectBackground(bone=bone, analysis="logistic", show.plot=FALSE)
      bone <- BP_DetectForeground(bone=bone, analysis="logistic", show.plot=FALSE)
      bone <- BP_DetectCenters(bone=bone, analysis="logistic", show.plot=FALSE)
      bone <- BP_EstimateCompactness(bone, analysis="logistic", 
                                     rotation.angle=rotation, 
                                     center=center, 
                                     cut.angle = angles,
                                     cut.distance = distances,
                                     show.plot=FALSE)

      bone <- BP_FitMLCompactness(bone, analysis="logistic", silent=TRUE, 
                                  fixed.parameters = c(K1=1, K2=1), twosteps=TRUE, 
                                  fitted.parameters = c(P=0.5, S=0.1, Max=0.99, Min=0.01))
      # fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
      # bone <- BP_FitMLCompactness(bone, analysis="logistic", silent=TRUE, 
      #                             fixed.parameters = c(K1=0, K2=0), 
      #                             fitted.parameters = c(fittedpar, Max=0.99, Min=0.01))
      # 
      fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
      bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
      bone <- BP_FitMLCompactness(bone, 
                                  fitted.parameters=c(fittedpar, K1=0.9, K2=0.9), 
                                  fixed.parameters=NULL, analysis="flexit", silent=TRUE, twosteps=TRUE)
      # bone <- BP_FitBayesianCompactness(bone, analysis="flexit")
      # mcmc <- RM_get(bone, RMname = "flexit", value="mcmc")
      # fittedpar <- as.parameters(mcmc)
      # bone <- BP_FitMLCompactness(bone, 
      #                             fitted.parameters=fittedpar, 
      #                             fixed.parameters=NULL, analysis="flexit", silent=TRUE)
      
      outAIC <- compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", alloptim=TRUE), 
                            Flexit=BP_GetFittedParameters(bone, analysis="flexit", alloptim=TRUE), silent = TRUE)
      if (outAIC$DeltaAIC[1]==0) {
        # Model Logistic
        selected.model <- "logistic"
        output$ResultOut1 <- renderText(paste0("The selected model based on AIC is the logistic model. The 
                                        probability that the logistic model is the best among the two 
                                        tested is ", specify_decimal(outAIC[1, 3], 3), "."))
      } else {
        # Model Flexit
        selected.model <- "flexit"
        output$ResultOut1 <- renderText(paste0("The selected model based on AIC is the flexit model. The 
                                        probability that the flexit model is the best among the two 
                                        tested is ", specify_decimal(outAIC[2, 3], 3), "."))
      }
      bone <- BP_FitBayesianCompactness(bone, analysis=selected.model)
      bone <- BP_FitMLRadialCompactness(bone, analysis=selected.model, silent=TRUE, twosteps=twosteps)
      
      
      output$DataOut <- renderTable(outAIC, rownames=TRUE, colnames = TRUE)
      output$TitleOut1 <- renderText("<h2><center>Analysis results</center></h2>")
      output$TitleOut2 <- renderText("<h3>Logistic and flexit model selection</h3>")
      
      output$TitleOut3  <- renderText("<h3>Global compacity</h3>")
      output$ResultOut2 <- renderText(paste0("<b>Observed compacity: </b>", 
                                             specify_decimal(RM_get(x=bone, RMname = selected.model, valuename = "global.compactness"), 3), 
                                             "<p><b>Modeled compacity by MCMC (2.5%, 50%, 97.5%): </b>", 
                                             specify_decimal(mean(RM_get(x=bone, RMname = selected.model, valuename = "mcmc")$quantiles["2.5%", ]), 3), ", ", 
                                             specify_decimal(mean(RM_get(x=bone, RMname = selected.model, valuename = "mcmc")$quantiles["50%", ]), 3), ", ", 
                                             specify_decimal(mean(RM_get(x=bone, RMname = selected.model, valuename = "mcmc")$quantiles["97.5%", ]), 3)
      ))
      
      output$GlobalOut <- renderTable(RM_get(x=bone, RMname=selected.model, valuename = "mcmc")$summary.table, 
                                      rownames=TRUE, colnames = TRUE, digits = 3)
      
      # nenv <- new.env()
      # assign("bone", bone, envir = nenv)
      # 
      # output$PlotModel <- renderPlot({
      #   plot(bone, analysis = selected.model,type="observations+model", CI = "MCMC")
      # }, env = nenv)
      
      bone <<- bone
      selected.model <<- selected.model
      
      output$PlotModel <- renderPlot({
        plot(bone, analysis = selected.model,type="observations+model", CI = "MCMC")
      })
      
      
      output$TitleOut4  <- renderText("<h3>Radial compacity</h3>")
      
      out1 <- RM_get(x=bone, RMname=selected.model, valuename = "optimRadial")
      ost <- out1$summary.table
      colnames(ost)[2] <- "sd"
      output$RadialOut1 <- renderTable(ost, 
                                       rownames=TRUE, colnames = TRUE, digits = 3)
      
      delta <- (out1$angles[2]-out1$angles[1])/2
      
      tbl1 <- cbind(data.frame(Angle=paste0(specify_decimal(out1$angles-delta, 3), ";", 
                                            specify_decimal(out1$angles+delta, 3))), 
                    out1$synthesis)
      tbl1 <- cbind(tbl1, data.frame(Modeled=out1$radial.modeled.compactness, 
                                     Observed=out1$observed.compactness, 
                                     Observed.modeled=out1$observed.modeled.compactness))
      
      tbl1 <<- tbl1
      
      output$RadialOut2 <- renderTable(tbl1, rownames=FALSE, 
                                       colnames = TRUE, hover = TRUE, digits=3)
      
      output$PlotModelRadial <- renderPlot({
        v <- NULL
        if (isolate(input$RadialVarMin)) v <- c(v, "Min")
        if (isolate(input$RadialVarMax)) v <- c(v, "Max")
        if (isolate(input$RadialVarP)) v <- c(v, "P")
        if (isolate(input$RadialVarS)) v <- c(v, "S")
        if (isolate(input$RadialVarK1)) v <- c(v, "K1")
        if (isolate(input$RadialVarK2)) v <- c(v, "K1")
        if (isolate(input$RadialVarTRC)) v <- c(v, "TRC")
        
        plot(bone, analysis = selected.model, 
             type="radial", radial.variable = v)
      })
      
      plot(bone, analysis = selected.model, show.grid=TRUE)
      
    }
  }
  )
  
  
  
  
  output$Plot <- renderPlot({
    
    outplotfn()
    
  })
  
  
  output$ExcelButton <- downloadHandler(
    
    filename = "Export.xlsx",
    content = function(file) {
      
      # writeLines(paste0(c("bonjour", "Hello"), collapse = "\n"), con=file)
      
      wb <- openxlsx::createWorkbook(creator = "author"
                                     , title = "title"
                                     , subject = "BoneProfileR report"
                                     , category = "")
      
      openxlsx::addWorksheet(
        wb=wb,
        sheetName="Global")
      
      openxlsx::addWorksheet(
        wb=wb,
        sheetName="Radial")
      
      out1 <- RM_get(x=bone, RMname=selected.model, valuename = "mcmc")
      
      out <- RM_list(x=bone, silent=TRUE)
      date <- out[[selected.model]]$timestamp
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x="Global model of compactness",
        startCol = 1,
        startRow = 1)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=date,
        startCol = 1,
        startRow = 2)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=selected.model,
        startCol = 1,
        startRow = 3)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x="Observed compactness",
        startCol = 1,
        startRow = 4)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=RM_get(x=bone, RMname = selected.model, valuename = "global.compactness"),
        startCol = 2,
        startRow = 4)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x="Modeled compactness by MCMC (2.5%, 50%, 97.5%)",
        startCol = 1,
        startRow = 5)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=mean(RM_get(x=bone, RMname = selected.model, valuename = "mcmc")$quantiles["2.5%", ]),
        startCol = 2,
        startRow = 5)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=mean(RM_get(x=bone, RMname = selected.model, valuename = "mcmc")$quantiles["50%", ]),
        startCol = 3,
        startRow = 5)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=mean(RM_get(x=bone, RMname = selected.model, valuename = "mcmc")$quantiles["97.5%", ]),
        startCol = 4,
        startRow = 5)
      
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=out1$summary.table,
        startCol = 2,
        startRow = 10)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=rownames(out1$summary.table),
        startCol = 1,
        startRow = 11)
      
      # Le modèle radial
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x="Radial model of compactness",
        startCol = 1,
        startRow = 1)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=date,
        startCol = 1,
        startRow = 2)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=selected.model,
        startCol = 1,
        startRow = 3)
      
      
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=RM_get(x=bone, RMname=selected.model, valuename = "optimRadial")$summary.table,
        startCol = 2,
        startRow = 4)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x="SD",
        startCol = 3,
        startRow = 4)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=rownames(RM_get(x=bone, RMname=selected.model, valuename = "optimRadial")$summary.table),
        startCol = 1,
        startRow = 5)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=tbl1,
        startCol = 1,
        startRow = 13)
      
      openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
      
    })
  
  
  
})
