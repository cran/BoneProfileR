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
    method <- isolate(input$Method)
    periodic <- isolate(input$Periodic)
    
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
      
      
      if (method == 1) mt <- "fast"
      if (method == 2) mt <- "accurate"
      if (method == 3) mt <- "fastconvex"
      if (method == 4) mt <- "accurateconvex"
      
      options(mc.cores = 2)
      options(forking = ifelse(.Platform$OS.type == "windows", FALSE, TRUE))
      
      bone <- BP_OpenImage(file=file$datapath, ijtiff = ijtiff)
      name <- attributes(bone)$name
      bone <- BP_DetectBackground(bone=bone, analysis="logistic", show.plot=FALSE)
      bone <- BP_DetectForeground(bone=bone, analysis="logistic", show.plot=FALSE)
      bone <- BP_DetectCenters(bone=bone, analysis="logistic", show.plot=FALSE, method=mt)
      
      bone <- BP_EstimateCompactness(bone, analysis="logistic", 
                                     rotation.angle=rotation, 
                                     center=center, 
                                     cut.angle = angles,
                                     cut.distance = distances,
                                     method = mt, 
                                     show.plot=FALSE)
      
      bone <- BP_FitMLCompactness(bone, analysis="logistic", silent=TRUE, 
                                  fixed.parameters = c(K1=1, K2=1), twosteps=TRUE, 
                                  fitted.parameters = c(P=0.5, S=0.1, Max=0.99, Min=0.01))
      # fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
      # bone <- BP_FitMLCompactness(bone, analysis="logistic", silent=TRUE, 
      #                             fixed.parameters = c(K1=0, K2=0), 
      #                             fitted.parameters = c(fittedpar, Max=0.99, Min=0.01))
      # 
      fittedpar <- BP_GetFittedParameters(bone, analysis="logistic", ML=TRUE, return.all=FALSE)[, "mean"]
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
      
      outAIC <- compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", ML=TRUE, return.all=TRUE), 
                            Flexit=BP_GetFittedParameters(bone, analysis="flexit", ML=TRUE, return.all=TRUE), silent = TRUE)
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
      
      if (periodic) {
        par <- BP_GetFittedParameters(bone, analysis=selected.model, ML=TRUE, return.all=FALSE)[, "mean"]
        options(mc.cores=2)
        bone <- BP_FitMLPeriodicCompactness(bone, analysis=selected.model, control.optim=list(trace=2), 
                                            fitted.parameters=c(par, PSin=0.001, PCos=0.001, 
                                                                SSin=0.001, SCos=0.001, MinSin=0.001, MinCos=0.001, 
                                                                MaxSin=0.001, MaxCos=0.001), replicates.CI=2000)
      }
      
      output$DataOut <- renderTable(outAIC, rownames=TRUE, colnames = TRUE)
      output$TitleOut1 <- renderText("<h2><center>Analysis results</center></h2>")
      output$TitleOut2 <- renderText("<h3>Logistic and flexit model selection</h3>")
      
      
      compactness.synthesis <- RM_get(x=bone, RMname=selected.model, valuename = "compactness.synthesis")
      m <- compactness.synthesis$mineralized
      nm <- compactness.synthesis$unmineralize
      outmcmc <- RM_get(x=bone, RMname=selected.model, valuename = "mcmc")
      
      output$TitleOut3  <- renderText("<h3>Global compacity</h3>")
      output$ResultOut2 <- renderText(paste0("<b>Observed compacity: </b>", 
                                             specify_decimal(RM_get(x=bone, RMname = selected.model, valuename = "global.compactness"), 3), 
                                             "<p><b>Modeled compacity by MCMC (2.5%, 50%, 97.5%): </b>", 
                                             specify_decimal(sum(((m+nm)*outmcmc$quantiles["2.5%", ]))/sum((m+nm)), 3), ", ", 
                                             specify_decimal(sum(((m+nm)*outmcmc$quantiles["50%", ]))/sum((m+nm)), 3), ", ", 
                                             specify_decimal(sum(((m+nm)*outmcmc$quantiles["97.5%", ]))/sum((m+nm)), 3)
      ))
      
      
      output$GlobalOut <- renderTable(outmcmc$summary.table, 
                                      rownames=TRUE, colnames = TRUE, digits = 3)
      
      # nenv <- new.env()
      # assign("bone", bone, envir = nenv)
      # 
      # output$PlotModel <- renderPlot({
      #   plot(bone, analysis = selected.model,type="observations+model", CI = "MCMC")
      # }, env = nenv)
      
      bone <<- bone # Pourquoi ? Reste d'un vieux debug ? 20/4/2025. Non utile
      # selected.model <<- selected.model
      
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
      
      tbl1 <- cbind(data.frame(Angle=paste0("]", specify_decimal(out1$angles-delta, 3), ";", 
                                            specify_decimal(out1$angles+delta, 3), "]")), 
                    out1$synthesis[, -1])
      
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
        if (isolate(input$RadialOC)) v <- c(v, "observed.compactness")
        if (isolate(input$RadialLOC)) v <- c(v, "linearized.observed.compactness")
        if (isolate(input$RadialMC)) v <- c(v, "modeled.compactness")
        if (isolate(input$RadialLMC)) v <- c(v, "linearized.modeled.compactness")
        
        plot(bone, analysis = selected.model, 
             type="radial", parameter.name = v)
      })

      if (periodic) {
        output$PlotModelPeriodic <- renderPlot({
          plot(bone, analysis = selected.model,
               type="periodic", parameter.name="compactness", col=rainbow(128))
        })
      }
      
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
      
      
      outAIC <- compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", ML=TRUE, return.all=TRUE), 
                            Flexit=BP_GetFittedParameters(bone, analysis="flexit", ML=TRUE, return.all=TRUE), silent = TRUE)
      if (outAIC$DeltaAIC[1]==0) {
        # Model Logistic
        selected.model <- "logistic"
      } else {
        # Model Flexit
        selected.model <- "flexit"
      }
      
      
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
      
      openxlsx::addWorksheet(
        wb=wb,
        sheetName="Periodic")
      
      
      
      analysis <- selected.model
      
      
      compactness.synthesis <- RM_get(x=bone, RMname=analysis, valuename = "compactness.synthesis")
      m <- compactness.synthesis$mineralized
      nm <- compactness.synthesis$unmineralize
      
      
      out1 <- RM_get(x=bone, RMname=analysis, valuename = "optim")
      
      out <- RM_list(x=bone, silent=TRUE)
      date <- out[[analysis]]$timestamp
      
      if (!is.null(out1)) {
        
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
          x=analysis,
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
          x=RM_get(x=bone, RMname = analysis, valuename = "global.compactness"),
          startCol = 2,
          startRow = 4)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x="ML 2.5%",
          startCol = 2,
          startRow = 5)
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x="ML 50%",
          startCol = 3,
          startRow = 5)
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x="ML 95%",
          startCol = 4,
          startRow = 5)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x="Modeled compactness",
          startCol = 1,
          startRow = 6)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x=sum(((m+nm)*out1$quantiles["2.5%", ]))/sum((m+nm)),
          startCol = 2,
          startRow = 6)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x=sum(((m+nm)*out1$quantiles["50%", ]))/sum((m+nm)),
          startCol = 3,
          startRow = 6)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x=sum(((m+nm)*out1$quantiles["97.5%", ]))/sum((m+nm)),
          startCol = 4,
          startRow = 6)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x="Linear compactnes",
          startCol = 1,
          startRow = 7)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x=mean(out1$quantiles["2.5%", ]),
          startCol = 2,
          startRow = 7)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x=mean(out1$quantiles["50%", ]),
          startCol = 3,
          startRow = 7)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x=mean(out1$quantiles["97.5%", ]),
          startCol = 4,
          startRow = 7)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x="The objective of estimation of modeled compacity is to verify that the fitted model represents well the observed compacity.",
          startCol = 1,
          startRow = 8)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x="The linear compacity represents the integration under the compacity curve then without taking into account that surface at the center is smaller than surface at periphery. The biological meaning of the linear compacity is not clear.",
          startCol = 1,
          startRow = 9)
        
        
        outmcmc <- RM_get(x=bone, RMname=analysis, valuename = "mcmc")
        
        if (!is.null(outmcmc)) {
          
          openxlsx::writeData(
            wb=wb,
            sheet="Global",
            x="MCMC 2.5%",
            startCol = 5,
            startRow = 5)
          openxlsx::writeData(
            wb=wb,
            sheet="Global",
            x="MCMC 50%",
            startCol = 6,
            startRow = 5)
          openxlsx::writeData(
            wb=wb,
            sheet="Global",
            x="MCMC 95%",
            startCol = 7,
            startRow = 5)
          
          
          # openxlsx::writeData(
          #   wb=wb,
          #   sheet="Global",
          #   x="Modeled compactness by MCMC (2.5%, 50%, 97.5%)",
          #   startCol = 1,
          #   startRow = 6)
          
          openxlsx::writeData(
            wb=wb,
            sheet="Global",
            x=sum(((m+nm)*outmcmc$quantiles["2.5%", ]))/sum((m+nm)),
            startCol = 5,
            startRow = 6)
          
          openxlsx::writeData(
            wb=wb,
            sheet="Global",
            x=sum(((m+nm)*outmcmc$quantiles["50%", ]))/sum((m+nm)),
            startCol = 6,
            startRow = 6)
          
          openxlsx::writeData(
            wb=wb,
            sheet="Global",
            x=sum(((m+nm)*outmcmc$quantiles["97.5%", ]))/sum((m+nm)),
            startCol = 7,
            startRow = 6)
          
          # openxlsx::writeData(
          #   wb=wb,
          #   sheet="Global",
          #   x="Modeled corrected compactness by MCMC (2.5%, 50%, 97.5%)",
          #   startCol = 1,
          #   startRow = 8)
          
          openxlsx::writeData(
            wb=wb,
            sheet="Global",
            x=mean(outmcmc$quantiles["2.5%", ]),
            startCol = 5,
            startRow = 7)
          
          openxlsx::writeData(
            wb=wb,
            sheet="Global",
            x=mean(outmcmc$quantiles["50%", ]),
            startCol = 6,
            startRow = 7)
          
          openxlsx::writeData(
            wb=wb,
            sheet="Global",
            x=mean(outmcmc$quantiles["97.5%", ]),
            startCol = 7,
            startRow = 7)
        }
        
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
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x="Quantiles",
          startCol = 1,
          startRow = 20)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Global",
          x=t(out1$quantiles),
          startCol = 1,
          startRow = 21)
      }
      
      out1 <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")
      
      if (!is.null(out1)) {
        
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
          x=analysis,
          startCol = 1,
          startRow = 3)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x=out1$summary.table,
          startCol = 2,
          startRow = 4)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x=rownames(out1$summary.table),
          startCol = 1,
          startRow = 5)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x=out1$synthesis,
          startCol = 1,
          startRow = 20)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x="Observed compacteness",
          startCol = 9,
          startRow = 20)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x="Linearized observed compactness",
          startCol = 10,
          startRow = 20)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x="Modeled compactness",
          startCol = 11,
          startRow = 20)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x="Linearized modeled compactness",
          startCol = 12,
          startRow = 20)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x="The 'Observed compactness' is the observed compactness of the portion.",
          startCol = 1,
          startRow = 12)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x="The 'Linearized observed compactness' is the observed compactness of the portion weighted to be linearized.",
          startCol = 1,
          startRow = 13)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x="The 'Modeled compactness' is the modeled compactness.",
          startCol = 1,
          startRow = 14)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Radial",
          x="The 'Linearized modeled compactness' is the modeled compactness of the portion weighted to be linearized.",
          startCol = 1,
          startRow = 15)
      }
      
      out1 <- RM_get(x=bone, RMname=analysis, valuename = "optimPeriodic")
      
      if (!is.null(out1)) {
        
        openxlsx::writeData(
          wb=wb,
          sheet="Periodic",
          x="Periodic model of compactness",
          startCol = 1,
          startRow = 1)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Periodic",
          x=rbind(Mean=out1$par, SE=out1$SE),
          startCol = 2,
          startRow = 3)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Periodic",
          x="Mean",
          startCol = 1,
          startRow = 4)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Periodic",
          x="SE",
          startCol = 1,
          startRow = 5)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Periodic",
          x=out1$GlobalCompactness,
          startCol = 2,
          startRow = 8)
        
        openxlsx::writeData(
          wb=wb,
          sheet="Periodic",
          x=rownames(out1$GlobalCompactness),
          startCol = 1,
          startRow = 9)
        
      }
      
      openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
      
    })
  
  
  
})
