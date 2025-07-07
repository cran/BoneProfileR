#' plot.BoneProfileR displays a bone section
#' @title Plot a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Nothing
#' @param x The bone image
#' @param message The message to be displayed
#' @param type The type of plot; see description
#' @param angle Which angle model to show
#' @param parameter.name The parameter to plot
#' @param options.mcmc The option to plot type mcmc output
#' @param show.all.angles For periodic type and partial section, should all angles been shown?
#' @param show.colors Should the background and foreground colors be shown?
#' @param show.centers Should the centers be shown?
#' @param show.grid Should the grid be shown?
#' @param show.legend Should a legend be shown? 
#' @param analysis Name or number of analysis to be plotted
#' @param CI Which confidence interval should be plotted: MCMC or ML
#' @param replicates.CI How many replicates to estimate CI?
#' @param restorePar If TRUE, restore the par parameter at the exit
#' @param mar The margin for type being "model" or "observations"
#' @param angle.3D The angle between x and y for 3Dcolors graph
#' @param ... Default parameters for some functions
#' @description Display a bone section.\cr
#' type value can be:\cr
#' Image plot: `original`, `mineralized`, `unmineralized`, `section`\cr
#' Original is the original image, mineralized is the mineral interpretation of the section, 
#' unmineralized is the unmineralized interpretation of the section, section is the interpretation of the section.\cr
#' `colors` shows the histograms of pixel information with foreground and background colors if they are defined.\cr
#' `3Dcolors` show the pixels colors in 3D\cr
#' Global analysis: `observations`, `model`, `observations+model`\cr
#' Radial analysis: `radial`\cr
#' Periodic analysis: `periodic`\cr
#' If angle is not null and a radial analysis exists, it will show the model for this angle.\cr
#' `mcmc`: It will show the posterior distribution of parameter.\cr
#' For periodic analysis, you can see a particular parameter with parameter.name being
#' P, S, Min, Max, K1, or K2 or the global median compactness using parameter.name="compactness". 
#' You can use col=rainbow(128) or hcl.colors(128) to see the region of transition. You can 
#' also plot the average compactness using parameter.name="averagemodel".\cr
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#'  bone <- BP_OpenImage()
#'  # or 
#'  path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_DetectBackground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectForeground(bone=bone, analysis="logistic")
#'  plot(bone, type="colors")
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic")
#'  plot(bone, type="3Dcolors")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic")
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic")
#'  plot(bone)
#'  ############################################
#'  # Example with comparison between two models
#'  ############################################
#'  path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_DetectBackground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectForeground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic")
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic")
#'  plot(bone)
#'  plot(bone, type="observations")
#'  plot(bone, type="observations+model", analysis=1)
#'  fittedpar <- BP_GetFittedParameters(bone, analysis="logistic", 
#'                                      ML=TRUE, return.all = FALSE)[, "mean"]
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  bone <- BP_FitMLCompactness(bone, 
#'                 fitted.parameters=c(fittedpar, K1=1, K2=1), 
#'                 fixed.parameters=NULL, analysis="flexit")
#'  compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", 
#'                                              ML=TRUE, return.all = TRUE), 
#'              Flexit=BP_GetFittedParameters(bone, analysis="flexit", 
#'                                            ML=TRUE, return.all = TRUE))
#'  out4p <- plot(bone, type="observations+model", analysis="logistic")
#'  out6p <- plot(bone, type="observations+model", analysis="flexit")
#'  
#'  ############################################
#'  # Fit distribution using Bayesian model
#'  ############################################
#'  bone <- BP_FitBayesianCompactness(bone, analysis="logistic", n.adapt=100)
#'  # Test the output - New in version 3.2
#'  plot(bone, type="mcmc", options.mcmc = list(what="LnL"))
#'  #########################################################################
#'  # Clearly the distribution is not stationary; the adaptation was too short
#'  ######################################################################### 
#'  bone <- BP_FitBayesianCompactness(bone, analysis="logistic", n.adapt=10000)
#'  # Now it is ok
#'  plot(bone, type="mcmc", options.mcmc = list(what="LnL"))
#'  #########################################################################
#'  # New in version 3.2
#'  #########################################################################
#'  plot(bone, type="mcmc", options.mcmc = list(what="Posterior", 
#'       xlim=c(0.025, 0.035), breaks=seq(from=0.025, to=0.035, by=0.001)), 
#'       parameter.name = "S")
#'  plot(bone, type="mcmc", options.mcmc = list(what="MarkovChain", 
#'                                             ylim=c(0.02, 0.04)), 
#'                                             parameter.name = "S")
#'  #########################################################################
#'  # Check the priors and the output
#'  #########################################################################
#'  mcmc <- RM_get(x=bone, RMname="logistic", valuename = "mcmc")
#'  priors <- mcmc$parametersMCMC$parameters
#'  parameters <- as.parameters(mcmc, index="median")
#'  #########################################################################
#'  # Now it is ok. It can be used
#'  #########################################################################
#'  plot(bone, type="observations+model", CI="MCMC")
#'  plot(bone, type="observations+model", CI="ML")
#'  #########################################################################
#'  
#'  #############################################
#'  # Radial compactness
#'  #############################################
#'  bone <- BP_FitMLRadialCompactness(bone, progressbar=TRUE)
#'  plot(bone, type="radial", parameter.name=c("P", "S"))
#'  plot(bone, type="radial", parameter.name=c("P", "S", "Min", "Max"))
#'  plot(bone, type="radial", parameter.name="observed.compactness")
#'  plot(bone, type="radial", parameter.name="linearized.observed.compactness")
#'  
#'  #############################################
#'  # Periodic analysis
#'  # This model can take 10 minutes to be fitted
#'  # And still more if you use large replicates.CI value
#'  #############################################
#'  bone <- BP_FitMLPeriodicCompactness(bone, analysis="logistic", control.optim=list(trace=2), 
#'                                      fitted.parameters=c(par, PSin=0.001, PCos=0.001, 
#'                                      SSin=0.001, SCos=0.001, MinSin=0.001, MinCos=0.001, 
#'                                      MaxSin=0.001, MaxCos=0.001), replicates.CI=2000)
#'  plot(bone, type="periodic", parameter.name="compactness", col=rainbow(128))
#'  plot(bone, type="periodic", parameter.name="compactness")
#'  plot(bone, type="periodic", parameter.name="P", ylim=c(0, 1), 
#'        col=rgb(red = 0.7, green = 0.7, blue = 0.7, alpha = 0.2))
#'  plot(bone, type="periodic", parameter.name="averagemodel")
#'  
#' }
#' @method plot BoneProfileR
#' @export


plot.BoneProfileR <- function(x                     , 
                              message=NULL          , 
                              type="original"       , 
                              angle=NULL            , 
                              show.all.angles=FALSE ,
                              show.centers=TRUE     , 
                              show.colors=TRUE      , 
                              show.grid=TRUE        , 
                              analysis=1            , 
                              parameter.name = "S"  , 
                              options.mcmc = list() , 
                              restorePar = TRUE     , 
                              mar=NULL              , 
                              angle.3D = 55         , 
                              CI="ML"               , 
                              replicates.CI=1000    , 
                              show.legend=TRUE      , 
                              ...                   ) {
  
  # message=NULL; type="original"; angle=NULL; parameter.name = "S"; options.mcmc = list(); restorePar=TRUE; mar=NULL; show.centers=TRUE; show.colors=TRUE; show.grid=TRUE; analysis=1; CI="ML"; show.legend=TRUE
  
  # type <- "observations"
  # type <- "radial"
  # type <- "model"
  # type <- "colors"
  
  if (is.null(analysis)) {
    stop("You must choose which analysis to report.")
  }
  
  if (type != "original") {
    out <- RM_list(x=x, silent=TRUE)
    if (is.numeric(analysis))
      if (analysis > length(out)) {
        stop("The analysis does no exist.")
      } else {
        analysis <- names(out)[analysis]
      }
    
    if (is.character(analysis) & (all(analysis != names(out)))) {
      stop(paste("The analysis", analysis, "does not exit. Check your data."))
    }
  }
  
  p3p <- tryCatch(list(...), error=function(e) list()) # p3p <- list()
  
  oldpar <- par(no.readonly = TRUE)    # code line i
  if (restorePar) on.exit(par(oldpar))            # code line i + 1
  
  type <- tolower(type)
  
  type <- match.arg(type, choices = c("original", "mineralized", 
                                      "unmineralized", "section", "radial", 
                                      "periodic", "mcmcperiodic", 
                                      "observations", "model", "observations+model", 
                                      "mcmc", "colors", "3dcolors"))
  
  out <- BP_ListAnalyses(bone=x)
  if (is.null(out[analysis][[1]]) & (type != "original") & (type != "colors")) {
    stop(paste0("The analysis ", analysis, " does not exist"))
  }
  
  if ((type == "periodic") | (type == "mcmcperiodic")) {
    
    parameter.name <- tolower(parameter.name)
    parameter.name <- match.arg(parameter.name, choices = c("p", "s", "min", "max", "k1", "k2", "compactness", "averagemodel"))
    
    if (type == "periodic") {
      o <- RM_get(x, RMname = analysis, valuename="optimPeriodic")
    } else {
      o <- RM_get(x, RMname = analysis, valuename="mcmcPeriodic")
    }
    
    if (is.null(o)) {
      stop("The periodic analysis has not been done !")
    } else {
      
      if (parameter.name %in% c("p", "s", "min", "max", "k1", "k2")) {
        fixed.parameters <- o$fixed.parameters
        par <- c(o$par, fixed.parameters)
        if (all(names(par) != "Min")) fixed.parameters <- c(fixed.parameters, Min = 1E-10)
        if (all(names(par) != "Max")) fixed.parameters <- c(fixed.parameters, Max = 1-1E-10)
        if (all(names(par) != "K1")) fixed.parameters <- c(fixed.parameters, K1=1)
        if (all(names(par) != "K2")) fixed.parameters <- c(fixed.parameters, K2=1)
        if (all(names(par) != "PSin")) fixed.parameters <- c(fixed.parameters, PSin=0)
        if (all(names(par) != "PCos")) fixed.parameters <- c(fixed.parameters, PCos=0)
        if (all(names(par) != "SSin")) fixed.parameters <- c(fixed.parameters, SSin=0)
        if (all(names(par) != "SCos")) fixed.parameters <- c(fixed.parameters, SCos=0)
        if (all(names(par) != "MinSin")) fixed.parameters <- c(fixed.parameters, MinSin=0)
        if (all(names(par) != "MinCos")) fixed.parameters <- c(fixed.parameters, MinCos=0)
        if (all(names(par) != "MaxSin")) fixed.parameters <- c(fixed.parameters, MaxSin=0)
        if (all(names(par) != "MaxCos")) fixed.parameters <- c(fixed.parameters, MaxCos=0)
        if (all(names(par) != "K1Sin")) fixed.parameters <- c(fixed.parameters, K1Sin=0)
        if (all(names(par) != "K1Cos")) fixed.parameters <- c(fixed.parameters, K1Cos=0)
        if (all(names(par) != "K2Sin")) fixed.parameters <- c(fixed.parameters, K2Sin=0)
        if (all(names(par) != "K2Cos")) fixed.parameters <- c(fixed.parameters, K2Cos=0)
        if (type == "periodic") {
          rd <- RandomFromHessianOrMCMC(method = "hessian", 
                                        Hessian = o$hessian, 
                                        fitted.parameters = o$par, 
                                        fixed.parameters = fixed.parameters,
                                        replicates = replicates.CI, 
                                        probs = NULL, silent = TRUE)$random
        } else {
          rd <- RandomFromHessianOrMCMC(method = "mcmc", 
                                        mcmc = o, 
                                        fitted.parameters = o$par, 
                                        fixed.parameters = fixed.parameters,
                                        replicates = replicates.CI, 
                                        probs = NULL, silent = TRUE)$random
        }
        
        ag <- RM_get(x=x, RMname = analysis, valuename = "cut.angle")
        # ag <- RM_get(x=x, RMname = analysis, valuename = "peripherie")$angle.center
        sinag <- sin(ag)
        cosag <- cos(ag)
      }
      
      if (any(parameter.name =="p")) {
        es <- sapply(1:replicates.CI, FUN = function(i) (rd[i, "P"]+ rd[i, "PSin"]*sinag + rd[i, "PCos"]*cosag))
        
        ppx <- modifyList(modifyList(list(xlim=c(min(ag), max(ag)), ylim=c(min(es), max(es)), las=1, bty="n", 
                                          ylab="P parameter", xlab="Angle in radians"), p3p[names(p3p) %in% c("ylim", "las", "bty", "xlab", "ylab")]), 
                          list(x=ag, y=rep(0, length(ag)), type="n"))
        do.call("plot", args = ppx)
        for (i in 1:replicates.CI) {
          ppx <- modifyList(modifyList(list(lwd=0.1, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.2)), 
                                       p3p[names(p3p) %in% c("lwd", "col", "xlab", "ylab")]), 
                            list(x=ag, y=es[, i]))
          do.call("lines", args = ppx)
        }
        qrd <- sapply(X = 1:nrow(es), FUN = function(i) quantile(es[i, ], probs=c(0.025, 0.5, 0.975)))
        lines(ag, qrd["50%", ], lty=2, lwd=0.5, col="black")
        lines(ag, qrd["2.5%", ], lty=3, lwd=0.5, col="black")
        lines(ag, qrd["97.5%", ], lty=3, lwd=0.5, col="black")
      }
      
      if (any(parameter.name =="s")) {
        es <- sapply(1:replicates.CI, FUN = function(i) (rd[i, "S"]+ rd[i, "SSin"]*sinag + rd[i, "SCos"]*cosag))
        ppx <- modifyList(modifyList(list(xlim=c(min(ag), max(ag)), ylim=c(min(es), max(es)), las=1, bty="n", 
                                          ylab="S parameter", xlab="Angle in radians"), p3p[names(p3p) %in% c("ylim", "las", "bty", "xlab", "ylab")]), 
                          list(x=ag, y=rep(0, length(ag)), type="n"))
        do.call("plot", args = ppx)
        for (i in 1:replicates.CI)  {
          ppx <- modifyList(modifyList(list(lwd=0.1, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.2)), p3p[names(p3p) %in% c("lwd", "col", "xlab", "ylab")]), 
                            list(x=ag, y=es[, i]))
          do.call("lines", args = ppx)
        }
        qrd <- sapply(X = 1:nrow(es), FUN = function(i) quantile(es[i, ], probs=c(0.025, 0.5, 0.975)))
        lines(ag, qrd["50%", ], lty=2, lwd=0.5, col="black")
        lines(ag, qrd["2.5%", ], lty=3, lwd=0.5, col="black")
        lines(ag, qrd["97.5%", ], lty=3, lwd=0.5, col="black")
      }
      
      if (any(parameter.name =="min")) {
        es <- sapply(1:replicates.CI, FUN = function(i) (rd[i, "Min"]+ rd[i, "MinSin"]*sinag + rd[i, "MinCos"]*cosag))
        ppx <- modifyList(modifyList(list(xlim=c(min(ag), max(ag)), ylim=c(min(es), max(es)), las=1, bty="n", 
                                          ylab="Min parameter", xlab="Angle in radians"), p3p[names(p3p) %in% c("ylim", "las", "bty", "xlab", "ylab")]), 
                          list(x=ag, y=rep(0, length(ag)), type="n"))
        do.call("plot", args = ppx)
        for (i in 1:replicates.CI)  {
          ppx <- modifyList(modifyList(list(lwd=0.1, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.2)), 
                                       p3p[names(p3p) %in% c("lwd", "col", "xlab", "ylab")]), 
                            list(x=ag, y=es[, i]))
          do.call("lines", args = ppx)
        }
        qrd <- sapply(X = 1:nrow(es), FUN = function(i) quantile(es[i, ], probs=c(0.025, 0.5, 0.975)))
        lines(ag, qrd["50%", ], lty=2, lwd=0.5, col="black")
        lines(ag, qrd["2.5%", ], lty=3, lwd=0.5, col="black")
        lines(ag, qrd["97.5%", ], lty=3, lwd=0.5, col="black")
      }
      
      if (any(parameter.name =="max")) {
        es <- sapply(1:replicates.CI, FUN = function(i) (rd[i, "Max"]+ rd[i, "MaxSin"]*sinag + rd[i, "MaxCos"]*cosag))
        ppx <- modifyList(modifyList(list(xlim=c(min(ag), max(ag)), ylim=c(min(es), max(es)), las=1, bty="n", 
                                          ylab="Max parameter", xlab="Angle in radians"), p3p[names(p3p) %in% c("ylim", "las", "bty", "xlab", "ylab")]), 
                          list(x=ag, y=rep(0, length(ag)), type="n"))
        do.call("plot", args = ppx)
        for (i in 1:replicates.CI)  {
          ppx <- modifyList(modifyList(list(lwd=0.1, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.2)), 
                                       p3p[names(p3p) %in% c("lwd", "col", "xlab", "ylab")]), 
                            list(x=ag, y=es[, i]))
          do.call("lines", args = ppx)
        }
        qrd <- sapply(X = 1:nrow(es), FUN = function(i) quantile(es[i, ], probs=c(0.025, 0.5, 0.975)))
        lines(ag, qrd["50%", ], lty=2, lwd=0.5, col="black")
        lines(ag, qrd["2.5%", ], lty=3, lwd=0.5, col="black")
        lines(ag, qrd["97.5%", ], lty=3, lwd=0.5, col="black")
      }
      
      if (any(parameter.name =="k1")) {
        es <- sapply(1:replicates.CI, FUN = function(i) (rd[i, "K1"]+ rd[i, "K1Sin"]*sinag + rd[i, "K1Cos"]*cosag))
        ppx <- modifyList(modifyList(list(xlim=c(min(ag), max(ag)), ylim=c(min(es), max(es)), las=1, bty="n", 
                                          ylab="K1 parameter", xlab="Angle in radians"), p3p[names(p3p) %in% c("ylim", "las", "bty", "xlab", "ylab")]), 
                          list(x=ag, y=rep(0, length(ag)), type="n"))
        do.call("plot", args = ppx)
        for (i in 1:replicates.CI)  {
          ppx <- modifyList(modifyList(list(lwd=0.1, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.2)), 
                                       p3p[names(p3p) %in% c("lwd", "col", "xlab", "ylab")]), 
                            list(x=ag, y=es[, i]))
          do.call("lines", args = ppx)
        }
        qrd <- sapply(X = 1:nrow(es), FUN = function(i) quantile(es[i, ], probs=c(0.025, 0.5, 0.975)))
        lines(ag, qrd["50%", ], lty=2, lwd=0.5, col="black")
        lines(ag, qrd["2.5%", ], lty=3, lwd=0.5, col="black")
        lines(ag, qrd["97.5%", ], lty=3, lwd=0.5, col="black")
      }
      
      if (any(parameter.name =="k2")) {
        es <- sapply(1:replicates.CI, FUN = function(i) (rd[i, "K2"]+ rd[i, "K2Sin"]*sinag + rd[i, "K2Cos"]*cosag))
        ppx <- modifyList(modifyList(list(xlim=c(min(ag), max(ag)), ylim=c(min(es), max(es)), las=1, bty="n", 
                                          ylab="K2 parameter", xlab="Angle in radians"), p3p[names(p3p) %in% c("ylim", "las", "bty", "xlab", "ylab")]), 
                          list(x=ag, y=rep(0, length(ag)), type="n"))
        do.call("plot", args = ppx)
        for (i in 1:replicates.CI)  {
          ppx <- modifyList(modifyList(list(lwd=0.1, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.2)), 
                                       p3p[names(p3p) %in% c("lwd", "col", "xlab", "ylab")]), 
                            list(x=ag, y=es[, i]))
          do.call("lines", args = ppx)
        }
        qrd <- sapply(X = 1:nrow(es), FUN = function(i) quantile(es[i, ], probs=c(0.025, 0.5, 0.975)))
        lines(ag, qrd["50%", ], lty=2, lwd=0.5, col="black")
        lines(ag, qrd["2.5%", ], lty=3, lwd=0.5, col="black")
        lines(ag, qrd["97.5%", ], lty=3, lwd=0.5, col="black")
      }
      
      
      if (any(parameter.name == "compactness")) {
        if (is.null(o$PeriodicCompactness)) {
          message("To plot the periodic compactness, you must run the analysis BP_FitMLPeriodicCompactness() with replicate.CI different from NULL.")
        } else {
          bc <- o$PeriodicCompactness[, , 2, drop=TRUE]
          
          p <- RM_get(x=x, RMname = analysis, valuename = "peripherie")
          pa <- p$angle.begin
          
          ag <- RM_get(x=x, RMname = analysis, valuename = "cut.angle") 
          if ((RM_get(x=x, RMname = analysis, valuename = "partial")) & (!show.all.angles)) {
            pos <- which(!(round((ag %% (2*pi)), 5) %in% round((pa %% (2*pi)), 5)))
            pos <- pos[-length(pos)]
            bc[pos, ] <- NA
          }
          ag <- (ag[-1] + ag[-length(ag)])/2
          
          ac <- RM_get(x=x, RMname = analysis, valuename = "cut.distance.center")
          ac <- (ac[-length(ac)]+ac[-1])/2
          
          ppx <- modifyList(modifyList(list(xlab="Angle in radians", ylab="Distance from center", legend.lab="Compacity"), 
                                       p3p[names(p3p) %in% c("col", "xlab", "ylab")]), list(x=bc, xaxt="n", yaxt="n"))
          
          do.call(fields::image.plot, args = ppx)
          axis(1, at=seq(from=0, to=1, length.out=dim(bc)[1]), labels = specify_decimal(ag, 2))
          axis(2, at=seq(from=0, to=1, length.out=dim(bc)[2]), labels = specify_decimal(ac, 2), las=1)
        }
      }
      
      if (any(parameter.name =="averagemodel")) {
        if (is.null(o$PeriodicCompactness)) {
          message("To plot the average model, you must run the analysis BP_FitMLPeriodicCompactness() with replicate.CI different from NULL.")
        } else {
          gc <- o$GlobalCompactness
          ac <- RM_get(x=x, RMname = analysis, valuename = "cut.distance.center")
          ac <- (ac[-length(ac)]+ac[-1])/2
          
          ppx <- modifyList(modifyList(list(xlab="Distance from center", ylab="Compacity", bty="n", las=1, ylim=c(0, 1)), 
                                       p3p[names(p3p) %in% c("col", "xlab", "ylab")]), 
                            list(x=ac, y=gc[, "50%"], type="l", lty=2))
          
          do.call("plot", args = ppx)
          lines(x=ac, y=gc[, "2.5%"], lty=3)
          lines(x=ac, y=gc[, "97.5%"], lty=3)
          if (show.legend) legend("topleft", legend = c("Median", "95% CI including periodic effect"), lty=c(2, 3))
        }
      }
    } 
  }
  
  if (type == "3dcolors") {
    
    threshold <- RM_get(x=x, RMname=analysis, valuename = "threshold")
    
    DF_background <- data.frame(Red=as.numeric(x[, , 1, 1])[as.vector(!threshold[])], 
                                Green=as.numeric(x[, , 1, 2])[as.vector(!threshold[])], 
                                Blue=as.numeric(x[, , 1, 3])[as.vector(!threshold[])])
    
    DF_foreground <- data.frame(Red=as.numeric(x[, , 1, 1])[as.vector(threshold[])], 
                                Green=as.numeric(x[, , 1, 2])[as.vector(threshold[])], 
                                Blue=as.numeric(x[, , 1, 3])[as.vector(threshold[])])
    DF <- rbind(DF_background, DF_foreground)
    
    # install.packages("scatterplot3d") # Install
    # library("scatterplot3d") # load
    getFromNamespace(x="scatterplot3d", ns="scatterplot3d")(DF[,1:3], xlim = c(0, 1), ylim = c(0, 1),
                                                            zlim = c(0, 1), xlab = "Red", ylab = "Green", zlab = "Blue", 
                                                            pch=c(rep(19, nrow(DF_background)), rep(21, nrow(DF_foreground))), 
                                                            color=rgb(red = DF[,1], green = DF[,2], blue = DF[,3], alpha = 1), 
                                                            angle = angle.3D, 
                                                            bg="black"
    )
    if (show.legend) {
      legend("topleft", legend=c("Background", "Foreground"), 
             col=c(rgb(red = mean(DF_background[,1]), green = mean(DF_background[,2]), blue = mean(DF_background[,3]), alpha = 1), 
                   "black"), 
             pch=c(19, 21), pt.bg=c(rgb(red = mean(DF_background[,1]), green = mean(DF_background[,2]), blue = mean(DF_background[,3]), alpha = 1), 
                                    rgb(red = mean(DF_foreground[,1]), green = mean(DF_foreground[,2]), blue = mean(DF_foreground[,3]), alpha = 1)))
    }
  }
  
  if (type == "colors") {
    
    DF <- data.frame(Red=as.numeric(x[, , 1, 1]), 
                     Green=as.numeric(x[, , 1, 2]), 
                     Blue=as.numeric(x[, , 1, 3]))
    if (!is.null(analysis)) {
      bg <- col2rgb(RM_get(x=x, RMname=analysis, valuename = "bg"))/255
      fg <- col2rgb(RM_get(x=x, RMname=analysis, valuename = "fg"))/255
    } else {
      bg <- NULL
      fg <- NULL
    }
    
    layout(1:3)
    ppar <- par(mar=c(2, 4, 2, 1))
    par(xpd=TRUE)
    hist(DF$Red, main="", xlab="", col="red")
    ymax <- ScalePreviousPlot()$ylim["end"]
    if (!is.null(bg)) {
      segments(x0=bg["red", 1], x1=bg["red", 1], y0=0, y1=ymax)
      text(x=bg["red", 1], y=ymax*1.1, labels = "Foreground")
    }
    if (!is.null(fg)) {
      segments(x0=fg["red", 1], x1=fg["red", 1], y0=0, y1=ymax)
      text(x=fg["red", 1], y=ymax*1.1, labels = "Background")
    }
    hist(DF$Green, main="", xlab="", col="green")
    ymax <- ScalePreviousPlot()$ylim["end"]
    if (!is.null(bg)) {
      segments(x0=bg["green", 1], x1=bg["green", 1], y0=0, y1=ymax)
      text(x=bg["green", 1], y=ymax*1.1, labels = "Foreground")
    }
    if (!is.null(fg)) {
      segments(x0=fg["green", 1], x1=fg["green", 1], y0=0, y1=ymax)
      text(x=fg["green", 1], y=ymax*1.1, labels = "Background")
    }
    
    hist(DF$Blue, main="", xlab="", col="blue")
    ymax <- ScalePreviousPlot()$ylim["end"]
    if (!is.null(bg)) {
      segments(x0=bg["blue", 1], x1=bg["blue", 1], y0=0, y1=ymax)
      text(x=bg["blue", 1], y=ymax*1.1, labels = "Foreground")
    }
    if (!is.null(fg)) {
      segments(x0=fg["blue", 1], x1=fg["blue", 1], y0=0, y1=ymax)
      text(x=fg["blue", 1], y=ymax*1.1, labels = "Background")
    }
    layout(1)
    par(mar=ppar$mar)
  }
  
  if (type == "mcmc") {
    
    parameter.name <- match.arg(parameter.name, choices = c("P", "S", "Min", "Max", "K1", "K2"))
    outMCMC <- RM_get(x = x, RM = "RM", RMname = analysis, valuename = "mcmc")
    par(xpd=FALSE)
    if (!is.null(outMCMC)) {
      options.mcmc <- modifyList(options.mcmc, list(x=outMCMC, parameters=parameter.name))
      out <- do.call(what = getFromNamespace("plot.mcmcComposite", ns="HelpersMG"), args=options.mcmc)
    } else {
      out <- NULL
    }
  } 
  
  if (type == "radial") {
    if (is.null(RM_get(x=x, RMname=analysis, valuename = "optimRadial"))) stop("Radial analysis has not been perfomed")
    
    if (is.numeric(analysis)) analysis <- names(RM_list(x=x, silent = TRUE))[analysis]
    
    par(mar=c(4, 4, 2, 1)+0.4)
    out <- RM_get(x=x, RMname=analysis, valuename = "optimRadial")$synthesis
    
    # parameter.name <- match.arg(parameter.name, choices = colnames(out))
    
    if (length(parameter.name) != 1) {
      layout(mat = 1:length(parameter.name))
    } else {
      layout(1)
    }
    angles <- RM_get(x=x, RMname=analysis, valuename = "optimRadial")$angles
    
    t <- "l"
    if (length(angles) == 1) t = "p"
    
    anglestot <- RM_get(x=x, RMname=analysis, valuename = "cut.angle")
    
    for (i in parameter.name) {
      ylim <- c(0, 1)
      if (i %in% c("S", "K1", "K2")) ylim=NULL
      plot(angles, out[, i], las=1, bty="n", xlab="Angle", 
           ylab=i, type=t, ylim=ylim, xaxt="n", xlim=c(-pi, pi), pch=19)
      axis(1, at=seq(from=-pi, to=pi, by=anglestot[2]-anglestot[1]), 
           labels = specify_decimal(seq(from=-pi, to=pi, by=anglestot[2]-anglestot[1]), decimals = 2), cex.axis=0.5, las=2)
      axis(1, at=seq(from=-pi, to=pi, by=2*(anglestot[2]-anglestot[1])), 
           labels = FALSE, lwd.ticks = 2)
      
    }
    
  } 
  
  if (type %in% c("observations", "model", "observations+model")) {
    
    if (!is.null(angle) & (!is.null(RM_get(x=x, RMname=analysis, valuename = "optimRadial")))) {
      # Je montre une tranche
      distance.center <- RM_get(x=x, RMname=analysis, valuename = "compactness.synthesis")$distance.center
      angles <- RM_get(x=x, RMname=analysis, valuename = "optimRadial")$angles
      indice.angle <- which.min(abs(angles-angle))[1]
      
      main <- paste0(" : Angle [", specify_decimal(angles[indice.angle], decimals = 3), 
                     ",", specify_decimal(angles[ifelse(indice.angle == length(angles), 1, indice.angle+1)], decimals = 3), "]")
      
      array.compactness <- RM_get(x=x, RMname=analysis, valuename = "array.compactness")
      # angles <- RM_get(x=x, RMname=analysis, valuename = "optimRadial")$angles
      # indice.angle <- which.min(abs(angles-angle))
      
      data_nm <- array.compactness[indice.angle, , "0"]
      data_m <- array.compactness[indice.angle, , "1"]
      
      
      compactness.synthesis <- data.frame(distance.center=distance.center, 
                                          mineralized=data_m, 
                                          unmineralized=data_nm, 
                                          compactness=data_m/(data_m+data_nm))
      
      if ((type=="model") | (type=="observations+model")) {
        p <- RM_get(x=x, RMname=analysis, valuename = "optimRadial")$synthesis[indice.angle, , drop=TRUE]
        p <- unlist(p[unlist(lapply(p, FUN = is.numeric))])
      }
    } else {
      
      # Je montre tout
      main <- ""
      # if ((type =="observations") | (type=="observations+model")) {
      compactness.synthesis <- RM_get(x=x, RMname=analysis, valuename = "compactness.synthesis")
      # }
      if ((type=="model") | (type=="observations+model")) {
        p <- c(RM_get(x=x, RMname=analysis, valuename = "optim")$par, RM_get(x=x, RMname=analysis, valuename = "optim")$fixed.parameters)
      }
    }
    
    if (is.null(compactness.synthesis)) {
      stop("Bone section has not still been analyzed. Use BP_EstimateCompactness().")
    }
    
    if (type == "observations") {
      par(xaxs="i", yaxs="r")
      if (is.null(mar)) {
        par(mar=c(4, 4, 2, 5)+0.4)
      } else {
        par(mar=mar)
      }
      m <- compactness.synthesis$mineralized
      nm <- compactness.synthesis$unmineralize
      plot(compactness.synthesis$distance.center, compactness.synthesis$compactness, xlim=c(0, 1), ylim=c(0, 1), 
           type="l", las=1, bty="n", xlab="Distance from the center", ylab="Compactness", lwd=2, main = main)
      lines(x=compactness.synthesis$distance.center[(m+nm) != 0], y=((m+nm)/max(m+nm))[(m+nm) != 0], col="blue")
      axis(side = 4, at=seq(from=0, to=1, by=0.2), labels = round(seq(from=0, to=1, by=0.2)*max(m+nm), 0), 
           las=1, col.axis="blue", col="blue")
      mtext("Number of pixels", side=4, line=3, col="blue")
      out <- data.frame(distance.center=compactness.synthesis$distance.center, 
                        observed.compactness=compactness.synthesis$compactness)
      
      if (show.legend) {
        legend("bottomright", legend=c("Number of pixels", "Observed compactness"), 
               lty=c(1, 1), lwd=c(1, 2), col=c("blue", "black"), cex=0.8)
      }
      
    }
    
    if (type == "model") {
      if (is.numeric(analysis)) analysis <- names(RM_list(x=x, silent = TRUE))[analysis]
      
      # Min <- p["Min"]
      # Max <- p["Max"]
      
      # 21/02/2020
      # p["S"] <- 1/(4*p["S"])
      
      
      c <- BP_flexit(x = compactness.synthesis$distance.center, 
                     par = p) # * (Max - Min) + Min
      
      out <- data.frame(distance.center=compactness.synthesis$distance.center, 
                        modeled.compactness=c)
      par(xaxs="i", yaxs="r")
      if (is.null(mar)) {
        par(mar=c(4, 4, 2, 5)+0.4)
      } else {
        par(mar=mar)
      }
      
      plot(x=compactness.synthesis$distance.center, y=c, xlim=c(0, 1), ylim=c(0, 1), 
           type="n", las=1, bty="n", xlab="Distance from the center", ylab="Compactness", lwd=2, lty=3, 
           main=paste(analysis, main))
      
      if ((CI == "MCMC") & (is.null(angle)) & (!is.null(RM_get(x=x, RMname=analysis, valuename = "mcmc")))) {
        polygon(x=c(compactness.synthesis$distance.center, rev(compactness.synthesis$distance.center)), 
                y=c(RM_get(x=x, RMname=analysis, valuename = "mcmc")$quantiles["2.5%", ], rev(RM_get(x=x, RMname=analysis, valuename = "mcmc")$quantiles["97.5%", ])), 
                col="lightgrey", border="lightgrey", lwd=3)
      }
      if ((CI == "ML") & (is.null(angle)) & (!is.null(RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles))) {
        polygon(x=c(compactness.synthesis$distance.center, rev(compactness.synthesis$distance.center)), 
                y=c(RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles["2.5%", ], rev(RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles["97.5%", ])), 
                col="lightgrey", border="lightgrey", lwd=3)
      }
      
      lines(x=compactness.synthesis$distance.center, y=c, lwd=2, lty=3)
      
      if (show.legend) {
        if ((CI == "MCMC") & (is.null(angle)) & (!is.null(RM_get(x=x, RMname=analysis, valuename = "optim")))) {
          legend("bottomright", legend=c("Model", "95% Credibility interval MCMC"), 
                 lty=c(3, 1), lwd=c(2, 6), col=c("black", "lightgrey"), cex=0.8)
        } else {
          if ((CI == "ML") & (is.null(angle)) & (!is.null(RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles))) {
            legend("bottomright", legend=c("Model", "95% Confidence interval ML"), 
                   lty=c(3, 1), lwd=c(2, 6), col=c("black", "lightgrey"), cex=0.8)
          } else {
            legend("bottomright", legend=c("Model"), 
                   lty=c(3), lwd=c(2), col=c("black"), cex=0.8)
          }
        }
      }
      
    }
    
    if (type == "observations+model") {
      if (is.numeric(analysis)) analysis <- names(RM_list(x=x, silent = TRUE))[analysis]
      
      par(xaxs="i", yaxs="r")
      if (is.null(mar)) {
        par(mar=c(4, 4, 2, 5)+0.4)
      } else {
        par(mar=mar)
      }
      m <- compactness.synthesis$mineralized
      nm <- compactness.synthesis$unmineralize
      plot(compactness.synthesis$distance.center, compactness.synthesis$compactness, xlim=c(0, 1), ylim=c(0, 1), 
           type="n", las=1, bty="n", xlab="Distance from the center", ylab="Compactness", lwd=2, 
           main=paste(analysis, main))
      out <- data.frame(distance.center=compactness.synthesis$distance.center, 
                        observed.compactness=compactness.synthesis$compactness)
      
      if ((CI == "MCMC") & (is.null(angle)) & (!is.null(RM_get(x=x, RMname=analysis, valuename = "mcmc")))) {
        polygon(x=c(compactness.synthesis$distance.center, rev(compactness.synthesis$distance.center)), 
                y=c(RM_get(x=x, RMname=analysis, valuename = "mcmc")$quantiles["2.5%", ], rev(RM_get(x=x, RMname=analysis, valuename = "mcmc")$quantiles["97.5%", ])), 
                col="lightgrey", border="lightgrey", lwd=3)
        lines(x=compactness.synthesis$distance.center, y=RM_get(x=x, RMname=analysis, valuename = "mcmc")$quantiles["50%", ], xlim=c(0, 1), lwd=2, lty=3)
        out <- cbind(out, median.modeled.compactness=RM_get(x=x, RMname=analysis, valuename = "mcmc")$quantiles["50%", ])
      }
      if ((CI == "ML") & (is.null(angle)) & (!is.null(RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles))) {
        polygon(x=c(compactness.synthesis$distance.center, rev(compactness.synthesis$distance.center)), 
                y=c(RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles["2.5%", ], rev(RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles["97.5%", ])), 
                col="lightgrey", border="lightgrey", lwd=3)
        lines(x=compactness.synthesis$distance.center, y=RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles["50%", ], xlim=c(0, 1), lwd=2, lty=3)
        out <- cbind(out, median.modeled.compactness=RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles["50%", ])
      }
      
      lines(compactness.synthesis$distance.center, compactness.synthesis$compactness, lwd=2)
      
      lines(x=compactness.synthesis$distance.center, y=(m+nm)/max(m+nm), col="blue")
      axis(side = 4, at=seq(from=0, to=1, by=0.2), labels = round(seq(from=0, to=1, by=0.2)*max(m+nm), 0), 
           las=1, col.axis="blue", col="blue")
      mtext("Number of pixels", side=4, line=3, col="blue")
      
      # p <- c(RM_get(x=x, RMname=analysis, valuename = "optim")$par, RM_get(x=x, RMname=analysis, valuename = "optim")$fixed.parameters)
      
      # Min <- p["Min"]
      # Max <- p["Max"]
      
      # 21/02/2020
      # p["S"] <- 1/(4*p["S"])
      
      c <- BP_flexit(x = compactness.synthesis$distance.center, 
                     par = p) # * (Max - Min) + Min
      
      out <- cbind(out, modeled.compactness=c)
      
      lines(x=compactness.synthesis$distance.center, y=c, xlim=c(0, 1), lwd=2, lty=4)
      
      
      if (show.legend) {
        if ((CI == "MCMC") & (is.null(angle)) & (!is.null(RM_get(x=x, RMname=analysis, valuename = "mcmc")))) {
          legend("bottomright", legend=c("Number of pixels", "Observed compactness", "MCMC median model", "95% Credibility interval MCMC", "ML model"), 
                 lty=c(1, 1, 3, 1, 4), lwd=c(1, 2, 2, 6, 2), col=c("blue", "black", "black", "lightgrey", "black"), cex=0.8)
        } else {
          if ((CI == "ML") & (is.null(angle)) & (!is.null(RM_get(x=x, RMname=analysis, valuename = "optim")$quantiles))) {
            legend("bottomright", legend=c("Number of pixels", "Observed compactness", "ML median model", "95% Confidence interval ML", "ML model"), 
                   lty=c(1, 1, 3, 1, 4), lwd=c(1, 2, 2, 6, 2), col=c("blue", "black", "black", "lightgrey", "black"), cex=0.8)
          } else {
            legend("bottomright", legend=c("Number of pixels", "Observed compactness", "Model"), 
                   lty=c(1, 1, 4), lwd=c(1, 2, 2), col=c("blue", "black", "black"), cex=0.8)
          }
        }
      }
    }
  }
  
  if (type %in% c("original", "mineralized", "unmineralized", "section")) {
    # layout(1)
    out <- NULL
    bone_x <- NULL
    
    threshold <- RM_get(x=x, RMname=analysis, valuename = "threshold")
    contour <- RM_get(x=x, RMname=analysis, valuename = "contour")
    
    
    # cbin <- ifelse(contour, 1, 0)
    
    bone <- NULL
    
    if (type != "original") {
      
      
      if ((type == "mineralized") & !is.null(threshold)) bone_x <- threshold
      if ((type == "unmineralized") & !is.null(threshold) & !is.null(contour)) bone_x <- contour & !threshold
      if ((type == "section") & !is.null(contour)) bone_x <- contour
      
      if (!is.null(bone_x)) {
        bone <- x
        bone[, , 1, 1] <- !bone_x
        if (dim(bone)[4]>1) bone[, , 1, 2] <- !bone_x
        if (dim(bone)[4]>2) bone[, , 1, 3] <- !bone_x
      }
    } else {
      bone <- x
    }
    
    if (!is.null(bone)) {
      
      
      par(xaxs="i", yaxs="i")
      
      if (is.null(mar)) {
        par(mar=c(4, 0, 0, 0))
      } else {
        par(mar=mar)
      }
      
      getFromNamespace("plot.cimg", ns="imager")(bone, bty="n", axes=FALSE, xlab="", ylab="", 
                                                 asp = 1)
      xl <- ScalePreviousPlot()$xlim
      yl <- ScalePreviousPlot()$ylim
      par(xpd=TRUE)
      
      if (!is.null(message)) {
        text(x=xl["begin"], 
             y=yl["begin"]-yl["range"]*0.05, 
             labels = message, pos=4)
      }
      
      if (show.colors) {
        bg <- RM_get(x=x, RMname=analysis, valuename = "bg")
        if (!is.null(bg)) {
          text(x=xl["begin"]+xl["range"]*0.08, 
               y=yl["begin"]-yl["range"]*0.09, 
               labels = "Background\ncolor", pos=4, cex=0.8)
          polygon(x=c(xl["begin"]+xl["range"]*0.02, 
                      xl["begin"]+xl["range"]*0.07, 
                      xl["begin"]+xl["range"]*0.07, 
                      xl["begin"]+xl["range"]*0.02), 
                  y=c(yl["begin"]-yl["range"]*0.12, 
                      yl["begin"]-yl["range"]*0.12, 
                      yl["begin"]-yl["range"]*0.07, 
                      yl["begin"]-yl["range"]*0.07), 
                  col=bg)
        }
        fg <- RM_get(x=x, RMname=analysis, valuename = "fg")
        if (!is.null(fg)) {
          text(x=xl["center"]/2+xl["range"]*0.08, 
               y=yl["begin"]-yl["range"]*0.09, 
               labels = "Foreground\ncolor", pos=4, cex=0.8)
          polygon(x=c(xl["center"]/2+xl["range"]*0.02, 
                      xl["center"]/2+xl["range"]*0.07, 
                      xl["center"]/2+xl["range"]*0.07, 
                      xl["center"]/2+xl["range"]*0.02), 
                  y=c(yl["begin"]-yl["range"]*0.12, 
                      yl["begin"]-yl["range"]*0.12, 
                      yl["begin"]-yl["range"]*0.07, 
                      yl["begin"]-yl["range"]*0.07), 
                  col=fg)
        }
      }
      
      
      
      if (show.grid & !is.null(RM_get(x=x, RMname=analysis, valuename = "peripherie"))) {
        # J'affiche la grille
        # compactness <- RM_get(x=x, RMname=analysis, valuename = "compactness")
        peripherie <- RM_get(x=x, RMname=analysis, valuename = "peripherie")
        # angles <- RM_get(x=x, RMname=analysis, valuename = "cut.angle")
        
        if (RM_get(x=x, RMname=analysis, valuename = "partial")) {
          # angles <- (angles %% (2*pi))
          
          angles_tot <- unique(c(peripherie[peripherie[, "include.begin"], "angle.begin"] %% (2*pi), 
                                 peripherie[peripherie[, "include.last"], "angle.last"] %% (2*pi)))
          
          posdeb <- which(peripherie[, "include.begin"])[1]
          poslast <- rev(which(peripherie[, "include.last"]))[1]
          
          angle_ec <- c(peripherie[posdeb:poslast, "angle.begin"], peripherie[poslast, "angle.last"])
          
          # angle_ec <- peripherie[(peripherie$angle.begin >= angles[1]) & (peripherie$angle.last <= angles[length(angles)]), "angle.center"]
          angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
          # angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
          
          x_peripherie <- cos(angle_ecp)*c(peripherie[posdeb:poslast, "peripherie.begin"], peripherie[poslast, "peripherie.last"])
          y_peripherie <- sin(angle_ecp)*c(peripherie[posdeb:poslast, "peripherie.begin"], peripherie[poslast, "peripherie.last"])
          
        } else {
          angles_tot <- peripherie[peripherie[, "include.begin"], "angle.begin"] %% (2*pi)
          angle_ec <- c(peripherie$angle.begin, peripherie[nrow(peripherie), "angle.last"])
          angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
          # angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
          x_peripherie <- cos(angle_ecp)*c(peripherie$peripherie.begin, peripherie[nrow(peripherie), "peripherie.last"])
          y_peripherie <- sin(angle_ecp)*c(peripherie$peripherie.begin, peripherie[nrow(peripherie), "peripherie.last"])
        }
        
        for (angle_ec in angles_tot) {
          
          angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
          # angle_ecp <- ifelse(angle_ecp >= pi, -(2*pi)+angle_ecp, angle_ecp)
          
          if (min(abs(peripherie$angle.begin-angle_ec)) < min(abs(peripherie$angle.last-angle_ec))) {
            md <- peripherie$peripherie.begin[which.min(abs(peripherie$angle.begin-angle_ec))][1]
          } else {
            md <- peripherie$peripherie.last[which.min(abs(peripherie$angle.last-angle_ec))][1]
          }
          
          segments(x0=RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"], 
                   x1=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md), 
                   y0=RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"], 
                   y1=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md), 
                   col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.8))
        }
        
        bg <- col2rgb(RM_get(x=x, RMname=analysis, valuename = "bg"))/255
        fg <- col2rgb(RM_get(x=x, RMname=analysis, valuename = "fg"))/255
        
        par(xpd=TRUE)
        # L'angle 0 deg
        angle_ec <- 0
        angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
        angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
        
        md <- peripherie$peripherie.peripherie[which.min(abs(peripherie$angle.center-angle_ec))]*0.9
        text(x=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md*1.2), 
             y=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md*1.2), 
             labels = "0", col=(fg + bg)/2)
        angle_ec <- pi
        angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
        angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
        
        md <- peripherie$peripherie.peripherie[which.min(abs(peripherie$angle.center-angle_ec))]*0.9
        text(x=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md*1.2), 
             y=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md*1.2), 
             labels = "pi", col=(fg + bg)/2)
        angle_ec <- pi/2
        angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
        angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
        
        md <- peripherie$peripherie.peripherie[which.min(abs(peripherie$angle.center-angle_ec))]*0.9
        text(x=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md*1.2), 
             y=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md*1.2), 
             labels = "pi/2", col=(fg + bg)/2)
        angle_ec <- -pi/2
        angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
        angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
        
        md <- peripherie$peripherie.peripherie[which.min(abs(peripherie$angle.center-angle_ec))]*0.9
        text(x=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md*1.2), 
             y=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md*1.2), 
             labels = "-pi/2", col=(fg + bg)/2)
        
        
        
        
        for (ratio in RM_get(x=x, RMname=analysis, valuename = "cut.distance.center")[-1]) {
          lines(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+x_peripherie*ratio, 
                RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+y_peripherie*ratio, 
                col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.8))
        }
        
      }
      
      if (show.centers) {
        centers <- RM_get(x=x, RMname=analysis, valuename = "centers")
        if (!is.null(centers)) {
          if (!is.na(centers["GC_cortex.x"])) {
            points(centers["GC_cortex.x"], centers["GC_cortex.y"], pch=4, col="red")
            points(centers["GC_bone.x"], centers["GC_bone.y"], pch=3, col="red")
            points(centers["GC_medula.x"], centers["GC_medula.y"], pch=1, col="red")
            points(centers["GC_ontogenic.x"], centers["GC_ontogenic.y"], pch=19, col="blue")
            text(x=xl["end"]-xl["range"]*0.40, 
                 y=yl["begin"]-yl["range"]*0.01, 
                 labels = "X Center of the mineralized part", cex=0.8, col="red", pos=4)
            text(x=xl["end"]-xl["range"]*0.40, 
                 y=yl["begin"]-yl["range"]*0.04, 
                 labels = "O Center of the non-mineralized part", cex=0.8, col="red", pos=4)
            text(x=xl["end"]-xl["range"]*0.40, 
                 y=yl["begin"]-yl["range"]*0.07, 
                 labels = "+ Center of the section", cex=0.8, col="red", pos=4)
            points(x=xl["end"]-xl["range"]*0.375, 
                   y=yl["begin"]-yl["range"]*0.1, 
                   pch=19, col="blue")
            text(x=xl["end"]-xl["range"]*0.38, 
                 y=yl["begin"]-yl["range"]*0.1, 
                 labels = "Ontogenetic center", cex=0.8, col="blue", pos=4)
          } else {
            points(centers["GC_user.x"], centers["GC_user.y"], pch=1, col="red")
            text(x=xl["end"]-xl["range"]*0.40, 
                 y=yl["begin"]-yl["range"]*0.02, 
                 labels = "O User-defined center", cex=0.8, col="red", pos=4)
          }
        }
      }
    } else {
      stop("The information is not available")
    }
  }
  
  
  return(invisible(out))
}


