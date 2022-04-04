#' plot.BoneProfileR displays a bone section
#' @title Plot a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Nothing
#' @param x The bone image
#' @param message The message to be displayed
#' @param type The type of plot; see description
#' @param angle Which angle model to show
#' @param parameter.mcmc The posterior parameter to show for type = "mcmc"
#' @param options.mcmc The option to plot type mcmc output
#' @param show.colors Should the background and foreground colors be shown?
#' @param show.centers Should the centers be shown?
#' @param show.grid Should the grid be shown?
#' @param show.legend Should a legend be shown? 
#' @param analysis Name or number of analysis to be plotted
#' @param radial.variable Name of the radial variable to plot
#' @param CI Which confidence interval should be plotted: MCMC or ML
#' @param restorePar If TRUE, restore the par parameter at the exit
#' @param mar The margin for type being "model" or "observations"
#' @param angle.3D The angle between x and y for 3Dcolors graph
#' @param ... Not used
#' @description Display a bone section.\cr
#' type value can be:\cr
#' Image plot: original, mineralized, unmineralized, section\cr
#' Original is the original image, mineralized is the mineral interpretation of the section, 
#' unmineralized is the unmineralized interpretation of the section, section is the interpretation of the section.\cr
#' colors show the histograms of pixel information with foreground and background colors if they are defined.\cr
#' 3Dcolors show the pixels colors in 3D\cr
#' Global analysis: observations, model, observations+model\cr
#' Radial analysis: radial\cr
#' If angle is not null and a radial analysis exists, it will show the model for this angle.\cr
#' mcmc: It will show the posterior distribution of parameter
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
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic", rotation.angle = 1)
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic")
#'  plot(bone)
#'  # 
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
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
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  bone <- BP_FitMLCompactness(bone, 
#'                 fitted.parameters=c(fittedpar, K1=1, K2=1), 
#'                 fixed.parameters=NULL, analysis="flexit")
#'  compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", alloptim=TRUE), 
#'              Flexit=BP_GetFittedParameters(bone, analysis="flexit", alloptim=TRUE))
#'  out4p <- plot(bone, type="observations+model", analysis="logistic")
#'  out6p <- plot(bone, type="observations+model", analysis="flexit")
#'  bone <- BP_FitBayesianCompactness(bone, analysis="logistic")
#'  plot(bone, type="observations+model", CI="MCMC")
#'  bone <- BP_FitMLRadialCompactness(bone)
#'  plot(bone, type="radial", radial.variable=c("P", "S"))
#'  plot(bone, type="radial", radial.variable=c("P", "S", "Min", "Max"))
#' }
#' @method plot BoneProfileR
#' @export


plot.BoneProfileR <- function(x, message=NULL, type="original", angle=NULL, 
                              show.centers=TRUE, 
                              show.colors=TRUE, show.grid=TRUE, analysis=1, 
                              parameter.mcmc = "S", 
                              options.mcmc = list(), 
                              restorePar = TRUE, 
                              mar=NULL, angle.3D = 55, 
                              CI="ML", radial.variable= "S", show.legend=TRUE, ...) {
  
  # message=NULL; type="original"; angle=NULL; parameter.mcmc = "S"; options.mcmc = list(); restorePar=TRUE; mar=NULL; show.centers=TRUE; show.colors=TRUE; show.grid=TRUE; analysis=1; CI="ML"; radial.variable= "S"; show.legend=TRUE
  
  # type <- "observations"
  # type <- "radial"
  # type <- "model"
  # type <- "colors"
  bone <- x
  
  oldpar <- par(no.readonly = TRUE)    # code line i
  if (restorePar) on.exit(par(oldpar))            # code line i + 1
  
  type <- match.arg(type, choices = c("original", "mineralized", 
                                      "unmineralized", "section", "radial", 
                                      "observations", "model", "observations+model", 
                                      "mcmc", "colors", "3Dcolors"))
  
  out <- BP_ListAnalyses(bone=x)
  if (is.null(out[analysis][[1]]) & (type != "original") & (type != "colors")) {
    stop(paste0("The analysis ", analysis, " does not exist"))
  }
  
  if (type == "3Dcolors") {

    threshold <- RM_get(x=bone, RMname=analysis, valuename = "threshold")
    
    DF_background <- data.frame(Red=as.numeric(bone[, , 1, 1])[as.vector(!threshold[])], 
                                Green=as.numeric(bone[, , 1, 2])[as.vector(!threshold[])], 
                                Blue=as.numeric(bone[, , 1, 3])[as.vector(!threshold[])])
    
    DF_foreground <- data.frame(Red=as.numeric(bone[, , 1, 1])[as.vector(threshold[])], 
                                Green=as.numeric(bone[, , 1, 2])[as.vector(threshold[])], 
                                Blue=as.numeric(bone[, , 1, 3])[as.vector(threshold[])])
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
    
    DF <- data.frame(Red=as.numeric(bone[, , 1, 1]), 
                     Green=as.numeric(bone[, , 1, 2]), 
                     Blue=as.numeric(bone[, , 1, 3]))
    if (!is.null(analysis)) {
      bg <- col2rgb(RM_get(x=bone, RMname=analysis, valuename = "bg"))/255
      fg <- col2rgb(RM_get(x=bone, RMname=analysis, valuename = "fg"))/255
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
    
    parameter.mcmc <- match.arg(parameter.mcmc, choices = c("P", "S", "Min", "Max", "K1", "K2"))
    outMCMC <- RM_get(x = x, RM = "RM", RMname = analysis, valuename = "mcmc")
    par(xpd=FALSE)
    if (!is.null(outMCMC)) {
      options.mcmc <- modifyList(options.mcmc, list(x=outMCMC, parameters=parameter.mcmc))
      out <- do.call(what = getFromNamespace("plot.mcmcComposite", ns="HelpersMG"), args=options.mcmc)
    } else {
      out <- NULL
    }
  } 
  
  if (type == "radial") {
    if (is.null(RM_get(x=bone, RMname=analysis, valuename = "optimRadial"))) stop("Radial analysis has not been perfomed")
    
    if (is.numeric(analysis)) analysis <- names(RM_list(x=bone, silent = TRUE))[analysis]
    
    par(mar=c(4, 4, 2, 1)+0.4)
    out <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")$synthesis
    
    if (length(radial.variable) != 1) layout(mat = 1:length(radial.variable))
    angles <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")$angles
    for (i in radial.variable) {
      ylim <- c(0, 1)
      if (i %in% c("S", "K1", "K2")) ylim=NULL
      plot(angles, out[, i], las=1, bty="n", xlab="Angle", 
           ylab=i, type="l", ylim=ylim, xaxt="n", xlim=c(-pi, pi))
      axis(1, at=seq(from=-pi, to=pi, by=angles[2]-angles[1]), 
           labels = specify_decimal(seq(from=-pi, to=pi, by=angles[2]-angles[1]), decimals = 2), cex.axis=0.5, las=2)
      axis(1, at=seq(from=-pi, to=pi, by=2*(angles[2]-angles[1])), 
           labels = FALSE, lwd.ticks = 2)
      
    }
    
  } 
  if (type %in% c("observations", "model", "observations+model")) {
    
    if (!is.null(angle) & (!is.null(RM_get(x=bone, RMname=analysis, valuename = "optimRadial")))) {
      # Je montre une tranche
      distance.center <- RM_get(x=bone, RMname=analysis, valuename = "compactness.synthesis")$distance.center
      angles <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")$angles
      indice.angle <- which.min(abs(angles-angle))
      
      main <- paste0(" : Angle [", specify_decimal(angles[indice.angle], decimals = 3), 
                     ",", specify_decimal(angles[ifelse(indice.angle == length(angles), 1, indice.angle+1)], decimals = 3), "]")
      
      array.compactness <- RM_get(x=bone, RMname=analysis, valuename = "array.compactness")
      angles <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")$angles
      indice.angle <- which.min(abs(angles-angle))
      
      data_nm <- array.compactness[indice.angle, , "0"]
      data_m <- array.compactness[indice.angle, , "1"]
      
      
      compactness.synthesis <- data.frame(distance.center=distance.center, 
                                          mineralized=data_m, 
                                          unmineralized=data_nm, 
                                          compactness=data_m/(data_m+data_nm))
      
      if ((type=="model") | (type=="observations+model")) {
        p <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")$synthesis[indice.angle, , drop=TRUE]
      }
    } else {
      
      # Je montre tout
      main <- ""
      # if ((type =="observations") | (type=="observations+model")) {
      compactness.synthesis <- RM_get(x=bone, RMname=analysis, valuename = "compactness.synthesis")
      # }
      if ((type=="model") | (type=="observations+model")) {
        p <- c(RM_get(x=bone, RMname=analysis, valuename = "optim")$par, RM_get(x=bone, RMname=analysis, valuename = "optim")$fixed.parameters)
      }
    }
    
    if (is.null(compactness.synthesis)) {
      stop("Bone section has not still been analyzed. Use BP_EstimateCompactness().")
    }
    
    if (type =="observations") {
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
      lines(x=compactness.synthesis$distance.center, y=(m+nm)/max(m+nm), col="blue")
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
      if (is.numeric(analysis)) analysis <- names(RM_list(x=bone, silent = TRUE))[analysis]
      
      Min <- p["Min"]
      Max <- p["Max"]
      
      # 21/02/2020
      p["S"] <- 1/(4*p["S"])
      
      
      c <- flexit(x = compactness.synthesis$distance.center, 
                  par = p) * (Max - Min) + Min
      
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
      
      if ((CI == "MCMC") & (is.null(angle)) & (!is.null(RM_get(x=bone, RMname=analysis, valuename = "mcmc")))) {
        polygon(x=c(compactness.synthesis$distance.center, rev(compactness.synthesis$distance.center)), 
                y=c(RM_get(x=bone, RMname=analysis, valuename = "mcmc")$quantiles["2.5%", ], rev(RM_get(x=bone, RMname=analysis, valuename = "mcmc")$quantiles["97.5%", ])), 
                col="lightgrey", border="lightgrey", lwd=3)
      }
      if ((CI == "ML") & (is.null(angle)) & (!is.null(RM_get(x=bone, RMname=analysis, valuename = "optim")$quantiles))) {
        polygon(x=c(compactness.synthesis$distance.center, rev(compactness.synthesis$distance.center)), 
                y=c(RM_get(x=bone, RMname=analysis, valuename = "optim")$quantiles["2.5%", ], rev(RM_get(x=bone, RMname=analysis, valuename = "optim")$quantiles["97.5%", ])), 
                col="lightgrey", border="lightgrey", lwd=3)
      }
      
      lines(x=compactness.synthesis$distance.center, y=c, lwd=2, lty=3)
      
      if (show.legend) {
        if ((CI == "MCMC") & (is.null(angle)) & (!is.null(RM_get(x=bone, RMname=analysis, valuename = "optim")))) {
          legend("bottomright", legend=c("Model", "95% Credibility interval MCMC"), 
                 lty=c(3, 1), lwd=c(2, 6), col=c("black", "lightgrey"), cex=0.8)
        } else {
          if ((CI == "ML") & (is.null(angle)) & (!is.null(RM_get(x=bone, RMname=analysis, valuename = "optim")$quantiles))) {
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
      if (is.numeric(analysis)) analysis <- names(RM_list(x=bone, silent = TRUE))[analysis]
      
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
      
      if ((CI == "MCMC") & (is.null(angle)) & (!is.null(RM_get(x=bone, RMname=analysis, valuename = "mcmc")))) {
        polygon(x=c(compactness.synthesis$distance.center, rev(compactness.synthesis$distance.center)), 
                y=c(RM_get(x=bone, RMname=analysis, valuename = "mcmc")$quantiles["2.5%", ], rev(RM_get(x=bone, RMname=analysis, valuename = "mcmc")$quantiles["97.5%", ])), 
                col="lightgrey", border="lightgrey", lwd=3)
      }
      if ((CI == "ML") & (is.null(angle)) & (!is.null(RM_get(x=bone, RMname=analysis, valuename = "optim")$quantiles))) {
        polygon(x=c(compactness.synthesis$distance.center, rev(compactness.synthesis$distance.center)), 
                y=c(RM_get(x=bone, RMname=analysis, valuename = "optim")$quantiles["2.5%", ], rev(RM_get(x=bone, RMname=analysis, valuename = "optim")$quantiles["97.5%", ])), 
                col="lightgrey", border="lightgrey", lwd=3)
      }
      
      lines(compactness.synthesis$distance.center, compactness.synthesis$compactness, lwd=2)
      
      lines(x=compactness.synthesis$distance.center, y=(m+nm)/max(m+nm), col="blue")
      axis(side = 4, at=seq(from=0, to=1, by=0.2), labels = round(seq(from=0, to=1, by=0.2)*max(m+nm), 0), 
           las=1, col.axis="blue", col="blue")
      mtext("Number of pixels", side=4, line=3, col="blue")
      out <- data.frame(distance.center=compactness.synthesis$distance.center, 
                        observed.compactness=compactness.synthesis$compactness)
      
      # p <- c(RM_get(x=bone, RMname=analysis, valuename = "optim")$par, RM_get(x=bone, RMname=analysis, valuename = "optim")$fixed.parameters)
      
      Min <- p["Min"]
      Max <- p["Max"]
      
      # 21/02/2020
      p["S"] <- 1/(4*p["S"])
      
      c <- flexit(x = compactness.synthesis$distance.center, 
                  par = p) * (Max - Min) + Min
      
      out <- cbind(out, modeled.compactness=c)
      
      lines(x=compactness.synthesis$distance.center, y=c, xlim=c(0, 1), lwd=2, lty=3)
      
      
      if (show.legend) {
        if ((CI == "MCMC") & (is.null(angle)) & (!is.null(RM_get(x=bone, RMname=analysis, valuename = "mcmc")))) {
          legend("bottomright", legend=c("Number of pixels", "Observed compactness", "Model", "95% Credibility interval MCMC"), 
                 lty=c(1, 1, 3, 1), lwd=c(1, 2, 2, 6), col=c("blue", "black", "black", "lightgrey"), cex=0.8)
        } else {
          if ((CI == "ML") & (is.null(angle)) & (!is.null(RM_get(x=bone, RMname=analysis, valuename = "optim")$quantiles))) {
            legend("bottomright", legend=c("Number of pixels", "Observed compactness", "Model", "95% Confidence interval ML"), 
                   lty=c(1, 1, 3, 1), lwd=c(1, 2, 2, 6), col=c("blue", "black", "black", "lightgrey"), cex=0.8)
          } else {
            legend("bottomright", legend=c("Number of pixels", "Observed compactness", "Model"), 
                   lty=c(1, 1, 3), lwd=c(1, 2, 2), col=c("blue", "black", "black"), cex=0.8)
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
      
      
      
      if (show.grid & !is.null(RM_get(x=x, RMname=analysis, valuename = "compactness"))) {
        # J'affiche la grille
        compactness <- RM_get(x=x, RMname=analysis, valuename = "compactness")
        
        peripherie <- RM_get(x=x, RMname=analysis, valuename = "peripherie")
        
        l.angle <- levels(compactness$cut.angle)
        angles <- RM_get(x=bone, RMname=analysis, valuename = "cut.angle")
        if (RM_get(x=bone, RMname=analysis, valuename = "partial")) angles <- peripherie[, "angle"]
        l.distance <- levels(compactness$cut.distance.center)
        max.distance <- NULL
        
        for (angle_ec in angles[-1]) {
          
          
          
          angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
          angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
          
          md <- peripherie$peripherie[which.min(abs(peripherie$angle-angle_ec))]
          
          segments(x0=RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"], 
                   x1=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md), 
                   y0=RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"], 
                   y1=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md), 
                   col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.8))
        }
        
        # L'angle 0 deg
        angle_ec <- 0
        angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
        angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
        
        md <- peripherie$peripherie[which.min(abs(peripherie$angle-angle_ec))]*0.9
        text(x=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md*1.2), 
             y=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md*1.2), 
             labels = "0")
        angle_ec <- pi
        angle_ecp <- (angle_ec - RM_get(x=bone, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
        angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
        
        md <- peripherie$peripherie[which.min(abs(peripherie$angle-angle_ec))]*0.9
        text(x=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md*1.2), 
             y=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md*1.2), 
             labels = "pi")
        angle_ec <- pi/2
        angle_ecp <- (angle_ec - RM_get(x=bone, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
        angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
        
        md <- peripherie$peripherie[which.min(abs(peripherie$angle-angle_ec))]*0.9
        text(x=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md*1.2), 
             y=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md*1.2), 
             labels = "pi/2")
        angle_ec <- -pi/2
        angle_ecp <- (angle_ec - RM_get(x=bone, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
        angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
        
        md <- peripherie$peripherie[which.min(abs(peripherie$angle-angle_ec))]*0.9
        text(x=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.x"]+cos(angle_ecp)*md*1.2), 
             y=(RM_get(x=x, RMname=analysis, valuename = "used.centers")["center.y"]+sin(angle_ecp)*md*1.2), 
             labels = "-pi/2")
        
        
        angle_ec <- peripherie$angle
        angle_ecp <- (angle_ec - RM_get(x=x, RMname=analysis, valuename = "rotation.angle" )) %% (2*pi)
        angle_ecp <- ifelse(angle_ecp > pi, -(2*pi)+angle_ecp, angle_ecp)
        
        x_peripherie <- cos(angle_ecp)*peripherie$peripherie
        if (!RM_get(x=bone, RMname=analysis, valuename = "partial")) x_peripherie <- c(x_peripherie, x_peripherie[1])
        y_peripherie <- sin(angle_ecp)*peripherie$peripherie
        if (!RM_get(x=bone, RMname=analysis, valuename = "partial")) y_peripherie <- c(y_peripherie, y_peripherie[1])
        
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


