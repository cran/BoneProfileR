#' A model for bone compactness.
#'
#' \tabular{ll}{
#'  Package: \tab BoneProfileR\cr
#'  Type: \tab Package\cr
#'  Version: \tab 2.5 build 773\cr
#'  Date: \tab 2024-04-10\cr
#'  License: \tab GPL (>= 2)\cr
#'  LazyLoad: \tab yes\cr
#'  }
#' @title A Model for Bone Compactness.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @docType package
#' @name BoneProfileR-package
#' @description A Model for Bone Compactness.\cr
#' The lastest version of this package can always been installed using:\cr
#' install.packages(c("imager", "tiff", "ijtiff", "HelpersMG", "knitr", "rmarkdown", "openxlsx", "shiny"))\cr
#' install.packages("https://hebergement.universite-paris-saclay.fr/marcgirondot/CRAN/HelpersMG.tar.gz", repos=NULL, type="source")\cr
#' install.packages("https://hebergement.universite-paris-saclay.fr/marcgirondot/CRAN/BoneProfileR.tar.gz", repos=NULL, type="source")\cr
#' BoneProfileR uses a new results management software that is developed as part of the HelpersMG 
#' package. Using this results management system (RM), all the results are stored as part of the 
#' analyzed image.\cr
#' This results management software has been developed to help users to maintain the results 
#' associated with the methodology used to obtain it. It is part of the large movement in 
#' science of replicative research.\cr
#' An analysis is then stored with the image in a single file with the following information:\cr
#' name, timestamp, bg, fg, threshold, contour, centers, peripherie, compactness, 
#' array.compactness, cut.distance.center, cut.angle, used.centers, compactness.synthesis, 
#' partial, rotation.angle, global.compactness, optim, optimRadial\cr
#' Several analyses can be stored within a single file.\cr
#' \if{html}{\figure{E.png}{options: alt="BoneProfileR logo"}}
#' \if{latex}{\figure{E.png}}
#' @references Girondot M, Laurin M (2003) Bone Profiler: a tool to quantify, model, and 
#' statistically compare bone-section compactness profiles. Journal of 
#' Vertebrate Paleontology 23: 458-461
#' @references Laurin M, Girondot M, Loth M-M (2004) The evolution of long bone microstructure 
#' and lifestyle in lissamphibians. Paleobiology 30: 589-613 
#' @references Gônet, Jordan, Michel Laurin, and Marc Girondot. 2022. BoneProfileR: 
#' The Next Step to Quantify, Model and Statistically Compare Bone 
#' Section Compactness Profiles. Paleontologica Electronica. 25(1): a12
#' @references Gônet, J., Bardin, J., Girondot, M., Hutchinson, J., Laurin, M., (2023). Deciphering 
#' locomotion in reptiles: application of elliptic Fourier transforms to femoral microanatomy. 
#' Zoological Journal of the Linnean Society 198, 1070-1091.
#' @references Gônet, J., Bardin, J., Girondot, M., Hutchinson, J.R., Laurin, M., (2023). Locomotor 
#' and postural diversity among reptiles viewed through the prism of femoral microanatomy: 
#' palaeobiological implications for some Permian and Mesozoic taxa. Journal of Anatomy 242, 891-916.
#' @references Gônet, J., Bardin, J., Girondot, M., Hutchinson, J.R., Laurin, M., (2023). Unravelling 
#' the postural diversity of mammals: contribution of humeral cross-sections to palaeobiological 
#' inferences. Journal of Mamalian Evolution 30, 321-337.
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  plot(bone, type="original")
#'  bone <- BP_DetectBackground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectForeground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic")
#'  plot(bone, type="original")
#'  plot(bone, type="mineralized")
#'  plot(bone, type="unmineralized")
#'  plot(bone, type="section")
#'  plot(bone, type="colors")
#'  plot(bone, type="3Dcolors")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic", center="ontogenetic")
#'  plot(bone, type="original")
#'  plot(bone, type="mineralized")
#'  plot(bone, type="observations")
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic")
#'  plot(bone, type="model", analysis=1)
#'  plot(bone, type="observations+model", analysis=1)
#'  fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  bone <- BP_FitMLCompactness(bone, 
#'                 fitted.parameters=c(fittedpar, K1=1, K2=1), 
#'                 fixed.parameters=NULL, analysis="flexit")
#'  compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", alloptim=TRUE), 
#'              Flexit=BP_GetFittedParameters(bone, analysis="flexit", alloptim=TRUE))
#' # pdf(file = "Figure 2.pdf", width = 8, height = 10, pointsize = 12)
#' layout(1:2)
#' plot(bone, type="observations+model", analysis="logistic", restorePar=FALSE, mar=c(4, 4, 2, 5))
#' plot(bone, type="observations+model", analysis="flexit", restorePar=FALSE, mar=c(4, 4, 2, 5))
#' layout(1)
#' # dev.off()
#'  out4p <- plot(bone, type="observations+model", analysis="logistic")
#'  out6p <- plot(bone, type="observations+model", analysis="flexit")
#'  bone <- BP_FitBayesianCompactness(bone, analysis="logistic")
#'  plot(bone, type="observations+model", CI="MCMC")
#'  bone <- BP_FitBayesianCompactness(bone, analysis="flexit")
#'  plot(bone, type="observations+model", CI="MCMC", analysis="flexit")
#'  plot(bone, type="mcmc", parameter="P", 
#'       options.mcmc=list(xlim=c(0.55, 0.57), breaks=seq(from=0, to=1, by=0.001)))
#'  plot(bone, type="mcmc", parameter="S", 
#'       options.mcmc=list(xlim=c(0.02, 0.05), breaks=seq(from=0.02, to=.05, by=0.001)))
#'  plot(bone, type="mcmc", parameter="Min", 
#'       options.mcmc=list(xlim=c(0.05, 0.08), breaks=seq(from=0, to=1, by=0.001)))
#'  plot(bone, type="mcmc", parameter="Max", 
#'       options.mcmc=list(xlim=c(0.95, 0.97), breaks=seq(from=0, to=1, by=0.001)))
#'  outMCMC <- RM_get(x = bone, RM = "RM", RMname = "logistic", valuename = "mcmc")
#'  summary(outMCMC)
#'  outMCMC <- RM_get(x = bone, RM = "RM", RMname = "flexit", valuename = "mcmc")
#'  summary(outMCMC)
#'  # pdf(file = "Figure 3.pdf", width = 8, height = 10, pointsize = 12)
#'  layout(1:2)
#'  plot(bone, type="mcmc", parameter="K1", analysis="flexit", 
#'       options.mcmc=list(xlim=c(-1, 3), ylim=c(0,10), 
#'       breaks=seq(from=-1, to=3, by=0.001), 
#'       legend = FALSE, show.prior = FALSE, mar=c(4, 4, 1, 6)), restorePar=FALSE)
#'  segments(x0=1, x1=1, 
#'          y0=0, y1=10, lty=4, lwd=3)
#'  text(x=ScalePreviousPlot(x=0.95, y=0.95)$x, 
#'       y=ScalePreviousPlot(x=0.95, y=0.95)$y, labels="A", cex=3)
#'  plot(bone, type="mcmc", parameter="K2", analysis="flexit", 
#'       options.mcmc=list(xlim=c(-1, 3), ylim=c(0,10), 
#'       breaks=seq(from=-1, to=3, by=0.001), 
#'       legend = FALSE, show.prior = FALSE, mar=c(4, 4, 1, 6)), restorePar=FALSE)
#'  segments(x0=1, x1=1, 
#'          y0=0, y1=10, lty=4, lwd=3)
#'  text(x=ScalePreviousPlot(x=0.95, y=0.95)$x, 
#'       y=ScalePreviousPlot(x=0.95, y=0.95)$y, labels="B", cex=3)
#'  # dev.off()
#'  
#'  bone <- BP_FitMLRadialCompactness(bone, analysis = "flexit")
#'  plot(bone, type="radial", radial.variable=c("P", "S"), analysis = "flexit")
#'  plot(bone, type="radial", radial.variable=c("P", "S", "Min", "Max"), analysis = "flexit")
#'  out <- RM_get(x=bone, RMname="flexit", valuename = "optimRadial")$synthesis
#'  mean(out[, "P"]); sd(out[, "P"])
#'  range(out[, "S"])
#'  quantile(out[, "S"])
#'  # pdf(file = "Figure 4.pdf", width=7, height = 9, pointsize = 12)
#'  layout(1:2)
#'  plot(bone, type="radial", radial.variable="P", analysis = "flexit", restorePar=FALSE)
#'  text(x=ScalePreviousPlot(x=0.95, y=0.95)$x, 
#'       y=ScalePreviousPlot(x=0.95, y=0.95)$y, labels="A", cex=3)
#'  plot(bone, type="radial", radial.variable="S", analysis = "flexit", restorePar=FALSE)
#'  text(x=ScalePreviousPlot(x=0.95, y=0.95)$x, 
#'       y=ScalePreviousPlot(x=0.95, y=0.95)$y, labels="B", cex=3)
#'  # dev.off()
#'  #' # How many times this package has been download
#' library(cranlogs)
#' BoneProfileR <- cran_downloads("BoneProfileR", from = "2021-10-07", 
#'                             to = Sys.Date() - 1) 
#' sum(BoneProfileR$count)
#' plot(BoneProfileR$date, BoneProfileR$count, type="l", bty="n", 
#'      xlab="Download date", ylab="Number of downloads")
#' }

NULL

