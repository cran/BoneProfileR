#' BP_AutoFit fits model automatically
#' @title Fit model automatically
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Characteristics of an image with all the fit information
#' @param file The file to be opened
#' @param xlsx TRUE, FALSE or the name and path of the report
#' @param rotation.angle The angle of rotation for analysis
#' @param center Which center to be used.
#' @description Open an image, fit a model and generate a report.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_AutoFit(file=path_Hedgehog, xlsx=TRUE)
#'  # or to open a dialog box
#'  bone <- BP_AutoFit()
#' }
#' @export


BP_AutoFit <- function(file=file.choose(), xlsx=TRUE, 
                       rotation.angle=0, center="ontogenic") {
  
  pdf=FALSE
  docx=FALSE
  
  pb <- txtProgressBar(min = 0, max = 20, initial = 0, style = 3)
  
  bone <-BP_OpenImage(file=file)
  setTxtProgressBar(pb, 1)
  name <- attributes(bone)$name
  
  if (isTRUE(pdf)) {
    pdf.name <- file.path(getwd(), paste0(gsub("\\..+$", "", name), ".pdf"))
  }
  if (isFALSE(pdf)) {
    pdf.name <- NULL
  }
  if (isTRUE(xlsx)) {
    xlsx.name <- file.path(getwd(), paste0(gsub("\\..+$", "", name), ".xlsx"))
  }
  if (isFALSE(xlsx)) {
    xlsx.name <- NULL
  }
  if (isTRUE(docx)) {
    docx.name <- file.path(getwd(), paste0(gsub("\\..+$", "", name), ".docx"))
  }
  if (isFALSE(docx)) {
    docx.name <- NULL
  }
  
  
  
  bone <- BP_DetectBackground(bone=bone, analysis="logistic", show.plot=FALSE)
  setTxtProgressBar(pb, 2)
  bone <- BP_DetectForeground(bone=bone, analysis="logistic", show.plot=FALSE)
  setTxtProgressBar(pb, 3)
  bone <- BP_DetectCenters(bone=bone, analysis="logistic", show.plot=FALSE)
  setTxtProgressBar(pb, 4)
  bone <- BP_EstimateCompactness(bone, analysis="logistic", 
                                 rotation.angle=rotation.angle, 
                                 center=center, show.plot=FALSE)
  setTxtProgressBar(pb, 5)
  bone <- BP_FitMLCompactness(bone, analysis="logistic", silent=TRUE, 
                              fixed.parameters = c(K1=1, K2=1, Max=3, Min=-3), 
                              fitted.parameters = c(P=0.5, S=0.1))
  fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
  bone <- BP_FitMLCompactness(bone, analysis="logistic", silent=TRUE, 
                              fixed.parameters = c(K1=1, K2=1), 
                              fitted.parameters = c(fittedpar, Max=2, Min=-2))
  setTxtProgressBar(pb, 8)
  fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
  setTxtProgressBar(pb, 9)
  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
  setTxtProgressBar(pb, 10)
  bone <- BP_FitMLCompactness(bone, 
                              fitted.parameters=c(fittedpar, K1=1, K2=1), 
                              fixed.parameters=NULL, analysis="flexit", silent=TRUE)
  setTxtProgressBar(pb, 14)
  outAIC <- compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", alloptim=TRUE), 
                        Flexit=BP_GetFittedParameters(bone, analysis="flexit", alloptim=TRUE), silent = TRUE)
  
  if (outAIC$DeltaAIC[1]==0) {
    # Model Logistic
    bone <- RM_delete(bone, RMname ="flexit")
  } else {
    # Model Flexit
    bone <- RM_delete(bone, RMname ="logistic")
  }
  setTxtProgressBar(pb, 15)
  bone <- BP_FitBayesianCompactness(bone, analysis=1)
  setTxtProgressBar(pb, 18)
  bone <- BP_FitMLRadialCompactness(bone, analysis=1, silent=TRUE)
  setTxtProgressBar(pb, 20)
  if (xlsx | docx | pdf) BP_Report(bone=bone, 
            analysis=1,
            control.plot = list(message = NULL, show.centers = TRUE, show.colors = TRUE,
                                show.grid = TRUE, CI = "MCMC", show.legend = TRUE), 
            docx=docx.name, 
            pdf= pdf.name, 
            xlsx=xlsx.name, 
            author="Marc Girondot", 
            title=name)
  
  return(invisible(bone))
}

