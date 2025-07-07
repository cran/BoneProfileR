#' BP_GetFittedParameters returns the fitted parameters
#' @title Return the fitted parameters
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The fitted parameters
#' @param bone The bone image to be used
#' @param analysis Name or rank of analysis
#' @param periodic If TRUE, the periodic model is used (depreciated, use type)
#' @param return.all If TRUE, return the complete object
#' @param ML If TRUE, return the ML estimate and the SE ; if FALSE, returns the MCMC estimate
#' @param type Can be "global", "radial", or "periodic"
#' @description Return the fitted parameters or complete object.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_DetectBackground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectForeground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic")
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic")
#'  BP_GetFittedParameters(bone, analysis="logistic")
#' }
#' @export


BP_GetFittedParameters <- function(bone, 
                                   analysis=1, 
                                   return.all=FALSE, 
                                   ML=TRUE, 
                                   periodic=FALSE, 
                                   type="global") {

  type <- tolower(type)
  type <- match.arg(type, choices = c("global", "radial", "periodic"))
  
  if (periodic & (type != "periodic")) stop("When periodic is TRUE, type must be equal to 'periodic'")
  
  if (periodic | (type == "periodic")) {
    if (ML) {
      out <- RM_get(x=bone, RMname=analysis, valuename = "optimPeriodic")
      if (is.null(out)) stop(paste0("The periodic analysis has not been done with analysis ", as.character(analysis),"."))
      if (!return.all) {
        out <- as.matrix(data.frame(mean=out$par, se=out$SE))
      }
    } else {
      out <- RM_get(x=bone, RMname = analysis, valuename="mcmcPeriodic")
      if (is.null(out)) stop(paste0("The MCMC periodic analysis has not been done with analysis ", as.character(analysis),"."))
      if (!return.all) {
        out <- as.matrix(cbind(out$summary.table, t(as.parameters(out, index = "quantile"))))
      }
    }
  }
  
  if (type == "global") {
    if (ML) {
      out <- RM_get(x=bone, RMname=analysis, valuename = "optim")
      if (!return.all) {
        out <- as.matrix(data.frame(mean=out$par, se=out$SE))
      }
    } else {
      out <- RM_get(x=bone, RMname=analysis, valuename = "mcmc")
      if (is.null(out)) stop(paste0("The MCMC analysis has not been done with analysis ", as.character(analysis),"."))
      if (!return.all) {
        out <- as.matrix(cbind(out$summary.table, t(as.parameters(out, index = "quantile"))))
      }
    }
  }
  
  if (type == "radial") {
    if (ML) {
      out <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")
      if (is.null(out)) stop(paste0("The radial analysis has not been done with analysis ", as.character(analysis),"."))
      if (!return.all) {
        out <- as.matrix(out$summary.table)
      }
    } else {
      stop("The Bayesian version of the radial analysis is not available. Use Bayesian periodic analysis instead.")
    }
  }
  return(out)
}

