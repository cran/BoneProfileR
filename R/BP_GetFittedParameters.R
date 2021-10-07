#' BP_GetFittedParameters returns the fitted parameters
#' @title Return the fitted parameters
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The fitted parameters
#' @param bone The bone image to be used
#' @param analysis Name or rank of analysis
#' @param alloptim If TRUE, return the complete object returned by optim
#' @description Return the fitted parameters.
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


BP_GetFittedParameters <- function(bone, analysis=1, alloptim=FALSE) {
  out <- RM_get(x=bone, RMname=analysis, valuename = "optim")
  if (!alloptim) {
    out <- out$par
  }
  return(out)
}

