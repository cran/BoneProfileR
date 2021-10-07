#' BP_ListAnalyses lists the analyses stored in an object
#' @title List the analyses stored in an object
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The list of analyses
#' @param bone The bone image to be used
#' @param silent Should the results be shown ?
#' @param max.level If TRUE, will return all list element of the objects
#' @description Get the analyses stored in an object.
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
#'  plot(bone)
#'  plot(bone, type="observations")
#'  plot(bone, type="observations+model", analysis=1)
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  BP_ListAnalyses(bone)
#' }
#' @export


BP_ListAnalyses <- function(bone, silent=TRUE, max.level = FALSE) {
  out <- RM_list(x=bone, silent=silent, max.level = max.level)
  return(out)
}

