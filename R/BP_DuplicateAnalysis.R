#' BP_DuplicateAnalysis duplicates an analysis stored in an object
#' @title Duplicates an analysis stored in an object
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The orignial bone object with a new analysis
#' @param bone The bone image to be used
#' @param from The name or rank of analysis to be duplicated
#' @param to The name or rank of analysis to be created
#' @description Duplicates an analysis stored in an object.
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
#' }
#' @export


BP_DuplicateAnalysis <- function(bone, from=1, to=2) {
  
  bone <- RM_duplicate(x=bone, RMnamefrom = from, RMnameto = to)
  return(bone)
}

