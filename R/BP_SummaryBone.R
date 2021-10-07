#' summary.BoneProfileR displays a bone section
#' @title Plot a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return An invisible list with recorded information
#' @param object The bone image
#' @param max.level	If TRUE, will return all list element of the objects
#' @param ... Not used
#' @description Display information of bone section
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
#'  summary(bone)
#' }
#' @method summary BoneProfileR
#' @export


summary.BoneProfileR <- function(object, max.level=FALSE, ...) {
  out <- dim(object)
  cat(paste0("The image is ", as.character(out[1]), " pixels width and ", as.character(out[2]), " pixels height.\n"))
  
  an <- BP_ListAnalyses(object, max.level=max.level, silent = TRUE)
  if (length(an) == 0) {
    cat("There are no recorded analysis still.\n")
  } else {
  cat("The current recorded analysis are:\n")
  an <- BP_ListAnalyses(object, max.level=max.level, silent = FALSE)
  }
  out <- list(dim=out, analysis=an)
  return(invisible(out))
}


