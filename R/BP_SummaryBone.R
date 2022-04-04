#' summary.BoneProfileR displays a bone section
#' @title Plot a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return An invisible list with recorded information
#' @param object The bone image
#' @param max.level	If TRUE, will return all list element of the objects
#' @param analysis The analysis to report the global compacteness
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


summary.BoneProfileR <- function(object, max.level=FALSE, analysis=1, ...) {
  out <- dim(object)
  cat(paste0("The image is ", as.character(out[1]), " pixels width and ", as.character(out[2]), " pixels height.\n"))
  
  if (!is.null(RM_get(x=object, RMname = analysis, valuename = "global.compactness")))
  cat(paste0("The observed global compactness for ", analysis, " analysis is ", 
             specify_decimal(RM_get(x=object, RMname = analysis, valuename = "global.compactness"), decimals = 3), ".\n"))
  if (!is.null(RM_get(x=object, RMname = analysis, valuename = "optim")$quantiles["50%", ])) {
    cat(paste0("The median value for ", analysis, " analysis for modeled global compactness is ", 
               specify_decimal(mean(RM_get(x=object, RMname = analysis, valuename = "optim")$quantiles["50%", ]), decimals = 3), ".\n"))
    cat(paste0("The 95% confidence interval for ", analysis, " analysis for modeled global compactness is between ", 
               specify_decimal(mean(RM_get(x=object, RMname = analysis, valuename = "optim")$quantiles["2.5%", ]), decimals = 3), " and ",
               specify_decimal(mean(RM_get(x=object, RMname = analysis, valuename = "optim")$quantiles["97.5%", ]), decimals = 3),".\n"))
  }
  
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


