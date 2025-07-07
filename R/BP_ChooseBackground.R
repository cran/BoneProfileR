#' BP_ChooseBackground lets the use to choose the background color of an image
#' @title Let the use to choose the background color of an image
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The orignial bone object with a new attribute for background color
#' @param bone The bone image to be used
#' @param analysis The name or rank of analysis
#' @description Let the user to choose the background color of an image.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#'  path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.tif", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_ChooseBackground(bone=bone)
#'  bone <- BP_ChooseForeground(bone=bone)
#'  plot(bone)
#'  }
#' @export


BP_ChooseBackground <- function(bone, analysis=1) {
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  
  plot(bone, 
       message="Please choose the background color", restorePar=FALSE, analysis = analysis)
  pos <- getFromNamespace(".BP_DetectClick", ns="BoneProfileR")(bone)
  bg <- bone[pos["x"], pos["y"], 1, 1:3]
  bg <- rgb(red=bg[1], green=bg[2], blue=bg[3])
  
  bone <- RM_add(x=bone, RMname=analysis, 
                 valuename = "bg", value=bg)
  
  bone <- RM_delete(x=bone, RMname = analysis, valuename="threshold")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="contour")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="centers")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="compactness")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="array.compactness") 
  bone <- RM_delete(x=bone, RMname = analysis, valuename="cut.distance.center") 
  bone <- RM_delete(x=bone, RMname = analysis, valuename="cut.angle")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="compactness.synthesis")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimPeriodic")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optim")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="used.centers")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimRadial")
  
  plot(bone, message="Do not forget to check thresholding", analysis = analysis)
  return(bone)
}

