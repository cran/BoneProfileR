#' BP_ChooseForeground let the user to choose the foreground color of an image
#' @title Let the user to choose the foreground color of an image
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return The orignial bone object with a new attribute for foreground color
#' @param bone The bone image to be used
#' @param analysis The name or rank of analysis
#' @description Let the user to choose the foreground color of an image.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#'  bone <- BP_OpenImage()
#'  bone <- BP_ChooseBackground(bone=bone)
#'  bone <- BP_ChooseForeground(bone=bone)
#'  plot(bone)
#' }
#' @export


BP_ChooseForeground <- function(bone, analysis=1) {
  
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  
  plot(bone, 
              message="Please choose the foreground color", restorePar=FALSE)
  pos <- getFromNamespace(".BP_DetectClick", ns="BoneProfileR")(bone)
  fg <- bone[pos["x"], pos["y"], 1, 1:3]
  fg <- rgb(red=fg[1], green=fg[2], blue=fg[3])
  bone <- RM_add(x=bone, RMname=analysis, 
                 valuename = "fg", value=fg)
  
  bone <- RM_delete(x=bone, RMname = analysis, valuename="threshold")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="contour")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="centers")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="compactness")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="array.compactness") 
  bone <- RM_delete(x=bone, RMname = analysis, valuename="cut.distance.center") 
  bone <- RM_delete(x=bone, RMname = analysis, valuename="cut.angle")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="compactness.synthesis")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optim")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="used.centers")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimRadial")
  
  
  plot(bone, message="Do not forget to check thresholding")
  return(bone)
}

