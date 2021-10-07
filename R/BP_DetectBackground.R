#' BP_DetectBackground detects the background color of an image
#' @title Detects the background color of an image
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The orignial bone object with a new attribute for background color
#' @param bone The bone image to be used
#' @param analysis The name or rank of analysis
#' @param show.plot should plot is shown ?
#' @description Detects the background color of an image.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#'  bone <- BP_OpenImage()
#'  bone <- BP_DetectBackground(bone=bone)
#'  bone <- BP_DetectForeground(bone=bone)
#'  plot(bone)
#' }
#' @export


BP_DetectBackground <- function(bone, analysis=1, show.plot=TRUE) {
  red <- bone[c(2:5, dim(bone)[1]-2:5), c(2:5, dim(bone)[2]-2:5), 1, 1]
  green <- bone[c(2:5, dim(bone)[1]-2:5), c(2:5, dim(bone)[2]-2:5), 1, 2]
  blue <- bone[c(2:5, dim(bone)[1]-2:5), c(2:5, dim(bone)[2]-2:5), 1, 3]
  bg <- rgb(red=median(red), green=median(green), blue=median(blue))
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
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optim")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="used.centers")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimRadial")
  
  if (show.plot) plot(bone)
  return(bone)
}

