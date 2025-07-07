#' BP_DetectForeground detects the foreground color of an image
#' @title Detects the foreground color of an image
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The orignial bone object with a new attribute for foreground color
#' @param bone The bone image to be used
#' @param analysis The name or rank of analysis
#' @param show.plot should plot is shown ?
#' @description Detects the foreground color of an image.
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


BP_DetectForeground <- function(bone, analysis=1, show.plot=TRUE) {
  if (is.null(RM_get(x=bone, RM="RM", RMname=analysis, valuename = "bg"))) {
    red <- bone[c(2:5, dim(bone)[1]-2:5), c(2:5, dim(bone)[2]-2:5), 1, 1]
    green <- bone[c(2:5, dim(bone)[1]-2:5), c(2:5, dim(bone)[2]-2:5), 1, 2]
    blue <- bone[c(2:5, dim(bone)[1]-2:5), c(2:5, dim(bone)[2]-2:5), 1, 3]
    bg <- rgb(red=median(red), green=median(green), blue=median(blue))
    red <- col2rgb(bg)["red", 1]
    green <- col2rgb(bg)["green", 1]
    blue <- col2rgb(bg)["blue", 1]
  } else {
  bg <- RM_get(x=bone, RM="RM", RMname=analysis, valuename = "bg")
  red <- col2rgb(bg)["red", 1]
  green <- col2rgb(bg)["green", 1]
  blue <- col2rgb(bg)["blue", 1]
}
  pos <- sqrt((bone[, , 1, 1]-red/255)^2+
                          (bone[, , 1, 2]-green/255)^2+
                          (bone[, , 1, 3]-blue/255)^2)
  pos2 <- which(pos==max(pos), arr.ind = TRUE)[1, ]
  
  fg <- rgb(red=bone[pos2[1], pos2[2], 1, 1], 
            green=bone[pos2[1], pos2[2], 1, 2], 
            blue=bone[pos2[1], pos2[2], 1, 3])
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
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimPeriodic")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="used.centers")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimRadial")
  
  if (show.plot) plot(bone, message="Do not forget to check thresholding", analysis = analysis)
  
    return(bone)
}

