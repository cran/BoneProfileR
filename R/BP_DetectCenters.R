#' BP_DetectCenters detects the centers of an image
#' @title Detect the centers of an image
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The orignial bone object with a new attribute for centers
#' @param bone The bone image to be used
#' @param analysis The name or rank of analysis
#' @param method Can be Fast, Accurate, FastConvex, or AccurateConvex
#' @param show.plot should plot is shown ?
#' @description Detects the centers of an image. Note that this function must not be used with partial bone section.\cr
#' The method Fast works well with the convex bone section while if the section is concave, Accurate is slower but works well in all circonstances.\cr
#' Fast method is maintained here only for compatibility with versions <3.1 of BoneProfileR.\cr
#' If the section is concave, the methods FastConvex and AccurateConvex return a minimum convex section.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#'  bone <- BP_OpenImage()
#'  # or 
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_DetectBackground(bone=bone)
#'  bone <- BP_DetectForeground(bone=bone)
#'  bone <- BP_DetectCenters(bone=bone)
#'  plot(bone, type="mineralized", show.grid=FALSE)
#'  plot(bone, type="unmineralized", show.grid=FALSE)
#'  plot(bone, type="section", show.grid=FALSE)
#'  # Note that some parts of the section are concave but it does not give problems in the analysis
#'  # For section with very strong concavity, it could be safer to use:
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic", method="AccurateConvex")
#'  plot(bone, type="mineralized", show.grid=FALSE)
#'  plot(bone, type="unmineralized", show.grid=FALSE)
#'  plot(bone, type="section", show.grid=FALSE)
#' }
#' @export


BP_DetectCenters <- function(bone, analysis=1, show.plot=TRUE, method="Accurate") {
  
  if (is.null(RM_get(x=bone, RMname=analysis, valuename = "bg")) | 
      is.null(RM_get(x=bone, RMname=analysis, valuename = "fg"))) {
    stop("You must first setup background and foreground colors") 
  }
  
  bg <- RM_get(x=bone, RMname=analysis, valuename = "bg")
  fg <- RM_get(x=bone, RMname=analysis, valuename = "fg")
  
  # Je formatte la coupe en threshold
  # if (is.null(RM_get(x=bone, RMname=analysis, valuename = "threshold"))) {
    threshold <- getFromNamespace(".BP_threshold", ns="BoneProfileR")(bone)
    bone <- RM_add(x=bone, RMname = analysis, valuename="threshold", 
                   value=threshold)
 # } else {
 #   threshold <- RM_get(x=bone, RMname = analysis, valuename="threshold")
#  }
  
#  contour <- RM_get(x=bone, RMname = analysis, valuename="contour")
 #  if (is.null(contour)){
  contour <- getFromNamespace(".BP_contour", ns="BoneProfileR")(bone, 
                                                                threshold=threshold, 
                                                                analysis=analysis, 
                                                                partial=FALSE, 
                                                               center.x=NA, 
                                                               center.y=NA, 
                                                               method=method)
  bone <- RM_add(x=bone, RMname = analysis, valuename="contour", 
                 value=contour)
 # }
  
  # essai <- array(data=as.numeric(contour), dim=c(dim(contour), 1, 1))
  # class(essai) <- c("BoneProfileR", "cimg", "imager_array", "numeric" )
  # plot(essai)
  
  GC_cortex.x <- mean(which(threshold, arr.ind = TRUE)[, 1])
  GC_cortex.y <- mean(which(threshold, arr.ind = TRUE)[, 2])
  
  GC_bone.x <- mean(which(contour, arr.ind = TRUE)[, 1])
  GC_bone.y <- mean(which(contour, arr.ind = TRUE)[, 2])
  
  GC_medula.x <- mean(which(contour & !threshold, arr.ind = TRUE)[, 1])
  GC_medula.y <- mean(which(contour & !threshold, arr.ind = TRUE)[, 2])
  
  bone <- RM_delete(x=bone, RMname = analysis, valuename="compactness")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="array.compactness") 
  bone <- RM_delete(x=bone, RMname = analysis, valuename="cut.distance.center") 
  bone <- RM_delete(x=bone, RMname = analysis, valuename="cut.angle")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="compactness.synthesis")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optim")
  # bone <- RM_delete(x=bone, RMname = analysis, valuename="contour")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="used.centers")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimRadial")
  
  # 
  
  GC_ontoCenter.x <- GC_bone.x - (GC_cortex.x - GC_bone.x) + (GC_medula.x - GC_bone.x)
  GC_ontoCenter.y <- GC_bone.y - (GC_cortex.y - GC_bone.y) + (GC_medula.y - GC_bone.y)
  # points(GC_ontoCenter.x, GC_ontoCenter.y, col="blue", pch=19)
  
  bone <- RM_add(x=bone, RMname = analysis, valuename="centers", value=c(GC_cortex.x=GC_cortex.x, 
                                                                         GC_cortex.y=GC_cortex.y, 
                                                                         GC_bone.x=GC_bone.x, 
                                                                         GC_bone.y=GC_bone.y, 
                                                                         GC_medula.x=GC_medula.x, 
                                                                         GC_medula.y=GC_medula.y, 
                                                                         GC_user.x=NA, 
                                                                         GC_user.y=NA, 
                                                                         GC_ontogenic.x=GC_ontoCenter.x, 
                                                                         GC_ontogenic.y=GC_ontoCenter.y))
  bone <- RM_add(x=bone, RMname = analysis, valuename="method", value=method)
  
  if (show.plot) plot(bone, message="Do not forget to check thresholding")
  return(bone)
}
