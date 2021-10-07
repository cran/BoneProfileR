#' BP_ChooseCenter lets the use to choose the center of the bone
#' @title Let the user to choose the center of the bone
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The orignal bone object with a new attribute for center
#' @param bone The bone image to be used
#' @param analysis The name or rank of analysis
#' @description Let the user to choose the center of the bone.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_DetectBackground(bone=bone)
#'  bone <- BP_DetectForeground(bone=bone)
#'  bone <- BP_ChooseCenter(bone=bone)
#'  # For partial section, only BP_ChooseCenter() must be used
#'  path_Dicynodon <- system.file("extdata", "Dicynodon_tibia_11.11.1.T_b_b-1.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Dicynodon)
#'  bone <- BP_DetectBackground(bone=bone)
#'  bone <- BP_DetectForeground(bone=bone)
#'  bone <- BP_ChooseCenter(bone=bone)
#'  bone <- BP_EstimateCompactness(bone, center="user", partial=TRUE)
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic")
#'  plot(bone, type="observations+model")
#' }
#' @export


BP_ChooseCenter <- function(bone, analysis=1) {
  
  if (is.null(RM_get(x=bone, RMname=analysis, valuename = "bg")) | 
      is.null(RM_get(x=bone, RMname=analysis, valuename = "fg"))) {
    stop("You must first setup background and foreground colors") 
  }
  
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  
  plot(bone, 
       message="Please choose the center of the section", restorePar=FALSE)
  pos <- getFromNamespace(".BP_DetectClick", ns="BoneProfileR")(bone)
  
  GC_cortex.x <- NA
  GC_cortex.y <- NA
  
  GC_bone.x <- NA
  GC_bone.y <- NA
  
  GC_medula.x <- NA
  GC_medula.y <- NA
  
  bone <- RM_delete(x=bone, RMname = analysis, valuename="compactness")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="array.compactness") 
  bone <- RM_delete(x=bone, RMname = analysis, valuename="cut.distance.center") 
  bone <- RM_delete(x=bone, RMname = analysis, valuename="cut.angle")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="compactness.synthesis")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optim")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="used.centers")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimRadial")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="contour")
  
  bone <- RM_add(x=bone, RMname = analysis, valuename="centers", value=c(GC_cortex.x=GC_cortex.x, 
                                                                         GC_cortex.y=GC_cortex.y, 
                                                                         GC_bone.x=GC_bone.x, 
                                                                         GC_bone.y=GC_bone.y, 
                                                                         GC_medula.x=GC_medula.x, 
                                                                         GC_medula.y=GC_medula.y, 
                                                                         GC_user.x=unname(pos["x"]), 
                                                                         GC_user.y=unname(pos["y"])))
  
  bone <- RM_add(x=bone, RMname = analysis, valuename="used.centers", value=c(center.x=unname(pos["x"]), center.y=unname(pos["y"])))
  
  plot(bone)
  return(bone)
}

