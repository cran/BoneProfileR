# #' .BP_threshold estimate a countour matrix

.BP_threshold <- function(bone, analysis=1) {
  
  if (is.null(RM_get(x=bone, RMname=analysis, valuename = "bg")) | 
      is.null(RM_get(x=bone, RMname=analysis, valuename = "fg"))) {
    stop("You must first setup background and foreground colors") 
  }
  
  bg <- RM_get(x=bone, RMname=analysis, valuename = "bg")
  fg <- RM_get(x=bone, RMname=analysis, valuename = "fg")
  
  Distance_bg <- sqrt((bone[, , 1, 1]-col2rgb(bg)["red", 1]/255)^2+
                        (bone[, , 1, 2]-col2rgb(bg)["green", 1]/255)^2+
                        (bone[, , 1, 3]-col2rgb(bg)["blue", 1]/255)^2)
  Distance_fg <- sqrt((bone[, , 1, 1]-col2rgb(fg)["red", 1]/255)^2+
                        (bone[, , 1, 2]-col2rgb(fg)["green", 1]/255)^2+
                        (bone[, , 1, 3]-col2rgb(fg)["blue", 1]/255)^2)
  
  
  return(Distance_bg>Distance_fg)
}
