# #' .BP_threshold estimate a countour matrix

.BP_threshold <- function(bone, analysis=1) {
  
  if (is.null(RM_get(x=bone, RMname=analysis, valuename = "bg")) | 
      is.null(RM_get(x=bone, RMname=analysis, valuename = "fg"))) {
    stop("You must first setup background and foreground colors") 
  }
  
  bg <- RM_get(x=bone, RMname=analysis, valuename = "bg")
  fg <- RM_get(x=bone, RMname=analysis, valuename = "fg")
  
  bg_red <- col2rgb(bg)["red", 1]/255
  bg_green <- col2rgb(bg)["green", 1]/255
  bg_blue <- col2rgb(bg)["blue", 1]/255
  
  fg_red <- col2rgb(fg)["red", 1]/255
  fg_green <- col2rgb(fg)["green", 1]/255
  fg_blue <- col2rgb(fg)["blue", 1]/255
  
  
  Distance_bg <- sqrt((bone[, , 1, 1]-bg_red)^2+
                        (bone[, , 1, 2]-bg_green)^2+
                        (bone[, , 1, 3]-bg_blue)^2)
  Distance_fg <- sqrt((bone[, , 1, 1]-fg_red)^2+
                        (bone[, , 1, 2]-fg_green)^2+
                        (bone[, , 1, 3]-fg_blue)^2)
  
  return(Distance_bg > Distance_fg)
}
