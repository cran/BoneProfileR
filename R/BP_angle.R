# #' .BP_contour estimate a distance matrix

.BP_angle <- function(bone, threshold, analysis=1, center.x=NA, center.y=NA) {
  
  # analysis=1; partial=FALSE; center.x=NA; center.y=NA
  if (is.null(threshold)) threshold <- RM_get(x=bone, RMname = analysis, valuename="threshold")
  
  if (is.null(threshold)) {
    threshold <- getFromNamespace(".BP_threshold", ns="BoneProfileR")(bone, analysis=analysis)
  }
  
  # Dans bone et threshold, le x est en dimension 1 et le y en dimension 2
  if (is.na(center.x) | is.na(center.y)) {
    center.x <- mean(which(threshold, arr.ind = TRUE)[, 1])
    center.y <- mean(which(threshold, arr.ind = TRUE)[, 2])
  }
  
  # Je dois travailler avec le center
  # Attention, dans l'objet bone, row c'est x et col c'est y
  
  # La colonne 1 est y; je fais varier d'abord les y
  # La colonne 2 est x
  m <- expand.grid(1:dim(threshold)[2], 1:dim(threshold)[1])
  a <- sapply(1:nrow(m), FUN = function(r) {
    x <- m[r, 1]
    y <- m[r, 2]
    
    angle <- atan2(y-center.y, x-center.x) %% (2*pi)
    # angle <- ((angle+pi+rotation.angle) %% (2*pi))-pi
    return(angle)
    # return((atan((x-center.x)/(y-center.y))+pi*(as.numeric((y-center.y)<0))) %% (2*pi))
  })
  
  a_mat <- matrix(a, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = TRUE)
  return(a_mat)

}
