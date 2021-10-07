# #' .BP_contour estimate a distance matrix

.BP_distance <- function(bone, threshold, analysis=1, center.x=NA, center.y=NA) {
  
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
  # La colonne 1 va de 1 à 34; C'est x
  # La colonne 2 va de 1 à 24; c'est y
  d <- sapply(1:nrow(m), FUN = function(r) sqrt((m[r, 1]-center.x)^2+(m[r, 2]-center.y)^2))
  
  # d[2] c'est x=1, y=2
  # d[3] c'est x=1, y=3
  
  d_mat <- matrix(d, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = TRUE)
  return(d_mat)

}
