# #' .BP_contour estimate a contour matrix

.BP_contour <- function(bone               , 
                        threshold=NULL     , 
                        analysis=1         , 
                        partial=FALSE      , 
                        center.x=NA        , 
                        center.y=NA        , 
                        method="accurate"  ) {
  
  # method can be fast, accurate, fastconvex, or accurateconvex
  
# bone_x <- threshold
# boneEssai <- bone
#  boneEssai[, , 1, 1] <- !bone_x
#  if (dim(boneEssai)[4]>1) boneEssai[, , 1, 2] <- !bone_x
#  if (dim(boneEssai)[4]>2) boneEssai[, , 1, 3] <- !bone_x
# getFromNamespace("plot.cimg", ns="imager")(boneEssai, bty="n", axes=FALSE, xlab="", ylab="", asp = 1)
  
  method <- tolower(method)
  method <- match.arg(method, choices = c("accurate", "fast", "fastconvex", "accurateconvex"))
  
  # analysis=1; threshold <- NULL; partial=FALSE; center.x=NA; center.y=NA
  if (is.null(threshold)) threshold <- RM_get(x=bone, RMname = analysis, valuename="threshold")
  
  if (is.null(threshold)) {
    threshold <- getFromNamespace(".BP_threshold", ns="BoneProfileR")(bone, analysis=analysis)
  }
  
  # Dans bone et threshold, le x est en dimension 1 et le y en dimension 2
  if (is.na(center.x) | is.na(center.y)) {
    center.y <- mean(which(threshold, arr.ind = TRUE)[, 1])
    center.x <- mean(which(threshold, arr.ind = TRUE)[, 2])
  } else {
    c. <- center.x
    center.x <- center.y
    center.y <- c.
    # points(center.y, center.x, col="red")
  }
  
  contourInitial <- function(bone, analysis=1, center.x, center.y) {
    d_mat <- getFromNamespace(".BP_distance", ns="BoneProfileR")(bone, 
                                                                 threshold=threshold, 
                                                                 analysis=analysis, 
                                                                 center.x=center.x, 
                                                                 center.y=center.y)
    
    d_Threshold_mat <- ifelse(threshold, d_mat, NA)
    d_Threshold <- as.vector(d_Threshold_mat)
    
    a_mat <- getFromNamespace(".BP_angle", ns="BoneProfileR")(bone, 
                                                              threshold=threshold, 
                                                              analysis=analysis, 
                                                              center.x=center.x, 
                                                              center.y=center.y)
    a <- as.vector(a_mat)
    cut.angle <- (dim(threshold)[1]+dim(threshold)[2])*2
    # if (cut.angle > 360*2) cut.angle <- 360*2
    
    la <- seq(from=0, to=2*pi, length.out=cut.angle+1)
    fa <- findInterval(a, la)
    
    dmax_a <- aggregate(d_Threshold, by=list(fa), FUN=function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
    dmax_a2 <- rep(NA, cut.angle)
    dmax_a2[dmax_a[, 1]] <- dmax_a[, 2]
    
    dmax_a2[dmax_a2 < median(dmax_a2, na.rm = TRUE)/4] <- NA
    
    contour <- (as.vector(d_mat) < dmax_a2[fa])
    contour <- matrix(contour, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = FALSE)
    
    contour <- ifelse(is.na(contour), FALSE, contour)
    
    return(contour)
  }
  
  if (grepl("fast", method)) {
    contour <- contourInitial(bone, analysis, center.x, center.y)
    # bone_x <- contour
    # boneEssai <- bone
    #  boneEssai[, , 1, 1] <- !bone_x
    #  if (dim(boneEssai)[4]>1) boneEssai[, , 1, 2] <- !bone_x
    #  if (dim(boneEssai)[4]>2) boneEssai[, , 1, 3] <- !bone_x
    # getFromNamespace("plot.cimg", ns="imager")(boneEssai, bty="n", axes=FALSE, xlab="", ylab="", asp = 1)
    
  }
  
  if (grepl("accurate", method)) {
    contourHL <- contourInitial(bone, analysis, center.x = 1, center.y = 1)
    contourHR <- contourInitial(bone, analysis, center.x = attributes(bone)$dim[1], center.y = 1)
    contourBR <- contourInitial(bone, analysis, center.x = attributes(bone)$dim[1], center.y = attributes(bone)$dim[2])
    contourBL <- contourInitial(bone, analysis, center.x = 1, center.y = attributes(bone)$dim[2])
    contour_x <- contourInitial(bone, analysis, center.x, center.y)
 
    contourML <- contourInitial(bone, analysis, center.x = attributes(bone)$dim[1]/2, center.y = 1)
    contourMR <- contourInitial(bone, analysis, center.x = attributes(bone)$dim[1], center.y = attributes(bone)$dim[2]/2)
    contourBM <- contourInitial(bone, analysis, center.x = attributes(bone)$dim[1]/2, center.y = attributes(bone)$dim[2])
    contourBM <- contourInitial(bone, analysis, center.x = 1, center.y = attributes(bone)$dim[2]/2)

    contour <- (contourHL + contourHR + contourBR + contourBL + contour_x + 
                  contourML + contourMR + contourBM + contourBM  == 9)
    
    # attributes(contour)$dim <- c(attributes(contour)$dim, 1, 1)
    # attributes(contour)$class <- c("cimg", "imager_array", "numeric") 
    # getFromNamespace("plot.cimg", ns="imager")(contour, bty="n", axes=FALSE, xlab="", ylab="", asp = 1)
  }
  
  if (grepl("convex", method)) {
    if (all(rownames(installed.packages()) != "sp")) stop("Install the sp package to use convex option.")
    if (all(rownames(installed.packages()) != "spatstat.geom")) stop("Install the spatstat.geom package to use convex option.")
    dfcontour <- which(contour, arr.ind=TRUE)
    dfcontour <- spatstat.geom::convexhull.xy(x=dfcontour[, "col"], 
                                y=dfcontour[, "row"])
    coord <- expand.grid(1:dim(bone)[2], 1:dim(bone)[1])
    pointinside <- sp::point.in.polygon(coord$Var1,coord$Var2, dfcontour$bdry[[1]]$x, dfcontour$bdry[[1]]$y)
    contour <- matrix(data = as.logical(pointinside), nrow = nrow(contour), ncol = ncol(contour), byrow = TRUE)
    # attributes(contour)$dim <- c(attributes(contour)$dim, 1, 1)
    # attributes(contour)$class <- c("cimg", "imager_array", "numeric") 
    # getFromNamespace("plot.cimg", ns="imager")(contour, bty="n", axes=FALSE, xlab="", ylab="", asp = 1)
  }
  
  return(contour)
  
  
  # # 
  # # library(fields)
  # # image.plot(a_mat)
  # 
  # # cut.angle <- (dim(threshold)[c(1)]*2+dim(threshold)[c(2)])*2
  # 
  # cut.angle <- 360*2
  # 
  # la <- seq(from=0, to=2*pi, length.out=cut.angle+1)
  # fa <- findInterval(a, la)
  # # fa_mat <- matrix(fa, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = FALSE)
  # # fa <- as.vector(fa)
  # 
  # dmax_a <- aggregate(d_Threshold, by=list(fa), FUN=function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
  # dmax_a2 <- rep(NA, cut.angle)
  # dmax_a2[dmax_a[, 1]] <- dmax_a[, 2]
  # 
  # # 23451
  # # 12345
  # # 51234
  # 
  # # Si deux valeurs sont séparées d'un NA, on remplace le NA par la moyenne des deux valeurs
  # # Est-ce vraiment nécessaire ? Pas sur
  # if (FALSE) {
  #   ld <- length(dmax_a2)
  #   kb_dmax_a2 <- c(tail(dmax_a2, n=1), dmax_a2[1:(ld-1)])
  #   ke_dmax_a2 <- c(dmax_a2[-1], dmax_a2[1])
  #   w <- which(is.na(dmax_a2) & !is.na(kb_dmax_a2) & !is.na(ke_dmax_a2))
  #   dmax_a2[w] <- sapply(w, FUN = function(x) mean(kb_dmax_a2[x], ke_dmax_a2[x]))
  # }
  # 
  # # Nouvelle stratégie
  # if (FALSE) {
  #   ld <- length(dmax_a2)
  #   dmax_a2_rep <- c(0, 0, rep(dmax_a2, 3), 0, 0)
  #   r <- which(!is.na(dmax_a2_rep))
  #   
  #   
  #   for (v in seq_along(r[-length(r)])) {
  #     if (r[v+1] != r[v]+1)
  #       dmax_a2_rep[(r[v]+1):(r[v+1]-1)] <- rep(mean(c(dmax_a2_rep[r[v]], dmax_a2_rep[r[v+1]])), (r[v+1]-1)-(r[v]+1)+1)
  #   }
  #   dmax_a2 <- dmax_a2_rep[(ld+3):(ld+3+ld-1)]
  # }
  # 
  # 
  # dmax_a2[dmax_a2 < median(dmax_a2, na.rm = TRUE)/4] <- NA
  # 
  # 
  # contour <- (as.vector(d_mat) < dmax_a2[fa])
  # contour <- matrix(contour, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = FALSE)
  # 
  # contour <- ifelse(is.na(contour), FALSE, contour)
  # 
  # # contour <- matrix(data = FALSE, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = FALSE)
  # # for (x in 1:dim(threshold)[2])
  # #   for (y in 1:dim(threshold)[1]) {
  # #     # dans fa_mat[y, x] j'ai la catégorie et dmax_a2[fa_mat[y, x]] la distance max
  # #     contour[y, x] <- (d_mat[y, x] < dmax_a2[fa_mat[y, x]])
  # #   }
  # # 
  # # contour <- ifelse(is.na(contour), FALSE, contour)
  # 
  # # contour <- t(contour)
  # # 
  # # 
  # # contourTF <- sapply(seq_along(d_Threshold), FUN=function(x) {
  # #   # Dans d j'ai la distance au centre de chaque pixel x
  # #   # Dans fa[x] j'ai la catgéorie d'angle du pixel x
  # #   # Dans dmax_a2 j'ai la distance max pour cet angle
  # #   return(d[x] <= dmax_a2[fa[x]])
  # # })
  # # 
  # # # distance <- matrix(data = d, ncol=dim(threshold)[1], nrow=dim(threshold)[2], byrow = TRUE)
  # # # angle <- matrix(data = fa, ncol=dim(threshold)[1], nrow=dim(threshold)[2], byrow = TRUE)
  # # contour <- matrix(data = contourTF, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = FALSE)
  # # 
  # # 
  # # # d_Threshold <- as.vector(d_Threshold_mat)
  # # # 
  # # # a_mat <- matrix(NA, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = TRUE)
  # # # 
  # # # for (x in 1:dim(threshold)[2]) 
  # # #   for (y in 1:dim(threshold)[1]) {
  # # #     a_mat[y, x] <- (atan((x-center.x)/(y-center.y))+pi*(as.numeric((y-center.y)<0))) %% (2*pi)
  # # #   }
  # # # 
  # # # invisible(sapply(1:nrow(m), FUN = function(r) a_mat[m[r, 2], m[r, 1]] <<- (atan((m[r, 1]-center.y)/(m[r, 2]-center.x))+pi*(as.numeric((m[r, 2]-center.x)<0))) %% (2*pi)))
  # # 
  # # 
  # # 
  # # 
  # # # Droite = 0
  # # # Gauche = -pi
  # # # Haut = -pi/2
  # # # bas = pi/2
  # # # a_mat <- matrix(NA, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = TRUE)
  # # # invisible(sapply(1:nrow(m), FUN = function(r) a_mat[m[r, 2], m[r, 1]] <<- (atan((m[r, 1]-center.y)/(m[r, 2]-center.x))+pi*(as.numeric((m[r, 2]-center.x)<0))) %% (2*pi)))
  # # 
  # # # a <- sapply(1:nrow(m), FUN = function(r) (atan((m[r, 1]-center.y)/(m[r, 2]-center.x))+pi*(as.numeric((m[r, 2]-center.x)<0))) %% (2*pi))
  # # # a_mat <- matrix(a, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = TRUE)
  # # 
  # # cut.angle <- 360
  # # 
  # # la <- seq(from=0, to=2*pi, length.out=cut.angle+1)
  # # fa <- findInterval(a, la)
  # # fa_mat <- matrix(fa, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = FALSE)
  # # 
  # # fa_mat2 <- fa_mat
  # # 
  # # for (r in 1:nrow(m)) {
  # #   va <- a_mat[m[r, 2], m[r, 1]]
  # #   fa_mat2[m[r, 2], m[r, 1]] <- max(which(va>la))
  # # }
  # # 
  # # dmax_a <- aggregate(d_Threshold, by=list(fa), FUN=function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
  # # dmax_a2 <- rep(NA, cut.angle)
  # # dmax_a2[dmax_a[, 1]] <- dmax_a[, 2]
  # # 
  # # fa_mat_2 <- matrix(NA, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = TRUE)
  # # fa_mat_2[] <- dmax_a2[fa_mat]
  # # 
  # # 
  # # contourTF <- matrix(NA, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = TRUE)
  # # invisible(sapply(1:nrow(m), FUN = function(r) {
  # #   # Dans d j'ai la distance au centre de chaque pixel x
  # #   # Dans fa[x] j'ai la catgéorie d'angle du pixel x
  # #   # Dans dmax_a2 j'ai la distance max pour cet angle
  # #   contourTF[m[r, 2], m[r, 1]] <<- dmax_a2[fa[x]]
  # # }))
  # # 
  # # 
  # # 
  # # # Dans fa j'ai la classe de chaque pixel pour l'angle
  # # contourTF <- sapply(seq_along(d_Threshold), FUN=function(x) {
  # #   # Dans d j'ai la distance au centre de chaque pixel x
  # #   # Dans fa[x] j'ai la catgéorie d'angle du pixel x
  # #   # Dans dmax_a2 j'ai la distance max pour cet angle
  # #   return(d[x] <= dmax_a2[fa[x]])
  # # })
  # # 
  # # # distance <- matrix(data = d, ncol=dim(threshold)[1], nrow=dim(threshold)[2], byrow = TRUE)
  # # # angle <- matrix(data = fa, ncol=dim(threshold)[1], nrow=dim(threshold)[2], byrow = TRUE)
  # # contour <- matrix(data = contourTF, nrow=dim(threshold)[1], ncol=dim(threshold)[2], byrow = FALSE)
  # 
  # 
  # # To check if it works
  # 
  # # x <- bone
  # # x[, , 1, 1] <- contour
  # # if (dim(x)[4]>1) x[, , 1, 2] <- !contour
  # # if (dim(x)[4]>2) x[, , 1, 3] <- !contour
  # # par(xaxs="i", yaxs="i")
  # # par(mar=c(4, 0, 0, 0))
  # # 
  # # getFromNamespace("plot.cimg", ns="imager")(x, bty="n", axes=FALSE, xlab="", ylab="", 
  # #                                            asp = 1)
  # 
  # return(contour)
}
