#' BP_EstimateCompactness estimates the compactness of a bone section
#' @title Estimation of the compactness of a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The orignial bone object with a new attribute for compactness
#' @param bone The bone image to be used
#' @param center Which center to be used: user, mineralized, unmineralized, section, ontogenetic
#' @param cut.angle Number of angles
#' @param cut.distance Number of distances
#' @param NbRemoveDistanceExterior How many exterior sectors should be removed from analysis?
#' @param partial Is the section partial?
#' @param NbRemoveEdgePartial How many radial section to remove at the edge of partial section?
#' @param rotation.angle The angle of rotation for analysis
#' @param analysis The name or rank of analysis
#' @param method Can be Fast, Accurate, FastConvex, or AccurateConvex
#' @param show.plot should plot is shown ?
#' @param cut.max The number of slices for the internal estimation
#' @description Estimation of the compactness of a bone section.\cr
#' The reference for radial estimation of compactness is the trigonometric circle for rotation.angle=0 in 
#' BP_EstimateCompactness():\cr
#' - The top of the section is located at -pi/2.\cr
#' - The left of the section is located at -pi and +pi.\cr
#' - The bottom of the section is located at pi/2.\cr
#' - The right of the section is 0.\cr
#' If rotation.angle is different from 0, the value of rotation.angle is added to the angle modulo 2.pi.\cr
#' The method Fast works well with the convex bone section while if the section is concave, Accurate is slower but works well in all circonstances.\cr
#' Fast method is maintained here only for compatibility with versions <3.1 of BoneProfileR.\cr
#' If the section is concave, the methods FastConvex and AccurateConvex return a minimum convex section.\cr
#' If the center has been automatically detected, the method parameter is ignored because it has already
#' been used with the function BP_DetectCenters().
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
#'  bone <- BP_EstimateCompactness(bone)
#'  plot(bone, type="original", show.grid=FALSE)
#'  plot(bone, type="mineralized", show.grid=FALSE)
#'  plot(bone, type="unmineralized", show.grid=FALSE)
#'  plot(bone, type="section", show.grid=FALSE)
#'  plot(bone, type="original", show.grid=TRUE)
#' }
#' @export


BP_EstimateCompactness <- function(bone                          , 
                                   center="ontogenetic"          , 
                                   partial=FALSE                 , 
                                   NbRemoveEdgePartial=1         ,
                                   cut.angle=60                  , 
                                   cut.distance=100              , 
                                   NbRemoveDistanceExterior = 1  ,
                                   rotation.angle=0              , 
                                   analysis=1                    , 
                                   method="Fast"                 ,
                                   show.plot=TRUE                , 
                                   cut.max = 360                 ) {
  
  # center="ontogenetic"; partial=FALSE; cut.angle=60; NbRemoveEdgePartial=1; cut.distance=100; NbRemoveDistanceExterior = 1; analysis=1; rotation.angle=0; show.plot=TRUE; method="Fast"; cut.max = 360
  # center="user"; partial=TRUE; cut.angle=60; cut.distance=100; analysis=1; rotation.angle=0; show.plot=TRUE
  
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))                 # code line i + 1
  center <- tolower(center)
  center <- match.arg(center, choices = c("user", "mineralized", "ontogenic", 
                                          "unmineralized", "section", "ontogenetic"))
  if (center == "ontogenic") center <- "ontogenetic"
  if ((center != "user") & partial) {
    stop("When analysis of partial bone section is performed, only user center must be used.")
  }
  
  if (is.null(RM_get(x=bone, RMname=analysis, valuename = "centers"))) {
    stop("You must first setup centers using BP_DetectCenters() or BP_ChooseCenter()")
  }
  
  
  # Je formate la coupe en threshold
  threshold <- RM_get(x=bone, RMname = analysis, valuename="threshold")
  if (is.null(threshold)) {
    threshold <- getFromNamespace(".BP_threshold", ns="BoneProfileR")(bone, analysis=analysis)
    bone <- RM_add(x=bone, RMname = analysis, valuename="threshold", 
                   value=threshold)
  }
  
  if ((isTRUE(partial)) & (method=="Accurate")) {
    warning("Accurate method is not compatible with partial section; fast method will be used rather.")
    method <- "fast"
  }
  
  if ((isTRUE(partial)) & (method=="AccurateConvex")) {
    warning("AccurateConvex method is not compatible with partial section; FastConvex method will be used rather.")
    method <- "fastconvex"
  }
  
  # bone_gs <- grayscale(bone)
  
  if (center == "user") {
    center.x <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_user.x"]
    center.y <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_user.y"]
  }
  if (center == "mineralized") {
    center.x <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_cortex.x"]
    center.y <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_cortex.y"]
  }
  if (center == "unmineralized") {
    center.x <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_medula.x"]
    center.y <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_medula.y"]
  }
  if (center == "section") {
    center.x <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_bone.x"]
    center.y <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_bone.y"]
  }
  
  if (center == "ontogenetic") {
    center.x <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_ontogenic.x"]
    center.y <- RM_get(x=bone, RMname=analysis, valuename = "centers")["GC_ontogenic.y"]
  }
  
  if (is.na(center.x) | is.na(center.y)) stop("The requested center is not available")
  
  # Le contour calculé ici est pour un centre de la section minéralisé
  contour <- RM_get(x=bone, RMname=analysis, valuename = "contour")
  if (is.null(contour))  {
    contour <- getFromNamespace(".BP_contour", ns="BoneProfileR")(bone, 
                                                                  analysis=analysis, 
                                                                  threshold=threshold, 
                                                                  partial=partial, 
                                                                  center.x=center.x, 
                                                                  center.y=center.y, 
                                                                  method=method)
    
    bone <- RM_add(x=bone, RMname = analysis, valuename="contour", 
                   value=contour)
    bone <- RM_add(x=bone, RMname = analysis, valuename="method", value=method)
    
    if (FALSE) {
      plot(bone)
      vct <- cbind(expand.grid(1:dim(contour)[1], 1:dim(contour)[2]), point=as.vector(contour))
      vct <- subset(vct, subset=point)
      points(x=vct$Var1, y=vct$Var2, pch=".")
    }
  } else {
    if (RM_get(x=bone, RMname = analysis, valuename="method") != method) {
      contour <- getFromNamespace(".BP_contour", ns="BoneProfileR")(bone, 
                                                                    analysis=analysis, 
                                                                    threshold=threshold, 
                                                                    partial=partial, 
                                                                    center.x=center.x, 
                                                                    center.y=center.y, 
                                                                    method=method)
      
      bone <- RM_add(x=bone, RMname = analysis, valuename="contour", 
                     value=contour)
      bone <- RM_add(x=bone, RMname = analysis, valuename="method", value=method)
    }
  }
  
  
  # sum(contour) est le nombre de pixels de la section
  
  compactness <- expand.grid(1:dim(bone)[2], 1:dim(bone)[1])
  colnames(compactness) <- c("y", "x")
  m1 <- compactness[, "x"] - center.x
  m2 <- compactness[, "y"] - center.y
  
  compactness <- cbind(compactness, distance.center=sqrt(m1^2 + m2^2), 
                       distance.external=NA, 
                       angle=((atan2(m2, m1) + pi + rotation.angle) %% (2*pi))-pi,
                       mineral=as.numeric(t(threshold)), 
                       # gradient=as.numeric(t(bone_gs[, , 1, 1])), 
                       gradient=as.numeric(t(threshold)), 
                       contour=as.vector(t(contour)))
  
  l360 <- seq(from=-pi, to=pi, length.out = cut.max+1)
  compactness <- cbind(compactness, 
                       cut.360=cut(compactness$angle, breaks = l360))
  
  langle <- seq(from=-pi, to=pi, length.out = cut.angle+1)
  compactness <- cbind(compactness, 
                       cut.angle=cut(compactness$angle, breaks = langle))
  
  
  
  # J'ai les informations pour chaque pixel de la coupe
  # getFromNamespace("plot.cimg", ns="imager")(bone, bty="n", axes=FALSE, xlab="", ylab="", asp = 1)
  # Chaque segment est passé
  levels.360 <- cut(l360, breaks = l360)[-1]
  levels.angle <- cut(langle, breaks = langle)[-1]
  nlevelparangle <- length(levels.360)/length(levels.angle)
  levels.360.by.angle <- cut(l360+1E-6, breaks = langle)[-length(l360)]
  angles <- langle[-length(langle)]+(langle[2]-langle[1])/2
  names(angles) <- levels.angle
  
  # Maintenant je ne garde que ceux de contour
  # C'est là où je ne dois garder 
  compactness <- compactness[compactness$contour, c("x", "y", "distance.center", 
                                                    "distance.external", "angle", 
                                                    "mineral", "gradient", "cut.360", "cut.angle")]
  
  
  if (FALSE) {
    plot(bone)
    points(x=compactness$x, y=compactness$y, pch=".")
  }
  
  
   
  # Cut angle est à 360 sections donc 360 bords
  # Je calcule la distance max sur chacune des portions à 360°
  peripherie <- NULL
  for (l in levels.360) {
    dc <- subset(compactness, subset = (compactness$cut.360 == l), select="distance.center")
    if (nrow(dc) != 0) {
      m <- max(dc[, "distance.center"])
      s <- which(compactness$cut.360 == l)
      compactness[s, "distance.external"] <- rep(m, times=length(s))
      peripherie <- c(peripherie, m)
    } else {
      peripherie <- c(peripherie, NA)
    }
  }
  
  peripherie <- data.frame(peripherie=peripherie, levels.360=levels.360, l360.by.angle=levels.360.by.angle)
  
  if (partial) {
    # Je détermine les angles où j'ai de l'information
    borne <- seq(from=-pi, to=pi, length.out = cut.max+1)
    borne.begin <- borne[-length(borne)]
    borne.last <- borne[-1]
    # Je les duplique 3 fois et je vais chercher le centre
    peripheriex <- data.frame(angle.begin=c(borne.begin-2*pi, 
                                      borne.begin, 
                                      borne.begin+2*pi), 
                              angle.center=(c(borne.last-2*pi, 
                                              borne.last, 
                                              borne.last+2*pi)+c(borne.begin-2*pi, 
                                                                 borne.begin, 
                                                                 borne.begin+2*pi))/2, 
                              angle.last=c(borne.last-2*pi, 
                                           borne.last, 
                                           borne.last+2*pi), 
                              peripherie=rbind(peripherie, peripherie, peripherie))
    peripheriex <- cbind(peripheriex, peripherie.begin=NA)
    peripheriex <- cbind(peripheriex, peripherie.last=NA)
    
    # C'est quoi ça ??
    
    peripheriex[, "peripherie.begin"] <- (peripheriex[c(nrow(peripheriex), 1:(nrow(peripheriex)-1)), "peripherie.peripherie"] + 
                                            peripheriex[, "peripherie.peripherie"])/2
    peripheriex[, "peripherie.last"] <- (peripheriex[, "peripherie.peripherie"] + 
                                           peripheriex[c(2:(nrow(peripheriex)), 1), "peripherie.peripherie"])/2
    mper <- max(peripheriex$peripherie.peripherie, na.rm = TRUE)
    pos <- (peripheriex$peripherie.peripherie/mper) > 0.75
    pos <- ifelse(is.na(pos), FALSE, pos)
    # plot(1:nrow(peripheriex), peripherie$peripherie.peripherie, col=ifelse(pos, "red", "grey"))
    # posT <- which(pos)
    
    # De posT à postT2 j'ai TRUE
    posT <- which(pos[-1] & !pos[-length(pos)])+1 # FALSE TRUE
    posT2 <- which(pos[-length(pos)] & !pos[-1]) # TRUE FALSE
    
    
    deb <- posT[1]
    fin <- posT2[posT2 > deb][1]
    
    postab <- (((deb-1):(fin-1)) %% cut.max)+1
    postab <- postab[(NbRemoveEdgePartial+1):(length(postab)-NbRemoveEdgePartial)]
    
    peripherie <- peripheriex[postab, ]
    peripherie[, "angle.begin"] <- (peripherie[, "angle.begin"] %% (2*pi))
    peripherie[, "angle.last"] <- (peripherie[, "angle.last"] %% (2*pi))
    peripherie[, "angle.center"] <- (peripherie[, "angle.center"] %% (2*pi))
    
    
    # Je ne dois garder que ceux qui sont entiers c'est à dire que il y a des informations pour tous les niveaux
    
    nlevels <- round(cut.max/cut.angle, 0)
    tbl <- table(peripherie$peripherie.l360.by.angle)
    peripherie <- subset(peripherie, subset=peripherie$peripherie.l360.by.angle %in% (names(tbl)[tbl == nlevels]))
    peripherie <- na.omit(peripherie)
    
  } else {
    peripherie <- data.frame(angle=levels.360, 
                             angle.begin=l360[-length(l360)] %% (2*pi), 
                             angle.last=l360[-1] %% (2*pi), 
                             angle.center=(l360[-length(l360)]+(l360[-1]-l360[-length(l360)])/2)%% (2*pi),  
                             peripherie=peripherie)
    
    peripherie <- cbind(peripherie, peripherie.begin=NA)
    peripherie <- cbind(peripherie, peripherie.last=NA)
    
    peripherie[, "peripherie.begin"] <- (peripherie[c(nrow(peripherie), 1:(nrow(peripherie)-1)), "peripherie.peripherie"] + 
                                            peripherie[, "peripherie.peripherie"])/2
    peripherie[, "peripherie.last"] <- (peripherie[, "peripherie.peripherie"] + 
                                           peripherie[c(2:(nrow(peripherie)), 1), "peripherie.peripherie"])/2
  }
  
  peripherie <- cbind(peripherie, include.begin = FALSE)
  peripherie <- cbind(peripherie, include.last = FALSE)
  
  
  bc <- seq(from=-pi, to=pi, length.out = cut.angle+1) %% (2*pi)
  for (p in 1:nrow(peripherie)) {
    if (length(which(abs(bc - peripherie[p, "angle.begin"]) < 1E-6)) != 0) {
      peripherie[p, "include.begin"] <- TRUE
    }
    if (length(which(abs(bc - peripherie[p, "angle.last"]) < 1E-6)) != 0) {
      peripherie[p, "include.last"] <- TRUE
    }
  }
  peripherie <- peripherie[min(which(peripherie[, "include.begin"])): max(which(peripherie[, "include.last"])), ]
  peripherie <- cbind(peripherie, angles=unname(angles[match(peripherie$peripherie.l360.by.angle, names(angles))]))
  peripherie$peripherie.l360.by.angle <- droplevels(peripherie$peripherie.l360.by.angle)
  
  bone <- RM_add(x=bone, RMname = analysis, valuename="peripherie", value=peripherie)
  
  
  
  # Je dois aussi tronquer les autres
  compactness$cut.angle <- droplevels(compactness$cut.angle)
  compactness <- cbind(compactness, ratio.center=compactness$distance.center/compactness$distance.external)
  compactness <- cbind(compactness, 
                       cut.distance.center=cut(compactness$ratio.center, breaks = seq(from=0, to=1, length.out = cut.distance+1)))
  
  if (NbRemoveDistanceExterior != 0) {
    keepDistances <- as.character(rev(cut(seq(from=0, to=1, length.out = cut.distance+1), breaks = seq(from=0, to=1, length.out = cut.distance+1))[-1])[(NbRemoveDistanceExterior + 1):cut.distance])
    compactness <- compactness[(as.character(compactness$cut.distance.center) %in% keepDistances), ]
    compactness$cut.distance.center <- droplevels(compactness$cut.distance.center)
  }
  
  lac <- levels(peripherie$peripherie.l360.by.angle)
  lacx <- gsub("^\\(", "", lac)
  lacx <- order(as.numeric(gsub(",.+", "", lacx)))
  lac <- lac[lacx]
  
  compactness$cut.angle <- droplevels(compactness$cut.angle)
  compactness$cut.distance.center <- droplevels(compactness$cut.distance.center)
  
  t <- table(compactness$cut.angle, compactness$cut.distance.center, compactness$mineral)
  
  # Je dois exclure les angles pas bons
  
  ldc <- levels(compactness$cut.distance.center)
  ldcx <- gsub("^\\(", "", ldc)
  ldcx <- order(as.numeric(gsub(",.+", "", ldcx)))
  ldc <- ldc[ldcx]
  
  ldc1 <- gsub("^\\(", "", ldc)
  ldc1 <- gsub("]$", "", ldc1)
  ldc2 <- as.numeric(gsub(".+,", "", ldc1))
  ldc1 <- as.numeric(gsub(",.+$", "", ldc1))
  ldc1 <- (ldc1 + ldc2)/2
  
  m <- NULL
  nm <- NULL
  for (l in ldc) {
    m <- c(m, sum(t[, l, "1"]))
    nm <- c(nm, sum(t[, l, "0"]))
  }
  
  gr <- aggregate(compactness$gradient, by=list(compactness$cut.angle, compactness$cut.distance.center), 
                  FUN=mean)
  gr <- with(gr, {
    out <- matrix(nrow=nlevels(Group.1), ncol=nlevels(Group.2),
                  dimnames=list(levels(Group.1), levels(Group.2)))
    out[cbind(Group.1, Group.2)] <- x
    out
  })
  
  gr <- gr[rownames(t), colnames(t)]
  
  compactness.synthesis <- data.frame(distance.center=ldc1, 
                                      mineralized=m, 
                                      unmineralize=nm,
                                      compactness=m/(m+nm))
  

  
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optim")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimRadial")
  
  bone <- RM_add(x=bone, RMname = analysis, valuename="compactness", value=compactness)
  bone <- RM_add(x=bone, RMname = analysis, valuename="array.compactness", value=t)
  bone <- RM_add(x=bone, RMname = analysis, valuename="array.gradient", value=gr)
  bone <- RM_add(x=bone, RMname = analysis, valuename="cut.distance.center", value=seq(from=0, to=1, length.out = cut.distance+1)[1:(cut.distance+1-NbRemoveDistanceExterior)])
  bone <- RM_add(x=bone, RMname = analysis, valuename="cut.angle", value=seq(from=-pi, to=pi, length.out = cut.angle+1))
  bone <- RM_add(x=bone, RMname = analysis, valuename="used.centers", value=c(center.x=unname(center.x), center.y=unname(center.y)))
  bone <- RM_add(x=bone, RMname = analysis, valuename="compactness.synthesis", value=compactness.synthesis)
  bone <- RM_add(x=bone, RMname = analysis, valuename="partial", value=partial)
  bone <- RM_add(x=bone, RMname = analysis, valuename="rotation.angle", value=rotation.angle)
  bone <- RM_add(x=bone, RMname = analysis, valuename="global.compactness", value=sum(compactness[, "mineral"])/(nrow(compactness)))
  
  if (show.plot) {
    #  bonex <<- bone
    plot(bone, type="observations")
    
    # par(xaxs="i", yaxs="r")
    # par(mar=c(4, 4, 2, 4)+0.4)
    # plot(compactness.synthesis$distance.center, m/(m+nm), xlim=c(0, 1), ylim=c(0, 1),
    #      type="l", las=1, bty="n", xlab="Distance from the center", ylab="Compactness", lwd=2)
    # lines(x=compactness.synthesis$distance.center, y=(m+nm)/max(m+nm), col="blue")
    # axis(side = 4, at=seq(from=0, to=1, by=0.2), labels = round(seq(from=0, to=1, by=0.2)*max(m+nm), 0), las=1, 
    #      col.axis="blue", col="blue")
    # mtext("Number of pixels", side=4, line=3, col="blue")
  }
  return(bone)
}


