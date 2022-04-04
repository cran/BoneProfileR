#' BP_EstimateCompactness estimates the compactness of a bone section
#' @title Estimation of the compactness of a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The orignial bone object with a new attribute for compactness
#' @param bone The bone image to be used
#' @param center Which center to be used: user, mineralized, unmineralized, section, ontogenetic
#' @param cut.angle Number of angles
#' @param cut.distance Number of distances
#' @param partial Is the section partial?
#' @param rotation.angle The angle of rotation for analysis
#' @param analysis The name or rank of analysis
#' @param show.plot should plot is shown ?
#' @description Estimation of the compactness of a bone section.\cr
#' The reference for radial estimation of compactness is the trigonometric circle for rotation.angle=0 in 
#' BP_EstimateCompactness():\cr
#' - The top of the section is located at -pi/2.\cr
#' - The left of the section is located at -pi and +pi.\cr
#' - The bottom of the section is located at pi/2.\cr
#' - The right of the section is 0.\cr
#' If rotation.angle is different from 0, the value of rotation.angle is added to the angle modulo 2.pi.
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
#' }
#' @export


BP_EstimateCompactness <- function(bone, center="ontogenetic", 
                                   partial=FALSE, 
                                   cut.angle=60, cut.distance=100, 
                                   rotation.angle=0, 
                                   analysis=1, show.plot=TRUE) {
  
  # center="ontogenetic"; partial=FALSE; cut.angle=60; cut.distance=100; analysis=1; rotation.angle=0; show.plot=TRUE
  # center="user"; partial=TRUE; cut.angle=60; cut.distance=100; analysis=1; rotation.angle=0; show.plot=TRUE
  
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  
  center <- match.arg(center, choices = c("user", "mineralized", "ontogenic", 
                                          "unmineralized", "section", "ontogenetic"))
  if (center == "ontogenic") center <- "ontogenetic"
  if ((center!="user") & partial) {
    stop("When partial analysis is done, only user center must be used.")
  }
  
  if (is.null(RM_get(x=bone, RMname=analysis, valuename = "centers"))) {
    stop("You must first setup centers using BP_DetectCenters() or BP_ChooseCenter()")
  }
  
  
  # Je formatte la coupe en threshold
  threshold <- RM_get(x=bone, RMname = analysis, valuename="threshold")
  if (is.null(threshold)) {
    threshold <- getFromNamespace(".BP_threshold", ns="BoneProfileR")(bone, analysis=analysis)
    bone <- RM_add(x=bone, RMname = analysis, valuename="threshold", 
                   value=threshold)
  }
  
  
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
  
 
  contour <- RM_get(x=bone, RMname=analysis, valuename = "contour")
  if (is.null(contour)) {
    contour <- getFromNamespace(".BP_contour", ns="BoneProfileR")(bone, analysis=analysis, 
                                                                  threshold=threshold, 
                                                                  partial=partial, 
                                                                  center.x=center.x, 
                                                                  center.y=center.y)
    
    bone <- RM_add(x=bone, RMname = analysis, valuename="contour", 
                   value=contour)
  }

  
  # sum(contour) est le nombre de pixels mineralisés
  
  compactness <- data.frame(x=rep(NA, times=sum(contour)), 
                            y=rep(NA, times=sum(contour)), 
                            distance.center=rep(NA, times=sum(contour)), 
                            distance.external=rep(NA, times=sum(contour)), 
                            angle=rep(NA, times=sum(contour)),
                            mineral=rep(NA, times=sum(contour)))
  
  
  # if (pg) pb <- progress_bar$new(total = ncol(contour))
  cpt <- 1
  # contour a les x en ligne et les y en colonne
  ytot <- 1:ncol(contour)
  
  for (x in 1:nrow(contour)) {
    
    lgn <- contour[x, ]
    if (any(lgn)) {
      pixel <- as.numeric(threshold[x, lgn])
      y <- ytot[lgn]
      dc <- sqrt((y-center.y)^2+(x-center.x)^2)
      angle <- atan2(y-center.y, x-center.x)
      angle <- ((angle+pi+rotation.angle) %% (2*pi))-pi
      
      
      compactness[cpt:(cpt+length(y)-1), c("x", "y", "distance.center", 
                                           "distance.external", "angle", "mineral")] <- c(rep(x, times=length(y)), y, dc, 
                                                                                          rep(NA, times=length(y)), angle, pixel)
      cpt <- cpt + length(y)
    }
    
    # if (pg) pb$tick()
  }
  
  compactness <- cbind(compactness, 
                       cut.360=cut(compactness$angle, breaks = seq(from=-pi, to=pi, length.out = 360+1)))
  compactness <- cbind(compactness, 
                       cut.angle=cut(compactness$angle, breaks = seq(from=-pi, to=pi, length.out = cut.angle+1)))
  
  # Cut angle est à 360 sections
  peripherie <- NULL
  for (l in levels(compactness$cut.360)) {
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
  peripherie <- data.frame(angle=seq(from=-pi, to=pi, length.out = 360), 
                           peripherie=peripherie)
  
  if (partial) {
    peripherie <- peripherie[c(2:nrow(peripherie)-1), ]
  }
  
  
  peripherie <- na.omit(peripherie)
  
  bone <- RM_add(x=bone, RMname = analysis, valuename="peripherie", value=peripherie)
  
  # Je dois aussi tronquer les autres
  
  
  
  compactness <- cbind(compactness, ratio.center=compactness$distance.center/compactness$distance.external)
  compactness <- cbind(compactness, 
                       cut.distance.center=cut(compactness$ratio.center, breaks = seq(from=0, to=1, length.out = cut.distance+1)))
  
  t <- table(compactness$cut.angle, compactness$cut.distance.center, compactness$mineral)
  
  m <- NULL
  nm <- NULL
  
  for (l in levels(compactness$cut.distance.center)) {
    m <- c(m, sum(t[, l, "1"]))
    nm <- c(nm, sum(t[, l, "0"]))
  }
  
  compactness.synthesis <- data.frame(distance.center=(seq(from=0, to=1, length.out = cut.distance+1)[-1]+
                                                         rev(rev(seq(from=0, to=1, length.out = cut.distance+1))[-1]))/2, 
                                      mineralized=m, 
                                      unmineralize=nm,
                                      compactness=m/(m+nm))
  
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optim")
  bone <- RM_delete(x=bone, RMname = analysis, valuename="optimRadial")
  
  bone <- RM_add(x=bone, RMname = analysis, valuename="compactness", value=compactness)
  bone <- RM_add(x=bone, RMname = analysis, valuename="array.compactness", value=t)
  bone <- RM_add(x=bone, RMname = analysis, valuename="cut.distance.center", value=seq(from=0, to=1, length.out = cut.distance+1))
  bone <- RM_add(x=bone, RMname = analysis, valuename="cut.angle", value=seq(from=-pi, to=pi, length.out = cut.angle+1))
  bone <- RM_add(x=bone, RMname = analysis, valuename="used.centers", value=c(center.x=unname(center.x), center.y=unname(center.y)))
  bone <- RM_add(x=bone, RMname = analysis, valuename="compactness.synthesis", value=compactness.synthesis)
  bone <- RM_add(x=bone, RMname = analysis, valuename="partial", value=partial)
  bone <- RM_add(x=bone, RMname = analysis, valuename="rotation.angle", value=rotation.angle)
  bone <- RM_add(x=bone, RMname = analysis, valuename="global.compactness", value=sum(compactness[, "mineral"])/(nrow(compactness)))
  
  
  
  if (show.plot) {
    #  bonex <<- bone
    #  plot(bone, type="observations")
    
    par(xaxs="i", yaxs="r")
    par(mar=c(4, 4, 2, 4)+0.4)
    plot(compactness.synthesis$distance.center, m/(m+nm), xlim=c(0, 1), ylim=c(0, 1),
         type="l", las=1, bty="n", xlab="Distance from the center", ylab="Compactness", lwd=2)
    lines(x=compactness.synthesis$distance.center, y=(m+nm)/max(m+nm), col="blue")
    axis(side = 4, at=seq(from=0, to=1, by=0.2), labels = round(seq(from=0, to=1, by=0.2)*max(m+nm), 0), las=1, 
         col.axis="blue", col="blue")
    mtext("Number of pixels", side=4, line=3, col="blue")
  }
  return(bone)
}


