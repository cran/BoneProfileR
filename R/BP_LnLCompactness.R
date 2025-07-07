#' BP_LnLCompactness estimates likelihood of model of a bone section
#' @title Estimation of the likelihood of a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The -Ln L
#' @param par Parameters of the model
#' @param bone The bone image to be used
#' @param fixed.parameters Fixed parameters of the model
#' @param data_m Number of mineralized pixels
#' @param data_nm Number of non-mineralized pixels
#' @param distance.center Distances to the center
#' @param analysis Name or rank of analysis
#' @param sign The likelihood if multiplied by sign (-1 or +1) to return -Ln L or Ln L
#' @description Estimation of the compactness of a bone section.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_DetectBackground(bone=bone)
#'  bone <- BP_DetectForeground(bone=bone)
#'  bone <- BP_DetectCenters(bone=bone)
#'  bone <- BP_EstimateCompactness(bone)
#'  plot(bone)
#' }
#' @export


BP_LnLCompactness <- function(par, bone=NULL, 
                              data_m=NULL, data_nm=NULL, 
                              distance.center=NULL,
                              fixed.parameters=NULL, 
                              analysis=1, 
                              sign = -1) {
  
  p <- c(par, fixed.parameters)
  # Min <- p["Min"]
  # Max <- p["Max"]
  
  # print(p)
  
  # 21/02/2020
  # p["S"] <- 1/(4*p["S"])
  
  
  if (inherits(bone, "BoneProfileR")) {
    # p <- c(P=0.5, S=0.1, K1=1, K2=1, Min=0.05, Max=0.99)
    data <- RM_get(x=bone, RMname=analysis, valuename = "compactness.synthesis")
    if (is.null(distance.center)) distance.center <- data$distance.center
    if (is.null(data_m)) data_m <- data$mineralized
    if (is.null(data_nm)) data_nm <- data$unmineralize
  }
  
  # data_m[length(data_m)] <- data_m[length(data_m)] + data_nm[length(data_nm)]
  # data_nm[length(data_nm)] <- 0
  
  c <- BP_flexit(x = distance.center, 
              par = p) # * (Max - Min) + Min
  # print(d(c(c, p)))
  c <- ifelse(c<1E-10, 1E-10, c)
  c <- ifelse(c>1-1E-10, 1-(1E-10), c)
  L <- dbinom(x=data_m, 
              size = data_nm + data_m, 
              prob=c, log = TRUE)
  
  LnL <- sign * sum(L)
  if (is.na(LnL)) LnL <- sign * -1E6 # < log(1E-250)
  # message(paste0("-LnL", as.character(LnL)))
  return(LnL)
}


