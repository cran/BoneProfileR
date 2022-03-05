#' BP_FitMLCompactness estimates likelihood of model of a bone section
#' @title Estimation of the likelihood of a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The -Ln L
#' @param fitted.parameters Parameters of the model to be fitted
#' @param bone The bone image to be used
#' @param fixed.parameters Fixed parameters of the model
#' @param priors Priors used for intermediate estimations
#' @param replicates.CI Number of replicates to estimate confidence interval
#' @param analysis Name or rank of analysis
#' @param twosteps Does a 2-steps analysis be performed?
#' @param silent Should information be shown?
#' @description Estimation of the model of compactness of a bone section.\cr
#' The two-steps analysis performs first a quasi-Newton method, then a Bayesian MCMC and finally again a quasi-Newton method. 
#' It generally ensures that global minimum is found. On the other hand, it doubles the time to complete.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#'  bone <- BP_OpenImage()
#'  # or, to use the package imager to open a tiff image
#'  bone <- BP_OpenImage(ijtiff=TRUE)
#' library(BoneProfileR)
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_DetectBackground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectForeground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic")
#'  plot(bone, type="mineralized", show.grid=FALSE)
#'  plot(bone, type="unmineralized", show.grid=FALSE)
#'  plot(bone, type="section", show.grid=FALSE)
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic", twosteps=TRUE)
#'  BP_GetFittedParameters(bone)
#'  plot(bone)
#'  plot(bone, type="observations")
#'  plot(bone, type="observations+model", analysis=1)
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  BP_ListAnalyses(bone)
#'  bone <- BP_FitMLCompactness(bone, 
#'                 fitted.parameters=c(fittedpar, K1=1, K2=1), 
#'                 fixed.parameters=NULL, analysis="flexit", twosteps=TRUE)
#'  compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", alloptim=TRUE), 
#'              Flexit=BP_GetFittedParameters(bone, analysis="flexit", alloptim=TRUE))
#'  out4p <- plot(bone, type="observations+model", analysis="logistic")
#'  out6p <- plot(bone, type="observations+model", analysis="flexit")
#' }
#' @export


BP_FitMLCompactness <- function(bone, fitted.parameters=c(P=0.5, S=0.05, Min=0.001, Max=0.999),
                                priors=NULL, 
                                fixed.parameters=c(K1=1, K2=1), twosteps=TRUE, 
                                replicates.CI=10000, analysis=1, silent=FALSE) {
  
  # fitted.parameters=c(P=0.5, S=0.05, Min=-2, Max=5); fixed.parameters=c(K1=1, K2=1); analysis=1
  
  # BP_LnLCompactness(par=fitted.parameters, bone, data_m=NULL, data_nm=NULL, distance.center=NULL, fixed.parameters=fixed.parameters, analysis=analysis)
  
  
  lower_limit <- c(P=0, S=-2, Min=0, Max=0.2, K1=-1000, k2=-1000)
  upper_limit <- c(P=1, S=2, Min=0.8, Max=1, K1=1000, k2=1000)
  lower <- lower_limit[names(fitted.parameters)]
  upper <- upper_limit[names(fitted.parameters)]
  
  o <- optim(par=fitted.parameters, fn=BP_LnLCompactness, 
             fixed.parameters=fixed.parameters, bone=bone, 
             upper = upper, lower = lower, 
             method = "L-BFGS-B", hessian = TRUE, analysis=analysis)
  
  if (twosteps) {
    # lancement en mcmc
    p <- o$par
    if (is.null(priors)) {
      priors <- data.frame(Density=character(), 
                           Prior1=numeric(), 
                           Prior2=numeric(), 
                           SDProp=numeric(), 
                           Min=numeric(), 
                           Max=numeric(), 
                           Init=numeric(), stringsAsFactors = FALSE)
      
      if (!is.na(p["P"])) {
        priors <- rbind(priors, data.frame(Density="dunif", 
                                           Prior1=0, 
                                           Prior2=1, 
                                           SDProp=0.005, 
                                           Min=0, 
                                           Max=1, 
                                           Init=unname(p["P"]), stringsAsFactors = FALSE, 
                                           row.names = "P"))
      }
      if (!is.na(p["S"])) {
        priors <- rbind(priors, data.frame(Density="dunif", 
                                           Prior1=0, 
                                           Prior2=max(c(unname(p["S"])*2, +10)), 
                                           SDProp=0.3, 
                                           Min=0, 
                                           Max=max(c(unname(p["S"])*2, +10)), 
                                           Init=unname(p["S"]), stringsAsFactors = FALSE, 
                                           row.names = "S"))
      }
      if (!is.na(p["Min"])) {
        priors <- rbind(priors, data.frame(Density="dunif", 
                                           Prior1=0, 
                                           Prior2=0.8, 
                                           SDProp=0.2, 
                                           Min=0, 
                                           Max=0.8, 
                                           Init=unname(p["Min"]), stringsAsFactors = FALSE, 
                                           row.names = "Min"))
      }
      if (!is.na(p["Max"])) {
        priors <- rbind(priors, data.frame(Density="dunif", 
                                           Prior1=0.2, 
                                           Prior2=1, 
                                           SDProp=0.2, 
                                           Min=0.2, 
                                           Max=1, 
                                           Init=unname(p["Max"]), stringsAsFactors = FALSE, 
                                           row.names = "Max"))
      }
      if (!is.na(p["K1"])) {
        priors <- rbind(priors, data.frame(Density="dunif", 
                                           Prior1=min(c(unname(p["K1"])*2, -10)), 
                                           Prior2=max(c(unname(p["K1"])*2, +10)), 
                                           SDProp=0.2, 
                                           Min=min(c(unname(p["K1"])*2, -10)), 
                                           Max=max(c(unname(p["K1"])*2, +10)), 
                                           Init=unname(p["K1"]), stringsAsFactors = FALSE, 
                                           row.names = "K1"))
      }
      if (!is.na(p["K2"])) {
        priors <- rbind(priors, data.frame(Density="dunif", 
                                           Prior1=min(c(unname(p["K2"])*2, -10)), 
                                           Prior2=max(c(unname(p["K2"])*2, +10)), 
                                           SDProp=0.2, 
                                           Min=min(c(unname(p["K2"])*2, -10)), 
                                           Max=max(c(unname(p["K2"])*2, +10)), 
                                           Init=unname(p["K2"]), stringsAsFactors = FALSE, 
                                           row.names = "K2"))
      }
    }
    mcmc <- HelpersMG::MHalgoGen(
      likelihood = BP_LnLCompactness,
      bone=bone, 
      fixed.parameters=fixed.parameters, 
      parameters_name = "par", 
      parameters = priors, n.iter = 10000,
      n.chains = 1,
      n.adapt = 100,
      thin = 1, adaptive = TRUE)
    
    fitted.parameters <- HelpersMG::as.parameters(mcmc)
    
    o <- optim(par=fitted.parameters, fn=BP_LnLCompactness, 
               fixed.parameters=fixed.parameters, bone=bone, 
               upper = upper, lower = lower, 
               method = "L-BFGS-B", hessian = TRUE, analysis=analysis)
  }
  
  o$fixed.parameters <- fixed.parameters
  o$AIC <- 2*o$value+2*length(o$par)
  o$SE <- SEfromHessian(o$hessian)
  
  rd <- RandomFromHessianOrMCMC(method = "hessian", 
                                Hessian = o$hessian, 
                                fitted.parameters = o$par, 
                                fixed.parameters = o$fixed.parameters,
                                replicates = replicates.CI, 
                                probs = NULL, silent = silent)
  
  Min <- rd$random[, "Min"]
  Max <- rd$random[, "Max"] 
  
  rd$random[, "Min"] <- Min
  rd$random[, "Max"] <- Max
  
  m <- matrix(data =c(apply(rd$random, MARGIN=2, FUN=mean), apply(rd$random, MARGIN=2, FUN=sd)), 
              ncol=2)
  
  colnames(m) <- c("Mean", "SE")
  rownames(m) <- colnames(rd$random)
  
  o$summary.table <- m
  
  # Calcul de l'intervalle de confiance avec Hessian
  
  data <- RM_get(x=bone, RMname=analysis, valuename = "compactness.synthesis")
  outHessian <- matrix(NA, ncol=nrow(data), nrow=replicates.CI)
  
  for (iter in 1:replicates.CI) {
    p <- unlist(rd$random[iter, , drop=TRUE])
    
    # 21/02/2020
    p["S"] <- 1/(4*p["S"])
    
    c <- flexit(x = data$distance.center, 
                par = p) * (p["Max"] - p["Min"]) + p["Min"]
    
    outHessian[iter, ] <- c
  }
  
  qHessian <- apply(X = outHessian, MARGIN = 2, FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
  colnames(qHessian) <- data$distance.center
  
  o$quantiles <- qHessian
  
  bone <- RM_add(x=bone, RMname = analysis, valuename = "optim", value=o)
  
  return(bone)
}


