#' BP_FitBayesianCompactness estimates Bayesian model of a bone section
#' @title Estimation of Bayesian model of a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The -Ln L
#' @param bone The bone image to be used
#' @param priors Priors 
#' @param n.iter Number of iterations
#' @param n.chains Number of chains
#' @param n.adapt Number of iteration to adapt
#' @param thin Thin parameter for analysis
#' @param analysis Name or rank of analysis
#' @param adaptive Should SDProp be changed during iterations
#' @param silent Should some information must me shown ?
#' @description Estimation of Bayesian model of a bone section./r
#' Get information using ?MHalgoGen.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#'  library(BoneProfileR)
#'  path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_DetectBackground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectForeground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic")
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic")
#'  plot(bone)
#'  plot(bone, type="observations")
#'  plot(bone, type="observations+model", analysis=1)
#'  fittedpar <- BP_GetFittedParameters(bone, analysis="logistic", ML=TRUE, return.all=FALSE)[, "mean"]
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  bone <- BP_FitMLCompactness(bone, 
#'                 fitted.parameters=c(fittedpar, K1=1, K2=1), 
#'                 fixed.parameters=NULL, analysis="flexit")
#'  compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", ML=TRUE, return.all=TRUE), 
#'              Flexit=BP_GetFittedParameters(bone, analysis="flexit", ML=TRUE, return.all=TRUE))
#'  out4p <- plot(bone, type="observations+model", analysis="logistic")
#'  out6p <- plot(bone, type="observations+model", analysis="flexit")
#'  bone <- BP_FitBayesianCompactness(bone, analysis="logistic")
#'  plot(bone, type="observations+model", CI="MCMC", analysis="logistic")
#'  bone <- BP_FitBayesianCompactness(bone, analysis="flexit")
#'  plot(bone, type="observations+model", CI="MCMC", analysis="flexit")
#' }
#' @export


BP_FitBayesianCompactness <- function(bone=stop("A result from BP_FitMLCompactness() must be provided") ,
                                      priors=NULL                                                       , 
                                      n.iter = 10000                                                    ,
                                      n.chains = 1                                                      ,
                                      n.adapt = 5000                                                    ,
                                      thin = 10                                                          , 
                                      analysis=1                                                        , 
                                      adaptive = TRUE                                                   , 
                                      silent=TRUE                                                       ) {
  
  # priors=NULL; n.iter = 10000; n.chains = 1; n.adapt = 100; thin = 1, analysis=1
  
  if (is.null(RM_get(x=bone, RMname=analysis, valuename = "optim"))) stop("The model must be first fitted with BP_FitMLCompactness()")
  
  if (is.null(priors)) {
    priors <- data.frame(Density=character(), 
                         Prior1=numeric(), 
                         Prior2=numeric(), 
                         SDProp=numeric(), 
                         Min=numeric(), 
                         Max=numeric(), 
                         Init=numeric(), stringsAsFactors = FALSE)
    p <- BP_GetFittedParameters(bone, analysis = analysis, ML=TRUE, return.all=FALSE)[, "mean"]

    if (!is.na(p["P"])) {
      priors <- rbind(priors, data.frame(Density="dunif", 
                                         Prior1=0, 
                                         Prior2=1, 
                                         SDProp=0.2, 
                                         Min=0, 
                                         Max=1, 
                                         Init=unname(p["P"]), stringsAsFactors = FALSE, 
                                         row.names = "P"))
    }
    if (!is.na(p["S"])) {
      priors <- rbind(priors, data.frame(Density="dunif", 
                                         Prior1=0, 
                                         Prior2=10, 
                                         SDProp=0.2, 
                                         Min=0, 
                                         Max=10, 
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
                                         Prior1=-10, 
                                         Prior2=10, 
                                         SDProp=0.2, 
                                         Min=-10, 
                                         Max=10, 
                                         Init=unname(p["K1"]), stringsAsFactors = FALSE, 
                                         row.names = "K1"))
    }
    if (!is.na(p["K2"])) {
      priors <- rbind(priors, data.frame(Density="dunif", 
                                         Prior1=-10, 
                                         Prior2=10, 
                                         SDProp=0.2, 
                                         Min=-10, 
                                         Max=10, 
                                         Init=unname(p["K2"]), stringsAsFactors = FALSE, 
                                         row.names = "K2"))
    }
    if (!silent) priors
  }
  
  fixedpar <- BP_GetFittedParameters(bone, analysis = analysis, ML=TRUE, return.all = TRUE)$fixed.parameters
  
  mcmc <- MHalgoGen(
    likelihood = BP_LnLCompactness,
    bone=bone, 
    fixed.parameters=fixedpar, 
    parameters_name = "par", 
    parameters = priors, n.iter = n.iter,
    n.chains = n.chains,
    n.adapt = n.adapt,
    thin = thin, adaptive = adaptive)
  
  data <- RM_get(x=bone, RMname=analysis, valuename = "compactness.synthesis")
  n.iter_x <- nrow(x = mcmc$resultMCMC[[1]])
  outmcmc <- matrix(NA, ncol=nrow(data), nrow=n.iter_x)
  
  for (iter in 1:n.iter_x) {
    p <- c(mcmc$resultMCMC[[1]][iter, ], 
           fixedpar)
    
    # Min <- p["Min"]
    # Max <- p["Max"]
    
    # 21/2/2020
    # p["S"] <- 1/(4*p["S"])
    
    c <- BP_flexit(x = data$distance.center, 
                par = p) # * (Max - Min) + Min
    
    outmcmc[iter, ] <- c
  }
  
  qmcmc <- apply(X = outmcmc, MARGIN = 2, FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
  colnames(qmcmc) <- data$distance.center
  
  mcmc <- modifyList(mcmc, list(quantiles=qmcmc))
  
  outmcmc <- apply(X = outmcmc, MARGIN = 1, FUN = mean)
  qmcmc <- quantile(outmcmc, probs = c(0.025, 0.5, 0.975))
  
  mcmc <- modifyList(mcmc, list(quantiles.global=qmcmc))
  
  mcmc$timestamp <- date()
  
  mcmc$resultMCMC[[1]][, "Min"] <- mcmc$resultMCMC[[1]][, "Min"]
  mcmc$resultMCMC[[1]][, "Max"] <- mcmc$resultMCMC[[1]][, "Max"]
  
  summary.table <- data.frame(mean=apply(mcmc$resultMCMC[[1]], MARGIN=2, FUN = mean), 
  se=apply(mcmc$resultMCMC[[1]], MARGIN=2, FUN = sd))
  mcmc <- modifyList(mcmc, list(summary.table=summary.table))
  
  # Je retire bone car ça prend trop de place
  mcmc$parametersMCMC$control$bone <- NULL
  
  bone <- RM_add(x=bone, RMname = analysis, valuename = "mcmc", value=mcmc)
  
  return(bone)
}


