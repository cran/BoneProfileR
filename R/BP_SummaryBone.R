#' summary.BoneProfileR displays a bone section
#' @title Plot a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return An invisible list with recorded information
#' @param object The bone image
#' @param analysis The analysis to report the compactness
#' @param periodic.angles A vector indicating which angle to report for periodic analysis
#' @param periodic.angles.replicate.CI Number of replicates to estimate CI
#' @param ... Not used
#' @description Display information of bone section
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#'  bone <- BP_OpenImage()
#'  # or 
#'  path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  summary(bone)
#'  
#'  bone <- BP_DetectBackground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectForeground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic", cut.angle = 60)
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic", twosteps=TRUE)
#'  plot(bone, type="observations+model", analysis="logistic")
#'  par <- BP_GetFittedParameters(bone, analysis="logistic", ML=TRUE, return.all=FALSE)[, "mean"]
#'  options(mc.cores=parallel::detectCores())
#'  
#'  #############################################
#'  # Periodic analysis
#'  #############################################
#'  bone <- BP_FitMLPeriodicCompactness(bone, analysis="logistic", control.optim=list(trace=2), 
#'                                      fitted.parameters=c(par, PSin=0.001, PCos=0.001, 
#'                                      SSin=0.001, SCos=0.001, MinSin=0.001, MinCos=0.001, 
#'                                      MaxSin=0.001, MaxCos=0.001), replicates.CI=2000)
#'  bone <- BP_FitBayesianPeriodicCompactness(bone, analysis="logistic", replicates.CI=2000)
#'  
#'  summary(object=bone, analysis="logistic")
#'  summary(object=bone, analysis="logistic", 
#'          periodic.angles=seq(from=-0.1, to=0.1, length.out=10))
#' }
#' @method summary BoneProfileR
#' @export


summary.BoneProfileR <- function(object, analysis=1, 
                                 periodic.angles="all", 
                                 periodic.angles.replicate.CI=2000, 
                                 ...) {
  mc.cores <- getOption("mc.cores", parallel::detectCores())
  forking <- getOption("forking", ifelse(.Platform$OS.type == "windows", FALSE, TRUE))
  
  if (all(is.numeric(periodic.angles)) & (length(periodic.angles) == 1)) {
    stop("You must give at least two numeric values for periodic.angles.")
  }
  
  out <- dim(object)
  cat(paste0("The image is ", as.character(out[1]), " pixels width and ", as.character(out[2]), " pixels height.\n"))
  cat("__________________________________\n")
  an <- BP_ListAnalyses(object, max.level=FALSE, silent = TRUE)
  if (length(an) == 0) {
    cat("There are no recorded analysis still.\n")
  } else {
    if (length(an) >1)
      cat("Recorded analyses:\n")
    else
      cat("Recorded analysis:\n")
    an <- BP_ListAnalyses(object, max.level=FALSE, silent = FALSE)
  }
  out <- list(dim=out, analysis=an)
  cat("__________________________________\n")
  if (is.numeric(analysis)) analysis <- names(BP_ListAnalyses(object, max.level=FALSE, silent = TRUE))[analysis]
  
  if (!is.null(RM_get(x=object, RMname = analysis, valuename = "global.compactness")))
    cat(paste0("The observed global compactness for '", analysis, "' analysis is ", 
               specify_decimal(RM_get(x=object, RMname = analysis, valuename = "global.compactness"), decimals = 3), ".\n"))
  
  if (!is.null(RM_get(x=object, RMname = analysis, valuename = "optim")$quantiles["50%", ])) {
    cat(paste0("The median value for '", analysis, "' analysis for ML modeled global compactness is ", 
               specify_decimal(mean(RM_get(x=object, RMname = analysis, valuename = "optim")$quantiles["50%", ]), decimals = 3), ".\n"))
    cat(paste0("The 95% confidence interval for '", analysis, "' analysis for ML modeled global compactness is between ", 
               specify_decimal(mean(RM_get(x=object, RMname = analysis, valuename = "optim")$quantiles["2.5%", ]), decimals = 3), " and ",
               specify_decimal(mean(RM_get(x=object, RMname = analysis, valuename = "optim")$quantiles["97.5%", ]), decimals = 3),".\n"))
  }
  
  if (!is.null(RM_get(x=object, RMname = analysis, valuename = "mcmc"))) {
    mcmc <- RM_get(x=object, RMname = analysis, valuename = "mcmc")$quantiles.global
    
    cat(paste0("The median value for '", analysis, "' analysis for MCMC modeled global compactness is ", 
               specify_decimal(mcmc["50%"], decimals = 3), ".\n"))
    cat(paste0("The 95% credible interval for '", analysis, "' analysis for MCMC modeled global compactness is between ", 
               specify_decimal(mcmc["2.5%"], decimals = 3), " and ",
               specify_decimal(mcmc["97.5%"], decimals = 3),".\n"))
  }
  
  if (!is.null(RM_get(x=object, RMname = analysis, valuename = "optimRadial"))) {
    
    radial <- RM_get(x=object, RMname = analysis, valuename = "optimRadial")$synthesis
    roc <- quantile(radial$observed.compactness, probs=c(0.025, 0.5, 0.975))
    cat(paste0("The median value for '", analysis, "' analysis for ML modeled radial compactness is ", 
               specify_decimal(roc["50%"], decimals = 3), ".\n"))
    cat(paste0("The 95% confidence interval for '", analysis, "' analysis for ML modeled radial compactness is between ", 
               specify_decimal(roc["2.5%"], decimals = 3), " and ",
               specify_decimal(roc["97.5%"], decimals = 3),".\n"))
    
    rlmc <- quantile(radial$linearized.modeled.compactness, probs=c(0.025, 0.5, 0.975))
    cat(paste0("The median value for '", analysis, "' analysis for ML linearized modeled radial compactness is ", 
               specify_decimal(rlmc["50%"], decimals = 3), ".\n"))
    cat(paste0("The 95% confidence interval for '", analysis, "' analysis for ML linearized modeled radial compactness is between ", 
               specify_decimal(rlmc["2.5%"], decimals = 3), " and ",
               specify_decimal(rlmc["97.5%"], decimals = 3),".\n"))
  }
  
  if (!is.null(RM_get(x=object, RMname = analysis, valuename = "optimPeriodic"))) {
    periodic <- RM_get(x=object, RMname = analysis, valuename = "optimPeriodic")
    
    if (!is.null(periodic$PeriodicGlobalCompactness)) {
      cat(paste0("The median value for '", analysis, "' analysis for ML modeled periodic compactness is ", 
                 specify_decimal(periodic$PeriodicGlobalCompactness["50%"], decimals = 3), ".\n"))
      cat(paste0("The 95% confidence interval for '", analysis, "' analysis for ML modeled periodic compactness is between ", 
                 specify_decimal(periodic$PeriodicGlobalCompactness["2.5%"], decimals = 3), " and ",
                 specify_decimal(periodic$PeriodicGlobalCompactness["97.5%"], decimals = 3),".\n"))
    }
    
      cat(paste0("The median value for '", analysis, "' analysis for ML linearized modeled periodic compactness is ", 
                 specify_decimal(periodic$PeriodicLinearCompactness["50%"], decimals = 3), ".\n"))
      cat(paste0("The 95% confidence interval for '", analysis, "' analysis for ML linearized modeled periodic compactness is between ", 
                 specify_decimal(periodic$PeriodicLinearCompactness["2.5%"], decimals = 3), " and ",
                 specify_decimal(periodic$PeriodicLinearCompactness["97.5%"], decimals = 3),".\n"))
    if (is.numeric(periodic.angles)) {
      cut.angle <- periodic.angles
      
      rd <- RandomFromHessianOrMCMC(method = "hessian", 
                                    Hessian = periodic$hessian, 
                                    fitted.parameters = periodic$par, 
                                    fixed.parameters = periodic$fixed.parameters,
                                    replicates = periodic.angles.replicate.CI, 
                                    probs = NULL, silent = TRUE)$random
      tablePeriodic <- getFromNamespace(".tablePeriodic", ns="BoneProfileR")
      cut.distance.center <- RM_get(x=object, RMname=analysis, valuename = "cut.distance.center")
      
      out_l <- universalmclapply(X=1:periodic.angles.replicate.CI, mc.cores =  mc.cores, forking = forking, 
                                 FUN = function(i) {
                                   p <- unlist(rd[i, , drop=TRUE])
                                   
                                   par_ec <- tablePeriodic(par=p                                    ,
                                                           cut.distance.center=cut.distance.center  , 
                                                           cut.angle=cut.angle                      , 
                                                           array.compactness=NA                     , 
                                                           analysis=NULL                            , 
                                                           fixed.parameters=NULL                    , 
                                                           sp=NULL                                  , 
                                                           cp=NULL                                  )
                                   # -sum(dbinom(x=par_ec[, "Os"], size=par_ec[, "Pixels"], prob=par_ec[, "Compactness"], log=TRUE))
                                   return(par_ec[, "Compactness"])
                                 }, 
                                 progressbar = TRUE, 
                                 clusterExport = list(tablePeriodic, rd, cut.angle, cut.distance.center, analysis), 
                                 clusterEvalQ=list(expr=as.expression("library(HelpersMG)")))
      
      out <- array(data = NA, dim=c( length(cut.angle)-1, 
                                     length(cut.distance.center)-1, 
                                     periodic.angles.replicate.CI), 
                   dimnames = list(paste0("]", specify_decimal(cut.angle[-length(cut.angle)], 2), ", ", specify_decimal(cut.angle[-1], 2), "]"), 
                                   paste0("]", specify_decimal(cut.distance.center[-length(cut.distance.center)], 2), ", ", specify_decimal(cut.distance.center[-1], 2), "]"), 
                                   as.character(1:periodic.angles.replicate.CI)))
      for (i in 1:length(out_l)) {
        out[, , i] <- matrix(out_l[[i]], nrow=length(cut.angle)-1, ncol=length(cut.distance.center)-1, byrow = TRUE)
      }
      
      out <- ifelse(out < 1E-10, 1E-10, out)
      out <- ifelse(out > 1-1E-10, 1-1E-10, out)
      out_l <- quantile(unlist(lapply(X = 1:dim(out)[3], FUN=function(i) mean(out[, , i]))), probs=c(0.025, 0.5, 0.975))
      
      cat(paste0("Periodic analysis for angles ", paste(specify_decimal(periodic.angles, decimals = 3), collapse = "; "), ".\n"))
      cat(paste0("The median value for '", analysis, "' analysis for ML linearized modeled periodic compactness is ", 
                 specify_decimal(out_l["50%"], decimals = 3), ".\n"))
      cat(paste0("The 95% confidence interval for '", analysis, "' analysis for ML linearized modeled periodic compactness is between ", 
                 specify_decimal(out_l["2.5%"], decimals = 3), " and ",
                 specify_decimal(out_l["97.5%"], decimals = 3),".\n"))
    }
  }
  
  if (!is.null(RM_get(x=object, RMname = analysis, valuename = "mcmcPeriodic"))) {
    periodic <- RM_get(x=object, RMname = analysis, valuename = "mcmcPeriodic")
    if (!is.null(periodic$PeriodicGlobalCompactness)) {
      cat(paste0("The median value for '", analysis, "' analysis for Bayesian modeled periodic compactness is ", 
                 specify_decimal(periodic$PeriodicGlobalCompactness["50%"], decimals = 3), ".\n"))
      cat(paste0("The 95% credible interval for '", analysis, "' analysis for Bayesian modeled periodic compactness is between ", 
                 specify_decimal(periodic$PeriodicGlobalCompactness["2.5%"], decimals = 3), " and ",
                 specify_decimal(periodic$PeriodicGlobalCompactness["97.5%"], decimals = 3),".\n"))
    }
    
    cat(paste0("The median value for '", analysis, "' analysis for Bayesian linearized modeled periodic compactness is ", 
               specify_decimal(periodic$PeriodicLinearCompactness["50%"], decimals = 3), ".\n"))
    cat(paste0("The 95% credible interval for '", analysis, "' analysis for Bayesian linearized modeled periodic compactness is between ", 
               specify_decimal(periodic$PeriodicLinearCompactness["2.5%"], decimals = 3), " and ",
               specify_decimal(periodic$PeriodicLinearCompactness["97.5%"], decimals = 3),".\n"))
    
    if (is.numeric(periodic.angles)) {
      cut.angle <- periodic.angles
      
      rd <- RandomFromHessianOrMCMC(method = "MCMC", 
                                    mcmc = periodic, 
                                    fixed.parameters = periodic$fixed.parameters,
                                    replicates = periodic.angles.replicate.CI, 
                                    probs = NULL, silent = FALSE)$random
      
      
      tablePeriodic <- getFromNamespace(".tablePeriodic", ns="BoneProfileR")
      cut.distance.center <- RM_get(x=object, RMname=analysis, valuename = "cut.distance.center")
      
      out_l <- universalmclapply(X=1:periodic.angles.replicate.CI, mc.cores =  mc.cores, forking = forking, 
                                 FUN = function(i) {
                                   p <- unlist(rd[i, , drop=TRUE])
                                   
                                   par_ec <- tablePeriodic(par=p                                    ,
                                                           cut.distance.center=cut.distance.center  , 
                                                           cut.angle=cut.angle                      , 
                                                           array.compactness=NA                     , 
                                                           analysis=NULL                            , 
                                                           fixed.parameters=NULL                    , 
                                                           sp=NULL                                  , 
                                                           cp=NULL                                  )
                                   # -sum(dbinom(x=par_ec[, "Os"], size=par_ec[, "Pixels"], prob=par_ec[, "Compactness"], log=TRUE))
                                   return(par_ec[, "Compactness"])
                                 }, 
                                 progressbar = TRUE, 
                                 clusterExport = list(tablePeriodic, rd, cut.angle, cut.distance.center, analysis), 
                                 clusterEvalQ=list(expr=as.expression("library(HelpersMG)")))
      
      out <- array(data = NA, dim=c( length(cut.angle)-1, 
                                     length(cut.distance.center)-1, 
                                     periodic.angles.replicate.CI), 
                   dimnames = list(paste0("]", specify_decimal(cut.angle[-length(cut.angle)], 2), ", ", specify_decimal(cut.angle[-1], 2), "]"), 
                                   paste0("]", specify_decimal(cut.distance.center[-length(cut.distance.center)], 2), ", ", specify_decimal(cut.distance.center[-1], 2), "]"), 
                                   as.character(1:periodic.angles.replicate.CI)))
      for (i in 1:length(out_l)) {
        out[, , i] <- matrix(out_l[[i]], nrow=length(cut.angle)-1, ncol=length(cut.distance.center)-1, byrow = TRUE)
      }
      
      out <- ifelse(out < 1E-10, 1E-10, out)
      out <- ifelse(out > 1-1E-10, 1-1E-10, out)
      out_l <- quantile(unlist(lapply(X = 1:dim(out)[3], FUN=function(i) mean(out[, , i]))), probs=c(0.025, 0.5, 0.975))
      
      cat(paste0("Periodic analysis for angles ", paste(specify_decimal(periodic.angles, decimals = 3), collapse = "; "), ".\n"))
      cat(paste0("The median value for '", analysis, "' analysis for MCMC linearized modeled periodic compactness is ", 
                 specify_decimal(out_l["50%"], decimals = 3), ".\n"))
      cat(paste0("The 95% confidence interval for '", analysis, "' analysis for MCMC linearized modeled periodic compactness is between ", 
                 specify_decimal(out_l["2.5%"], decimals = 3), " and ",
                 specify_decimal(out_l["97.5%"], decimals = 3),".\n"))
    }
    
  }
  
  
  return(invisible(out))
}


