#' BP_FitBayesianPeriodicCompactness estimates likelihood of global model of a bone section
#' @title Estimation of the likelihood of a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The -Ln L
#' @param bone The bone image to be used
#' @param fitted.parameters Parameters of the model to be fitted
#' @param fixed.parameters Fixed parameters of the model
#' @param analysis Name or rank of analysis
#' @param priors The priors of Bayesian analysis
#' @param replicates.CI Number of replicates to estimate confidence interval using Hessian
#' @param amplitude.max The maximum allowed amplitude for each parameter
#' @param control.MHalgoGen The control parameters of MHalgoGen()
#' @param silent Should the function displays some information?
#' @description Estimation of the compactness of a bone section using Bayesian periodic model.\cr
#' To control the parallel computing, use: \cr
#' options(mc.cores = \[put here the number of cores you want use\])\cr
#' options(forking = FALSE) or options(forking = TRUE)\cr
#' The maximum number of cores is obtained by: parallel::detectCores()\cr
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run
#' library(BoneProfileR)
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
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
#'  mcmc <- RM_get(bone, RMname="logistic", valuename="mcmcPeriodic")
#'  plot(mcmc, parameters="P", what="MarkovChain", ylim=c(0.555, 0.565), main="P parameter")
#'  
#'  plot(bone, type="mcmcPeriodic", parameter.name="compactness", col=rainbow(128))
#'  plot(bone, type="mcmcPeriodic", parameter.name="compactness", 
#'                col=hcl.colors(12, "YlOrRd", rev = TRUE))
#'  plot(bone, type="mcmcPeriodic", parameter.name="averagemodel")
#'  plot(bone, type="mcmcPeriodic", parameter.name="P", 
#'                rgb(red = 0.7, green = 0.7, blue = 0.7, alpha = 0.2))
#'  plot(bone, type="mcmcPeriodic", parameter.name="P", ylim=c(0, 1), 
#'                rgb(red = 0.7, green = 0.7, blue = 0.7, alpha = 0.2))
#'  
#' }
#' @export


BP_FitBayesianPeriodicCompactness <- function(bone                                      , 
                                              fitted.parameters=NULL                    , 
                                              priors=NULL                               , 
                                              fixed.parameters=NULL                     , 
                                              analysis=1                                , 
                                              silent=FALSE                              , 
                                              replicates.CI=2000                        , 
                                              amplitude.max=0.1                         ,
                                              control.MHalgoGen=list(n.iter = 10000,    
                                                                     n.chains = 1,      
                                                                     trace = TRUE,      
                                                                     n.adapt = 5000,    
                                                                     thin = 1,          
                                                                     adaptive = TRUE)     
                                              ) {
  
  
  if (is.null(analysis)) {
    stop("You must choose which analysis to analyse.")
  }
  
  mc.cores <- getOption("mc.cores", parallel::detectCores())
  forking <- getOption("forking", ifelse(.Platform$OS.type == "windows", FALSE, TRUE))
  
  
  out <- RM_list(x=bone, silent=TRUE)
  if (is.numeric(analysis))
    if (analysis > length(out)) {
      stop("The analysis does no exist.")
    } else {
      analysis <- names(out)[analysis]
    }
  
  if (is.character(analysis) & (all(analysis != names(out)))) {
    stop(paste("The analysis", analysis, "does not exit. Check your data."))
  }
  
  o <- RM_get(x=bone, RMname = analysis, valuename="optimPeriodic")
  
  if (is.null(o)) {
    stop("You must do first a periodic analysis with BP_FitMLPeriodicCompactness() !")
  }
  
  fitted.parameters <- BP_GetFittedParameters(bone, analysis=analysis, type="periodic", ML=TRUE, return.all=FALSE)[, "mean"]
  fixed.parameters <- BP_GetFittedParameters(bone, analysis=analysis, type="periodic", ML=TRUE, return.all = TRUE)$fixed.parameters
  
  lnLG <- getFromNamespace(".lnLG", ns="BoneProfileR")
  tablePeriodic <- getFromNamespace(".tablePeriodic", ns="BoneProfileR")
  
  array.compactness <- RM_get(x=bone, RMname=analysis, valuename = "array.compactness")
  Os <- as.vector(t(array.compactness[, , 2]))
  Pixels <- as.vector(t(array.compactness[, , 1])) + Os
  
  cut.distance.center <- RM_get(x=bone, RMname=analysis)$cut.distance.center
  cut.angle <- RM_get(x=bone, RMname=analysis)$cut.angle
  
  # sin((jour/365.25)*2*pi)+ cos((jour/365.25)*2*pi)
  dc <- (cut.distance.center[-length(cut.distance.center)] + cut.distance.center[-1])/2
  da <- (cut.angle[-length(cut.angle)] + cut.angle[-1])/2
  par_ec <- expand.grid(dc, da)
  
  sp <- sin(par_ec$Var2)
  cp <- cos(par_ec$Var2)
  V1 <- par_ec$Var1
  
  rules <- rbind(data.frame(Name="P", Min=0, Max=1),
                 data.frame(Name="PSin", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="PCos", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="S", Min=-3, Max=3),
                 data.frame(Name="SSin", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="SCos", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="Min", Min=0, Max=0.7),
                 data.frame(Name="MinSin", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="MinCos", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="Max", Min=0.7, Max=1),
                 data.frame(Name="MaxSin", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="MaxCos", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="K1", Min=-1000, Max=1000),
                 data.frame(Name="K1Sin", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="K1Cos", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="K2", Min=-1000, Max=1000),
                 data.frame(Name="K2Sin", Min=-amplitude.max, Max=amplitude.max),
                 data.frame(Name="K2Cos", Min=-amplitude.max, Max=amplitude.max))
  
  
  priors <- setPriors(
    par = fitted.parameters,
    se = NULL,
    density = "dunif",
    rules = rules,
    silent = FALSE
  )
  
  priors[, "SDProp"] <- priors[, "SDProp"]/10
  
  message("MCMC estimation")
  
  o <- do.call("MHalgoGen", args = modifyList(list(likelihood = lnLG,
                                                      fixed.parameters=fixed.parameters,
                                                      Os=Os, Pixels=Pixels, 
                                                      sp=sp, 
                                                      cp=cp, 
                                                      V1=V1, 
                                                      sign = -1,
                                                      parameters_name = "par",
                                                      parameters = priors), control.MHalgoGen))

  fitted.parameters <- as.parameters(o, index="median")
  par_ec <- tablePeriodic(par=fitted.parameters, bone=bone, analysis=analysis, 
                          fixed.parameters=fixed.parameters)
  
  o$TableCompactness <- par_ec
  
  o$fixed.parameters <- fixed.parameters
  par <- c(o$par, fixed.parameters)
  
  o$summary.table <- data.frame(mean=apply(o$resultMCMC[[1]], MARGIN=2, FUN = mean), 
                              se=apply(o$resultMCMC[[1]], MARGIN=2, FUN = sd))
  
  o$par <- fitted.parameters
  
  rd <- RandomFromHessianOrMCMC(method = "MCMC", 
                                mcmc = o, 
                                fixed.parameters = o$fixed.parameters,
                                replicates = replicates.CI, 
                                probs = NULL, silent = silent)$random
  
  out_l <- universalmclapply(X=1:replicates.CI, mc.cores =  mc.cores, forking = forking, 
                             FUN = function(i) {
    p <- unlist(rd[i, , drop=TRUE])
    
    par_ec <- tablePeriodic(par=p,
                            cut.distance.center=cut.distance.center, 
                            cut.angle=cut.angle, 
                            array.compactness=array.compactness, 
                            analysis=analysis, 
                            fixed.parameters=fixed.parameters, 
                            sp=sp, cp=cp)
    # -sum(dbinom(x=par_ec[, "Os"], size=par_ec[, "Pixels"], prob=par_ec[, "Compactness"], log=TRUE))
    return(par_ec[, "Compactness"])
  }, 
  progressbar = TRUE, 
  clusterExport = list(tablePeriodic, rd, cut.angle, cut.distance.center, array.compactness, analysis, fixed.parameters), 
  clusterEvalQ=list(expr=as.expression("library(HelpersMG)"))
  )
  
  out <- array(data = NA, dim=c( length(cut.angle)-1, 
                                 length(cut.distance.center)-1, 
                                 replicates.CI), 
               dimnames = list(paste0("]", specify_decimal(cut.angle[-length(cut.angle)], 2), ", ", specify_decimal(cut.angle[-1], 2), "]"), 
                               paste0("]", specify_decimal(cut.distance.center[-length(cut.distance.center)], 2), ", ", specify_decimal(cut.distance.center[-1], 2), "]"), 
                               as.character(1:replicates.CI)))
  for (i in 1:length(out_l)) {
    out[, , i] <- matrix(out_l[[i]], nrow=length(cut.angle)-1, ncol=length(cut.distance.center)-1, byrow = TRUE)
  }
  
  out <- ifelse(out < 1E-10, 1E-10, out)
  out <- ifelse(out > 1-1E-10, 1-1E-10, out)
  
  rm(out_l)
  
  o$PeriodicLinearCompactness <- quantile(unlist(lapply(X = 1:dim(out)[3], FUN=function(i) mean(out[, , i]))), probs=c(0.025, 0.5, 0.975))
  
  nbpixels <- array.compactness[, , 1]+array.compactness[,,2]
  nbpixels <- nbpixels/(sum(nbpixels))
  
  if (all(dim(out)[1:2] == dim(nbpixels))) o$PeriodicGlobalCompactness <- quantile(unlist(lapply(X = 1:dim(out)[3], FUN=function(i) sum(nbpixels*(out[, , i])))), probs=c(0.025, 0.5, 0.975))
  
  
  out_q <- array(data = NA, dim=c(length(cut.angle)-1, 
                                  length(cut.distance.center)-1, 
                                  3), 
                 dimnames = list(paste0("]", specify_decimal(cut.angle[-length(cut.angle)], 2), ", ", specify_decimal(cut.angle[-1], 2), "]"), 
                                 paste0("]", specify_decimal(cut.distance.center[-length(cut.distance.center)], 2), ", ", specify_decimal(cut.distance.center[-1], 2), "]"), 
                                 c("2.5%", "50%", "97.5%")))
  
  for (nr in 1:dim(out_q)[1]) for(nc in 1:dim(out_q)[2]) out_q[nr, nc, 1:3] <- quantile(out[nr, nc, ], probs=c(0.025, 0.50, 0.975))
  
  o$PeriodicCompactness <- out_q
  
  out_gc <- array(data = NA, dim=c(length(cut.distance.center)-1, 
                                   3), 
                  dimnames = list(paste0("]", specify_decimal(cut.distance.center[-length(cut.distance.center)], 2), ", ", specify_decimal(cut.distance.center[-1], 2), "]"), 
                                  c("2.5%", "50%", "97.5%")))
  
  for(nc in 1:dim(out_gc)[1]) out_gc[nc, 1:3] <- quantile(out[, nc, ], probs=c(0.025, 0.50, 0.975))
  
  o$GlobalCompactness <- out_gc
  
  bone <- RM_add(x=bone, RMname = analysis, valuename = "mcmcPeriodic", value=o)
  
  return(bone)
}


