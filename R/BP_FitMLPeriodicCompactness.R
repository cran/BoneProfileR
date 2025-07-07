#' BP_FitMLPeriodicCompactness estimates likelihood of global model of a bone section
#' @title Estimation of the likelihood of a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The -Ln L
#' @param bone The bone image to be used
#' @param fitted.parameters Parameters of the model to be fitted
#' @param fixed.parameters Fixed parameters of the model
#' @param analysis Name or rank of analysis
#' @param twosteps Should a 2-steps analysis be performed? It can be sometimes useful.
#' @param replicates.CI Number of replicates to estimate confidence interval using Hessian
#' @param amplitude.max The maximum allowed amplitude for each parameter
#' @param control.optim The list of options for optim.
#' @param silent Should the function displays some information?
#' @description Estimation of the compactness of a bone section using radial model.\cr
#' If the fitted.parameters and fixed.parameters are NULL and the analysis includes a 
#' BP_FitMLCompactness() result, the values of this result is used as a reference for 
#' fitted.parameters and fixed.parameters.\cr
#' If no BP_FitMLCompactness() result is available, it will use:\cr
#' fitted.parameters=c(P=0.5, S=0.05, Min=-2, Max=5); fixed.parameters=c(K1=1, K2=1).\cr
#' The reference for radial estimation of compactness is the trigonometric circle for rotation.angle=0 in 
#' BP_EstimateCompactness():\cr
#' - The top of the section is located at -pi/2.\cr
#' - The left of the section is located at -pi and +pi.\cr
#' - The bottom of the section is located at pi/2.\cr
#' - The right of the section is 0.\cr
#' If rotation.angle is different from 0, the value of rotation.angle is added to the angle modulo 2.pi.\cr
#' The two-steps analysis performs first a quasi-Newton method, then a Bayesian MCMC and finally again a quasi-Newton method. 
#' It generally ensures that global minimum is found. On the other hand, it doubles the time to complete for each angle.\cr
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
#'  analysisP <- BP_GetFittedParameters(bone, analysis="logistic", type="periodic", 
#'                                      ML=TRUE, return.all=FALSE)[, "mean"]
#'  analysisP$par                                    
#'  plot(bone, type="periodic", parameter.name="compactness", col=rainbow(128))
#'  plot(bone, type="periodic", parameter.name="compactness", 
#'                col=hcl.colors(12, "YlOrRd", rev = TRUE))
#'  plot(bone, type="periodic", parameter.name="averagemodel")
#'  plot(bone, type="periodic", parameter.name="P", 
#'                rgb(red = 0.7, green = 0.7, blue = 0.7, alpha = 0.2))
#'  plot(bone, type="periodic", parameter.name="P", ylim=c(0, 1), 
#'                rgb(red = 0.7, green = 0.7, blue = 0.7, alpha = 0.2))
#'  boneNoPeriodic <- BP_FitMLPeriodicCompactness(bone, analysis="logistic", 
#'                                                control.optim=list(trace=2), 
#'                                      fitted.parameters=par, replicates.CI=2000)
#'  analysisNP <- BP_GetFittedParameters(boneNoPeriodic, analysis="logistic", ML=TRUE, 
#'                                       return.all=TRUE, type="periodic")
#'  analysisNP$par
#'  compare_AIC(PeriodicModel=analysisP, 
#'              NoPeriodicModel=analysisNP)
#'  
#'  #############################################
#'  
#'  # Note that the absolute likelihood is dependent on the number of angle cut
#'  # Only models analyzed with the same number of angle cuts can be compared
#'  
#'  dbinom(5, 10, prob=0.4, log=TRUE); 
#'        dbinom(2, 5, prob=0.4, log=TRUE)+dbinom(3, 5, prob=0.4, log=TRUE)
#'  # But the likelihood difference between two models are not:
#'  dbinom(5, 10, prob=0.4, log=TRUE)-dbinom(5, 10, prob=0.3, log=TRUE)
#'  dbinom(2, 5, prob=0.4, log=TRUE)+dbinom(3, 5, prob=0.4, log=TRUE)- 
#'       dbinom(2, 5, prob=0.3, log=TRUE)-dbinom(3, 5, prob=0.3, log=TRUE)
#' }
#' @export


BP_FitMLPeriodicCompactness <- function(bone                        , 
                                        fitted.parameters=NULL      , 
                                        # priors=NULL               , 
                                        fixed.parameters=NULL       , 
                                        analysis=1                  , 
                                        silent=FALSE                , 
                                        replicates.CI=NULL          , 
                                        twosteps=FALSE              , 
                                        amplitude.max=0.1           ,
                                        control.optim=list(trace=1) ) {
  
  # fitted.parameters=c(P=0.5, S=0.05, Min=0.001, Max=0.999); fixed.parameters=c(K1=1, K2=1); analysis=NULL; silent=FALSE; twosteps=TRUE
  # fitted.parameters=NULL; priors=NULL; fixed.parameters=NULL; analysis=1; silent=FALSE; twosteps=TRUE
  
  mc.cores <- getOption("mc.cores", parallel::detectCores())
  forking <- getOption("forking", ifelse(.Platform$OS.type == "windows", FALSE, TRUE))
  
  
  if (is.null(analysis)) {
    stop("You must choose which analysis to report.")
  }
  
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
  
  if (is.null(fitted.parameters)) {
    if (is.null(BP_GetFittedParameters(bone, analysis=analysis, type = "global", ML=TRUE, return.all = FALSE))) {
      fitted.parameters=c(P=0.5, S=0.05, Min=0.001, Max=0.999, 
                          PSin=0, PCos=0, 
                          SSin=0, SCos=0, 
                          MinSin=0, MinCos=0, 
                          MaxSin=0, MaxCos=0)
      fixed.parameters=c(K1=1, K1Sin=0, K1Cos=0, K2=1, K2Sin=0, K2Cos=0)
    } else {
      fitted.parameters <- BP_GetFittedParameters(bone, analysis=analysis, type = "global", ML=TRUE, return.all = FALSE)[, "mean"]
      fixed.parameters <- BP_GetFittedParameters(bone, analysis=analysis, type = "global",  ML=TRUE, return.all = TRUE)$fixed.parameters
    }
  }
  
  if (is.null(fixed.parameters)) {
    fixed.parameters <- BP_GetFittedParameters(bone, analysis=analysis, type = "global",  ML=TRUE, return.all = TRUE)$fixed.parameters
  }
  
  lower_limit <- c(P=0, PSin=-amplitude.max, PCos=-amplitude.max, 
                   S=-3, SSin=-amplitude.max, SCos=-amplitude.max, 
                   Min=0, MinSin=-amplitude.max, MinCos=-amplitude.max, 
                   Max=0.2, MaxSin=-amplitude.max, MaxCos=-amplitude.max, 
                   K1=-1000, K1Sin=-amplitude.max, K1Cos=-amplitude.max, 
                   K2=-1000, K2Sin=-amplitude.max, K2Cos=-amplitude.max)
  upper_limit <- c(P=1, PSin=amplitude.max, PCos=amplitude.max, 
                   S=3, SSin=amplitude.max, SCos=amplitude.max, 
                   Min=0.8, MinSin=amplitude.max, MinCos=amplitude.max, 
                   Max=1, MaxSin=amplitude.max, MaxCos=amplitude.max, 
                   K1=1000, K1Sin=amplitude.max, K1Cos=amplitude.max, 
                   K2=1000, K2Sin=amplitude.max, K2Cos=amplitude.max)
  
  lower <- lower_limit[names(fitted.parameters)]
  upper <- upper_limit[names(fitted.parameters)]
  
  fitted.parameters[names(upper)] <- ifelse(fitted.parameters[names(upper)] >= upper, upper-0.001*upper, fitted.parameters[names(upper)])
  fitted.parameters[names(lower)] <- ifelse(fitted.parameters[names(lower)] <= lower, lower+0.001*lower, fitted.parameters[names(lower)])
  
  # lnLG(par=fitted.parameters, bone=bone, analysis=analysis, fixed.parameters=NULL, sign=1)
  
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
  
  
  if (twosteps) {
    
    o <- optim(par=fitted.parameters, fn=lnLG, 
               fixed.parameters=fixed.parameters, 
               Os=Os, Pixels=Pixels, 
               sp=sp, 
               cp=cp, 
               V1=V1, 
               upper = upper, lower = lower, 
               # bone=bone, analysis=analysis, 
               control=modifyList(list(maxit=1000), control.optim), 
               sign=-1, 
               method = "L-BFGS-B", hessian = FALSE)
    
    fitted.parameters <- o$par
    
    fitted.parameters[names(upper)] <- ifelse(fitted.parameters[names(upper)] >= upper, upper-0.01*upper, fitted.parameters[names(upper)])
    fitted.parameters[names(lower)] <- ifelse(fitted.parameters[names(lower)] <= lower, lower+0.01*lower, fitted.parameters[names(lower)])
    
    message("Second round and Hessian estimation")
    o <- optim(par=fitted.parameters, fn=lnLG, 
               fixed.parameters=fixed.parameters, 
               Os=Os, Pixels=Pixels, 
               sp=sp, 
               cp=cp, 
               V1=V1, 
               upper = upper, lower = lower, 
               # bone=bone, analysis=analysis, 
               control=modifyList(list(maxit=1000), control.optim), 
               sign=-1, 
               method = "L-BFGS-B", hessian = TRUE)
  } else {
    message("Fit and Hessian estimation")
    o <- optim(par=fitted.parameters, fn=lnLG, 
               fixed.parameters=fixed.parameters, 
               Os=Os, Pixels=Pixels, 
               sp=sp, 
               cp=cp, 
               V1=V1, 
               upper = upper, lower = lower, 
               # bone=bone, analysis=analysis, 
               control=modifyList(list(maxit=1000), control.optim), 
               sign=-1, 
               method = "L-BFGS-B", hessian = TRUE)
  }
  
  par_ec <- tablePeriodic(par=o$par, bone=bone, analysis=analysis, 
                          fixed.parameters=fixed.parameters)
  
  o$TableCompactness <- par_ec
  
  o$fixed.parameters <- fixed.parameters
  par <- c(o$par, fixed.parameters)
  o$AIC <- 2*o$value+2*length(o$par)
  o$SE <- SEfromHessian(o$hessian)
  
  
  rd <- RandomFromHessianOrMCMC(method = "hessian", 
                                Hessian = o$hessian, 
                                fitted.parameters = o$par, 
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
                             clusterExport = list(tablePeriodic, rd, cut.angle, cut.distance.center, 
                                                  array.compactness, analysis, fixed.parameters, sp, cp), 
                             clusterEvalQ=list(expr=as.expression("library(HelpersMG)")))
  
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
  
  bone <- RM_add(x=bone, RMname = analysis, valuename = "optimPeriodic", value=o)
  
  return(bone)
}

.tablePeriodic <- function(par                      , 
                           bone=NULL                , 
                           analysis=1               , 
                           cut.distance.center=NULL , 
                           cut.angle=NULL           , 
                           array.compactness=NULL   , 
                           fixed.parameters=NULL    , 
                           sp=NULL                  , 
                           cp=NULL                  ) {
  
  para <- c(par, fixed.parameters)
  
  P <- para["P"]
  PSin <- para["PSin"]; if (is.na(PSin)) PSin <- 0
  PCos <- para["PCos"]; if (is.na(PCos)) PCos <- 0
  S <- para["S"]
  SSin <- para["SSin"]; if (is.na(SSin)) SSin <- 0
  SCos <- para["SCos"]; if (is.na(SCos)) SCos <- 0
  Min <- para["Min"]; if (is.na(Min)) Min <- 1E-10
  MinSin <- para["MinSin"]; if (is.na(MinSin)) MinSin <- 0
  MinCos <- para["MinCos"]; if (is.na(MinCos)) MinCos <- 0
  Max <- para["Max"]; if (is.na(Max)) Max <- 1-1E-10
  MaxSin <- para["MaxSin"]; if (is.na(MaxSin)) MaxSin <- 0
  MaxCos <- para["MaxCos"]; if (is.na(MaxCos)) MaxCos <- 0
  K1 <- para["K1"]; if (is.na(K1)) K1 <- 1
  K1Sin <- para["K1Sin"]; if (is.na(K1Sin)) K1Sin <- 0
  K1Cos <- para["K1Cos"]; if (is.na(K1Cos)) K1Cos <- 0
  K2 <- para["K2"]; if (is.na(K2)) K2 <- 1
  K2Sin <- para["K2Sin"]; if (is.na(K2Sin)) K2Sin <- 0
  K2Cos <- para["K2Cos"]; if (is.na(K2Cos)) K2Cos <- 0
  
  if (is.null(cut.distance.center)) cut.distance.center <- HelpersMG::RM_get(x=bone, RMname=analysis)$cut.distance.center
  if (is.null(cut.angle)) cut.angle <- HelpersMG::RM_get(x=bone, RMname=analysis)$cut.angle
  
  # sin((jour/365.25)*2*pi)+ cos((jour/365.25)*2*pi)
  dc <- (cut.distance.center[-length(cut.distance.center)] + cut.distance.center[-1])/2
  da <- (cut.angle[-length(cut.angle)] + cut.angle[-1])/2
  par_ec <- expand.grid(dc, da)
  par_ec <- cbind(par_ec, expand.grid(1:length(dc), 1:length(da)))
  colnames(par_ec) <- c("Var1", "Var2", "dc", "da")
  
  if (is.null(sp)) sp <- sin(par_ec$Var2)
  if (is.null(cp)) cp <- cos(par_ec$Var2)
  
  par_ec <- cbind(par_ec, P=P + PSin*sp + PCos*cp)
  par_ec <- cbind(par_ec, S=S + SSin*sp + SCos*cp)
  par_ec <- cbind(par_ec, Min=Min + MinSin*sp + MinCos*cp)
  par_ec <- cbind(par_ec, Max=Max + MaxSin*sp + MaxCos*cp)
  par_ec <- cbind(par_ec, K1=K1 + K1Sin*sp + K1Cos*cp)
  par_ec <- cbind(par_ec, K2=K2 + K2Sin*sp + K2Cos*cp)
  
  par_ec <- cbind(par_ec, Compactness=BoneProfileR::BP_flexit(x = par_ec$Var1, P=par_ec$P, S=par_ec$S, K1=par_ec$K1, 
                                                              K2=par_ec$K2, Min=par_ec$Min, Max=par_ec$Max))
  
  if (length(array.compactness) != 1) {
  if (is.null(array.compactness)) array.compactness <- HelpersMG::RM_get(x=bone, RMname=analysis, valuename = "array.compactness")
  
  # par_ec <- cbind(par_ec, Os=array.compactness[par_ec[, "da"], par_ec[, "dc"], 1])
  par_ec <- cbind(par_ec, Os=as.vector(t(array.compactness[, , 2])))
  par_ec <- cbind(par_ec, Pixels=as.vector(t(array.compactness[, , 1])) + par_ec[, "Os"])
  }
  return(par_ec)
}


.lnLG <- function(par                  , 
                  fixed.parameters=NULL, 
                  Os=NULL              ,
                  Pixels = NULL        ,
                  sp=NULL              , 
                  cp=NULL              , 
                  V1=NULL              ,
                  sign = -1) {
  
  # bone=NULL; analysis=1; cut.distance.center=NULL; cut.angle=NULL; array.compactness=NULL; fixed.parameters=NULL; sign = -1
  # par <- c(P=0.3, S=1, Min=0.1, Max=0.9); fixed.parameters=NULL
  
  para <- c(par, fixed.parameters)
  
  P <- para["P"]
  PSin <- para["PSin"]; if (is.na(PSin)) PSin <- 0
  PCos <- para["PCos"]; if (is.na(PCos)) PCos <- 0
  S <- para["S"]
  SSin <- para["SSin"]; if (is.na(SSin)) SSin <- 0
  SCos <- para["SCos"]; if (is.na(SCos)) SCos <- 0
  Min <- para["Min"]; if (is.na(Min)) Min <- 1E-10
  MinSin <- para["MinSin"]; if (is.na(MinSin)) MinSin <- 0
  MinCos <- para["MinCos"]; if (is.na(MinCos)) MinCos <- 0
  Max <- para["Max"]; if (is.na(Max)) Max <- 1-1E-10
  MaxSin <- para["MaxSin"]; if (is.na(MaxSin)) MaxSin <- 0
  MaxCos <- para["MaxCos"]; if (is.na(MaxCos)) MaxCos <- 0
  K1 <- para["K1"]; if (is.na(K1)) K1 <- 1
  K1Sin <- para["K1Sin"]; if (is.na(K1Sin)) K1Sin <- 0
  K1Cos <- para["K1Cos"]; if (is.na(K1Cos)) K1Cos <- 0
  K2 <- para["K2"]; if (is.na(K2)) K2 <- 1
  K2Sin <- para["K2Sin"]; if (is.na(K2Sin)) K2Sin <- 0
  K2Cos <- para["K2Cos"]; if (is.na(K2Cos)) K2Cos <- 0
  
  par_ec <- (Pixels != 0)
  
  P <- (P + PSin*sp + PCos*cp)[par_ec]
  S <- (S + SSin*sp + SCos*cp)[par_ec]
  Min <- (Min + MinSin*sp + MinCos*cp)[par_ec]
  Max <- (Max + MaxSin*sp + MaxCos*cp)[par_ec]
  K1 <- (K1 + K1Sin*sp + K1Cos*cp)[par_ec]
  K2 <- (K2 + K2Sin*sp + K2Cos*cp)[par_ec]
  V1 <- V1[par_ec]
  
  # Compactness <- unlist(parallel::mclapply(X = seq_along(V1), function(i) BP_flexit(x = V1[i], P=P[i], S=S[i], K1=K1[i], K2=K2[i], Min=Min[i], Max=Max[i])))
  
  Compactness <- BP_flexit(x = V1, P=P, S=S, K1=K1, K2=K2, Min=Min, Max=Max)
  
  Compactness <- ifelse(Compactness < 1E-10, 1E-10, Compactness)
  Compactness <- ifelse(Compactness > 1-1E-10, 1-(1E-10), Compactness)
  
  LnL <- sign * sum(dbinom(x=Os[par_ec], size=Pixels[par_ec], prob = Compactness, log=TRUE))
  if (is.infinite(LnL)) LnL <- sign * -1E6
  # if (!silent) print(LnL)
  # if (!silent) print(d(par))
  return(LnL)
}

# rules <- rbind(data.frame(Name="P", Min=0, Max=1), 
#                data.frame(Name="PSin", Min=-0.1, Max=0.1), 
#                data.frame(Name="PCos", Min=-0.1, Max=0.1), 
#                data.frame(Name="S", Min=-3, Max=3), 
#                data.frame(Name="SSin", Min=-0.1, Max=0.1), 
#                data.frame(Name="SCos", Min=-0.1, Max=0.1), 
#                data.frame(Name="Min", Min=0, Max=0.7), 
#                data.frame(Name="MinSin", Min=-0.1, Max=0.1), 
#                data.frame(Name="MinCos", Min=-0.1, Max=0.1), 
#                data.frame(Name="Max", Min=0.7, Max=1), 
#                data.frame(Name="MaxSin", Min=-0.1, Max=0.1), 
#                data.frame(Name="MaxCos", Min=-0.1, Max=0.1), 
#                data.frame(Name="K1", Min=-1000, Max=1000), 
#                data.frame(Name="K1Sin", Min=-0.1, Max=0.1), 
#                data.frame(Name="K1Cos", Min=-0.1, Max=0.1), 
#                data.frame(Name="K2", Min=-1000, Max=1000), 
#                data.frame(Name="K2Sin", Min=-0.1, Max=0.1), 
#                data.frame(Name="K2Cos", Min=-0.1, Max=0.1))
# 
# 
# priors <- setPriors(
#   par = o$par,
#   se = NULL,
#   density = "dunif",
#   rules = rules,
#   silent = FALSE
# )
# priors[, "SDProp"] <- priors[, "SDProp"]/10
# message("MCMC estimation")
# mcmc <- HelpersMG::MHalgoGen(
#   likelihood = lnLG,
#   bone=bone, 
#   fixed.parameters=fixed.parameters, 
#   sign = -1, 
#   parameters_name = "par", 
#   parameters = priors, 
#   n.iter = 10000,
#   n.chains = 1,
#   trace = TRUE, 
#   n.adapt = 5000,
#   thin = 1, adaptive = TRUE)
