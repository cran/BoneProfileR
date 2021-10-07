#' BP_FitMLRadialCompactness estimates likelihood of model of a bone section
#' @title Estimation of the likelihood of a bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The -Ln L
#' @param bone The bone image to be used
#' @param fitted.parameters Parameters of the model to be fitted
#' @param fixed.parameters Fixed parameters of the model
#' @param analysis Name or rank of analysis
#' @param twosteps Should a 2-steps analysis be performed?
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
#' It generally ensures that global minimum is found. On the other hand, it doubles the time to complete for each angle.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run
#' library(BoneProfileR)
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  # or 
#'  bone <- BP_OpenImage(ijtiff=TRUE)
#'  bone <- BP_DetectBackground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectForeground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic", cut.angle=30)
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic")
#'  plot(bone)
#'  plot(bone, type="observations")
#'  plot(bone, type="observations+model", analysis=1)
#'  fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  bone <- BP_FitMLCompactness(bone, 
#'                 fitted.parameters=c(fittedpar, K1=1.01, K2=1.01), 
#'                 fixed.parameters=NULL, analysis="flexit")
#'  bone <- BP_FitBayesianCompactness(bone, analysis="flexit")
#'  mcmc <- RM_get(bone, RMname = "flexit", value="mcmc")
#'  fittedpar <- as.parameters(mcmc)
#'  bone <- BP_FitMLCompactness(bone, 
#'                 fitted.parameters=fittedpar, 
#'                 fixed.parameters=NULL, analysis="flexit")
#'  compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", alloptim=TRUE), 
#'              Flexit=BP_GetFittedParameters(bone, analysis="flexit", alloptim=TRUE))
#'  out4p <- plot(bone, type="observations+model", analysis="logistic")
#'  out6p <- plot(bone, type="observations+model", analysis="flexit")
#'  # The twosteps fit is more acurate but is around 100 times slower
#'  bone <- BP_FitMLRadialCompactness(bone, analysis="logistic", twosteps=TRUE)
#'  bone <- BP_FitMLRadialCompactness(bone, analysis="logistic", twosteps=FALSE)
#'  plot(bone, type="observations", angle=0)
#'  plot(bone, type="model", analysis="logistic", angle=0)
#'  plot(bone, type="observations+model", angle=0)
#'  plot(bone, type="observations+model", angle=pi)
#'  plot(bone, type="radial", radial.variable=c("P", "S"), analysis="logistic")
#'  plot(bone, type="radial", radial.variable=c("P", "S", "Min", "Max"), analysis="logistic")
#'  plot(bone, type="radial", radial.variable=c("TRC"), analysis="logistic")
#'  # Test using the change of orientation using default.angle from BP_EstimateCompactness():
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="logistic_rotation_pi")
#'  # With a pi rotation, the top moves to the bottom and the left moves to the right
#'  bone <- BP_EstimateCompactness(bone, rotation.angle=pi, analysis="logistic_rotation_pi")
#'  bone <- BP_FitMLRadialCompactness(bone, analysis="logistic_rotation_pi")
#'  plot(bone, type="radial", radial.variable=c("P", "S"), analysis="logistic")
#'  plot(bone, type="radial", radial.variable=c("P", "S"), analysis="logistic_rotation_pi")
#'  BP_Report(bone=bone, 
#'            analysis=1,
#'            docx=NULL, 
#'            pdf=NULL, 
#'            xlsx=file.path(getwd(), "report.xlsx"), 
#'            author="Marc Girondot", 
#'            title=attributes(bone)$name)
#' }
#' @export


BP_FitMLRadialCompactness <- function(bone, fitted.parameters=NULL,
                                      fixed.parameters=NULL, analysis=1, silent=FALSE, twosteps=TRUE) {
  
  # fitted.parameters=c(P=0.5, S=0.05, Min=0.001, Max=0.999); fixed.parameters=c(K1=1, K2=1); analysis=NULL; silent=FALSE; twosteps=TRUE
  # fitted.parameters=NULL; fixed.parameters=NULL; analysis=1; silent=FALSE; twosteps=TRUE
  
  
  if (is.null(fitted.parameters)) {
    if (is.null(BP_GetFittedParameters(bone, analysis=analysis))) {
      fitted.parameters=c(P=0.5, S=0.05, Min=0.001, Max=0.999)
      fixed.parameters=c(K1=1, K2=1)
    } else {
      fitted.parameters <- BP_GetFittedParameters(bone, analysis=analysis)
      fixed.parameters <- BP_GetFittedParameters(bone, analysis=analysis, alloptim = TRUE)$fixed.parameters
      
    }
  }
  
  array.compactness <- RM_get(x=bone, RMname=analysis, valuename = "array.compactness")
  partial <- RM_get(x=bone, RMname=analysis, valuename = "partial")
  
  result.radial <- matrix(data = NA, ncol=length(fitted.parameters), nrow=dim(array.compactness)[1])
  colnames(result.radial) <- names(fitted.parameters)
  # LnL <- 0
  # anglefait <- NULL
  # radial.modeled.compactness <- NULL
  # observed.modeled.compactness <- NULL
  # observed.compactness <- NULL
  
  distance.center <- RM_get(x=bone, RMname=analysis, valuename = "compactness.synthesis")$distance.center
  
  outmcl <- universalmclapply(1:dim(array.compactness)[1], FUN=function(angle) {
    
    # for (angle in 1:dim(array.compactness)[1]) {
    data_nm <- array.compactness[angle, , "0"]
    data_m <- array.compactness[angle, , "1"]
    if ((!partial) | (any(data_nm+data_m != 0))) {
      # o <- optim(par=fitted.parameters, fn=BP_LnLCompactness, 
      #            fixed.parameters=fixed.parameters, 
      #            data_nm=data_nm, data_m=data_m, 
      #            distance.center = RM_get(x=bone, RMname=analysis, valuename = "compactness.synthesis")$distance.center, 
      #            method = "Nelder-Mead")
      lower_limit <- c(P=0, S=-2, Min=0, Max=0.2, K1=-10, k2=-10)
      upper_limit <- c(P=1, S=2, Min=0.8, Max=1, K1=10, k2=10)
      lower <- lower_limit[names(fitted.parameters)]
      upper <- upper_limit[names(fitted.parameters)]
      
      
      o <- optim(par=fitted.parameters, fn=BP_LnLCompactness, 
                 fixed.parameters=fixed.parameters, 
                 data_nm=data_nm, data_m=data_m, 
                 upper = upper, lower = lower, 
                 distance.center = distance.center, 
                 control=list(maxit=1000), 
                 method = "L-BFGS-B", hessian = FALSE)
      
      if (twosteps) {
        # lancement en mcmc
        p <- o$par
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
                                             Prior2=10, 
                                             SDProp=0.3, 
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
        
        mcmc <- HelpersMG::MHalgoGen(
          likelihood = BP_LnLCompactness,
          bone=NULL, 
          data_nm=data_nm, data_m=data_m, 
          distance.center = distance.center, 
          fixed.parameters=fixed.parameters, 
          parameters_name = "par", 
          parameters = priors, n.iter = 10000,
          n.chains = 1,
          n.adapt = 100,
          thin = 1, adaptive = TRUE)
        
        fitted.parameters <- HelpersMG::as.parameters(mcmc)
        
        o <- optim(par=fitted.parameters, fn=BP_LnLCompactness, 
                   fixed.parameters=fixed.parameters, 
                   data_nm=data_nm, data_m=data_m, 
                   upper = upper, lower = lower, 
                   distance.center = distance.center, 
                   control=list(maxit=1000), 
                   method = "L-BFGS-B", hessian = FALSE)
      }
      
      LnL <- o$value
      result.radial <- o$par
      anglefait <- angle
      
      p <- c(o$par, fixed.parameters)
      
      Min <- p["Min"]
      Max <- p["Max"]
      
      
      # 21/02/2020
      p["S"] <- 1/(4*p["S"])
      cp <- HelpersMG::flexit(x = distance.center, par = p)
      radial.modeled.compactness <- mean(cp * (Max - Min) + Min)
      observed.modeled.compactness <- mean(data_m/(data_m+data_nm), na.rm = TRUE)
      observed.compactness <- sum(data_m)/sum(data_m+data_nm)
      limitlow <- distance.center[which.min(abs(cp-0.025))[1]]
      limithigh <- distance.center[which.min(abs(cp-0.975))[1]]
      TRC <- abs(limithigh-limitlow)
      
    } else {
      result.radial <- NA
      radial.modeled.compactness <- NA
      observed.modeled.compactness <- NA
      observed.compactness <- NA
      anglefait <- NA
      LnL <- 0
      TRC <- NA
    }
    
    return(list(result.radial=result.radial, 
                radial.modeled.compactness=radial.modeled.compactness, 
                observed.modeled.compactness=observed.modeled.compactness, 
                observed.compactness=observed.compactness, 
                anglefait=anglefait, 
                TRC=TRC, 
                LnL=LnL))
  }, mc.cores = parallel::detectCores(), 
  clusterExport=list(varlist=c("array.compactness", "twosteps", "fixed.parameters", 
                               "fitted.parameters", "distance.center"), envir=environment()))
  
  result.radial <- sapply(X = outmcl, FUN = function(x) x["result.radial"])
  result.radial <- t(as.data.frame(result.radial))
  rownames(result.radial) <- NULL
  
  radial.modeled.compactness <- unname(unlist(sapply(X = outmcl, FUN = function(x) x["radial.modeled.compactness"])))
  observed.modeled.compactness <- unname(unlist(sapply(X = outmcl, FUN = function(x) x["observed.modeled.compactness"])))
  observed.compactness <- unname(unlist(sapply(X = outmcl, FUN = function(x) x["observed.compactness"])))
  anglefait <- na.omit(unname(unlist(sapply(X = outmcl, FUN = function(x) x["anglefait"]))))
  LnL <- sum(unname(unlist(sapply(X = outmcl, FUN = function(x) x["LnL"]))), na.rm = TRUE)
  TRC <- unname(unlist(sapply(X = outmcl, FUN = function(x) x["TRC"])))
  
  o <- list()
  o$par <- result.radial
  o$value <- LnL
  o$counts <- NULL
  o$convergence <- NULL
  o$message <- NULL
  o$fixed.parameters <- fixed.parameters
  o$angles <- (RM_get(x=bone, RMname=analysis, valuename = "cut.angle")[-1]+
                 rev(rev(RM_get(x=bone, RMname=analysis, valuename = "cut.angle"))[-1]))/2
  o$radial.modeled.compactness <- radial.modeled.compactness
  o$observed.modeled.compactness <- observed.modeled.compactness
  o$observed.compactness <- observed.compactness
  par <- result.radial
  
  if (!is.null(fixed.parameters)) {
    fpmat <- matrix(rep(fixed.parameters, nrow(par)), nrow=nrow(par), byrow = TRUE)
    colnames(fpmat) <- names(fixed.parameters)
    par <- cbind(par, fpmat)
  }
  
  
  # Min <- par[, "Min"]
  # Max <- par[, "Max"]
  
  # par[, "Min"] <- Min
  # par[, "Max"] <- Max
  
  tablestat <- data.frame(mean=numeric(ncol(par)), 
                          sd=numeric(ncol(par)), stringsAsFactors = FALSE)
  for (i in 1:ncol(par)) {
    tablestat[i, "mean"] <- mean(par[, i], na.rm = TRUE)
    tablestat[i, "sd"] <- sd(par[, i], na.rm = TRUE)
  }
  
  rownames(tablestat) <- colnames(par)
  
  
  o$summary.table <- tablestat
  
  # S <- par[ ,"S"]
  # S <- 1/(4*S)
  # 
  # K1 <- par[ ,"K1"]
  # K2 <- par[ ,"K2"]
  # P <- par[ , "P"]
  # 
  # l <- 0.025
  # K1 <- ifelse(K1==0, 1E-9, K1)
  # K2 <- ifelse(K2==0, 1E-9, K2)
  # 
  # K1 <- ifelse(is.infinite(2^(K1)), sign(K1)*500, K1)
  # K2 <- ifelse(is.infinite(2^(K2)), sign(K2)*500, K2)
  # 
  # S1 <- (2^(K1 - 1)*K1*S)/(2^(K1) - 1)
  # S2 <- (2^(K2 - 1)*K2*S)/(2^(K2) - 1)
  # 
  # limit.low.TRC <- P -  log(((1/(1-l)) ^ K1 - 1)/(2^K1 - 1))/(4 * S1)
  # limit.high.TRC <- 1/(4 * S2) * log(((1/(1-l))^K2 - 1)/(2^K2 - 1)) + P
  # 
  # TRC <- abs(unname(limit.high.TRC-limit.low.TRC))
  # 
  par <- cbind(par, TRC=TRC)
  
  o$synthesis <- cbind(par, angles=attributes(bone)$optimRadial[[analysis]]$angles)
  
  bone <- RM_add(x=bone, RMname = analysis, valuename = "optimRadial", value=o)
  
  if (!silent) print(tablestat)
  
  return(bone)
}


