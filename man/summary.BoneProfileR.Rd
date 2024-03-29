% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BP_SummaryBone.R
\name{summary.BoneProfileR}
\alias{summary.BoneProfileR}
\title{Plot a bone section}
\usage{
\method{summary}{BoneProfileR}(object, max.level = FALSE, analysis = 1, ...)
}
\arguments{
\item{object}{The bone image}

\item{max.level}{If TRUE, will return all list element of the objects}

\item{analysis}{The analysis to report the global compacteness}

\item{...}{Not used}
}
\value{
An invisible list with recorded information
}
\description{
Display information of bone section
}
\details{
summary.BoneProfileR displays a bone section
}
\examples{
\dontrun{
# Not run:
library(BoneProfileR)
 bone <- BP_OpenImage()
 # or 
 path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
                             package = "BoneProfileR")
 bone <- BP_OpenImage(file=path_Hedgehog)
 summary(bone)
}
}
\seealso{
Other BoneProfileR: 
\code{\link{BP_AutoFit}()},
\code{\link{BP_ChooseBackground}()},
\code{\link{BP_ChooseCenter}()},
\code{\link{BP_ChooseForeground}()},
\code{\link{BP_DetectBackground}()},
\code{\link{BP_DetectCenters}()},
\code{\link{BP_DetectForeground}()},
\code{\link{BP_DuplicateAnalysis}()},
\code{\link{BP_EstimateCompactness}()},
\code{\link{BP_FitBayesianCompactness}()},
\code{\link{BP_FitMLCompactness}()},
\code{\link{BP_FitMLRadialCompactness}()},
\code{\link{BP_GetFittedParameters}()},
\code{\link{BP_ListAnalyses}()},
\code{\link{BP_LnLCompactness}()},
\code{\link{BP_OpenImage}()},
\code{\link{BP_Report}()},
\code{\link{Erinaceus_europaeus}},
\code{\link{plot.BoneProfileR}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{BoneProfileR}
