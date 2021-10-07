#' BP_Report save a pdf report for the analyzed bone
#' @title Generate a pdf report for the analyzed bone
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Nothing
#' @param bone The bone image
#' @param control.plot A list with the parameters used for plot
#' @param analysis Indicate analysis name or rank that you want report
#' @param pdf Name of pdf file
#' @param docx Name of Word file
#' @param xlsx Name of Excel file
#' @param author Name indicated in the report
#' @param title Title of the report
#' @description Generate a docx, xlsx, or pdf report.\cr
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  bone <- BP_DetectBackground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectForeground(bone=bone, analysis="logistic")
#'  bone <- BP_DetectCenters(bone=bone, analysis="logistic")
#'  bone <- BP_EstimateCompactness(bone, analysis="logistic")
#'  bone <- BP_FitMLCompactness(bone, analysis="logistic")
#'  fittedpar <- BP_GetFittedParameters(bone, analysis="logistic")
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="flexit")
#'  bone <- BP_FitMLCompactness(bone, 
#'                 fitted.parameters=c(fittedpar, K1=1, K2=1), 
#'                 fixed.parameters=NULL, analysis="flexit")
#'  compare_AIC(Logistic=BP_GetFittedParameters(bone, analysis="logistic", alloptim=TRUE), 
#'              Flexit=BP_GetFittedParameters(bone, analysis="flexit", alloptim=TRUE))
#'  bone <- BP_FitMLRadialCompactness(bone, analysis="logistic")
#'  # Test using the change of orientation using default.angle from BP_EstimateCompactness():
#'  bone <- BP_DuplicateAnalysis(bone, from="logistic", to="logistic_rotation_pi")
#'  # With a pi rotation, the top moves to the bottom and the left moves to the right
#'  bone <- BP_EstimateCompactness(bone, rotation.angle=pi, analysis="logistic_rotation_pi")
#'  bone <- BP_FitMLRadialCompactness(bone, analysis="logistic_rotation_pi")
#'  BP_Report(bone=bone, 
#'            analysis=1,
#'            docx=NULL, 
#'            pdf=NULL, 
#'            xlsx=file.path(getwd(), "report.xlsx"), 
#'            author="Marc Girondot", 
#'            title=attributes(bone)$name)
#'            
#'  BP_Report(bone=bone, 
#'            analysis=1,
#'            docx=NULL, 
#'            pdf=file.path(getwd(), "report.pdf"), 
#'            xlsx=NULL, 
#'            author="Marc Girondot", 
#'            title=attributes(bone)$name)
#'            
#'  BP_Report(bone=bone, 
#'            analysis=1,
#'            docx=file.path(getwd(), "report.docx"), 
#'            pdf=NULL, 
#'            xlsx=NULL, 
#'            author="Marc Girondot", 
#'            title=attributes(bone)$name)
#' }
#' @export


BP_Report <- function(bone=stop("A bone section must be provided"), 
                      control.plot=list(message=NULL, show.centers=TRUE, 
                                        show.colors=TRUE, show.grid=TRUE, 
                                        CI="ML", show.legend=TRUE), 
                      analysis=1, 
                      docx=file.path(getwd(), "report.docx"), 
                      pdf=file.path(getwd(), "report.pdf"), 
                      xlsx=file.path(getwd(), "report.xlsx"), 
                      author=NULL, title=attributes(bone)$name) {
  
  # control.plot=list(message=NULL, show.centers=TRUE,  show.colors=TRUE, show.grid=TRUE, CI="ML", show.legend=TRUE); analysis=1; pdf=NULL; docx=file.path(getwd(), "report.docx"); xlsx=file.path(getwd(), "report.xlsx"); author=NULL; title=attributes(bone)$name
  
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
    stop("The analysis does no exist.")
  }
  
  texte <- NULL
  
  date <- out[[analysis]]$timestamp
  
  if (!is.null(pdf)) {
    control.output.pdf=list(Title=title, 
                            Author=author, 
                            Date=date, 
                            filename=pdf, 
                            dirname=NULL)
    
    texte.pdf <- c('---', paste0('title: "', control.output.pdf$Title, '"'), 
                   paste0('author: "', control.output.pdf$Author, '"'), 
                   paste0('date: "', control.output.pdf$Date, '"'), 
                   'output: "pdf_document"', 
                   '---'      )
    
    texte.pdf <- c(texte.pdf, "")
    texte.pdf <- c(texte.pdf, paste0("# ", analysis))
    texte.pdf <- c(texte.pdf, "")
    
    texte.pdf <- c(texte.pdf,"```{r echo=FALSE}" , 
                   paste0("control.plot <- ", paste0(capture.output(dput(control.plot)), collapse = "")), 
                   "do.call(getFromNamespace('plot.BoneProfileR', ns='BoneProfileR'), 
                          modifyList(control.plot, list(x=env$bone)))",
                   "```")
    
    texte.pdf <- c(texte.pdf, "")
    
    out1 <- RM_get(x=bone, RMname=analysis, valuename = "optim")
    
    if (!is.null(out1)) {
      
      texte.pdf <- c(texte.pdf, "## Global model of compactness")
      
      texte.pdf <- c(texte.pdf, 
                     "```{r echo=FALSE}" , 
                     paste0("do.call(getFromNamespace('plot.BoneProfileR', ns='BoneProfileR'), 
                          modifyList(control.plot, list(x=env$bone, type='observations+model', analysis='", analysis,"')))"),
                     "```")
      texte.pdf <- c(texte.pdf, "")
      
      texte.pdf <- c(texte.pdf, knitr::kable(out1$summary.table)) 
      texte.pdf <- c(texte.pdf, "")
    }
    
    out1 <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")
    
    if (!is.null(out1)) {
      
      texte.pdf <- c(texte.pdf, "## Radial model of compactness")
      
      texte.pdf <- c(texte.pdf, 
                     "```{r echo=FALSE}" , 
                     paste0("do.call(getFromNamespace('plot.BoneProfileR', ns='BoneProfileR'), 
                          modifyList(control.plot, list(x=env$bone, type='radial', analysis='", analysis,"')))"),
                     "```")
      texte.pdf <- c(texte.pdf, "")
      
      texte.pdf <- c(texte.pdf, knitr::kable(out1$summary.table)) 
      texte.pdf <- c(texte.pdf, "")
    }
    
    
    texte.pdf <- c(texte.pdf, "", 
                   "### This software is provided by [Marc Girondot](https://max2.ese.u-psud.fr/epc/conservation/Marc.html), Ecologie, Syst\u00E9matique, Evolution, CNRS, Universit\u00E9 Paris Saclay, AgroParisTech.")
    
    texte.pdf <- iconv(texte.pdf, from = "", to="UTF-8")
    tmp <- tempdir()
    Rmd <- file.path(tmp, "temppdf.Rmd")
    writeLines(texte.pdf, con=Rmd)

    env <- new.env()
    assign("bone", bone, envir = env)

    rmarkdown::render(input=Rmd, 
                      output_file = pdf, 
                      output_dir=NULL, 
                      envir = env
    )
    
  }
  if (!is.null(docx)) {
    control.output.docx=list(Title=title, 
                             Author=author, 
                             Date=date, 
                             filename=docx, 
                             dirname=NULL)
    
    texte.docx <- c('---', paste0('title: "', control.output.docx$Title, '"'), 
                    paste0('author: "', control.output.docx$Author, '"'), 
                    paste0('date: "', control.output.docx$Date, '"'), 
                    'output: "word_document"', 
                    #            ':reference_docx: template.docx', 
                    '---'      )
    
    texte.docx <- c(texte.docx, "")
    texte.docx <- c(texte.docx, paste0("# ", analysis))
    texte.docx <- c(texte.docx, "")
    
    texte.docx <- c(texte.docx,"```{r echo=FALSE}" , 
                    paste0("control.plot <- ", paste0(capture.output(dput(control.plot)), collapse = "")), 
                    "do.call(getFromNamespace('plot.BoneProfileR', ns='BoneProfileR'), 
                          modifyList(control.plot, list(x=bone)))",
                    "```")
    
    texte.docx <- c(texte.docx, "")
    out1 <- RM_get(x=bone, RMname=analysis, valuename = "optim")
    
    if (!is.null(out1)) {
      
      texte.docx <- c(texte.docx, "## Global model of compactness")
      
      texte.docx <- c(texte.docx, 
                      "```{r echo=FALSE}" , 
                      paste0("do.call(getFromNamespace('plot.BoneProfileR', ns='BoneProfileR'), 
                          modifyList(control.plot, list(x=bone, type='observations+model', analysis='", analysis,"')))"),
                      "```")
      texte.docx <- c(texte.docx, "")
      
      texte.docx <- c(texte.docx, knitr::kable(out1$summary.table)) 
      texte.docx <- c(texte.docx, "")
    }
    
    out1 <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")
    
    if (!is.null(out1)) {
      
      texte.docx <- c(texte.docx, "## Radial model of compactness")
      
      texte.docx <- c(texte.docx, 
                      "```{r echo=FALSE}" , 
                      paste0("do.call(getFromNamespace('plot.BoneProfileR', ns='BoneProfileR'), 
                          modifyList(control.plot, list(x=bone, type='radial', analysis='", analysis,"')))"),
                      "```")
      texte.docx <- c(texte.docx, "")
      
      texte.docx <- c(texte.docx, knitr::kable(out1$summary.table)) 
      texte.docx <- c(texte.docx, "")
    }
    
    texte.docx <- c(texte.docx, "", 
                    "### This software is provided by [Marc Girondot](https://max2.ese.u-psud.fr/epc/conservation/Marc.html), Ecologie, Syst\u00E9matique, Evolution, CNRS, Universit\u00E9 Paris Saclay, AgroParisTech.")
    
    texte.docx <- iconv(texte.docx, from = "", to="UTF-8")
 
    tmp <- tempdir()
    Rmd <- file.path(tmp, "tempdocx.Rmd")
    writeLines(texte.docx, con=Rmd)
    
    env <- new.env()
    assign("bone", bone, envir = env)
    
    rmarkdown::render(input=Rmd, 
                      output_file = docx, 
                      output_dir=NULL, 
                      envir = env
    )
  }
  
  
  if (!is.null(xlsx)) {
    
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("openxlsx package is absent; Please install it first to export Excel file")
    }
    

    wb <- openxlsx::createWorkbook(creator = author
                                   , title = title
                                   , subject = "BoneProfileR report"
                                   , category = "")
    
    openxlsx::addWorksheet(
      wb=wb,
      sheetName="Global")
    
    openxlsx::addWorksheet(
      wb=wb,
      sheetName="Radial")
    
    out1 <- RM_get(x=bone, RMname=analysis, valuename = "optim")
    
    if (!is.null(out1)) {
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x="Global model of compactness",
        startCol = 1,
        startRow = 1)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=date,
        startCol = 1,
        startRow = 2)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=analysis,
        startCol = 1,
        startRow = 3)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x="Observed compactness",
        startCol = 1,
        startRow = 4)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=RM_get(x=bone, RMname = analysis, valuename = "global.compactness"),
        startCol = 2,
        startRow = 4)
   
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x="Modeled compactness by ML (2.5%, 50%, 97.5%)",
        startCol = 1,
        startRow = 5)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=mean(RM_get(x=bone, RMname = analysis, valuename = "optim")$quantiles["2.5%", ]),
        startCol = 2,
        startRow = 5)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=mean(RM_get(x=bone, RMname = analysis, valuename = "optim")$quantiles["50%", ]),
        startCol = 3,
        startRow = 5)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=mean(RM_get(x=bone, RMname = analysis, valuename = "optim")$quantiles["97.5%", ]),
        startCol = 4,
        startRow = 5)
      
      if (!is.null(RM_get(x=bone, RMname = analysis, valuename = "mcmc"))) {
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x="Modeled compactness by MCMC (2.5%, 50%, 97.5%)",
        startCol = 1,
        startRow = 6)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=mean(RM_get(x=bone, RMname = analysis, valuename = "mcmc")$quantiles["2.5%", ]),
        startCol = 2,
        startRow = 6)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=mean(RM_get(x=bone, RMname = analysis, valuename = "mcmc")$quantiles["50%", ]),
        startCol = 3,
        startRow = 6)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=mean(RM_get(x=bone, RMname = analysis, valuename = "mcmc")$quantiles["97.5%", ]),
        startCol = 4,
        startRow = 6)
      }
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=out1$summary.table,
        startCol = 2,
        startRow = 10)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=rownames(out1$summary.table),
        startCol = 1,
        startRow = 11)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x="Quantiles",
        startCol = 1,
        startRow = 20)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Global",
        x=t(out1$quantiles),
        startCol = 1,
        startRow = 21)
    }
    
    out1 <- RM_get(x=bone, RMname=analysis, valuename = "optimRadial")
    
    if (!is.null(out1)) {
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x="Radial model of compactness",
        startCol = 1,
        startRow = 1)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=date,
        startCol = 1,
        startRow = 2)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=analysis,
        startCol = 1,
        startRow = 3)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=out1$summary.table,
        startCol = 2,
        startRow = 4)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=rownames(out1$summary.table),
        startCol = 1,
        startRow = 5)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=out1$synthesis,
        startCol = 2,
        startRow = 13)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x="Angle",
        startCol = 1,
        startRow = 13)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=out1$angles,
        startCol = 1,
        startRow = 14)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x="Radial modeled compactness",
        startCol = 10,
        startRow = 13)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=out1$radial.modeled.compactness,
        startCol = 10,
        startRow = 14)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x="Radial observed modeled compactness",
        startCol = 11,
        startRow = 13)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=out1$observed.modeled.compactness,
        startCol = 11,
        startRow = 14)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x="Radial observed compactness",
        startCol = 12,
        startRow = 13)
      
      openxlsx::writeData(
        wb=wb,
        sheet="Radial",
        x=out1$observed.compactness,
        startCol = 12,
        startRow = 14)
    }
    
    openxlsx::saveWorkbook(wb, file = xlsx, overwrite = TRUE)
  }
}


