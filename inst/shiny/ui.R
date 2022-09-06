# library(shiny); runApp("/Users/marcgirondot/Documents/Espace_de_travail_R/_shiny/BoneProfileR")
# library(shiny); runApp("http://134.158.74.46/BoneProfileR/")


library(shiny)
package.BoneProfileR <- require('BoneProfileR')
version <- "2.4 build 766"

# 
# mycss <- "
# #plot-container {
# position: relative;
# }
# #loading-spinner {
# position: absolute;
# left: 50%;
# top: 50%;
# z-index: -1;
# margin-top: -33px;  /* half of the spinner's height */
# margin-left: -33px; /* half of the spinner's width */
# }
# #plot.recalculating {
# z-index: -2;
# }
# "

splitLayout <- function (..., cellWidths = NULL, cellArgs = list()) 
{
  children <- list2(...)
  childIdx <- !nzchar(names(children) %||% character(length(children)))
  attribs <- children[!childIdx]
  children <- children[childIdx]
  count <- length(children)
  if (length(cellWidths) == 0 || any(is.na(cellWidths))) {
    cellWidths <- sprintf("%.3f%%", 100/count)
  }
  cellWidths <- rep(cellWidths, length.out = count)
  cellWidths <- sapply(cellWidths, validateCssUnit)
  do.call(tags$div, c(list(class = "shiny-split-layout"), attribs, 
                      mapply(children, cellWidths, FUN = function(x, w) {
                        do.call(tags$div, c(list(style = sprintf("width: %s;", 
                                                                 w)), cellArgs, list(x)))
                      }, SIMPLIFY = FALSE)))
}

# Define UI for application that draws a histogram
fluidPage(
  titlePanel(h1("Bone Profile",
             img(src="Rlogo.png", height=40, width=40), align = "center"), 
             windowTitle = "Bone ProfileR"), # h1("Bone Profile", align = "center"), HTML('<h1><center>Bone Profile<img scr="Rlogo.png"></center></h1>')), # h1("Bone Profile", align = "center")), 
  p(HTML("<b><a href=\"https://max2.ese.u-psud.fr/epc/conservation/index.html\">Marc Girondot</a></b> - Laboratoire Ecologie, Systématique, Evolution"), align = "center"),
  p(HTML("Université Paris-Saclay, CNRS, AgroParisTech, France."), align = "center"), 
  
  wellPanel(
 
  p(HTML("<strong>BoneProfileR is a scientific method and a software used 
  to model bone section for paleontological and ecological studies.</strong>"), align = "left"),
  p(paste0("This web server version v. ",version, " is a simplified version of the complete tools available as an R package.")), 
  
  p(HTML("Open a bone section image, choose the options and 
         click 'Run the analysis' button."), align = "left")
  ), 

  # Show a plot of the generated distribution
  wellPanel(
    splitLayout(
      fileInput(inputId="FileOpen", label="Choose a local file with an image of a bone section", 
                multiple = FALSE, accept = c(".tif", ".png", ".jpg")), 
      checkboxInput(inputId="ijtiff", label = "Use IJTiff package to import image?", value = FALSE)
      , cellWidths=c("70%", "30%")
    ), 
    splitLayout(
      selectInput("center", label="Choose the center to be used "
                  , choices=list("Ontogenic center"="ontogenic", 
                                 "Section center"="section", 
                                 "Mineralized center"="mineralized", 
                                 "Unmineralized center"= "unmineralized")
                  , selected = "ontogenic", multiple = FALSE,
                  selectize = FALSE, size = NULL),
      textInput(inputId="rotation", label="Rotation angle", value="0")
      , cellWidths=c("50%", "50%")
    ), 
    splitLayout(
      textInput(inputId="angles", label="Number of angles", value="60")
      , textInput(inputId="distances", label="Number of ribbons from center", value="100")
      , cellWidths=c("50%", "50%")
    ), 
    splitLayout(
      checkboxInput(inputId="twosteps", label = "Use a two-steps fit? The time requires for a two-steps fit is around 5 minutes.", value = TRUE)
      , cellWidths=c("100%")
    ), 
    
    p(HTML("<b>Choose the variables to plot for radial analysis</b>"), align = "left"), 
    splitLayout(
      checkboxInput(inputId="RadialVarP", label="P", value=TRUE)
      , checkboxInput(inputId="RadialVarS", label="S", value=FALSE)
      , checkboxInput(inputId="RadialVarMin", label="Min", value=FALSE)
      , checkboxInput(inputId="RadialVarMax", label="Max", value=FALSE)
      , checkboxInput(inputId="RadialVarK1", label="K1", value=FALSE)
      , checkboxInput(inputId="RadialVarK2", label="K2", value=FALSE)
      , checkboxInput(inputId="RadialVarTRC", label="Transitional Range of Compacity", value=TRUE)
      , cellWidths=c("10%", "10%", "10%", "10%", "10%", "10%", "40%")
    ), 
    
    p("The Transitional Range of Compacity is the range of distances from the center where the compacity is between 2.5% to 97.5% of the compacity between Min and Max.")
    , 
    actionButton(inputId="goButton", label="Run the analysis", width="30%", 
                 icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    p(HTML("Be patient, to be completed, the analysis requires between 40 seconds or more than 5 minutes if you use 2-steps analysis."), align = "left")
  ), 
  p(""), 
  wellPanel(

    
    htmlOutput(outputId="TitleOut1"), 
    plotOutput("Plot"), 
    # tags$head(tags$style(HTML(mycss)))
    # , div(id = "plot-container",
    #       tags$img(src = "spinner.gif",
    #                id = "loading-spinner"),
    #       plotOutput("Plot")
    # ), 
    htmlOutput(outputId="TitleOut2"), 
    tableOutput(outputId="DataOut"), 
    htmlOutput(outputId="ResultOut1"), 
    htmlOutput(outputId="TitleOut3"), 
    htmlOutput(outputId="ResultOut2"), 
    tableOutput(outputId="GlobalOut"), 
    plotOutput("PlotModel"), 
    htmlOutput(outputId="TitleOut4"), 
    tableOutput(outputId="RadialOut1"), 
    tableOutput(outputId="RadialOut2"), 
    plotOutput("PlotModelRadial"),
    downloadButton("ExcelButton", "Export formated data for Excel"), 
  ), 
  HTML("<small><i><font color='#006699'>The Virtual Data initiative, run by LABEX P2IO and supported by Université Paris-Saclay, is thanked for providing computing resources on its cloud infrastructure.</font></i></small>")

)
