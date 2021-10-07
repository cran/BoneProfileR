#' BP runs a shiny application to fit bone section
#' @title Run a shiny application to fit bone section
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Nothing
#' @description Run a shiny application to fit bone section
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#' BP()
#' }
#' @export


BP <- function() {
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is absent; Please install it first")
  }
  
  if(interactive()){
    getFromNamespace("runApp", ns="shiny")(appDir = system.file("shiny", package="BoneProfileR"), 
                                         launch.browser =TRUE)
  }
  
}
