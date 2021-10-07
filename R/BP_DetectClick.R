# #' BP_DetectClick waits for a click in the image
# #' @title Wait for a click in the image
# #' @author marc.girondot@@u-psud.fr
# #' @return a vector with x and y
# #' @param bone The bone image
# #' @description Wait for a click in the image.
# #' @family BoneProfileR
# #' @examples
# #' \dontrun{
# #' bone <- BP_OpenImage()
# #' BP_PlotBone(bone)
# #' pos <- BP_DetectClick(bone)
# #' }
# #' @export


.BP_DetectClick <- function(bone) {
  
  
  
  click.locn <- locator(n = 1, type = "n")
  
  # click.locn$x <- click.locn$x * dim(bone)[1]
  # click.locn$y <- click.locn$y * dim(bone)[2]
  
  if (is.null(click.locn)) {
    click.locn <- list(x=NA, y=NA)
  } else {
    
    if (click.locn$x<1) click.locn$x <- 1
    if (click.locn$x>dim(bone)[1]) click.locn$x <- dim(bone)[1]
    if (click.locn$y<1) click.locn$y <- 1
    if (click.locn$y>dim(bone)[2]) click.locn$y <- dim(bone)[2]
  }
  
  return(c(x=click.locn$x, y=click.locn$y))
}
