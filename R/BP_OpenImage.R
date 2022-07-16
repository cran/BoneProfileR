#' BP_OpenImage opens an image
#' @title Open an image
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Characteristics of an image
#' @param file The file to be opened
#' @param name Name of this slice
#' @param ijtiff Should the ijtiff must be used to read tiff image
#' @description Open an image.
#' @family BoneProfileR
#' @examples
#' \dontrun{
#' # Not run:
#' library(BoneProfileR)
#' path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  plot(bone)
#'  path_Hedgehog <- system.file("extdata", "Erinaceus_europaeus_fem_2-1_small.tif", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Hedgehog)
#'  plot(bone)
#'  bone <- BP_OpenImage(file=path_Hedgehog, ijtiff=TRUE)
#'  plot(bone)
#'  # A partial section
#'  path_Dicynodon <- system.file("extdata", "Dicynodon_tibia_11.11.1.T_b_b-1.png", 
#'                              package = "BoneProfileR")
#'  bone <- BP_OpenImage(file=path_Dicynodon)
#'  plot(bone)
#'  # To open a file with a dialog:
#'  bone <- BP_OpenImage()
#' }
#' @export


BP_OpenImage <- function(file=file.choose(), name=NULL, ijtiff=FALSE) {
  
  if (grepl("\\.tif", file) | grepl("\\.TIF", file)) {
    
    if (isFALSE(ijtiff)) {
      
      if (!requireNamespace("tiff", quietly = TRUE)) {
        stop("tiff package is absent; Please install it first to read tiff image")
      }
 

      bone <- suppressWarnings(tiff::readTIFF(file))
      if (length(dim(bone))==3) {
        bone <- aperm(bone, c(2, 1, 3))
      } else {
        bone_pre <- array(data = NA, dim=c(dim(bone)[2], dim(bone)[1], 1, 1))
        bone_pre[, , 1, 1] <- bone
        bone <- bone_pre
      }
      bone <- suppressWarnings(as.cimg(bone))
    } else {
      
      if (!requireNamespace("ijtiff", quietly = TRUE)) {
        stop("ijtiff package is absent; Please install it first to read tiff image with option ijtiff")
      }
      
      bone <- ijtiff::read_tif(file)
      bone <- aperm(bone, c(2, 1, 4, 3))
      bone <- suppressWarnings(as.cimg(bone))
      if (max(bone[, , , 1]) > 1) bone[, , , 1] <- bone[, , , 1] / 255
      if (dim(bone)[4] > 1) if (max(bone[, , , 2]) > 1) bone[, , , 2] <- bone[, , , 2] / 255
      if (dim(bone)[4] > 2) if (max(bone[, , , 3]) > 1) bone[, , , 3] <- bone[, , , 3] / 255
    }
    
  } else {
    bone <- suppressWarnings(load.image(file))
  }
  
  if (dim(bone)[4] != 3) {
    
    bone_pre <- array(data = 0, dim=c(dim(bone)[1], dim(bone)[2], 1, 3))
    bone_pre <- suppressWarnings(as.cimg(bone_pre))
    
    if (dim(bone)[4] == 1) {
      # bone_pre[, , 1, 1] <- bone
      bone_pre <- add.colour(bone)
      # class(bone_pre) <- c("cimg", "imager_array", "numeric" )
    }
    if (dim(bone)[4] == 2) {
      bone_pre[, , 1, 1:dim(bone)[4]] <- bone
      # bone_pre[, , 1, 3] <- 1
      # class(bone_pre) <- c("cimg", "imager_array", "numeric" )
    }
    
    if (dim(bone)[4] == 4) {
      bone_pre <- flatten.alpha(bone)
      # bone_pre[, , 1, 3] <- 1
      # class(bone_pre) <- c("cimg", "imager_array", "numeric" )
    }
    
    
    bone <- bone_pre
  }
  
  
  if (is.null(name)) {
    name <- basename(file)
  }
  
  bone <- addS3Class(bone, "BoneProfileR")
  attributes(bone) <- modifyList(attributes(bone), list(name=name))
  
  return(bone)
}

