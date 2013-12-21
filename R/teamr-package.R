#' TEAM Remote Sensing Processing Tools
#'
#' A set of scripts for processing remote sensing data to support the work of 
#' the Tropical Ecology Assessment and Monitoring (TEAM) Network. Intended to 
#' be used in conjunction with the set of Python scripts available in the 
#' \href{https://github.com/azvoleff/teampy}{teampy repository} on github.
#'
#' @name teamr-package
#' @aliases teamr
#' @docType package
#' @title TEAM Remote Sensing Processing Tools
#' @author Alex Zvoleff, \email{azvoleff@@conservation.org}
#' @keywords package
#' @useDynLib teamr
NULL
#' Landsat 5 Surface Reflectance Image from February 6, 1986 (path 15, row 53)
#' 
#' Portion of Landsat 5 Surface Reflectance image from the Landsat Climate Data 
#' Record archive. This subset of the image includes only bands 1-4.
#'
#' @docType data
#' @name L5TSR_1986
#' @seealso L5TSR_2001
NULL
#' Landsat 5 Surface Reflectance Image from January 14, 2001 (path 15, row 53)
#' 
#' Portion of Landsat 5 Surface Reflectance image from the Landsat Climate Data 
#' Record archive. This subset of the image includes only bands 1-4.
#'
#' @docType data
#' @name L5TSR_2001
#' @seealso L5TSR_1986
NULL
#' Subset of ASTER Digital Elevation Model V002
#' 
#' @docType data
#' @name ASTER_V002_WEST
#' @seealso ASTER_V002_EAST
NULL
#' Training polygons for 1986 and 2001 Landsat 5 Surface Reflectance images
#' 
#' Polygons digitized from 1986 and 2001 Landsat 5 Surface Reflectance image 
#' from the Landsat Climate Data Record archive. The training polygons can be 
#' used for testing classification algorithms.
#'
#' There are three columns in the dataset. "class_1986" is the cover class for 
#' the pixels in the polygon from the 1986 image. "class_2001" is the cover class 
#' for the pixels in the polygon from the 2001 image.
#'
#' @docType data
#' @name L5TSR_1986_2001_training
NULL
#' Subset of ASTER Digital Elevation Model V002
#' 
#' @docType data
#' @name ASTER_V002_WEST
#' @seealso ASTER_V002_EAST
NULL
#' Subset of ASTER Digital Elevation Model V002
#' 
#' @docType data
#' @name ASTER_V002_EAST
#' @seealso ASTER_V002_WEST
NULL
#' Randomly generated 100x100 test image
#'
#' Used in testing the output from the GLCM texture statistics C++ code.
#' 
#' @docType data
#' @name glcm_test_raster
#' @examples
#' # The image was generated with the following code:
#' set.seed(0)
#' test_matrix <- matrix(runif(100)*32, nrow=10)
#' glcm_test_raster <- raster(test_matrix, crs='+init=EPSG:4326')
#' glcm_test_raster <- cut(glcm_test_raster, seq(0, 32))
NULL
#' GLCM textures calculated in EXELIS ENVI
#'
#' This is the output from running a "co-occurrence measures" calculation to 
#' calculate GLCM textures in EXELIS ENVI from the \code{glcm_test_raster} 
#' included in the \code{teamr} package. The following settings were used:
#'     window size: 3x3
#' 	   co-occurrence shift: 1 (x), 1 (y)
#' 	   greyscale quantization levels: 32
#'     textures to compute: mean, variance, homogeneity, contrast, dissimilarity, 
#'         entropy, second moment, correlation
#' 
#' @docType data
#' @name glcm_test_raster_ENVI_textures
NULL
#' Example output from classify_image.R
#'
#' This is the output from running the classification example from the
#' \code{classify_image} function.
#' 
#' @docType data
#' @name classified_LT5SR_1986
NULL
