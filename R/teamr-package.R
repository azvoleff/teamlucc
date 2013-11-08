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
#' Training polygons for 1986 and 2001 Landsat 5 Surface Reflectance images
#' 
#' Polygons digitized from 1986 and 2001 Landsat 5 Surface Reflectance image 
#' from the Landsat Climate Data Record archive. The training polygons can be 
#' used for testing classification algorithms.
#'
#' There are three columns in the dataset. "t1_class" is the cover class for 
#' the pixels in the polygon from the 1986 image. "t2_class" is the cover class 
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
