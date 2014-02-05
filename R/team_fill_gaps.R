#' Perform SLC-off gap fill of Landsat 7 ETM+ image
#'
#' Calls GNSPI.pro IDL script by Zhu Xiaolin to fill gaps in SLC-off Landsat 7 
#' ETM+ image.
#'
#' @param sample_size the sample size of sample pixels
#' @param size_wind the maximum window size
#' @param class_num the estimated number of classes
#' @param DN_min the minimum DN value of the image
#' @param DN_max the maximum DN value of the image
#' @param patch_long the size of block, to process whole ETM scene, set to 1000
#' @export
#' @references Zhu, X., Liu, D., Chen, J., 2012. A new geostatistical approach 
#' for filling gaps in Landsat ETM+ SLC-off images. Remote Sensing of 
#' Environment 124, 49-60.
team_fill_gaps <- function(sample_size=20, size_wind=12, class_num=4, 
                           DN_min=0.0, DN_max=1.0, patch_long=500) {

    #TODO: calculate num_series based on length of input image timeseries

}
