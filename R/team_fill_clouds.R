#' Perform heavy cloud filling
#'
#' Calls CLOUD_REMOVE.pro IDL script by Zhu Xiaolin to fill heavy clouds in 
#' Landsat image.
#'
#' @export
#' @importFrom tools file_path_sans_ext
#' @param img_cloudy the clear image (base image)
#' @param img_clear the cloudy image
#' @param img_cloud_mask cloud mask, with clouded areas set to
#' @param num_class set the estimated number of classes in image
#' @param extent1 set the range of cloud neighborhood
#' @param DN_min minimum of DN
#' @param DN_max maximum of DN
#' @param patch_long block size when process large image
#' @references Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified 
#' neighborhood similar pixel interpolator approach for removing thick clouds 
#' in Landsat images. Geoscience and Remote Sensing Letters, IEEE 9, 521-525.
team_fill_clouds <- function(img_cloudy, img_clear, img_cloud_mask, 
                             script_path='CLOUD_REMOVE.sav', num_class=1, 
                             extent1=1, DN_min=0, DN_max=255,
                             idl="C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe") {
   idl <- "C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe"
   script_path <- 'C:/Users/azvoleff/Code/TEAM/team_IDL_code/cloud_removal/CLOUD_REMOVE.pro'

   #idl_cmd <- paste('CD,', script_dir, ';\n.run CLOUD_REMOVE;\nexit\n;', sep='')

   out <- system(paste(shQuote(idl), shQuote(script_path)), invisible=FALSE, intern=TRUE)
   print(out)

}
