.onLoad <- function(libname, pkgname) {
    if(is.null(getOption("teamr.zoifile.suffix.geo")))
        options(teamr.zoifile.suffix.geo='_ZOI_GEO.shp')
    
    if(is.null(getOption("teamr.zoifile.suffix.geosimple")))
        options(teamr.zoifile.suffix.geo='_ZOI_GEO_simple.shp')
    
    if(is.null(getOption("teamr.zoifile.suffix.utm")))
        options(teamr.zoifile.suffix.geo='_ZOI_UTM.shp')
}
