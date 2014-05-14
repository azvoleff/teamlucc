# teamlucc

[![Build Status](https://travis-ci.org/azvoleff/teamlucc.png)](https://travis-ci.org/azvoleff/teamlucc)

## Overview

The `teamlucc` package is designed to facilitate analysis of land use and cover 
change (LUCC) around the monitoring sites of the Tropical Ecology Assessment 
and Monitoring (TEAM) Network. The [TEAM Network](http://www.teamnetwork.org/) 
is a global network of sites in tropical forests focused on collecting 
standardized real-time data measuring the response of tropical forests to 
changing climate, land cover and land use, and population.

`teamlucc` assists with processing and analysis of TEAM data and remote sensing 
imagery of TEAM sites for measuring change in ecosystems from the 
plot-landscape scales. `teamlucc` supports a range of analyses using satellite 
imagery with TEAM monitoring data, including:

* Image selection from USGS archive
    * Parsing metadata files from USGS EarthExplorer
    * Plotting available imagery for an area of interest (AOI), including 
      AOIs that cover more than one path/row
    * Formatting an image order for upload to ESPA system
    * Downloading images from a USGS ESPA order

* Preprocessing
    * Extraction and file conversion of surface reflectance imagery from the 
      Landsat Climate Data Record (CDR) archive
    * Topographic correction using parallel processing (Goslee, 2011)
    * Cloud fill and gap fill (for SLC-off Landsat 7 scenes) with R scripts 
      wrapping the modified Neighborhood Similar Pixel Interpolator (NSPI) and
      Geostatistical Neighborhood Similar Pixel Interpolator (GNSPI) by Zhu et 
      al. (2012a, 2012b)
    * Calculation of vegetation indices and image texture measures from 
      grey-level co-occurrence matrices (GLCMs)

* Image classification using support vector machine (SVM)

* Change detection: implementation of Change Vector Analysis in Posterior 
  Probability Space (CVAPS) and Double Window Flexible Pace Search (DFPS) 
  algorithms (Chen et al. 2011)

* Accuracy assessment using quantity agreement and disagreement (Pontius and 
  Millones, 2011)

The toolkit is still under active development. Follow the [TEAM 
website](http://www.teamnetwork.org/) for news, and the [toolkit project page
on github](https://github.com/azvoleff/teamlucc) for the latest version of the 
code.


## Package installation

### Installing `teamlucc`

As `teamlucc` is still under development, it is not yet listed on 
[CRAN](http://cran.r-project.org).  The easiest way to install the (beta 
version) of the `teamlucc` package is using the 
[`devtools`](http://cran.r-project.org/web/packages/devtools/index.html) 
package by Hadley Wickham.

To install `devtools` type:

```R
install.packages('devtools')
```

at the R command prompt. This will fecth the latest version of `devtools` from 
CRAN. After installing `devtools` type:

```R
install_github('teamlucc', username='azvoleff')
```

at the R prompt to install the latest version of `teamlucc`. Typing the above 
command will also work if you already have `teamlucc` installed and want to 
install an updated version of the package.

*NOTE:* If you are installing on Windows, you will first need to install the 
appropriate version of [Rtools](http://cran.r-project.org/bin/windows/Rtools/) 
for your version of R (as `teamlucc` contains C++ code).

### (optional) Install IDL
[IDL](http://www.exelisvis.com/ProductsServices/IDL.aspx) is required for 
running the cloud fill and Landsat 7 SLC-off gap fill routines in `teamlucc`. If 
you do not have an IDL license, you can also use the (free) IDL virtual 
machine (VM). See [this 
page](http://www.exelisvis.com/Support/HelpArticlesDetail/TabId/219/ArtMID/900/ArticleID/12395/The-IDL-Virtual-Machine.aspx) 
for details on the IDL VM.

### (optional) Install GDAL

This step is only required if you want to use the `unstack_ledapscdr` 
function in `teamlucc`. All the other functions in `teamlucc` will work without 
installing GDAL.

#### Windows:

Download the [32bit](http://download.osgeo.org/osgeo4w/osgeo4w-setup-x86.exe) 
or [
64bit](http://download.osgeo.org/osgeo4w/osgeo4w-setup-x86_64.exe) [OSGeo4W](http://trac.osgeo.org/osgeo4w/) installer.

Run the installer. Choose the "Express Desktop Install".  On the "Select 
Packages" screen, ensure the GDAL screen package is checked. You can uncheck 
the boxes for QGIS and GRASS GIS if you don't want them installed (though I 
highly recommend QGIS).

[Edit your environment variables](http://support.microsoft.com/kb/310519):

1. Add "C:\OSGeo4W\bin" (or "C:\OSGeo4W64\bin" if you installed the 64bit 
version) to the "PATH" environment variable.
2. Add a new "GDAL_DATA" environment variable equal to "C:\OSGeo4W\share\gdal" 
(or "C:\OSGeo4W64\share\gdal" for the 64bit version).

#### Linux (ubuntu):

At a shell prompt, type:

``` sh
sudo apt-get install gdal-bin libgdal-dev
```

## Using teamlucc

For more information on using `teamlucc`, see the online help in R, and the 
[`teamlucc` webpage](http://www.azvoleff.com/teamlucc). The webpage includes 
examples of a number of specific applications of `teamlucc`, including:

* [Filtering and downloading Landsat 
  scenes](/articles/filtering-landsat-with-teamlucc)

* [Preprocessing imagery and 
  DEMS](/articles/preprocessing-imagery-with-teamlucc)

* [Cloud removal](/articles/cloud-removal-with-teamlucc)

* [Image classification](/articles/image-classification-with-teamlucc)

## Author Contact Information

[Alex Zvoleff](mailto:azvoleff@conservation.org)  
Postdoctoral Associate  
Tropical Ecology Assessment and Monitoring (TEAM) Network  
Conservation International  
2011 Crystal Dr. Suite 500  
Arlington, VA 22202  
USA

## References
Chen, J., Chen, X., Cui, X., Chen, J., 2011. Change vector analysis in 
posterior probability space: a new method for land cover change detection. IEEE 
Geoscience and Remote Sensing Letters 8, 317--321.

Goslee, S.C., 2011. Analyzing remote sensing data in R: the landsat package. 
Journal of Statistical Software 43, 1--25.

Pontius, R.G., Millones, M., 2011. Death to Kappa: birth of quantity 
disagreement and allocation disagreement for accuracy assessment. International 
Journal of Remote Sensing 32, 4407--4429.

Zhu, X., Gao, F., Liu, D., Chen, J., 2012a. A modified neighborhood similar 
pixel interpolator approach for removing thick clouds in Landsat images. 
Geoscience and Remote Sensing Letters, IEEE 9, 521--525.

Zhu, X., Liu, D., Chen, J., 2012b. A new geostatistical approach for filling 
gaps in Landsat ETM+ SLC-off images. Remote Sensing of Environment 124, 49--60.
