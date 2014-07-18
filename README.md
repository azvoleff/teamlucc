# teamlucc

[![Build Status](https://travis-ci.org/azvoleff/teamlucc.png)](https://travis-ci.org/azvoleff/teamlucc)

## Overview

The `teamlucc` package is designed to facilitate analysis of land use and cover 
change (LUCC) around the monitoring sites of the Tropical Ecology Assessment 
and Monitoring (TEAM) Network. The [TEAM Network](http://www.teamnetwork.org/) 
is a global network of sites in tropical forests wth standardized real-time 
data collection designed to measure tropical forest responses to climate 
variability and change, land cover and land use change, and other threats.

`teamlucc` assists with processing and analysis of remote sensing imagery. 
`teamlucc` supports a range of preprocessing steps and analyses, including: 

* Image selection from USGS archive
    * Parsing metadata files from USGS EarthExplorer
    * Plotting available imagery for an area of interest (AOI), including 
      AOIs that cover more than one path/row
    * Formatting an image order for upload to ESPA system
    * ~~Downloading images from a USGS ESPA order~~ Not working as of 7/1/2014 
      due to changes in the ESPA system.

* Preprocessing
    * Extraction and file conversion of surface reflectance imagery from the 
      Landsat Climate Data Record (CDR) archive
    * Topographic correction using parallel processing (Goslee, 2011)
    * Cloud fill and gap fill (for SLC-off Landsat 7 scenes), including support 
      for the modified Neighborhood Similar Pixel Interpolator (NSPI) and 
      Geostatistical Neighborhood Similar Pixel Interpolator (GNSPI) by Zhu et 
      al.  (2012a, 2012b)
    * Image normalization
    
* Calculation of vegetation indices and image texture measures from grey-level 
  co-occurrence matrices (GLCMs)

* Image classification using random forests

* Change detection using the Change Vector Analysis in Posterior Probability 
  Space (CVAPS) and Double Window Flexible Pace Search (DFPS) algorithms (Chen 
  et al. 2011)

* Accuracy assessment using user's, producer's and overall accuracies, in 
  addition to quantity agreement and disagreement (Pontius and Millones, 2011)

The toolkit is under active development. Follow the [TEAM 
website](http://www.teamnetwork.org/) for news, and the [toolkit project page
on github](https://github.com/azvoleff/teamlucc) for the latest updates.

## Package installation

### Installing `teamlucc`

**NOTE: If you are installing on Windows, you will need to install the 
appropriate version of [Rtools](http://cran.r-project.org/bin/windows/Rtools/) 
for your version of R (as `teamlucc` contains C++ code) before you follow the 
below steps.**

As `teamlucc` is still under development, it is not yet listed on 
[CRAN](http://cran.r-project.org).  The easiest way to install the `teamlucc` 
package is using the 
[`devtools`](http://cran.r-project.org/web/packages/devtools/index.html) 
package by Hadley Wickham.

To install `devtools` type:

```R
install.packages('devtools')
```

at the R command prompt. This will fecth the latest version of `devtools` from 
CRAN. After installing `devtools` type:

```R
library(devtools)
install_github('azvoleff/teamlucc')
```

at the R prompt to install the latest version of `teamlucc`. Typing the above 
command will also work if you already have `teamlucc` installed and want to 
install an updated version of the package.

### Install GDAL

`teamlucc` uses the `gdalUtils` package to facilitate fast image reprojection 
and mosaicking. `gdalUtils` requires having a local GDAL installation. Follow 
the below steps to install GDAL on your system:

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

### (optional) Install IDL
[IDL](http://www.exelisvis.com/ProductsServices/IDL.aspx) is required for 
running the IDL cloud fill and Landsat 7 SLC-off gap fill routines in 
`teamlucc`. There are two native R cloud fill routines that can be used without 
an IDL license.

## Using teamlucc

For more information on using `teamlucc`, see the online help in R, and the 
[`teamlucc` webpage](http://www.azvoleff.com/teamlucc). The webpage includes 
examples of a number of specific applications of `teamlucc`, including:

* [Filtering and downloading Landsat 
  scenes](http://www.azvoleff.com/articles/filtering-landsat-with-teamlucc)

* [Preprocessing imagery and 
  DEMS](http://www.azvoleff.com/articles/preprocessing-imagery-with-teamlucc)

* [Cloud removal](http://www.azvoleff.com/articles/cloud-removal-with-teamlucc)

* [Image 
  classification](http://www.azvoleff.com/articles/image-classification-with-teamlucc)

## Installing `teamlucc` development version

If you want the very latest version of `teamlucc` (though be aware this version 
might not install as it is not as well tested as the stable version), you can 
install from the development branch by typing:

```R
library(devtools)
install_github('azvoleff/teamlucc', ref="development")
```

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
