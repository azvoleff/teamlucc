# TEAM Data Processing Tools (R)

[![Build Status](https://travis-ci.org/azvoleff/teamr.png)](https://travis-ci.org/azvoleff/teamr)

## Overview

The `teamr` package is designed to facilitate the work of the Tropical Ecology 
Assessment and Monitoring (TEAM) Network. The [TEAM 
Network](http://www.teamnetwork.org/) is a global network of sites in tropical 
forests focused on collecting standardized real-time data measuring 
the response of tropical forests to changing climate, land cover and land use, 
and population.

The TEAM Data Processing Tools are a collection of scripts written in R and 
Python to facilitate the processing and analysis of TEAM data and remote 
sensing imagery of TEAM sites for measuring change in ecosystems from the 
plot-landscape scales. This package includes the R scripts in the toolkit. The 
[teampy package](https://github.com/azvoleff/teampy) includes the Python 
scripts that are part of the toolkit. Both packages are required in order to 
use the toolkit.

The toolkit is still under active development. Follow the [TEAM 
website](http://www.teamnetwork.org/) for news, and the [toolkit project page
on github](https://github.com/azvoleff/teamr) for the latest version of the 
code.

## Package Installation

The easiest way to install the (beta version) of the `teamr` package is using 
the [devtools](http://cran.r-project.org/web/packages/devtools/index.html) 
package from Hadley Wickham. After installing `devtools` from CRAN, type:

```R
install_github('teamr', username='azvoleff')
```

at the R prompt to install `teamr`. As `teamr` contains C++ code, you will need 
to have the the appropriate version of 
[Rtools](http://cran.r-project.org/bin/windows/Rtools/) installed for your 
version of R.

## Author Contact Information

[Alex Zvoleff](mailto:azvoleff@conservation.org)  
Postdoctoral Associate  
Tropical Ecology Assessment and Monitoring (TEAM) Network  
Conservation International  
2011 Crystal Dr. Suite 500  
Arlington, VA 22202  
USA
