## Submission

This address the following from Prof Brian Ripley

- you are shipping a log file SpatialEpi/vignettes/Manuscript/Manuscript.log
- You really must not recommend the use of MacPorts.  Proper, static builds of gdal and geos are available from mac.r-project.org, and binary macOS packages for rgeos and rgdal are available,  And you do not need to use a terminal to download and install packages: it can and probably should be done from inside R, either via install.packages() or from the GUI menus.
- Similarly, Windows users can just use binary packages. OSGeo4W is not needed even when installing from source.
- That vignette contains only misinformation and should be removed.
- Why are you including dirs Manuscript and Old_Vignette in the sources?


## Test environments

* local macOS install, R 4.0.3
* win-builder (release, devel, oldrelease)
* GitHub Actions
    + macOS: latest
* Rhub via devtools::check_rhub(env_vars=c(R_COMPILE_AND_INSTALL_PACKAGES = "always"))
    + 


## R CMD check results

* Warnings
    + I got the following warning: Found the following (possibly) invalid URLs: URL: https://developer.apple.com/download/all/. However I'm able to access the link without issue.
* Notes
    + Version contains large components (1.2.5.9000)

