## Resubmission

This is a resubmission of SpatialEpi v1.2.4 which I originally submitted on 2021-07-27. In this version I have addressed all NOTES and ERRORS in https://cran.r-project.org/web/checks/check_results_SpatialEpi.html:

Version: 1.2.4
Check: dependencies in R code
Result: NOTE
    Namespaces in Imports field not imported from:
     ‘dplyr’ ‘ggplot2’ ‘grDevices’ ‘graphics’ ‘methods’ ‘sf’ ‘stats’
     All declared Imports should be used.
Flavors: r-devel-linux-x86_64-fedora-clang, r-devel-linux-x86_64-fedora-gcc, r-devel-windows-x86_64-gcc10-UCRT, r-patched-solaris-x86, r-release-macos-arm64, r-release-macos-x86_64, r-oldrel-macos-x86_64


Version: 1.2.4
Check: examples
Result: ERROR
    Running examples in ‘SpatialEpi-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: NYleukemia_sf
    > ### Title: Upstate New York Leukemia
    > ### Aliases: NYleukemia_sf
    > ### Keywords: datasets
    >
    > ### ** Examples
    >
    >
    > # Static map of NY Leukemia rate per county
    > library(ggplot2)
    > ggplot(NYleukemia_sf) +
    + geom_sf(aes(fill= cases/population)) +
    + scale_fill_gradient(low = "white", high = "red")
    Warning in CPL_transform(x, crs, aoi, pipeline, reverse, desired_accuracy, :
     GDAL Error 1: No PROJ.4 translation for source SRS, coordinate transformation initialization has failed.
    Error in CPL_transform(x, crs, aoi, pipeline, reverse, desired_accuracy, :
     OGRCreateCoordinateTransformation() returned NULL: PROJ available?
    Calls: <Anonymous> ... st_transform.sfc -> st_sfc -> structure -> CPL_transform
    Execution halted
Flavor: r-patched-solaris-x86




## Test environments

* local macOS install, R 4.1.0
* win-builder (release, devel, oldrelease)
* GitHub Actions
    + ubuntu-16.04: latest
    + windows: latest
    + macOS: latest, devel
* Rhub via devtools::check_rhub(env_vars=c(R_COMPILE_AND_INSTALL_PACKAGES = "always"))
    + Ubuntu Linux 20.04.1 LTS, R-release, GCC
    + Windows Server 2008 R2 SP1, R-devel, 32/64 bit
    + Oracle Solaris 10, x86, 32 bit, R-release


## R CMD check results

* Warnings
    + I got the following warning: Found the following (possibly) invalid URLs: URL: https://developer.apple.com/download/all/. However I'm able to access the link without issue.
* NOTES only for Oracle Solaris 10, x86, 32 bit, R-release
    + Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.
    + Compilation used the following non-portable flag(s): '-march=pentiumpro'



## Comments

I am switching the package maintainer from albert@stat.washington.edu to albert.ys.kim@gmail.com