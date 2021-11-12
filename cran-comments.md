## Submission

This fixes the errors found at https://cran.r-project.org/web/checks/check_results_SpatialEpi.html

Version: 1.2.5
Check: examples
Result: ERROR
    Running examples in 'SpatialEpi-Ex.R' failed
    The error most likely occurred in:
    
    > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
    > ### Name: EBpostthresh
    > ### Title: Produce the probabilities of exceeding a threshold given a
    > ### posterior gamma distribution.
    > ### Aliases: EBpostthresh
    >
    > ### ** Examples
    >
    > data(scotland)
    > Y <- scotland$data$cases
    > E <- scotland$data$expected
    > ebresults <- eBayes(Y,E)
    > #Find probabilities of exceedence of 3
    > thresh3 <- EBpostthresh(Y, E, alpha=ebresults$alpha, beta=ebresults$beta, rrthresh=3)
    > mapvariable(thresh3, scotland$spatial.polygon)
    Warning in wkt(obj) : CRS object has no comment
    Error in rgdal::OSRIsProjected(obj) : Can't parse user input string
    Calls: mapvariable ... is.projected -> is.projected -> is.projected -> <Anonymous>
    Execution halted
Flavors: r-devel-linux-x86_64-debian-clang, r-devel-linux-x86_64-debian-gcc


## Test environments

* local macOS install, R 4.0.3
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


