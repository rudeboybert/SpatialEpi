## Test environments

* local macOS install, R 4.1.0
* win-builder (release, devel, oldrelease)
* GitHub Actions
    + ubuntu-16.04: latest
    + windows: latest
    + macOS: latest, devel
* Rhub via devtools::check_rhub(env_vars=c(R_COMPILE_AND_INSTALL_PACKAGES = "always"))
    + Fedora Linux, R-devel, clang, gfortran
    + Ubuntu Linux 20.04.1 LTS, R-release, GCC
    + Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results

There were no ERRORs or WARNINGs. 

## Comments

I am switching the package maintainer from albert@stat.washington.edu to albert.ys.kim@gmail.com