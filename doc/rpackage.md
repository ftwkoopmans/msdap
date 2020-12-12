
# Installing the MS-DAP R package

This guide helps you install the MS-DAP R package and all of its
software dependencies on your computer. If you are looking for the
*Dockerized* MS-DAP, [use this guide](docker.md).

## prequisites

**Windows**

    None, you may immediately go to the next section.

You’ll need to install the following system packages on Linux and macOS;
netcdf and poppler

**Ubuntu** (for Ubuntu 18.04 or older, see further
<https://github.com/ropensci/pdftools> ):

    sudo apt-get install libnetcdf-dev netcdf-bin libpoppler-cpp-dev

**Fedora**

    sudo dnf install netcdf-devel netcdf poppler-cpp-devel

**MacOS** first, install Homebrew if you haven’t: <https://brew.sh>

    brew install netcdf poppler automake
    export PKG_CONFIG_PATH=/usr/local/Cellar/zlib/1.2.8/lib/pkgconfig:/usr/local/lib/pkgconfig:/opt/X11/lib/pkgconfig

**troubleshooting**

*poppler* is a requirement for pdftools to function, which is used by
MS-DAP to create PDF reports, if you have any trouble installing please
check these resources:

  - <https://github.com/ropensci/pdftools>
  - <https://github.com/politza/pdf-tools>

## installing R (if you have not already)

An installation of R version 3.6 or newer is required, version 3.6.3 is
recommended (4.0 is not supported yet, some of the packages we rely on
are not compatible with 4.0 yet)

**Windows**

  - base R 3.6.3 @
    <https://cloud.r-project.org/bin/windows/base/old/3.6.3/>
  - Rtools35 @
    <https://cran.r-project.org/bin/windows/Rtools/history.html>

**Ubuntu**

    sudo apt-get install r-base r-base-dev

**Fedora**

    sudo yum install R

**MacOS** - download and install *R-3.6.3.nn.pkg* @
<https://cran.r-project.org/bin/macosx/> - make sure to select Tcl/Tk
and Textinfo during installation

## install RStudio Desktop

Download and install RStudio Desktop (free edition) from
<https://rstudio.com/products/rstudio/download/>

## install the MS-DAP R package

Start RStudio (restart if it was already up-and-running to get a clean
session) and run below code line-by-line so you can keep an eye on the
results of each intermediate step in the RStudio console.

Notes on windows;

1)  When prompted to optionally compile packages from source you can
    simply select ‘no’ (to install a readily available binary instead).
2)  As noted during the TinyTex installation, the two popup errors
    stating *“The code execution cannot proceed because luatex.dll …”*
    can be safely dismissed.

<!-- end list -->

``` r
### 1) setup required R packages
install.packages(c("devtools", "tidyverse", "tinytex", "BiocManager"))
# install LaTeX which is required to create the pdf report (you can safely dismiss both popup errors)
tinytex::install_tinytex()
# install packages from bioconductor
BiocManager::install(c('ProtGenerics', 'MSnbase', 'limma'), update=T, ask=F)

### 2) install MS-DAP R package directly from github
# first, make sure we don't halt on minor issues that can be ignored (eg; your R installation is a minor version behind)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# install MS-DAP R package and all required dependencies (but don't check for updates on all packages, to minimize complexity)
devtools::install_github("ftwkoopmans/msdap", upgrade = "never")

# 3) if all went well, you can now load the msdap package.
library(msdap)

# If you get any warnings, no problem; you only need to be on the lookout for errors
```

**troubleshooting**

  - if you encounter “installation path not writable” errors the user
    running R (/RStudio) may not have write access to the directory
    where the R packages are stored. Use the R command `.libPaths()` to
    print those locations.
      - For windows users, the most robust workaround is to start a
        command prompt as administrator (start \> type `cmd` \>
        right-click `run as administrator`), open a commandline R
        session (eg; on the cmd prompt, type `"C:\Program
        Files\R\R-3.6.3\bin\R.exe"`) and retry your package
        installation.
  - Further documentation for TinyTex is available at
    <https://yihui.org/tinytex/>
