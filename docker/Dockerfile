# Docker container based on R 3.6.3 - note that some packages have moved to R 4.x since and have broken dependency links for R 3.6.x so we here hardcode package installation order and Bioconductor version accordingly
FROM rocker/verse:3.6.3
LABEL MAINTAINER="github.com/ftwkoopmans/msdap"

# libpoppler is a pdftools requirement
# libnetcdf is required downstream by R package ncdf4, which is an upstream dependency
RUN apt-get update && apt-get install -y --no-install-recommends \
  libnetcdf-dev \
  netcdf-bin \
  libpoppler-cpp-dev

# install missing tinytex packages
RUN /opt/TinyTeX/bin/*/tlmgr install ifxetex ifluatex oberdiek graphics graphics-cfg graphics-def

# hardcoded install of CRAN dependencies; bugfix for failure to perform dependency checks/installs @ Bioconductor Msnbase  &  due to caching in `docker build` this saves a lot of time when updating the MS-DAP package (since dependencies are all cached after this line)
RUN R -e "install.packages(c('devtools', 'tidyverse', 'data.table', 'arrangements', 'bit64', 'crayon', 'DBI', 'doParallel', 'forcats', 'foreach', 'formatR', 'ggpubr', 'ggrepel', 'patchwork', 'gtools', 'knitr', 'xtable', 'lme4', 'MASS', 'Matrix', 'matrixStats', 'nloptr', 'openxlsx', 'pdftools', 'pROC', 'purrr', 'colorspace', 'readr', 'rmarkdown', 'RSQLite', 'scales', 'stringr', 'styler', 'tinytex', 'viridis', 'testthat', 'iq'), repos = 'https://cloud.r-project.org')"

# bugfix; directly link to proper R 3.6.3 version of this package. If we add MALDIquant to above install.packages() there is an error; package ‘MALDIquant’ is not available (for R version 3.6.3)
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/MALDIquant/MALDIquant_1.19.3.tar.gz')"

# explicitly use BiocManager for Bioconductor dependencies
# (at the time of testing, seems more robust than relying on devtools to install our Bioconductor dependencies  @  ProtGenerics/MSnbase)
RUN R -e "BiocManager::install(c('ProtGenerics', 'MSnbase', 'limma', 'vsn', 'pcaMethods', 'DEqMS'), update=F, ask=F, version = '3.10')"

# dependencies from github
RUN R -e "devtools::install_github('zimmerlab/MS-EmpiRe', upgrade = 'never')"
RUN R -e "devtools::install_github('vdemichev/diann-rpackage', upgrade = 'never')"

# install msdap R package
###RUN R -e "devtools::install_github('ftwkoopmans/msdap@1.0', upgrade = 'never')"
RUN mkdir /msdap
COPY temp/msdap_1.0.tar.gz /msdap/
RUN R -e "devtools::install_local('/msdap/msdap_1.0.tar.gz', upgrade = 'never')"

# make sure the rstudio user can install and upgrade packages
RUN chown -R rstudio:rstudio /usr/local/lib/R/library
RUN chown -R rstudio:rstudio /usr/local/lib/R/site-library

# optionally, install example data
COPY temp/exampledata /exampledata/
RUN chown -R rstudio:rstudio /exampledata/*

# our standard mount point: shared directory between host and container where we can store data
RUN mkdir -m 777 /data

# change default working directory for RStudio server to /data
RUN echo "session-default-working-dir=/data" >> /etc/rstudio/rsession.conf
RUN echo "session-default-new-project-dir=/data" >> /etc/rstudio/rsession.conf

