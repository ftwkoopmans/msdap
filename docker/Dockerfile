FROM rocker/verse:4.4.0
LABEL MAINTAINER="github.com/ftwkoopmans/msdap"


### system dependencies
# libpoppler is a pdftools requirement
# libnetcdf is required downstream by R package ncdf4, which is an upstream dependency
RUN apt-get update && apt-get install -y --no-install-recommends \
  libnetcdf-dev \
  netcdf-bin \
  libpoppler-cpp-dev


### R dependencies; hardcode installation of main dependencies not in base image yet
### so we can reuse intermediate images in future updates where only MS-DAP will change

# tinytex packages (tinytex is already installed)
RUN R -e "tinytex::tlmgr_install(c('ifxetex', 'ifluatex', 'oberdiek', 'graphics', 'graphics-cfg', 'graphics-def', 'amsmath', 'latex-amsmath-dev', 'euenc', 'fontspec', 'tipa', 'unicode-math', 'xunicode', 'ltxcmds', 'kvsetkeys', 'etoolbox', 'xcolor', 'geometry', 'fancyvrb', 'framed', 'bigintcalc', 'bitset', 'etexcmds', 'gettitlestring', 'hycolor', 'hyperref', 'intcalc', 'kvdefinekeys', 'letltxmacro', 'pdfescape', 'refcount', 'rerunfilecheck', 'stringenc', 'uniquecounter', 'zapfding'))"

# CRAN R packages we need that aren't in this rocker/verse
RUN R -e "install.packages(c('archive', 'styler', 'formatR', 'pdftools', 'matrixStats', 'lme4', 'nloptr', 'pROC', 'iq', 'doParallel', 'foreach', 'missForest', 'ggpubr', 'ggrepel', 'patchwork', 'openxlsx'), repos = 'https://cloud.r-project.org')"

# BioConductor R packages
RUN R -e "BiocManager::install(c('ProtGenerics', 'MSnbase', 'limma', 'vsn', 'pcaMethods', 'DEqMS', 'BiocParallel', 'variancePartition'), update=F, ask=F)"

# GitHub R packages
RUN R -e "devtools::install_github('zimmerlab/MS-EmpiRe', upgrade = 'never')"
RUN R -e "devtools::install_github('vdemichev/diann-rpackage', upgrade = 'never')"

# 'variancePartition' package now requires 'remaCor' package version >= 0.0.15
RUN R -e "devtools::install_version('remaCor', '0.0.18', upgrade = 'never', repos = 'https://cloud.r-project.org')"

### MS-DAP
# from local devtools::build() file
RUN mkdir /msdap
COPY temp/msdap_1.0.7.tar.gz /msdap/
RUN R -e "devtools::install_local('/msdap/msdap_1.0.7.tar.gz', upgrade = 'never')"
# alternatively, install from github; RUN R -e "devtools::install_github('ftwkoopmans/msdap@1.0.7', upgrade = 'never')"


### finish setting up
# make sure the rstudio user can install and upgrade packages
RUN chown -R rstudio:rstudio /usr/local/lib/R/library
RUN chown -R rstudio:rstudio /usr/local/lib/R/site-library

# install example data
COPY temp/exampledata /exampledata/
RUN chown -R rstudio:rstudio /exampledata/*

# our standard mount point: shared directory between host and container where we can store data
RUN mkdir -m 777 /data

# change default working directory for RStudio server to /data
RUN echo "session-default-working-dir=/data" >> /etc/rstudio/rsession.conf
RUN echo "session-default-new-project-dir=/data" >> /etc/rstudio/rsession.conf
