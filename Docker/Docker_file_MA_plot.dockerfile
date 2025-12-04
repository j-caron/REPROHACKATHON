FROM rocker/r-base:4.4.0

## Install system dependencies required by Bioconductor packages
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean

## Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

## Install required Bioconductor packages
RUN R -e "BiocManager::install(c('KEGGREST', 'EnrichmentBrowser'))"

## Install CRAN packages
RUN R -e "install.packages(c('ggplot2', 'ggrepel', 'stringr'), repos='https://cloud.r-project.org')"

CMD ["/bin/bash"]
