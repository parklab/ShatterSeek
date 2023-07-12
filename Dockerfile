# docker build -t shahcompbio/shatterseek .
# docker push shahcompbio/shatterseek

FROM rocker/r-ver:4.3.1

RUN apt-get update

# for R install
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev zlib1g-dev curl git

# R install
RUN apt-get install -y libbz2-dev liblzma-dev 
RUN R -e "install.packages(c('devtools'))"
RUN R -e "install.packages(c('BiocManager'))"
RUN R -e "install.packages(c('MASS'))"
RUN R -e "install.packages(c('argparse'))"
RUN R -e "install.packages(c('tidyverse'))"
RUN R -e "install.packages(c('gridExtra'))"
RUN R -e "BiocManager::install('BiocGenerics')"
RUN R -e "BiocManager::install('graph')"
RUN R -e "BiocManager::install('S4Vectors')"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('IRanges')"

RUN R -e "devtools::install_github('parklab/ShatterSeek')"

RUN R -e "install.packages(c('optparse'))"
