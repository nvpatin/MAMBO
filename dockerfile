# Use the R studio image as a base
FROM rocker/r-ver:4.3.2

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=2023.12.0+369
ENV DEFAULT_USER=rstudio
ENV PANDOC_VERSION=default
ENV QUARTO_VERSION=default
ENV PASSWORD=mamboPW

RUN /rocker_scripts/install_rstudio.sh
RUN /rocker_scripts/install_pandoc.sh
RUN /rocker_scripts/install_quarto.sh

# Install JAGS dependency (system-wide)
RUN apt-get update && . /etc/environment \
  && wget sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-4.3.2.tar.gz -O jags.tar.gz \
  && tar -xf jags.tar.gz \
  && cd JAGS* && ./configure && make -j4 && make install

# Install other dependencies (system-wide)
RUN apt-get update && apt-get install -y libudunits2-dev  # for R "units" library 
RUN apt-get update && apt-get install -y cmake  # for R "car" library
RUN apt-get update && apt-get install -y build-essential  # for R "devtools" library
RUN apt-get update && apt-get install -y libcurl4-gnutls-dev libssl-dev # for R "devtools" library
RUN apt-get update && apt-get install -y libfontconfig1-dev # for R "devtools" library
RUN apt-get update && apt-get install -y libgit2-dev # for R "devtools" library
RUN apt-get update && apt-get install -y libxml2-dev # for R "devtools" library
RUN apt-get update && apt-get install -y libharfbuzz-dev libfribidi-dev # for R "devtools" library
RUN apt-get update && apt-get install -y libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev # for R "devtools" library
RUN apt-get update && apt-get install -y libgdal-dev # for R "sf" library


# Install R package dependencies
RUN R -e "install.packages('udunits2', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('httpuv')"
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages(c('rjags'), dependencies=TRUE)"
RUN R -e "install.packages(c('abind', 'runjags', 'dplyr', 'purrr', 'ggplot2'), dependencies=TRUE)"
RUN R -e "install.packages(c('car'),  dependencies=TRUE)"
RUN R -e "install.packages(c('swfscMisc'), repos='http://cran.rstudio.com/', dependencies=TRUE)"
RUN R -e "devtools::install_github('nvpatin/mambo/mambo')"
RUN echo "auth-timeout-minutes=20000" >> /etc/rstudio/rserver.conf


# Copy your R script or analysis files into the container
WORKDIR /home/app
COPY . .


EXPOSE 8787
CMD ["/init"]