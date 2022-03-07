FROM rocker/shiny:latest
#FROM bioconductor/bioconductor_docker:devel

LABEL maintainer="Chieh-Yu Lee <chiehyulee97@gmail.com>"

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

RUN R -e "install.packages(pkgs=c('shiny','shinybusy','shinythemes','shinydashboard','shinyjs','shinyWidgets','shinycssloaders','dplyr','Seurat','plotly','reshape2','ggplot2','gridExtra','RColorBrewer'),dependencies=TRUE, repos='http://cran.rstudio.com/')"

#RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install("biomaRt")'


RUN addgroup --system app \
    && adduser --system --ingroup app app
WORKDIR /app
COPY app .
RUN chown app:app -R /app
USER app
EXPOSE 3939
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3939)"]
