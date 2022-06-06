FROM satijalab/seurat

LABEL maintainer="Chieh-Yu Lee <chiehyulee97@gmail.com>"

# system libraries of general use
## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

RUN R -e "install.packages(pkgs=c('shiny','shinybusy','shinythemes','shinydashboard','shinyjs','shinyWidgets','shinycssloaders','dplyr','plotly','reshape2','ggplot2','gridExtra','RColorBrewer'),dependencies=TRUE, repos='http://cran.rstudio.com/')"



RUN addgroup --system app \
    && adduser --system --ingroup app app

WORKDIR /app

COPY app .
RUN chown app:app -R /app
USER app
EXPOSE 3939
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3939)"]
