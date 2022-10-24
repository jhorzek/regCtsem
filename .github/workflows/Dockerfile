FROM rocker/tidyverse

RUN R -e "install.packages('numDeriv',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN apt-get update && apt-get install -y qpdf
