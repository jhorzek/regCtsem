FROM rocker/tidyverse

RUN R -e "install.packages('ctsemOMX',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Rcpp',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('optimx',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Rsolnp',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('numDeriv',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('testthat',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('knitr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('rmarkdown',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN apt-get update && apt-get install -y qpdf
