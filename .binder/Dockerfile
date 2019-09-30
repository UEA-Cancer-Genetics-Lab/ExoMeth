FROM rocker/binder:3.5.3
LABEL maintainer='Shea Connell'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}



RUN wget https://github.com/Shedimus/ExoMeth/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-29/')); devtools::install_deps()"

RUN rm DESCRIPTION.1; exit 0
