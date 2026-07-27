FROM rocker/r-ver:4.4.1
RUN apt-get update && apt-get install -y --no-install-recommends \
      zlib1g-dev git && rm -rf /var/lib/apt/lists/*
RUN R -e 'install.packages(c("Matrix","remotes"), repos="https://cloud.r-project.org")' \
 && R -e 'install.packages("data.table", type="source", repos="https://cloud.r-project.org")'
# git clone https://github.com/gjhunt/hspe.git hspe_src
COPY hspe_src/lib_hspe /tmp/hspe_pkg
RUN Rscript -e 'remotes::install_local("/tmp/hspe_pkg"); \
                library(hspe); library(Matrix); library(data.table); \
                cat("ALL THREE LOAD OK\n")'