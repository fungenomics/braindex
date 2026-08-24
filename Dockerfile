FROM rocker/shiny:3.6.3

ENV RENV_PATHS_CACHE=/renv/cache
ENV RENV_PATHS_LIBRARY=/renv/library
ENV RENV_ACTIVATE_PROJECT=FALSE

RUN sed -i 's|deb.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    sed -i 's|security.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    sed -i '/buster-updates/d' /etc/apt/sources.list

RUN apt-get update && apt-get install -y \
    curl \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libuv1-dev \
    build-essential

COPY renv.lock renv.lock

RUN Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
RUN R -e 'remotes::install_version("renv", version = "0.15.4", repos="https://cloud.r-project.org")'

RUN R -e 'renv::init(bare=TRUE)'
# Restore the R packages that need to be installed from the pre-created renv.lock
RUN R -e "renv::restore(lockfile = 'renv.lock', prompt = FALSE)"
RUN R -e "renv::install('R.utils', version = '2.13.0', repos = 'https://cloud.r-project.org')"

# Register /renv/library/R-3.6/x86_64-pc-linux-gnu/ on every R session's search path, regardless of user or cwd
RUN echo '.libPaths(c("/renv/library/R-3.6/x86_64-pc-linux-gnu/", .libPaths()))' >> $(R RHOME)/etc/Rprofile.site

RUN chmod -R a+rX /renv

EXPOSE 3838
