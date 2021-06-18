# braindex: brain developmental expression portal

![](www/README_img.png)


## About

This repository hosts the code for a directory of shiny apps created to visualize and interrogate
brain developmental gene expression data. Each app is designed for a different type of data.


## Directory organization

In this directory, several resources pertain globally to the portal:

* `www`: Special directory to expose resources for the website, including the CSS
specification for the index page and individual apps, the portal and lab logos
* `www/layout`: Files containing the HTML for various page parts (e.g. header, footer)
inserted into each app for custom styling and navigation, using the helper functions
in `www/ui_functions.R`
* `index.html`: Contains the HTML for the homepage, which links to the apps hosted in the directory
* `style.R`: An R script which can be sourced by each app, containing styles shared
across apps. Currently contains a ggplot2 theme, `theme_min()`.


There is one directory for each app, containing at minimum two R files `server.R` and `ui.R`,
and a data directory `data`.

* `lifespan`: Containts the R code for the shiny app for exploring gene expression in the
brain across the lifespan, based on data from the BrainSpan project
* `clusters`: Contains the R code for the shiny app for exploring gene expression
based on clusters from the scRNAseq atlas of the developing mouse forebrain and pons
generated our lab


## Datasets which can be explored with the app

* Miller, J.A. et al. (2014) Transcriptional landscape of the prenatal human brain, Nature 508: 199-206. [doi:10.1038/nature13185](https://doi.org/10.1038/nature13185.) (Data is © 2010 Allen Institute for Brain Science. Allen Human Brain Atlas. Available from: https://www.brainspan.org/)

* Jessa, S. et al. (2019) Stalled developmental programs at the root of pediatric brain tumors, Nature Genetics 51: 1702-1713. [doi:10.1038/s41588-019-0531-7](https://doi.org/10.1038/s41588-019-0531-7)

