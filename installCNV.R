local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)
})


install.packages('devtools')

library(devtools)

devtools::install_github(repo = "honzee/RNAseqCNV")