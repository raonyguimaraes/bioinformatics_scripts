# sudo apt-get install libxml2-dev
# sudo apt-get install r-cran-xml
# 
# sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
# sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
# 
# sudo apt-get update
# sudo apt-get install r-base
# 
# sudo apt-get install libcurl4-gnutls-dev

source("https://bioconductor.org/biocLite.R")
biocLite("PureCN")

library(PureCN)

browseVignettes("PureCN")

#sudo apt-get install libssl-dev

install.packages("devtools")
library(devtools)
install_github("Bioconductor-mirror/PureCN")

library(PureCN)
