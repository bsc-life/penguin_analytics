
# this is a fundtion that install packages if not existing already and loads them to the R sesssion.
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}


# my required packages
usePackage("dplyr")
usePackage("tidyr")
usePackage("tidyverse")
usePackage("ggplot2")
usePackage("gridExtra")
usePackage("purrr")
usePackage("jsonlite")
usePackage("Hmisc")
usePackage("igraph")
usePackage("ggpubr")
usePackage("stringr")
usePackage("scales")
usePackage("sqldf")
usePackage("grid")
usePackage("data.table")
