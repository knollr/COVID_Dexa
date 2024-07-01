library(ggplot2)
library(ggfortify)
library(ggsignif)
library(ggExtra)
library(RColorBrewer)
library(minfi)
library(limma)
library(knitr)
library(bookdown)
library(Gviz)
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(stringr)
library(readxl)
# not working installing from source?library(DMRcate)
library(devtools)
library(factoextra)
library(ChIPseeker)
library(rtracklayer)
library(MethylSeekR)
library(RnBeads)
library(fuzzyjoin)
# additional pipes, like %$%
library(magrittr)
library(dendextend) # more options for dendograms
library(pheatmap)
library(wateRmelon)
library(tidyverse)
library(Biobase)
library(here)
library(GEOquery)
library(sva)
library(ggrepel)
library(fgsea)
library(methylGSA)
library(DMRcate)
select <- dplyr::select
barplot <- graphics::barplot
theme_set(theme_classic())
set.seed(123)

R_path <- file.path("/home/rstudio/analysis/Code/R_projects/R-general-functions")
source(paste0(R_path, "/plot.R"))
source(paste0(R_path, "/function.R"))

# color palette for the following analysis
pal <- brewer.pal(8,"Dark2")

# color palette for a variable with more factors
pal20 <- colorRampPalette(brewer.pal(8, "Dark2"))(20)

# colors used for Dexa paper
pal_Dexa <- c("Non-Responder" = "#5B7674", "Responder" = "#57B391")

# paths to relevant directories
codeDirectory <-file.path("/home/rstudio/analysis/Code/R_projects/EPIC_Analysis/Dexa-repo")
dataDirectory <-file.path("/home/rstudio/analysis/Data_s_drive/Epic/Dexa/idat_files")
general_dataDirectory <-file.path("/home/rstudio/analysis/Data_s_drive/Epic/general_data")
outputDirectory <-file.path("/home/rstudio/analysis/Output_Analysis/EPIC_Output/Dexa_Output")