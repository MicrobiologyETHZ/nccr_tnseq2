# Title     : TODO
# Objective : TODO
# Created by: ansintsova
# Created on: 03/06/2021

# Title     : TODO
# Objective : TODO
# Created by: ansintsova
# Created on: 23.01.21
#
#install.packages(c('BiocManager','Hmisc'),
#                 dependencies='Depends',
#                repo = "http://cran.us.r-project.org")
#BiocManager::install(c('DESeq2','apeglm'))
#


library(DESeq2)

args <- commandArgs(trailingOnly=TRUE)

sdf_file <-args[1]
edf_file <- args[2]

experiment <- unlist(strsplit(basename(sdf_file), '_sdf'))[[1]]
prefix <- unlist(strsplit(sdf_file, '_sdf'))[[1]]
sdf <- read.csv(sdf_file, row.names=1)
edf <- read.csv(edf_file, row.names=1)

calculate_vst_transform <- function(sdf, edf, prefix){
    print('Design: day')
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=edf, colData=sdf, design=~day)
    dds$day <- relevel(dds$day, ref = "d0")
    dds <- DESeq2::DESeq(dds)

    vsd <- DESeq2::vst(dds, blind=TRUE)
    write.csv(assay(vsd), paste0(prefix, "_vst.csv"))
}

calculate_vst_transform(sdf, edf, prefix)

