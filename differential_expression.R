library(tidyr)
library(tibble)
library("DEP")
library("dplyr")
library(SummarizedExperiment)
library(ggplot2)
library(collections)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("/Users/user/Documents/proteomics")

AYA<-read.csv("proteomics_gene_names.csv")