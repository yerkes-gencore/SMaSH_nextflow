#!/usr/bin/env Rscript

# Plot heatmap of the exponents of SMaSH pvalues (because the pvalue floats are too small for R to handle for many datasets)

## Pass in args for plot size and text size
args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one arguments: if not, set default plot width (in cm)
if (length(args)==0) {
  args[1] <- 30
} else if (length(args)==1) {
  # default plot test size
  args[2] = 4
}

# plotWidth <- args[1]
# plotTextSize <- args[2]
plotWidth <- 30
plotTextSize <- 8

## Setup
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(forcats)
library(ggplot2)

## Read in data
## I removed the file paths from the pval_out.txt file in "smash_out" and save it here.
pval_out <- read_tsv("pval_out.txt", col_types = paste(rep("c",360), collapse=''))

pval_out <- pval_out %>%
  mutate(across(ends_with(".bam"), function(x) str_extract(x, pattern = "E-[0-9]+")),
         across(ends_with(".bam"), function(x) as.numeric(str_extract(x, pattern = "[0-9]+"))))

to_plot <- pval_out %>%
  rename(sample_A = ".") %>%
  pivot_longer(cols = ends_with(".bam"),
               names_to = "sample_B",
               values_to = "pval_neg_exp") %>%
  # mutate(sample_A = str_replace_all(sample_A, "-", "_"),
  #        sample_B = str_replace_all(sample_B, "-", "_")) %>%
  mutate(sample_A = str_replace(str_replace(sample_A, "-", "_"),"-","_"),
         sample_B = str_replace(str_replace(sample_B, "-", "_"),"-","_")) %>%
  separate_wider_delim(cols = c("sample_A"),
                       delim = "_",
                       names = c("projID_A", "sampID_A", "indivID_A", "repID_A"),
                       too_many = "drop") %>%
  separate_wider_delim(cols = c("sample_B"),
                       delim = "_",
                       names = c("projID_B", "sampID_B", "indivID_B", "repID_B"),
                       too_many = "drop") %>%
  mutate(self = indivID_A == indivID_B,
         proj_sampID_A = paste0(projID_A, sampID_A),
         proj_sampID_B = paste0(projID_B, sampID_B),
         indiv_repID_A = paste0(indivID_A, "_", repID_A),
         indiv_repID_B = paste0(indivID_B, "_", repID_B)) %>%
  select(proj_sampID_A, proj_sampID_B, everything())

p <- to_plot %>%
  ggplot(data=., aes(x=fct_reorder(indiv_repID_A, projID_A), 
                     y=fct_reorder(indiv_repID_B, projID_B), 
                     fill=pval_neg_exp)) +
  geom_raster() +
  scale_fill_viridis_c(option="magma") +
  theme(axis.text = element_text(size=plotTextSize,face = "bold"),
        axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5),
        aspect.ratio = 1) 

ggsave(filename = "SMaSH_pval_heatmap.png", plot = p, device = "png", width = plotWidth, height = plotWidth, units = "cm")

