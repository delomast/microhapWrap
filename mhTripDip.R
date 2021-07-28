#!/usr/bin/env Rscript
# ploidy inference (triploid vs diploid)

cmdArgs <- commandArgs(trailingOnly=TRUE)
if(length(cmdArgs) != 6) stop("Internal error: wrong number of arguments to mhTripDip.R")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tripsAndDipR))

prefix <- cmdArgs[6]
cmdArgs <- as.numeric(cmdArgs[1:5])
if(cmdArgs[5] < cmdArgs[4]) stop("Thresholds not consistent.")

readCounts <- read_tsv("snp_counts.txt", col_types = "ccdd") %>%
	filter(!grepl("NTCTray", Indiv))

ref <- readCounts %>% select(-Count2) %>% spread(Locus, Count1) %>% arrange(Indiv) %>% select(Indiv, sort(colnames(.)))
alt <- readCounts %>% select(-Count1) %>% spread(Locus, Count2) %>% arrange(Indiv) %>% select(Indiv, sort(colnames(.)))

ref <- as.data.frame(ref)
rownames(ref) <- ref$Indiv
ref <- ref %>% select(-Indiv) %>% as.matrix()

alt <- as.data.frame(alt)
rownames(alt) <- alt$Indiv
alt <- alt %>% select(-Indiv) %>% as.matrix()

ploid <- tripsAndDip(counts = ref, counts_alt = alt, h = rep(1, ncol(ref)), eps = rep(.01, ncol(ref)),
							min_reads = cmdArgs[1], min_loci = cmdArgs[2], binom_p_value = cmdArgs[3])
ploid$genPloidy <- "U"
ploid$genPloidy[!is.na(ploid$LLR) & ploid$LLR < cmdArgs[4]] <- "2n"
ploid$genPloidy[!is.na(ploid$LLR) & ploid$LLR > cmdArgs[5]] <- "3n"

ploid <- ploid %>% mutate(fullName = sample_name, 
					  sample_name = gsub("^initial|^f[0-9]|^rr?[0-9]|^qc", "", sample_name)) %>%
	select(fullName, sample_name, genPloidy, everything()) %>% rename(`Individual Name` = sample_name)

write.csv(ploid, paste0(prefix, "_genPloidy.txt"), row.names = FALSE, quote = FALSE)
