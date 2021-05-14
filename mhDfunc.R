#!/usr/bin/env Rscript
# duplicate search with microhaps using duplicate search function
# in rubias

cmdArgs <- commandArgs(trailingOnly=TRUE)
if(length(cmdArgs) != 3) stop("Three (and only three) command line arguments should be present. ",
										"First is minimum proportion of genotypes present in both individuals. ", 
										"Second is minimum proportion of genotypes matching in both individuals. ",
										"Third is prefix for the output file")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rubias))


# minimum proportion of genotypes present in both individuals
minComparable <- as.numeric(cmdArgs[1])
# minimum proportion of genotypes that are identical
propMatch <- as.numeric(cmdArgs[2])
cat(minComparable, "minimum proportion of genotypes present in both individuals\n")
cat(propMatch, "minimum proportion of genotypes matching in both individuals\n")

prefix <- cmdArgs[3]

# get all DfuncItUp.txt files
inputs <- dir(pattern = "*_DfuncItUp.txt")

genos <- tibble()
for(i in inputs){
	genos <- genos %>% bind_rows(read_tsv(i, col_types = "cccccccc"))
}
genos <- genos %>% mutate(Allele1 = na_if(Allele1, "000"), 
								  Allele2 = na_if(Allele2, "000"),
								  sampleName = paste0(Type, Indiv))
md <- genos %>% select(sampleName, Indiv, Type, Status, Plate, Well) %>% distinct() %>%
	rename(Sample = Indiv)

rIn <- genos %>% select(sampleName, Locus, Allele1, Allele2) %>%
	gather(a, c, Allele1:Allele2) %>% 
	mutate(a = gsub("llele", "", a),
			 Locus = paste0(Locus, ".", a)) %>% select(-a) %>%
	spread(Locus, c) %>% select(sampleName, sort(colnames(.))) %>%
	rename(indiv = sampleName) %>%
	mutate(sample_type = "reference", repunit = "placeholder", 
			 collection = "placeholder") %>% 
	select(indiv, sample_type, repunit, collection, everything())

# duplicate search but don't print anything to the terminal
quiet <- capture.output(suppressMessages(suppressWarnings(dupTable <- close_matching_samples(rIn, gen_start_col = 5, 
											  min_frac_non_miss = minComparable,
											  min_frac_matching = propMatch))))

# add metadata
dupTable <- dupTable %>% select(num_non_miss, num_match, indiv_1, indiv_2)
colnames(md)[colnames(md) != "sampleName"] <- paste0(colnames(md)[colnames(md) != "sampleName"], "_1")
dupTable <- dupTable %>% left_join(md, by = c("indiv_1" = "sampleName"))
colnames(md)[colnames(md) != "sampleName"] <- gsub("1$", "2", colnames(md)[colnames(md) != "sampleName"])
dupTable <- dupTable %>% left_join(md, by = c("indiv_2" = "sampleName"))

# calculate additional columns
dupTable <- dupTable %>% mutate(propAlike = num_match / num_non_miss, numUnalike = num_non_miss - num_match,
						  concordance = Sample_1 == Sample_2, 
						  samNumApart = abs(as.numeric(gsub("^.+_", "", Sample_1)) - as.numeric(gsub("^.+_", "", Sample_2))))

write.csv(dupTable, paste0(prefix, "_dupSearchResults.csv"), row.names = FALSE, quote = FALSE)
