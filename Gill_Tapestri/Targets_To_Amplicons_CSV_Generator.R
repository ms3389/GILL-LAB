# Generate Amplicon list for use in Tapestri Designer

# Load base dependencies
library(tidyverse)
library(magrittr)
#Read in data
Targets <- read.delim("~/Gill_Tapestri/Targets.txt")
# View(Targets)

#Subset Data by sites with 2 or less mismatches
Targets <- tibble(subset(Targets, subset = Targets$Mismatches <= 2))

#remove columns unnecessary for Tapestri Amplicon Panel Designer, and match the downstream query for merging
Targets %<>%
  select_('Chromosome', 'Position', 'Direction', 'DNA', 'crRNA', 'Mismatches') %>%
  mutate_all(~gsub("chr", "", .)) %>%
   rename_(chromosome_name = 'Chromosome') %>% 
  rename_(strand = 'Direction') 
Targets$Position <- as.integer(Targets$Position)
Targets$strand <- gsub(x = Targets$strand, pattern = "+", replacement = "1", fixed = TRUE)
Targets$strand <- gsub(x = Targets$strand, pattern = "-", replacement = "-1", fixed = TRUE)
Targets$strand <- as.integer(Targets$strand)
#View(Targets)

#The amplicon designer requires a specific format to the CSV, described here:https://support.missionbio.com/hc/en-us/articles/360038037893
#Needed columns are: gene, region, dbSNP, COSMIC, and HGVS

#Use biomaRt to generate metadata
require(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
attr <- c("chromosome_name", "strand", "start_position", "end_position", "hgnc_symbol", "percentage_gene_gc_content",  "description")
bm <- getBM(attributes = attr, mart = ensembl)

#create function for finding matches for gene locations and gene ranges
require(IRanges)
source("~/Gill_Tapestri/is_X_between_Y_n_Z.R")
is_X_between_Y_n_Z(Targets$Position, bm$start_position, bm$end_position) -> Pairs

#Format Pairs Table and remove unnecessary columns
Pairs %>% 
  select_('first.start', 'second.start', 'second.end') -> Pairs

#rename columns
Pairs %>% 
  rename(first.start = 'Position') %>%
  rename(second.start = 'start_position') %>% 
  rename(second.end = 'end_position') -> Pairs

#Combine Pair info with biomaRt annotations
left_join(Pairs, bm, by = c('start_position', 'end_position')) -> Pairs
#Combine Target info with Pairs Info
left_join(x = Targets, Pairs, by = c('Position', 'strand', 'chromosome_name')) -> Targets
View(Targets)

#Write to CSV
#write_csv(Targets, "~/Gill_Tapestri/Amplicons.csv")

