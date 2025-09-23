# Select 2023 NEON Samples to use

#### 1. Setup #####
# Libraries
library(tidyverse)
library(ggside)
library(ozmaps)
library(sf)
library(scales)
library(arules)
library(phyloNEON)

# Functions
`%notin%` <- Negate(`%in%`)

# Working directory
setwd("~/Documents/GitHub/Oxygen/")

#### 2. Explore ####
# Metadata sent by Hugh Cross Jan 29th, 2025 (wasn't publicly available yet)
# Total without filtering: 263 from 27 sites
# Non-wetland: 242 from 27 sites
# Non-wetland and depth starting at 0: 233 from 27 sites
d <- read.delim("data/neon_soil_metagenome_samples_2023.tsv") %>%
  drop_na(Genome.Size....assembled) %>% # Remove those not sequenced yet (n = 9)
  separate(depth..meters, remove = F, sep = " - ", into = c("topDepth", "botDepth")) %>%
  mutate(topDepth = as.numeric(topDepth),
         botDepth = as.numeric(botDepth)) %>%
  filter(ecosystem_subtype != "Wetlands") %>% # Remove wetlands (n = 21)
  filter(topDepth == 0) %>% # Remove those not starting at 0 (n = 9)
  mutate(Site = substr(sample.name, start = 0, stop = 4))
sites <- as.data.frame(table(d$Site))
nrow(sites) # 27 sites
range(sites$Freq) # 3 - 13 samples per site
range(d$topDepth) # 0 to 0
range(d$botDepth) # 0.007 to 0.36
hist(d$botDepth)



#### 3. MetaG ####
# Use Hugh's R package to find the MetaG at JGI
View(neon.metaDB)
openIMG('BONA_013-O-20230710-COMP-DNA1')

# Example - get BONA samples with > 10 bins
query1 <- neon.metaDB %>%
  dplyr::filter(siteID == 'BONA') %>%
  dplyr::filter(metaBATbinCount > 10)
dim(query1) 
# [1] 11 32
# you need to input a list, so specify the dnaSampleID of your query:
openIMG(query1$dnaSampleID)

# Example - Wind River, 2021
query2 <- neon.metaDB %>%
  dplyr::filter(siteID == 'WREF') %>%
  dplyr::filter(as.Date(collectDate, "%Y%m%d") > '2021-01-01' & as.Date(collectDate, "%Y%m%d") < '2021-12-31') %>%
  dplyr::pull(dnaSampleID)
# check how many samples (using length() this time)
length(query2) 
# [1]  7
# open the pages for each site. Note here the query is a vector so we do not have to input query2$dnaSampleID
openIMG(query2)

# Example - JGI samples (better depth)
query3 <- neon.metaDB %>%
  dplyr::filter(`Sequencing Center` == 'DOE Joint Genome Institute  (JGI)')
# or by ITS Proposal ID
query3a <- neon.metaDB %>%
  dplyr::filter(`ITS Proposal ID` %in% c('509938','509462'))
