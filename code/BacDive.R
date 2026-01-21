# Assemble BacDive Dataset and Genome Download List
# By Cliff Bueno de Mesquita, Fierer Lab, May 2025

#### Setup ####
# Libraries
library(tidyverse)
#install.packages("BacDive", repos="http://R-Forge.R-project.org")
#library(BacDive)
library(ggrepel)
library(ggbeeswarm)
library(cowplot)
`%notin%` <- Negate(`%in%`)
setwd("~/Documents/GitHub/Oxygen/")

# Test BacDive
bacdive <- open_bacdive("cliff.buenodemesquita@colorado.edu", "enterpassword")
example(fetch)
one_strain <- fetch(object = bacdive, ids = 5621)
one_strain
View(one_strain$results)
one_strain$results$`5621`$`Physiology and metabolism`$`oxygen tolerance`[2]$`oxygen tolerance`
one_strain$results$`5621`$`Sequence information`$`Genome sequences`[[1]]$`NCBI tax ID`
one_strain$results$`5621`$`Name and taxonomic classification`$phylum
one_strain$results[[1]]$`Name and taxonomic classification`$genus
two_strain <- fetch(object = bacdive, ids = 5621, 139709)

# Working Directory
setwd("~/Documents/GitHub/Oxygen/")

# Download all strains with oxygen tolerance data (May 20, 2025)
# Pie online says total count 24162
# But in the Advanced Search, it says hits 23255
# The downloaded file has 23255 with ID and O2 info



#### GTDB ####
# Bacterial metadata from GTDB r226
bacGT <- read_tsv("~/Desktop/Fierer/AEGIS/Oxygen/bac120_metadata_r226.tsv") %>% # 715230 genomes
  mutate(gb_accession = substr(ncbi_genbank_assembly_accession,
                               start = 1,
                               stop = nchar(ncbi_genbank_assembly_accession) - 2))
names(bacGT)
bacGT_rep <- bacGT %>%
  filter(gtdb_representative == TRUE) # 136646 representative genomes
length(unique(bacGT$ncbi_species_taxid)) # 75260
length(unique(bacGT$ncbi_taxid)) # 100433
length(unique(bacGT_rep$ncbi_species_taxid)) # 40675
length(unique(bacGT_rep$ncbi_taxid)) # 41876

bacGT_useful <- bacGT %>%
  dplyr::select(accession, checkm2_completeness, checkm2_contamination,
                checkm_completeness, checkm_contamination, gc_percentage, genome_size,
                gtdb_representative, gtdb_taxonomy, ncbi_assembly_level, ncbi_country,
                ncbi_genbank_assembly_accession, ncbi_genome_category, ncbi_isolation_source, 
                ncbi_refseq_category, ncbi_species_taxid, ncbi_taxid)



#### BacDive Download ####
# First Pass
# Keep main classes. 23161
# Keep type strains. 11108
d <- read.csv("data/advsearch_bacdive_2025-05-20.csv") %>%
  filter(!is.na(ID)) %>%
  filter(Oxygen.tolerance %in% c("obligate aerobe", "aerobe", "facultative anaerobe",
                                 "microaerophile", "anaerobe", "obligate anaerobe")) %>%
  filter(is_type_strain_header == 1)
table(d$Oxygen.tolerance)
split_list <- split(d, ceiling(seq_len(nrow(d)) / 100))
info <- list()



# Second Pass - don't filter by type strains. only keep 4 categories. only Bacteria. With genomes.
# Will filter to GTDB representatives later
# n = 8673
# n = 7124 using the 4 main categories
d_ncbi <- read.csv("data/advsearch_bacdive_2025-06-03_Dom_Phy_O2_GCA_ID.csv") %>%
  filter(!is.na(ID)) %>%
  filter(Oxygen.tolerance %in% c("obligate aerobe", "aerobe", "anaerobe", "obligate anaerobe")) %>%
  mutate(Genome.Sequence.database = "NCBI")

# Add img? But not linked to GTDB
d_img_phy <- read.csv("data/advsearch_bacdive_2025-06-03_Dom_Phy_O2_notGCA_ID.csv") %>%
  filter(!is.na(ID)) %>%
  filter(Oxygen.tolerance %in% c("obligate aerobe", "aerobe", "anaerobe", "obligate anaerobe")) %>%
  dplyr::select(ID, phylum)
d_img <- read.csv("data/advsearch_bacdive_2025-06-03_Dom_O2_notGCA_ID_IMG.csv") %>%
  filter(!is.na(ID)) %>%
  filter(Oxygen.tolerance %in% c("obligate aerobe", "aerobe", "anaerobe", "obligate anaerobe")) %>%
  left_join(., d_img_phy, by = "ID") %>%
  filter(ID %notin% d_ncbi$ID) %>%
  filter(Genome.Sequence.associated.NCBI.tax.ID %notin% d$Genome.Sequence.associated.NCBI.tax.ID) %>%
  dplyr::select(all_of(names(d_ncbi))) # 212 genomes
table(d_img$phylum)
sum(d_img$is_type_strain_header) # 192
sum(bacGT$accession %in% d_img$Genome.seq..accession.number) # 0

d <- d_ncbi %>%
  #rbind(., d_img) %>%
  mutate(phylum = recode(phylum,
                         "Proteobacteria" = "Pseudomonadota",
                         "Thermomicrobiota" = "Chloroflexota")) %>%
  mutate(Oxygen.tolerance = factor(Oxygen.tolerance,
                                   levels = c("obligate aerobe", "aerobe", 
                                              "anaerobe", "obligate anaerobe"))) %>%
  mutate(Oxygen2 = case_match(Oxygen.tolerance,
                              c("obligate aerobe", "aerobe") ~ "aerobic",
                              c("anaerobe", "obligate anaerobe") ~ "anaerobic"))

table(d$phylum)
sum(d$is_type_strain_header) # 6596
phylum_counts <- as.data.frame(table(d$phylum)) %>%
  arrange(desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = Var1))
phylum_counts_clear <- d %>%
  group_by(phylum, Oxygen2) %>%
  summarise(Freq = n()) %>%
  ungroup() %>%
  mutate(phylum = factor(phylum, levels = phylum_counts$Var1))
ggplot(phylum_counts_clear, aes(phylum, Freq, fill = Oxygen2)) +
  geom_bar(stat = "identity") +
  geom_text(data = phylum_counts,
            aes(x = Var1, y = Freq + 150, label = Freq), size = 2.5, angle = 45, 
            vjust = 1, inherit.aes = F) +
  labs(x = "Phylum", y = "# genomes w O2 info") +
  scale_fill_manual(values = c("red", "blue")) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 3000)) +
  ggtitle("BacDive May 2025", subtitle = "n bacteria with clear O2 and genomes = 7124") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(1,1),
        legend.justification = c(1,1))

# Merge with GTDB
sum(bacGT$gb_accession %in% d$Genome.seq..accession.number) # 7034 of 7124
sum(d$Genome.seq..accession.number %in% bacGT$gb_accession) # 7034 of 7124
missings <- d %>%
  filter(Genome.seq..accession.number %notin% bacGT$gb_accession)
# Missing ones don't appear to be a formatting issue. They just aren't there.
bacdive_gtdb <- bacGT %>%
  filter(gb_accession %in% d$Genome.seq..accession.number)
table(bacdive_gtdb$gtdb_representative)
bacdive_gtdb_rep <- bacdive_gtdb %>%
  filter(gtdb_representative == TRUE)

# Download these with list of NCBI accessions
# write.table(bacdive_gtdb_rep$ncbi_genbank_assembly_accession,
#             "taxon_ids_5521.txt", sep = "\t", row.names = F, col.names = F)
# perl -pi -e 's/\r\n/\n/g' taxon_ids_5521.txt
# sed 's/"//g' taxon_ids_5521.txt > clean.txt
# First pass downloaded 5520 of 5521. Which missed?
d_gtdb <- bacdive_gtdb_rep %>%
  left_join(., d, by = c("gb_accession" = "Genome.seq..accession.number")) %>%
  dplyr::select(accession, checkm2_completeness, checkm2_contamination,
                checkm_completeness, checkm_contamination, gc_percentage, genome_size,
                gtdb_representative, gtdb_taxonomy, ncbi_assembly_level, ncbi_country,
                ncbi_genbank_assembly_accession, ncbi_genome_category, ncbi_isolation_source, 
                ncbi_refseq_category, ncbi_species_taxid, ncbi_taxid, Oxygen.tolerance,
                Oxygen2)
dwn <- read.delim("data/genome_filenames.txt", header = F) %>%
  separate(V1, into = c("GenomeID", "Junk1", "Junk2", "Junk3"), sep = "_") %>%
  mutate(GenomeID = paste(GenomeID, Junk1, sep = "_"))
miss <- d_gtdb %>%
  filter(ncbi_genbank_assembly_accession %notin% dwn$GenomeID)
# This one GCA_005925325.1 is suppressed! So final number is 5520.
d_gtdb <- d_gtdb %>%
  filter(ncbi_genbank_assembly_accession %notin% miss$ncbi_genbank_assembly_accession) %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, remove = F, into = c("Domain", "Phylum", "Class", "Order", "Family",
                                               "Genus", "Species"), sep = ";")
#write.csv(d_gtdb, "bacdive_gtdb_metadata_5520.csv")
d_gtdb <- read.csv("data/bacdive_gtdb_metadata_5520.csv")
phylum_counts <- as.data.frame(table(d_gtdb$Phylum)) %>%
  arrange(desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = Var1))
phylum_counts_clear <- d_gtdb %>%
  group_by(Phylum, Oxygen2) %>%
  summarise(Freq = n()) %>%
  ungroup() %>%
  mutate(Phylum = factor(Phylum, levels = phylum_counts$Var1))
phy_fig <- ggplot(phylum_counts_clear, aes(Phylum, Freq, fill = Oxygen2)) +
  geom_bar(stat = "identity") +
  geom_text(data = phylum_counts,
            aes(x = Var1, y = Freq + 150, label = Freq), size = 2.5, angle = 45, 
            vjust = 1, inherit.aes = F) +
  labs(x = "Phylum", y = "# genomes w O2 info") +
  scale_fill_manual(values = c("red", "blue")) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 3000)) +
  ggtitle("BacDive May 2025", subtitle = "n GTDB reps with clear O2 and genomes = 5520") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(1,1),
        legend.justification = c(1,1))
phy_fig
pdf("InitialFigs/BacDive_BacterialPhylumNumbers_Binary_Reps.pdf", width = 8, height = 6)
phy_fig
dev.off()

phy_fig <- ggplot(phylum_counts_clear, aes(Phylum, Freq, fill = Oxygen2)) +
  geom_bar(stat = "identity") +
  geom_text(data = phylum_counts,
            aes(x = Var1, y = Freq + 150, label = Freq), size = 3, angle = 45, 
            vjust = 1, inherit.aes = F) +
  labs(x = "Phylum", y = "Genomes") +
  scale_fill_manual(values = c("red", "blue")) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 3000)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(1,1),
        legend.justification = c(1,1))
phy_fig
pdf("FinalFigs/FigureS1.pdf", width = 8, height = 6)
phy_fig
dev.off()
png("FinalFigs/FigureS1.png", width = 8, height = 6, units = "in", res = 300)
phy_fig
dev.off()



#### Get Info Pass 1 ####
# Loop through 112 dataframe with 100 rows each!
for (j in 112:length(split_list)) {
  info[[j]] <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 10)) %>%
    set_names(c("ID", "Oxygen", "NCBI_ID", "Domain", "Phylum", "Class", "Order", 
                "Family", "Genus", "Species"))
  for (i in 1:nrow(split_list[[j]])) {
    info[[j]]$ID[i] <- split_list[[j]]$ID[i]
    info[[j]]$Oxygen[i] <- split_list[[j]]$Oxygen.tolerance[i]
    f <- fetch(object = bacdive, ids = split_list[[j]]$ID[i])
    if ("Genome sequences" %in% names(f$results[[1]]$`Sequence information`)) {
      if (is.integer(f$results[[1]]$`Sequence information`$`Genome sequences`[[1]])) {
        info[[j]]$NCBI_ID[i] <- f$results[[1]]$`Sequence information`$`Genome sequences`[[length(f$results[[1]]$`Sequence information`$`Genome sequences`)]]
      } 
      if (is.list(f$results[[1]]$`Sequence information`$`Genome sequences`[[1]])) {
        if (!is.null(f$results[[1]]$`Sequence information`$`Genome sequences`[[1]]$`NCBI tax ID`)) {
          info[[j]]$NCBI_ID[i] <- f$results[[1]]$`Sequence information`$`Genome sequences`[[1]]$`NCBI tax ID`
        }
      }
    }
    if ("Genome sequences" %notin% names(f$results[[1]]$`Sequence information`)) {
      info[[j]]$NCBI_ID[i] <- "NoGenome"
    }
    if ("domain" %in% names(f$results[[1]]$`Name and taxonomic classification`)) {
      info[[j]]$Domain[i] <- f$results[[1]]$`Name and taxonomic classification`$domain
    }
    info[[j]]$Phylum[i] <- f$results[[1]]$`Name and taxonomic classification`$phylum
    info[[j]]$Class[i] <- f$results[[1]]$`Name and taxonomic classification`$class
    info[[j]]$Order[i] <- f$results[[1]]$`Name and taxonomic classification`$order
    info[[j]]$Family[i] <- f$results[[1]]$`Name and taxonomic classification`$family
    info[[j]]$Genus[i] <- f$results[[1]]$`Name and taxonomic classification`$genus
    info[[j]]$Species[i] <- f$results[[1]]$`Name and taxonomic classification`$species
  }
}

info_merged <- info[[1]]
for (i in 2:length(info)) {
  info_merged <- rbind(info_merged, info[[i]])
}
info_merged <- info_merged %>%
  filter(!is.na(ID))
info_merged_wGenome <- info_merged %>%
  filter(NCBI_ID != "NoGenome") # 8349

#saveRDS(info_merged_wGenome, "info_merged_wGenome.rds")
#write.table(info_merged_wGenome$NCBI_ID, "taxon_ids.txt", sep = "\t")



#### Info ####
info_merged_wGenome <- readRDS("info_merged_wGenome.rds")
table(info_merged_wGenome$Domain)
table(info_merged_wGenome$Phylum)
table(info_merged_wGenome$Oxygen)
phylum_counts <- as.data.frame(table(info_merged_wGenome$Phylum)) %>%
  arrange(desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = Var1))
pdf("BacDive_PhylumNumbers.pdf", width = 8, height = 6)
ggplot(phylum_counts, aes(Var1, Freq)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x = Var1, y = Freq + 150, label = Freq), size = 2.5, angle = 45, vjust = 1) +
  labs(x = "Phylum", y = "# genomes w O2 info") +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 3250)) +
  ggtitle("BacDive May 2025", subtitle = "n type strains with O2 and genomes = 8349") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
table(info_merged_wGenome$Class)
table(info_merged_wGenome$Order)
table(info_merged_wGenome$Family)
table(info_merged_wGenome$Genus)

# Remake with just the main 4 categories
# Note: This removes 1251 genomes
# And just bacteria
# Note: This removes 256 genomes
# Also fix BacDive phyla names
info_merged_wGenome_clear <- info_merged_wGenome %>%
  filter(Oxygen %in% c("obligate aerobe", "aerobe", "anaerobe", "obligate anaerobe")) %>%
  filter(Domain == "Bacteria") %>%
  mutate(Oxygen = factor(Oxygen,
                         levels = c("obligate aerobe", "aerobe", "anaerobe", "obligate anaerobe"))) %>%
  mutate(Oxygen2 = case_match(Oxygen,
                              c("obligate aerobe", "aerobe") ~ "aerobic",
                              c("anaerobe", "obligate anaerobe") ~ "anaerobic")) %>%
  mutate(Phylum2 = recode(Phylum,
                              "Proteobacteria" = "Pseudomonadota",
                              "Actinobacteria" = "Actinomycetota",
                              "Firmicutes" = "Bacillota",
                              "Bacteroidetes" = "Bacteroidota",
                              "Deinococcus-Thermus" = "Deinococcota",
                              "Thermotogae" = "Thermotogota",
                              "Spirochaetes" = "Spirochaetota",
                              "Chloroflexi" = "Chloroflexota",
                              "Planctomycetes" = "Planctomycetota",
                              "Acidobacteria" = "Acidobacteriota",
                              "Fusobacteria" = "Fusobacteriota",
                              "Not assigned to order" = "Unassigned",
                              "Verrucomicrobia" = "Verrucomicrobiota",
                              "Synergistetes" = "Synergistota",
                              "Aquificae" = "Aquificota",
                              "Thermodesulfobacteria" = "Desulfobacteriota",
                              "Chlorobi" = "Bacteroidota",
                              "Deferribacteres" = "Chrysiogenota",
                              "Balneolaeota" = "Bacteroidota",
                              "Tenericutes" = "Bacillota",
                              "Nitrospirae" = "Nitrospirota",
                              "Rhodothermaeota" = "Bacteroidota",
                              "Chrysiogenetes" = "Chrysiogenota",
                              "Dictyoglomi" = "Dictyoglomota",
                              "Fibrobacteres" = "Fibrobacterota",
                              "Lentisphaerae" = "Verrucomicrobiota",
                              "Armatimonadetes" = "Armatimonadota",
                              "Balneolota" = "Bacteroidota",
                              "Caldiserica" = "Caldisericota",
                              "Calditrichaeota" = "Electryoneota",
                              "Gemmatimonadetes" = "Gemmatimonadota",
                              "Thermomicrobia" = "Chloroflexota"))
phylum_counts <- as.data.frame(table(info_merged_wGenome_clear$Phylum2)) %>%
  arrange(desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = Var1))
phylum_counts_clear <- info_merged_wGenome_clear %>%
  group_by(Phylum2, Oxygen2) %>%
  summarise(Freq = n()) %>%
  ungroup() %>%
  mutate(Phylum2 = factor(Phylum2, levels = phylum_counts$Var1))
pdf("BacDive_BacterialPhylumNumbers_Binary.pdf", width = 8, height = 6)
ggplot(phylum_counts_clear, aes(Phylum2, Freq, fill = Oxygen2)) +
  geom_bar(stat = "identity") +
  geom_text(data = phylum_counts,
            aes(x = Var1, y = Freq + 150, label = Freq), size = 2.5, angle = 45, 
            vjust = 1, inherit.aes = F) +
  labs(x = "Phylum", y = "# genomes w O2 info") +
  scale_fill_manual(values = c("red", "blue")) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 3000)) +
  ggtitle("BacDive May 2025", subtitle = "n bacterial type strains with clear O2 and genomes = 6846") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(1,1),
        legend.justification = c(1,1))
dev.off()

length(unique(info_merged_wGenome_clear$NCBI_ID)) # 6824 of 6846
sum(is.na(info_merged_wGenome_clear$NCBI_ID))
taxid_count <- info_merged_wGenome_clear %>%
  group_by(NCBI_ID) %>%
  summarise(count = n())
sum(taxid_count$NCBI_ID %in% bacGT_useful$ncbi_species_taxid)
write.table(taxid_count$NCBI_ID, "taxon_ids.txt", sep = "\t", 
            row.names = F, col.names = F) # 6824 unique IDs
# Remember need to fix file too with 2 lines of code
#dos2unix taxon_ids.txt
#sed -i 's/"//g' taxon_ids.txt



#### Prevalence ####
# Start here
d <- read.csv("data/bacdive_gtdb_metadata_5520.csv") %>%
  mutate(Oxygen2 = as.factor(Oxygen2))
d_Aer <- d %>%
  filter(Oxygen2 == "aerobic")
d_An <- d %>%
  filter(Oxygen2 == "anaerobic")
pfams_merged_pa <- readRDS("data/pfams_merged_pa.rds")
pfams_aer <- pfams_merged_pa %>%
  select(all_of(d_Aer$ncbi_genbank_assembly_accession))
pfams_an <- pfams_merged_pa %>%
  select(all_of(d_An$ncbi_genbank_assembly_accession))

# Check prevalence of Pfams
pfams_prev <- data.frame(Pfam = rownames(pfams_merged_pa),
                         Prev = rowSums(pfams_merged_pa)) %>%
  mutate(Perc = Prev/5520*100) %>%
  filter(Prev > 0)
hist(pfams_prev$Prev)
sum(pfams_prev$Prev == 0) # 0 (duh)
sum(pfams_prev$Prev == 5520) # 73 (these need to be removed, not informative)
pfams_all <- pfams_prev %>%
  filter(Prev == 5520)

pfams_aer_prev <- data.frame(Pfam = rownames(pfams_aer),
                             Aer_Prev = rowSums(pfams_aer)) %>%
  mutate(Aer_Perc = Aer_Prev/4227*100) %>%
  filter(Aer_Prev > 0)
hist(pfams_aer_prev$Aer_Prev)
sum(pfams_aer_prev$Aer_Prev == 0) # 907 just in anaerobes
sum(pfams_aer_prev$Aer_Prev == 4227) # 90 in all
sum(pfams_aer_prev$Aer_Perc >= 80) # 1357 in ≥ 80%
pfams_aer_80 <- pfams_aer_prev %>%
  filter(Aer_Perc >= 80)
pfams_aer_50 <- pfams_aer_prev %>%
  filter(Aer_Perc >= 50)

pfams_an_prev <- data.frame(Pfam = rownames(pfams_an),
                            An_Prev = rowSums(pfams_an)) %>%
  mutate(An_Perc = An_Prev/1293*100) %>%
  filter(An_Prev > 0)
hist(pfams_an_prev$An_Prev)
sum(pfams_an_prev$An_Prev == 0) #  4182 just in aerobes
sum(pfams_an_prev$An_Prev == 1293) # 142 in all
sum(pfams_an_prev$An_Perc >= 80) # 1096 in ≥ 80%
pfams_an_80 <- pfams_an_prev %>%
  filter(An_Perc >= 80)
pfams_an_50 <- pfams_an_prev %>%
  filter(An_Perc >= 50)

# Elias model results top 20
pfam_info <- read.csv("data/selected_important_features.csv") %>%
  dplyr::select(Feature, Aggregated_Score) %>%
  set_names(c("Pfam", "SHAP")) %>%
  left_join(., pfams_aer_prev, by = "Pfam") %>%
  left_join(., pfams_an_prev, by = "Pfam") %>%
  mutate(Indicator = ifelse(Aer_Perc > 50, 
                            "Aerobic indicator", 
                            "Anaerobic indicator")) %>%
  mutate(PfamCode = substr(Pfam, 1, 7))

pdf("FinalFigs/Figure2.pdf", width = 7, height = 5)
ggplot(pfam_info, aes(x = Aer_Perc, y = An_Perc, colour = Indicator)) +
  geom_point(size = 3, pch = 16) + 
  geom_text_repel(aes(x = Aer_Perc, y = An_Perc, label = PfamCode),
                  inherit.aes = FALSE, size = 2) +
  scale_colour_manual(values = c("grey80", "black")) +
  labs(x = "Prevalence in aerobes (%)",
       y = "Prevalence in anaerobes (%)") +
  xlim(0, 100) +
  ylim(0, 100) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(1,1),
        legend.justification.inside = c(1,1),
        legend.background = element_blank())
dev.off()
png("FinalFigs/Figure2.png", width = 7, height = 5, units = "in", res = 300)
ggplot(pfam_info, aes(x = Aer_Perc, y = An_Perc, colour = Indicator)) +
  geom_point(size = 3, pch = 16) + 
  geom_text_repel(aes(x = Aer_Perc, y = An_Perc, label = PfamCode),
                  inherit.aes = FALSE, size = 2) +
  scale_colour_manual(values = c("grey80", "black")) +
  labs(x = "Prevalence in aerobes (%)",
       y = "Prevalence in anaerobes (%)") +
  xlim(0, 100) +
  ylim(0, 100) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(1,1),
        legend.justification.inside = c(1,1),
        legend.background = element_blank())
dev.off()



#### SHAP ####
# Plot the SHAP scores for the top 20 genes
# For supplementary figure
# pfam_info <- read.csv("data/selected_important_features.csv") %>%
#   dplyr::select(Feature, Aggregated_Score) %>%
#   set_names(c("Pfam", "SHAP")) %>%
#   left_join(., pfams_aer_prev, by = "Pfam") %>%
#   left_join(., pfams_an_prev, by = "Pfam") %>%
#   mutate(Indicator = ifelse(Aer_Perc > 50, 
#                             "Aerobic indicator", 
#                             "Anaerobic indicator")) %>%
#   mutate(PfamCode = substr(Pfam, 1, 7))
# pfam_info_sorted <- pfam_info %>%
#   arrange(Indicator, SHAP) %>%
#   mutate(Ind = gsub(" indicator", "", Indicator))
# saveRDS(pfam_info_sorted, "data/pfam_info_sorted.rds")
pfam_info_sorted <- readRDS("data/pfam_info_sorted.rds")
s <- read.csv("data/ensemble_shap_values_Class_1.csv")
# These are the 1104 genomes used for testing (20% of 5520)
table(s$true_class) # 845 0 (aerobe) and 259 1 (anaerobe)
# Note - the figure is colored by gene presence/absence NOT genome aerobe vs anaerobe!
# So, need to pull in and merge those data.
# Note - Careful with the index matching!
g <- read.csv("data/genes_testset.csv") %>%
  dplyr::select(-Oxygen2) %>%
  set_names(substr(names(.), 1, 7)) %>%
  pivot_longer(cols = 2:21) %>%
  mutate(Genome_Pfam = paste(X, name, sep = "_")) %>%
  rename(PA = value) %>%
  dplyr::select(Genome_Pfam, PA)
s <- read.csv("data/ensemble_shap_values_Class_1.csv") %>%
  pivot_longer(cols = c(4:23)) %>%
  mutate(true_class = as.factor(true_class),
         name = substr(name, 1, 7)) %>%
  left_join(., pfam_info_sorted, by = c("name" = "PfamCode")) %>%
  mutate(name = factor(name,
                       levels = pfam_info_sorted$PfamCode)) %>%
  mutate(Genome_Pfam = paste(X, name, sep = "_")) %>%
  left_join(., g, by = "Genome_Pfam") %>%
  mutate(PA = as.factor(PA))
figs2a <- ggplot(s, aes(name, value, colour = PA)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_beeswarm(method = "swarm", cex = 0.07, size = 1, priority = "density") +
  labs(x = NULL, y = "SHAP value", colour = NULL) +
  scale_colour_manual(values = c("blue", "red"),
                      labels = c("Absent", "Present")) +
  facet_grid(Ind ~ 1, space = "free", scales = "free_y") +
  coord_flip() +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.4),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
figs2a

# Make the confusion matrix for figs2b
TClass <- factor(c("Aerobe", "Aerobe", "Anaerobe", "Anaerobe"))
PClass <- factor(c("Aerobe", "Anaerobe", "Aerobe", "Anaerobe"))
Y      <- c("98%", "2%", "11%", "89%")
df <- data.frame(TClass, PClass, Y)
figs2b <- ggplot(df, aes(x = TClass, y = PClass)) +
  geom_tile(fill = "white", colour = "black", linewidth = 1) +
  geom_text(aes(label = Y), size = 7) +
  labs(x = "True", y = "Predicted") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14, vjust = -3),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, vjust = 5),
        axis.title.x = element_text(size = 14, vjust = 5),
        plot.margin = margin(0.1,0.1,0.1,0.1))
figs2b

pdf("FinalFigs/FigureS2.pdf", width = 8, height = 7)
plot_grid(figs2a, figs2b, labels = "auto", label_x = 0.2, rel_widths = c(0.55, 0.45))
dev.off()
png("FinalFigs/FigureS2.png", width = 8, height = 7, units = "in", res = 300)
plot_grid(figs2a, figs2b, labels = "auto", label_x = 0.2, rel_widths = c(0.55, 0.45))
dev.off()



#### Testing ####
# For supplementary table, make different phyla subsets to test and get metrics
# Elías already built the "common" phyla and "uncommon" phyla
# Here build an "all" phyla, and also each of the big 4 individually
oxygen_pfams <- read.csv("data/Oxygen_pfams.csv")
tax <- read.csv("data/bacdive_gtdb_metadata_5520.csv") %>%
  rename(Accession = ncbi_genbank_assembly_accession) %>%
  dplyr::select(Accession, Phylum) %>%
  rename(phylum = Phylum)
all <- read.csv("data/oxygen_pfams_data.csv") %>%
  rename("Accession" = "X") %>%
  left_join(., tax, by = "Accession") %>%
  #dplyr::select(Accession, phylum, all_of(oxygen_pfams$Pfam_full), Oxygen2)
  dplyr::select(Accession, phylum, everything()) %>%
  mutate(Oxygen2 = gsub("anaerobic", 1, Oxygen2)) %>%
  mutate(Oxygen2 = gsub("aerobic", 0, Oxygen2)) %>%
  mutate(Oxygen2 = as.integer(Oxygen2))
common <- all %>%
  filter(phylum %in% c("Pseudomonadota", "Actinomycetota", "Bacillota", "Bacteroidota"))
pseudo <- all %>%
  filter(phylum %in% c("Pseudomonadota"))
actino <- all %>%
  filter(phylum %in% c("Actinomycetota"))
bacill <- all %>%
  filter(phylum %in% c("Bacillota"))
bacter <- all %>%
  filter(phylum %in% c("Bacteroidota"))
uncommon <- all %>%
  filter(phylum %notin% c("Pseudomonadota", "Actinomycetota", "Bacillota", "Bacteroidota"))

# write.csv(all, "data/test_all.csv", row.names = FALSE)
# write.csv(common, "data/test_common.csv", row.names = FALSE)
# write.csv(pseudo, "data/test_pseudo.csv", row.names = FALSE)
# write.csv(actino, "data/test_actino.csv", row.names = FALSE)
# write.csv(bacill, "data/test_bacill.csv", row.names = FALSE)
# write.csv(bacter, "data/test_bacter.csv", row.names = FALSE)
# write.csv(uncommon, "data/test_uncommon.csv", row.names = FALSE)

# Make Table S5 from validation results
v <- read.csv("data/validation_results.csv") %>%
  group_by(phylum) %>%
  summarize(nT = sum(correct == TRUE),
            nF = sum(correct == FALSE)) %>%
  ungroup() %>%
  mutate(n = nT + nF) %>%
  mutate(Accuracy = round(nT/n*100, digits = 2)) %>%
  dplyr::select(phylum, Accuracy) %>%
  arrange(desc(Accuracy)) %>%
  rename(Phylum = phylum)
# write_xlsx(v, "~/Desktop/Fierer/AEGIS/Oxygen/Manuscript/TableS5.xlsx",
#            format_headers = FALSE)


#### End Script ####