# Calculate % aerobe for ancient metagenomes
# Use the model from the simulations

#### 1. Setup ####
# Libraries
library(tidyverse)
library(mgcv)
library(readxl)
library(car)
library(segmented)
library(ggbeeswarm)
`%notin%` <- Negate(`%in%`)

# Working directory
setwd("~/Documents/GitHub/Oxygen/")

# The 20 Pfams
oxygen_pfams <- read.csv("data/Oxygen_pfams.csv")
aerobic_pfams <- oxygen_pfams %>%
  filter(Oxygen == "aerobic")
anaerobic_pfams <- oxygen_pfams %>%
  filter(Oxygen == "anaerobic")

# Gene IDs for all of the variants of each Pfam
map <- read.table("data/pfam_headers_table.txt", sep = "\t", header = TRUE, 
                  stringsAsFactors = FALSE, quote = "") %>%
  separate(Header, into = c("Header", "Junk"), sep = " ") %>%
  select(-Junk)
table(map$Pfam) # Number of variants per Pfam

# Mean gene lengths for the 20 Pfams
pfam_gene_length <- read.delim("data/pfam_mean_lengths.tsv") %>%
  separate(Pfam, sep = "-", into = c("Junk1", "Junk2", "Pfam")) %>%
  mutate(Gene.length = MeanLength * 3) %>%
  dplyr::select(Pfam, Gene.length)

# Model (from simulations)
gam_model <- readRDS("data/gam_model.rds")



#### 2. Parse Diamond ####
# Loop through all the ancient samples
# First try bacterial + archaeal reads, then bacteria only
# Note: S3, S9, S11, S15, S16, S17, S19, S20 had no Bact + Arch. (no S10)
# Note: S3, S9, S11, S12, S15, S16, S17, S19, S20 had no Bact.
# Note: Need to tweak the Diamond cutoffs for ancient DNA!
# Note: mean Bacterial read lengths range from 44 to 72, min 30, max 121 bp

#### _Bact + Arch ####
setwd("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_ancientBA/")
files <- list.files()
length(files)
results <- as.data.frame(matrix(nrow = length(files), ncol = 5, NA)) %>%
  set_names(c("filename", "ratio", "Pfams", "Aerobe_Pfams", "Anaerobe_Pfams"))


for (i in 1:length(files)) {
  
  # Add filename to the dataframe
  results$filename[i] <- files[i]
  
  # Read in and filter the data
  d <- read.table(files[i]) %>%
    set_names(c("qseqid",	"sseqid",	"pident",	"length",	"qstart",	"qend",	
                "sstart",	"send",	"evalue",	"bitscore")) %>%
    filter(pident >= 45) %>% # Changed from 60 to 45
    filter(evalue < 0.1) %>% # Changed from 0.001 to 0.1
    filter(bitscore >= 25) %>% # Changed from 50 to 25
    mutate(Gene.length = abs(send - sstart)) %>%
    left_join(., map, by = c("sseqid" = "Header"))
  
  # Number of the 20 Pfams found
  pf_count <- as.data.frame(table(d$Pfam)) 
  results$Pfams[i] <- nrow(pf_count)
  results$Aerobe_Pfams[i] <- sum(pf_count$Var1 %in% aerobic_pfams$Pfam)
  results$Anaerobe_Pfams[i] <- sum(pf_count$Var1 %in% anaerobic_pfams$Pfam)
  
  # Gene hits
  gene.hits = d %>% 
    group_by(Pfam) %>% 
    summarise(total_count=n())
  
  # Correct for gene length (reads per kilobase)
  # gene.hit.length.correction <- merge(mean.gene.lengths, gene.hits, by = "Pfam") %>%
  #   mutate(RPK = total_count / (1000*Gene.length)) %>%
  #   left_join(., oxygen_pfams, by = "Pfam")
  gene.hit.length.correction <- gene.hits %>%
    left_join(., pfam_gene_length, by = "Pfam") %>%
    mutate(RPK = total_count / (1000*Gene.length)) %>%
    left_join(., oxygen_pfams, by = "Pfam")
  
  # Now sum by aerobe indicator vs anaerobe indicator
  oxygen_rpk <- gene.hit.length.correction %>%
    group_by(Oxygen) %>%
    summarize(RPKsum = sum(RPK))
  
  # Calculate the ratio and add it to the dataframe
  results$ratio[i] <- oxygen_rpk$RPKsum[1] / oxygen_rpk$RPKsum[2]
  
  # Status
  message("Done processing table ", i)
  
}
warnings()
setwd("~/Documents/GitHub/Oxygen/")
new_data <- data.frame(ratio = results$ratio)
results$Per_aerobe <- predict(gam_model, newdata = new_data)
results <- results %>%
  mutate(Per_aerobe = ifelse(Per_aerobe > 100, 100, Per_aerobe)) %>%
  mutate(Per_aerobe = ifelse(Per_aerobe < 0, 0, Per_aerobe)) %>%
  mutate(sampleID = gsub("_R1_ann.tsv", "", filename))
#saveRDS(results, "data/results_ancientBA.rds")



#### _Bact ####
setwd("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_ancientB/")
files <- list.files()
length(files)
results <- as.data.frame(matrix(nrow = length(files), ncol = 5, NA)) %>%
  set_names(c("filename", "ratio", "Pfams", "Aerobe_Pfams", "Anaerobe_Pfams"))

for (i in 1:length(files)) {
  
  # Add filename to the dataframe
  results$filename[i] <- files[i]
  
  # Read in and filter the data
  d <- read.table(files[i]) %>%
    set_names(c("qseqid",	"sseqid",	"pident",	"length",	"qstart",	"qend",	
                "sstart",	"send",	"evalue",	"bitscore")) %>%
    filter(pident >= 45) %>% # Changed from 60 to 45
    filter(evalue < 0.1) %>% # Changed from 0.001 to 0.1
    filter(bitscore >= 25) %>% # Changed from 50 to 25
    mutate(Gene.length = abs(send - sstart)) %>%
    left_join(., map, by = c("sseqid" = "Header"))
  
  # Number of the 20 Pfams found
  pf_count <- as.data.frame(table(d$Pfam)) 
  results$Pfams[i] <- nrow(pf_count)
  results$Aerobe_Pfams[i] <- sum(pf_count$Var1 %in% aerobic_pfams$Pfam)
  results$Anaerobe_Pfams[i] <- sum(pf_count$Var1 %in% anaerobic_pfams$Pfam)
  
  # Gene hits
  gene.hits = d %>% 
    group_by(Pfam) %>% 
    summarise(total_count=n())
  
  # Correct for gene length (reads per kilobase)
  # gene.hit.length.correction <- merge(mean.gene.lengths, gene.hits, by = "Pfam") %>%
  #   mutate(RPK = total_count / (1000*Gene.length)) %>%
  #   left_join(., oxygen_pfams, by = "Pfam")
  gene.hit.length.correction <- gene.hits %>%
    left_join(., pfam_gene_length, by = "Pfam") %>%
    mutate(RPK = total_count / (1000*Gene.length)) %>%
    left_join(., oxygen_pfams, by = "Pfam")
  
  # Now sum by aerobe indicator vs anaerobe indicator
  oxygen_rpk <- gene.hit.length.correction %>%
    group_by(Oxygen) %>%
    summarize(RPKsum = sum(RPK))
  
  # Calculate the ratio and add it to the dataframe
  results$ratio[i] <- oxygen_rpk$RPKsum[1] / oxygen_rpk$RPKsum[2]
  
  # Status
  message("Done processing table ", i)
  
}
warnings()
setwd("~/Documents/GitHub/Oxygen/")
new_data <- data.frame(ratio = results$ratio)
results$Per_aerobe <- predict(gam_model, newdata = new_data)
results <- results %>%
  mutate(Per_aerobe = ifelse(Per_aerobe > 100, 100, Per_aerobe)) %>%
  mutate(Per_aerobe = ifelse(Per_aerobe < 0, 0, Per_aerobe)) %>%
  mutate(sampleID = gsub("_R1_ann.tsv", "", filename))
#saveRDS(results, "data/results_ancientB.rds")



#### 3. Analyze ####
# Now analyze the predicted % aerobe by habitat and depth
read_stats <- read.table("data/ancient_fastq_length_stats.tsv", 
                         sep = "\t",
                         header = F) %>%
  separate(V2, into = c("reads", "mean", "min", "max"), sep = " ") %>%
  separate(V1, into = c("sampleID", "filename"), sep = "_") %>%
  mutate(reads = as.numeric(reads),
         mean = as.numeric(mean),
         min = as.numeric(min),
         max = as.numeric(max))
meta1 <- read.table("data/per-sample-ancient-subglacial-only-nmds-loading.txt",
                    header = TRUE) %>%
  mutate(Cluster = ifelse(MDS1 < -0.5, "Left", "Right"))
table(meta1$Cluster)
ggplot(meta1, aes(MDS1, MDS2)) +
  geom_point() +
  theme_bw()
meta2 <- read_xlsx("data/subglacial_metadata_300325_AH_BDS.xlsx") %>%
  filter(!is.na(Filename))

resultsB <- readRDS("data/results_ancientB.rds") %>%
  separate(filename, into = c("sampleID", "filename"), sep = "_") %>%
  filter(sampleID %notin% c("left", "right")) %>%
  filter(Pfams >= 2) %>% # Have to have at least 2/20 Pfams detected
  left_join(., read_stats, by = "sampleID") %>%
  left_join(., meta2, by = c("sampleID" = "Filename")) %>%
  left_join(., meta1, by = c("Paper code" = "Sample")) %>%
  mutate(Method = ifelse(is.na(ratio), "Presence/Absence", "Ratio")) %>%
  dplyr::select(sampleID, Cluster, ratio, Pfams, Aerobe_Pfams, Anaerobe_Pfams,
                Per_aerobe, Method, reads, mean, min, max, Age, AgeBin, Composition,
                Lat, Lon, GeoCluster.x, δ15N, δ18OVPDB, δ13Ccarb, δ13Corg, `C:N`, 
                `%C`, `%N`, `87Sr/86Sr`, `P/Ca`, `234/238i`, `Color/Descrption`,
                `Collection Environment`, MDS1, MDS2) %>%
  mutate(Per_aerobe = ifelse(is.na(ratio), 
                             ifelse(Aerobe_Pfams == 0, 0, 100),
                             Per_aerobe))
table(resultsB$Cluster)

# Reads and Pfams
lm_fit <- lm(Pfams ~ reads, data = resultsB)
seg_fit <- segmented(lm_fit, seg.Z = ~reads, psi = list(reads = 30000))
seg_fit$psi
resultsB$fit <- fitted(seg_fit)
ggplot(resultsB, aes(reads, Pfams)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "# bacterial reads",
       y = "# Pfams detected") +
  scale_y_continuous(breaks = seq(1, 12)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
range(resultsB$reads)
pdf("InitialFigs/Bianca_ReadsPfams.pdf", width = 6, height = 4)
ggplot(resultsB, aes(mean, Pfams)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Mean bacterial read length (bp)",
       y = "# Pfams detected") +
  scale_y_continuous(breaks = seq(1, 12)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
dev.off()

# Cluster and % aerobes
pdf("InitialFigs/Bianca_Cluster_Aerobes.pdf", width = 6, height = 4)
ggplot(resultsB, aes(Cluster, Per_aerobe)) +
  geom_beeswarm(size = 4, pch = 21, color = "white", fill = "black") +
  labs(x = "NMDS Cluster",
       y = "Predicted % abundance of aerobic bacteria") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
dev.off()

ggplot(resultsB, aes(Cluster, Aerobe_Pfams)) +
  geom_beeswarm(size = 4, pch = 21, color = "white", fill = "black") +
  labs(x = "NMDS Cluster",
       y = "# Aerobic Pfams") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))

pdf("InitialFigs/Bianca_Cluster_AnaerobePfams.pdf", width = 6, height = 4)
ggplot(resultsB, aes(Cluster, Anaerobe_Pfams)) +
  geom_beeswarm(size = 4, pch = 21, color = "white", fill = "black") +
  labs(x = "NMDS Cluster",
       y = "# Anaerobic Pfams") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
dev.off()

# Other variables
ggplot(resultsB, aes(δ15N, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "δ15N",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(δ18OVPDB, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "δ18OVPDB",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(δ13Ccarb, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "δ13Ccarb",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(δ13Corg, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "δ13Corg",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(`C:N`, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "C:N",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(`%C`, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "%C",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(`%N`, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "%N",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(`87Sr/86Sr`, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "87Sr/86Sr",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(`P/Ca`, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "P/Ca",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(`234/238i`, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "234/238i",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
ggplot(resultsB, aes(Age, Per_aerobe, fill = Cluster)) +
  geom_point(size = 4, pch = 21, color = "white") +
  labs(x = "Age (kya)",
       y = "Predicted % abundance of aerobic bacteria") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12))
