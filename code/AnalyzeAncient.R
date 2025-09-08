# Calculate % aerobe for ancient metagenomes
# Use the model from the simulations

#### 1. Setup ####
# Libraries
library(tidyverse)
library(mgcv)
library(readxl)
library(car)

# Working directory
setwd("~/Desktop/Fierer/AEGIS/Oxygen")

# The 20 Pfams
oxygen_pfams <- read.csv("Oxygen_pfams.csv")
aerobic_pfams <- oxygen_pfams %>%
  filter(Oxygen == "aerobic")
anaerobic_pfams <- oxygen_pfams %>%
  filter(Oxygen == "anaerobic")

# Gene IDs for all of the variants of each Pfam
map <- read.table("pfam_headers_table_fixed.txt", sep = "\t", header = TRUE, 
                  stringsAsFactors = FALSE, quote = "") %>%
  separate(Header, into = c("Header", "Junk"), sep = " ") %>%
  select(-Junk)
table(map$Pfam) # Number of variants per Pfam

# Model (from simulations)
gam_model <- readRDS("gam_model.rds")



#### 2. Parse Diamond ####
# Loop through all the ancient samples
# First try bacterial + archaeal reads, then bacteria only
# Note: S3, S9, S11, S15, S16, S17, S19, S20 had no diamond output.
setwd("diamond_output_ancientBA/")
files <- list.files()
length(files)
results <- as.data.frame(matrix(nrow = length(files), ncol = 2, NA)) %>%
  set_names(c("filename", "ratio"))


for (i in 1:length(files)) {
  
  # Add filename to the dataframe
  results$filename[i] <- files[i]
  
  # Read in and filter the data
  d <- read.table(files[i]) %>%
    set_names(c("qseqid",	"sseqid",	"pident",	"length",	"qstart",	"qend",	
                "sstart",	"send",	"evalue",	"bitscore")) %>%
    filter(pident >= 60) %>%
    filter(evalue < 0.001) %>%
    filter(bitscore >= 50) %>%
    mutate(Gene.length = abs(send - sstart)) %>%
    left_join(., map, by = c("sseqid" = "Header"))
  
  # Gene hits
  gene.hits = d %>% 
    group_by(Pfam) %>% 
    summarise(total_count=n())
  
  # Correct for gene length (reads per kilobase)
  gene.hit.length.correction <- merge(mean.gene.lengths, gene.hits, by = "Pfam") %>%
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
setwd("../")
new_data <- data.frame(ratio = results$ratio)  # or a sequence of values
results$Per_aerobe <- predict(gam_model, newdata = new_data)
results <- results %>%
  mutate(Per_aerobe = ifelse(Per_aerobe > 100, 100, Per_aerobe)) %>%
  mutate(Per_aerobe = ifelse(Per_aerobe < 0, 0, Per_aerobe)) %>%
  mutate(sampleID = gsub("_R1_ann.tsv", "", filename))
saveRDS(results, "results_modern.rds")



#### 3. Analyze ####
# Now analyze the predicted % aerobe by habitat and depth
meta <- read_xlsx("AncientMetaG.xlsx")
length(unique(meta$sampleID)) # n = 278
results <- readRDS("results_modern.rds")
length(unique(results$sampleID)) # 278
read_counts <- read.delim("read_counts_table.tsv") %>%
  mutate(Per_Bact = Bact_Reads/QC_reads)
length(unique(read_counts$sampleID))
sum(read_counts$Bact_Reads > 1000000) # 109 with > 1 mil bact reads
d <- read_xlsx("ModernMetaG.xlsx") %>%
  left_join(., results, by = "sampleID") %>%
  left_join(., read_counts, by = "sampleID") %>%
  mutate(Depth_cm = as.numeric(Depth_cm),
         Oxygen = as.numeric(Oxygen),
         O2_H2S = as.numeric(O2_H2S),
         Temp = as.numeric(Temp),
         pH = as.numeric(pH),
         Habitat = as.factor(Habitat),
         Habitat2 = as.factor(Habitat2),
         Site = as.factor(Site)) %>%
  filter(is.na(Per_aerobe) == FALSE) %>%
  droplevels()
levels(d$Habitat)
levels(d$Habitat2)
range(d$Bact_Reads)
m <- aov(Bact_Reads ~ Habitat, data = d)
summary(m) # 8 1.617e+16 2.022e+15   23.49 <2e-16 ***
# From Paul Tol muted palette
pdf("InitialFigs/ReadCounts.pdf", width = 7, height = 5)
ggplot(d, aes(reorder(Habitat, Bact_Reads, mean), Bact_Reads, colour = Habitat)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 3, alpha = 0.8, width = 0.25) +
  scale_colour_manual(values = c("#882255", "brown1", "#882255", 
                                 "#88CCEE", 
                                 "#882255",
                                 "#DDCC77", "#44AA99", 
                                 #"#882255", 
                                 "#DDCC77", "#88CCEE", 
                                 "#44AA99", "#88CCEE", "#DDCC77", "black", "black")) +
  labs(x = NULL, y = "# Bacterial reads") +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12))
dev.off()



#### _Habitat ####
table(d$Habitat)
m <- aov(Per_aerobe ~ Habitat, data = d)
summary(m) # 11 254127   23102   57.94 <2e-16 ***
ggplot(d, aes(reorder(Habitat, Per_aerobe, mean), Per_aerobe, colour = Habitat)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 3, alpha = 0.8, width = 0.25) +
  scale_colour_manual(values = c("#882255", "brown1", "#88CCEE", "#882255", 
                                 "#88CCEE", "#882255", "#DDCC77", "#44AA99", 
                                 "#DDCC77", "#882255", "#DDCC77", "#88CCEE", 
                                 "#44AA99", "#DDCC77")) +
  labs(x = NULL, y = "Aerobe relative abundance (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12))

table(d$Habitat2)
m <- aov(Per_aerobe ~ Habitat2, data = d)
summary(m) # 16 300316   18770   97.64 <2e-16 ***
ggplot(d, aes(reorder(Habitat2, Per_aerobe, mean), Per_aerobe)) +
  geom_boxplot(outliers = F, aes(colour = Habitat2)) +
  geom_jitter(pch = 21, size = 3, alpha = 0.8, width = 0.25, aes(fill = Habitat2)) +
  scale_colour_manual(values = c("brown1", "#DDCC77", "#332288", "#332288", "#88CCEE", 
                                 "#88CCEE", "#DDCC77", "#44AA99", "#DDCC77", "#DDCC77", 
                                 "#332288", "#332288", "#44AA99", "#882255", "#882255", 
                                 "#882255","#882255", "#882255", "#882255", "#DDCC77")) +
  scale_fill_manual(values = c("brown1", "#DDCC77", "#332288", "#332288", "#88CCEE", 
                                 "#88CCEE", "#DDCC77", "#44AA99", "#DDCC77", "#DDCC77", 
                                 "#332288", "#332288", "#44AA99", "#882255", "#882255", 
                                 "#882255","#882255", "#882255", "#882255", "#DDCC77")) +
  labs(x = NULL, y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12))

d_plot <- d %>%
  filter(Plot == "Yes") %>%
  droplevels()
table(d_plot$Habitat3)
table(d_plot$Type)
m <- aov(Per_aerobe ~ Habitat3, data = d_plot)
summary(m) # 15 281245   18750   235.7 <2e-16 ***
pdf("InitialFigs/Habitat_Aerobe.pdf", width = 7, height = 5)
ggplot(d_plot, aes(reorder(Habitat3, Per_aerobe, mean), Per_aerobe)) +
  geom_boxplot(outliers = F, aes(colour = Type)) +
  geom_jitter(pch = 21, size = 3, alpha = 0.8, width = 0.25, aes(fill = Type)) +
  scale_colour_manual(values = c("#88CCEE", "#44AA99", "#332288", "#DDCC77", "brown1", 
                                 "#882255")) +
  scale_fill_manual(values = c("#88CCEE", "#44AA99", "#332288", "#DDCC77", "brown1", 
                               "#882255")) +
  labs(x = NULL, y = "Predicted aerobe relative abundance (%)") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 0),
        legend.justification.inside = c(1, 0),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12))
dev.off()



# End script