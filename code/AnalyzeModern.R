# Calculate % aerobic bacteria in modern metagenomes
# by Cliff Bueno de Mesquita, Fierer Lab, 2025
# for AEGIS
# Use the GAM model from the simulated metagenomes

# Contents:
# 1. Setup
# 2. Parse Diamond
# 3. Analyze
# 4. Test Seq Depth
# 5. Application



#### 1. Setup ####
# Libraries
library(tidyverse)
library(mgcv)
library(readxl)
library(writexl)
library(car)
library(bestglm)
library(cowplot)
library(neonUtilities)
library(neonOS)
library(respirometry)
library(phyloNEON)
library(scales)
library(FSA)
library(ggside)
library(Z10)
library(httr)
library(jsonlite)
library(broom)
library(gamlss)
library(mctoolsr)
library(tidytext)
library(ozmaps)
library(sf)
library(ggspatial)
library(urbnmapr)
library(khroma) # According to Paul Tol’s technical note, the bright, contrast, vibrant and muted color schemes are color-blind safe.
highcontrast <- color("high contrast")
highcontrast(3)[1:3]

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
  select(-Junk) %>%
  filter(!duplicated(Header))
table(map$Pfam) # Number of variants per Pfam

# Mean gene lengths for the 20 Pfams
pfam_gene_length <- read.delim("data/pfam_mean_lengths.tsv") %>%
  separate(Pfam, sep = "-", into = c("Junk1", "Junk2", "Pfam")) %>%
  mutate(Gene.length = MeanLength * 3) %>%
  dplyr::select(Pfam, Gene.length)
# write.table(pfam_gene_length, 
#             "data/pfam_lengths.tsv", 
#             sep = "\t",
#             row.names = F)

# Note: RPK = (Number of reads mapped to the gene) / (Gene length in base pairs / 1000)

# Model (from simulations)
gam_model <- readRDS("data/gam_model.rds")

# Tran et al. 2021 Oxygen
# (metaG don't have O2, but can infer from depth)
# tran <- read.csv("data/TranDO.csv")
# ggplot(tran, aes(Depth, DO)) +
#   geom_point() +
#   geom_smooth() +
#   geom_vline(xintercept = 122, color = "red") +
#   geom_hline(yintercept = 12, color = "red") +
#   geom_vline(xintercept = 1, color = "red") +
#   geom_hline(yintercept = 84, color = "red") +
#   labs(x = "Depth (m)",
#        y = "DO") +
#   theme_bw()
# tran_model <- gam(DO ~ s(Depth), data = tran)
# summary(tran_model)
# plot(tran_model)
# formula(tran_model)
# new_tran_data <- data.frame(Depth = c(0, 120, 150, 200, 250, 300, 400, 1200))
# predicted_tran_o2 <- predict(tran_model, newdata = new_tran_data)
# predicted_tran_o2
# Add these values as the oxygen numbers for Freshwater Lake in ModernMetaG.xlsx
# But use 1.3 for the minimum instead of negative numbers



#### 2. Parse Diamond ####
#### _Test ####
# Trial run on the first table 
d <- read.table("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_modern/19485_R1_ann.tsv") %>%
  set_names(c("qseqid",	"sseqid",	"pident",	"length",	"qstart",	"qend",	
              "sstart",	"send",	"evalue",	"bitscore")) %>%
  filter(pident >= 60) %>%
  filter(evalue < 0.001) %>%
  filter(bitscore >= 50) %>%
  left_join(., map, by = c("sseqid" = "Header"))
table(d$Pfam) # PF05425 (aerobic) and PF13597 (anaerobic) are missing!

gene.hits = d %>% 
  group_by(Pfam) %>% 
  summarise(total_count=n())

gene.hit.length.correction <- gene.hits %>%
  left_join(., pfam_gene_length, by = "Pfam") %>%
  mutate(RPK = total_count / (Gene.length/1000)) %>%
  left_join(., oxygen_pfams, by = "Pfam")

oxygen_rpk <- gene.hit.length.correction %>%
  group_by(Oxygen) %>%
  summarize(RPKsum = sum(RPK))

ratio <- oxygen_rpk$RPKsum[1] / oxygen_rpk$RPKsum[2]

new_data <- data.frame(ratio = ratio)  # or a sequence of values
predicted_y <- predict(gam_model, newdata = new_data)
predicted_y # 81.93% aerobe (oldold RPK was 79.56%, old map was 82.32%, old RPK was 81.74%)

# Test more individual samples
# Need to diagnose outliers and unexpected results
# One question is how many of the 20 genes are present in certain samples

# Hot Spring Sediment
t <- read.table("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_modern/SRR25522566_R1_ann.tsv") %>%
  set_names(c("qseqid",	"sseqid",	"pident",	"length",	"qstart",	"qend",	
              "sstart",	"send",	"evalue",	"bitscore")) %>%
  filter(pident >= 60) %>%
  filter(evalue < 0.001) %>%
  filter(bitscore >= 50) %>%
  left_join(., map, by = c("sseqid" = "Header"))
t_df <- as.data.frame(table(t$Pfam)) 
nrow(t_df) # 15/20 genes present
sum(t_df$Var1 %in% aerobic_pfams$Pfam) # 10 aerobic pfams
sum(t_df$Var1 %in% anaerobic_pfams$Pfam) # 5 anaerobic pfams



#### _Modern ####
# Now loop through all the modern samples
setwd("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_modern/")
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
    filter(pident >= 60) %>%
    filter(evalue < 0.001) %>%
    filter(bitscore >= 50) %>%
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
  
  gene.hit.length.correction <- gene.hits %>%
    left_join(., pfam_gene_length, by = "Pfam") %>%
    mutate(RPK = total_count / (Gene.length/1000)) %>%
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
  mutate(Per_aerobe = ifelse(ratio > 35, 100, Per_aerobe)) %>%
  mutate(sampleID = gsub("_R1_ann.tsv", "", filename))
#saveRDS(results, "data/results_modern.rds")



#### 3. Analyze ####
# Now analyze the predicted % aerobe by habitat and depth
meta <- read_xlsx("data/ModernMetaG.xlsx")
length(unique(meta$sampleID)) # n = 278
results <- readRDS("data/results_modern.rds")
hist(results$ratio)
sum(results$ratio > 35)
length(unique(results$sampleID)) # 278
read_counts <- read.delim("data/read_counts_table.tsv") %>%
  mutate(Per_Bact = Bact_Reads/QC_reads)
length(unique(read_counts$sampleID)) # 278
range(read_counts$Bact_Reads) # 28616 to 99898651
mean(read_counts$Bact_Reads) # 9696925
sd(read_counts$Bact_Reads) # 13165046
sum(read_counts$Bact_Reads > 1000000) # 215 with > 1 mil bact reads
d <- read_xlsx("data/ModernMetaG.xlsx") %>%
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



#### _Habitat ####
# Tried different versions. Use this. Habitat3 variable.
d_plot <- d %>%
  filter(Plot == "Yes") %>%
  droplevels()
TableS1 <- d %>%
  filter(Plot == "Yes") %>%
  droplevels() %>%
  dplyr::select(-pH, -Temp, -filename, -Study, -Plot, -Habitat, -Habitat2) %>%
  rename(Habitat = Habitat3)
# writexl::write_xlsx(TableS1, 
#                     "~/Desktop/Fierer/AEGIS/Oxygen/Manuscript/TableS1.xlsx",
#                     format_headers = FALSE)
  

# Read stats for this subset
range(d_plot$Bact_Reads) # 28616 to 99898651
mean(d_plot$Bact_Reads) # 10919808
sd(d_plot$Bact_Reads) # 14609572

table(d_plot$Habitat3)
table(d_plot$Type)
m <- aov(Per_aerobe ~ Habitat3, data = d_plot)
summary(m) # 15 302114   20141   212.7 <2e-16 ***
ggplot(d_plot, aes(Habitat, Per_aerobe)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
ggplot(d_plot, aes(Habitat2, Per_aerobe)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
fig3 <- ggplot(d_plot, aes(reorder(Habitat3, Per_aerobe, mean), Per_aerobe)) +
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
fig3
pdf("FinalFigs/Figure3.pdf", width = 7, height = 5)
fig3
dev.off()
png("FinalFigs/Figure3.png", width = 7, height = 5, units = "in", res = 300)
fig3
dev.off()

# Reads
hist(d$QC_reads)
hist(d$Bact_Reads)
hist(d$Per_Bact)
m <- aov(Bact_Reads ~ Habitat3, data = d_plot)
summary(m) # 15 2.852e+16 1.901e+15   24.36 <2e-16 ***
pdf("InitialFigs/Habitat_Bact_Reads.pdf", width = 7, height = 5)
ggplot(d_plot, aes(reorder(Habitat3, Bact_Reads, mean), Bact_Reads)) +
  geom_boxplot(outliers = F, aes(colour = Type)) +
  geom_jitter(pch = 21, size = 3, alpha = 0.8, width = 0.25, aes(fill = Type)) +
  scale_colour_manual(values = c("#88CCEE", "#44AA99", "#332288", "#DDCC77", "brown1", 
                                 "#882255")) +
  scale_fill_manual(values = c("#88CCEE", "#44AA99", "#332288", "#DDCC77", "brown1", 
                               "#882255")) +
  scale_y_log10() +
  labs(x = NULL, 
       y = "Bacterial reads") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 0),
        legend.justification.inside = c(1, 0),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        plot.margin = margin(5,5,5,20))
dev.off()



#### _Black Sea ####
bs <- d %>%
  filter(Study == "CaseStudy1") %>%
  droplevels()
table(bs$Habitat)
m <- lm(Per_aerobe ~ Oxygen, data = bs)
summary(m) # p = 0.16
ggplot(bs, aes(Oxygen, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Oxygen", y = "Aerobe relative abundance (%)") +
  ylim(0, 100) +
  scale_x_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

m <- lm(Per_aerobe ~ poly(Depth_cm, 3, raw = T), data = bs)
summary(m) # p = 0.006
ggplot(bs, aes(Depth_cm, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3, raw = TRUE), color = "blue") +
  labs(x = "Depth (km)", y = "Aerobe relative abundance (%)") +
  #ylim(0, 100) +
  scale_x_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

m <- lm(Per_aerobe ~ log(O2_H2S), data = bs)
summary(m) # p < 0.001
m <- lm(Per_aerobe ~ poly(O2_H2S, 3, raw = T), data = bs)
summary(m) # p < 0.001
pdf("InitialFigs/BlackSea_Aerobe_O2_H2S.pdf", width = 7, height = 5)
ggplot(bs, aes(O2_H2S, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(se = F) +
  labs(x = "[O2] / [H2S]", y = "Predicted aerobe relative abundance (%)") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(-5, 110)) +
  scale_x_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()

range(bs$O2_H2S)
bs_gam <- gam(Per_aerobe ~ s(O2_H2S), data = bs, method = "REML", bs = "cs")
summary(bs_gam)
plot(bs_gam)
formula(bs_gam)
new_data <- data.frame(O2_H2S = seq(min(bs$O2_H2S), max(bs$O2_H2S, length.out = 100)))
new_data$Per_aerobe <- predict(bs_gam, newdata = new_data)

p1 <- ggplot(bs, aes(O2_H2S, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(se = F) +
  labs(x = expression("[" * O[2] * "] / [" * H[2] * S * "] ratio"),
       y = "Predicted aerobes (%)") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(-5, 110)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")) +
  ggtitle("a) Black Sea water") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
p1



#### _Baltic ####
# Baltic Sea sediments, oxygen gradient
baltic <- d %>%
  filter(Study == "CaseStudy5") %>%
  droplevels()
table(baltic$Habitat)
m <- lm(Per_aerobe ~ Oxygen, data = baltic)
summary(m) # p = 0.01, R2 = 0.38
pdf("InitialFigs/Baltic_Aerobe_O2.pdf", width = 7, height = 5)
ggplot(baltic, aes(Oxygen, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "Dissolved oxygen (mg/L)", 
       y = "Predicted aerobe relative abundance (%)") +
  ylim(0, 100) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()

p2 <- ggplot(baltic, aes(Oxygen, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Dissolved oxygen (mg/L)", 
       y = "Predicted aerobes (%)") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(-5, 110)) +
  ggtitle("b) Baltic Sea sediment") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
p2



#### _Freshwater Lake ####
# Have aerobic and anaerobic Lake Tanganyika samples
lake <- d %>%
  filter(Habitat == "Freshwater lake") %>%
  droplevels()
table(lake$Habitat)
table(lake$Habitat2)
m <- t.test(Per_aerobe ~ Habitat2, data = lake)
m # p = 0.0005
pdf("InitialFigs/Tanganyika_Aerobe_Depth.pdf", width = 7, height = 5)
ggplot(lake, aes(Depth_cat, Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 3, alpha = 0.75, width = 0.25, height = 0) +
  labs(x = "Depth", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12))
dev.off()

lab <- data.frame(x = c("< 121 m", "> 149 m"),
                  y = c(107, 85),
                  label = c("a", "b"))
p3 <- ggplot(lake, aes(Depth_cat, Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 3, alpha = 0.75, width = 0.25, height = 0) +
  geom_text(data = lab, aes(x = x, y = y, label = label),
            inherit.aes = F, size = 5) +
  labs(x = "Depth", y = "Predicted aerobes (%)") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(-5, 110)) +
  ggtitle("c) Lake Tanganyika water") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
p3

ggplot(lake, aes(Depth_cm/100, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "solid") +
  labs(x = "Depth (m)", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggplot(lake, aes(Oxygen, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "Oxygen (mg/L)", y = "Predicted aerobes (%)") +
  ylim(0, 100) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

fig2 <- plot_grid(p1, p2, p3, ncol = 3, align = "h",
                  rel_widths = c(0.35, 0.3, 0.3))
pdf("FinalFigs/Figure2.pdf", width = 8, height = 4)
fig2
dev.off()
png("FinalFigs/Figure2.png", width = 8, height = 4, units = "in", res = 300)
fig2
dev.off()



#### 4. Seq Depth ####
# Want to investigate the effect of sequencing depth
# Choose 1 high read depth clearly aerobic sample to downsample
# Choose 1 high read depth clearly anaerobic sample to downsample
# Choose 1 high read depth with ~50% aerobes to downsample
# SRR7053051 (aerobic lake) has 52,283,874 and 100% aerobe
# SRR7008372 (deeper lake) has 32,901,414 and 61% aerobe
# SRR10680549 (human gut) has 20,981,616 and 0% aerobe
# Plot reads versus predicted % aerobes and total pfams
# Should be stable for most seq depths except for a few thousand reads
# Do 1k 2k 4k 8k reads etc




#### _Parse ####
# Loop through the seqdepth diamond outputs
setwd("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_seqdepth/")
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
    filter(pident >= 60) %>%
    filter(evalue < 0.001) %>%
    filter(bitscore >= 50) %>%
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
  gene.hit.length.correction <- gene.hits %>%
    left_join(., pfam_gene_length, by = "Pfam") %>%
    mutate(RPK = total_count / (Gene.length/1000)) %>%
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
  mutate(Per_aerobe = ifelse(ratio > 35, 100, Per_aerobe)) %>%
  mutate(sampleID = gsub("_R1_ann.tsv", "", filename))
#saveRDS(results, "data/results_seqdepth.rds")



#### _Analyze ####
seqdepth <- readRDS("data/results_seqdepth.rds") %>%
  separate(filename, into = c("sampleID", "SeqDepth", "Junk")) %>%
  select(-Junk) %>%
  mutate(SeqDepth = as.numeric(SeqDepth)) %>%
  mutate(Habitat = gsub("SRR10680549", "Gut", sampleID)) %>%
  mutate(Habitat = gsub("SRR7053051", "Surface Lake", Habitat)) %>%
  mutate(Habitat = gsub("SRR7008372", "Deep Lake", Habitat)) %>%
  mutate(Habitat = factor(Habitat, levels = c("Surface Lake", "Deep Lake", "Gut"))) %>%
  arrange(Habitat, SeqDepth)
seqdepth2 <- read.table("~/Desktop/per_aerobe_predictions.tsv", header = T) %>%
  separate(SampleID, into = c("sampleID", "SeqDepth"), sep = "_") %>%
  mutate(SeqDepth = as.numeric(SeqDepth)) %>%
  mutate(Habitat = gsub("SRR10680549", "Gut", sampleID)) %>%
  mutate(Habitat = gsub("SRR7053051", "Surface Lake", Habitat)) %>%
  mutate(Habitat = gsub("SRR7008372", "Deep Lake", Habitat)) %>%
  mutate(Habitat = factor(Habitat, levels = c("Surface Lake", "Deep Lake", "Gut"))) %>%
  arrange(Habitat, SeqDepth)
plot(seqdepth$ratio, seqdepth2$ratio)
  
pdf("InitialFigs/SeqDepth_Aerobes.pdf", width = 8, height = 4)
ggplot(seqdepth, aes(SeqDepth, Per_aerobe, colour = Habitat)) +
  geom_point(size = 4) +
  geom_line() +
  scale_x_continuous(trans = "log2",
                     breaks = c(1000, 2000, 4000, 8000, 16000, 32000, 64000, 
                                128000, 256000, 512000, 1024000, 2048000, 4096000, 
                                8192000, 16384000, 32768000),
                     labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Reads",
       y = "Predicted % aerobes") +
  scale_colour_manual(values = highcontrast(3)[1:3]) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.3))
dev.off()

pdf("InitialFigs/SeqDepth_Ratio.pdf", width = 8, height = 4)
ggplot(seqdepth, aes(SeqDepth, ratio, colour = Habitat)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(trans = "log2",
                     breaks = c(1000, 2000, 4000, 8000, 16000, 32000, 64000, 
                                128000, 256000, 512000, 1024000, 2048000, 4096000, 
                                8192000, 16384000, 32768000),
                     labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Reads",
       y = "Ratio") +
  scale_colour_manual(values = highcontrast(3)[1:3]) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.3))
dev.off()

pdf("InitialFigs/SeqDepth_Pfams.pdf", width = 8, height = 4)
ggplot(seqdepth, aes(SeqDepth, Pfams, colour = Habitat)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(trans = "log2",
                     breaks = c(1000, 2000, 4000, 8000, 16000, 32000, 64000, 
                                128000, 256000, 512000, 1024000, 2048000, 4096000, 
                                8192000, 16384000, 32768000),
                     labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Reads",
       y = "Pfams") +
  scale_colour_manual(values = highcontrast(3)[1:3]) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.3))
dev.off()

pdf("InitialFigs/SeqDepth_AerobePfams.pdf", width = 8, height = 4)
ggplot(seqdepth, aes(SeqDepth, Aerobe_Pfams, colour = Habitat)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(trans = "log2",
                     breaks = c(1000, 2000, 4000, 8000, 16000, 32000, 64000, 
                                128000, 256000, 512000, 1024000, 2048000, 4096000, 
                                8192000, 16384000, 32768000),
                     labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Reads",
       y = "Aerobic Pfams") +
  scale_colour_manual(values = highcontrast(3)[1:3]) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.3))
dev.off()

pdf("InitialFigs/SeqDepth_AnaerobePfams.pdf", width = 8, height = 4)
ggplot(seqdepth, aes(SeqDepth, Anaerobe_Pfams, colour = Habitat)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(trans = "log2",
                     breaks = c(1000, 2000, 4000, 8000, 16000, 32000, 64000, 
                                128000, 256000, 512000, 1024000, 2048000, 4096000, 
                                8192000, 16384000, 32768000),
                     labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Reads",
       y = "Anaerobic Pfams") +
  scale_colour_manual(values = highcontrast(3)[1:3]) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.3))
dev.off()



#### 5. Application ####
# Use the tool on a big dataset to answer a question
# Which variables (precip, aridity, moisture, carbon, texture) are most closely
# associated with predicted % aerobic bacteria relative abundance?
# And - is this consistent across 3 independent datasets?
# Test across Australia, Panama, NEON metagenomes and metadata



#### _Australia ####
# Analyze all 331 natural surface Australian metagenomes

# Loop through the diamond output for the 331 samples
setwd("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_soil/")
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
    filter(pident >= 60) %>%
    filter(evalue < 0.001) %>%
    filter(bitscore >= 50) %>%
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
  
  gene.hit.length.correction <- gene.hits %>%
    left_join(., pfam_gene_length, by = "Pfam") %>%
    mutate(RPK = total_count / (Gene.length/1000)) %>%
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
  mutate(Per_aerobe = ifelse(ratio > 35, 100, Per_aerobe)) %>%
  mutate(sampleID = gsub("_R1_ann.tsv", "", filename)) %>%
  separate(filename, into = c("sampleID", "Junk1", "Junk2", "Junk3"),
           sep = "_") %>%
  dplyr::select(-Junk1, -Junk2, -Junk3)
#saveRDS(results, "data/results_soil.rds")

# Analyze
results <- readRDS("data/results_soil.rds")
range(results$Per_aerobe) # 64% to 100%
hist(results$Per_aerobe)
mean(results$Per_aerobe) # 97.16
median(results$Per_aerobe) # 100
sd(results$Per_aerobe) # 6.26
sum(results$Per_aerobe == 100) # 243 with 100, 88 with anaerobes
aus <- read_delim("~/Documents/GitHub/AussieStrains/data/metadata_331.txt") %>%
  mutate(sampleID = as.character(sampleID)) %>%
  left_join(., results, by = "sampleID") %>%
  mutate("ClimateClass" = ifelse(AI < 0.5, "Arid to semi-arid",
                                  ifelse(AI >= 0.5 & AI < 0.65, "Dry sub-humid",
                                         "Humid"))) %>%
  separate(utc_date_sampled, into = c("Month", "Day", "Year"), sep = "\\/") %>%
  mutate(Month = as.character(Month)) %>%
  mutate(Month = na_if(Month, "NA")) %>%
  mutate(Month = factor(Month,
                        levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                   "9", "10", "11", "12")))

# Get relevant environmental data with no NA
input <- readRDS("~/Documents/GitHub/AussieStrains/data/input_sylph.rds")
# Check environmental variables
env_vars <- c("boron_hot_cacl2", "clay", "conductivity",
              "dtpa_copper", "dtpa_iron", "dtpa_manganese", "dtpa_zinc",
              "elev",
              "exc_aluminium", "exc_calcium", "exc_magnesium",
              "exc_potassium", "exc_sodium", "latitude",
              "longitude", "nitrate_nitrogen", "organic_carbon",
              "ph", "phosphorus_colwell",
              "sand", "silt", "sulphur", 
              "water_content", # 51 NA
              "bio1", "bio12", "AI")
d_env <- input$map_loaded %>%
  dplyr::select(sampleID, all_of(env_vars))
n_na <- c()
for (i in 1:ncol(d_env)) {
  n_na[i] <- sum(is.na(d_env[,i]))
}
n_na # Most NA is 67, which means 201 have, which is pretty good.
# envfit can handle NA anyway
# remove boron_hot_cacl2, exc_sodium, sulphur, bio12 because too correlated
# Also sand and iron
# For this though keep bio12 because it's relevant
d_env <- input$map_loaded %>%
  dplyr::select(clay, conductivity, dtpa_copper, dtpa_manganese, 
                dtpa_zinc, elev, exc_aluminium, exc_calcium, exc_magnesium, exc_potassium,
                latitude, longitude, nitrate_nitrogen, organic_carbon, ph, 
                phosphorus_colwell, silt, water_content, bio1, bio12, AI)
d_env2 <- d_env %>%
  dplyr::select(-latitude, -longitude) %>%
  scale() %>%
  as.data.frame() %>%
  set_names(c("Clay", "Conductivity", "Cu", "Mn", "Zn", "Elevation", "Al", "Ca", "Mg", "K",
              "NO3", "C", "pH", "P", "Silt", "H2O", "Temp.", "Precip.", "Aridity")) %>%
  mutate(ClaySilt = Clay + Silt)
X <- d_env2 %>%
  dplyr::select(ClaySilt, Conductivity, NO3, C, pH, P, H2O, Elevation,
                Temp., Precip., Aridity) %>%
  filter_at(vars(1:11), all_vars(!is.na(.)))
y <- aus %>%
  filter(sampleID %in% rownames(X)) %>%
  dplyr::select(Per_aerobe)
Xy <- cbind(X, y) # n = 152
bestM <- bestglm(Xy,
                 family = gaussian,
                 IC = "AIC",
                 t = "default",
                 CVArgs = "default", 
                 qLevel = 0.99, 
                 TopModels = 10, 
                 method = "exhaustive", 
                 intercept = TRUE, 
                 weights = NULL, 
                 nvmax = "default", 
                 RequireFullEnumerationQ = FALSE)
bestM # P, Elevation, Temp.
# Estimate Std. Error    t value      Pr(>|t|)
# (Intercept) 97.537932  0.4848727 201.161946 2.366554e-182
# P            1.615514  0.9565573   1.688883  9.334754e-02
# Elevation    1.662556  0.7100120   2.341589  2.053527e-02
# Temp.        1.603956  0.6270567   2.557912  1.153551e-02
bestM$BestModels

# Linear Models
aus_uni <- aus %>%
  mutate(ClaySilt = clay + silt) %>%
  dplyr::select(Per_aerobe, ClaySilt, conductivity, nitrate_nitrogen, organic_carbon, 
                ph, phosphorus_colwell, water_content, elev, bio1, bio12,
                AI) %>%
  set_names(c("Per_aerobe", "ClaySilt", "Conductivity", "NO3", "C", "pH", "P", 
              "H2O", "Elevation", "Temp.", "Precip.", "Aridity")) %>%
  dplyr::select(ClaySilt, Conductivity, NO3, C, pH, P, H2O, Elevation, Temp., 
                Precip., Aridity, Per_aerobe)
models <- as.data.frame(matrix(NA, 11, 3)) %>%
  set_names("Variable", "R", "p")
for (i in 1:11) {
  models$Variable[i] <- names(aus_uni)[i]
  m <- lm(Per_aerobe ~ aus_uni[[i]], data = aus_uni)
  s <- summary(m)
  models$R[i] <- s$r.squared
  models$p[i] <- s$coefficients[2,4]
}
# ClaySilt, H2O, and C are significant

# Univariate
# Continuous
pdf("InitialFigs/Aus_Aerobe_Texture.pdf", width = 7, height = 5)
ggplot(aus_uni, aes(ClaySilt, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "% Clay + Silt",
       y = "Predicted % aerobic bacteria") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()
ggplot(aus, aes(water_content, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "% water content",
       y = "Predicted % aerobic bacteria") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggplot(aus, aes(organic_carbon, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "% C",
       y = "Predicted % aerobic bacteria") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

ggplot(aus, aes(AI, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "drier <==  Aridity index  ==> wetter",
       y = "Predicted % aerobic bacteria") +
  theme_classic()
ggplot(aus, aes(bio1, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth() +
  labs(x = "MAT",
       y = "Predicted % aerobic bacteria") +
  theme_classic()
ggplot(aus, aes(bio12, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth() +
  labs(x = "MAP",
       y = "Predicted % aerobic bacteria") +
  theme_classic()
ggplot(aus, aes(nitrate_nitrogen, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth() +
  labs(x = "N",
       y = "Predicted % aerobic bacteria") +
  theme_classic()
ggplot(aus, aes(ph, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth() +
  theme_classic()
ggplot(aus, aes(phosphorus_colwell, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth() +
  theme_classic()
ggplot(aus, aes(conductivity, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth() +
  labs(x = "Conductivity",
       y = "Predicted % aerobic bacteria") +
  theme_classic()
ggplot(aus, aes(phosphorus_colwell, Per_aerobe)) +
  geom_point(alpha = 0.75) +
  geom_smooth() +
  theme_classic()

# Categorical
m <- aov(Per_aerobe ~ vegetation_type, data = aus)
summary(m)
TukeyHSD(m)
ggplot(aus, aes(reorder(vegetation_type, Per_aerobe, median), Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_jitter(alpha = 0.75, width = 0.25) +
  labs(x = "Vegetation",
       y = "Predicted % aerobic bacteria") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
m <- aov(Per_aerobe ~ ClimateClass, data = aus)
summary(m)
TukeyHSD(m)
ggplot(aus, aes(reorder(ClimateClass, Per_aerobe, median), Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_jitter(alpha = 0.75, width = 0.25) +
  labs(x = "Aridity class",
       y = "Predicted % aerobic bacteria") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
m <- aov(Per_aerobe ~ Month, data = aus)
summary(m)
TukeyHSD(m)
ggplot(aus, aes(Month, Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_jitter(alpha = 0.75, width = 0.25) +
  labs(x = "Month",
       y = "Predicted % aerobic bacteria") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# What about interactions or multivariate
m <- lm(Per_aerobe ~ vegetation_type*sand, data = aus)
summary(m)
m <- lm(Per_aerobe ~ sand + water_content + organic_carbon, data = aus)
summary(m)

# What about correlations with known anaerobic taxa (e.g. methanogens)
# Methanogens aren't in the model, so independent test
# Load in mTAGs data from Bradyrhizobium manuscript GitHub repo
input_mtags <- readRDS("~/Documents/GitHub/AussieStrains/data/input_filt_rare_mtags331.rds")
input_methano <- filter_taxa_from_input(input = input_mtags,
                                        taxa_to_keep = "Methano")
# Make sure no non-methanogens with "Methano" in name like Methanoperedens
input_methano <- filter_taxa_from_input(input = input_methano,
                                        taxa_to_remove = "Methanoperedenaceae")
nrow(input_methano$data_loaded) # 18 OTUs
View(input_methano$taxonomy_loaded)
tax_sum_otu <- summarize_taxonomy(input = input_mtags,
                                  level = 8,
                                  relative = TRUE,
                                  report_higher_tax = FALSE) %>%
  filter(rownames(.) %in% input_methano$taxonomy_loaded$taxonomy8)
methanos <- data.frame(SampleID = colnames(tax_sum_otu),
                       Methano = colSums(tax_sum_otu)) %>%
  mutate(sampleID = gsub("X", "", SampleID)) %>%
  mutate(Methano_Cat = ifelse(Methano > 0, "Methano", "No Methano"))
range(methanos$Methano) # 0 to 0.03
mean(methanos$Methano) # 0.00027
sum(methanos$Methano > 0) # Present in 29/331 samples
results <- results %>%
  left_join(., methanos, by = "sampleID")
ggplot(results, aes(Per_aerobe, Methano)) +
  geom_point(size = 3, pch = 16, alpha = 0.75) +
  labs(x = "Predicted % aerobes",
       y = "% methanogens") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
                axis.text = element_text(size = 12))
ggplot(results, aes(Methano_Cat, Per_aerobe)) +
  geom_violin() +
  geom_jitter(size = 3, pch = 16, alpha = 0.75, height = 0, width = 0.1) +
  labs(y = "Predicted % aerobes",
       x = "Methanogen presence") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))



#### _NEON ####
# Use 2023 metagenomes and metadata
# Downloaded 247 metagenomes from JGI Data Portal (should be 248, but one file broken)



#### __Metadata ####
# Need to acquire and merge metadata for those samples
# Get all NEON sites metadata
# res <- httr::GET("https://data.neonscience.org/api/v0/sites")
# stop_for_status(res)
# sites_meta <- content(res, as = "text", encoding = "UTF-8") %>%
#   jsonlite::fromJSON(simplifyVector = TRUE)
# site_df <- sites_meta$data %>% 
#   as_tibble() %>%
#   dplyr::select(1:10)

# Or use the file provided here https://www.neonscience.org/field-sites/explore-field-sites
# This was updated Oct 8, 2025, use this for MAT and MAP
explore_sites <- read.csv("data/NEON_Field_Site_Metadata_20251008.csv") %>%
  filter(site_id %in% site_list) %>%
  mutate(Elevation = mean_evelation_m,
         MAT = as.numeric(gsub("°C", "", mean_annual_temperature_C)),
         MAP = mean_annual_precipitation_mm) %>%
  dplyr::select(site_id, site_name, latitude, longitude, site_state,
                Elevation, MAT, MAP,
                mean_canopy_height_m, avg_number_of_green_days)
# Also export coords for QGIS - need to confer Aridity index
# Remember: Global-AI values need to be multiplied by 0.0001
# Loaded 2.5 min BIOCLIM 12 and eto into QGIS
# Point sampling tool
# Aridity = MAP/PET
# coords <- explore_sites %>%
#   dplyr::select(site_id, latitude, longitude)
# write.csv(coords, "data/neon_coords.csv", row.names = FALSE)
neon_ai <- read.csv("data/neon_ai.csv") %>%
  dplyr::select(-latitude, -longitude)
  
# Run a bit of oxymetag code first to get the site_list, then get metadata for those
# Soil physical and chemical properties, periodic
# d15N, orgd13C, %N, %C, C:N
soilChem <- loadByProduct(dpID = 'DP1.10086.001',
                          startdate = "2023-01",
                          enddate = "2023-12",
                          check.size = FALSE,
                          package = 'expanded',
                          site = site_list)
View(soilChem$sls_soilChemistry)
View(soilChem$sls_metagenomicsPooling)
View(soilChem$sls_soilCoreCollection)
View(soilChem$sls_soilMoisture)
View(soilChem$sls_soilpH)
soilCN <- soilChem$sls_soilChemistry %>%
  dplyr::select(uid, siteID, plotID, sampleID, d15N, organicd13C, nitrogenPercent,
                organicCPercent, CNratio) %>%
  mutate(Date = substr(sampleID, nchar(sampleID) - 7, nchar(sampleID))) %>%
  mutate(SitePlotH = substr(sampleID, 1, 10)) %>%
  mutate(SitePlotHDate = paste(SitePlotH, Date, sep = "-")) %>%
  filter(SitePlotHDate %in% neon_per_aerobe$SitePlotHDate) %>%
  group_by(SitePlotHDate) %>%
  summarize(d15N = mean(d15N),
            organicd13C = mean(organicd13C),
            C = mean(organicCPercent),
            N = mean(nitrogenPercent),
            CN = mean(CNratio)) %>%
  ungroup() # n = 78
length(unique(soilCN$siteID)) # Only 9 of 26 sites
soilMoisture <- soilChem$sls_soilMoisture %>%
  dplyr::select(uid, siteID, plotID, sampleID, horizon, soilMoisture) %>%
  mutate(Date = substr(sampleID, nchar(sampleID) - 7, nchar(sampleID))) %>%
  mutate(SitePlotH = paste(plotID, horizon, sep = "-")) %>%
  mutate(SitePlotHDate = paste(SitePlotH, Date, sep = "-")) %>%
  filter(SitePlotHDate %in% neon_per_aerobe$SitePlotHDate) %>%
  group_by(SitePlotHDate) %>%
  summarize(SoilMoisture = mean(soilMoisture)) %>%
  ungroup() %>%
  mutate(SoilMoisture = SoilMoisture * 100)
soilChemMG <- soilChem$sls_metagenomicsPooling

# Soil physical and chemical properties, Megapit
# These were done once per site between 2012 and 2018
# These will have to just get one value for the site/horizon and then duplicate values
soilMega <- loadByProduct(dpID = 'DP1.00096.001',
                          startdate = "2012-01",
                          enddate = "2018-12",
                          check.size = FALSE,
                          package = 'expanded', 
                          site = site_list)
View(soilMega$mgp_perarchivesample)
View(soilMega$mgp_perbiogeosample) # Texture is here!
View(soilMega$mgp_perbulksample)
View(soilMega$mgp_perhorizon)
View(soilMega$mgp_permegapit) # Soil types
table(soilMega$mgp_permegapit$soilOrder)
soilType <- soilMega$mgp_permegapit %>%
  dplyr::select(siteID, nlcdClass, soilOrder)
soilTexture <- soilMega$mgp_perbiogeosample %>%
  filter(horizonName %in% c("Oa", "Oe", "Oi", 
                            "A", "A1", "A2", "A3", "Ac", "AE", "AO", "Ap")) %>%
  filter(!is.na(sandTotal)) %>%
  dplyr::select(siteID, horizonName, sandTotal, siltTotal, clayTotal) %>%
  group_by(siteID) %>%
  summarize(MegaSand = mean(sandTotal),
            MegaSilt = mean(siltTotal),
            MegaClay = mean(clayTotal)) %>%
  ungroup() # Have for 23/26 sites 

# Soil physical and chemical properties, distributed initial characterization
# These were done once per site between 2015 and 2021
# These will have to just get mean by site and then duplicate values
soilInit <- loadByProduct(dpID = 'DP1.10047.001',
                          startdate = "2015-01",
                          enddate = "2021-12",
                          check.size = FALSE,
                          package = 'expanded',
                          site = site_list)
View(soilInit$spc_biogeochem)
View(soilInit$spc_bulkdensity)
View(soilInit$spc_externalLabSummary) # Methods
View(soilInit$spc_particlesize) # Texture is here!
View(soilInit$spc_perhorizon)
View(soilInit$spc_perplot) # Soil types
soilTexture2 <- soilInit$spc_particlesize %>%
  dplyr::select(uid, siteID, plotID, horizonName, biogeoTopDepth, biogeoBottomDepth,
                sandTotal, siltTotal, clayTotal) %>%
  filter(plotID %in% neon_per_aerobe$SitePlot) %>%
  filter(biogeoBottomDepth <= 36) %>%
  filter(!is.na(sandTotal)) %>%
  mutate(Horizon = ifelse(horizonName == "Oa", "O", "M")) %>%
  mutate(SitePlotH = paste(plotID, Horizon, sep = "-")) %>%
  group_by(SitePlotH) %>%
  summarize(Sand = mean(sandTotal),
            Silt = mean(siltTotal),
            Clay = mean(clayTotal)) %>%
  ungroup() # Have for 57 samples matched by plot and horizon
  
metaDB <- neon.metaDB %>%
  filter(`ITS Proposal ID` == 509938) %>%
  filter(`Ecosystem Category` == "Terrestrial") %>%
  filter(grepl(2023, dnaSampleID))
table(metaDB$`Ecosystem Subtype`)

# Note some samples are not composited, need to join by Site_Plot_Horizon combo
metaDB_Hugh <- read.delim("data/neon_soil_metagenome_samples_2023.tsv") %>%
  drop_na(Genome.Size....assembled) %>% # Remove those not sequenced yet (n = 9)
  filter(grepl(2023, dnaSampleID)) %>%
  separate(depth..meters, remove = F, sep = " - ", into = c("topDepth", "botDepth")) %>%
  rename("Elevation" = `elevation..meters`) %>%
  rename("Ecosystem" = "ecosystem_subtype") %>%
  mutate(topDepth = as.numeric(topDepth),
         botDepth = as.numeric(botDepth)) %>%
  mutate(MAT = as.numeric(gsub(" Cel", "", `mean.annual.temperature`))) %>%
  mutate(MAP = as.numeric(gsub(" mm", "", `mean.annual.precipitation`))) %>%
  dplyr::select(dnaSampleID, Ecosystem, 
                #Elevation, MAP, MAT, # Get these from explore_sites
                pH, topDepth, botDepth,
                soil.horizon) %>%
  mutate(SitePlotH = substr(dnaSampleID, 1, 10))

neon_files <- read.csv("data/Globus_Download_521076_File_Manifest.csv") %>%
  mutate(MetaG.Name = gsub("-DNA1", "", Genome.Metagenome.Name))

neon_reads <- read.delim("data/neon_read_stats.tsv", header = FALSE) %>%
  separate(V2, into = c("Reads", "Mean", "Min", "Max"), sep = " ") %>%
  mutate(SampleID = gsub("BactReads/", "", V1)) %>%
  mutate(SampleID = gsub("_R1_bacterial.fastq.gz", "", SampleID)) %>%
  dplyr::select(-V1) %>%
  mutate(Reads = as.integer(Reads))
range(neon_reads$Reads)
mean(neon_reads$Reads)
se(neon_reads$Reads)
sd(neon_reads$Reads)



#### __oxymetag ####
# Merge everything!!
# Note there were 38 samples that were the individual samples
# Those also contained composite samples, so just use composite samples
neon_per_aerobe <- read.delim("data/neon_per_aerobe_predictions.tsv") %>%
  mutate(Per_anaerobe = 100 - Per_aerobe) %>%
  left_join(., neon_reads, by = "SampleID") %>% # Reads
  left_join(., neon_files, by = c("SampleID" = "Short.Organism.Name")) %>% # IDs
  mutate(Composite = ifelse(grepl("GEN", Genome.Metagenome.Name), "No", "Yes")) %>%
  filter(Composite == "Yes") %>%
  mutate(Site = substr(Genome.Metagenome.Name, 1, 4)) %>%
  mutate(SitePlot = substr(Genome.Metagenome.Name, 1, 8)) %>%
  mutate(SitePlotH = substr(Genome.Metagenome.Name, 1, 10)) %>%
  mutate(SitePlotHDate = substr(Genome.Metagenome.Name, 1, 19)) %>%
  left_join(., metaDB_Hugh, by = "SitePlotH") %>% # Ecosystem, pH, depth
  left_join(., soilCN, by = "SitePlotHDate") %>% # Soil C, N, n = 78
  left_join(., soilMoisture, by = "SitePlotHDate") %>% # Soil moisture
  left_join(., explore_sites, by = c("Site" = "site_id")) %>% # MAP, MAT, Elevation
  left_join(., neon_ai, by = c("Site" = "site_id")) %>% # Aridity index
  left_join(., soilType, by = c("Site" = "siteID")) %>% # NLCD class and soil order
  left_join(., soilTexture, by = c("Site" = "siteID")) %>% # Megapit
  left_join(., soilTexture2, by = "SitePlotH") %>% # Initial characterization
  mutate(Ecosystem2 = case_match(nlcdClass,
                                 "deciduousForest" ~ "Deciduous forest",
                                 "dwarfScrub" ~ "Dwarf scrub",
                                 "evergreenForest" ~ "Evergreen forest",
                                 "grasslandHerbaceous" ~ "Grassland",
                                 "shrubScrub" ~ "Shrub scrub",
                                 "woodyWetlands" ~ "Wetland")) %>%
  mutate(Dataset = "NEON (n = 247)")
names(neon_per_aerobe)
table(neon_per_aerobe$Composite)
TableS3 <- neon_per_aerobe %>%
  mutate(ClimateClass = ifelse(AI < 0.03, "Hyper arid",
                               ifelse(AI >= 0.03 & AI < 0.2, "Arid",
                                      ifelse(AI >= 0.2 & AI < 0.5, "Semi-arid",
                                             ifelse(AI >= 0.5 & AI < 0.65, "Dry sub-humid",
                                                    "Humid"))))) %>%
  dplyr::select(-Ecosystem) %>%
  rename(Ecosystem = Ecosystem2) %>%
  mutate(Ecosystem2 = case_match(Ecosystem,
                                 "Dwarf scrub" ~ "Shrubland",
                                 "Shrub scrub" ~ "Shrubland",
                                 "Evergreen forest" ~ "Forest",
                                 "Deciduous forest" ~ "Forest",
                                 .default = Ecosystem)) %>%
  dplyr::select(-Dataset, -Composite, -md5.checksum)
# write_xlsx(TableS3, "~/Desktop/Fierer/AEGIS/Oxygen/Manuscript/TableS3.xlsx",
#            format_headers = FALSE)
neon_per_aerobe <- read_xlsx("~/Desktop/Fierer/AEGIS/Oxygen/Manuscript/TableS3.xlsx",
                             skip = 1)

# Basic info
site_list <- unique(neon_per_aerobe$Site)
table(neon_per_aerobe$Site)
table(neon_per_aerobe$`Ecosystem Subtype`)

hist(neon_per_aerobe$Per_aerobe)
range(neon_per_aerobe$Per_aerobe)
mean(neon_per_aerobe$Per_aerobe)
se(neon_per_aerobe$Per_aerobe)
sd(neon_per_aerobe$Per_aerobe)
median(neon_per_aerobe$Per_aerobe)
sum(neon_per_aerobe$Per_aerobe == 100) # 164 100, 45 < 100

range(neon_per_aerobe$Reads)
mean(neon_per_aerobe$Reads)
se(neon_per_aerobe$Reads)
sd(neon_per_aerobe$Reads)

# Distribution
pdf("InitialFigs/NEON_PerAerobe_histogram.pdf", width = 7, height = 5)
ggplot(neon_per_aerobe, aes(Per_aerobe)) +
  geom_histogram() +
  labs(x = "Predicted aerobes (%)", y = "Count") +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
dev.off()

ggplot(neon_per_aerobe, aes(Dataset, Per_aerobe)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 2, pch = 16, alpha = 0.5) +
  geom_ysidedensity(aes(x = after_stat(density))) +
  labs(x = NULL, y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme_ggside_minimal()

# Predictor variables
pdf("InitialFigs/NEON_PerAerobe_Sites.pdf", width = 7, height = 5)
ggplot(neon_per_aerobe, aes(reorder(Site, Per_aerobe, mean), Per_aerobe)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 2, pch = 16, alpha = 0.75) +
  labs(x = "Site", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))
dev.off()

ggplot(neon_per_aerobe, aes(reorder(Ecosystem, Per_aerobe, mean), Per_aerobe)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 2, pch = 16, alpha = 0.75, width = 0.25) +
  labs(x = "Ecosystem", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(reorder(Ecosystem2, Per_aerobe, mean), Per_aerobe)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 2, pch = 16, alpha = 0.75, width = 0.25) +
  labs(x = "Ecosystem", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(reorder(soilOrder, Per_aerobe, mean), Per_aerobe)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 2, pch = 16, alpha = 0.75) +
  labs(x = "Soil order", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(MegaSand, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "% Sand", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(Sand, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "% Sand", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(SoilMoisture, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "GWC (%)", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(data = subset(neon_per_aerobe, Ecosystem2 != "Wetland"),
       aes(SoilMoisture, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "GWC (%)", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(MAP, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "MAP (mm)", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(MAT, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "MAT (C)", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(AI, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "AI", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(pH, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "pH", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(Elevation, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "Elevation (m)", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(mean_canopy_height_m, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "Canopy height", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

ggplot(neon_per_aerobe, aes(avg_number_of_green_days, Per_aerobe)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth() +
  labs(x = "Green days", y = "Predicted aerobes (%)") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14))

an <- neon_per_aerobe %>%
  filter(Per_aerobe < 100) %>%
  droplevels()
table(an$Site) # 13 of 26 sites
an_sites <- neon_per_aerobe %>%
  group_by(Site) %>%
  summarize(mean = mean(Per_aerobe)) %>%
  ungroup() %>%
  filter(mean < 100) %>%
  left_join(., site_df, by = c("Site" = "siteCode")) %>%
  dplyr::select(Site, siteName, everything()) %>%
  arrange(mean)



#### _Categorical, Fig 6/S4 ####
# For Figure 6, we need to plot by habitat and add histogram as an inset
# Make climate classes
# < 0.03 Hyper Arid
# 0.03 – 0.2 Arid
# 0.2 – 0.5 Semi-Arid
# 0.5 – 0.65 Dry sub-humid
# > 0.65 Humid
# Panel for Australia and panel for U.S.
base6 <- aus %>%
  dplyr::select(vegetation_type, Per_aerobe, AI, water_content) %>%
  mutate(Per_anaerobe = (100 - Per_aerobe) / 100) %>%
  mutate(ClimateClass = ifelse(AI < 0.03, "Hyper arid",
                               ifelse(AI >= 0.03 & AI < 0.2, "Arid",
                                      ifelse(AI >= 0.2 & AI < 0.5, "Semi-arid",
                                             ifelse(AI >= 0.5 & AI < 0.65, "Dry sub-humid",
                                                    "Humid"))))) %>%
  rename(Ecosystem = vegetation_type,
         SoilMoisture = water_content) %>%
  mutate(Ecosystem2 = case_match(Ecosystem,
                                 "Savannah" ~ "Grassland",
                                 "Heathland" ~ "Shrubland",
                                 .default = Ecosystem)) %>%
  mutate(Dataset = "a) Australia (n = 331)",
         Dataset2 = "a) Australia (n = 247)")
table(base6$ClimateClass)
table(base6$Ecosystem2)
neon6 <- neon_per_aerobe %>%
  dplyr::select(Ecosystem2, Per_aerobe, AI, SoilMoisture) %>%
  mutate(Per_anaerobe = (100 - Per_aerobe) / 100) %>%
  mutate(ClimateClass = ifelse(AI < 0.03, "Hyper arid",
                               ifelse(AI >= 0.03 & AI < 0.2, "Arid",
                                      ifelse(AI >= 0.2 & AI < 0.5, "Semi-arid",
                                             ifelse(AI >= 0.5 & AI < 0.65, "Dry sub-humid",
                                                    "Humid"))))) %>%
  rename(Ecosystem = Ecosystem2) %>%
  mutate(Ecosystem2 = case_match(Ecosystem,
                                 "Dwarf scrub" ~ "Shrubland",
                                 "Shrub scrub" ~ "Shrubland",
                                 "Evergreen forest" ~ "Forest",
                                 "Deciduous forest" ~ "Forest",
                                 .default = Ecosystem)) %>%
  mutate(Dataset = "b) U.S. (n = 209)",
         Dataset2 = "b) U.S. (n = 209)")
table(neon6$ClimateClass)
table(neon6$Ecosystem2)


# Stats
base6 <- base6 %>%
  dplyr::select(Per_anaerobe, Ecosystem, Ecosystem2, ClimateClass, Dataset)
b0 <- gamlss(base6$Per_anaerobe ~ 1, family = BEZI, trace = F)
b1 <- gamlss(base6$Per_anaerobe ~ base6$Ecosystem, family = BEZI, trace = F)
LR.test(b0, b1) # p = 0.20
summary(b1)
ll_null <- logLik(gamlss(Per_anaerobe ~ 1, data = base6, family = BEZI))
ll_full <- logLik(b1)
n <- nrow(base6)
r2_nagelkerke <- as.numeric((1 - exp((2/n)*(ll_null - ll_full))) / (1 - exp(2*ll_null/n)))
r2_nagelkerke
summary_b1 <- summary(b1)
summary_b1

neon6 <- neon6 %>%
  dplyr::select(Per_anaerobe, Ecosystem, Ecosystem2, ClimateClass, Dataset)
b0 <- gamlss(neon6$Per_anaerobe ~ 1, family = BEZI, trace = F)
b1 <- gamlss(neon6$Per_anaerobe ~ neon6$Ecosystem, family = BEZI, trace = F)
LR.test(b0, b1) # p = 0.50
summary(b1)
ll_null <- logLik(gamlss(Per_anaerobe ~ 1, data = neon6, family = BEZI))
ll_full <- logLik(b1)
n <- nrow(neon6)
r2_nagelkerke <- as.numeric((1 - exp((2/n)*(ll_null - ll_full))) / (1 - exp(2*ll_null/n)))
r2_nagelkerke
summary_b1 <- summary(b1)
summary_b1



b0 <- gamlss(base6$Per_anaerobe ~ 1, family = BEZI, trace = F)
b1 <- gamlss(base6$Per_anaerobe ~ base6$Ecosystem2, family = BEZI, trace = F)
LR.test(b0, b1) # p = 0.25
summary(b1)
ll_null <- logLik(gamlss(Per_anaerobe ~ 1, data = base6, family = BEZI))
ll_full <- logLik(b1)
n <- nrow(base6)
r2_nagelkerke <- as.numeric((1 - exp((2/n)*(ll_null - ll_full))) / (1 - exp(2*ll_null/n)))
r2_nagelkerke
summary_b1 <- summary(b1)
summary_b1

b0 <- gamlss(neon6$Per_anaerobe ~ 1, family = BEZI, trace = F)
b1 <- gamlss(neon6$Per_anaerobe ~ neon6$Ecosystem2, family = BEZI, trace = F)
LR.test(b0, b1) # p = 0.23
summary(b1)
ll_null <- logLik(gamlss(Per_anaerobe ~ 1, data = neon6, family = BEZI))
ll_full <- logLik(b1)
n <- nrow(neon6)
r2_nagelkerke <- as.numeric((1 - exp((2/n)*(ll_null - ll_full))) / (1 - exp(2*ll_null/n)))
r2_nagelkerke
summary_b1 <- summary(b1)
summary_b1


b0 <- gamlss(base6$Per_anaerobe ~ 1, family = BEZI, trace = F)
b1 <- gamlss(base6$Per_anaerobe ~ base6$ClimateClass, family = BEZI, trace = F)
LR.test(b0, b1) # p = 0.80
summary(b1)
ll_null <- logLik(gamlss(Per_anaerobe ~ 1, data = base6, family = BEZI))
ll_full <- logLik(b1)
n <- nrow(base6)
r2_nagelkerke <- as.numeric((1 - exp((2/n)*(ll_null - ll_full))) / (1 - exp(2*ll_null/n)))
r2_nagelkerke # 0.008
summary_b1 <- summary(b1)
summary_b1

b0 <- gamlss(neon6$Per_anaerobe ~ 1, family = BEZI, trace = F)
b1 <- gamlss(neon6$Per_anaerobe ~ neon6$ClimateClass, family = BEZI, trace = F)
LR.test(b0, b1) # p = 0.24
summary(b1)
ll_null <- logLik(gamlss(Per_anaerobe ~ 1, data = neon6, family = BEZI))
ll_full <- logLik(b1)
n <- nrow(neon6)
r2_nagelkerke <- as.numeric((1 - exp((2/n)*(ll_null - ll_full))) / (1 - exp(2*ll_null/n)))
r2_nagelkerke # 0.042
summary_b1 <- summary(b1)
summary_b1


# Plot
fig6_df <- rbind(base6, neon6) %>%
  mutate(ClimateClass = factor(ClimateClass,
                               levels = c("Arid", "Semi-arid", 
                                          "Dry sub-humid", "Humid")))

figS4 <- ggplot(fig6_df, aes(ClimateClass, Per_aerobe)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 2, pch = 21, alpha = 0.75, width = 0.25, height = 0, 
              colour = "black", aes(fill = ClimateClass)) +
  scale_fill_manual(values = c("red", "orange", "yellow", "blue")) +
  labs(x = "Climate class", 
       y = "Predicted aerobes (%)") +
  facet_wrap(~ Dataset, scales = "free_x") +
  scale_x_reordered() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "none")
pdf("FinalFigs/FigureS4.pdf", width = 7, height = 5)
figS4
dev.off()

# Draft1 with Ecosystem
main <- ggplot(fig6_df, aes(reorder_within(Ecosystem, Per_aerobe, Dataset), Per_aerobe)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 2, pch = 21, alpha = 0.75, width = 0.25, height = 0, 
              colour = "white", fill = "black") +
  labs(x = "Ecosystem", 
       y = "Predicted aerobes (%)") +
  facet_wrap(~ Dataset, scales = "free_x") +
  scale_x_reordered() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14))
main

# Histogram insets
i1 <- ggplot(aus, aes(Per_aerobe)) +
  geom_histogram() +
  labs(x = "Predicted aerobes (%)", y = "Count") +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        #axis.title = element_text(size = 10),
        axis.title = element_blank())

i2 <- ggplot(neon_per_aerobe, aes(Per_aerobe)) +
  geom_histogram() +
  labs(x = "Predicted aerobes (%)", y = "Count") +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        #axis.title = element_text(size = 10),
        axis.title = element_blank())

fig6 <- ggdraw() +
  draw_plot(main) +
  draw_plot(i1, x = 0.33, y = 0.23, width = 0.2, height = 0.2) +
  draw_plot(i2, x = 0.78, y = 0.23, width = 0.2, height = 0.2)
pdf("FinalFigs/Figure6.pdf", width = 7, height = 5)
fig6
dev.off()

# Draft 2 with Ecosystem 2
main <- ggplot(fig6_df, aes(reorder_within(Ecosystem2, Per_aerobe, Dataset), Per_aerobe)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 2, pch = 21, alpha = 1, width = 0.25, height = 0, 
              colour = "black", aes(fill = Ecosystem2)) +
  scale_fill_manual(values = c("antiquewhite", "darkgreen", "gold", "chartreuse3",
                               "blue", "brown")) +
  labs(x = "Ecosystem", 
       y = "Predicted aerobes (%)") +
  facet_wrap(~ Dataset, scales = "free_x") +
  scale_x_reordered() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "none")
main

# Histogram insets
i1 <- ggplot(aus, aes(Per_aerobe)) +
  geom_histogram() +
  labs(x = "Predicted aerobes (%)", y = "Count") +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        #axis.title = element_text(size = 10),
        axis.title = element_blank())

i2 <- ggplot(neon_per_aerobe, aes(Per_aerobe)) +
  geom_histogram() +
  labs(x = "Predicted aerobes (%)", y = "Count") +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        #axis.title = element_text(size = 10),
        axis.title = element_blank())

fig6 <- ggdraw() +
  draw_plot(main) +
  draw_plot(i1, x = 0.33, y = 0.23, width = 0.2, height = 0.2) +
  draw_plot(i2, x = 0.78, y = 0.23, width = 0.2, height = 0.2)
pdf("FinalFigs/Figure6.pdf", width = 7, height = 5)
fig6
dev.off()
png("FinalFigs/Figure6.png", width = 7, height = 5, units = "in", res = 300)
fig6
dev.off()



#### _Continuous, Fig S5 ####
# Australia and U.S., run zero-inflated beta regression for continuous variables
# Both datasets have MAP, MAT, AI for all samples
# Then, each dataset has a suite of soil variable, with varying amounts of NAs
# Report these results as a table

# Australia
aus_uni <- aus_uni %>%
  mutate(ClaySilt = 100 - ClaySilt) %>%
  rename(Sand = ClaySilt) %>%
  mutate(Per_anaerobe = (100 - aus_uni$Per_aerobe) / 100)
m_aus <- as.data.frame(matrix(NA, 11, 4)) %>%
  set_names("Variable", "R2", "p", "n")
for (i in 1:11) {
  # Record the variable tested
  m_aus$Variable[i] <- names(aus_uni)[i]
  
  # Get a dataframe with no NA for that specific variable
  d <- aus_uni %>%
    filter(!is.na(.[[i]])) %>%
    dplyr::select(i, Per_anaerobe)
  
  # Run zero-inflated beta regression, get metrics
  b0 <- gamlss(Per_anaerobe ~ 1, family = BEZI, trace = F, data = d)
  b1 <- gamlss(Per_anaerobe ~ d[[1]], family = BEZI, trace = F, data = d)
  lr <- LR.test(b0, b1, print = FALSE)
  ll_null <- logLik(b0)
  ll_full <- logLik(b1)
  n <- nrow(d)
  r2_nagelkerke <- as.numeric((1-exp((2/n)*(ll_null-ll_full)))/(1-exp(2*ll_null/n)))
  
  # Add Nagelkerke pseudo-R2 and p
  m_aus$R2[i] <- round(r2_nagelkerke, 3)
  m_aus$p[i] <- round(lr$p.val, 3)
  m_aus$n[i] <- n
}


# U.S.
us_uni <- neon_per_aerobe %>%
  dplyr::select(Sand, pH, C, N, CN, SoilMoisture, Elevation, MAT, MAP, AI,
                Per_anaerobe) %>%
  mutate(Per_anaerobe = Per_anaerobe/100)
m_us <- as.data.frame(matrix(NA, 10, 4)) %>%
  set_names("Variable", "R2", "p", "n")
for (i in 1:10) {
  # Record the variable tested
  m_us$Variable[i] <- names(us_uni)[i]
  
  # Get a dataframe with no NA for that specific variable
  d <- us_uni %>%
    filter(!is.na(.[[i]])) %>%
    dplyr::select(i, Per_anaerobe)
  
  # Run zero-inflated beta regression, get metrics
  b0 <- gamlss(Per_anaerobe ~ 1, family = BEZI, trace = F, data = d)
  b1 <- gamlss(Per_anaerobe ~ d[[1]], family = BEZI, trace = F, data = d)
  lr <- LR.test(b0, b1, print = FALSE)
  ll_null <- logLik(b0)
  ll_full <- logLik(b1)
  n <- nrow(d)
  r2_nagelkerke <- as.numeric((1-exp((2/n)*(ll_null-ll_full)))/(1-exp(2*ll_null/n)))
  
  # Add Nagelkerke pseudo-R2 and p
  m_us$R2[i] <- round(r2_nagelkerke, 3)
  m_us$p[i] <- round(lr$p.val, 3)
  m_us$n[i] <- n
}

m_aus
m_us



# Plot
d_aus <- base6 %>%
  filter(!is.na(SoilMoisture))
b1 <- gamlss(Per_anaerobe ~ SoilMoisture, family = BEZI, trace = F, data = d_aus)
newdat <- data.frame(SoilMoisture = seq(min(d_aus$SoilMoisture), max(d_aus$SoilMoisture),
                                        length.out = 100))
pred <- predict(b1, newdata = newdat, type = "response")
newdat$fit <- pred
newdat$Dataset2 <- "a) Australia (n = 247)"
ggplot(d_aus, aes(x = SoilMoisture, y = Per_anaerobe, colour = Ecosystem2)) +
  geom_point(size = 2, pch = 16, alpha = 0.5) +
  geom_line(data = newdat, aes(y = fit), color = "blue", linewidth = 1.2) +
  labs(x = "Gravimetric water content (g water/g dry soil)", 
       y = "Predicted proportion anaerobes") +
  theme_bw()

d_us <- neon_per_aerobe %>%
  filter(!is.na(Sand)) %>%
  mutate(Per_anaerobe = Per_anaerobe/100) %>%
  dplyr::select(Sand, Per_anaerobe, Ecosystem2)
b2 <- gamlss(Per_anaerobe ~ Sand, family = BEZI, trace = F, data = d_us)
newdat2 <- data.frame(Sand = seq(min(d_us$Sand), max(d_us$Sand),
                                 length.out = 100))
pred2 <- predict(b2, newdata = newdat2, type = "response")
newdat2$fit <- pred2
ggplot(d_us, aes(x = Sand, y = Per_anaerobe, colour = Ecosystem2)) +
  geom_point(size = 2, pch = 16, alpha = 0.5) +
  geom_line(data = newdat2, aes(y = fit), color = "blue", linewidth = 1.2) +
  labs(x = "% Sand", 
       y = "Predicted proportion anaerobes") +
  theme_bw()

d_us <- neon6 %>%
  filter(!is.na(SoilMoisture))
b2 <- gamlss(Per_anaerobe ~ SoilMoisture, family = BEZI, trace = F, data = d_us)
newdat2 <- data.frame(SoilMoisture = seq(min(d_us$SoilMoisture), max(d_us$SoilMoisture),
                                         length.out = 100))
pred2 <- predict(b2, newdata = newdat2, type = "response")
newdat2$fit <- pred2
newdat2$Dataset2 <- "b) U.S. (n = 209)"
ggplot(d_us, aes(x = SoilMoisture, y = Per_anaerobe, colour = Ecosystem2)) +
  geom_point(size = 2, pch = 16, alpha = 0.5) +
  geom_line(data = newdat2, aes(y = fit), color = "blue", linewidth = 1.2) +
  labs(x = "Gravimetric water content (% dry mass)", 
       y = "Predicted proportion anaerobes") +
  theme_bw()



# Combine
fig7 <- rbind(d_aus, d_us)
nd <- rbind(newdat, newdat2)
lab <- data.frame("SoilMoisture" = c(20, 250),
                  "Per_anaerobe" = c(0.5, 0.5),
                  "label" = c("atop(R^2 == 0.059, italic(p) == 0.026)",
                              "atop(R^2 == 0.08, italic(p) == 0.005)"),
                  "Dataset2" = c("a) Australia (n = 247)",
                                 "b) U.S. (n = 209)"))
fig7 <- ggplot(fig7, aes(x = SoilMoisture, y = Per_anaerobe)) +
  geom_point(size = 3, pch = 21, alpha = 1, aes(fill = Ecosystem2)) +
  geom_line(data = nd, aes(y = fit), color = "blue", linewidth = 1.2) +
  geom_text(data = lab, aes(label = label), parse = TRUE) +
  scale_fill_manual(values = c("antiquewhite", "darkgreen", "gold", "chartreuse3",
                               "blue", "brown")) +
  labs(x = "Gravimetric water content (% dry mass)", 
       y = "Predicted proportion anaerobes",
       fill = "Ecosystem") +
  facet_wrap(~ Dataset2, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "right")
pdf("FinalFigs/FigureS5.pdf", width = 7, height = 5)
figS5
dev.off()
png("FinalFigs/FigureS5.png", width = 7, height = 5, units = "in", res = 300)
figS5
dev.off()



#### 6. Map ####
# Make a map to include as a supplemental figure
# Australia map was already made for Bueno de Mesquita et al. 2025
# Make similar U.S. map
d_331 <- read.delim("~/Documents/GitHub/AussieStrains/data/metadata_331.txt") %>%
  mutate("ClimateClass" = ifelse(AI < 0.5, "Arid to semi-arid",
                                 ifelse(AI >= 0.5 & AI < 0.65, "Dry sub-humid",
                                        "Humid")))
coords_trans <- st_as_sf(d_331,
                         coords = c('longitude', 'latitude'),
                         crs=4326)
sf_oz <- ozmap("states")
map <- ggplot(data = sf_oz) +
  geom_sf(fill = "grey90", color = "white") +
  geom_sf(data = coords_trans,
          aes(fill = AI, shape = ClimateClass),
          size = 3, alpha = 1, color = "black", stroke = 0.3) +
  scale_shape_manual(values = c(21, 24, 22)) +
  annotation_scale() +
  annotation_north_arrow(pad_x = unit(1, "cm"), pad_y = unit(1, "cm"),
                         height = unit(1, "cm"), width = unit(1, "cm"),) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  guides(shape = guide_legend(order = 2),
         fill = guide_colorbar(order = 1)) +
  xlim(111, 155) +
  ylim(43, 12) +
  labs(fill = "Aridity\nindex") +
  theme_minimal()
map

neon_per_aerobe <- neon_per_aerobe %>%
  mutate(ClimateClass = ifelse(AI < 0.03, "Hyper arid",
                               ifelse(AI >= 0.03 & AI < 0.2, "Arid",
                                      ifelse(AI >= 0.2 & AI < 0.5, "Semi-arid",
                                             ifelse(AI >= 0.5 & AI < 0.65, "Dry sub-humid",
                                                    "Humid")))))
coords_trans2 <- st_as_sf(neon_per_aerobe,
                          coords = c('longitude', 'latitude'),
                          crs=4326)
states <- get_urbn_map("states", sf = TRUE)
ggplot(states) +
  geom_sf(fill = "gray90", color = "white") +
  geom_sf(data = coords_trans2,
          aes(fill = AI, shape = ClimateClass),
          size = 3, alpha = 1, color = "black", stroke = 0.3) +
  scale_shape_manual(values = c(21, 24, 22, 23)) +
  annotation_scale() +
  annotation_north_arrow(pad_x = unit(1, "cm"), pad_y = unit(1, "cm"),
                         height = unit(1, "cm"), width = unit(1, "cm"),) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  guides(shape = guide_legend(order = 2),
         fill = guide_colorbar(order = 1)) +
  #xlim(111, 155) +
  #ylim(43, 12) +
  labs(fill = "Aridity\nindex") +
  theme_minimal()


library(usmap)
library(ggplot2)

# Example NEON-like sites
sites <- data.frame(
  siteID = c("BONA", "DELA", "CPER"),
  lon = c(-147.5, -122.0, -104.7),
  lat = c(65.0, 46.0, 40.8)
)

# rename columns to x/y so usmap recognises them
sites_usmap <- rename(sites, x = lon, y = lat)

plot_usmap(regions = "states", fill = "gray90", color = "white") +
  geom_point(
    data = sites_usmap,
    aes(x = x, y = y),
    color = "red",
    size = 2
  ) +
  theme_void()



# End Script