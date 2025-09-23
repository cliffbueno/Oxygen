# Calculate % aerobic bacteria in modern metagenomes
# by Cliff Bueno de Mesquita, Fierer Lab, 2025
# for AEGIS
# Use the model from the simulations



#### 1. Setup ####
# Libraries
library(tidyverse)
library(mgcv)
library(readxl)
library(car)
library(bestglm)

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

# Tran et al. 2021 Oxygen
# (metaG don't have O2, but can infer from depth)
tran <- read.csv("data/TranDO.csv")
ggplot(tran, aes(Depth, DO)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = 122, color = "red") +
  geom_hline(yintercept = 12, color = "red") +
  geom_vline(xintercept = 1, color = "red") +
  geom_hline(yintercept = 84, color = "red") +
  labs(x = "Depth (m)",
       y = "DO") +
  theme_bw()
tran_model <- gam(DO ~ s(Depth), data = tran)
summary(tran_model)
plot(tran_model)
formula(tran_model)
new_tran_data <- data.frame(Depth = c(0, 120, 150, 200, 250, 300, 400, 1200))
predicted_tran_o2 <- predict(tran_model, newdata = new_tran_data)
predicted_tran_o2
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
  mutate(Gene.length = abs(send - sstart)) %>%
  left_join(., map, by = c("sseqid" = "Header"))
hist(d$Gene.length) # Gene length of the hits
d$sseqid[1365] # 2 pfams (cytochrome c)
table(d$Pfam) # PF05425 (aerobic) and PF13597 (anaerobic) are missing!

gene.hits = d %>% 
  group_by(Pfam) %>% 
  summarise(total_count=n())

# mean.gene.lengths <- d %>% group_by(Pfam) %>% 
#   summarise_at(vars(Gene.length), list(name = mean)) %>%
#   set_names(c("Pfam", "Gene.length"))
# gene.hit.length.correction <- merge(mean.gene.lengths, gene.hits, by = "Pfam") %>%
#   mutate(RPK = total_count / (1000*Gene.length)) %>%
#   left_join(., oxygen_pfams, by = "Pfam")
gene.hit.length.correction <- gene.hits %>%
  left_join(., pfam_gene_length, by = "Pfam") %>%
  mutate(RPK = total_count / (1000*Gene.length)) %>%
  left_join(., oxygen_pfams, by = "Pfam")

oxygen_rpk <- gene.hit.length.correction %>%
  group_by(Oxygen) %>%
  summarize(RPKsum = sum(RPK))

ratio <- oxygen_rpk$RPKsum[1] / oxygen_rpk$RPKsum[2]

new_data <- data.frame(ratio = ratio)  # or a sequence of values
predicted_y <- predict(gam_model, newdata = new_data)
predicted_y # 82.32% aerobe (old RPK was 79.56% aerobe)

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
  mutate(Gene.length = abs(send - sstart)) %>%
  left_join(., map, by = c("sseqid" = "Header"))
hist(t$Gene.length)
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
#saveRDS(results, "data/results_modern.rds")



#### 3. Analyze ####
# Now analyze the predicted % aerobe by habitat and depth
meta <- read_xlsx("data/ModernMetaG.xlsx")
length(unique(meta$sampleID)) # n = 278
results <- readRDS("data/results_modern.rds")
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
table(d_plot$Habitat3)
table(d_plot$Type)
m <- aov(Per_aerobe ~ Habitat3, data = d_plot)
summary(m) # 15 281245   18750   235.7 <2e-16 ***
fig2 <- ggplot(d_plot, aes(reorder(Habitat3, Per_aerobe, mean), Per_aerobe)) +
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
fig2
pdf("FinalFigs/Figure2.pdf", width = 7, height = 5)
fig2
dev.off()
png("FinalFigs/Figure2.png", width = 7, height = 5, units = "in", res = 300)
fig2
dev.off()

# Reads
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
  labs(x = NULL, y = "Predicted aerobe relative abundance (%)") +
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
summary(m) # p = 0.14
ggplot(bs, aes(Oxygen, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "Oxygen", y = "Aerobe relative abundance (%)") +
  ylim(0, 100) +
  scale_x_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

m <- lm(Per_aerobe ~ poly(Depth_cm, 3, raw = T), data = bs)
summary(m) # p = 0.016
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
summary(m) # p = 0.14
m <- lm(Per_aerobe ~ poly(O2_H2S, 3, raw = T), data = bs)
summary(m) # p< 0.001
pdf("InitialFigs/BlackSea_Aerobe_O2_H2S.pdf", width = 7, height = 5)
ggplot(bs, aes(O2_H2S, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth() +
  labs(x = "[O2] / [H2S]", y = "Predicted aerobe relative abundance (%)") +
  scale_x_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()



#### _Baltic ####
# Baltic Sea sediments, oxygen gradient
baltic <- d %>%
  filter(Study == "CaseStudy5") %>%
  droplevels()
table(baltic$Habitat)
m <- lm(Per_aerobe ~ Oxygen, data = baltic)
summary(m) # p = 0.009, R2 = 0.42
pdf("InitialFigs/Baltic_Aerobe_O2.pdf", width = 7, height = 5)
ggplot(baltic, aes(Oxygen, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "Oxygen (mg/L)", y = "Predicted aerobe relative abundance (%)") +
  ylim(0, 100) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()



#### _Freshwater Lake ####
# Have aerobic and anaerobic Lake Tanganyika samples
lake <- d %>%
  filter(Habitat == "Freshwater lake") %>%
  droplevels()
table(lake$Habitat)
table(lake$Habitat2)
m <- t.test(Per_aerobe ~ Habitat2, data = lake)
m # p = 0.0007
pdf("InitialFigs/Tanganyika_Aerobe_Depth.pdf", width = 7, height = 5)
ggplot(lake, aes(Depth_cat, Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 3, alpha = 0.75, width = 0.25) +
  labs(x = "Depth", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12))
dev.off()

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
  labs(x = "Oxygen (mg/L)", y = "Predicted aerobe relative abundance (%)") +
  ylim(0, 100) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))



#### _OWC ####
owc <- d %>%
  filter(Study == "CaseStudy2") %>%
  droplevels()
table(owc$Habitat)
m <- lm(Per_aerobe ~ Oxygen, data = owc)
summary(m) # p = 0.06
m <- lm(Per_aerobe ~ log(Oxygen+1), data = owc)
summary(m) # p = 0.08
pdf("InitialFigs/OWC_Aerobe_O2.pdf", width = 7, height = 5)
ggplot(owc, aes(Oxygen, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "dashed") +
  labs(x = "Oxygen (µM)", y = "Predicted aerobe relative abundance (%)") +
  #ylim(0, 100) +
  #scale_x_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()



#### _CZO ####
# 1 site on microbe, 4 sites on MG-RAST (merged, qc'd by them)
czo <- d %>%
  filter(Study == "CaseStudy3") %>%
  mutate(Depth_cm = as.numeric(Depth_cm)) %>%
  mutate(Site = factor(Site,
                       levels = c("Shale Hills",
                                  "Boulder",
                                  "Catalina-Jemez",
                                  "IML",
                                  "Luquillo"))) %>%
  droplevels()
table(czo$Habitat)
table(czo$Site)
m <- lm(Per_aerobe ~ Site * Depth_cm, data = czo)
summary(m)
Anova(m)
Anova(m, type = "III")
pdf("InitialFigs/CZO_Aerobe_Depth.pdf", width = 7, height = 5)
ggplot(czo, aes(Depth_cm, Per_aerobe, fill = Site, colour = Site)) +
  geom_point(size = 3, pch = 21, alpha = 1, position = position_dodge(width = 4),
             colour = "black") +
  geom_line(position = position_dodge(width = 4)) +
  #geom_smooth(method = "lm", linewidth = 0.5, se = F) +
  labs(x = "Soil depth (cm)", y = "Predicted aerobe relative abundance (%)") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()



#### _Rhizo ####
# Compare bulk soil vs. rhizosphere soil
rhizo <- d %>%
  filter(Study == "CaseStudy4") %>%
  droplevels()
table(rhizo$Treatment)
table(rhizo$SoilType)
m <- t.test(Per_aerobe ~ SoilType, data = rhizo)
m # p = 0.08
ggplot(rhizo, aes(SoilType, Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_point(size = 3, alpha = 0.75) +
  labs(x = NULL, y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12))







#### _Hot Spring ####
hot <- d %>%
  filter(Habitat == "Hot spring sediment") %>%
  droplevels()
table(hot$Habitat)
m <- lm(Per_aerobe ~ Oxygen, data = hot)
summary(m) # p = 0.84
m <- lm(Per_aerobe ~ Temp, data = hot)
summary(m) # p = 0.85
m <- lm(Per_aerobe ~ pH, data = hot)
summary(m) # p = 0.27
m <- aov(Per_aerobe ~ Site, data = hot)
summary(m) # p = 0.175
pdf("InitialFigs/HotSpring_Aerobe_O2.pdf", width = 7, height = 5)
ggplot(hot, aes(Oxygen, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "dashed") +
  labs(x = "Oxygen (mg/L)", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()

ggplot(hot, aes(Temp, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "dashed") +
  labs(x = "Temperature (˚C)", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggplot(hot, aes(pH, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "dashed") +
  labs(x = "pH", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggplot(hot, aes(Site, Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_point(size = 3, alpha = 0.75) +
  labs(x = "Region", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12))



#### _Ocean Water ####
# Have aerobic and anaerobic ocean water samples
# But the low O2 ones are still not very deep, and hypoxic, not anoxic
# Deep chlorophyll maximum - Malaspina Expedition
# Surface - Tara Oceans Expedition
ocean <- d %>%
  filter(Habitat == "Ocean water") %>%
  droplevels()
table(ocean$Habitat)
m <- lm(Per_aerobe ~ Oxygen, data = ocean)
summary(m) # p = 0.06
pdf("InitialFigs/Ocean_Aerobe_O2.pdf", width = 7, height = 5)
ggplot(ocean, aes(Oxygen, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "dashed") +
  labs(x = "Oxygen (µmol/kg)", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()

ggplot(ocean, aes(Habitat2, Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 3, alpha = 0.75, width = 0.25) +
  labs(x = NULL, y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))



#### _Hydrothermal ####
# Brothers Volcano, "Hydrothermal deposit"
vent <- d %>%
  filter(Habitat == "Hydrothermal deposit") %>%
  droplevels()
table(vent$Habitat)
m <- lm(Per_aerobe ~ Temp, data = vent)
summary(m) # p = 0.002
m <- lm(Per_aerobe ~ pH, data = vent)
summary(m) # p = 0.1
m <- aov(Per_aerobe ~ Site, data = vent)
summary(m) # p = 0.193

ggplot(vent, aes(Temp, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "solid") +
  labs(x = "Temperature (˚C)", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggplot(vent, aes(pH, Per_aerobe)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "dashed") +
  labs(x = "pH", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggplot(vent, aes(Site, Per_aerobe)) +
  geom_boxplot(outliers = F) +
  geom_point(size = 3, alpha = 0.75) +
  labs(x = "Region", y = "Predicted aerobe relative abundance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12))



#### 4. Application ####
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
#saveRDS(results, "data/results_soil.rds")

# Analyze
results <- readRDS("data/results_soil.rds") %>%
  separate(filename, into = c("sampleID", "Junk1", "Junk2", "Junk3"),
           sep = "_") %>%
  dplyr::select(-Junk1, -Junk2, -Junk3)
range(results$Per_aerobe)
hist(results$Per_aerobe)
mean(results$Per_aerobe) # 96.3
median(results$Per_aerobe) # 98.7
sd(results$Per_aerobe) # 5.9
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
                dtpa_zinc, exc_aluminium, exc_calcium, exc_magnesium, exc_potassium,
                latitude, longitude, nitrate_nitrogen, organic_carbon, ph, 
                phosphorus_colwell, silt, water_content, bio1, bio12, AI)
d_env2 <- d_env %>%
  dplyr::select(-latitude, -longitude) %>%
  scale() %>%
  as.data.frame() %>%
  set_names(c("Clay", "Conductivity", "Cu", "Mn", "Zn", "Al", "Ca", "Mg", "K",
              "NO3", "C", "pH", "P", "Silt", "H2O", "Temp.", "Precip.", "Aridity")) %>%
  mutate(ClaySilt = Clay + Silt)
X <- d_env2 %>%
  dplyr::select(ClaySilt, Conductivity, NO3, C, pH, P, H2O, Temp., Precip., Aridity) %>%
  filter_at(vars(1:10), all_vars(!is.na(.)))
y <- aus %>%
  filter(sampleID %in% rownames(X)) %>%
  dplyr::select(Per_aerobe)
Xy <- cbind(X, y)
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
bestM # P, Precip.
# Estimate Std. Error    t value      Pr(>|t|)
# (Intercept) 96.4930205  0.4520363 213.462992 3.491420e-187
# P            1.3331709  0.8069191   1.652174  1.006048e-01
# Precip.      0.6272267  0.4347038   1.442883  1.511524e-01
bestM$BestModels

# Linear Models
aus_uni <- aus %>%
  mutate(ClaySilt = clay + silt) %>%
  dplyr::select(Per_aerobe, ClaySilt, conductivity, nitrate_nitrogen, organic_carbon, 
                ph, phosphorus_colwell, water_content, bio1, bio12,
                AI) %>%
  set_names(c("Per_aerobe", "ClaySilt", "Conductivity", "NO3", "C", "pH", "P", 
              "H2O", "Temp.", "Precip.", "Aridity")) %>%
  dplyr::select(ClaySilt, Conductivity, NO3, C, pH, P, H2O, Temp., 
                Precip., Aridity, Per_aerobe)
models <- as.data.frame(matrix(NA, 10, 3)) %>%
  set_names("Variable", "R", "p")
for (i in 1:10) {
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



#### _Panama ####
# Analyze all Panama metagenomes 



#### _NEON ####
# Use 2023 metagenomes and metadata



# End Script