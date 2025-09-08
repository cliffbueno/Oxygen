# Calculate % aerobe for modern metagenomes
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

# Mean gene lengths for the 20 Pfams
pfam_gene_length <- read.delim("pfam_prot_rev/pfam_mean_lengths.tsv") %>%
  separate(Pfam, sep = "-", into = c("Junk1", "Junk2", "Pfam")) %>%
  mutate(Gene.length = MeanLength * 3) %>%
  dplyr::select(Pfam, Gene.length)

# Model (from simulations)
gam_model <- readRDS("gam_model.rds")

# Tran et al. 2021 Oxygen
# (metaG don't have O2, but can infer from depth)
tran <- read.csv("SampleInfo/TranDO.csv")
ggplot(tran, aes(Depth, DO)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = 120, color = "red") +
  geom_hline(yintercept = 12.5, color = "red") +
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
# Also make a categorical plot



#### 2. Parse Diamond ####
#### _Test ####
# Trial run on the first table 
d <- read.table("diamond_output_modern/19485_R1_ann.tsv") %>%
  set_names(c("qseqid",	"sseqid",	"pident",	"length",	"qstart",	"qend",	
              "sstart",	"send",	"evalue",	"bitscore")) %>%
  filter(pident >= 60) %>%
  filter(evalue < 0.001) %>%
  filter(bitscore >= 50) %>%
  mutate(Gene.length = abs(send - sstart)) %>%
  left_join(., map, by = c("sseqid" = "Header"))
hist(d$Gene.length)
d$sseqid[1244] # 2 pfams (cytochrome c)
table(d$Pfam) # PF05425 PF05721 are missing! (aerobic)

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
predicted_y # 79.56% aerobe

# Test more individual samples
# Need to diagnose outliers and unexpected results
# One question is how many of the 20 genes are present in certain samples

# Hot Spring Sediment
t <- read.table("diamond_output_modern/SRR25522566_R1_ann.tsv") %>%
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
setwd("diamond_output_modern/")
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
meta <- read_xlsx("ModernMetaG.xlsx")
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
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()



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



#### 4. Application ####
# Use the tool on a big dataset to answer a question
# Which variables (precip, aridity, moisture, carbon, texture) are most closely
# associated with predicted % aerobe relative abundance?
# And - is this consistent across 3 independent datasets?
# Test across Australia, Panama, NEON metagenomes and metadata


#### _Australia ####
# Analyze all 331 natural surface Australian metagenomes

#### _Panama ####
# Analyze all Panama metagenomes 

#### _NEON ####
# Use 2023 metagenomes and metadata


# End Script