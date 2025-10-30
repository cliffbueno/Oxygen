# Parse Diamond Output, Calculate Ratio, Regress with Oxygen

#### 1. Setup ####
# Libraries
library(tidyverse)

# Working directory
setwd("~/Documents/GitHub/Oxygen/")

# The 20 Pfams
oxygen_pfams <- read.csv("data/Oxygen_pfams.csv")

# Gene IDs for all of the variants of each Pfam
map <- read.table("data/pfam_headers_table.txt", sep = "\t", header = TRUE, 
                  stringsAsFactors = FALSE, quote = "") %>%
  separate(Header, into = c("Header", "Junk"), sep = " ") %>%
  select(-Junk) %>%
  filter(!duplicated(Header))
table(map$Pfam) # Number of variants per Pfam

# Mean gene lengths for the 20 Pfams
pfam_gene_length <- read.delim("data/pfam_lengths.tsv")



#### 2. Parse Diamond ####
# Trial run on the first table 
d <- read.table("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_sim/sample0_01_trimmed_ann") %>%
  set_names(c("qseqid",	"sseqid",	"pident",	"length",	"qstart",	"qend",	
              "sstart",	"send",	"evalue",	"bitscore")) %>%
  filter(pident >= 60) %>%
  filter(evalue < 0.001) %>%
  filter(bitscore >= 50) %>%
  left_join(., map, by = c("sseqid" = "Header"))
table(d$Pfam) # PF05425 PF05721 are missing! (aerobic)

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
ratio



#### _5 mil ####
# Now loop through all 18 simulations of 5 million seq depth
setwd("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_sim")
files <- list.files()
length(files)
results <- as.data.frame(matrix(nrow = length(files), ncol = 2, NA)) %>%
  set_names(c("filename", "ratio")) %>%
  mutate(perc_aerobe = c(0, 0, 0, 100, 100, 100, 20, 20, 20, 40, 40, 40,
                         60, 60, 60, 80, 80, 80)) %>%
  mutate(replicate = c(rep(c(1, 2, 3), 6)))
  
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



#### _25 mil ####
# Now loop through all 18 simulations of 25 million seq depth
setwd("~/Desktop/Fierer/AEGIS/Oxygen/diamond_output_sim25/")
files <- list.files()
length(files)
results25 <- as.data.frame(matrix(nrow = length(files), ncol = 2, NA)) %>%
  set_names(c("filename", "ratio")) %>%
  mutate(perc_aerobe = c(0, 0, 0, 100, 100, 100, 20, 20, 20, 40, 40, 40,
                         60, 60, 60, 80, 80, 80)) %>%
  mutate(replicate = c(rep(c(1, 2, 3), 6)))

for (i in 1:length(files)) {
  
  # Add filename to the dataframe
  results25$filename[i] <- files[i]
  
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
  
  gene.hit.length.correction <- gene.hits %>%
    left_join(., pfam_gene_length, by = "Pfam") %>%
    mutate(RPK = total_count / (Gene.length/1000)) %>%
    left_join(., oxygen_pfams, by = "Pfam")
  
  # Now sum by aerobe indicator vs anaerobe indicator
  oxygen_rpk <- gene.hit.length.correction %>%
    group_by(Oxygen) %>%
    summarize(RPKsum = sum(RPK))
  
  # Calculate the ratio and add it to the dataframe
  results25$ratio[i] <- oxygen_rpk$RPKsum[1] / oxygen_rpk$RPKsum[2]
  
  # Status
  message("Done processing table ", i)
  
}

warnings()
setwd("~/Documents/GitHub/Oxygen/")



#### 3. Regress and Plot ####
# Relationship between % aerobes and the aerobic:anaerobic PFAM abundances

ggplot(results, aes(perc_aerobe, ratio)) +
  geom_point(size = 3, alpha = 0.75) +
  #geom_smooth(method = "lm", formula = y ~ I(x) + I(x^2)) +
  geom_smooth() +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  labs(x = "% aerobe genomes",
       y = "Aerobe:Anaerobe gene ratio (RPK)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(results25, aes(perc_aerobe, ratio)) +
  geom_point(size = 3, alpha = 0.75) +
  #geom_smooth(method = "lm", formula = y ~ I(x) + I(x^2)) +
  geom_smooth() +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  labs(x = "% aerobe genomes",
       y = "Aerobe:Anaerobe gene ratio (RPK)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

results$Depth <- "5 to 8 million reads"
results25$Depth <- "26 to 41 million reads"
c <- rbind(results, results25) %>%
  mutate(Depth = factor(Depth, levels = c("5 to 8 million reads", "26 to 41 million reads")))

pdf("FinalFigs/FigureS4.pdf", width = 7, height = 5)
ggplot(c, aes(perc_aerobe, ratio)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth() +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  labs(x = "% aerobe genomes",
       y = "Aerobe:Anaerobe gene ratio (RPK)") +
  facet_wrap(~ Depth) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))
dev.off()

png("FinalFigs/FigureS4.png", width = 7, height = 5, units = "in", res = 300)
ggplot(c, aes(perc_aerobe, ratio)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth() +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  labs(x = "% aerobe genomes",
       y = "Aerobe:Anaerobe gene ratio (RPK)") +
  facet_wrap(~ Depth) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))
dev.off()
range(c$ratio)

# Relationship is not linear and exponential doesn't quite fit right either
# Need GAM
# Use the data with greater sequencing depth because its more precise
library(mgcv)
gam_model <- gam(perc_aerobe ~ s(ratio), data = results25)
summary(gam_model)
plot(gam_model)
formula(gam_model)
new_data <- data.frame(ratio = c(5, 10, 15))  # or a sequence of values
predicted_y <- predict(gam_model, newdata = new_data)
predicted_y

# Save the GAM model
saveRDS(gam_model, "data/gam_model.rds")
