# Run 66 total simulations and make ternary plot
# Assess the role of facultative bacteria
# This was done in response to a reviewer comment and the manuscript has been updated
# This helps interpret the different scenarios that can lead to a given OxyMetaG output
# The OxyMetaG output is now changed from "Per_aerobe" to "Oxygen_status"
# This is more accurate becasue 100% aerobe could have actually had some facultative taxa
# Oxygen status is just a unitless value scaled from 0 to 100 that gives you an idea
# of the oxygen availability in your sample, with 100 being high O2 and 0 being low O2
# Based on the analysis below, there are a few key points
# OS > 63 means you have at least some aerobes. 63 is the max with 100% facultative
# OS = 0 means you have all anaerobes
# OS = 100 means you have all aerobe
# Intermediate OS, could be various combos of all 3 - see ternary plot!



#### 1. Setup ####
# Libraries
library(tidyverse)
library(ggtern)
library(plotly)
library(interp)
library(geometry)

# Working directory (GitHub repo)
setwd("~/Documents/GitHub/Oxygen/")

# Data
# Note: There were 3 batches done
# Together they are all 66 possible combos at 10% increments
d1 <- read.delim("data/per_aerobe_predictions_n5.tsv")
oldd2 <- read.delim("data/per_aerobe_predictions_fac_n6.tsv") %>%
  mutate(Batch = as.integer(gsub("batch_|_trimmed", "", SampleID))) %>%
  arrange(Batch) %>%
  mutate(aerobic = c(50,50,50,40,40,40,30,30,30,20,20,20,10,10,10,0,0,0),
         anaerobic = c(0,0,0,20,20,20,40,40,40,60,60,60,80,80,80,100,100,100),
         facultative = c(50,50,50,40,40,40,30,30,30,20,20,20,10,10,10,0,0,0)) %>%
  group_by(aerobic, anaerobic, facultative) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(Batch)
d2 <- read.delim("data/per_aerobe_predictions_fac_n6_redo.tsv") %>%
  dplyr::select(-Bad) %>%
  mutate(Batch = as.integer(gsub("batch_|_trimmed", "", SampleID))) %>%
  arrange(Batch) %>%
  mutate(aerobic = c(50,40,30,20,10,0),
         anaerobic = c(0,20,40,60,80,100),
         facultative = c(50,40,30,20,10,0))
d3_info <- read.csv("data/sim_info_n55.csv") %>%
  mutate(Batch = as.integer(rownames(.)))
d3 <- read.delim("data/per_aerobe_predictions_fac_n55.tsv") %>%
  mutate(Batch = as.integer(gsub("batch_|_trimmed", "", SampleID))) %>%
  arrange(Batch) %>%
  left_join(., d3_info, by = "Batch")
d <- rbind(d1, d2, d3) %>%
  arrange(aerobic, anaerobic, facultative)



#### 2. Plot ####
# Better names
d <- d %>%
  rename(aerobes = aerobic,
         anaerobes = anaerobic)

# Basic Ternary plot
ggplot(data = d, aes(x = aerobes, y = anaerobes, z = facultative)) +
  geom_mask() +
  geom_point(aes(color = Per_aerobe), size = 3) +
  scale_color_viridis_c(name = "Oxygen status") +
  coord_tern() +
  scale_T_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_L_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_R_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_minimal() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85,0.7),
        tern.axis.title.L = element_text(hjust = 0.2),
        tern.axis.title.R = element_text(hjust = 0.65))

# ggtern examples for interpolation and text
data(Feldspar)
ggtern(Feldspar,aes(Ab,An,Or,value=T.C)) + 
  stat_interpolate_tern(geom="polygon",
                        formula=value~x+y,
                        method=lm,n=100,
                        breaks=seq(0,1000,by=100),
                        aes(fill=..level..),expand=1) +
  geom_point()

ggtern() + 
  annotate(geom  = 'text',
           x     = c(0.5,1/3,0.0),
           y     = c(0.5,1/3,0.0),
           z     = c(0.0,1/3,1.0),
           angle = c(0,30,60),
           vjust = c(1.5,0.5,-0.5),
           label = paste("Point",c("A","B","C")),
           color = c("green","red",'blue')) +
  theme_dark() + 
  theme_nomask()

# Contour with text
pdf("InitialFigs/OxyMetaG_Ternary_Contour.pdf", width = 7, height = 5)
ggtern(d, aes(x = aerobes, y = anaerobes, z = facultative, value = Per_aerobe)) + 
  geom_mask() +
  stat_interpolate_tern(geom = "polygon",
                        formula = value ~ x + y,
                        method = lm, n = 100,
                        breaks = seq(0, 100, by = 0.5),
                        aes(fill=..level..), expand = 2) +
  annotate(geom = "text",
           x = d$aerobes, y = d$anaerobes, z = d$facultative,
           label = round(d$Per_aerobe, 0), size = 3, color = "grey40") +
  scale_fill_viridis_c(name = "Oxygen status") +
  scale_T_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_L_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_R_continuous(breaks = seq(0, 1, by = 0.1)) +
  guides(T = "none", L = "none", R = "none") +
  theme_minimal() +
  theme_hidegrid() +
  theme_showarrows() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9,0.65),
        tern.axis.title.L = element_text(hjust = 0.2),
        tern.axis.title.R = element_text(hjust = 0.65),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        tern.panel.grid.major = element_blank(),
        tern.panel.grid.minor = element_blank(),
        tern.axis.line = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(colour = "white"))
dev.off()

# The interpolation does not fill the triangle and is inaccurate by the axes
# Need to try to interpolate outside of plot, then plot
interp_ternary <- function(df, n = 150) {
  tri <- delaunayn(cbind(df$aerobes, df$anaerobes))
  
  grid <- expand.grid(aerobes = seq(0, 100, length.out = n),
                      anaerobes = seq(0, 100, length.out = n))
  
  grid$facultative <- 100 - grid$aerobes - grid$anaerobes
  
  grid <- grid[grid$facultative >= 0 & grid$facultative <= 100, ]
  
  simplex <- tsearchn(cbind(df$aerobes, df$anaerobes), tri,
                      cbind(grid$aerobes, grid$anaerobes))
  
  interp_vals <- rep(NA, nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    s <- simplex$idx[i]
    if (is.na(s)) next
    
    verts <- tri[s, ]
    bary <- simplex$p[i, ]
    
    interp_vals[i] <- sum(bary * df$Per_aerobe[verts])
  }
  
  grid$Per_aerobe <- interp_vals
  grid
}
interp_df <- interp_ternary(d, n = 200)
plot_ly(data = interp_df,
        type = "scatterternary",
        mode = "markers",
        a = ~aerobes,
        b = ~anaerobes,
        c = ~facultative,
        marker = list(color = ~Per_aerobe,
                      colorscale = "Viridis",
                      size = 3,
                      showscale = TRUE))
e <- d %>%
  mutate(across(c(aerobes, anaerobes, facultative), ~ ifelse(. == 0, 2, .)))

# Use this - write as new Figure 2!
g <- ggtern(interp_df, 
            aes(x = aerobes, y = anaerobes, z = facultative, value = Per_aerobe)) + 
  geom_mask() +
  geom_point(aes(color = Per_aerobe), size = 0.5) +
  annotate(geom = "text",
           x = e$aerobes, y = e$anaerobes, z = e$facultative,
           label = round(e$Per_aerobe, 0), size = 3, color = "white") +
  annotate("text", x = 1, y = 0, z = 0, label = "aerobes", vjust = 1.1, hjust = 1) +
  annotate("text", x = 0, y = 1, z = 0, label = "anaerobes", vjust = -0.9) +
  annotate("text", x = 0, y = 0, z = 1, label = "facultative", vjust = 1.1, hjust = 0) +
  scale_color_viridis_c(option = "B",
                        name = "Predicted\noxygen level (%)") +
  scale_T_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_L_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_R_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs("T" = "% anaerobe genomes",
       "L" = "% aerobe genomes",
       "R" = "% facultative genomes") +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme_minimal() +
  theme_hidegrid() +
  theme_showarrows() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.2,0.7),
        tern.axis.title.T = element_blank(),
        tern.axis.title.L = element_blank(),
        tern.axis.title.R = element_blank(),
        axis.ticks = element_line(colour = "white"),
        plot.margin = margin(0, 0, 0, 0))
g
pdf("FinalFigs/Figure2.pdf", width = 7, height = 5)
g
dev.off()
png("FinalFigs/Figure2.png", width = 7, height = 5, res = 300, units = "in")
g
dev.off()

# Or points with text
low <- d %>%
  filter(Per_aerobe < 30)
high <- d %>%
  filter(Per_aerobe >= 30)
pdf("InitialFigs/OxyMetaG_Ternary.pdf", width = 7, height = 5)
ggtern(d, aes(x = aerobes, y = anaerobes, z = facultative, value = Per_aerobe)) + 
  geom_mask() +
  geom_point(aes(color = Per_aerobe), size = 4.6, pch = 16) +
  annotate(geom = "text",
           x = d$aerobes, y = d$anaerobes, z = d$facultative, 
           label = round(d$Per_aerobe, 0), size = 2.25, colour = "white") +
  scale_color_viridis_c(name = "Oxygen status") +
  scale_T_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_L_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_R_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_minimal() +
  theme_showarrows() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9,0.65),
        legend.title = element_text(size = 12, hjust = 0.5),
        tern.axis.title.L = element_text(hjust = 0.2),
        tern.axis.title.R = element_text(hjust = 0.65))
dev.off()

# Interactive - plot_ly
plot_ly(data = d,
        type = "scatterternary",
        mode = "markers",
        a = ~anaerobes,
        b = ~aerobes,
        c = ~facultative,
        marker = list(color = ~Per_aerobe,
                      colorscale = "Viridis",
                      showscale = TRUE),
        hovertemplate = paste("Anerobes: %{a}<br>",
                              "Aerobes: %{b}<br>",
                              "Facultative: %{c}<br>",
                              "Predicted: %{marker.color:.2f}<extra></extra>")) %>%
  layout(ternary = list(aaxis = list(title = "Anaerobes"),
                        baxis = list(title = "Aerobes"),
                        caxis = list(title = "Facultative")))

# Binary plots for inspection, keeping the 3rd variable consistent
fac0 <- d %>%
  filter(facultative == 0)
ggplot(fac0, aes(aerobic, Per_aerobe)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 5, shape = 16) +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
  scale_y_continuous(limits = c(0,100),
                     breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
  labs(x = "% Aerobes",
       y = "Oxygen status") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank())
aer0 <- d %>%
  filter(aerobic == 0)
ggplot(aer0, aes(facultative, Per_aerobe)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 5, shape = 16) +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
  scale_y_continuous(limits = c(0,100),
                     breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
  labs(x = "% Facultative",
       y = "Oxygen status") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank())
an0 <- d %>%
  filter(anaerobic == 0)
ggplot(an0, aes(aerobic, Per_aerobe)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 5, shape = 16) +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
  scale_y_continuous(limits = c(0,100),
                     breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
  labs(x = "% Aerobic",
       y = "Oxygen status") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank())
