# Interactive script to generate the WGS coverage maps for the NGS Unit's manuscripts
# comparing coverage between the various Illumina library prep kits.
# Note R uses different file path syntax for Windows and Linux.

# load packages
library(tidyverse)
library(egg)
library(runner)
library(scales)

# read in mpileup data Taylor generated.
# everyone agreed to use the "proper pairs" data sets.
# grabbing the table first, then adding headers
mpileups_2015C_3865 <- read_tsv(file = file.path("..", "Taylors_materials", "samtools_mpileup", "samtoools_mpileup_proper_2015C-3865.tsv"),
                          col_names = FALSE,
                          trim_ws = TRUE,
                          skip_empty_rows = TRUE)

# no column names in this file, so clean it up and name them manually. 
# Column names extracted from bam-inputs.txt files in Taylor's directory.
mpileups_2015C_3865 <- dplyr::select(mpileups_2015C_3865, X2, X4, X7, X10, X13)
colnames(mpileups_2015C_3865) <- c("Position", "2015C-3865-M3235-16-039", "2015C-3865-M3235-18-002", "2015C-3865-M3235-18-007", "2015C-3865-M3235-19-034")

# lazy programmer repeats code for second file.
mpileups_2015EL_1671a <- read_tsv(file = file.path("..", "Taylors_materials", "samtools_mpileup", "samtoools_mpileup_proper_2015EL-1671a.tsv"),
                                col_names = FALSE,
                                trim_ws = TRUE,
                                skip_empty_rows = TRUE)

# no column names in this file, so clean it up and name them manually. 
# Column names extracted from bam-inputs.txt files in Taylor's directory.
mpileups_2015EL_1671a <- dplyr::select(mpileups_2015EL_1671a, X2, X4, X7, X10, X13)
colnames(mpileups_2015EL_1671a) <- c("Position", "2015EL-1671a-M3235-16-001", "2015EL-1671a-M347-18-007", "2015EL-1671a-M347-18-009", "2015EL-1671a-M347-20-006")

# Let's make this tidy.
# first turn into a 3 column table with SeqID, Coverage, and Position, then add a column for the isolate ID.
mpileups_2015C_3865 <- gather(mpileups_2015C_3865, 
                              SeqID, 
                              Coverage, 
                              "2015C-3865-M3235-16-039":"2015C-3865-M3235-19-034") %>%
                       mutate(Isolate = "2015C-3865", 
                              SeqID = str_replace(SeqID, "2015C-3865-", ""))
mpileups_2015EL_1671a <- gather(mpileups_2015EL_1671a, 
                                SeqID, 
                                Coverage, 
                                "2015EL-1671a-M3235-16-001":"2015EL-1671a-M347-20-006") %>%
                         mutate(Isolate = "2015EL-1671a", 
                                SeqID = str_replace(SeqID, "2015EL-1671a-", ""))

# combine 2 mpileup tables
mpileupsWGS <- rows_insert(mpileups_2015C_3865, 
                           mpileups_2015EL_1671a, 
                           by = c("SeqID", "Isolate", "Position"))

# these are big data sets, so clear out tables that you're done with
mpileups_2015C_3865 <- NULL
mpileups_2015EL_1671a <- NULL
gc()

# need to take the mean of 10,000 nt windows because can't plot full resolution data
mpileupWGS_grouped <- group_by(mpileupsWGS, SeqID)


# Runner introduces a bunch of garbage NAs between the end of the real mpileup data and the 5 million nt end
mean_mpileupWGS <- summarise(mpileupWGS_grouped, 
                          Iso = runner(at = seq(from = 10000, to = 5000000, by = 10000), # at provides the end point for the windows that runner is calculating on
                                             x = Isolate, # the value to use in the calculations
                                             f = function(x) x[1], # the calculation to do
                                             k = 10000, # the size of each window
                                             lag = 0), # don't offset
                          Pos_Window = runner(at = seq(from = 10000, to = 5000000, by = 10000),
                                              x = Position,
                                              f = function(x) max(x),
                                              k = 10000,
                                              lag = 0),
                          Mean_Cov = runner(at = seq(from = 10000, to = 5000000, by = 10000),
                                            x = Coverage,
                                            f = function(x) mean(x),
                                            k = 10000,
                                            lag = 0),
                          STDEV = runner(at = seq(from = 10000, to = 5000000, by = 10000),
                                         x = Coverage,
                                         f = function(x) sd(x),
                                         k = 10000,
                                         lag = 0))

# Remove garbage NA lines
mean_mpileupWGS <- filter(mean_mpileupWGS, 
                          is.na(Iso) == FALSE,
                          Pos_Window != -Inf,
                          is.na(STDEV) == FALSE)
# now add metadata for each row 
mean_mpileupWGS <- mutate(mean_mpileupWGS,
                    LibKit = case_when(
                              SeqID == "M3235-16-039" | SeqID == "M3235-16-001" | SeqID == "M3235-19-034" | SeqID == "M347-20-006" ~ 'Nextera XT',
                              SeqID == "M3235-18-007" | SeqID == "M347-18-009" | SeqID == "M3235-18-002" | SeqID == "M347-18-007" ~ 'DNA Prep'),
                    Cycles = case_when(
                            SeqID == "M3235-16-039" | SeqID == "M3235-16-001" | SeqID == "M3235-18-002" | SeqID == "M347-18-007" ~ '500',
                            SeqID == "M3235-18-007" | SeqID == "M347-18-009" | SeqID == "M3235-19-034" | SeqID == "M347-20-006" ~ '300'))

mean_mpileupWGS <- mutate(mean_mpileupWGS,
                          FacetVar = paste(LibKit, Cycles, "Cycles"))

# plot!

# How am I ending up with an NA in isolate???
#Don't forget to truncate yaxis
wgsPlot_Prep <- ggplot(filter(mean_mpileupWGS, LibKit == "DNA Prep"), aes(x = Pos_Window, y = Mean_Cov)) +
           geom_ribbon(aes(ymin = Mean_Cov - STDEV,
                           ymax = Mean_Cov + STDEV,
                           fill = Iso),
                       alpha = 0.5) +
           geom_line(aes(x = Pos_Window, y = Mean_Cov, linetype = Iso)) +
           facet_wrap(vars(Cycles), nrow = 2) +
           labs(x = "Alignment Position", y = "Mean Coverage / 10k Nucleotides", linetype = "Isolate", fill = "Isolate") +
           scale_fill_discrete(type = c("#e66101", "#5e3c99")) +
           scale_y_continuous(limits = c(15, 150), breaks = function(x) seq(from = 20, to = 140, by = 20), minor_breaks = NULL) +
           scale_x_continuous(limits = c(0, 4750000), labels = comma, expand = c(0, 0)) +
           geom_hline(yintercept = 40, size = 1) +
           theme(strip.text.x = element_text(face = "bold", size = 12), 
                 axis.title = element_text(face = "bold", size = 15), 
                 legend.title = element_text(face = "bold", size = 15),
                 axis.text = element_text(size = 12, colour = "black"),
                 legend.text = element_text(size = 12),
                 legend.key.size = unit(1, "cm"),
                 #legend.direction = "horizontal",
                 legend.position = c(0.095, 0.88))
wgsPlot_Prep <- tag_facet(wgsPlot_Prep) +
           theme(strip.text = element_text(), 
                 strip.background = element_rect())
wgsPlot_Prep

pdf("prep_wgsCoverage_Angela.pdf",
    width = 14,
    height = 7)
wgsPlot_Prep
dev.off()

wgsPlot_XT <- ggplot(filter(mean_mpileupWGS, LibKit == "Nextera XT"), aes(x = Pos_Window, y = Mean_Cov)) +
  geom_ribbon(aes(ymin = Mean_Cov - STDEV,
                  ymax = Mean_Cov + STDEV,
                  fill = Iso),
              alpha = 0.5) +
  geom_line(aes(x = Pos_Window, y = Mean_Cov, linetype = Iso)) +
  facet_wrap(vars(Cycles), nrow = 2) +
  labs(x = "Alignment Position", y = "Mean Coverage / 10k Nucleotides", linetype = "Isolate", fill = "Isolate") +
  scale_fill_discrete(type = c("#e66101", "#5e3c99")) +
  scale_y_continuous(limits = c(15, 180), breaks = function(x) seq(from = 20, to = 180, by = 20), minor_breaks = NULL) +
  scale_x_continuous(limits = c(0, 4750000), labels = comma, expand = c(0, 0)) +
  geom_hline(yintercept = 40, size = 1) +
  theme(strip.text.x = element_text(face = "bold", size = 12), 
        axis.title = element_text(face = "bold", size = 15), 
        legend.title = element_text(face = "bold", size = 15),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        legend.direction = "horizontal",
        legend.position = c(0.18, 0.94))
wgsPlot_XT <- tag_facet(wgsPlot_XT) +
  theme(strip.text = element_text(), 
        strip.background = element_rect())
wgsPlot_XT

pdf("xt_wgsCoverage_Angela.pdf",
    width = 14,
    height = 7)
wgsPlot_XT
dev.off()

# Joint plot
wgsPlot <- ggplot(filter(mean_mpileupWGS, Pos_Window < 4750001), aes(x = Pos_Window, y = Mean_Cov)) +
  geom_ribbon(aes(ymin = Mean_Cov - STDEV,
                  ymax = Mean_Cov + STDEV,
                  fill = Iso),
              alpha = 0.5) +
  geom_line(aes(x = Pos_Window, y = Mean_Cov, linetype = Iso)) +
  facet_wrap(vars(FacetVar), nrow = 4, scales = "free_y") +
  labs(x = "Alignment Position", y = "Mean Coverage / 10k Nucleotides", linetype = "Isolate", fill = "Isolate") +
  scale_fill_discrete(type = c("#e66101", "#5e3c99")) +
  #scale_y_continuous(limits = c(15, 180), breaks = function(x) seq(from = 20, to = 180, by = 20), minor_breaks = NULL) +
  scale_x_continuous(limits = c(0, 4750000), labels = comma, expand = c(0, 0)) +
  geom_hline(yintercept = 40, size = 1) +
  theme(strip.text.x = element_text(face = "bold", size = 12), 
        axis.title = element_text(face = "bold", size = 12), 
        legend.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.position = c(0.12, 0.98))
wgsPlot <- tag_facet(wgsPlot) +
  theme(strip.text = element_text(), 
        strip.background = element_rect())
wgsPlot

pdf("wgsCoverage_Angela.pdf",
    width = 21,
    height = 14)
wgsPlot
dev.off()

# The journal insisted. EPS not an option because you can't do transparency.
ggsave(filename = "wgsCoverage_Angela.tiff", 
       plot = wgsPlot, 
       device = "tiff",
       dpi = 600,
       width = 21,
       height = 14, 
       units = "in")
