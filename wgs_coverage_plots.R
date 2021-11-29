# Interactive script to generate the WGS coverage maps for the NGS Unit's manuscripts
# comparing coverage between the various Illumina library prep kits.
# Note R uses different file path syntax for Windows and Linux.

# load packages
library(tidyverse)
library(egg)
library(runner)

# read in mpileup data Taylor generated.
# everyone agreed to use the "proper pairs" data sets.
# grabbing the table first, then adding headers
mpileups_2015C_3865 <- read_tsv(file = file.path("..", "Taylors_materials", "samtoools_mpileup_proper_2015C-3865.tsv"),
                          col_names = FALSE,
                          trim_ws = TRUE,
                          skip_empty_rows = TRUE)

# no column names in this file, so clean it up and name them manually. 
# Column names extracted from bam-inputs.txt files in Taylor's directory.
mpileups_2015C_3865 <- dplyr::select(mpileups_2015C_3865, X2, X4, X7, X10, X13)
colnames(mpileups_2015C_3865) <- c("Position", "2015C-3865-M3235-16-039", "2015C-3865-M3235-18-002", "2015C-3865-M3235-18-007", "2015C-3865-M3235-19-034")

# lazy programmer repeats code for second file.
mpileups_2015EL_1671a <- read_tsv(file = file.path("..", "Taylors_materials", "samtoools_mpileup_proper_2015EL-1671a.tsv"),
                                col_names = FALSE,
                                trim_ws = TRUE,
                                skip_empty_rows = TRUE)

# no column names in this file, so clean it up and name them manually. 
# Column names extracted from bam-inputs.txt files in Taylor's directory.
mpileups_2015EL_1671a <- dplyr::select(mpileups_2015EL_1671a, X2, X4, X7, X10, X13)
colnames(mpileups_2015EL_1671a) <- c("Position", "2015EL-1671a-M3235-16-001", "2015EL-1671a-M347-18-007", "2015EL-1671a-M347-18-009", "2015EL-1671a-M347-20-006")

# Let's make this tidy.
# first turn into a 3 column table with SeqID, Coverage, and Position, then add a column for the reference.
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

mpileups_2015C_3865 <- NULL
mpileups_2015EL_1671a <- NULL
gc()

# now add metadata for each row 
mpileupsWGS <- mutate(mpileupsWGS, 
                      LibKit = case_when(
                        SeqID == "M3235-16-039" | SeqID == "M3235-16-001" | SeqID == "M3235-19-034" | SeqID == "M347-20-006" ~ 'Nextera XT',
                        SeqID == "M3235-18-007" | SeqID == "M347-18-009" | SeqID == "M3235-18-002" | SeqID == "M347-18-007" ~ 'DNA Prep'),
                      Cycles = case_when(
                        SeqID == "M3235-16-039" | SeqID == "M3235-16-001" | SeqID == "M3235-18-002" | SeqID == "M347-18-007" ~ '500',
                        SeqID == "M3235-18-007" | SeqID == "M347-18-009" | SeqID == "M3235-19-034" | SeqID == "M347-20-006" ~ '300'))

mpileupsWGS <- mutate(mpileupsWGS, 
                      Position = as.numeric(Position),
                      Coverage = as.numeric(Coverage))

# need to take the mean of 1000 nt windows because can't plot full resolution data
mpileupWGS_grouped <- group_by(mpileupsWGS, SeqID)
mpileupsWGS <- NULL
gc()
# This step is crapping out - using >10GB of memory despite clearing everything I can
mean_mpileupWGS <- summarise(mpileupWGS_grouped, 
                          Mean_Cov = runner(idx = Position,
                                            x = Coverage,
                                            f = function(x) mean(x),
                                            k = 1000,
                                            lag = 0),
                          STDEV = runner(idx = Position,
                                         x = Coverage,
                                         f = function(x) sd(x),
                                         k = 1000,
                                         lag = 0))


# plot!
wgsPlot <- ggplot(mpileupsWGS, aes(x = Position, y = Coverage)) +
           geom_line(aes(linetype = Cycles, colour = Isolate)) +
           facet_wrap(vars(LibKit), nrow = 2) +
           xlab("Alignment Position") +
           scale_color_manual(values = c("#e66101", "#5e3c99")) +
           scale_x_continuous() +
           geom_hline(yintercept = 40, size = 1) +
           theme(strip.text.x = element_text(face = "bold"), 
                 axis.title = element_text(face = "bold"), 
                 legend.title = element_text(face = "bold"))
wgsPlot <- tag_facet(wgsPlot) +
           theme(strip.text = element_text(), 
                 strip.background = element_rect())
wgsPlot

pdf("wgsCoverage_Angela.pdf",
    width = 14,
    height = 7)
wgsPlot
dev.off()