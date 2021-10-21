# Interactive script to generate the wzx gene coverage maps for the NGS Unit's manuscripts
# comparing coverage between the various Illumina library prep kits. Both manuscripts are 
# covered here.
# Note R uses different file path syntax

# load packages
library(tidyverse)
library(egg)

# read in mpileup data Taylor generated before leave.
# everyone agreed to use the "proper pairs" data sets.
# grabbing the table first, then adding headers
mpileups_O130 <- read_tsv(file = paste("..", "Taylors_materials", "samtoools_mpileup_proper_O130.tsv", sep = "//"),
                     col_names = FALSE,
                     trim_ws = TRUE,
                     skip = 10,
                     skip_empty_rows = TRUE)
# now get the header rows and munge into column names.
O130Header <- read_lines(file = paste("..", "Taylors_materials", "samtoools_mpileup_proper_O130.tsv", sep = "//"),
                         n_max = 10)
O130Header <- lapply(O130Header, function(x) str_extract(x, "2015EL-1671a-.*0\\."))
O130Header <- lapply(O130Header, function(x) str_replace(x, "\\.", ""))
colnames(mpileups_O130) <- c("position", O130Header)
# repeat for second file. A less lazy person might write a single function, but there's only 2 files.
mpileups_O181 <- read_tsv(file = paste("..", "Taylors_materials", "samtoools_mpileup_proper_O181.tsv", sep = "//"),
                          col_names = FALSE,
                          trim_ws = TRUE,
                          skip = 10,
                          skip_empty_rows = TRUE)
O181Header <- read_lines(file = paste("..", "Taylors_materials", "samtoools_mpileup_proper_O181.tsv", sep = "//"),
                         n_max = 10)
O181Header <- lapply(O181Header, function(x) str_extract(x, "2015C-3865-.*0\\."))
O181Header <- lapply(O181Header, function(x) str_replace(x, "\\.", ""))
colnames(mpileups_O181) <- c("position", O181Header)

# combine tables
mpileups <- left_join(mpileups_O130, mpileups_O181, by = "position")

# these column names are terrible and have way too much data embedded. Let's make this tidy.
# first turn into a 3 column table with SeqID, Coverage, and Position.
mpileups <- gather(mpileups, SeqID, Coverage, "2015EL-1671a-M347-18-009_Flex_300":"2015C-3865-M3235-16-039_XT_500")
# now split SeqID into multiple columns
mpileups <- mutate(mpileups, 
                   Isolate = str_replace(str_extract(SeqID, "^.*-M"), "-M", ""),
                   LibKit = str_extract(SeqID, "(XT|NEB|KAPA|Qia|Flex)"),
                   Cycles = str_replace(str_extract(SeqID, "_(300|500)"), "_", ""))
# set correct data types
mpileups <- mutate(mpileups,
                   position = as.numeric(position),
                   SeqID = as.factor(SeqID),
                   Coverage = as.numeric(Coverage),
                   Isolate = as.factor(Isolate),
                   LibKit = as.factor(LibKit),
                   Cycles = as.factor(Cycles))

# now we plot! But first, we split up the data for different manuscripts
angelaData <- filter(mpileups, LibKit == "Flex" | LibKit == "XT")
angelaData$LibKit <- droplevels(angelaData$LibKit)
angelaData$LibKit <- recode(angelaData$LibKit,
                            Flex = "DNA Prep",
                            XT = "Nextera XT")
jennyData <- filter(mpileups, LibKit != "XT")
jennyData$LibKit <- droplevels(jennyData$LibKit)
jennyData$LibKit <- recode(jennyData$LibKit, 
                          Flex = "Illumina DNA Prep",
                          KAPA = "Roche KAPA HyperPlus",
                          NEB = "NEBNext Ultra II FS",
                          Qia = "Qiagen QIAseq FX")

angelaPlot <- ggplot(angelaData, aes(x = position, y = Coverage)) +
              geom_line(aes(linetype = Cycles, colour = Isolate)) +
              facet_wrap(vars(LibKit), nrow = 2) +
              xlab("wzx Alignment Position") +
              scale_color_manual(values = c("#e66101", "#5e3c99")) +
              scale_x_continuous(breaks = seq(0, 1300, 100)) +
              geom_hline(yintercept = 40, size = 1) +
              theme(strip.text.x = element_text(face = "bold"), axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
angelaPlot <- tag_facet(angelaPlot) +
              theme(strip.text = element_text(), strip.background = element_rect())
angelaPlot

# in panels by lib kit
jennyPlot <- ggplot(jennyData, aes(x = position, y = Coverage)) +
             geom_line(aes(linetype = Cycles, colour = Isolate)) +
             facet_wrap(vars(LibKit), nrow = 2, ncol = 2, scales = "free_y") +
             xlab("wzx Alignment Position") +
             scale_color_manual(values = c("#e66101", "#5e3c99")) +
             scale_x_continuous(breaks = seq(0, 1300, 100)) +
             geom_hline(yintercept = 40, size = 1)
jennyPlot

# in panels by cycle and isolate
jennyPlot <- ggplot(jennyData, aes(x = position, y = Coverage)) +
             geom_line(aes(colour = LibKit)) +
             facet_grid(vars(Isolate), vars(Cycles), scales = "free_y") +
             xlab("wzx Alignment Position") +
             scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
                                name = "Library Kit") +
             scale_x_continuous(breaks = seq(0, 1300, 100)) +
             geom_hline(yintercept = 40, size = 1)
jennyPlot <- tag_facet(jennyPlot) +
             theme(strip.text = element_text(), strip.background = element_rect())
jennyPlot


# export to pdf
pdf("wzxCoverage_Angela.pdf",
    width = 7,
    height = 7)
angelaPlot
dev.off()

pdf("wzxCoverage_Jenny.pdf",
    width = 14,
    height = 7)
jennyPlot
dev.off()