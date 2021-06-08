# Interactive script to generate the wzx gene coverage maps for the NGS Unit's manuscripts
# comparing coverage between the various Illumina library prep kits. Both manuscripts are 
# covered here.
# Note R uses different file path syntax

# load packages
library(tidyverse)

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
