library(data.table)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

infile <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/prediction_all_markers/GEBVs_ranks.txt"
# infile <- "analysis/test_prediction/multi_env/no_gxe/100qtns_SVs_equal_0.5h2_pop1/prediction_all_markers/GEBVs_ranks.txt"

ranks <- fread(infile, header = TRUE, data.table = FALSE)
ranks$family <- sapply(ranks$genotype, function(x) gsub("*", "x", unlist(strsplit(x, split = "-"))[1], fixed = TRUE))
ranks <- ranks[order(ranks$env1), ]
ranks$genotype <- factor(ranks$genotype, levels = ranks$genotype)
ranks <- pivot_longer(ranks, cols = -c(genotype, family), names_to = "env", values_to = "rank")
ranks$env <- factor(ranks$env)

# # line plot with changes over environments
# ggplot(ranks) + 
#   geom_line(aes(x = env, y = rank, group = genotype, color = family), alpha = 0.5)


# heatmap showing color change from lowest rank to highest rank
ggplot(ranks) +
  geom_tile(aes(x = env, y = genotype, fill = rank)) +
  scale_fill_distiller(type = "div", palette = "RdBu", aesthetics = "fill")
