library(data.table)
library(ggplot2)
library(grid)
library(dplyr)


#### fig 1 ----

##### a) ld decay ----
script <- "scripts/plot_ld_decay.R"
arg1 <- "analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.window-1000kb.filter-0.25.ld.gz"
arg2 <- "../../../publications/2021_simulations/fig1a"
opts <- "--unequal-windows"

system(paste("Rscript", script, arg1, arg2, opts))


##### b) ld distribution SNP-SV -----
ld_results <- data.frame()
for (chr in 1:10) {
  ld_file_chr <- paste0("analysis/ld/window-1kb_filter-0.25/ld_usda_rils_snp-sv_only.chr", chr, ".window-1kb.filter-0.25.highest-ld.ld")
  ld_file_chr <- fread(ld_file_chr, header = TRUE, data.table = FALSE)
  ld_results <- rbind(ld_results, ld_file_chr)
}
rm(ld_file_chr, chr)

# ggplot(ld_results, aes(x = R2)) +
#   geom_histogram(fill = "#900721", binwidth = 0.01) +
#   labs(title = paste0("LD between SVs and SNPs"),
#        x = bquote("LD"~(r^2)),
#        y = "Count") +
#   coord_cartesian(xlim = c(0, 1)) +
#   theme(title = element_text(size = 15),
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 15))

plot_data <- ggplot_build(ggplot(ld_results, aes(x = R2)) + geom_histogram(binwidth = 0.01))
plot_data <- data.frame(plot_data$data)[c("x", "ymin", "ymax")]
plot_data[NROW(plot_data), "ymin"] <- 2000
plot_data <- rbind(plot_data, data.frame(x = plot_data[NROW(plot_data), "x"], ymin = 0, ymax = 100))
plot_data[plot_data$ymin > 100, "y_axis_break"] <- "above"
plot_data[plot_data$ymin <= 100, "y_axis_break"] <- "below"

dist_plot <- ggplot(plot_data, aes(x = x, y = ymin, xend = x, yend = ymax)) +
  geom_segment(size = 5) +
  facet_grid(y_axis_break ~ ., scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = bquote("LD"~(r^2)),
       y = "Count") +
  theme_minimal() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        strip.text.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line()) + 
  annotation_custom(grob = linesGrob(gp = gpar(lwd = 2)), xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = +Inf) +
  annotation_custom(grob = linesGrob(gp = gpar(lwd = 2)), xmin = +Inf, xmax = +Inf, ymin = -Inf, ymax = +Inf) +
  annotation_custom(grob = linesGrob(gp = gpar(lwd = 2)), xmin = -Inf, xmax = +Inf, ymin = -5, ymax = -5) +
  annotation_custom(grob = linesGrob(gp = gpar(lwd = 2)), xmin = -Inf, xmax = +Inf, ymin = 4231, ymax = 4231)

ggsave(filename = "../../../publications/2021_simulations/fig1b.png", dist_plot, device = "png")

# inspired by:
#   https://stackoverflow.com/a/65242880
#   https://www.j4s8.de/post/2018-01-15-broken-axis-with-ggplot2/


# plot MAF per LD -- way 1
ld_results %>%
  mutate(sv_maf = if_else(grepl("^del|ins|inv|dup|tra", SNP_A, perl = TRUE),
                          true = MAF_A, false = MAF_B),
         ld_category = case_when(R2 >= 0.8 ~ "High LD",
                                 R2 < 0.8 & R2 >= 0.5 ~ "Moderate LD",
                                 R2 < 0.5 ~ "Low LD")) %>%
  ggplot(aes(x = sv_maf)) +
  geom_density(aes(color = ld_category))

# plot MAF per LD -- way 2 (rils)
site_summary <- fread("analysis/ld/window-1kb_filter-0.25/SVs_with_LD_info_SiteSummary.txt",
                      header = TRUE, data.table = FALSE)
site_summary <- site_summary %>% 
  select(`Site Name`, `Major Allele`, `Major Allele Frequency`, `Minor Allele Frequency`) %>%
  mutate(sv_frequency = if_else(condition = `Major Allele` == "T", 
                                true = `Major Allele Frequency`,
                                false = `Minor Allele Frequency`))

ld_summary <- ld_results %>% 
  mutate(sv = if_else(grepl("^del|ins|inv|dup|tra", SNP_A, perl = TRUE),
                      true = SNP_A, false = SNP_B),
         ld_category = case_when(R2 >= 0.8 ~ "High LD",
                                 R2 < 0.8 & R2 >= 0.5 ~ "Moderate LD",
                                 R2 < 0.5 ~ "Low LD")) %>% 
  select(sv, R2, ld_category)
ld_summary <- ld_summary[!duplicated(ld_summary$sv), ]

if (all(ld_summary$sv == site_summary$`Site Name`)) {
  ld_summary <- cbind(ld_summary, sv_frequency = site_summary$sv_frequency)
  # plot
  ld_summary %>%
    ggplot(aes(x = sv_frequency, color = ld_category)) +
    geom_density()
}

# plot MAF per LD -- way 3 (parents)
site_summary <- fread("analysis/ld/window-1kb_filter-0.25/SVs_with_LD_info_SiteSummary.parents-only.txt",
                      header = TRUE, data.table = FALSE)
site_summary <- site_summary %>% 
  select(`Site Name`, `Major Allele`, `Major Allele Frequency`, `Minor Allele Frequency`) %>%
  mutate(sv_frequency = if_else(condition = `Major Allele` == "T", 
                                true = `Major Allele Frequency`,
                                false = `Minor Allele Frequency`))

ld_summary <- ld_results %>% 
  mutate(sv = if_else(grepl("^del|ins|inv|dup|tra", SNP_A, perl = TRUE),
                      true = SNP_A, false = SNP_B),
         ld_category = case_when(R2 >= 0.8 ~ "High LD",
                                 R2 < 0.8 & R2 >= 0.5 ~ "Moderate LD",
                                 R2 < 0.5 ~ "Low LD")) %>% 
  select(sv, R2, ld_category)
ld_summary <- ld_summary[!duplicated(ld_summary$sv), ]

if (all(ld_summary$sv == site_summary$`Site Name`)) {
  ld_summary <- cbind(ld_summary, sv_frequency = site_summary$sv_frequency)
  # plot
  ld_summary %>%
    ggplot(aes(x = sv_frequency, color = ld_category)) +
    geom_density()
}







#### supp fig X ----

# projection rate and accuracy

proj_summary_sv <- fread("analysis/projection_svs-snps/projection_svs_summary.txt",
                         header = TRUE, data.table = FALSE)
proj_summary_snp <- fread("analysis/projection_svs-snps/projection_reseq-snps_summary.txt",
                          header = TRUE, data.table = FALSE)
proj_summary <- rbind(data.frame(proj_summary_sv[, c(1, 5, 7)], marker = "SVs"),
                      data.frame(proj_summary_snp[, c(1, 3, 4)], marker = "SNPs"))
rm(proj_summary_sv, proj_summary_snp)

ggplot(proj_summary, aes(x = family, y = avg_percent_projected, fill = marker)) +
  geom_col(position = "dodge", width = 0.6) +
  labs(x = "Family", y = "Projection rate (%)", fill = "Marker type") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 1))+
  scale_x_discrete(labels = function(cross) gsub("x", "\nx\n", cross)) +
  scale_fill_manual(values = c("grey60", "grey20")) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        title = element_text(size = 20),
        text = element_text(size = 15))

ggplot(proj_summary, aes(x = family, y = proj_accuracy, fill = marker)) +
  geom_col(position = "dodge", width = 0.6) +
  labs(x = "Family", y = "Projection accuracy (%)", fill = "Marker type") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 1))+
  scale_x_discrete(labels = function(cross) gsub("x", "\nx\n", cross)) +
  scale_fill_manual(values = c("grey60", "grey20")) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        title = element_text(size = 20),
        text = element_text(size = 15))
