library(data.table)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)

usage <- function() {
  cat("
description: plot karyotypes

usage: Rscript karyoplots.R [chr_info] [centromere_info] [infile_name] [outfile_name] [...]

positional arguments:
  chr_info                tab-delimited file with chromosome coordinates (format: chr, length)
  centromere_info         tab-delimited file with centromere coordinates (format: chr, start_pos, end_pos)
  infile_name             tab-delimited input file (format: id, chr, pos, attribute)
  outfile_name            name of karyoplot (with '.pdf' extension)

optional argument:
  --help                  show this helpful message
  --line-width=VALUE      control line width (default: 6)
  --alpha=VALUE           control line transparency (default: 0.3)
  --add-legend            add legend to plot


credits: this code was adapted from Carles Hernandez-Ferrer's blog:
         https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/
         
"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(as.numeric(arg[2]))
  
}



#### command line options ----

# set default
line_width <- 6
alpha <- 0.3
add_legend <- FALSE

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 4) stop(usage(), "missing positional argument(s)")

if (length(args) > 4) {
  
  opt_args <- args[-1:-4]
  opt_args_allowed <- c("--line-width", "--alpha", "--add-legend")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")
  
  # change default based on the argument provided
  for (argument in opt_args_allowed) {
    if (any(grepl(argument, opt_args_requested))) {
      arg_name <- gsub("-", "_", gsub("--", "", argument))
      arg_value <- getArgValue(opt_args[grep(argument, opt_args)])
      assign(arg_name, arg_value)
    }
  }
  
}


#### plot karyotypes ----

# get arguments
chr_info <- args[1]
centromere_info <- args[2]
infile_name <- args[3]
outfile_name <- args[4]
# chr_info <- "data/B73_RefGen_V4_chrm_info.txt"
# centromere_info <- "data/centromeres_Schneider-2016-pnas_v4.bed"
# infile_name <- "karyoplot_toy.txt"
# outfile_name <- "karyoplot_toy.pdf"
# line_width <- 6
# alpha <- 0.3
# add_legend <- FALSE

# chromosomes
chrms <- fread(chr_info, header = TRUE, data.table = FALSE)
chrms <- data.frame(chr = chrms$chr, start_pos = 0, end_pos = chrms$length)

# centromeres
centros <- read.delim(centromere_info, sep = "\t", header = TRUE)
centros <- data.frame(chr = centros$chr, start_pos = centros$start_pos, end_pos = centros$end_pos)

# load input file
infile <- fread(infile_name, header = TRUE, data.table = FALSE)
colnames(infile) <- c("id", "chr", "pos", "att")

# make sure attributes are factors
infile[, "att"] <- as.factor(infile[, "att"])
att_levels <- levels(infile[, "att"])

# plot blank karyotypes 
plot <- ggplot() +
  geom_segment(data = chrms,
               aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
               lineend = "round", color = "Gainsboro", size = 5) +
  geom_segment(data = centros,
               aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
               lineend = "round", color = "DimGray", size = 5) +
  facet_grid(~chr, switch = "y") +
  scale_y_reverse(breaks = seq(0, 3.5e8, 0.50e8), labels = c(1, seq(50, 350, 50))) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(size = rel(1.1), color = "DimGray"),
        axis.text = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        strip.text.x = element_text(size = rel(2)),
        legend.position = "bottom") +
  labs(x = "Chromosomes",
       y = "Genomic positions (Mb)")

# get colors for each attribute
att_colors <- brewer.pal(n = length(att_levels), name = 'Dark2')

for (att in 1:length(att_levels)) {
  
  # filter data to have only specific attribute
  infile_filtered <- infile[which(infile[, "att"] == att_levels[att]), ]
  att_color <- att_colors[att]
  
  # add attributes one at a time
  plot <- plot +
    geom_segment(data = infile_filtered,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 10^line_width),
                 lineend = "butt", color = att_color, size = 5, alpha = alpha)
  
}

if (add_legend) {
  
  # open graphic device
  pdf(outfile_name)
  
  grid.newpage()
  grid.draw(arrangeGrob(plot))
  
  for (att in 1:length(att_levels)) {
    grid.rect(x = 0.7, y = 0.25 - (att / 20), width = 0.03, height = 0.01,
              gp = gpar(col = NA, fill = att_colors[att]))
    grid.text(x = 0.75, y = 0.25 - (att / 20), label = att_levels[att],
              gp = gpar(fontsize = 15))
  }
  
  # close graphic device for saving plot
  dev.off()
  
} else {
  
  # save plot
  ggsave(filename = outfile_name, plot = plot, device = "pdf")
  
}
