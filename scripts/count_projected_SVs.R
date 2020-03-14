#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script summarizes SV projection and plots the amount projected
                   and projection accuracy per family.
      
      Usage: Rscript count_projected_SVs.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 3) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript count_projected_svs.R [...]
       ")
}

# assign arguments to variables
sv.filename <- args[1]
cross.info <- args[2]
proj.folder <- args[3]


# sv.filename <- "data/usda_SVs_7parents.sorted.hmp.txt"
# cross.info <- "data/usda_biparental-crosses.txt"
# proj.folder <- "analysis/projection"


#### libraries used ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



# open file with SV positions
sv.file <- fread(sv.filename, header = TRUE, data.table = FALSE)
sv.info <- sv.file[, c("chrom", "pos")]

# open file with family name
df.crosses <- fread(cross.info, header = TRUE, data.table = FALSE)

# create an empty df to store summary of SVs for that population
summary.projection <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(summary.projection) <- c("family", "total_SVs", "missing_SVs", "polymorphic_SVs", "avg_percent_projected", "avg_percent_projected_poly", "proj_accuracy")

# create empty vector to store the percentage of SVs projected for all individuals
proj.SVs.all.rils <- c()
proj.poly.SVs.all.rils <- c()


# generate summary for each cross
for (row in 1:NROW(df.crosses)) {
  
  # get cross name, and change * by x in the name
  cross <- df.crosses[row, "cross"]
  cross <- gsub(pattern = "*", replacement = "x", x = cross, fixed = TRUE)

  filename.after.proj <- list.files(path = proj.folder,
                                     pattern = paste0(cross, "_RILs.projected.hmp.txt"),
                                     full.names = TRUE)

  
  # make sure cross was use for projection
  if (length(filename.after.proj) > 0) {
    
    # open file before and after projection
    hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
    
    # filter file to have only SVs (transform below code in function...)
    hmp.after.filtered <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp.after)))
    colnames(hmp.after.filtered) <- colnames(hmp.after)
    
    for (chr in unique(sv.info[, "chrom"])) {
      
      sv.info.chr <- subset(sv.info, chrom == chr)
      positions.with.SVs <- sv.info.chr[, "pos"]
      
      hmp.after.chr <- subset(hmp.after, chrom == chr)
      hmp.after.chr.SVs <- hmp.after.chr[which(hmp.after.chr[, "pos"] %in% positions.with.SVs), ]
      hmp.after.filtered <- rbind(hmp.after.filtered, hmp.after.chr.SVs)
    }
    
    # count how many SVs were projected per RIL
    SVs.projected <- list()
    for (RIL in colnames(hmp.after.filtered)[12:NCOL(hmp.after.filtered)]) {
      SVs.projected[[RIL]] <- sum(hmp.after.filtered[, RIL] != "NN")
    }
    
    # count number total SVs
    number.total <- NROW(sv.file)
    # filter file to SV file to have only polymorphic SVs
    parent1 <- unlist(strsplit(cross, split = "x"))[1]
    parent2 <- unlist(strsplit(cross, split = "x"))[2]
    sv.file.poly <- sv.file[, c("chrom", "pos", parent1, parent2)]
    # filter missing data (if any)
    number.missing <- sum(sv.file.poly[, parent1] == "NN" & sv.file.poly[, parent2] == "NN")
    missing.SVs <- sv.file.poly[, parent1] == "NN" | sv.file.poly[, parent2] == "NN"
    sv.file.poly <- sv.file.poly[!missing.SVs, ]
    # keep only polymorphic SVs
    sv.file.poly <- sv.file.poly[which(sv.file.poly[, parent1] != sv.file.poly[, parent2]), ]
    number.poly.SVs <- NROW(sv.file.poly)
    
    # filter hapmap frame to have only projected polymorphic SVs
    hmp.after.filtered.poly <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp.after.filtered)))
    colnames(hmp.after.filtered.poly) <- colnames(hmp.after.filtered)
    for (chr in unique(sv.file.poly[, "chrom"])) {
      # subset data
      sv.file.poly.chr <- subset(sv.file.poly, chrom == chr)
      hmp.after.filtered.chr <- subset(hmp.after.filtered, chrom == chr)
      # select positions
      poly.SVs.positions <- sv.file.poly.chr[, "pos"]
      hmp.after.filtered.chr <- subset(hmp.after.filtered.chr, pos %in% poly.SVs.positions)
      # append to new df
      hmp.after.filtered.poly <- rbind(hmp.after.filtered.poly, hmp.after.filtered.chr)
    }
    # count how many polymorphic SVs were projected per RIL
    SVs.projected.poly <- list()
    for (RIL in colnames(hmp.after.filtered.poly)[12:NCOL(hmp.after.filtered.poly)]) {
      SVs.projected.poly[[RIL]] <- sum(hmp.after.filtered.poly[, RIL] != "NN")
    }
    
    # plot distribution
    proj.distribution <- as.numeric(unlist(SVs.projected))
    proj.distribution.plot <- ggplot(data.frame(dist = proj.distribution), aes(x = dist)) +
      geom_histogram() +
      labs(x = "SVs projected",
           y = "Number of RILs")
    ggsave(plot = proj.distribution.plot, filename = paste0(proj.folder, "/", cross, "_projection_distribution.png"), device = "png")

    # plot distribution polymorphic
    proj.distribution.poly <- as.numeric(unlist(SVs.projected.poly))
    proj.distribution.poly.plot <- ggplot(data.frame(dist = proj.distribution.poly), aes(x = dist)) +
      geom_histogram() +
      labs(x = "Polymorphic SVs projected",
           y = "Number of RILs")
    ggsave(plot = proj.distribution.poly.plot, filename = paste0(proj.folder, "/", cross, "_projection_distribution_poly-SVs.png"), device = "png")
    
    # add percent SVs projected of all individuals to vector
    proj.SVs.all.rils <- append(proj.SVs.all.rils, proj.distribution/(number.total - number.missing))
    proj.poly.SVs.all.rils <- append(proj.poly.SVs.all.rils, proj.distribution.poly/number.poly.SVs)
    
    # get accuracy
    # filename.accuracy <- list.files(path = "analysis/projection",
    #                                 pattern = paste0(cross, "_RILs.projected.hmp.Accuracy.txt"),
    #                                 full.names = TRUE)
    filename.accuracy <- list.files(path = proj.folder,
                                    pattern = paste0(cross, "_RILs.projected.hmp.Accuracy.txt"),
                                    full.names = TRUE)
    file.accuracy <- scan(file = filename.accuracy, what = "character", sep = "\n")
    proj.accuracy <- unlist(strsplit(file.accuracy[2], split = "\t"))
    proj.accuracy <- round(as.numeric(proj.accuracy[length(proj.accuracy)]), digits = 4)
    
    
    # get summary
    average.SVs <- round(mean(proj.distribution), digits = 1)
    percent.SVs <- round(average.SVs / (number.total - number.missing), digits = 4)
    average.SVs.poly <- round(mean(proj.distribution.poly), digits = 1)
    percent.SVs.poly <- round(average.SVs.poly / number.poly.SVs, digits = 4)
    
    summary.projection <- rbind(summary.projection,
                                list(family = cross,
                                     total_SVs = number.total,
                                     missing_SVs = number.missing,
                                     polymorphic_SVs = number.poly.SVs,
                                     avg_percent_projected = percent.SVs,
                                     avg_percent_projected_poly = percent.SVs.poly,
                                     proj_accuracy = proj.accuracy),
                                stringsAsFactors = FALSE)
    
    
    cat(cross, ": ", average.SVs, " (", round(percent.SVs * 100, digits = 2), "%) of SVs projected on average\n", sep = "")
    
  }
  
}

cat("Average SV projection: ", round(mean(summary.projection$avg_percent_projected) * 100, digits = 2), "% (",
    round(mean(summary.projection$proj_accuracy) * 100, digits = 2), "% accuracy)", sep = "")

outfile.sum <- paste0(proj.folder, "/projection_summary.txt")
fwrite(summary.projection, file = outfile.sum, sep = "\t", quote = FALSE, row.names = FALSE, na = NA)

proj.sum.plot <- ggplot(summary.projection, aes(x = family, y = avg_percent_projected)) +
  geom_col(fill = "#900721") + 
  labs(x = "Cross", y = "Projected SVs (%)", title = "SVs projected") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 1))+
  scale_x_discrete(labels = function(cross) gsub("x", "\nx\n", cross)) +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        text = element_text(size = 15))

ggsave(plot = proj.sum.plot, filename = paste0(proj.folder, "/projection_summary.png"), device = "png")


proj.sum.poly.plot <- ggplot(summary.projection, aes(x = family, y = avg_percent_projected_poly)) +
  geom_col(fill = "#900721") + 
  labs(x = "Cross", y = "Projected polymorphic SVs (%)", title = "Polymorphic SVs projected") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 1))+
  scale_x_discrete(labels = function(cross) gsub("x", "\nx\n", cross)) +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        text = element_text(size = 15))

ggsave(plot = proj.sum.poly.plot, filename = paste0(proj.folder, "/projection_poly-SVs_summary.png"), device = "png")


proj.accu.plot <- ggplot(summary.projection, aes(x = family, y = proj_accuracy)) +
  geom_col(fill = "#900721") + 
  labs(x = "Cross", y = "Projection accuracy (%)", title = "Projection accuracy") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 1))+
  scale_x_discrete(labels = function(cross) gsub("x", "\nx\n", cross)) +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        text = element_text(size = 15))

ggsave(plot = proj.accu.plot, filename = paste0(proj.folder, "/projection_accuracy_summary.png"), device = "png")
