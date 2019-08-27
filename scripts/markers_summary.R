#### description ----

# this script generates a summary of genotypic data from parents or RILs



#### library ----

library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)



#### function ----

CountMarkers <- function(list_individuals, genotype_data) {
  
  df.results <- data.frame(line = as.character(), total = as.numeric(),
                           one_allele_missing = as.numeric(), both_alleles_missing = as.numeric(),
                           percent_missing = as.numeric(), available = as.numeric(),
                           homo = as.numeric(), het = as.numeric(),  percent_het = as.numeric(),
                           stringsAsFactors = FALSE)
  
  for (individual in list_individuals) {
    
    count.two.missing <- 0
    count.one.missing <- 0
    count.homo <- 0
    count.het <- 0
    
    for (marker in 1:NROW(genotype_data)) {
      
      alleles <- as.character(genotype_data[marker, individual])
      
      # if there are missing allele(s)...
      if (grepl("N", alleles)) {
        if (str_count(alleles, pattern = "N") == 2) {
          # ...check if both are missing...
          count.two.missing <- count.two.missing + 1
        } else {
          # ...or just one
          count.one.missing <- count.one.missing + 1
        }
      } else { 
        # if there both alleles are present, split...
        alleles <- unlist(strsplit(alleles, split = ""))
        if (alleles[1] == alleles[2]) {
          # ...and check if they are homozygous...
          count.homo <- count.homo + 1
        } else {
          # ...or het
          count.het <- count.het + 1
        }
      }
      
    }
    
    total.missing <- count.one.missing + count.two.missing
    total.markers <- count.homo + count.het + total.missing
    
    df.results <- rbind(df.results,
                        list(line = individual,
                             total = total.markers,
                             one_allele_missing = count.one.missing,
                             both_alleles_missing = count.two.missing,
                             percent_missing = round(total.missing/total.markers, digits = 3),
                             available = total.markers - total.missing,
                             homo = count.homo,
                             het = count.het,
                             percent_het = round(count.het / (total.markers - total.missing),
                                                 digits = 3)),
                        stringsAsFactors = FALSE)
    
  }
  
  return(df.results)
  
}




#### summarizing parents ----

# load data
geno.parents <- fread("data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt",
                      header = TRUE, data.table = FALSE)

# get names of parents
parent.lines <- colnames(geno.parents)[12:NCOL(geno.parents)]
# get genotypes of parents
parent.genotypes <- geno.parents[, which(colnames(geno.parents) %in% parent.lines)]
# get summary of markers
parent.summary <- CountMarkers(list_individuals = parent.lines,
                               genotype_data = parent.genotypes)

# write summary
fwrite(parent.summary, file = "analysis/qc/summary_markers_7parents.txt", sep = "\t", quote = FALSE)

# PHG35 has much more 10x hets than other parents (0.03 vs 0.003)...

# boxplots of missing and het data
ggplot(parent.summary, aes(x = line, y = percent_missing)) +
  geom_col() +
  labs(x = "Parent",
       y = "% missing markers") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 0.1))

ggsave(filename = "analysis/qc/summary_markers_7parents_missing.png", device = "png")

ggplot(parent.summary, aes(x = line, y = percent_het)) +
  geom_col() +
  labs(x = "Parent",
       y = "% het markers") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 0.1))

ggsave(filename = "analysis/qc/summary_markers_7parents_het.png", device = "png")





#### summarizing RILs per cross ----

# load data
geno.rils <- fread("data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt",
                   header = TRUE, data.table = FALSE)
# also need to load file telling which RILs belong to which crosses
crosses <- fread("data/usda_biparental-crosses.txt",
                 header = TRUE, data.table = FALSE)
# create empty df to store results
RILs.summary <- data.frame(stringsAsFactors = FALSE)

# get summary of markers for RILs in each cross
for (row in 1:NROW(crosses)) {
  
  # get cross name
  cross <- gsub(pattern = "*", replacement = "x", x = crosses[row, "cross"], fixed = TRUE)
  # get RILs in that cross
  RILs <- unlist(strsplit(crosses[row, "RILs"], split = ","))
  # get genotypes for those RILs
  geno.rils.by.cross <- geno.rils[, which(colnames(geno.rils) %in% RILs)]
  # make sure to use only RILs that are genotype
  RILs <- colnames(geno.rils.by.cross)
  
  # make sure that the cross have genotyped RILs
  if (length(RILs) > 0) {
    
    cat("Analyzing ",cross,"...\n")
    
    # get summary of markers
    markers.summary <- CountMarkers(list_individuals = RILs, genotype_data = geno.rils.by.cross)
    # add cross name into df
    cross.summary <- cbind(cross, markers.summary, stringsAsFactors = FALSE)
    # add results into final table
    RILs.summary <- rbind(RILs.summary, cross.summary, stringsAsFactors = FALSE)
    
    cat("Done!\n")
    
  }
  
}

# write summary
fwrite(RILs.summary, file = "analysis/qc/summary_markers_325rils.txt", sep = "\t", quote = FALSE)

# boxplots of missing and het data
ggplot(RILs.summary, aes(x = cross, y = percent_missing)) +
  geom_boxplot() +
  labs(x = "RIL population",
       y = "% missing markers") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "analysis/qc/summary_markers_325rils_missing-per-cross.png", device = "png")

ggplot(RILs.summary, aes(x = cross, y = percent_het)) +
  geom_boxplot() +
  labs(x = "RIL population",
       y = "% het markers") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 0.3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "analysis/qc/summary_markers_325rils_het-per-cross.png", device = "png")




#### summarizing RILs per parent ----

# load data
geno.rils <- fread("data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt",
                   header = TRUE, data.table = FALSE)

geno.parents <- fread("data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt",
                      header = TRUE, data.table = FALSE)

# also need to load file telling which RILs belong to which crosses
crosses <- fread("data/usda_biparental-crosses.txt",
                 header = TRUE, data.table = FALSE)

# get parents names
parents.names <- colnames(geno.parents)[12:NCOL(geno.parents)]

# create a list to store RILs belonging to the same parent
RILs.per.parent.list <- vector(mode = "list", length = 7)
names(RILs.per.parent.list) <- parents.names

for (row in 1:NROW(crosses)) {
  
  # get parent names for that cross
  parent1 <- unlist(strsplit(crosses[row, "cross"], split = "*", fixed = TRUE))[1]
  parent2 <- unlist(strsplit(crosses[row, "cross"], split = "*", fixed = TRUE))[2]
  
  # if both parents are in the list
  if (parent1 %in% parents.names & parent2 %in% parents.names) {
    
    # get RIL names
    RILs <- unlist(strsplit(crosses[row, "RILs"], split = ","))
    
    # append RIL names to each parent
    RILs.per.parent.list[[parent1]] <- append(RILs.per.parent.list[[parent1]], RILs)
    RILs.per.parent.list[[parent2]] <- append(RILs.per.parent.list[[parent2]], RILs)
  }
  
}

# create empty df to store results
RILs.summary.per.parent <- data.frame(stringsAsFactors = FALSE)

# now that i have all RIL names that share the same parent, I need to get the summary for each member of the list
for (parent in names(RILs.per.parent.list)) {
  
  
  
  # get all RIL names from the parent
  RILs.per.parent <- RILs.per.parent.list[[parent]]
  # get genotypes for each RIL
  geno.rils.by.parent <- geno.rils[, which(colnames(geno.rils) %in% RILs.per.parent)]
  # keep only genotyped rils in the list
  RILs.per.parent <- RILs.per.parent[RILs.per.parent %in% colnames(geno.rils.by.parent)]
  
  if (length(RILs.per.parent) > 0) {
    
    cat("Analyzing",parent,"...\n")
    
    # get summary of markers
    markers.summary <- CountMarkers(list_individuals = RILs.per.parent,
                                    genotype_data = geno.rils.by.parent)
    # add cross name into df
    parent.summary <- cbind(parent, markers.summary, stringsAsFactors = FALSE)
    # add results into final table
    RILs.summary.per.parent <- rbind(RILs.summary.per.parent, parent.summary, stringsAsFactors = FALSE)
    
    cat("Done!\n")
    
  }
  
}

# boxplots of missing and het data
ggplot(RILs.summary.per.parent, aes(x = parent, y = percent_missing)) +
  geom_boxplot() +
  labs(x = "Parent",
       y = "% missing markers") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 0.1))

ggsave(filename = "analysis/qc/summary_markers_325rils_missing-per-parent.png", device = "png")

ggplot(RILs.summary.per.parent, aes(x = parent, y = percent_het)) +
  geom_boxplot() +
  labs(x = "Parent",
       y = "% het markers") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 0.3))

ggsave(filename = "analysis/qc/summary_markers_325rils_het-per-parent.png", device = "png")




#### qc of results ----

# each line should have 20,139 markers
all(parent.summary$total == 20139)
all(RILs.summary$total == 20139)

# see how many missing data is due to one allele only
sum(parent.summary$one_allele_missing > 0)
sum(RILs.summary$one_allele_missing > 0)

# make sure the categories sum up to total
attach(parent.summary)
all(one_allele_missing + both_alleles_missing + homo + het == total)
detach()

attach(RILs.summary)
all(one_allele_missing + both_alleles_missing + homo + het == total)
detach()

# check which RILs have heterozygosity bigger than 3% (expected heterozigosity for F6s)
RILs.summary[which(RILs.summary[, "percent_het"] > 0.03), c("cross", "line")]

# get range/mean/median of missing data and het per RIL population
RILs.summary %>% 
  group_by(cross) %>% 
  summarize(min = min(percent_missing),
            max = max(percent_missing),
            mean = round(mean(percent_missing), digits = 3),
            median = median(percent_missing))


RILs.summary %>% 
  group_by(cross) %>% 
  summarize(min = min(percent_het),
            max = max(percent_het),
            mean = round(mean(percent_het), digits = 3),
            median = median(percent_het))

