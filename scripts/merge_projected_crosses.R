#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script merges the results of SV projection of RILs into one file

      
      Usage: Rscript merge_projected_crosses.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 3) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript merge_projected_crosses.R [...]
       ")
}

# assign arguments to variables
hmp.file.before <- args[1]
proj.folder <- args[2]
cross.info <- args[3]


# hmp.file.before <- "data/usda_SNPs-SVs_325rils.not-in-SVs.hmp.txt"
# proj.folder <-"analysis/projection"
# cross.info <- "data/usda_biparental-crosses.txt"


#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### merge files ----


# load hmp before projection
hmp.before.projection <- fread(hmp.file.before, header = TRUE, data.table = FALSE)
# transform NA to NN in ril data
hmp.before.projection[, 12:NCOL(hmp.before.projection)] <- apply(hmp.before.projection[, 12:NCOL(hmp.before.projection)],
                                                                 MARGIN = 2, function(x) {
                                                                   x[which(is.na(x))] <- "NN"
                                                                   return(x)
                                                                 })

# open file with family name
df.crosses <- fread(cross.info, header = TRUE, data.table = FALSE)

# create main data frame to store results
hmp.proj.merged <- hmp.before.projection

# parse each cross
for (row in 1:NROW(df.crosses)) {
  
  # get cross name, and change * by x in the name
  cross <- df.crosses[row, "cross"]
  cross <- gsub(pattern = "*", replacement = "x", x = cross, fixed = TRUE)
  # get ril names for that cross
  ril.names <- unlist(strsplit(df.crosses[row, "RILs"], split = ","))
  
  # find file name for that cross after projection
  filename.after.proj <- list.files(path = proj.folder,
                                    pattern = paste0(cross, "_RILs.projected.hmp.txt"),
                                    full.names = TRUE)
  
  # make sure cross was use for projection
  if (length(filename.after.proj) > 0) {
    
    # open hmp after projection
    hmp.after.proj.cross <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
    # get only genotyped rils
    ril.names <- ril.names[which(ril.names %in% colnames(hmp.after.proj.cross)[12:NCOL(hmp.after.proj.cross)])]
    # add projected data into main data frame
    for (ril in ril.names) {
      hmp.proj.merged[, ril] <- hmp.after.proj.cross[, ril]
    }
    
  }
}


# write results
outfile <- gsub(".hmp.txt", ".projected.hmp.txt", hmp.file.before, fixed = TRUE)
fwrite(hmp.proj.merged, outfile, sep = "\t", quote = FALSE, na = NA)


# fix allele columns using TASSEL
commands.hapdip <- paste0("/home/hirschc1/della028/software/tassel-5-standalone/run_pipeline.pl",
                          " -importGuess ", outfile,
                          " -export ", outfile,
                          " -exportType HapmapDiploid")
system(commands.hapdip)
