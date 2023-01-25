# arguments for command line ----
args <- commandArgs(trailingOnly = TRUE)

# assign arguments to variables
filename <- args[1]

# change sv names
library(data.table)
map <- fread(filename, header = FALSE, data.table = FALSE)
new.names <- sapply(1:NROW(map), function(x) paste0("sv_", x))
map[, 2] <- new.names
fwrite(map, filename, sep = "\t", na = NA, row.names = FALSE, col.names = FALSE)
