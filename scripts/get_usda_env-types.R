library(data.table)
library(EnvRtype)

usage <- function() {
  cat("
description: get environment similarity from a file with geographic coordinates.

usage: Rscript get_usda_env-types.R [coords_filename] [outfile] [...]

positional arguments:
  coords_filename         file containing geographic coordinates (4 columns:
                          Location, Symbol, Latitude, Longitude)
  outfile                 name of output correlation matrix

optional argument:
  --help                  show this helpful message
  --country=VALUE         3-letter ISO code for country (default: USA)
  --start=VALUE           planting date (default: 2020-04-01)
  --end=VALUE             harvesting date (default: 2020-10-31)
  --heatmap               create a heatmap of environmental similarity
  
  
"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])
  
}



#### command line options ----

# set default
country <- "USA"
start <- "2020-04-01"
end <- "2020-10-31"
heatmap <- FALSE

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {
  
  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--country", "--start", "--end", "--heatmap")
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

# get positional arguments
coords_filename <- args[1]
outfile <- args[2]



#### get environmental similarity ----

# load file with location coordinates
sites_coord <- fread(coords_filename, header = TRUE, data.table = FALSE)

# get weather data
env_data <- get_weather(env.id = sites_coord$Symbol,
                        lat = sites_coord$Latitude,
                        lon = sites_coord$Longitude,
                        country = country,
                        start.day = rep(start, NROW(sites_coord)),
                        end.day = rep(end, NROW(sites_coord)))
# remove ALT columns
env_data <- env_data[, colnames(env_data) != "ALT"]

# # summarize data
# summaryWTH(env.data = env_data)

# process data
env_data = processWTH(env.data = env_data)

# T2M         Temperature at 2 Meters
# T2M_MAX     Maximum Temperature at 2 Meters
# T2M_MIN     Minimum Temperature at 2 Meters
# PRECTOT     Precipitation
# WS2M        Wind Speed at 2 Meters
# RH2M        Relative Humidity at 2 Meters
# T2MDEW      Dew/Frost Point at 2 Meters
# n           Actual duration of sunshine (hour)
# N           Daylight hours (hour)
# RTA         Extraterrestrial radiation (MJ/m^2/day)
# SRAD        Solar radiation (MJ/m^2/day)
# SPV         Slope of saturation vapour pressure curve (kPa.Celsius)
# VPD         Vapour pressure deficit (kPa)
# ETP         Potential Evapotranspiration (mm.day)
# PETP        Deficit by Precipitation (mm.day)
# GDD         Growing Degree Day (oC/day)
# FRUE        Effect of temperature on radiation use efficiency (from 0 to 1)
# T2M_RANGE   Daily Temperature Range (oC day)

# select environmental variables to use
env_vars <- c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOT", "WS2M", "RH2M", "T2MDEW",
              "n", "N", "RTA", "SRAD", "SPV", "VPD", "ETP", "PETP", "GDD","FRUE",
              "T2M_RANGE")

# # get environmental types using data-driven limits (instead of known cardinals)
# ET <- env_typing(env.data = env_data, env.id = 'env', var.id = env_vars)
# plot_panel(ET, title = 'Panel of Environmental Types')

# get environmental covariables
EC <- W_matrix(env.data = env_data, env.id = 'env', var.id = env_vars,
               statistic = "mean", by.interval = TRUE)
# plot_panel(EC, title = 'Panel of Environmental Covariables')

# get environmental similarity
K_E <- env_kernel(env.data = EC)[[2]]
K_E <- cov2cor(K_E)

if (heatmap) {
  library(pheatmap)
  heatmap_out <- rev(unlist(strsplit(outfile, ".", fixed = TRUE)))[-1]
  heatmap_out <- paste0(rev(heatmap_out)[1], ".pdf")
  pdf(heatmap_out)
  pheatmap(mat = K_E, breaks = seq(-1, 1, length.out = 100), angle_col = 0, main = 'Environmental Similarity')
  dev.off()
}

# write correlations
fwrite(K_E, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
