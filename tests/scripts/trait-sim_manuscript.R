#### description ----

# Simulate different genetetic architectures for a trait in a single or multiple environments
# controlled by SNPs and/or SVs and run k-fold GS cross-validation on a genotypic dataset of a
# diversity panel.

# This script is actually a wrapper of functions/scripts written by Samuel Fernandes and Alex Lipka
# and adjusted by Rafael Della Coletta to generate the data for the simulation manuscript on SVs.




#### things to improve in following versions of the script ----

# 1. For the "simulate_trait function, need to find a way to sample equally from SNPs and from SVs
#    to be the source of variation, since there are much more SNPs in the dataset (instead of 
#    randomly selecting variants between SNPs and SVs).

# 2. For the "filter_SNPs_in_LD" function, I will need to how to select SNPs in LD by actually
#    calculating LD and not using an arbitrary number of bases up and downstream the SV.

# 3. Add functionality to write a log with the output produced in R console.

# 4. Run kfold validation for all 50 replicates: will need to uncomment a line in the function,
#    and find a way to optimize performance.

# 5. Use the latest GAPIT version available on their website (instead of the one Alex gave me)



#### BUG DETECTED ----

# If I run the simulate_trait() function right away, it doesn't work. I have to run the
# create_simulated_data() function outside the simulate_trait() function first, and then
# i can run simulate_trait() without any problems. I think there is a problem in setting up
# a connection to write the phenotypic files.



#### load libraries and functions ----

# the working directory will be set automatically when loading the R project

library(data.table)
library(lqmm)
library(mvtnorm)
library(rrBLUP)
library(MASS)
library(multtest)
library(gplots)

source("scripts/simulation_auxiliar-functions.R")
source("scripts/GAPIT_Code_from_Internet_20120411_Allelic_Effect.R")





#### functions ----

numericalize_hapmap <- function(geno_data, SNP_effect = "Add", SNP_impute = "Major") {
  hm = GAPIT.HapMap(geno_data, SNP.effect = SNP_effect, SNP.impute = SNP_impute, heading=TRUE)
  # create a new data frame for the numeric format
  num_hmp <- data.frame(matrix(data = NA, nrow = NROW(hm$GI), ncol = length(hm$GT) + 5))
  colnames(num_hmp)[1:5] <- c("Snp", "allele", "chr", "pos", "cm")
  colnames(num_hmp)[6:NCOL(num_hmp)] <- as.character(hm$GT)
  # add values to df
  num_hmp[,1] <- hm$GI[["SNP"]]
  num_hmp[,3] <- hm$GI[["Chromosome"]]
  num_hmp[,4] <- hm$GI[["Position"]]
  num_hmp[,6:NCOL(num_hmp)] <- t(hm$GD)
  return(num_hmp)
}



read_input_files <- function(snp_file = NULL, output_num = FALSE, sv_file = NULL, merge_SNP_SV = FALSE) {
  # if a snp file is provided
  if (!is.null(snp_file)) {
    # check if file is in standard or numeric hapmap format
    hmp <- fread(snp_file, header = TRUE, data.table = F)
    
    # check if file is in standard or numeric hapmap format by
    # looking at the elements in the the first genotype column
    if (is.numeric(hmp[,12]) == FALSE) {
      # if it's hapmap is not numeric, read in again but make sure to use "header = FALSE", for
      # compatibility with the rest of the script
      hmp <- fread(snp_file, header = FALSE, data.table = F)
      # if not numeric, transform to numeric
      num_hmp <- numericalize_hapmap(hmp, SNP_effect = "Add", SNP_impute = "Major")
      # write numeric hapmap if specified by the user
      if (output_num == TRUE) {
        outfile_name <- sub(pattern = ".txt", replacement = "", x = paste0(snp_file))
        write.table(num_hmp, file = paste0(outfile_name, "_NUM-test.txt"), quote = FALSE, sep = "\t",
                    row.names = FALSE)
        print("Numeric hapmap created")
        print(getwd())
      }
    }
    else {
      # file was in numeric format
      num_hmp <- hmp
    }
  }
  
  # if a sv file is provided
  if (!is.null(sv_file)) {
    # read in file with SVs only -- in the next version, this can be a simple list of SVs (if the
    # SVs are already in the middle of SNPs)
    geno_SVs <- fread(sv_file, header = TRUE, data.table = F)
    # merge snp and sv data, if user requested
    if (merge_SNP_SV == TRUE) {
      
      # # merge SVs into numeric hapmap SNP dataset by changing values to actual copy number
      # SVs <- geno_SVs[,1]
      # num_hmp_sv <- num_hmp
      # for (sv in SVs) {
      #   num_hmp_sv[which(num_hmp_sv[, 1] == sv),] <- geno_SVs[which(geno_SVs[, 1] == sv),]
      # }
      
      # merge by position
      geno_SVs <- geno_SVs[order(geno_SVs[,3], geno_SVs[,4]),]
      num_hmp_sv <- rbind(num_hmp, geno_SVs)
      num_hmp_sv <- num_hmp_sv[order(num_hmp_sv[,3], num_hmp_sv[,4]),]
    }
  }
  
  # return variables according to user input
  if (!is.null(snp_file) & is.null(sv_file) & merge_SNP_SV == FALSE) {
    return(num_hmp)
  }
  if (is.null(snp_file) & !is.null(sv_file) & merge_SNP_SV == FALSE) {
    return(geno_SVs)
  }
  if (!is.null(snp_file) & !is.null(sv_file) & merge_SNP_SV == FALSE) {
    return(list("num_hmp" = num_hmp, "geno_SVs" = geno_SVs))
  }
  if (!is.null(snp_file) & !is.null(sv_file) & merge_SNP_SV == TRUE) {
    return(list("num_hmp" = num_hmp, "geno_SVs" = geno_SVs, "num_hmp_sv" = num_hmp_sv))
  }
}



filter_SNPs_in_LD <- function(geno_data, sv_info, LD_type = "linked") {
  # since i still have to develop a way to calculate LD, i will use just the position of the SNPs
  # relative to a SV to "simulate" SNPs linked to SV (e.g. +/- 10kb) or in varying degrees of LD
  # (e.g. +/- 50kb)
  
  # 'geno_data' should be a data frame with SNP data (not a file name!)
  # 'sv_info' should be a tab-delimied data frame with 3 columns: sv_id, chrm, position
  # 'LD_type' should be "linked" or "varying"
  
  # create an empty data frame to store results
  geno_data_LD <- data.frame(matrix(nrow = 0, ncol = NCOL(geno_data)))
  colnames(geno_data_LD) <- colnames(geno_data)
  
  if (LD_type == "linked") print("Keeping only SNPs in LD with SVs...")
  if (LD_type == "varying") print("Keeping only SNPs in varying LD with SVs...")
  
  # for each chromosome in sv_info
  for (chrm in unique(sv_info[,2])) {
    # for each position in this subset
    pos_in_curr_chrm <- which(sv_info[,2] == chrm)
    for (pos in sv_info[pos_in_curr_chrm,3]) {
      # add 10kb up and downstream (if "linked") or 50kb (if "varying") to create a range of
      # possible values to select snps
      linked_SNPs_distance <- 10000
      if (LD_type == "linked") {
        # select the range of bp with SNPs linked to SV by adding 10kb up and downstream the SV position 
        range_LD <- seq(from = pos - linked_SNPs_distance, to = pos + linked_SNPs_distance, by = 1)
      }
      if (LD_type == "varying") {
        # select the range of bp with SNPs linked to SV by adding 50kb up and downstream the SV position
        # and by removing the bp within 10kb of the SV
        not_linked_SNPs_distance <- 50000
        range_linked_SNPs <- seq(from = pos - linked_SNPs_distance, to = pos + linked_SNPs_distance, by = 1)
        range_not_linked_SNPs <- seq(from = pos - not_linked_SNPs_distance, to = pos + not_linked_SNPs_distance, by = 1)
        # select only SNPs that are not linked to SV but have some degree of LD
        range_LD <- setdiff(range_not_linked_SNPs, range_linked_SNPs)
      }
      # go to geno_data, subset data by the appropriate chromosome
      geno_data_subset <- geno_data[which(geno_data[,3] == chrm),]
      # grab all rows for which the position column is in the possible range created
      geno_data_subset <- geno_data_subset[which(geno_data_subset[,4] %in% range_LD),]
      # append SNPs in LD to data frame
      geno_data_LD <- rbind(geno_data_LD, geno_data_subset)
    }
  }
  # remove any duplicates that may have been selected in case there were two SVs very close to each other
  geno_data_LD <- geno_data_LD[!duplicated(geno_data_LD),]
  # return a data frame only with SNPs in LD
  return(geno_data_LD)
}   



simulate_trait <- function(snp_data = NULL, sv_data = NULL, snp_with_sv_data = NULL,
                           list_of_SV_IDs = NULL, source_trait_variation = "SVs",
                           QTN_number = c(3,25,75), large_QTN_effect = 0.7,
                           small_QTN_effect = 0.3, heritability = c(0.2, 0.5, 0.9),
                           replicates = 50, output_folder_name = "analysis/one-trait_one-env",
                           seed_number = 2019) {
  
  # set home directory
  home_dir <- getwd()
  
  # try to find errors in the arguments provide
  if (is.null(snp_data) & is.null(sv_data) & is.null(snp_with_sv_data)) {
    stop("No genotypic data found! Please provide the name of a file with SNP and/or SV data.")
  }
  
  # if trait is explained by SVs -- select QTNs from genotypic data only with SVs
  if (source_trait_variation == "SVs") {
    
    # assert
    if (!is.null(snp_data) & is.null(sv_data) & is.null(snp_with_sv_data)) {
      stop("Wrong argument used! Please provide the name of a file with SV data on the 'sv_data' or 'snp_with_sv_data' arguments.")
    }
    # read in file depending on user input
    if (!is.null(sv_data)) {
      # get the correct genotypic data
      geno2trait_sim <- read_input_files(snp_file = NULL, sv_file = sv_data, merge_SNP_SV = FALSE)
    }
    if (is.null(sv_data) & !is.null(snp_with_sv_data)) {
      # get the correct genotypic data
      geno_data <- fread(snp_with_sv_data, header = TRUE, data.table = F)
      # make sure list of SVs is provided
      if (is.null(list_of_SV_IDs)) {
        stop("No list of SV IDs provided!")
      }
      # keep only SVs
      geno2trait_sim <- geno_data[which(geno_data[,1] %in% list_of_SV_IDs),]
    }
    
  }
  
  # if trait is explained by SNPs -- select QTNs from genotypic data only with SNPs
  if (source_trait_variation == "SNPs") {
    
    # assert
    if (is.null(snp_data) & !is.null(sv_data) & is.null(snp_with_sv_data)) {
      stop("Wrong argument used! Please provide the name of a file with SNP data on the 'snp_data' or 'snp_with_sv_data' arguments.")
    }
    # read in file depending on user input
    if (!is.null(snp_data)) {
      # get the correct genotypic data
      geno2trait_sim <- read_input_files(snp_file = snp_data, sv_file = NULL, merge_SNP_SV = FALSE)
    }
    if (is.null(snp_data) & !is.null(snp_with_sv_data)) {
      # get the correct genotypic data
      geno_data <- fread(snp_with_sv_data, header = TRUE, data.table = F)
      # make sure list of SVs is provided
      if (is.null(list_of_SV_IDs)) {
        stop("No list of SV IDs provided!")
      }
      # remove all SVs
      geno2trait_sim <- geno_data[which(!geno_data[,1] %in% list_of_SV_IDs),]
    }
    
  }
  
  # if trait is associated with both SNPs and SVs -- select QTNs randomly from both subsets
  # (SVs only and SNPs only; still need to find a way to sample equally from SNPs and from SVs,
  # since there are much more SNPs in the dataset below)
  if (source_trait_variation == "both") {
    
    # assert
    if (!is.null(snp_data) & is.null(sv_data) & is.null(snp_with_sv_data)) {
      stop("Missing SV file! Please provide the name of a file with SV data.")
    }
    if (is.null(snp_data) & !is.null(sv_data) & is.null(snp_with_sv_data)) {
      stop("Missing SNP file! Please provide the name of a file with SNP data.")
    }
    # read in file depending on user input
    if (!is.null(snp_data) & !is.null(sv_data) & is.null(snp_with_sv_data)) {
      print("Merging SNP and SV datasets...")
      geno_data <- read_input_files(snp_file = snp_data, sv_file = sv_data, merge_SNP_SV = TRUE)
      geno2trait_sim <- geno_data$num_hmp_sv
    }
    if (!is.null(snp_with_sv_data)) {
      geno2trait_sim <- fread(snp_with_sv_data, header = TRUE, data.table = F)
    }
    
    # will need to add more functionalities if i want to have equal numbers of SNPs and
    # SVs controlling the trait
    
  }
  
  # loop through different combinations of genetic architecture
  for (QTNs in QTN_number) {
    
    # create specific output folder names
    folder_name <- paste0(output_folder_name, "/", QTNs, "-QTNs_from_",
                          source_trait_variation, "/")
    # check if directory already exists, create a new one if it doesn't
    if (dir.exists(folder_name) == FALSE) {
      dir.create(folder_name, recursive = TRUE)
    }
    
    # simulate trait for single environment
    create.simulated.data(
      genotypes = geno2trait_sim,
      output.dir = folder_name,
      Additive.QTN.number = QTNs,
      additive.effect = small_QTN_effect,
      big.additive.QTN.effect = large_QTN_effect,
      rep = replicates,
      h2 = heritability,
      seed = seed_number
    )
    
    # go back to home directory
    setwd(home_dir)
  }
}



kfold_validation_on_sim_traits <- function(snp_data = NULL, sv_data = NULL, snp_with_sv_data = NULL,
                                           list_of_SV_IDs = NULL, sv_info = NULL,
                                           marker_data_type = 1, number_of_folds = 5,
                                           use_all_SNPs = TRUE, SNPs_to_sample = NULL,
                                           dir_with_sim_traits = "analysis/one-trait_one-env",
                                           seed_number = -673994, testing = FALSE) {
  
  # possible values for "marker_data_type" argument:
  # 1 = "all SNPs", 2 = "SNPs LD", 3 = "SNPs varying LD", 4 = "only SVs", 5 ="SNPs and SVs"
  possible_markers_for_GS <- c("all-SNPs", "SNPs-LD", "snps-varying-LD", "only-SVs", "SNPs-and-SVs")
  
  # try to find errors in the arguments provide
  if (is.null(snp_data) & is.null(sv_data) & is.null(snp_with_sv_data)) {
    stop("No genotypic data found! Please provide the name of a file with SNP and/or SV data.")
  }
  if (!marker_data_type %in% 1:5) {
    stop("Invalid 'marker_data_type' value! Possible values for this argument are '1' (all SNPs), '2' (SNPs in LD with SVs), '3' (SNPs in varying LD with SVs), '4' (only SVs), or '5' (SNPs and SVs).")
  }
  
  # read in the genotypic data according to the type of file provided by the user
  if (!is.null(snp_data) & is.null(sv_data) & is.null(snp_with_sv_data)) {
    # get the correct genotypic data
    geno_data <- read_input_files(snp_file = snp_data, sv_file = NULL, merge_SNP_SV = FALSE)
    print(paste("Total number of markers:", NROW(geno_data)))
    
    ### REMOVE MONOMORPHIC MARKERS
    print("Keeping only polymorphic markers...")
    monomorphic.markers <- c()
    for (row in 1:NROW(geno_data)) {
      # get all genotypes from a marker in the population
      marker.genotypes <- as.numeric(geno_data[row, 6:NCOL(geno_data)])
      # return TRUE if there is all individuals have the same genotype
      if (length(unique(marker.genotypes)) == 1) {
        # get the row numbers of polymorphic markers
        monomorphic.markers <- append(monomorphic.markers, row)
      }
    }
    print(paste("Monomorphic markers removed:", length(monomorphic.markers)))
    geno_data <- geno_data[-monomorphic.markers, ]
    
    
    # get list of SNP IDs
    SNPs_list <- geno_data[,1]
    
    # filter data if marker type is 2 (i.e., keep only SNPs in LD with SV)
    if (marker_data_type == 2) {
      if (is.null(sv_info)) {
        stop("No information about SV position found! Please provide a tab-delimied data frame with 3 columns: SV_id, chrm, position")
      }
      geno_data <- filter_SNPs_in_LD(geno_data, sv_info, LD_type = "linked")
      print(paste("Number of SNPs in LD kept:", NROW(geno_data)))
      # get list of SNP IDs
      SNPs_list <- geno_data[,1]
    }
    
    # filter data if marker type is 3 (i.e., keep only SNPs in varying LD with SV)
    if (marker_data_type == 3) {
      if (is.null(sv_info)) {
        stop("No information about SV position found! Please provide a tab-delimied data frame with 3 columns: SV_id, chrm, position")
      }
      geno_data <- filter_SNPs_in_LD(geno_data, sv_info, LD_type = "varying")
      print(paste("Number of SNPs in LD kept:", NROW(geno_data)))
      # get list of SNP IDs
      SNPs_list <- geno_data[,1]
    }
    
    # make sure that only the correct type of marker data can be used
    if (marker_data_type == 4 || marker_data_type == 5) {
      stop(paste("Cannot use", possible_markers_for_GS[marker_data_type], "markers with this dataset. Please, set 'marker_data_type' argument to '1' (all SNPs), '2' (SNPs in LD with SVs), or '3' (SNPs in varying LD with SVs)"))
    }
  }
  
  if (is.null(snp_data) & !is.null(sv_data) & is.null(snp_with_sv_data)) {
    # get the correct genotypic data
    geno_data <- read_input_files(snp_file = NULL, sv_file = sv_data, merge_SNP_SV = FALSE)
    print(paste("Total number of markers:", NROW(geno_data)))
    # make sure these options are set to default
    use_all_SNPs <- TRUE
    SNPs_to_sample <- NULL
    
    # make sure that only the correct type of marker data can be used
    if (marker_data_type != 4) {
      stop(paste("Cannot use", possible_markers_for_GS[marker_data_type], "markers with this dataset. Please, set 'marker_data_type' argument to '4' (only SVs)"))
    }
  }
  
  if (!is.null(snp_data) & !is.null(sv_data) & is.null(snp_with_sv_data)) {
    # get the correct genotypic data for marker types 1, 2 or 3
    if (marker_data_type == 1 ||  marker_data_type == 2 || marker_data_type == 3) {
      geno_data <- read_input_files(snp_file = snp_data, sv_file = NULL, merge_SNP_SV = FALSE)
      print(paste("Total number of markers:", NROW(geno_data)))
      
      ### REMOVE MONOMORPHIC MARKERS
      print("Keeping only polymorphic markers...")
      monomorphic.markers <- c()
      for (row in 1:NROW(geno_data)) {
        # get all genotypes from a marker in the population
        marker.genotypes <- as.numeric(geno_data[row, 6:NCOL(geno_data)])
        # return TRUE if there is all individuals have the same genotype
        if (length(unique(marker.genotypes)) == 1) {
          # get the row numbers of polymorphic markers
          monomorphic.markers <- append(monomorphic.markers, row)
        }
      }
      print(paste("Monomorphic markers removed:", length(monomorphic.markers)))
      geno_data <- geno_data[-monomorphic.markers, ]
      
      # get list of SNP IDs
      SNPs_list <- geno_data[,1]
      
      # filter data if marker type is 2 (i.e., keep only SNPs in LD with SV)
      if (marker_data_type == 2) {
        if (is.null(sv_info)) {
          stop("No information about SV position found! Please provide a tab-delimied data frame with 3 columns: SV_id, chrm, position")
        }
        geno_data <- filter_SNPs_in_LD(geno_data, sv_info, LD_type = "linked")
        print(paste("Number of SNPs in LD kept:", NROW(geno_data)))
        # get list of SNP IDs
        SNPs_list <- geno_data[,1]
      }
      
      # filter data if marker type is 3 (i.e., keep only SNPs in varying LD with SV)
      if (marker_data_type == 3) {
        if (is.null(sv_info)) {
          stop("No information about SV position found! Please provide a tab-delimied data frame with 3 columns: SV_id, chrm, position")
        }
        geno_data <- filter_SNPs_in_LD(geno_data, sv_info, LD_type = "varying")
        print(paste("Number of SNPs in LD kept:", NROW(geno_data)))
        # get list of SNP IDs
        SNPs_list <- geno_data[,1]
      }
    }
    
    # get the correct genotypic data for marker type 4
    if (marker_data_type == 4) {
      geno_data <- read_input_files(snp_file = NULL, sv_file = sv_data, merge_SNP_SV = FALSE)
      print(paste("Total number of markers:", NROW(geno_data)))
      # make sure these options are set to default
      use_all_SNPs <- TRUE
      SNPs_to_sample <- NULL
    }
    
    # get the correct genotypic data for marker type 5
    if (marker_data_type == 5) {
      print("Merging SNP and SV datasets...")
      geno_data <- read_input_files(snp_file = snp_data, sv_file = sv_data, merge_SNP_SV = TRUE)
      print(paste("Total number of markers:", NROW(geno_data)))
      # get sv info
      sv_info <- geno_data$geno_SVs
      sv_info <- sv_info[,c(1,3,4)]
      sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
      # get combined SNP + SV info
      geno_data <- geno_data$num_hmp_sv
      
      ### REMOVE MONOMORPHIC MARKERS
      print("Keeping only polymorphic markers...")
      monomorphic.markers <- c()
      for (row in 1:NROW(geno_data)) {
        # get all genotypes from a marker in the population
        marker.genotypes <- as.numeric(geno_data[row, 6:NCOL(geno_data)])
        # return TRUE if there is all individuals have the same genotype
        if (length(unique(marker.genotypes)) == 1) {
          # get the row numbers of polymorphic markers
          monomorphic.markers <- append(monomorphic.markers, row)
        }
      }
      print(paste("Monomorphic markers removed:", length(monomorphic.markers)))
      geno_data <- geno_data[-monomorphic.markers, ]
      
      # get list of SNP IDs
      SVs_IDs <- as.character(sv_info[,1])
      SNPs_list <- geno_data[which(geno_data[,1] %in% SVs_IDs), 1]
    }
  }
  
  if (!is.null(snp_with_sv_data)) {
    # get the correct genotypic data
    geno_data <- fread(snp_with_sv_data, header = TRUE, data.table = F)
    print(paste("Total number of markers:", NROW(geno_data)))
    
    ### REMOVE MONOMORPHIC MARKERS
    print("Keeping only polymorphic markers...")
    monomorphic.markers <- c()
    for (row in 1:NROW(geno_data)) {
      # get all genotypes from a marker in the population
      marker.genotypes <- as.numeric(geno_data[row, 6:NCOL(geno_data)])
      # return TRUE if there is all individuals have the same genotype
      if (length(unique(marker.genotypes)) == 1) {
        # get the row numbers of polymorphic markers
        monomorphic.markers <- append(monomorphic.markers, row)
      }
    }
    print(paste("Monomorphic markers removed:", length(monomorphic.markers)))
    geno_data <- geno_data[-monomorphic.markers, ]
    
    # get list of SNP IDs
    SNPs_list <- geno_data[which(!geno_data[,1] %in% list_of_SV_IDs),1]
    
    # use SV information from "list_of_SV_IDs" argument to filter SVs according to marker data type
    if (marker_data_type == 1 ||  marker_data_type == 2 || marker_data_type == 3) {
      # make sure list of SVs is provided
      if (is.null(list_of_SV_IDs)) {
        stop("No list of SV IDs provided!")
      }
      # remove all SVs
      geno_data <- geno_data[which(!geno_data[,1] %in% list_of_SV_IDs),]
      # get list of SNP IDs
      SNPs_list <- geno_data[,1]
      
      # filter data if marker type is 2 (i.e., keep only SNPs in LD with SV)
      if (marker_data_type == 2) {
        if (is.null(sv_info)) {
          stop("No information about SV position found! Please provide a tab-delimied data frame with 3 columns: SV_id, chrm, position")
        }
        geno_data <- filter_SNPs_in_LD(geno_data, sv_info, LD_type = "linked")
        print(paste("Number of SNPs in LD kept:", NROW(geno_data)))
        # get list of SNP IDs
        SNPs_list <- geno_data[,1]
      }
      
      # filter data if marker type is 3 (i.e., keep only SNPs in varying LD with SV)
      if (marker_data_type == 3) {
        if (is.null(sv_info)) {
          stop("No information about SV position found! Please provide a tab-delimied data frame with 3 columns: SV_id, chrm, position")
        }
        geno_data <- filter_SNPs_in_LD(geno_data, sv_info, LD_type = "varying")
        print(paste("Number of SNPs in LD kept:", NROW(geno_data)))
        # get list of SNP IDs
        SNPs_list <- geno_data[,1]
      }
    }
    
    if (marker_data_type == 4) {
      # make sure list of SVs is provided
      if (is.null(list_of_SV_IDs)) {
        stop("No list of SV IDs provided!")
      }
      # keep only SVs
      geno_data <- geno_data[which(geno_data[,1] %in% list_of_SV_IDs),]
      # make sure these options are set to default
      use_all_SNPs <- TRUE
      SNPs_to_sample <- NULL
    }
    # if "marker_data_type = 5" don't need to do anything
  }
  
  # assign all folders with simulated traits to a variable
  trait_dirs <- list.dirs(path = dir_with_sim_traits, recursive = FALSE)
  
  
  
  ### USE THE NUMBER OF MARKERS PROVIDED BY USER -- TEST EACH FILE COMBINATION!
  if (use_all_SNPs == FALSE) {
    # assert if SNPs_to_sample is NULL
    if (is.null(SNPs_to_sample)) {
      stop("Please provide the number of SNPs to be sampled for the prediction model.")
    }
    print(paste("Sampling", SNPs_to_sample, "markers..."))
    # sample only SNPs without replacement
    set.seed(seed_number)
    SNPs.sampled <- sample(SNPs_list, SNPs_to_sample, replace = FALSE)
    geno_data <- geno_data[which(geno_data[,1] %in% SNPs.sampled), ]
  }

  # loop though each folder
  for (curr_dir in trait_dirs) {
    # once inside a folder, open each simulated trait file at a time
    # select only files that start (^) with "Simulated.Data"
    trait_files <- list.files(path = curr_dir, pattern = "^Simulated.Data")
    
    # for each simulated trait file, perform k-fold cross validation for each of the 50 replicates
    # (or should i take average from 50 reps and do only one cross validation?)
    for (curr_file in trait_files) {
      # debug
      print("------------------------------------")
      print(paste0("file being parsed: ", curr_dir, curr_file))
      print(paste0("dimensions geno_data: ", NROW(geno_data), " x ", NCOL(geno_data)))
      print(paste0("marker type: ", marker_data_type))
      print("------------------------------------")
      # read in file with simulated trait replicated 50 times
      sim_trait_file <- fread(paste0(curr_dir, "/", curr_file), data.table = F)
      
      # select columns with trait values for each replicate
      columns_of_replicates <- 2:NCOL(sim_trait_file)
      
      if (testing == TRUE) {
        # validate only the first replicate (column 2)
        columns_of_replicates <- 2:4
      }
      
      # perform k-fold validation in each replicate
      for (i in columns_of_replicates) {
        # select one replicate at a time
        myY <- sim_trait_file[,c(1,i)]
        # retrieve name of trait being parsed (heritability and replicate number)
        trait_name <- names(sim_trait_file[i])
        
        # create new directory to store validation results
        dir_heritability <- sub(pattern = "_Rep_[0-9]+", replacement = "", x = trait_name)
        if (use_all_SNPs == TRUE) {
          new_dir <- paste0(curr_dir, "/gs_without_subsetting_markers/k-fold_validation_",
                            possible_markers_for_GS[marker_data_type], "/", dir_heritability, "/")
        }
        if (use_all_SNPs == FALSE) {
          new_dir <- paste0(curr_dir, "/gs_with_", SNPs_to_sample, "_SNPs/k-fold_validation_",
                            possible_markers_for_GS[marker_data_type], "/", dir_heritability, "/")
        }
        # check if directory already exists, create a new one if it doesn't
        if (dir.exists(new_dir) == FALSE) {
          dir.create(new_dir, recursive = TRUE)
        }
        
        # since the genotypic data will already be in the numeric format, i need to format
        # the genotypic like GAPIT.HapMap() would before running the validation. I will also
        # have to using an edited version of the kfold script, where i removed the numericalization
        # process
        
        hm <- list()
        
        GT <- matrix(colnames(geno_data)[-(1:5)], ncol = 1)  
        GI <- geno_data[,c(1,3,4)]
        colnames(GT)="taxa"
        colnames(GI)=c("SNP","Chromosome","Position")
        GD <- geno_data[,-(1:5)]
        
        hm[["GT"]] <- GT
        hm[["GD"]] <- apply(GD, 1, function(x) as.numeric(x))
        hm[["GI"]] <- GI
        
        # run k-fold cross validation -- use my slightly modified version of Alex's script
        rrblup.kfoldfoldCV.numHapMap.input(Y = myY,
                                           Geno = hm,
                                           traitname = trait_name,
                                           path.for.results = new_dir,
                                           number.of.folds = number_of_folds,
                                           seed.number = seed_number)
      }
    }
  }
}



#### simulate trait for single environment ----

# run this first so my simulate_trait() function can run (do this until i debug)
# this will even throw an error saying it couldn't find "num_hmp". It's ok. I just
# need to run this so the function following this code works. There must be some
# problem in the sink() behavior.

{
  # write a text file just saying to ignore this folder
  dir.create("analysis/test_toy/test_function/", recursive = TRUE)
  cat("Please ignore this entire folder. If you are reading this, you can delete it.",
      file = "analysis/test_toy/test_function/README")
  
  # test run
  geno_SVs <- fread("data/Structural_Variation_demo.txt", header = TRUE, data.table = F)
  SVs <- geno_SVs[,1]
  QTN_number <- c(3,15,30)
  origin_trait_variance <- c("SVs", "SNPs", "both")
  for (qtn in QTN_number) {
    for (origin in origin_trait_variance) {
      # create genetic architecture
      Additive.QTN.number <- qtn
      big.additive.QTN.effect <- 0.7
      additive.effect <- 0.3
      heritabilities.vector <- c(0.2, 0.5, 0.9)
      replicates <- 50
      origin_trait_variation <- origin   # choose between "SNPs", "SVs", "both"
      
      if (origin_trait_variation == "SVs") {
        # if trait is explained by SVs -- select QTNs from genotypic data only with SVs
        geno2trait_sim <- geno_SVs
      }
      
      if (origin_trait_variation == "SNPs") {
        # if trait is explained by SNPs -- select QTNs from genotypic data only with SNPs
        geno2trait_sim <- num_hmp[which(!num_hmp[,1] %in% SVs),]
        # colnames(geno2trait_sim) <- geno2trait_sim[1,]
        # geno2trait_sim <- geno2trait_sim[-1,]
      }
      
      if (origin_trait_variation == "both") {
        # if trait is associated with both SNPs and SVs -- select QTNs randomly from both subsets
        # (SVs only and SNPs only; still need to find a way to sample equally from SNPs and from SVs,
        # since there are much more SNPs in the dataset below)
        geno2trait_sim <- num_hmp_sv
      }
      
      # check if directory already exists, create a new one if it doesn't
      folder_name <- paste0("analysis/test_toy/test_function/", Additive.QTN.number, "-QTNs_from_",
                            origin_trait_variation, "/")
      if (dir.exists(folder_name) == FALSE) {
        dir.create(folder_name, recursive = TRUE)
      }
      
      # simulate trait for single environment
      create.simulated.data(
        genotypes = geno2trait_sim,
        output.dir = folder_name,
        Additive.QTN.number = Additive.QTN.number,
        additive.effect = additive.effect,
        big.additive.QTN.effect = big.additive.QTN.effect,
        rep = replicates,
        h2 = heritabilities.vector,
        seed = 2019
      )
    }
  }
}


# now that the code above was ran, the simulation using my function simulate_trait() will work.

# get the SV IDs list
geno_SVs <- fread("data/Structural_Variation_demo.txt", header = TRUE, data.table = F)
SVs <- geno_SVs[,1]

# simulate data for toy dataset
for (source_var in c("SNPs", "SVs", "both")) {
  simulate_trait(snp_data = NULL,
                 sv_data = NULL,
                 snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                 list_of_SV_IDs = SVs,
                 source_trait_variation = source_var,
                 QTN_number = c(3,25,75),
                 large_QTN_effect = 0.7,
                 small_QTN_effect = 0.3,
                 heritability = c(0.2, 0.5, 0.9),
                 replicates = 50,
                 output_folder_name = "analysis/test_toy",
                 seed_number = 2019)
}

# simulate data for USDA dataset (only SNPs, since I don't have SV info yet)
simulate_trait(snp_data = "data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt",
               sv_data = NULL,
               snp_with_sv_data = NULL,
               list_of_SV_IDs = NULL,
               source_trait_variation = "SNPs",
               QTN_number = c(3,25,75),
               large_QTN_effect = 0.7,
               small_QTN_effect = 0.3,
               heritability = c(0.2, 0.5, 0.9),
               replicates = 50,
               output_folder_name = "analysis/test_usda-parents",
               seed_number = 2019)

simulate_trait(snp_data = "data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt",
               sv_data = NULL,
               snp_with_sv_data = NULL,
               list_of_SV_IDs = NULL,
               source_trait_variation = "SNPs",
               QTN_number = c(3,25,75),
               large_QTN_effect = 0.7,
               small_QTN_effect = 0.3,
               heritability = c(0.2, 0.5, 0.9),
               replicates = 50,
               output_folder_name = "analysis/test_usda-rils",
               seed_number = 2019)


#### simulate trait for multiple environment ----

# later...



#### k-fold cross validation (toy dataset) ----

# get SV information
geno_SVs <- fread("data/Structural_Variation_demo.txt", header = TRUE, data.table = F)
SVs <- geno_SVs[,1]
sv_info <- geno_SVs[,c(1,3,4)]
sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]

# set a seed number -- this will ensure that the same folds are being used.
# this.seed.number <- sample(-1000000:1000000,1)
this.seed.number <- -673994


# genotypic data should be one of the following: (1) all SNPs, (2) only SNPs in LD with SV, (3) SNPs
# in varying LD with SV, (4) the actual SV, or (5) all SNPs and SVs. So far i have only data with
# all SNPs. Will have to select SNPs in LD later.


#-------------- using all markers --------------#

# (1) all SNPs as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 1,
                               number_of_folds = 5,
                               use_all_SNPs = TRUE,
                               SNPs_to_sample = NULL,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 26min

# (2) only SNPs in LD with SV as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 2,
                               number_of_folds = 5,
                               use_all_SNPs = TRUE,
                               SNPs_to_sample = NULL,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 5min

# (3) SNPs in varying LD with SV as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 3,
                               number_of_folds = 5,
                               use_all_SNPs = TRUE,
                               SNPs_to_sample = NULL,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 5min

# (4) only actual SVs as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 4,
                               number_of_folds = 5,
                               use_all_SNPs = TRUE,
                               SNPs_to_sample = NULL,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 2min

# (5) all SNPs and SVs as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 5,
                               number_of_folds = 5,
                               use_all_SNPs = TRUE,
                               SNPs_to_sample = NULL,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 35min


#-------------- using 1000 markers --------------#

############## (1) all SNPs as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 1,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 1000,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 2min

# (2) only SNPs in LD with SV as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 2,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 1000,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 3min

# (3) SNPs in varying LD with SV as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 3,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 1000,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 3min

# (5) all SNPs and SVs as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 5,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 1000,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 2min


#-------------- using 50 markers --------------#

# (1) all SNPs as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 1,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 50,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 2min

# (2) only SNPs in LD with SV as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 2,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 50,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 3min

# (3) SNPs in varying LD with SV as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 3,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 50,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 3min

# (5) all SNPs and SVs as markers
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = NULL,
                               sv_data = NULL,
                               snp_with_sv_data = "data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
                               list_of_SV_IDs = SVs,
                               sv_info = sv_info,
                               dir_with_sim_traits = "analysis/test_toy/",
                               marker_data_type = 5,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 50,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 2min



#### k-fold cross validation (USDA dataset) ----

#-------------- using all markers --------------#

# (1) all SNPs as markers - parental data
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = "data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt",
                               sv_data = NULL,
                               snp_with_sv_data = NULL,
                               list_of_SV_IDs = NULL,
                               sv_info = NULL,
                               dir_with_sim_traits = "analysis/test_usda-parents",
                               marker_data_type = 1,
                               number_of_folds = 5,
                               use_all_SNPs = TRUE,
                               SNPs_to_sample = NULL,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 18s

# (1) all SNPs as markers - ril data
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = "data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt",
                               sv_data = NULL,
                               snp_with_sv_data = NULL,
                               list_of_SV_IDs = NULL,
                               sv_info = NULL,
                               dir_with_sim_traits = "analysis/test_usda-rils",
                               marker_data_type = 1,
                               number_of_folds = 5,
                               use_all_SNPs = TRUE,
                               SNPs_to_sample = NULL,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 4min


#-------------- using 1000 markers --------------#

# (1) all SNPs as markers - parental data
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = "data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt",
                               sv_data = NULL,
                               snp_with_sv_data = NULL,
                               list_of_SV_IDs = NULL,
                               sv_info = NULL,
                               dir_with_sim_traits = "analysis/test_usda-parents",
                               marker_data_type = 1,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 1000,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 10s

# (1) all SNPs as markers - ril data
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = "data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt",
                               sv_data = NULL,
                               snp_with_sv_data = NULL,
                               list_of_SV_IDs = NULL,
                               sv_info = NULL,
                               dir_with_sim_traits = "analysis/test_usda-rils",
                               marker_data_type = 1,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 1000,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 1min


#-------------- using 50 markers --------------#

# (1) all SNPs as markers - parental data
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = "data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt",
                               sv_data = NULL,
                               snp_with_sv_data = NULL,
                               list_of_SV_IDs = NULL,
                               sv_info = NULL,
                               dir_with_sim_traits = "analysis/test_usda-parents",
                               marker_data_type = 1,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 50,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 8s

# (1) all SNPs as markers - ril data
start_time <- proc.time()
kfold_validation_on_sim_traits(snp_data = "data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt",
                               sv_data = NULL,
                               snp_with_sv_data = NULL,
                               list_of_SV_IDs = NULL,
                               sv_info = NULL,
                               dir_with_sim_traits = "analysis/test_usda-rils",
                               marker_data_type = 1,
                               number_of_folds = 5,
                               use_all_SNPs = FALSE,
                               SNPs_to_sample = 50,
                               seed_number = -673994,
                               testing = TRUE)
timetaken(start_time)  # 1min


#### test functions  ----

############ read in datasets

# # manually
# 
# hmp <- fread("SNP55K_maize282_AGPv2_20100513_1.hmp.txt", header = FALSE, data.table = F)
# num_hmp <- fread("SNP55K_maize282_AGPv2_20100513_NUM.txt", header = TRUE, data.table = F)
# geno_SVs <- fread("Structural_Variation_demo.txt", header = TRUE, data.table = F)
# 
# SVs <- geno_SVs[,1]
# num_hmp_sv <- num_hmp
# for (sv in SVs) {
#   num_hmp_sv[which(num_hmp_sv[, 1] == sv),] <- geno_SVs[which(geno_SVs[, 1] == sv),]
# }
# 
# # or using function
# 
# num_hmp <- read_input_files(snp_file = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",
#                             output_num = FALSE, sv_file = NULL, merge_SNP_SV = FALSE)
# 
# geno_SVs <- read_input_files(snp_file = NULL, output_num = FALSE,
#                              sv_file = "Structural_Variation_demo.txt", merge_SNP_SV = FALSE)
# 
# geno_list <- read_input_files(snp_file = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",
#                               output_num = FALSE, sv_file = "Structural_Variation_demo.txt",
#                               merge_SNP_SV = TRUE)
# 
# num_hmp_sv <- geno_list$num_hmp_sv
# 
# 
# # get SV info data
# sv_info <- geno_SVs[,c(1,3,4)]
# # sort by chrm and position
# sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
# # since SVs were simulated from SNPs, I will add 1 to the each position so they don't overlap with SNPs
# sv_info[,3] <- sv_info[,3] + 1



############ trait simulation function

# # no dataset
# simulate_trait(snp_data = NULL,
#                sv_data = NULL,
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SNPs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# # source variation: SNPs
# simulate_trait(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",  #### this should work
#                sv_data = NULL,
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SNPs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = NULL,
#                sv_data = "Structural_Variation_demo.txt",  #### this should NOT work
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SNPs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt",  #### this should work
#                sv_data = NULL,
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SNPs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",  #### this should work
#                sv_data = "Structural_Variation_demo.txt",
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SNPs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = NULL,
#                sv_data = NULL,
#                snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",  #### this should NOT work
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SNPs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = NULL,
#                sv_data = NULL,
#                snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",  #### this should work
#                list_of_SV_IDs = SVs,
#                source_trait_variation = "SNPs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# 
# # source variation: SVs
# simulate_trait(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",  #### this should NOT work
#                sv_data = NULL,
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SVs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = NULL,
#                sv_data = "Structural_Variation_demo.txt",  #### this should work
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SVs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt",  #### this should NOT work
#                sv_data = NULL,
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SVs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",
#                sv_data = "Structural_Variation_demo.txt",  #### this should work
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SVs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = NULL,
#                sv_data = NULL,
#                snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",  #### this should NOT work
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "SVs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = NULL,
#                sv_data = NULL,
#                snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",  #### this should work
#                list_of_SV_IDs = SVs,
#                source_trait_variation = "SVs",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# 
# # source variation: both
# simulate_trait(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",  #### this should NOT work
#                sv_data = NULL,
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "both",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = NULL,
#                sv_data = "Structural_Variation_demo.txt",  #### this should NOT work
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "both",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt",  #### this should NOT work
#                sv_data = NULL,
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "both",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",
#                sv_data = "Structural_Variation_demo.txt",  #### this should work
#                snp_with_sv_data = NULL,
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "both",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = NULL,
#                sv_data = NULL,
#                snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",  #### this should work
#                list_of_SV_IDs = NULL,
#                source_trait_variation = "both",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)
# 
# simulate_trait(snp_data = NULL,
#                sv_data = NULL,
#                snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",  #### this should work
#                list_of_SV_IDs = SVs,
#                source_trait_variation = "both",
#                QTN_number = c(3,15,30),
#                large_QTN_effect = 0.7,
#                small_QTN_effect = 0.3,
#                heritability = c(0.2, 0.5, 0.9),
#                replicates = 50,
#                output_folder_name = "test_function",
#                seed_number = 2019)



############ kfold validation function

# # no dataset
# kfold_validation_on_sim_traits(snp_data = NULL,
#                                sv_data = NULL,
#                                snp_with_sv_data = NULL,
#                                list_of_SV_IDs = NULL,
#                                sv_info = NULL,
#                                marker_data_type = 1,
#                                number_of_folds = 5,
#                                seed_number = -673994,
#                                testing = TRUE)
# # RESULT: worked as expected
# 
# 
# # wrong marker data type (0)
# kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",
#                                sv_data = NULL,
#                                snp_with_sv_data = NULL,
#                                list_of_SV_IDs = NULL,
#                                sv_info = NULL,
#                                marker_data_type = 0,
#                                number_of_folds = 5,
#                                seed_number = -673994,
#                                testing = TRUE)
# # RESULT: worked as expected
# 
# 
# # geno data: SNPs -- format: hmp -- markers: all SNPs
# kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",
#                                sv_data = NULL,
#                                snp_with_sv_data = NULL,
#                                list_of_SV_IDs = NULL,
#                                sv_info = NULL,
#                                marker_data_type = 5,  # test all 5 possible marker types manually
#                                number_of_folds = 5,
#                                seed_number = -673994,
#                                testing = TRUE)
# # RESULT: worked as expected -- only marker type 1 works
# 
# 
# # geno data: SNPs -- format: hmp -- markers: all SNPs, SNPs in LD, SNPs in varying LD
# for (i in 1:5) {
#   kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",
#                                  sv_data = NULL,
#                                  snp_with_sv_data = NULL,
#                                  list_of_SV_IDs = NULL,
#                                  sv_info = sv_info,
#                                  marker_data_type = i,
#                                  number_of_folds = 5,
#                                  seed_number = -673994,
#                                  testing = TRUE)
# }
# # RESULT: worked as expected -- only marker types 1, 2 and 3 worked.
# 
# 
# # geno data: SNPs -- format: numeric -- markers: all SNPs
# kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt",
#                                sv_data = NULL,
#                                snp_with_sv_data = NULL,
#                                list_of_SV_IDs = NULL,
#                                sv_info = NULL,
#                                marker_data_type = 5,  # test all 5 possible marker types manually
#                                number_of_folds = 5,
#                                seed_number = -673994,
#                                testing = TRUE)
# # RESULT: worked as expected -- only marker type 1 works
# 
# 
# # geno data: SNPs -- format: numeric -- markers: all SNPs, SNPs in LD, SNPs in varying LD
# for (i in 1:5) {
#   kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt",
#                                  sv_data = NULL,
#                                  snp_with_sv_data = NULL,
#                                  list_of_SV_IDs = NULL,
#                                  sv_info = sv_info,
#                                  marker_data_type = i,
#                                  number_of_folds = 5,
#                                  seed_number = -673994,
#                                  testing = TRUE)
# }
# # RESULT: worked as expected -- only marker types 1, 2 and 3 worked.
# 
# 
# # geno data: SVs -- format: numeric -- markers: only SVs
# kfold_validation_on_sim_traits(snp_data = NULL,
#                                sv_data = "Structural_Variation_demo.txt",  # no way to be on standard hmp format
#                                snp_with_sv_data = NULL,
#                                list_of_SV_IDs = NULL,
#                                sv_info = NULL,
#                                marker_data_type = 5, # test all 5 possible marker types manually
#                                number_of_folds = 5,
#                                seed_number = -673994,
#                                testing = TRUE)
# # RESULT: worked as expected -- only marker type 4 works
# 
# 
# # geno data: SNPs and SVs (separate files) -- format: hmp -- markers: all SNPs, only SVs, SNPs and SVs
# kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",
#                                sv_data = "Structural_Variation_demo.txt",
#                                snp_with_sv_data = NULL,
#                                list_of_SV_IDs = NULL,
#                                sv_info = NULL,
#                                marker_data_type = 5,  # test all 5 possible marker types manually 
#                                number_of_folds = 5,
#                                seed_number = -673994,
#                                testing = TRUE)
# # RESULT: worked as expected -- only marker types 1, 4 and 5 worked.
# 
# 
# # geno data: SNPs and SVs (separate files) -- format: hmp -- markers: all 5 types
# for (i in 1:5) {
#   kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_1.hmp.txt",
#                                  sv_data = "Structural_Variation_demo.txt",
#                                  snp_with_sv_data = NULL,
#                                  list_of_SV_IDs = NULL,
#                                  sv_info = sv_info,
#                                  marker_data_type = i,  # test all 5 possible marker types manually 
#                                  number_of_folds = 5,
#                                  seed_number = -673994,
#                                  testing = TRUE)
# }
# # RESULT: all marker types worked
# 
# 
# # geno data: SNPs and SVs (separate files) -- format: hmp -- markers: all SNPs, only SVs, SNPs and SVs
# kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt",
#                                sv_data = "Structural_Variation_demo.txt",
#                                snp_with_sv_data = NULL,
#                                list_of_SV_IDs = NULL,
#                                sv_info = NULL,
#                                marker_data_type = 5,  # test all 5 possible marker types manually 
#                                number_of_folds = 5,
#                                seed_number = -673994,
#                                testing = TRUE)
# # RESULT: worked as expected -- only marker types 1, 4 and 5 worked.
# 
# 
# # geno data: SNPs and SVs (separate files) -- format: numeric -- markers: all 5 possibilities
# for (i in 1:5) {
#   kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt",
#                                  sv_data = "Structural_Variation_demo.txt",
#                                  snp_with_sv_data = NULL,
#                                  list_of_SV_IDs = NULL,
#                                  sv_info = sv_info,
#                                  marker_data_type = i,
#                                  number_of_folds = 5,
#                                  seed_number = -673994,
#                                  testing = TRUE)
# }
# # RESULT: worked for all markers
# 
# 
# # geno data: SNPs and SVs (same file) -- format: num_hmp -- markers: SNPs and SVs
# kfold_validation_on_sim_traits(snp_data = NULL,
#                                sv_data = NULL,
#                                snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
#                                list_of_SV_IDs = NULL,
#                                sv_info = NULL,
#                                marker_data_type = 5, # test all 5 possible marker types manually 
#                                number_of_folds = 5,
#                                seed_number = -673994,
#                                testing = TRUE)
# # RESULT: worked as expected -- only marker type 5 worked
# 
# 
# 
# # geno data: SNPs and SVs (same file) -- format: num_hmp -- markers: all SNPs, only SVs, SNPs and SVs
# for (i in 1:5) {
#   print("################")
#   print(paste("Marker type:", i))
#   print("################")
#   try(
#     kfold_validation_on_sim_traits(snp_data = NULL,
#                                    sv_data = NULL,
#                                    snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
#                                    list_of_SV_IDs = SVs,
#                                    sv_info = NULL,
#                                    marker_data_type = i,
#                                    number_of_folds = 5,
#                                    seed_number = -673994,
#                                    testing = TRUE)
#   )
# }
# # RESULT: worked as expected -- only marker types 1, 4 and 5 worked
# 
# 
# # geno data: SNPs and SVs (same file) -- format: num_hmp -- markers: all 5 possibilities
# for (i in 1:5) {
#   kfold_validation_on_sim_traits(snp_data = NULL,
#                                  sv_data = NULL,
#                                  snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
#                                  list_of_SV_IDs = SVs,
#                                  sv_info = sv_info,
#                                  marker_data_type = i,
#                                  number_of_folds = 5,
#                                  seed_number = -673994,
#                                  testing = TRUE) 
# }
# # RESULT: worked for all marker types
# 
# 
# # geno data: SNPs and SVs (same file) -- format: num_hmp -- markers: SNPs and SVs
# for (i in 1:5) {
#   print("################")
#   print(paste("Marker type:", i))
#   print("################")
#   try(
#   kfold_validation_on_sim_traits(snp_data = NULL,
#                                  sv_data = NULL,
#                                  snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
#                                  list_of_SV_IDs = NULL,
#                                  sv_info = sv_info,
#                                  marker_data_type = i,
#                                  number_of_folds = 5,
#                                  seed_number = -673994,
#                                  testing = TRUE)
#   )
# }
# # RESULT: worked as expected -- only marker 5 worked
# 
# 
# # should work for all markers
# # geno data: SNPs and SVs (same file) -- format: num_hmp -- markers: all 5 possibilities
# # note: snps in LD should be like all snps for now, since i still need to find a way to filter them
# for (i in 1:5) {
#   kfold_validation_on_sim_traits(snp_data = NULL,
#                                  sv_data = "Structural_Variation_demo.txt",
#                                  snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
#                                  list_of_SV_IDs = SVs,
#                                  sv_info = sv_info,
#                                  marker_data_type = i,
#                                  number_of_folds = 5,
#                                  seed_number = -673994,
#                                  testing = TRUE)
# }
# # RESULT: all markers worked
# 
# 
# # should work for all markers
# # geno data: SNPs and SVs (same file) -- format: num_hmp -- markers: all 5 possibilities
# for (i in 1:5) {
#   kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt",
#                                  sv_data = NULL,
#                                  snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
#                                  list_of_SV_IDs = SVs,
#                                  sv_info = sv_info,
#                                  marker_data_type = i,
#                                  number_of_folds = 5,
#                                  seed_number = -673994,
#                                  testing = TRUE)
# }
# # RESULT: all markers worked
# 
# 
# # should work for all markers -- and number 5 should not merge
# # geno data: SNPs and SVs (same file) -- format: num_hmp -- markers: all 5 possibilities
# for (i in 1:5) {
#   kfold_validation_on_sim_traits(snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt",
#                                  sv_data = "Structural_Variation_demo.txt",
#                                  snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt",
#                                  list_of_SV_IDs = SVs,
#                                  sv_info = sv_info,
#                                  marker_data_type = i,
#                                  number_of_folds = 5,
#                                  seed_number = -673994,
#                                  testing = TRUE)
# }
# # RESULT: all markers worked



#### testing k-fold validation by filtering SNPs ----

# the following tests will only go until the filtering part (i.e., i did not run the k-fold
# validation per se). Just went line by line of the function until the sampling loop.

# test sampling condition
{
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = NULL
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = NULL
  marker_data_type = 1
  number_of_folds = 5
  use_all_SNPs = TRUE    # this skips sampling markers
  SNPs_to_sample = NULL
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

{
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = NULL
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = NULL
  marker_data_type = 1
  number_of_folds = 5
  use_all_SNPs = FALSE   # this gives an error because SNPs_to_sample needs to have a value
  SNPs_to_sample = NULL
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}



# snp_data only -- this should reduce geno_data to 1000 markers
{
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = NULL
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = NULL
  marker_data_type = 1
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

{
  geno_SVs <- fread("Structural_Variation_demo.txt", header = TRUE, data.table = F)
  SVs <- geno_SVs[,1]
  sv_info <- geno_SVs[,c(1,3,4)]
  sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
  
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = NULL
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = sv_info
  marker_data_type = 2
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

{
  geno_SVs <- fread("Structural_Variation_demo.txt", header = TRUE, data.table = F)
  SVs <- geno_SVs[,1]
  sv_info <- geno_SVs[,c(1,3,4)]
  sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
  
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = NULL
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = sv_info
  marker_data_type = 3
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}


# sv_data only -- this shouldn't do anything to geno_data
{
  snp_data = NULL
  sv_data = "Structural_Variation_demo.txt"
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = NULL
  marker_data_type = 4
  number_of_folds = 5
  use_all_SNPs = TRUE
  SNPs_to_sample = NULL
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

{
  snp_data = NULL
  sv_data = "Structural_Variation_demo.txt"
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = NULL
  marker_data_type = 4
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

# snp_data + sv_data

{
  geno_SVs <- fread("Structural_Variation_demo.txt", header = TRUE, data.table = F)
  SVs <- geno_SVs[,1]
  sv_info <- geno_SVs[,c(1,3,4)]
  sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
  
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = "Structural_Variation_demo.txt"
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = sv_info
  marker_data_type = 1
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

{
  geno_SVs <- fread("Structural_Variation_demo.txt", header = TRUE, data.table = F)
  SVs <- geno_SVs[,1]
  sv_info <- geno_SVs[,c(1,3,4)]
  sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
  
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = "Structural_Variation_demo.txt"
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = sv_info
  marker_data_type = 2
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

{
  geno_SVs <- fread("Structural_Variation_demo.txt", header = TRUE, data.table = F)
  SVs <- geno_SVs[,1]
  sv_info <- geno_SVs[,c(1,3,4)]
  sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
  
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = "Structural_Variation_demo.txt"
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = sv_info
  marker_data_type = 3
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

{
  geno_SVs <- fread("Structural_Variation_demo.txt", header = TRUE, data.table = F)
  SVs <- geno_SVs[,1]
  sv_info <- geno_SVs[,c(1,3,4)]
  sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
  
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = "Structural_Variation_demo.txt"
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = sv_info
  marker_data_type = 4
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

{
  geno_SVs <- fread("Structural_Variation_demo.txt", header = TRUE, data.table = F)
  SVs <- geno_SVs[,1]
  sv_info <- geno_SVs[,c(1,3,4)]
  sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
  
  snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
  sv_data = "Structural_Variation_demo.txt"
  snp_with_sv_data = NULL
  list_of_SV_IDs = NULL
  sv_info = sv_info
  marker_data_type = 5
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  dir_with_sim_traits = "folder_name"
  seed_number = -673994
  testing = TRUE
}

# snp_with_sv_data

{
  geno_SVs <- fread("Structural_Variation_demo.txt", header = TRUE, data.table = F)
  SVs <- geno_SVs[,1]
  sv_info <- geno_SVs[,c(1,3,4)]
  sv_info <- sv_info[order(sv_info[,2], sv_info[,3]),]
  
  snp_data = NULL
  sv_data = NULL
  snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt"
  list_of_SV_IDs = SVs
  sv_info = sv_info
  dir_with_sim_traits = "folder_name"
  marker_data_type = 5
  number_of_folds = 5
  use_all_SNPs = FALSE
  SNPs_to_sample = 1000
  seed_number = -673994
  testing = FALSE
}


snp_data = "SNP55K_maize282_AGPv2_20100513_NUM.txt"
sv_data = "Structural_Variation_demo.txt"
snp_with_sv_data = "SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt"
