#### description ----

# this script plots some results of k-fold validations on simulated traits



#### libraries used ----

library(ggplot2)
library(data.table)



#### plots toy dataset ----

setwd("/Users/Della/Documents/Rafael/OneDrive/University of Minnesota/PhD/hirsch_lab/projects/genomic_prediction/lipka_lab_visit/Simulation_SV/trait-sim_manuscript/one-trait_one-env/")
setwd("/Users/rafael/OneDrive/University of Minnesota/PhD/hirsch_lab/projects/genomic_prediction/lipka_lab_visit/Simulation_SV/trait-sim_manuscript/one-trait_one-env/")

trait_dirs <- list.dirs(recursive = TRUE)

for (qtn_number in c("3", "15", "30")) {
  
  for (heritability in c("0.2", "0.5" ,"0.9")) {
    
    # select one number of QTN to look at a time
    curr_qtn <- grep(paste0(qtn_number, "-QTNs_from_*"), x = trait_dirs)
    qtn_dirs <- trait_dirs[curr_qtn]
    
    # select one type of heritability to look at a time
    curr_herit <- grep(paste0("Heritability_", heritability, "$"), x = qtn_dirs)
    herit_dirs <- qtn_dirs[curr_herit]
    
    # create an empty data frame to store values for plotting
    df_plot <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(df_plot) <- c("var_source", "marker_type", "avg_performance", "avg_SE")
    
    # loop though each folder
    for (curr_dir in herit_dirs) {

      # print current source of variation
      var_source <- strsplit(curr_dir, "/")[[1]][2]
      var_source <- strsplit(var_source, "_from_")[[1]][2]
      print(paste("Variation source:", var_source))
      # print current marker type used in GS
      marker_type <- strsplit(curr_dir, "/")[[1]][3]
      marker_type <- strsplit(marker_type, "k-fold_validation_")[[1]][2]
      print(paste("Marker used for GS:", marker_type))
      
      # retrieve files with predition performances
      trait_files <- list.files(path = curr_dir, pattern = "5-fold_CV_Results_Heritability_", full.names = TRUE)
      replicate_files <- trait_files[grep("txt$", x = trait_files)]
      
      # create empty vectors to store results
      r_avg <- c()
      r_std <- c()
      
      for (curr_file in replicate_files) {
        
        # read in file
        kfold_result <- fread(curr_file, header = TRUE, data.table = F)
        # get average correlation
        r_avg <- append(r_avg, kfold_result[,6])
        # get std error. Attention that i need to transform from std dev to std error, which is
        # std_dev/sqrt(n)), where n is the number of fold validations (5 in this case)
        r_SE <- append(r_std, kfold_result[,7]/sqrt(5))
        
      }
      
      # calculate average
      mean_r_avg <- mean(r_avg)
      # calculate std error
      mean_r_SE <- mean(r_SE)
      
      # append values to data frame
      df_plot[NROW(df_plot)+1, ] <- list(var_source, marker_type, mean_r_avg, mean_r_SE)

    }
    
    # reorder levels of marker type for better visualization
    df_plot$marker_type <- factor(df_plot$marker_type, levels = c("all-SNPs", "SNPs-LD", "snps-varying-LD",
                                                                  "only-SVs", "SNPs-and-SVs"))
    
    # create the plot with the info on data frame
    ggplot(data = df_plot, aes(x = marker_type, y = avg_performance, fill = marker_type)) + 
      geom_bar(stat = "identity", show.legend = FALSE) +
      facet_grid(cols = vars(var_source)) +
      geom_errorbar(aes(ymin = avg_performance - avg_SE, ymax = avg_performance + avg_SE),
                    width=.2) +
      ylim(0,1) + 
      xlab("Marker type used in prediciton") +
      ylab("Average performance") +
      ggtitle(paste(qtn_number, "QTNs,", heritability, "heritability")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      # scale_fill_brewer(type = "div", palette = "RdYlBu")
      scale_colour_viridis_d(option = "B", begin = 0.3, end = 0.8, aesthetics = "fill")
    
    # save plot
    plot_name <- paste0("plot_avg-perform_", heritability, "-herit_", qtn_number, "-QTNs.png")
    ggsave(filename = plot_name)
  }
  
}



#### plots usda dataset (no SV info yet! only GS with SNPs) ----

# this is just for the preliminary data
# once I have all SV data, I will run the code above

setwd("/Users/Della/Documents/Rafael/OneDrive/University of Minnesota/PhD/hirsch_lab/projects/genomic_prediction/lipka_lab_visit/Simulation_SV/trait-sim_manuscript/test_usda-rils/")
setwd("/Users/rafael/OneDrive/University of Minnesota/PhD/hirsch_lab/projects/genomic_prediction/lipka_lab_visit/Simulation_SV/trait-sim_manuscript/test_usda-rils/")

trait_dirs <- list.dirs(recursive = TRUE)

# create an empty data frame to store values for plotting
df_plot <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df_plot) <- c("var_source", "marker_type", "qtn_number", "heritability", "avg_performance", "avg_SE")

for (qtn_number in c("3", "15", "30")) {
  
  for (heritability in c("0.2", "0.5" ,"0.9")) {
    
    # select one number of QTN to look at a time
    curr_qtn <- grep(paste0(qtn_number, "-QTNs_from_*"), x = trait_dirs)
    qtn_dirs <- trait_dirs[curr_qtn]
    
    # select one type of heritability to look at a time
    curr_herit <- grep(paste0("Heritability_", heritability, "$"), x = qtn_dirs)
    herit_dirs <- qtn_dirs[curr_herit]
    
    # loop though each folder
    for (curr_dir in herit_dirs) {
      
      # print current source of variation
      var_source <- strsplit(curr_dir, "/")[[1]][2]
      var_source <- strsplit(var_source, "_from_")[[1]][2]
      print(paste("Variation source:", var_source))
      # print current marker type used in GS
      marker_type <- strsplit(curr_dir, "/")[[1]][3]
      marker_type <- strsplit(marker_type, "k-fold_validation_")[[1]][2]
      print(paste("Marker used for GS:", marker_type))
      
      # retrieve files with predition performances
      trait_files <- list.files(path = curr_dir, pattern = "5-fold_CV_Results_Heritability_", full.names = TRUE)
      replicate_files <- trait_files[grep("txt$", x = trait_files)]
      
      # create empty vectors to store results
      r_avg <- c()
      r_std <- c()
      
      for (curr_file in replicate_files) {
        
        # read in file
        kfold_result <- fread(curr_file, header = TRUE, data.table = F)
        # get average correlation
        r_avg <- append(r_avg, kfold_result[,6])
        # get std error. Attention that i need to transform from std dev to std error, which is
        # std_dev/sqrt(n)), where n is the number of fold validations (5 in this case)
        r_SE <- append(r_std, kfold_result[,7]/sqrt(5))
        
      }
      
      # calculate average
      mean_r_avg <- mean(r_avg)
      # calculate std error
      mean_r_SE <- mean(r_SE)
      
      # append values to data frame
      df_plot[NROW(df_plot)+1, ] <- list(var_source, marker_type, qtn_number, heritability, mean_r_avg, mean_r_SE)
      
    }
    
  }
  
}

# reorder levels of marker type for better visualization
df_plot$qtn_number <- factor(df_plot$qtn_number, levels = c("3", "15", "30"))

# create the plot with the info on data frame
ggplot(data = df_plot, aes(x = heritability, y = avg_performance, fill = heritability)) + 
  facet_grid(cols = vars(qtn_number)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_errorbar(aes(ymin = avg_performance - avg_SE, ymax = avg_performance + avg_SE),
                width=.2) +
  ylim(0,1) + 
  xlab("Heritability") +
  ylab("Average performance") +
  scale_fill_manual(values = c("#810f7c", "#8856a7", "#8c96c6"))

# save plot
plot_name <- paste0("plot_avg-perform_usda-rils_SNPs-only.png")
ggsave(filename = plot_name)
