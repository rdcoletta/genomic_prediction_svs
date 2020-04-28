#### description ----

# this script plots some results of k-fold validations on simulated traits



#### libraries used ----

library(ggplot2)
library(data.table)




#### function ----

PlotPredictionAccuracies <- function(subset.snps = c("all", "1000", "50"),
                                     qtn.number = c("3", "25", "75"),
                                     heritability = c("0.2", "0.5" ,"0.9"),
                                     dir.with.results = "analysis/test_toy/") {
  
  # save path to home directory
  home.dir <- getwd()
  # change to directory with all predictions
  setwd(dir.with.results)
  
  for (subset.snps in c("all", "1000", "50")) {
    
    # get all folders in home directory
    trait.dirs <- list.dirs(recursive = TRUE)
    
    # use different folder names depending if the all SNPs or only a subset was used for GS
    if (subset.snps == "all") {
      folder.name <- "gs_without_subsetting_markers"
    } else {
      folder.name <- paste0("gs_with_", subset.snps, "_SNPs")
    }
    
    # get only the directories correspondent to the subset above
    curr.subset <- grep(folder.name, x = trait.dirs)
    subset.dirs <- trait.dirs[curr.subset]
    
    for (qtn.number in c("3", "25", "75")) {
      
      for (heritability in c("0.2", "0.5" ,"0.9")) {
        
        # select one number of QTN to look at a time
        curr.qtn <- grep(paste0(qtn.number, "-QTNs_from_*"), x = subset.dirs)
        qtn.dirs <- subset.dirs[curr.qtn]
        
        # select one type of heritability to look at a time
        curr.herit <- grep(paste0("Heritability_", heritability, "$"), x = qtn.dirs)
        herit.dirs <- qtn.dirs[curr.herit]
        
        # create an empty data frame to store values for plotting
        df.plot <- data.frame(matrix(ncol = 4, nrow = 0))
        colnames(df.plot) <- c("var_source", "marker_type", "avg_performance", "avg_SE")
        
        # loop though each folder
        for (curr.dir in herit.dirs) {
          
          # print current source of variation
          var.source <- strsplit(curr.dir, "/")[[1]][2]
          var.source <- strsplit(var.source, "_from_")[[1]][2]
          print(paste("Variation source:", var.source))
          # print current marker type used in GS
          marker.type <- strsplit(curr.dir, "/")[[1]][4]
          marker.type <- strsplit(marker.type, "k-fold_validation_")[[1]][2]
          print(paste("Marker used for GS:", marker.type))
          
          # retrieve files with predition performances
          trait.files <- list.files(path = curr.dir, pattern = "5-fold_CV_Results_Heritability_", full.names = TRUE)
          replicate.files <- trait.files[grep("txt$", x = trait.files)]
          
          # create empty vectors to store results
          r.avg <- c()
          r.SE <- c()
          
          for (curr.file in replicate.files) {
            
            # read in file
            kfold.result <- fread(curr.file, header = TRUE, data.table = F)
            # get average correlation
            r.avg <- append(r.avg, kfold.result[,6])
            # get std error. Attention that i need to transform from std dev to std error, which is
            # std_dev/sqrt(n)), where n is the number of fold validations (5 in this case)
            r.SE <- append(r.SE, kfold.result[,7]/sqrt(5))
            
          }
          
          # calculate average
          mean.r.avg <- mean(r.avg)
          # calculate std error
          mean.r.SE <- mean(r.SE)
          
          # append values to data frame
          df.plot[NROW(df.plot)+1, ] <- list(var.source, marker.type, mean.r.avg, mean.r.SE)
          
        }
        
        # reorder levels of marker type for better visualization
        df.plot$marker_type <- factor(df.plot$marker_type, levels = c("all-SNPs", "SNPs-LD",
                                                                      "snps-varying-LD", "only-SVs",
                                                                      "SNPs-and-SVs"))
        
        # add information to the plot title about which kind of subset markers were used
        if (subset.snps == "all") {
          subset.plot.title <- ""
        } else {
          subset.plot.title <- paste0(", ", subset.snps, "SNPs in GS")
        }
        
        # create the plot with the info on data frame
        ggplot(data = df.plot, aes(x = marker_type, y = avg_performance, fill = marker_type)) + 
          geom_bar(stat = "identity", show.legend = FALSE) +
          facet_grid(cols = vars(var_source)) +
          geom_errorbar(aes(ymin = avg_performance - avg_SE, ymax = avg_performance + avg_SE),
                        width=.2) +
          scale_x_discrete(drop = FALSE) +
          ylim(0,1) + 
          xlab("Marker type used in prediciton") +
          ylab("Average performance") +
          ggtitle(paste0(qtn.number, "QTNs, ", heritability, " heritability", subset.plot.title)) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          # scale_fill_brewer(type = "div", palette = "RdYlBu")
          scale_colour_viridis_d(option = "B", begin = 0.3, end = 0.8, aesthetics = "fill",
                                 drop = FALSE)
        
        # save plot
        if (subset.snps == "all") {
          plot.name <- paste0("plot_avg-perform_without-subsetting_", heritability, "-herit_",
                              qtn.number, "-QTNs.png")
        } else {
          plot.name <- paste0("plot_avg-perform_with-", subset.snps, "-SNPs_", heritability, "-herit_",
                              qtn.number, "-QTNs.png")
        }
        ggsave(filename = plot.name)
      }
      
    }
    
  }
  
  # change back to home directory
  setwd(home.dir)
}


#### plots toy dataset ----

# PlotPredictionAccuracies(subset.snps = c("all", "1000", "50"),
#                          qtn.number = c("3", "25", "75"),
#                          heritability = c("0.2", "0.5" ,"0.9"),
#                          dir.with.results = "tests/analysis/test_toy/")



#### plots usda dataset (no SV info yet! only GS with SNPs) ----

PlotPredictionAccuracies(subset.snps = c("all", "1000", "50"),
                         qtn.number = c("3", "25", "75"),
                         heritability = c("0.2", "0.5" ,"0.9"),
                         dir.with.results = "analysis/snps-svs_rils")

# this is just for the preliminary data
# once I have all SV data, I will run the PlotPredictionAccuracies() function above
# for now,  i will run modified commands from the function above


# save path to home directory
home.dir <- getwd()
# change to directory with all predictions
setwd("analysis/snps-svs_rils")

for (subset.snps in c("all", "1000", "50")) {
  
  # get all folders in home directory
  trait.dirs <- list.dirs(recursive = TRUE)
  
  # use different folder names depending if the all SNPs or only a subset was used for GS
  if (subset.snps == "all") {
    folder.name <- "gs_without_subsetting_markers"
  } else {
    folder.name <- paste0("gs_with_", subset.snps, "_SNPs")
  }
  
  # get only the directories correspondent to the subset above
  curr.subset <- grep(folder.name, x = trait.dirs)
  subset.dirs <- trait.dirs[curr.subset]
  
  # create an empty data frame to store values for plotting
  df.plot <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(df.plot) <- c("var_source", "marker_type", "qtn_number", "heritability", "avg_performance", "avg_SE")
  
  for (qtn.number in c("3", "25", "75")) {
    
    for (heritability in c("0.2", "0.5" ,"0.9")) {
      
      # select one number of QTN to look at a time
      curr.qtn <- grep(paste0(qtn.number, "-QTNs_from_*"), x = subset.dirs)
      qtn.dirs <- subset.dirs[curr.qtn]
      
      # select one type of heritability to look at a time
      curr.herit <- grep(paste0("Heritability_", heritability, "$"), x = qtn.dirs)
      herit.dirs <- qtn.dirs[curr.herit]
      
      # loop though each folder
      for (curr.dir in herit.dirs) {
        
        # print current source of variation
        var.source <- strsplit(curr.dir, "/")[[1]][2]
        var.source <- strsplit(var.source, "_from_")[[1]][2]
        print(paste("Variation source:", var.source))
        # print current marker type used in GS
        marker.type <- strsplit(curr.dir, "/")[[1]][4]
        marker.type <- strsplit(marker.type, "k-fold_validation_")[[1]][2]
        print(paste("Marker used for GS:", marker.type))
        
        # retrieve files with predition performances
        trait.files <- list.files(path = curr.dir, pattern = "5-fold_CV_Results_Heritability_", full.names = TRUE)
        replicate.files <- trait.files[grep("txt$", x = trait.files)]
        
        # create empty vectors to store results
        r.avg <- c()
        r.SE <- c()
        
        for (curr.file in replicate.files) {
          
          # read in file
          kfold.result <- fread(curr.file, header = TRUE, data.table = F)
          # get average correlation
          r.avg <- append(r.avg, kfold.result[,6])
          # get std error. Attention that i need to transform from std dev to std error, which is
          # std_dev/sqrt(n)), where n is the number of fold validations (5 in this case)
          r.SE <- append(r.SE, kfold.result[,7]/sqrt(5))
          
        }
        
        # calculate average
        mean.r.avg <- mean(r.avg)
        # calculate std error
        mean.r.SE <- mean(r.SE)
        
        # append values to data frame
        df.plot[NROW(df.plot)+1, ] <- list(var.source, marker.type, qtn.number, heritability, mean.r.avg, mean.r.SE)
        
      }
      
    }
    
  }
  
  # reorder levels of marker type for better visualization
  df.plot$qtn_number <- factor(df.plot$qtn_number, levels = c("3", "25", "75"))
  
  # add plot title
  if (subset.snps == "all") {
    subset.plot.title <- "All SNPs used in GS"
  } else {
    subset.plot.title <- paste("Only", subset.snps, "SNPs used in GS")
  }
  
  # create the plot with the info on data frame
  ggplot(data = df.plot, aes(x = heritability, y = avg_performance, fill = heritability)) + 
    facet_grid(cols = vars(qtn_number)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_errorbar(aes(ymin = avg_performance - avg_SE, ymax = avg_performance + avg_SE),
                  width=.2) +
    ylim(0,1) + 
    xlab("Heritability") +
    ylab("Average performance") +
    ggtitle(subset.plot.title) +
    scale_fill_manual(values = c("#810f7c", "#8856a7", "#8c96c6"))
  
  # save plot
  if (subset.snps == "all") {
    plot.name <- paste0("plot_avg-perform_without-subsetting_usda-rils_SNPs-only.png")
  } else {
    plot.name <- paste0("plot_avg-perform_with-", subset.snps, "-SNPs_usda-rils_SNPs-only.png")
  }
  ggsave(filename = plot.name)
  
}

# change back to home directory
setwd(home.dir)


