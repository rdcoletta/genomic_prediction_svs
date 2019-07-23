# Functions provided by Samuel Fernandes to necessary to simulate traits.

#### Numericalization v2 ----
`GAPIT.Numericalization_v2` <-
  function(x,
           bit = NULL,
           effect = "Add",
           impute = "None",
           Create.indicator = FALSE,
           Major.allele.zero = FALSE,
           byRow = TRUE) {
    #Object: To convert character SNP genotpe to numerical
    #Output: Coresponding numerical value
    #Authors: Feng Tian and Zhiwu Zhang (Modified by Samuel Fernandes)
    # Last update: April 11, 2019
    #---------------------------------------------------------------------------
    if (bit == 1)  {
      x[x == "X" | x == "-" | x == "+" | x == "/"] = "N"
      #K (for GT genotype)is replaced by Z to ensure heterozygose has the largest value
      x[x == "K"] = "Z"
      if (class(x) != "matrix") {
        x <- as.matrix(x)
      }
      #Genotype counts
      count = table(x)
      lev = setdiff(names(count), "N")
      len = length(lev)
      max.c <- which.max(count)
      min.c <- which.min(count)
      
      if (Major.allele.zero) {
        if (len > 1 & len <= 3) {
          #One bit: Make sure that the SNP with the major allele is on the top,
          #and the SNP with the minor allele is on the second position
          order <- c(max.c, min.c, setdiff(1:len, c(max.c, min.c)))
          count = count[order]
          lev = lev[order]
        }
      } #End  if(Major.allele.zero)
      
      #make two  bit order genotype as AA,AT and TT, one bit as A(AA),T(TT) and X(AT)
      if (len == 3) {
        temp = count[2]
        count[2] = count[3]
        count[3] = temp
      }
      position = order(count)
    }
    
    if (bit == 2)  {
      x[x == "XX" | x == "--" | x == "++" | x == "//" | x == "NN"] = "N"
      
      #Genotype counts
      count = table(x)
      lev = setdiff(names(count), "N")
      len = length(lev)
      
      if (Major.allele.zero & (len > 1 & len <= 3)) {
        max.c <- which.max(count)
        min.c <- which.min(count)
        #Two bit: Make sure that the SNP with the major allele is on the top,
        #and the SNP with the minor allele is on the third position
        order <- c(max.c,  setdiff(1:len, c(max.c, min.c)), min.c)
        count = count[order]
        lev = lev[order]
      } #End  if(Major.allele.zero)
      position = order(count)
    }
    
    #1 status other than 2 or 3
    if (len <= 1 | len > 3) {
      x1 <- rep(0, times = length(x))
    }
    
    #2 status
    if (len == 2) {
      x1 <- 1:length(x)
      x1[x == "N"] <- NA
      x1[x == lev[1]] <- 0
      x1[x != lev[1]] <- 2
    }
    
    #3 status
    if (len == 3) {
      if (bit == 1) {
        x1 <- 1:length(x)
        x1[x == "N"] <- NA
        x1[x != lev[1] & x != lev[3]] <- 2
        x1[x == lev[1]] <- 0
        x1[x == lev[3]] <- 1
      } else{
        x1 <- 1:length(x)
        x1[x == "N"] <- NA
        x1[x != lev[1] & x != lev[3]] <- 1
        x1[x == lev[1]] <- 0
        x1[x == lev[3]] <- 2
      }
    }
    
    #missing data imputation
    if (impute == "Middle") {
      x1[is.na(x1)] = 1
    }
    
    if (len == 3) {
      if (impute == "Minor")  {
        x1[is.na(x1)] = position[1]  - 1
      }
      if (impute == "Major")  {
        x1[is.na(x1)] = position[len] - 1
      }
      
    } else{
      if (impute == "Minor")  {
        x1[is.na(x1)] = 2 * (position[1]  - 1)
      }
      if (impute == "Major")  {
        x1[is.na(x1)] = 2 * (position[len] - 1)
      }
    }
    
    #alternative genetic models
    if (effect == "Dom") {
      x1[x1 == 1] = 1
      x1[x1 != 1] = 0
    }
    if (effect == "Left")
      x1[x1 == 1] = 0
    if (effect == "Right")
      x1[x1 == 1] = 2
    
    if (byRow) {
      result = matrix(x1, length(x1), 1)
    } else{
      result = matrix(x1, 1, length(x1))
    }
    return(result)
  }#end of GAPIT.Numericalization_v2 function

#### HapMap ----
`GAPIT.HapMap` <-
  function(G,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           heading = TRUE,
           Create.indicator = FALSE,
           Major.allele.zero = FALSE) {
    #Object: To convert character SNP genotpe to numerical (calling GAPIT.Numericalization_v2)
    #Output: Coresponding numerical value
    #Authors: Feng Tian and Zhiwu Zhang
    # Last update: May 30, 2011
    #---------------------------------------------------------------------------
    print(paste0(
      "Converting HapMap format to numerical under model of ",
      SNP.impute
    ))
    if (heading) {
      GT = t(G[1, -(1:11)])
      GI = G[-1, c(1, 3, 4)]
    } else{
      GT = NULL
      GI = G[, c(1, 3, 4)]
    }
    
    #Set column names
    if (heading)
      colnames(GT) = "taxa"
    colnames(GI) = c("SNP", "Chromosome", "Position")
    
    #Initial GD
    GD = NULL
    #to determine number of bits of genotype
    bit = nchar(as.character(G[2, 12]))
    print("Performing numericalization")
    if (heading) {
      if (!Create.indicator)
        GD = apply(G[-1, -(1:11)], 1, function(one)
          GAPIT.Numericalization_v2(
            one,
            bit = bit,
            effect = SNP.effect,
            impute = SNP.impute,
            Major.allele.zero = Major.allele.zero
          ))
      if (Create.indicator)
        GD = t(G[-1, -(1:11)])
    } else{
      if (!Create.indicator)
        GD = apply(G[, -(1:11)], 1, function(one)
          GAPIT.Numericalization_v2(
            one,
            bit = bit,
            effect = SNP.effect,
            impute = SNP.impute,
            Major.allele.zero = Major.allele.zero
          ))
      if (Create.indicator)
        GD = t(G[, -(1:11)])
    }
    
    #set GT and GI to NULL in case of null GD
    if (is.null(GD)) {
      GT = NULL
      GI = NULL
    }
    
    if (!Create.indicator) {
      print(paste0("Succesfuly finished converting HapMap which has bits of ",
                   bit))
    }
    return(list(GT = GT, GD = GD, GI = GI))
  }#end of GAPIT.HapMap function

#### Fragment_v2 ----
`GAPIT.Fragment_v2` <-
  function(file.path = NULL,
           file.from = NULL,
           file.to = NULL,
           file.G = NULL,
           file.Ext.G = NULL,
           seed = NULL,
           SNP.fraction = 1,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           genoFormat = NULL,
           file.GD = NULL,
           file.Ext.GD = NULL,
           file.GM = NULL,
           file.Ext.GM = NULL,
           file.fragment = Inf,
           file = 1,
           frag = 1,
           LD.chromosome = NULL,
           LD.location = NULL,
           LD.range = NULL,
           Create.indicator = FALSE,
           Major.allele.zero = FALSE) {
    #Object: To load SNPs on a (frag)ment in file (this is to replace sampler)
    #Output: genotype data sampled
    #Authors: Alex Lipka and Zhiwu Zhang (Modified by Samuel Fernandes)
    # Last update: April 11, 2019
    #---------------------------------------------------------------------------
    # Using function fread from package data.table for fast reading HMP files
    require(data.table)
    genoFormat = "hapmap"
    if (!is.null(file.GD) & is.null(file.G)) {
      genoFormat = "EMMA"
    }
    
    if (genoFormat == "hapmap") {
      G = NULL
      if (frag == 1) {
        skip.1 = 0
        G <-
          try(fread(
            paste0(file.path, file.G, file, ".", file.Ext.G),
            head = FALSE,
            skip = skip.1,
            nrows = file.fragment + 1,
            na.strings = "NA",
            data.table = F
          ),
          silent = TRUE)
      } else{
        skip.1 <- (frag - 1) * file.fragment + 1
        G <-
          try(fread(
            paste0(file.path, file.G, file, ".", file.Ext.G),
            head = FALSE,
            skip = skip.1,
            nrows = file.fragment
          ),
          silent = TRUE)
      }
      
      if (inherits(G, "try-error"))  {
        G = NULL
        return(list(
          GD = NULL,
          GI = NULL,
          GT = NULL,
          linesRead = NULL,
          GLD = NULL,
          heading = NULL
        ))
      }
      
      #print("Calling hapmap...")
      heading = (frag == 1)
      
      #Recording number of lineas read
      if (heading) {
        n = nrow(G) - 1
      } else{
        n = nrow(G)
      }
      
      linesRead = n
      
      #Sampling
      if (SNP.fraction < 1) {
        #print("Number of SNP in this pragment:")
        #print(n)
        
        #set.seed(seed+(file*1000)+frag)
        if (!is.null(seed)) {
          set.seed(seed)
        }
        #mySample=sample(1:n,max(2,floor(n*as.numeric(as.vector(SNP.fraction)))))
        mySample = sample(1:n, max(2, floor(n * SNP.fraction)))
        
        #print(length(mySample))
        if (heading) {
          G = G[c(1, (1 + mySample)),]
        } else{
          G = G[mySample,]
        }
      } #end of if(SNP.fraction<1)
      
      hm = GAPIT.HapMap(
        G,
        SNP.effect = SNP.effect,
        SNP.impute = SNP.impute,
        heading = heading,
        Create.indicator = Create.indicator,
        Major.allele.zero = Major.allele.zero
      )
      
      #print("Extracting snps for LD plot...")
      #Extract SNPs for LD plot
      if (!is.null(LD.chromosome) & !is.null(hm$GD)) {
        index = (G[, 3] == LD.chromosome[1]) &
          abs((as.numeric(G[, 4]) - as.numeric(LD.location[1])) <
                (as.numeric(LD.range[1]) / 2))
        GLD = G[index,]
      } else{
        GLD = NULL
      }
      
      #rm(G)
      #gc()
      print("hapmap called successfuly from fragment")
      
      return(
        list(
          GD = hm$GD,
          GI = hm$GI,
          GT = hm$GT,
          linesRead = linesRead,
          GLD = GLD,
          heading = heading,
          G = G
        )
      )
      
      print("ERROR: It should not get here!!!")
    } #end of "hapmap"
    
    
    if (genoFormat == "EMMA") {
      #Initial GD
      GD = NULL
      skip.1 <- (frag - 1) * file.fragment
      #Skip the remaining columns
      GD.temp <-
        try(fread(
          paste0(file.path, file.GD, file, ".", file.Ext.GD),
          head = TRUE,
          nrows = 1,
          na.strings = "NA",
          data.table = F
        ),
        silent = TRUE)
      num.SNP <- ncol(GD.temp) - 1
      rm(GD.temp)
      read.in <- min(file.fragment, (num.SNP - skip.1))
      skip.2 <- max((num.SNP - (skip.1 + read.in)), 0)
      
      GD <-
        try(fread(
          paste0(file.path, file.GD, file, ".", file.Ext.GD),
          head = TRUE,
          na.strings = "NA",
          data.table = F,
          colClasses = c(
            "factor",
            rep("NULL", skip.1),
            rep("numeric", read.in),
            rep("NULL", skip.2)
          )
        ) ,
        silent = TRUE)
      GI <-
        try(fread(
          paste0(file.path, file.GM, file, ".", file.Ext.GM),
          head = TRUE,
          na.strings = "NA",
          data.table = F,
          skip = skip.1,
          nrows = file.fragment
        ) ,
        silent = TRUE)
      
      if (inherits(GD, "try-error"))  {
        GD = NULL
        print("File end reached for GD!!!")
      }
      if (inherits(GI, "try-error"))  {
        GI = NULL
        print("File end reached for GI!!!")
      }
      
      if (is.null(GD))
        return(list(
          GD = NULL,
          GI = NULL,
          GT = NULL,
          linesRead = NULL,
          GLD = NULL
        ))
      
      GT = GD[, 1]  #Extract infividual names
      
      GD = GD[,-1] #Remove individual names
      #print("Numerical file read sucesfuly from fragment")
      linesRead = ncol(GD)
      if (SNP.fraction == 1)
        return(list(
          GD = GD,
          GI = GI,
          GT = GT,
          linesRead = linesRead,
          GLD = NULL
        ))
      
      if (SNP.fraction < 1) {
        n = ncol(GD)
        #set.seed(seed+file)
        if (!is.null(seed)) {
          set.seed(seed)
        }
        sample = sample(1:n, floor(n * SNP.fraction))
        return(list(
          GD = GD[, sample],
          GI = GI[sample,],
          GT = GT,
          linesRead = linesRead,
          GLD = NULL
        ))
      }
    } # end of the "EMMA"
    #print("fragment ended succesfully!")
  }#End of GAPIT.Fragment_v2 function

#### Genotypes ----
`Genotypes` <-
  function(file.path = NULL,
           maf_cutoff = NULL,
           seed = 123,
           file.G = NULL,
           file.Ext.G = NULL,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           file.fragment = Inf,
           Create.indicator = FALSE,
           Major.allele.zero = FALSE,
           file.from = 1,
           file.to = 1) {
    #Object: Generate a numeric (dosaje) HapMap file
    #Output: A numeric HapMap
    #Authors: Alex lipka and Samuel Fernandes
    # Last update: April 19, 2019
    #---------------------------------------------------------------------------
    count.filenum <- 0
    for (filenum in file.from:file.to) {
      
      if(!any(grepl(paste0(file.G, filenum,".",file.Ext.G) , dir()))){
        stop(paste("File" , paste0("\'",file.G, filenum,".",file.Ext.G,"\'"), "not in the directory!","\n"))
      }
      
      myFRG = GAPIT.Fragment_v2(
        file.path = NULL,
        file.from = file.from,
        file.to = file.to,
        file.G = file.G,
        file.Ext.G = file.Ext.G,
        seed = seed,
        SNP.effect = SNP.effect,
        SNP.impute = SNP.impute,
        genoFormat = NULL,
        file.GD = NULL,
        file.Ext.GD = NULL,
        file.GM = NULL,
        file.Ext.GM = NULL,
        file = filenum,
        file.fragment = file.fragment,
        LD.chromosome = NULL,
        LD.location = NULL,
        LD.range = NULL,
        Create.indicator = Create.indicator,
        Major.allele.zero = Major.allele.zero
      )
      
      if (count.filenum == 0) {
        all.FRGGs <- myFRG
      } else{
        all.FRGGs$GD <- cbind(all.FRGGs$GD, myFRG$GD)
        all.FRGGs$GI <- rbind(all.FRGGs$GI, myFRG$GI)
      }#End if(count.filenum == 0)
      count.filenum <- count.filenum + 1
    }#end for(filenum in file.from:file.to)
    
    if (!is.null(maf_cutoff)) {
      hm <- list(GT = all.FRGGs$GT,
                 GD = all.FRGGs$GD,
                 GI = all.FRGGs$GI)
      
      #Obtain the mafs of all SNPs
      #-------------------------------------------------------------------------
      #Total number of lines
      ns <- nrow(hm$GD)
      
      #Sum of the allele scores for each SNP
      ss <- apply(hm$GD, 2, sum)
      
      #Combine two situations: one where the allele coded as "2" is major;
      #one where "0" is coded as major.
      maf.matrix <- rbind((.5 * ss / ns), (1 - (0.5 * ss / ns)))
      
      #Copy the minor allele frequencies for all SNPs
      maf <- apply(maf.matrix, 2, min)
      
      #Find out which SNPs have MAF < maf_cutoff
      snps.below.maf <- which(maf < maf_cutoff)
      
      # Remove these SNPs from hm$GD
      
      hm.GD.without.snps.below.maf <- hm$GD[,-snps.below.maf]
      
      genotypes <-
        data.frame(hm$GI[, 1],
                   rep(NA, nrow(hm$GI)),
                   hm$GI[, 2:3],
                   rep(NA, nrow(hm$GI)),
                   t(hm.GD.without.snps.below.maf))
      
      colnames(genotypes) <-
        c("Snp", "allele", "chr", "pos", "cm", t(as.character(hm$GT)))
      
    } else {
      genotypes <-
        data.frame(all.FRGGs$GI[, 1],
                   rep(NA, nrow(all.FRGGs$GI)),
                   all.FRGGs$GI[, 2:3],
                   rep(NA, nrow(all.FRGGs$GI)),
                   t(all.FRGGs$GD))
      colnames(genotypes) <-
        c("Snp", "allele", "chr", "pos", "cm", t(as.character(all.FRGGs$GT)))
    }
    return(genotypes)
  }#End of Genotypes function

#### QTN_pleiotropic ----
`QTN_pleiotropic` <-
  function(genotypes = NULL,
           seed = NULL,
           Additive.QTN.number = NULL,
           Epistatic.QTN.number = NULL) {
    #---------------------------------------------------------------------------
    #Object: Select SNPs to be assigned as QTNs
    #Output: Genotype of selected SNPs
    #Authors: Alex lipka and Samuel Fernandes
    # Last update: April 19, 2019
    #---------------------------------------------------------------------------
    #Randomly select (without replacement) k additive QTN, and assign an effect size
    # using fwrite from data.table package
    require(data.table)
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    vector.of.add.QTN <-
      sample(1:nrow(genotypes), Additive.QTN.number, replace = FALSE)
    Add.QTN.genotypic.information <- genotypes[vector.of.add.QTN, ]
    
    #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
    write.table(
      ifelse(
        !is.null(seed),
        paste0("Here_is_the_seed_number: ", seed),
        "set.seed not assigned"
      ),
      paste0("seed.number.for.",
             Additive.QTN.number,
             "Add.QTN",
             ".txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    fwrite(
      Add.QTN.genotypic.information,
      paste0(
        "Genotypic.information.for.",
        Additive.QTN.number,
        ".Additive.QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    
    #Create a "base line" trait, which is basically just the additive effects;
    #this is what we would see if the heritability of the simulated trait were 1
    additive.effect.trait.object <-
      t(Add.QTN.genotypic.information[, -c(1:5)])
    
    colnames(additive.effect.trait.object) <-
      paste0("Chr_",
             Add.QTN.genotypic.information[, 3],
             "_",
             Add.QTN.genotypic.information[, 4])
    
    if ( !is.null(Epistatic.QTN.number) ){
      #Randomly select (without replacement) 2*k epistatic QTN, and assign an effect size
      if (!is.null(seed)) {
        set.seed(seed * seed)
      }
      
      vector.of.epi.QTN <-
        sample(1:nrow(genotypes), (2 * Epistatic.QTN.number), replace = FALSE)
      Epi.QTN.genotypic.information <- genotypes[vector.of.epi.QTN, ]
      
      #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
      write.table(
        ifelse(
          !is.null(seed),
          paste0("Here_is_the_seed_number: ", seed * seed),
          "set.seed not assigned"
        ),
        paste0("seed.number.for.",
               Epistatic.QTN.number,
               "Epi.QTN",
               ".txt"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      fwrite(
        Epi.QTN.genotypic.information,
        paste0(
          "Genotypic.information.for.",
          Epistatic.QTN.number,
          ".Epistatic.QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
      
      epistatic.effect.trait.object <-
        t(Epi.QTN.genotypic.information[, -c(1:5)])
      
      colnames(epistatic.effect.trait.object) <-
        paste0("Chr_",
               Epi.QTN.genotypic.information[, 3],
               "_",
               Epi.QTN.genotypic.information[, 4])
      return(
        list(
          additive.effect.trait.object = list(additive.effect.trait.object),
          epistatic.effect.trait.object = list(epistatic.effect.trait.object)
        )
      )
    } else {
      return(
        list(
          additive.effect.trait.object = list(additive.effect.trait.object)
        )
      )
    }
  }

#### QTN_partially_pleiotropic ----
`QTN_partially_pleiotropic` <-
  function(genotypes = NULL,
           seed = NULL,
           overlap = NULL,
           overlapE = NULL,
           specific.QTN.number = NULL,
           specific.E.QTN.number = NULL, 
           ntraits=NULL) {
    #---------------------------------------------------------------------------
    #Object: Select SNPs to be assigned as QTNs
    #Output: Genotype of selected SNPs
    #Authors: Alex lipka and Samuel Fernandes
    # Last update: April 19, 2019
    #---------------------------------------------------------------------------
    #Randomly select (without replacement) k additive QTN, and assign an effect size
    # using fwrite from data.table package
    require(data.table)
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    vector.of.pleiotropic.add.QTN <-
      sample(1:nrow(genotypes), overlap, replace = FALSE)
    
    Add.pleiotropic.QTN.genotypic.info <- genotypes[vector.of.pleiotropic.add.QTN, ]
    
    fwrite(
      Add.pleiotropic.QTN.genotypic.info,
      paste0(
        "Genotypic.information.for.",
        overlap,
        ".pleiotropic.Additive.QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    
    snps <- setdiff(1:nrow(genotypes),vector.of.pleiotropic.add.QTN )
    vector.of.specific.add.QTN <- list()
    Add.specific.QTN.genotypic.info <- list()
    ss <- c()
    for(i in 1:ntraits){
      if (!is.null(seed)) {
        ss[i] <- seed+i
        set.seed(seed+i)
      }
      vector.of.specific.add.QTN[[i]] <-
        sample(snps, specific.QTN.number[i], replace = FALSE)
      
      snps <- setdiff(snps,vector.of.specific.add.QTN[[i]])
      
      Add.specific.QTN.genotypic.info[[i]] <- genotypes[vector.of.specific.add.QTN[[i]],]
      
      fwrite(
        Add.specific.QTN.genotypic.info[[i]],
        paste0(
          "Trait", i, ".Additive.Genotypic.info.for.",specific.QTN.number[i],
          ".specific.QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
    }
    
    #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
    write.table(
      c(seed, ss),
      paste0("seed.number.for.",
             paste0(specific.QTN.number+overlap, collapse = "_"),
             "Add.QTN",
             ".txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    
    if(!is.null(overlapE) | !is.null(specific.E.QTN.number)){
      #Randomly select (without replacement) 2*k epistatic QTN, and assign an effect size
      if (!is.null(seed)) {
        set.seed(seed * seed)
      }
      
      vector.of.pleiotropic.epi.QTN <-
        sample(1:nrow(genotypes), (2 * overlapE), replace = FALSE)
      Epi.pleiotropic.QTN.genotypic.info <- genotypes[vector.of.pleiotropic.epi.QTN, ]
      
      #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
      
      fwrite(
        Epi.pleiotropic.QTN.genotypic.info,
        paste0(
          "Genotypic.information.for.",
          overlapE,
          ".pleiotropic.Epistatic.QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
      
      snpse <- setdiff(1:nrow(genotypes),vector.of.pleiotropic.epi.QTN )
      vector.of.specific.epi.QTN <- list()
      Epi.specific.QTN.genotypic.info <- list()
      sse <- c()
      for(i in 1:ntraits){
        if (!is.null(seed)) {
          sse[i] <- (seed+i)*seed
          set.seed((seed+i)*seed)
        }
        vector.of.specific.epi.QTN[[i]] <-
          sample(snps, (2 * specific.E.QTN.number[i]) , replace = FALSE)
        
        snps <- setdiff(snps,vector.of.specific.epi.QTN[[i]])
        
        Epi.specific.QTN.genotypic.info[[i]] <- genotypes[vector.of.specific.epi.QTN[[i]],]
        
        fwrite(
          Epi.specific.QTN.genotypic.info[[i]],
          paste0(
            "Trait", i, ".Epistatic.Genotypic.info.for.",specific.E.QTN.number[i],
            ".specific.QTN",
            ".txt"
          ),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
      
      write.table(
        c(seed * seed,sse),
        paste0("seed.number.for.",
               paste0(specific.E.QTN.number+overlapE, collapse = "_"),
               "pleiotropic.Epi.QTN",
               ".txt"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      
      #Create a "base line" trait, which is basically just the additive effects;
      #this is what we would see if the heritability of the simulated trait were 1
      additive.effect.trait.object <- list()
      epistatic.effect.trait.object <- list()
      for(i in 1:ntraits){
        additive.effect.trait.object[[i]] <- t(rbind(Add.pleiotropic.QTN.genotypic.info, Add.specific.QTN.genotypic.info[[i]])[, -c(1:5)])
        epistatic.effect.trait.object[[i]] <- t(rbind(Epi.pleiotropic.QTN.genotypic.info, Epi.specific.QTN.genotypic.info[[i]])[, -c(1:5)])
        
        colnames(additive.effect.trait.object[[i]]) <-
          c(paste0("Chr_",
                   Add.pleiotropic.QTN.genotypic.info[, 3],
                   "_",
                   Add.pleiotropic.QTN.genotypic.info[, 4]),
            paste0("Chr_",
                   Add.specific.QTN.genotypic.info[[i]][, 3],
                   "_",
                   Add.specific.QTN.genotypic.info[[i]][, 4]) )
        
        colnames(epistatic.effect.trait.object[[i]]) <-
          c(paste0("Chr_",
                   Epi.pleiotropic.QTN.genotypic.info[, 3],
                   "_",
                   Epi.pleiotropic.QTN.genotypic.info[, 4]),
            paste0("Chr_",
                   Epi.specific.QTN.genotypic.info[[i]][, 3],
                   "_",
                   Epi.specific.QTN.genotypic.info[[i]][, 4]) )
      }
      
      return(
        list(
          additive.effect.trait.object = additive.effect.trait.object,
          epistatic.effect.trait.object = epistatic.effect.trait.object
        )
      )
    } else{
      
      #Create a "base line" trait, which is basically just the additive effects;
      #this is what we would see if the heritability of the simulated trait were 1
      additive.effect.trait.object <- list()
      
      for(i in 1:ntraits){
        additive.effect.trait.object[[i]] <- t(rbind(Add.pleiotropic.QTN.genotypic.info, Add.specific.QTN.genotypic.info[[i]])[, -c(1:5)])
        
        colnames(additive.effect.trait.object[[i]]) <-
          c(paste0("Chr_",
                   Add.pleiotropic.QTN.genotypic.info[, 3],
                   "_",
                   Add.pleiotropic.QTN.genotypic.info[, 4]),
            paste0("Chr_",
                   Add.specific.QTN.genotypic.info[[i]][, 3],
                   "_",
                   Add.specific.QTN.genotypic.info[[i]][, 4]) )
      }
      
      return(
        list(
          additive.effect.trait.object = additive.effect.trait.object
        )
      )
    }  
  }

#### Base_line_single_trait ----
`Base_line_single_trait` <-
  function(additive.object = NULL,
           epistatic.object = NULL,
           additive.effect = NULL,
           epistatic.effect = NULL,
           big.additive.QTN.effect = NULL,
           seed = NULL) {
    #Object: Calculate genetic value based on QTN objects
    #Output: A vector of Genetic values
    #Authors: Alex Lipka and  Samuel Fernandes
    # Last update: April 19, 2019
    #---------------------------------------------------------------------------
    additive.object <- as.data.frame(additive.object)
    
    #make base.line.trait additive.component and epistatic.component
    additive.component <-
      as.data.frame(matrix(0,
                           nrow = nrow(additive.object),
                           ncol = 1))
    
    additive.component <- additive.component +
      (additive.object[, 1] * (big.additive.QTN.effect))
    
    addNumber <- ncol(additive.object)
    if (addNumber >= 2) {
      for (i in 2:addNumber) {
        additive.component <- additive.component +
          (additive.object[, i] * (additive.effect[1] ^ (i - 1)))
      }
    }#end if(addNumber >= 2)
    
    rownames(additive.component) <- rownames(additive.object)
    colnames(additive.component) <- "Additive.effect"
    additive.genetic.variance <- var(additive.component)
    
    if(!is.null(epistatic.object)){
      epistatic.object <- as.data.frame(epistatic.object)
      
      epistatic.component <-
        as.data.frame(matrix(0,
                             nrow = nrow(epistatic.object),
                             ncol = 1))
      
      eNumber <- ncol(epistatic.component)
      for (i in 0:(eNumber - 1)) {
        epistatic.component <-
          epistatic.component +
          ((epistatic.object[, ((2 * i) + 1)] *
              epistatic.object[, ((2 * i) + 2)]) *
             (epistatic.effect ^ (i + 1)))
      }
      rownames(epistatic.component) <- rownames(epistatic.object)
      colnames(epistatic.component) <- "Epistatic.effect"
      epistatic.genetic.variance <- var(epistatic.component)
      
      base.line.trait <- additive.component + epistatic.component
      
      return(list(
        base.line = base.line.trait,
        VA = c(additive.genetic.variance),
        VE = c(epistatic.genetic.variance)
      ))
    } else {
      
      base.line.trait <- additive.component
      
      return(list(
        base.line = base.line.trait,
        VA = c(additive.genetic.variance)
      ))
    }
  }

#### Base_line_multi_traits ----
`Base_line_multi_traits` <-
  function(additive.object = NULL,
           epistatic.object = NULL,
           additive.effect = NULL,
           epistatic.effect = NULL,
           big.additive.QTN.effect = NULL,
           seed = seed,
           set.cor = TRUE,
           ntraits = NULL,
           correlation = NULL,
           model = "pleiotropic"){
    #Object: Calculate genetic value based on QTN objects
    #Output: A matrix of Genetic values for multiple traits
    #Authors: Samuel Fernandes
    # Last update: May 15, 2019
    #---------------------------------------------------------------------------
    if (set.cor) {
      if(model=="pleiotropic"){
        Genetic_value <- Base_line_single_trait(
          additive.object = additive.object[[1]],
          epistatic.object = epistatic.object[[1]],
          additive.effect = additive.effect[1],
          epistatic.effect = epistatic.effect[1],
          big.additive.QTN.effect = big.additive.QTN.effect[1],
          seed = seed
        )
        
        n <- nrow(Genetic_value$base.line)
        varg <- var(Genetic_value$base.line)
        traits <-  matrix(NA, n, ncol = ntraits)
        traits[, 1] <- Genetic_value$base.line$Additive.effect
        
        for(i in 2:ntraits){
          traits[,i] <- Genetic_value$base.line$Additive.effect +  rnorm(n, 0,sqrt(varg*0.01))
        }
        
        Genetic_s<- apply(traits, 2, scale)
        
        # whiting transformation 
        L <- t(chol(cov(Genetic_s)))
        G_white <- t(solve(L)%*%t(Genetic_s))
        
        # coloring transformation 
        #approximation to make it positive definite
        if(!is.positive.definite(correlation)){
          cat("Correlation matrix not positive definite! Using make.positive.definite \n")
          correlation <- lqmm::make.positive.definite(correlation)
        }
        L <- t(chol(correlation))
        traits <- t(L%*%t(G_white))
        
        rownames(traits) <- rownames(additive.object[[1]])
        
        #bring it back to the original scale
        cor_original_trait <- c()
        for(i in 1:ncol(traits)){
          traits[,i] <- (traits[,i]*sd(Genetic_value$base.line$Additive.effect))+mean(Genetic_value$base.line$Additive.effect)
          
          cor_original_trait[i] <- cor( traits[,i], Genetic_value$base.line$Additive.effect)
        }
        sample.cor <- cor(traits)
        
        results <- if(!is.null(epistatic.object)){list(
          base.line = traits,
          VA = c(Genetic_value$VA),
          VE = c(Genetic_value$VE),
          sample.cor = sample.cor
        )} else {
          list(
            base.line = traits,
            VA = c(Genetic_value$VA),
            sample.cor = sample.cor
          )
        }
        return(results)
        
      }else{
        Genetic_value <-
          matrix(NA, nrow(additive.object[[1]]), ncol = ntraits)
        VA <- c()
        VE <- c()
        for (j in 1:ntraits) {
          trait_temp <-  Base_line_single_trait(
            additive.object = additive.object[[j]],
            epistatic.object = epistatic.object[[j]],
            additive.effect = additive.effect[j],
            epistatic.effect = epistatic.effect[j],
            big.additive.QTN.effect = big.additive.QTN.effect[j],
            seed = seed
          )
          Genetic_value[, j] <- trait_temp$base.line$Additive.effect
          VA[j] <- trait_temp$VA
          if(!is.null(epistatic.object)){VE[j] <- trait_temp$VE}
        }
        sdg <- apply(Genetic_value, 2, sd)
        meang <- apply(Genetic_value, 2, mean)
        Genetic_s<- apply(Genetic_value, 2, scale)
        
        # whiting transformation 
        L <- t(chol(cov(Genetic_s)))
        G_white <- t(solve(L)%*%t(Genetic_s))
        
        # coloring transformation 
        if(!is.positive.definite(correlation)){
          cat("Correlation matrix not positive definite! \nUsing lqmm::make.positive.definite \n")
          correlation <- lqmm::make.positive.definite(correlation)
        }
        L <- t(chol(correlation))
        traits <- t(L%*%t(G_white))
        
        rownames(traits) <- rownames(additive.object[[1]])
        #bring it back to the original scale
        cor_original_trait <- c()
        for(i in 1:ncol(traits)){
          traits[,i] <- (traits[,i]*sdg[i])+meang[i]
          cor_original_trait[i] <- cor( traits[,i], Genetic_value[,i])
        }
        
        sample.cor <- cor(traits)
        results <-  if(!is.null(epistatic.object)){list(
          base.line = traits,
          VA = VA,
          VE = VE,
          sample.cor = sample.cor,
          cor_original_trait = cor_original_trait
        )}else{
          list(
            base.line = traits,
            VA = VA,
            sample.cor = sample.cor,
            cor_original_trait = cor_original_trait
          )
        }
        return(results)
      }
      
    } else {
      traits <-
        matrix(NA, nrow(additive.object[[1]]), ncol = ntraits)
      VA <- c()
      VE <- c()
      for (i in 1:ntraits) {
        ifelse(model=="pleiotropic",j<- 1 , j <- i)
        trait_temp <-  Base_line_single_trait(
          additive.object = additive.object[[j]],
          epistatic.object = epistatic.object[[j]],
          additive.effect = additive.effect[i],
          epistatic.effect = epistatic.effect[i],
          big.additive.QTN.effect = big.additive.QTN.effect[i],
          seed = seed
        )
        traits[, i] <- trait_temp$base.line[, 1]
        VA[i] <- trait_temp$VA
        if(!is.null(epistatic.object)){VE[i] <- trait_temp$VE}
      }
      rownames(traits) <- rownames(additive.object[[1]])
      sample.cor <- cor(traits)
      
      results <- if(!is.null(epistatic.object)){
        list(
          base.line = traits,
          VA = VA,
          VE = VE, 
          sample.cor = sample.cor
        )} else {
          list(
            base.line = traits,
            VA = VA,
            sample.cor = sample.cor
          )
        }
      return(results)
    }
  }


#### Phenotypes ----
`Phenotypes` <-
  function(base.line.trait = NULL,
           h2 = NULL,
           rep = NULL,
           seed = NULL,
           ntraits = NULL,
           h2_MT = NULL,
           format = "multi-file") {
    #Object: Generate environmental effects based on a given heritability
    #Output: Phenotypes for ntraits traits
    #Authors: Alex Lipka and  Samuel Fernandes
    # Last update: April 19, 2019
    #---------------------------------------------------------------------------
    require(data.table)
    if (ntraits > 1) {
      if (format == "multi-file") {
        #Create a working directory for the output results:
        dir.create("Phenotypes")
        
        #Set the working directory
        setwd("./Phenotypes")
      }
      #For loop through the vector of heritabilities
      for (i in h2) {
        simulated.data <- vector("list", rep)
        ss <- c()
        
        #If heritability is zero
        if (i == 0) {
          #Simulate rep experiments with a given heritability
          for (z in 1:rep) {
            simulated.data[[z]] <-
              data.frame(rownames(base.line.trait),
                         matrix(NA, nrow(base.line.trait), ntraits),
                         rep = z)
            colnames(simulated.data[[z]]) <-
              c("<Taxa>",
                "Target",
                paste0("Trait_", 2:ntraits),
                "Rep")
            
            if (!is.null(seed)) {
              ss <- c(ss, seed + z)
              set.seed(seed + z)
            }
            #using multivariate normal with ntrait variances and 0 covariances
            #for independent residuals
            simulated.data[[z]][, 2:(ntraits + 1)] <-
              mvtnorm::rmvnorm(
                n = nrow(base.line.trait),
                mean = rep(0, ntraits),
                sigma = diag(1, ntraits)
              )
          }
          
          if (format == "multi-file") {
            
            invisible(lapply(1:rep, function (x) {
              fwrite(
                simulated.data[[x]][-(ntraits + 2)],
                paste0("Simulated.Data.", ".Rep", x, ".Herit.", i, ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }))
          } else if (format == "long") {
            simulated.data <-  do.call(rbind, simulated.data)
            fwrite(
              simulated.data,
              paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else {
            x <- simulated.data[[1]][-(ntraits + 2)]
            for (j in 2:rep) {
              x <- cbind(x, simulated.data[[j]][-c(1, (ntraits + 2))])
            }
            fwrite(
              x,
              paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
          
          #Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            ifelse(format == "multi-file", 
                   paste0("../seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt"),
                   paste0("seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt")),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          
        } else {
          #Calcualte V_e, the residual variance
          residual.cov <-
            diag((apply(base.line.trait, 2, var) / c(i, h2_MT)
            ) - apply(base.line.trait, 2, var))
          #Simulate rep experiments with a given heritability
          for (z in 1:rep) {
            simulated.data[[z]] <-
              data.frame(rownames(base.line.trait),
                         matrix(NA, nrow(base.line.trait), ntraits),
                         rep = z)
            
            colnames(simulated.data[[z]]) <-
              c("<Taxa>", "Target", c(paste0("Trait_", 2:ntraits, "_H2_", i), "Rep"))
            
            if (!is.null(seed)) {
              set.seed(round((seed * z * z) * i))
              ss <- c(ss, round((seed * z * z) * i))
            }
            #using multivariate normal with ntrait variances and 0 covariances
            #for independent residuals
            simulated.data[[z]][, 2:(ntraits + 1)] <-
              base.line.trait +
              mvtnorm::rmvnorm(
                n = nrow(base.line.trait),
                mean = rep(0, ntraits),
                sigma = residual.cov
              )
          }
          
          H2 <- sapply(1:rep, function (x) {
            apply(base.line.trait, 2, var) /
              apply(simulated.data[[x]][1:ntraits+1], 2, var)
          })
          
          H2 <- apply(H2, 1, mean)
          cat("\nPopulational heritability: \n")
          HH <- c(i, h2_MT); names(HH) <- names(H2)
          print(HH)
          
          cat(paste("Sample heritability (Average of",rep," replications): \n"))
          print(H2)
          
          
          if (format == "multi-file") {
            invisible(lapply(1:rep, function (x) {
              fwrite(
                simulated.data[[x]][-(ntraits + 2)],
                paste0("Simulated.Data.", ".Rep", x, ".Herit.", i, ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }))
          } else if (format == "long") {
            simulated.data <-  do.call(rbind, simulated.data)
            fwrite(
              simulated.data,
              paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else {
            x <- simulated.data[[1]][-(ntraits + 2)]
            for (j in 2:rep) {
              x <- cbind(x, simulated.data[[j]][-c(1, (ntraits + 2))])
            }
            fwrite(
              x,
              paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
          
          #Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            ifelse(format == "multi-file", 
                   paste0("../seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt"),
                   paste0("seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt")),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          
        }
      }#End for(i in h2)
      return(paste("Files saved in:", getwd()))
      
    } else{
      #For loop through the vector of heritabilities
      for (i in h2) {
        #If heritability is zero
        ss <- c()
        if (i == 0) {
          simulated.data <-
            data.frame(rownames(base.line.trait), matrix(NA, nrow(base.line.trait), rep))
          #Format the output file for the simulated phenotypes
          colnames(simulated.data) <-
            c("<Taxa>",
              paste0("the.normal.random.variables.", 1:rep))
          
          #Simulate rep experiments with a given heritability
          for (j in 1:rep) {
            if (!is.null(seed)) {
              set.seed(seed + j)
              ss <- c(ss, seed + j)
            }
            simulated.data[, j + 1] <-
              rnorm(nrow(base.line.trait),
                    mean = 0,
                    sd = 1)
          }
          
          #Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            paste0("seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          fwrite(
            simulated.data,
            paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          
        } else{
          #Calcualte V_e, the residual variance
          residual.variance <- (var(base.line.trait) / i) - var(base.line.trait)
          
          simulated.data <-
            data.frame(rownames(base.line.trait), matrix(NA, nrow(base.line.trait), rep))
          #Format the output file for the simulated phenotypes
          colnames(simulated.data) <-
            c("<Taxa>", c(paste0(
              "Heritability_", i, "_Rep_", 1:rep
            )))
          
          #Simulate rep experiments with a given heritability
          for (j in 1:rep) {
            if (!is.null(seed)) {
              ss <- c(ss, round((seed * j * j) * i))
              set.seed(round((seed * j * j) * i))
            }
            the.normal.random.variables <-
              rnorm(
                nrow(base.line.trait),
                mean = 0,
                sd = sqrt(residual.variance)
              )
            
            simulated.data[, j + 1] <-
              base.line.trait + the.normal.random.variables
          }
          
          H2 <- mean(as.vector(var(base.line.trait))/
                       apply(simulated.data[,1:rep+1], 2, var))
          cat("\nPopulational heritability: \n")
          #HH <- c(i); names(HH) <- names(H2)
          print(i)
          
          cat(paste("Sample heritability (Average of",rep," replications): \n"))
          print(H2)
          
          #Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            paste0("seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          fwrite(
            simulated.data,
            paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
      }#End for(i in h2)
      
      return(paste("Files saved in:", getwd()))
    }
  }

#### Create simulated data ----
create.simulated.data <-
  function(genotypes = NULL,
           file.G = NULL,
           file.Ext.G = NULL,
           #default for numericalization#
           file.fragment = Inf,
           file.from = 1,
           file.to = 1,
           maf_cutoff = NULL,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           Create.indicator = FALSE,
           Major.allele.zero = FALSE,
           #-------------------------#
           Additive.QTN.number = NULL,
           Epistatic.QTN.number = NULL,
           additive.effect = NULL,
           epistatic.effect = NULL,
           big.additive.QTN.effect = NULL,
           model = c("pleiotropic", "partially", "LD"),
           overlap = NULL,
           overlapE = NULL,
           specific.QTN.number = NULL,
           specific.E.QTN.number = NULL,
           rep = NULL,
           ntraits = 1,
           h2 = NULL,
           h2_MT = NULL,
           set.cor = TRUE,
           correlation = NULL,
           seed = NULL,
           home.dir=getwd(),
           output.dir = NULL,
           format = "multi-file",
           out.geno = FALSE
  ) {
    #Object: Main function for simulation of ntraits phenotypes based on a SNP file
    #Output: Numeric hapmap, selected QTNs, phenotypes for ntraits traits
    #Authors: Alex Lipka and  Samuel Fernandes
    # Last update: April 19, 2019
    #-----------------------------------------------------------------------------
    try({
      setwd(home.dir)
      if (is.null(genotypes)) {
        genotypes <- Genotypes(
          file.path = paste0(home.dir, "/"),
          maf_cutoff = maf_cutoff,
          seed = 123,
          file.G = file.G,
          file.Ext.G = file.Ext.G,
          SNP.effect = SNP.effect,
          SNP.impute = SNP.impute,
          file.fragment = file.fragment,
          Create.indicator = Create.indicator,
          Major.allele.zero = Major.allele.zero,
          file.from = file.from,
          file.to = file.to
        )
      }
      
      #Create a working directory for the output results:
      if(!is.null(output.dir)) {dir.create(paste0(home.dir, "/", output.dir))}
      
      #Set the working directory
      setwd(paste0(home.dir, "/", output.dir, "/"))
      sink("Log_Sim.txt", type = "output")
      
      if(ntraits > 1){
        mm <- ifelse(model=="pleiotropic", "pleiotropic", 
                     ifelse(model=="partially","partially pleiotropic",
                            ifelse(model=="LD", "Linkage Desiquilibrium",
                                   {stop("model used is not valid!", call. = F); sink()}) ))
        
        cat(paste("Simulation of a", mm,"genetic model \n"))
      }
      if(is.null(file.G)){cat(paste("genotype object:", deparse(substitute(genotypes))))
      } else{cat(paste("genotype file:",file.G ,"*",file.Ext.G))}
      
      if(ntraits > 1){
        if(model=="pleiotropic"){
          cat(paste("\nNumber of additive QTNs:", Additive.QTN.number,
                    "\nNumber of epistatic QTNs:", Epistatic.QTN.number, "\n"))
        }
        
        if(model=="partially"){
          if(is.null(overlapE) | is.null(specific.E.QTN.number)){
            cat(paste("\nNumber of pleiotropic additive QTNs:", overlap,
                      "\nNumber of trait specific additive QTNs:", paste(specific.QTN.number, collapse = ", "), "\n \n"))
          }else{
            cat(paste("\nNumber of pleiotropic additive QTNs:", overlap,
                      "\nNumber of trait specific additive QTNs:", paste(specific.QTN.number, collapse = ", "),
                      "\nNumber of pleiotropic epistatic QTNs:", overlapE,
                      "\nNumber of trait specific epistatic QTNs:", paste(specific.E.QTN.number, collapse = ", "), "\n \n"))
          }     
        }
      }
      if(out.geno){
        fwrite(
          genotypes,
          paste0(ifelse(is.null(output.dir),"", "../"), file.G, "NUM.txt"),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        print(paste0("Numeric Genotypes saved at:",  home.dir))
      }
      if(ntraits > 1){
        if(set.cor & !all(ntraits==dim(correlation))){stop("ntraits and Correlation matrix do not match!", call. = F); sink()}
        
        if(ntraits != (length(h2_MT)+1) ){
          if(length(h2_MT)>(ntraits-1)){
            cat(paste("Length of h2_MT > number of correlated traits! using the first",ntraits-1,"(",paste(h2_MT[1:(ntraits-1)], collapse = ", "),")","values\n"))
            h2_MT <- h2_MT[1:(ntraits-1)]
          } else{
            cat(paste("Heritability not assigned for all correlated traits! \nSetting all traits with the same heritability (",h2[1],")\n"))
            h2_MT <- rep(h2[1], (ntraits-1))
          }
        }
        
        if(model == "partially" & !is.null(specific.QTN.number) & ntraits != length(specific.QTN.number) ){
          if(length(specific.QTN.number)>ntraits){
            cat(paste("Length of additive.effect > ntraits! using the first",ntraits,"(",paste(specific.QTN.number[1:ntraits], collapse = ", "),")","values\n"))
            specific.QTN.number <- specific.QTN.number[1:ntraits]
          } else{
            stop("Please set up [additive] specific.QTN.number", call. = F)
            sink()
          }
        }
        
        if(model == "partially" & !is.null(specific.E.QTN.number) &  ntraits != length(specific.E.QTN.number)){
          if(length(specific.E.QTN.number)>ntraits){
            cat(paste("Length of specific.E.QTN.number > ntraits! using the first",ntraits,"(",paste(specific.E.QTN.number[1:ntraits], collapse = ", "),")","values\n"))
            specific.E.QTN.number <- specific.E.QTN.number[1:ntraits]
          } else{
            stop("Please set up specific.E.QTN.number", call. = F)
            sink()
          }
        }
        
        if(ntraits != length(additive.effect)){
          if(length(additive.effect)>ntraits){
            cat(paste("Length of additive.effect > ntraits! using the first",ntraits,"(",paste(additive.effect[1:ntraits], collapse = ", "),")","values\n"))
            additive.effect <- additive.effect[1:ntraits]
          } else{
            cat(paste("Length additive.effect < ntraits! \nSetting all traits with the additive.effect = ",additive.effect[1],"\n"))
            additive.effect <- rep(additive.effect[1], ntraits)
          }
        }
        
        if(!is.null(epistatic.effect) & ntraits != length(epistatic.effect)){
          if(length(epistatic.effect)>ntraits){
            cat(paste("Length of epistatic.effect > ntraits! using the first",ntraits,"(",paste(epistatic.effect[1:ntraits], collapse = ", "),")","values\n"))
            epistatic.effect <- epistatic.effect[1:ntraits]
          } else{
            cat(paste("Length epistatic.effect < ntraits! \nSetting all traits with the epistatic.effect = ",epistatic.effect[1],"\n"))
            epistatic.effect <- rep(epistatic.effect[1], ntraits)
          }
        }
        
        if( ntraits != length(big.additive.QTN.effect)){
          if(length(big.additive.QTN.effect)>ntraits){
            cat(paste("Length of big.additive.QTN.effect > ntraits! using the first",ntraits,"(",paste(big.additive.QTN.effect[1:ntraits], collapse = ", "),")","values\n"))
            big.additive.QTN.effect <- big.additive.QTN.effect[1:ntraits]
          } else{
            cat(paste("Length big.additive.QTN.effect < ntraits! \nSetting all traits with the big.additive.QTN.effect = ",big.additive.QTN.effect[1],"\n"))
            big.additive.QTN.effect <- rep(big.additive.QTN.effect[1], ntraits)
          }
        }
      }
      
      if(ntraits == 1 | !any(model!="pleiotropic") ){
        QTN <- QTN_pleiotropic(
          genotypes = genotypes,
          seed = seed,
          Additive.QTN.number = Additive.QTN.number,
          Epistatic.QTN.number = Epistatic.QTN.number
        )
      } 
      
      if(ntraits > 1 & !any(model!="partially")){
        QTN <- QTN_partially_pleiotropic(
          genotypes = genotypes,
          seed = seed,
          overlap = overlap,
          overlapE = overlapE,
          specific.QTN.number = specific.QTN.number,
          specific.E.QTN.number = specific.E.QTN.number,
          ntraits =ntraits
        ) 
      }
      
      if (ntraits == 1) {
        Genetic_value <- Base_line_single_trait(
          seed = seed,
          additive.object = QTN$additive.effect.trait.object,
          epistatic.object = QTN$epistatic.effect.trait.object,
          additive.effect = additive.effect,
          epistatic.effect = epistatic.effect,
          big.additive.QTN.effect = big.additive.QTN.effect
        )
      } else {
        Genetic_value <- Base_line_multi_traits(
          seed = seed,
          set.cor = set.cor,
          ntraits = ntraits,
          correlation = correlation,
          additive.object = QTN$additive.effect.trait.object,
          epistatic.object = QTN$epistatic.effect.trait.object,
          additive.effect = additive.effect,
          epistatic.effect = epistatic.effect,
          big.additive.QTN.effect = big.additive.QTN.effect
        )
        cat("\nDIAGNOSTICS: \n")
        if(set.cor){
          cat("Populational Correlation \n")
          colnames(correlation) <-  c("Target", c(paste0("Trait_", 2:ntraits)))
          rownames(correlation) <- c("Target", c(paste0("Trait_", 2:ntraits)))
          print(correlation)
        }
        cat("\nSample Correlation \n")
        sample.cor <- Genetic_value$sample.cor
        colnames(sample.cor) <-  c("Target", c(paste0("Trait_", 2:ntraits)))
        rownames(sample.cor) <- c("Target", c(paste0("Trait_", 2:ntraits)))
        print(sample.cor)
      }
      
      Phenotypes(
        seed = seed,
        base.line.trait = Genetic_value$base.line,
        h2 = heritabilities.vector,
        rep = replicates,
        ntraits = ntraits,
        h2_MT =  h2_MT,
        format = format
      )
      
      closeAllConnections()
      setwd(home.dir)
      file.show(paste0(output.dir,"/Log_Sim.txt"))
    }, silent = TRUE)
    
  }#end "create.simluated.data()"




# Function modified by Rafael Della Coletta to make it compatible to simulation script.

#### K-fold validation using numeric hapmap as input ----
rrblup.kfoldfoldCV.numHapMap.input <- function(Y = NULL, Geno = NULL, traitname = "Name_of_Trait",
                                               path.for.results = path.for.results,
                                               number.of.folds = 10, seed.number = 999){
  #######################################################################################
  #Uses the GAPIT and rrBLUP R packages
  #Written by Alex Lipka on February 7, 2012
  #Read in a vector of phenotypes, and genotypic data in HapMap format (note that "head = FALSE" for reading in these data)
  #Process the data
  
  CV=Y[,1:2]
  CV[,2]=1
  colnames(CV)=c("taxa","overall")
  
  #### the line below was commented by Rafael Della Coletta
  #hm=GAPIT.HapMap(G = Geno,SNP.effect="Add",SNP.impute="Major")
  #### the line below was added by Rafael Della Coletta
  hm=Geno
  
  # the two modifications above were necessary to use kfold validation using a numeric HapMap input, which is the only
  # input that can accomodate different copy number of markers (structural variants). All the rest of the script
  # reamins the same as written by Alex Lipka.
  
  
  #####################################
  #Obtain the mafs of all SNPs
  
  #Total number of lines
  ns <- nrow(hm$GD)
  
  #Sum of the allele scores for each SNP
  ss <- apply(hm$GD, 2, sum)
  
  #Combine two situations: one where the allele coded as "2" is major; one where "0" is coded as major.
  maf.matrix <- rbind((.5*ss/ns), (1-(0.5*ss/ns)))
  
  #Copy the minor allele frequencies for all SNPs
  maf <- apply(maf.matrix, 2, min)
  
  #Find out which SNPs have MAF < 0.05
  snps.below.0.05.maf <- which(maf < 0.05)
  
  # Remove these SNPs from hm$GD
  #### the line below was commented by Rafael Della Coletta
  #hm.GD.without.snps.below.0.05.maf <- hm$GD[,-snps.below.0.05.maf]
  #### the lines below were added by Rafael Della Coletta
  if (length(snps.below.0.05.maf) > 0) {
    hm.GD.without.snps.below.0.05.maf <- hm$GD[,-snps.below.0.05.maf]
  } else {
    hm.GD.without.snps.below.0.05.maf <- hm$GD
  }
  
  # this code was added because if there is no SNP with maf < 0.05, R will create a vector with
  # "integer(0)", which will give an error when running "hm$GD[,-snps.below.0.05.maf]". I added
  # an "if else" statement to avoid this behavior, since if there is no SNP with maf < 0.05, there
  # is nothing to exclude from "hm$GD".
  
  ###############################
  
  GK <- cbind(hm$GT, hm.GD.without.snps.below.0.05.maf)
  
  qc=GAPIT.QC(Y = Y, GT = hm$GT, CV = CV, GK = GK)
  
  y <- as.matrix(qc$Y[-1])
  
  G <- as.numeric(qc$GK[,-1])
  
  G <- matrix(G, nrow(y), ncol(qc$GK[,-1]))
  
  G <- G - 1
  
  cv <- (as.matrix(qc$CV[,-1]))
  
  taxa.names <- qc$CV[,1]
  
  #Calculate the kinship matrix in rrBLUP
  A1 <- A.mat(G,shrink=TRUE)
  
  
  ############################################################################################
  #Run k-fold cross validation
  
  sample.size <- nrow(y)
  
  #Randomly sort the number of lines, and subdivide them into ten subgroups
  set.seed(seed.number)
  sequence.sample <- rep(1:sample.size)
  random.sample <- sample(1:sample.size, replace = FALSE)
  increment <- ceiling(length(random.sample)/number.of.folds) 
  r.gy <- NULL 
  
  #have a "for" loop, start it at 0, and end it at 9
  #I am setting up "k" to denote the nubmer of folds - 1. This is done
  # so that the for loop will work correctly.
  k <- number.of.folds - 1
  
  vector.of.observed.values <- NULL
  vector.of.predicted.values <- NULL
  vector.of.taxa.names <- NULL
  for(i in 0:k){
    pred <- random.sample[((increment*i)+1):min(((increment*i)+increment) , sample.size)]
    train <- random.sample[-(((increment*i)+1):min(((increment*i)+increment) , sample.size))] 
    
    yNA <- y
    yNA[pred] <- NA
    
    data1 <- data.frame(y=yNA,gid=1:length(y), cv = cv)
    the.cv.names <- NULL
    for(j in 1:ncol(cv)) the.cv.names <- c(the.cv.names, paste("CV_",j,sep = ""))
    
    colnames(data1) <- c("y","gid", the.cv.names)
    
    rownames(A1) <- 1:nrow(A1)
    ans1 <- kin.blup(data1,K=A1,geno="gid",pheno="y", covariate = the.cv.names)
    r.gy <- c(r.gy, cor(ans1$g[pred], y[pred]) )
    
    
    #ans <- kinship.BLUP(y=y[train],G.train=G[train,],G.pred=G[pred,],K.method="RR")
    
    #r.gy <- c(r.gy, cor(ans$g.pred,y[pred]))
    
    vector.of.observed.values <- c(vector.of.observed.values, y[pred])
    vector.of.predicted.values <- c(vector.of.predicted.values, ans1$g[pred])
    vector.of.taxa.names <- c(vector.of.taxa.names, as.character(taxa.names[pred]))
  }
  
  
  #Fit a regression model, where:
  # Observed.Phenotype_i = beta_0 + beta_1*Predicrted.Phenotype_i + epsilon_i
  SLR.model <- lm(vector.of.observed.values ~ vector.of.predicted.values)
  
  pdf(paste0(path.for.results,number.of.folds,"-fold_CV_Results_",traitname,"Obs.vs.Pred.pdf"))
  plot(vector.of.observed.values ~ vector.of.predicted.values, col = "blue", xlab = "Predicted values", ylab = "Observed Values")
  abline(SLR.model)
  legend("topleft", paste("Intercept = ", round(coef(SLR.model)["(Intercept)"],2), ", Slope = ", round(coef(SLR.model)["vector.of.predicted.values"],2), sep = ""))
  dev.off()
  # Once the loop is over, output the values of the correlation coefficients, as well as their means and standard deviations  
  
  #Create a vector of column names for the k-fold cross validation output
  colname.r.gy <- NULL
  for(i in 1:number.of.folds) colname.r.gy <- c(colname.r.gy, paste("r_CV",i,sep = ""))
  
  r.gy <- c(r.gy, mean(r.gy), sd(r.gy))
  
  r.gy.output <- t(as.matrix(r.gy))
  
  
  colnames(r.gy.output) <- c(colname.r.gy, "r_avg", "r_std")
  
  #
  
  write.table(r.gy.output, paste0(path.for.results,number.of.folds,"-fold_CV_Results_",traitname,".txt"), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  
  #Print out the observed and predicted values
  the.observed.predicted.and.taxa.names <- cbind(as.character(vector.of.taxa.names), vector.of.observed.values, vector.of.predicted.values)
  
  colnames(the.observed.predicted.and.taxa.names) <- c("SampleID", paste("Observed_", traitname, sep = ""), paste("Predicted_", traitname, sep = ""))
  
  write.table(the.observed.predicted.and.taxa.names, 
              paste0(path.for.results,"Obs.and.Predicted.Trait.values",traitname,".txt"), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
} #End rrblup.tenfoldCV