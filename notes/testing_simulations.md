# Testing simulations

Since my visit to Alex Lipka's lab, I've been working to implement the trait simulation code written by Samuel Fernandes (Alex's postdoc) and the k-fold validation of genomic prediction accuracies written by Alex in one single script according to our USDA dataset and our goals for the manuscript(s). Our first goal is to test the effects of structural variants in genomic prediction on the USDA RIL population only (i.e., simulating a breeding program for inbred lines). This will be our first manuscript. We will also add the effects of different environmnents later. Only in a follow-up manuscript we will simulate hybrid genotypes and test the effects of SVs in genomic prediction in hybrids.

A more detailed overview of the goals for this project can be seen in the power point presentation `notes/simulation_meeting_2019-06-19.pptx`.

The script `tests/script/trait-sim_manuscript.R` was written to first simulates trait using the `simulate_trait()` function, and then use k-fold validation to get the genomic prediction accuracies with the `kfold_validation_on_sim_traits()` function. The script `tests/script/results_validations.R` was written to plot results of accuracies for different traits. The scripts are heavily commented, so I will just highlight the main options/results/issues.




## `simulate_trait()` function

I simulated traits using the function `simulate_trait()` that I wrote. I basically added more functionalities to the function `create_simulated_data()` written by Samuel to make sure the input data with or without structural variants are read correcly, and ensure that all options for simulating a trait defined in our meetings (QTN number, heritability, etc.) are also taken into consideration. Basically, with this function, you can:

* Use a genotypic dataset containing only SNPs, only SVs, or both SNPs and SVs (for the latter, you should also provide a list of SVs IDs to distinguish between the two types of variants).
* Choose the source of variation: if only SNPs are responsible for the trait, only SVs, or both SNPs and SVs contribute.
* Select the number of QTNs controlling the trait, the effect of the largest and smallest QTNs, the heritability and the number of replicates per simulated trait.

> There is currently a bug in the code. If I run my `simulate_trait()` function right away, it generates an incomplete output, but it doesn't show the simulated traits. However, if I run the `create_simulated_data()` function (from Samuel) outside `simulate_trait()` and then run `simulate_trait()`, my function works just fine (i.e., the simulated traits are generated). I think there is a problem in setting up a connection to write the phenotypic files, but will need to investigate this more.




## `kfold_validation_on_sim_traits()` function

I used k-fold validation to get genomic prediction accuracies for the simulated traits described above. This functions builds up from the `rrblup.kfoldfoldCV()` function written by Alex Lipka, and accounts for SNP and/or SV data. Briefly, with this function, you can:

* Use a genotypic dataset containing only SNPs, only SVs, or both SNPs and SVs (for the latter, you should also provide a list or file with SVs IDs to distinguish between the two types of variants).
* Choose what types of markers you want to use for genomic prediction: all SNPs (1), SNPs in LD with SV (2), SNPs with varying LD with SV (3), only SVs (4), or both SNPs and SVs information (5).
* Define how many folds to perform validation (default is 5-fold validation).
* Decide if you want to use all SNPs in the genotypic dataset, or randomly sample a certain number of SNPs (e.g., 1000).

> I still have to implement the calculation of LD in this script. I've been used fixed distances to determine SNP in LD with SV (+/- 10kb) or varying LD (+/- 50 kb and excluding the first +/- 10kb).

It should be noted that I did small modifications to the original `rrblup.kfoldfoldCV()` function from Alex to adapt it to my needs. Thus, I'm using my version of this function, which is called `rrblup.kfoldfoldCV.numHapMap.input()`. Basically, the main difference is that my version requires a numeric hapmap format as input instead of the standard hapmap.




## Toy dataset

Samuel created the dataset `tests/data/Structural_Variation_demo.txt` with the goal of simulating copy number variants (CNVs). He selected 2,000 random SNPs from the `tests/data/SNP55K_maize282_AGPv2_20100513_1.hmp.txt` dataset (55K SNPs from 282 diverse maize lines) and converted the 0,1,2 scale (i.e., homozygous recessive, heterozygote, homozygous dominant) to 1,3,5 scale (i.e., 1, 3 or 5 copies of a marker).

I generated other related datasets using the function `read_input_files()` from script `tests/script/trait-sim_manuscript.R`. The file `tests/data/SNP55K_maize282_AGPv2_20100513_NUM.txt` is just the numeric version of the standard hapmap file, and `tests/data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt` is also a numeric version of the hapmap but with the dummy SVs included. These files were created just to speed up testing (use them instead of doing the conversion every time).

After running the `simulate_trait()` function with the toy dataset, the simulated traits can be found in `tests/analysis/test_toy` folder. I generated traits using the following combination of arguments:

* SNPs, SVs, or both were used as source of variation for the traits;
* QTN numbers were 3, 25 or 75;
* the large QTN effect was 0.7 and the small was 0.3;
* heritability was either 0.2, 0.5 or 0.9;
* number of replicates was set to 50;
* seed number was 2019.

After running the `kfold_validation_on_sim_traits()` function, the genomic prediction results and validations can be found in sub-folders for each simulated trait in `tests/analysis/test_toy`. The options I used for k-fold validation are a combination of the following:

* 5-fold validation without subsetting the marker data, using only 1000 SNPs for genomic prediction, or using only 50 SNPs;
* all marker types were used when not subsetting marker data (1 = all SNPs, 2 = SNPs in LD with SV, 3 = SNPs with varying LD with SV, 4 = only SVs, or 5 = both SNPs and SVs), and all but `marker_data_type = 4` were used when subsetting SNPs to 1000 or 50;
* testing arguments was set to `TRUE` (i.e., only 3 replicates of the simulated trait was analyzed).

> Check script to know which dataset file was used for each type of marker used.





## USDA dataset

The original USDA dataset information and the transformations I made are described in more details in `notes/gs_simulation.md`. Briefly, it's a 22k SNP chip dataset of 525 RILs generated from 7 inbred parents, which I had to transform from `.csv` to `.hmp.txt` format and then reduce the numbers of RILs (from 525 to 325) to match the ones we planted last year to generate hybrids. Therefore, the datasets that will be testing here are `usda_22kSNPs_7parents.hmp.txt` and `usda_22kSNPs_325rils.hmp.txt`.

After running the `simulate_trait()` function with the toy dataset, the simulated traits can be found in `tests/analysis/test_usda-parents` or `tests/analysis/test_usda-rils` folders for traits of parents and RILs, respectively.I generated traits using the following combination of arguments:

* SNPs, SVs, or both as source of variation for the traits;
* QTN numbers were 3, 25 or 75;
* the large QTN effect was 0.7 and the small was 0.3;
* heritability was either 0.2, 0.5 or 0.9;
* number of replicates was set to 50;
* seed number was 2019.

After running the `kfold_validation_on_sim_traits()` function, the genomic prediction results and validations can be found in sub-folders for each simulated trait in `tests/analysis/test_usda-parents` or `tests/analysis/test_usda-rils`. The options I used for k-fold validation are a combination of the following:

* 5-fold validation without subsetting the marker data, using only 1000 SNPs for genomic prediction, or using only 50 SNPs;
* SNPs were the only marker type used (`marker_data_type = 1`), because I still don't have the SV information for the USDA dataset;
* testing arguments was set to `TRUE` (i.e., only 3 replicates of the simulated trait was analyzed).

> Once I get the SV data, I will be able to run the other options (2 to 5), similarly to the toy dataset.




## Plot partial results

I wrote the script `results_validations.R` to visualize prediction accuracies for different simulated traits. Plots can be found at `tests/analysis/test_toy` or `tests/analysis/test_usda-rils`. Note that I didn't plot results from parents of the USDA dataset since there were only 7 parents and predictions are not accurates at all (this dataset was used basically to test if the scripts were working, since it had very small data).

Overall, for both toy and USDA datasets, there is very little difference in prediction accuracies among simulated traits controlled by 3, 25 or 75 QTNs. In addition, increasing trait heritability consistently increased prediction accuracy. Although increasing the number of markers used (50, 1000 or all) also increased prediction accuracy, it's interesting to see that there is not much difference between using 1000 or all SNPs, especially in the USDA dataset.

Looking at the plots from the toy dataset, you can see that using only SVs to predict traits has poor performance compared to other types of markers, even when the source of trait variation is due to SVs. This might have to do with how these SVs were simulated, and it will be interesting to see what happens when we add the real SVs from USDA data into the simulation.

It will be interesting to test different parameters, such as the size of the large and small effects.




# Preliminary SV calls dataset

On August 9, 2019, Patrick sent me 5 `.vcf` files containing structural variation calls from the software Lumpy. Each file is a SV call of 100 lines against one of the following reference genomes: B73, Mo17, PH207, PHB47, or W22. **This results are preliminary** and there are a lot of false-positives in there. However, they will be useful for me to write scripts to select the 8 lines I need, project the SVs from parents to RILs, and incorporate these "real" SV data in my simulation scripts (instead of the fake toy dataset).

The files are located at `tests/data/`, and I will use only the SVs called against the B73 reference genome `B73v4_2019-08-09.ls.RT.vcf` for testing.



# Simulations with SVs included

`tests/scripts/trait-sim_manuscript_SVs.R`


<mark>TO DO:</mark>
* Test different leveles of QTN effects. Look at Sam's code and find out what really is going on with the code: only large effect go for geometric series? Or does the geometric series start from the smallest effect QTLs.

* Find out what happens with missing data when doing numericalization! --> TRANSFORM TO MAJOR ALLELE BASED ON WHICH ALLELES FLANK IT!





<mark>TO DO!!!!
* CORRECT PLINK LD CALCULATION IN SIMULATION SCRIPT!!!
