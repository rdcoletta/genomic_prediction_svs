# Testing simulations

Since my visit to Alex Lipka's lab, I've been working to implement the trait simulation code written by Samuel Fernandes (Alex's postdoc) and the k-fold validation of genomic prediction accuracies written by Alex in one single script according to our USDA dataset and our goals for the manuscript(s). Our first goal is to test the effects of structural variants in genomic prediction on the USDA RIL population only (i.e., simulating a breeding program for inbred lines). This will be our first manuscript. We will also add the effects of different environmnents later. Only in a follow-up manuscript we will simulate hybrid genotypes and test the effects of SVs in genomic prediction in hybrids.

A more detailed overview of the goals for this project can be seen in the power point presentation `notes/simulation_meeting_2019-06-19.pptx`.

The script `tests/script/trait-sim_manuscript.R` was written to first simulates trait using the `simulate_trait()` function, and then use k-fold validation to get the genomic prediction accuracies with the `kfold_validation_on_sim_traits()` function. The script `tests/script/results_validations.R` was written to plot results of accuracies for different traits. The scripts are heavily commented, so I will just highlight the main options/results/issues.

## `simulate_trait()` function

I simulated traits using the function `simulate_trait()` that I wrote. I basically added more functionalities to the function `create_simulated_data()` written by Samuel to make sure the input data with or without structural variants are read correcly, and ensure that all options for simulating a trait defined in our meetings (QTN number, heritability, etc.) are also taken into consideration. Basically, with this function, you can:

* Use a genotypic dataset containing only SNPs, only SVs, or both SNPs and SVs (for the latter, you should also provide a list of SVs IDs to distinguish between the two types of variants).
* Choose the source of variation: if only SNPs are responsible for the trait, only SVs, or both SNPs and SVs contribute.
* Select the number of QTNs controlling the trait, the effect of the largest and smallest QTNs, the heritability and the number of replicates per simulated trait.

> **Note:** there is currently a bug in the code. If I run my `simulate_trait()` function right away, it generates an incomplete output, but it doesn't show the simulated traits. However, if I run the `create_simulated_data()` function (from Samuel) outside `simulate_trait()` and then run `simulate_trait()`, my function works just fine (i.e., the simulated traits are generated). I think there is a problem in setting up a connection to write the phenotypic files, but will need to investigate this more.


## `kfold_validation_on_sim_traits()` function

I used k-fold validation to get genomic prediction accuracies for the simulated traits described above. This functions builds up from the `rrblup.kfoldfoldCV()` function written by Alex Lipka, and accounts for SNP and/or SV data. Briefly, with this function, you can:

* Use a genotypic dataset containing only SNPs, only SVs, or both SNPs and SVs (for the latter, you should also provide a list or file with SVs IDs to distinguish between the two types of variants).
* Choose what types of markers you want to use for genomic prediction: all SNPs (1), SNPs in LD with SV (2), SNPs with varying LD with SV (3), only SVs (4), or both SNPs and SVs information (5).
* Define how many folds to perform validation (default is 5-fold validation).
* Decide if you want to use all markers in the genotypic dataset, or randomly sample a certain number of SNPs (e.g., 1000).


## Toy dataset

Samuel created the dataset `tests/data/Structural_Variation_demo.txt` with the goal of simulating copy number variants (CNVs). He selected 2,000 random SNPs from the `tests/data/SNP55K_maize282_AGPv2_20100513_1.hmp.txt` dataset (55K SNPs from 282 diverse maize lines) and converted the 0,1,2 scale (i.e., homozygous recessive, heterozygote, homozygous dominant) to 1,3,5 scale (i.e., 1, 3 or 5 copies of a marker).

I generated other related datasets using the function `read_input_files()` from script `tests/script/trait-sim_manuscript.R`. The file `tests/data/SNP55K_maize282_AGPv2_20100513_NUM.txt` is just the numeric version of the standard hapmap file, and `tests/data/SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt` is also a numeric version of the hapmap but with the dummy SVs included. These files were created just to speed up testing (use them instead of doing the conversion every time).

After running the `simulate_trait()` function with the toy dataset, the simulated traits can be found in `tests/analysis/test_toy` folder. I generated traits using SNPs, SVs, or both as source of variation for the traits; QTN numbers were 3, 25 or 75; the large QTN effect was 0.7 and the small was 0.3; heritability was either 0.2, 0.5 or 0.9; number of replicates was set to 50; and seed number was 2019.

After running the `kfold_validation_on_sim_traits()` function with the toy dataset, five folders containing the prediction accuracies from the k-fold validation was generated (`k-fold_validation_all-SNPs`, `k-fold_validation_only-SVs`, `k-fold_validation_SNPs-and-SVs`, `k-fold_validation_SNPs-LD`, and `k-fold_validation_snps-varying-LD`) inside each of the folders containing the simulated traits above. I used 5-fold validation to get genomic prediction accuracies, using all SNPs (1), SNPs in LD with SV (2), SNPs with varying LD with SV (3), only SVs (4), or both SNPs and SVs information (5). Check script for the other parameters, such as dataset file used, since it varies depending on type of marker used.




## USDA dataset

<mark>TO DO:</mark> describe dataset and where i got it from...

After running the `simulate_trait()` function with the toy dataset, the simulated traits can be found in `tests/analysis/test_usda-parents` or `tests/analysis/test_usda-rils` folders for traits of parents and RILs, respectively. I generated traits using SNPs, SVs, or both as source of variation for the traits; QTN numbers were 3, 25 or 75; the large QTN effect was 0.7 and the small was 0.3; heritability was either 0.2, 0.5 or 0.9; number of replicates was set to 50; and seed number was 2019.

<mark>TO DO:</mark> run k-fold cross validation on USDA dataset







<mark>TO DO:</mark>
1. Add all auxiliar functions in the script `simulation_auxiliar-functions.R` (except for the GAPIT script, which will also copied to the scripts folder).
2. Make sure code follows Google style guide for R.
3. Move `results_validations.R` to this simulation folder, change the name of the script to `...`, and then run it to generate plots for different number of markers.
4. Move `trait-sim_manuscript.R` to this simulation folder, change the name of the script to `...`, and start version control.
