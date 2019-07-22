# Testing simulations

Since my visit to Alex Lipka's lab, I've been working to implement the trait simulation code written by Samuel Fernandes (Alex's postdoc) and the k-fold validation of genomic prediction accuracies written by Alex in one single script according to our USDA dataset and our goals for the manuscript(s). Our first goal is to test the effects of structural variants in genomic prediction on the USDA RIL population only (i.e., simulating a breeding program for inbred lines). This will be our first manuscript. We will also add the effects of different environmnents later. Only in a follow-up manuscript we will simulate hybrid genotypes and test the effects of SVs in genomic prediction in hybrids.

A more detailed overview of the goals for this project can be seen in the power point presentation `notes/simulation_meeting_2019-06-19.pptx`.


## Toy dataset

Samuel created the dataset `Structural_Variation_demo.txt` with the goal of simulate copy number variants (CNVs). He selected 2,000 random SNPs from the `SNP55K_maize282_AGPv2_20100513_1.hmp.txt` dataset (55K SNPs from 282 diverse maize lines) and converted the 0,1,2 scale (i.e., homozygous recessive, heterozygote, homozygous dominant) to 1,3,5 scale (i.e., 1, 3 or 5 copies of a marker).

```
SNP55K_maize282_AGPv2_20100513_1.hmp_NUM-test.txt
SNP55K_maize282_AGPv2_20100513_NUM.txt
SNP55K_maize282_AGPv2_20100513_SNP-SV-merged_test.txt
```

## RILs only, one location

Right now, I'm focusing on the first manuscript. My first objective is to adjust the script to simulate 27 traits in one location (replicated 50 times), run predictions models, and the use k-fold validation to determine prediction accuracies.







<mark>TO DO:</mark>
1. Add all auxiliar functions in the script `simulation_auxiliar-functions.R` (except for the GAPIT script, which will also copied to the scripts folder).
2. Make sure code follows Google style guide for R.
3. Move `results_validations.R` to this simulation folder, change the name of the script to `...`, and then run it to generate plots for different number of markers.
4. Move `trait-sim_manuscript.R` to this simulation folder, change the name of the script to `...`, and start version control.
