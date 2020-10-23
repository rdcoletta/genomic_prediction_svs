# Structural variation and genomic prediction: a simulation approach

by Rafael Della Coletta, Alex Lipka, Martin Bohn, and Candice Hirsch (June 2019 - February 2020).

> The objective of this project is to simulate traits in one or more environments and analyze the effects of structural variants in genome prediction accuracy. The overall workflow for the simulations involves simulating traits, running genomic prediction models and validating results.

<!-- TOC START min:1 max:3 link:true asterisk:false update:true -->
- [Structural variation and genomic prediction: a simulation approach](#structural-variation-and-genomic-prediction-a-simulation-approach)
  - [Project folder](#project-folder)
  - [SNP chip dataset](#snp-chip-dataset)
    - [HapMap format](#hapmap-format)
    - [Data QC](#data-qc)
    - [Remove low quality SNPs](#remove-low-quality-snps)
    - [Check reference genome version of SNP chip](#check-reference-genome-version-of-snp-chip)
    - [Convert SNP coordinates from B73v2 to B73v4](#convert-snp-coordinates-from-b73v2-to-b73v4)
    - [Divide dataset by cross](#divide-dataset-by-cross)
  - [SV dataset](#sv-dataset)
    - [Hapmap format](#hapmap-format-1)
    - [Divide dataset by cross](#divide-dataset-by-cross-1)
  - [Remove chip SNPs inside deletions](#remove-chip-snps-inside-deletions)
  - [Correct SNP chip miscalls with sliding window approach](#correct-snp-chip-miscalls-with-sliding-window-approach)
  - [Resequencing SNPs](#resequencing-snps)
    - [Divide dataset by cross](#divide-dataset-by-cross-2)
    - [Remove resequecing SNPs inside deletions and keep only polymorphic SNPs](#remove-resequecing-snps-inside-deletions-and-keep-only-polymorphic-snps)
  - [Merge SNP and SV data](#merge-snp-and-sv-data)
  - [Project SVs and resequencing SNPs from parents to RILs](#project-svs-and-resequencing-snps-from-parents-to-rils)
    - [Sliding window approach](#sliding-window-approach)
    - [Add monomorphic resequencing SNPs back](#add-monomorphic-resequencing-snps-back)
    - [Add back SNPs and SVs not present in a certain cross](#add-back-snps-and-svs-not-present-in-a-certain-cross)
    - [Merge projected markers from different families](#merge-projected-markers-from-different-families)
  - [LD structure between SNPs and SVs](#ld-structure-between-snps-and-svs)
    - [Distribution](#distribution)
    - [Decay](#decay)
    - [Closest SNPs to SVs](#closest-snps-to-svs)
    - [LD by common parent](#ld-by-common-parent)
    - [LD by family](#ld-by-family)
    - [PCA](#pca)
  - [Trait Simulation](#trait-simulation)
<!-- TOC END -->




## Project folder

All data, scripts, analyses, notes and other things related to this project is located on my MSI account:

```bash
cd /home/hirschc1/della028/projects/genomic_prediction/simulation/
mkdir {analysis,data,scripts}
```




## SNP chip dataset

The USDA project contains 525 RILs generated from 7 inbred parents (B73, PHJ40, PHG39, PHG47, PH207, PHG35, LH82), which are all ex-PVPs. In addition, 400 F1 hybrids were generated from a partial diallel cross of those RILs.

There is genotypic information (22k SNP chip) for all inbred parents and their RILs. The dataset was provided by Martin Bohn from a shared folder on DropBox ("03_DispensibleGenome/EMAMP_SNP_2018"). The `.zip` files were dowloaded and decompressed on May 28, 2019. They were then transferred to `data/SNP_chip/` via FileZilla.

```
Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.csv
Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part1_122318.csv
Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part2_122318.csv
```
> There are more genotypes (parents and RILs) in those files besides the ones used in the USDA project.

After some preliminary analysis, we found some errors in the names of RILs used for genotyping:
* 19 RILs have duplicated genotypic information but actually had different source IDs. Martin sent me the file `data/snp-chip_duplicated-missing_rils_MB.xlsx` with the correct sources to use. I manually deleted the duplicated RIL with the wrong source ID in the above excel files.
* 17 RILs used to create hybrids didn't have genotypic information in the SNP chip files. In reality, they were genotyped but under different names. Martin sent me the file `data/B73_PH207_PHG39_Issue.xlsx` wtih the correct name of RILs. I manually changed the RIL names in the above excel files.
* 1 RIL (LH82\*PHG47-B-B-31-1-1-B-B) was genotyped a generation earlier (LH82\*PHG47-B-B-31-1-1-B). I manually changed the RIL name in the above excel files.

Therefore, the final SNP chip files used for analysis were:

```
Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.csv
Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part1_122318_corrected.csv
Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part2_122318_corrected.csv
```



### HapMap format

Before doing any analysis, it's important to transform the genotypic data that comes from the SNP chip into the HapMap format (see <https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load> for more info about the format). Here's the basic format:

| rs#   | alleles | chrom | pos | strand | assembly  | center   | protLSID | assayLSID | QCcode      | Line 1 | Line 2 | ... |
| ----- | ------- | ----- | --- | ------ | --------- | -------- | -------- | --------- | ----------- | ------ | ------ | --- |
| SNP 1 | A/C     | 1     | 20  | +      | refgen_v4 | MaizeGDB | NA       | NA        | maize_panel | CC     | AA     | ... |

The first 11 columns are required, but not all are required to have information. Usually, only the fields `rs#`, `alleles`, `chrom`, `pos` and the lines' genotypes will be used in downstream analysis (the other fields can be NA).

I used the `scripts/usda_geno2hmp.R` script to transform all parental and RIL data to hapmap format:

```bash
Rscript scripts/usda_geno2hmp.R data/SNP_chip/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.csv \
                                data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part1_122318_corrected.csv \
                                data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part2_122318_corrected.csv \
                                data
```

This script generated these files on `data/` folder:
* Genotypic data on hapmap format:
  ```
  Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318_corrected.hmp.txt
  Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318_corrected.hmp.txt
  ```
* Lists relating genotype name with the ID of the source used to generate that RIL:
  ```
  id_table_Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318_corrected.txt
  id_table_Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318_corrected.txt
  ```

> RILs have their genotype ID instead of name because some RILs had multiple IDs. That's why I generated a list relating names to IDs. Also, I converted `--` to `NN` to represent missing data for compatibility with TASSEL 5 software.

Since there are genotypic data for more parents and RILs than used in the USDA project, I need to keep only the genotypes described in the USDA project that will have phenotypic data collected. To do this, Candy sent the file `data/2018_field_planning.xlsx` that contains all the RILs crossesd to generate the 400 hybrids. They can be found under **USDA crosses** on `X10_Nursery_Book` and `X9_Nursery_Book` sheets. I manually copied all this information, removed duplicates (in Excel) and saved as file `data/usda_RILs_2018.txt`.

> Before creating this file, I noticed that there were 328 unique RILs and one hybrid (`LH82*PHG47`). That was a mistake, and I changed the name of the  of this hybrid to `LH82*PHG47-B-B-4-1-1-B-B`, `LH82*PHG47-B-B-5-1-1-B-B`, `LH82*PHG47-B-B-6-1-1-B-B`, `LH82*PHG47-B-B-7-1-1-B-B`, and `LH82*PHG47-B-B-8-1-1-B-B`. Therefore, the total number of RILs are 333.

The `scripts/remove_extra_geno-data.R` script removes extra genotypic data that will not be used in the USDA project from hapmap files:

```bash
Rscript scripts/remove_extra_geno-data.R data/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.hmp.txt \
                                         data/usda_22kSNPs_parents.hmp.txt \
                                         data/Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318_corrected.hmp.txt \
                                         data/usda_22kSNPs_rils.hmp.txt \
                                         data/usda_RILs_2018.txt
```

I will use the software [Tassel 5](https://www.maizegenetics.net/tassel) for some basic QC, because this software can generate summaries (including allele frequency) pretty quickly.

Tassel requires that the genotypic data is sorted by position and chromosome number. The genotypic data `data/usda_22kSNPs_rils.hmp.txt`has SNPs from chromosome 10 comes after chromosome 1, so it's not sorted correctly. Thus, I used `SortGenotypeFilePlugin` from Tassel to quickly sort the genotypic data and create another hapmap file called `data/usda_22kSNPs_rils.sorted.hmp.txt`. Then I ran Tassel again to transform the sorted data into diploid hapmap format. This file will be used for the rest of QC analysis. I also did the same procedure with the parental data.

```bash
# sort ril data
run_pipeline.pl -Xmx2g -SortGenotypeFilePlugin \
                       -inputFile data/usda_22kSNPs_rils.hmp.txt \
                       -outputFile data/usda_22kSNPs_rils.sorted.hmp.txt \
                       -fileType Hapmap
run_pipeline.pl -Xmx2g -importGuess data/usda_22kSNPs_rils.sorted.hmp.txt \
                       -export data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                       -exportType HapmapDiploid

# sort parents
run_pipeline.pl -Xmx2g -SortGenotypeFilePlugin \
                       -inputFile data/usda_22kSNPs_parents.hmp.txt \
                       -outputFile data/usda_22kSNPs_parents.sorted.hmp.txt \
                       -fileType Hapmap
run_pipeline.pl -Xmx2g -importGuess data/usda_22kSNPs_parents.sorted.hmp.txt \
                       -export data/usda_22kSNPs_parents.sorted.diploid.hmp.txt \
                       -exportType HapmapDiploid
```

These commands generated the following files on `data/` folder, and will be used in further downstream analysis:

```bash
# parental data
usda_22kSNPs_parents.sorted.diploid.hmp.txt
# RIL data (315 out of 333 RILs were genotyped)
usda_22kSNPs_rils.sorted.diploid.hmp.txt
```


### Data QC

The QC of parental and RIL data will be done on a family basis. Thus, I have to find out which RILs belong to the same biparental cross. I wrote the `scripts/create_list_biparental-crosses.py` script to generate a table with a biparental cross in a column and all the RIL IDs in another, and ran this in the command line to create the table:

```bash
python scripts/create_list_biparental-crosses.py data/id_table_Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318_corrected.txt \
                                                 data/usda_biparental-crosses.txt
```

I wrote `scripts/markers_summary.R` to summarize the number of missing and heterozygous markers for the parents and RILs and save them at the `analysis/qc` folder.

For parents, it uses information from `data/usda_22kSNPs_parents.sorted.diploid.hmp.txt`, and generates the table `summary_markers_7parents.txt` and the boxplots `summary_markers_7parents_missing.png` and `summary_markers_7parents_het.png`.

For RILs, the script takes information from `data/usda_22kSNPs_rils.sorted.diploid.hmp.txt` and `data/usda_biparental-crosses.txt`, and generates the table `summary_markers_rils.txt` and the boxplots `summary_markers_rils_missing-per-cross.png` and `summary_markers_rils_hets-per-cross.png` at the `analysis/qc` folder. In addition, the script also checks the number of missing and heterozygous genotypes in RILs that have the same parent to help identify any problems with a parent during generation of RILs (`summary_markers_325rils_missing-per-parent.png` and `summary_markers_325rils_het-per-parent.png`).

```bash
# create directory to store qc results
mkdir -p analysis/qc
# run qc
Rscript scripts/markers_summary.R data/usda_22kSNPs_parents.sorted.diploid.hmp.txt \
                                  analysis/qc/summary_markers_parents.txt \
                                  data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                                  analysis/qc/summary_markers_rils.txt \
                                  data/usda_biparental-crosses.txt
```

#### Parental data

It's very important to make sure the parental lines also have good quality data, especially because parental genotypic data will be used for projection of SVs into RILs. By looking at the bar plots of parents genrated earlier, it's very clear that **PHG35 has higher missing and het markers than other parents** (~4% compared to ~0.3% of other parents).

Another way to inspect parents is by plotting the karyotypes with homo, het, and missing markers to see if there are some regions of the genome with higher proportion of hets and missing data.  To plot karyotypes I need information from chromosome length and centromere positions. The chromosome info was manually extracted from `GCA_000005005.6_B73_RefGen_v4_assembly_stats.txt` of B73_RefGen_v4 assembly at [MaizeGDB](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4) (accessed 07/19/2019) and saved in a tab-delimited text file (`data/B73_RefGen_V4_chrm_info.txt`).

Centromeres positions were extracted from Table S2 of [Schneider et al 2016](https://doi.org/10.1073/pnas.1522008113) and stored at `data/centromeres_Schneider-2016-pnas_v2-v3.bed`. However, coordinates of centromeres from chromosomes 1 to 9 are from AGPv2 assembly, and chromosome 10 is from AGPv3. Thus, I used the [Assembly Converter](http://ensembl.gramene.org/Zea_mays/Tools/AssemblyConverter?db=core) from Gramene  with the data above to convert coordinates from v2 to v4 (`data/centromeres_v2-to-v4.bed`) and v3 to v4 (`data/centromeres_v3-to-v4.bed`). Then, I wrote the R script `script/get_centromeres_v4_coord.R` and filtered the previous files to get the centromeres coordinates in the v4 assembly, which can be seen in `data/centromeres_Schneider-2016-pnas_v4.bed` and in the table below:

```bash
Rscript scripts/get_centromeres_v4_coord.R data/centromeres_v2-to-v4.bed \
                                           data/centromeres_v3-to-v4.bed \
                                           data/centromeres_Schneider-2016-pnas_v4.bed
```

| chr | cent_start | cent_end  | origin_coordinates |
| --- | ---------- | --------- | ------------------ |
| 1   | 136822299  | 137672564 | RefGen_v4          |
| 2   | 95520346   | 97475434  | RefGen_v4          |
| 3   | 112459593  | 113664599 | RefGen_v4          |
| 4   | 107740595  | 110268603 | RefGen_v4          |
| 5   | 104656279  | 106318158 | RefGen_v4          |
| 6   | 51342080   | 60057564  | RefGen_v4          |
| 7   | 63706708   | 64300427  | RefGen_v4          |
| 8   | 49979865   | 51601476  | RefGen_v4          |
| 9   | 53908684   | 55031062  | RefGen_v4          |
| 10  | 51516139   | 52771682  | RefGen_v4          |


Finally, to plot the parental karyotypes, I wrote `scripts/plot_parents_karyotypes.R` and saved them at `analysis/qc/parents`.

```bash
Rscript scripts/plot_parents_karyotypes.R data/B73_RefGen_V4_chrm_info.txt \
                                          data/centromeres_Schneider-2016-pnas_v4.bed \
                                          data/usda_22kSNPs_parents.sorted.diploid.hmp.txt \
                                          analysis/qc/parents
```

Looking at the karyotypes, the distribution of both missing and het markers seem to occur more further away from the centromeres but in a bit random way. Except for PHG35, I don't see any big clusters of missing and/or het markers. Since PHG35 had much more missing and het data than other parents, it's easier to see some clusters of such markers.

After talking to Candy, it's likely that the PHG35 souce used to be genotyped had issues (either some backcross or cross-contamination). Thus, in order to correctly project the SVs from this parent to its RILs, we need to make sure the genotype from PHG35 agrees with the RILs. To do that, we will reconstruct the PHG35 genome based on the markers available from the RILs.



#### Reconstruct PHG35 genotype

Based on previous QC, the PHG35 parent has much more heterozygotes than expected for a fully inbred line. This might very likely be due to polen contamination in the seeds used to do the genotyping. Now that we have SNP chip data in the refgen v4 assembly, I can use the marker genotypes in the RIL data from all PHG35 progeny and figure out the PHG35 haplotype based on the genotype of the other parent used to develop that progeny. To do that, I wrote `scripts/reconstruct_PHG35_from_RIL_data.R`.

```bash
Rscript scripts/reconstruct_PHG35_from_RIL_data.R data/usda_biparental-crosses.txt \
                                                  data/usda_22kSNPs_parents.sorted.diploid.hmp.txt \
                                                  data/usda_22kSNPs_rils.sorted.diploid.hmp.txt
```

To make sure the reconstruction worked, I ran again `scripts/markers_summary.R`. Plots at `analysis/qc` folder show that, despite higher mising data relatively to other parents, the number of heterozygotes is now zero.

```bash
Rscript scripts/markers_summary.R data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.hmp.txt \
                                  analysis/qc/summary_markers_parents_PHG35-reconstructed.txt \
                                  data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                                  analysis/qc/summary_markers_rils.txt \
                                  data/usda_biparental-crosses.txt
```



#### RIL data

By looking at the boxplot of missing data of RILs, what stands out is that cross `B73xPHG47` had the highest variation and amount of missing data than other populations, and few populations showed very few RILs with more extreme missing data. Although the maximum amount of missing data for a RIL was about 9%, the median amount of missing data per cross was below 2.5%. For heterozygous markers, they were also mostly between 0 and 3% (which would be the expected for an F6 population), however population `B73xPHG47` also displayed a lot of variation and higher values of heterozygous markers. Very few RILs displayed high levels of heterozygosity, but those values were much higher than the expected (between 10 to 30%). Finally, I cannot see any parent that contributes with more heterozygosity or missing data in the RILs.

**Allele frequency of each marker**

With RILs from a biparental cross, we have the expectation that each segregating locus will have only two alleles with frequency ~50% each. Large deviations from this expectation may indicate some errors during genotyping and the locus should be removed from analysis.

Since I had 59 biparental crosses in `data/usda_biparental-crosses.txt`, it would be very annoying to perform the analysis one-by-one. Thus, I decided to use the command-line version of Tassel by using the following commands:

```bash
# save current field delimiter
curr_IFS=$IFS
{
  # skip header of "usda_biparental-crosses.txt" file
  read
  # set delimeter to tab
  IFS="\t"
  # read file line by line
  while read -r line; do
    # split line and assign each column to a variable
    # also change "*" by "x" in the cross name; otherwise my OneDrive won't upload that folder into the cloud
    cross=$(echo $line | cut -f1 | tr "*" "x")
    ril_list=$(echo $line | cut -f2)
    # check if directory exists; if it doesn't, create one to store results
    [[ -d analysis/qc/$cross ]] || mkdir -p analysis/qc/$cross
    # run tassel (added "\" at the end of the line just to improve readability)
    run_pipeline.pl -Xmx6g -importGuess data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                    -FilterTaxaBuilderPlugin -taxaList $ril_list -endPlugin \
                    -GenotypeSummaryPlugin -endPlugin \
                    -export analysis/qc/$cross/$cross\_OverallSummary,analysis/qc/$cross/$cross\_AlleleSummary,analysis/qc/$cross/$cross\_SiteSummary,analysis/qc/$cross/$cross\_TaxaSummary
  done
} < "data/usda_biparental-crosses.txt" > "analysis/qc/tassel_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir analysis/qc/*
```

> **Note**: this bash script also outputs a log file at `analysis/qc/tassel_log.txt`. I ran the TASSEL's summary to double check the results I got from `scripts/markers_summary.R`, and they are, indeed, the same.



**Recombination frequency**

I will use the sorted diploid hapmap files from parents (`data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.hmp.txt`) and RILs (`data/usda_22kSNPs_rils.sorted.diploid.hmp.txt`) to estimate recombination frequency for each biparental population. For this purpose, I wrote the `scripts/recomb-freq_biparental-crosses.R` script. Specifically, it does three things:

1. Filter the sorted diploid hapmap files to have genotypic information only from parents and RILs that make up a specific biparental cross. It takes information from the file `data/usda_biparental-crosses.txt` to do that, and output files on `data/biparental-crosses` (two files per cross: one hapmap for the parents, and another for RILs).

  > The `alleles` column of the filtered hapmap files are not correct, because it shows the alleles present in all parents and RILs (not only for that particular biparental cross). This might be corrected in future versions of the script, if needed.

2. Combine the hapmap files for each biparental population and converts them into the correct input format required by [rqtl](http://www.rqtl.org/), the R package that estimates the recombination frequency.

3. Run rqtl to estimate recombination frequencies in each biparental population, and plot the genetic distance (cM) by physical distance (Mb) for each cross as well. In fact, here is the complete list of things the function `EstimateRecombinationFreq()` does:

  * Read biparental cross information (genotypic data formated to rqtl).

  * Remove missing data (remove individual if half of the markers or more are missing, remove markers missing in half of individuals or more).

    > These cut-offs were arbitrary.

  * Remove duplicated individuals.

  * Remove duplicated markers.

  * Remove markers with segregation distortion (FDR adjusted p-value < 0.05).

  * Estimate recombination frequency.

  * **Plot** genotype frequencies by individual, and **write** genotypic data per RIL after all filtering.

  * **Plot** genetic and physical positions of markers per chromosome.

  * **Write summaries** of each cross before and after filtering.

  * **Write table** with genetic and physical positions of markers.

    > All ouput goes to the the respective **cross folder** in `analysis/qc/`. It also generates a log file called `rqtl_log.txt`.

```bash
Rscript scripts/recomb-freq_biparental-crosses.R data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.hmp.txt \
                                                 data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                                                 data/usda_biparental-crosses.txt \
                                                 analysis/qc \
                                                 data/biparental-crosses
```


**Allele frequencies of filtered markers**

I also wrote the file `scripts/usda_allele-freq_dist.R`. This script extracts the marker names from the recombination frequency tables for each cross, which were filtered to be polymorphic between parents and not show segregation distortion, to filter Tassel's `_SiteSummary` tables. This filtered table was used to make plots of the distributions of allele frequencies for each biparental population. The output was stored their respective **cross folder** in `analysis/qc`.

> Some of my first plots of the distributions had a weird shape. This was solved when I increased `binwidth` of histograms from `0.05` to `0.07`.

Then, as Candy suggested, I checked if markers that have allele frequency < 0.25 or > 0.75 were the same across populations. The majority of such markers (910) are unique of a population, while only 65 of these markers have extreme allele frequencies in more than one population.

```bash
Rscript scripts/usda_allele-freq_dist.R analysis/qc
```

**Karyotype of RILs**

I wrote `scripts/plot_ril_karyotypes.R`to see how the parental haplotypes are distributed in 5 randomnly selected RILs per population. For example, after looking at the karyotypes from RIL 2 (genotype frequency ~50/50) and RIL 3 (genotype frequency ~25/75) of cross B73 x LH82, I was not able to see any striking difference or something weird that can be causing the deviation on allele frequency. Overall, no major issues were observed in other "karyoplots".

```bash
Rscript scripts/plot_ril_karyotypes.R data/B73_RefGen_V4_chrm_info.txt \
                                      data/centromeres_Schneider-2016-pnas_v4.bed \
                                      analysis/qc
```

> It's important to check the number of recombinations per cromosome. At each generation you advance in RILs, you expect about 1-2 recombination events (although after later generations, F5-6, you should expect less because then you start recommbining blocks that have the same genotypes). For those RILs that have extreme genotypes, you should expect to see recombination more towards the end of chromosomes, that's why you would see a overrepressentation of one genotype (recombination is not always in the same place, it's a distribution).



### Remove low quality SNPs

The rqtl program detected 1,183 SNPs with segregation distortion (i.e. allele frequencies are too different from what is expected in RILs), and they should not be used as they can introduce bias to our downstream analysis.

```bash
# create file with problematic snps from all crosses
cat analysis/qc/*/markers-seg* analysis/qc/*/markers-missing* | sort | uniq > analysis/qc/rqtl_snps-to-remove_all-crosses.txt

# filter parental file
run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.hmp.txt \
                -excludeSiteNamesInFile analysis/qc/rqtl_snps-to-remove_all-crosses.txt \
                -export data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.hmp.txt \
                -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                -excludeSiteNamesInFile analysis/qc/rqtl_snps-to-remove_all-crosses.txt \
                -export data/usda_22kSNPs_rils.sorted.diploid.filtered.hmp.txt \
                -exportType HapmapDiploid
```



### Check reference genome version of SNP chip

One important thing that Candy remind me is to check whether the version of the reference genome used to make the probes in the SNP chip is the same as the one used in the SV calls (i.e., refgen v4). Since I will merge the two datasets, it's important that the coordinates are aligned.

To do that, I first downloaded the `.fasta` files of the chromosome 1 for the 4 refgen versions:

```bash
# create directory to save files
mkdir data/check_refgen_SNPchip
# go to that directory
cd data/check_refgen_SNPchip

# download chr1 from v1 to v4
wget http://ftp.maizesequence.org/release-4a.53/assembly/ZmB73_AGPv1_chr1.fasta.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/Zea_mays.AGPv2.dna.chromosome.1.fa
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-22/fasta/zea_mays/dna/Zea_mays.AGPv3.22.dna.chromosome.1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr1.fna.gz

# de-compress files
gunzip *.gz

# change names for easier identification
mv ZmB73_AGPv1_chr1.fasta v1_chr1.fasta
mv Zea_mays.AGPv2.dna.chromosome.1.fa v2_chr1.fasta
mv Zea_mays.AGPv3.22.dna.chromosome.1.fa v3_chr1.fasta
mv chr1.fna v4_chr1.fasta

# return to project's home directory
cd ~/projects/genomic_prediction/simulation
```

> Links above were obtained by searching the [MaizeGDB website](https://www.maizegdb.org/genome)

Then, I extracted the genotypic data of B73 (chromosome 1) from the SNP chip dataset using UNIX commands, to generate the file `data/check_refgen_SNPchip/markers_b73_chr1.txt`:

```bash
# create a temporary file with the hapmap header
head -n 1 data/usda_22kSNPs_parents.sorted.diploid.hmp.txt > data/check_refgen_SNPchip/tmp.hmp.txt
# keep only markers in the first chromosome
awk '$3 == 1' data/usda_22kSNPs_parents.sorted.diploid.hmp.txt >> data/check_refgen_SNPchip/tmp.hmp.txt
# select only position and B73 columns and create a new file
cut -f 4,12 data/check_refgen_SNPchip/tmp.hmp.txt > data/check_refgen_SNPchip/markers_b73_chr1.txt
# remove tmp file
rm data/check_refgen_SNPchip/tmp.hmp.txt
```

To find out which reference genome the SNP chip dataset was derived from, I wrote the `scripts/check_refgen_SNPchip.py`. This script extracts the positions from `data/check_refgen_SNPchip/markers_b73_chr1.txt` file and writes the respective positions for each refgen version. Additionally, it counts how many matches there are between each refgen and the SNP chip.

```bash
python scripts/check_refgen_SNPchip.py data/check_refgen_SNPchip/markers_b73_chr1.txt data/check_refgen_SNPchip
```

| refgen     | matches          |
| ---------- | ---------------- |
| B73_v1     | 774 (24.2%)      |
| **B73_v2** | **2605 (81.3%)** |
| B73_v3     | 786 (24.5%)      |
| B73_v4     | 788 (24.6%)      |


The results suggest that refgen version 2 was used to design the SNP chip, since there is ~80% match. The rest ~20% may be due to the probes used in the chip be tagging the opposite strand of the DNA. To test that, I wrote a short script called `scripts/check_refgen_SNPchip_strands.R` and noticed that **all markers** that didn't match v2 ref gen was actually the complementary nucleotide (while for other reference genome versions it was only ~30% of markers matching with complementary base):

```bash
Rscript scripts/check_refgen_SNPchip_strands.R data/check_refgen_SNPchip/ref-markers_chr1.txt

# Same alleles between marker and v1 is 0.2416485
# Missing alleles removed: 35
# Reverse complement alleles between marker and v1 is 0.3187135

# Same alleles between marker and v2 is 0.8133
# Missing alleles removed: 35
# Reverse complement alleles between marker and v2 is 1

# Same alleles between marker and v3 is 0.2453949
# Missing alleles removed: 35
# Reverse complement alleles between marker and v3 is 0.3375315

# Same alleles between marker and v4 is 0.2460194
# Missing alleles removed: 34
# Reverse complement alleles between marker and v4 is 0.3175136
```

At this point, we are pretty confident that the probes designed for the SNP chip were based on the **reference genome assembly v2**.




### Convert SNP coordinates from B73v2 to B73v4

Now, I need to convert the SNP chip coordinates from refgen v2 to v4 before merging this SNP data with the SV data. The first thing I have to do is download the refgen v2 assembly and extract 100bp sequences around SNPchip probe positions.

```bash
# make directory to store refgen sequence, and go to that directory
mkdir -p data/refgen/B73v2
cd data/refgen/B73v2

# download v2 sequences for all chromosomes
for chr in {1..10}; do
  wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/Zea_mays.AGPv2.dna.chromosome.$chr.fa
done
# download README
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/README

# return to project's home directory
cd ~/projects/genomic_prediction/simulation
# extract probes sequences
python scripts/extract_probe_seqs.py data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.hmp.txt \
                                     data/refgen/B73v2/ \
                                     100 \
                                     data/probes-100bp_22kSNPchip_B73v2.fa
```

Once I have a fasta file with all probe sequences, I can align them to the refgen v4 assembly using `bowtie`. But first, I need to download the v4 assembly and build an index (by running `scripts/bowtie_index_refgenv4.sh`).

```bash
# make directory to store refgen sequence, and go to that directory
mkdir -p data/refgen/B73v4
cd data/refgen/B73v4

# download sequences for all chromosomes
for chr in {1..10}; do
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr$chr.fna.gz
done
# decompress them
gunzip *.gz
# change extension name .fna to .fa for bowtie compatibility
for file in *.fna; do
    mv -- "$file" "${file%.fna}.fa"
done

# return to project's home directory
cd ~/projects/genomic_prediction/simulation

# build index for refgen v4
qsub scripts/bowtie_index_refgenv4.sh
```

Once the index is built, mapping reads with bowtie is pretty straightforward:

```bash
# load bowtie
module load bowtie/1.1.2

# map probes to refgen v4
# bowtie [options] [index_basename] [unpaired_reads] [output_file]
bowtie -f -v 0 -m 1 --un data/probes-100bp_22kSNPchip_not-aligned-to-B73v4.fa -S \
       data/refgen/B73v4/index \
       data/probes-100bp_22kSNPchip_B73v2.fa \
       data/probes-100bp_22kSNPchip_aligned-to-B73v4.sam 2> data/probes-100bp_22kSNPchip_aligned-to-B73v4.stats.txt
```

> Bowtie options:
> * `-f`: input is a fasta file
> * `-v 0`: map end-to-end with 0 mismatches allowed
> * `-m 1`: don't report alignments if read map in more than 1 place
> * `--un data/probes-100bp_22kSNPchip_not-aligned-to-B73v4.fa`: unmapped reads are written to this file
> * `-S`: output is SAM
> * `data/refgen/B73v4/index`: basename of the ref gen v4 index
> * `data/probes-100bp_22kSNPchip_B73v2.fa`: reads to map
> * `data/probes-100bp_22kSNPchip_aligned-to-B73v4.sam`: output name
> * `2> data/probes-100bp_22kSNPchip_aligned-to-B73v4.stats.txt`: this will print alignment statistics from std err to this file

Here are the stats of the alignment:

| Reads                   | Stats  |
| ----------------------- | ------ |
| Processed               | 18,956 |
| Uniquely mapped         | 18,683 |
| Multimapped (supressed) | 64     |
| Unmapped                | 209    |

Now, I need to correct the positions of the SNPs in the hapmap files `data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.hmp.txt` and `data/usda_22kSNPs_rils.sorted.diploid.filtered.hmp.txt` based on new coordinates from SAM file. I wrote a python script called `scripts/convert_hmp_v2-to-v4.py` to do that and create the files `usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt` and `usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt`, and another file with SNP coordinates in both v2 and v4. This script also removes reads that uniquely mapped with no mismatch to a different chromosome in v4, since it's very unlikely that there would be such major differences between v2 and v4 assemblies (i.e., these are probably mismapped).

```bash
python scripts/convert_hmp_v2-to-v4.py data/probes-100bp_22kSNPchip_aligned-to-B73v4.sam \
                                       data/SNP_positions_v2-to-v4_probes-100bp.txt \
                                       data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.hmp.txt,data/usda_22kSNPs_rils.sorted.diploid.filtered.hmp.txt
# 22 reads discarded due to mapping into different chromosome in v4
```

I performed a quick QC (using only chromosome 1 data) to see if SNPs from v4 were actually the same as in v2 by running `scripts/check_refgen_SNPchip_v2tov4_probes.R`. The results showed that some markers that didn't match between v2 and v4, and that these mismatches were due to mapping to complementary strand:

```bash
Rscript scripts/check_refgen_SNPchip_v2tov4_probes.R data/check_refgen_SNPchip/ref-markers_chr1.txt data/SNP_positions_v2-to-v4_probes-100bp.txt
# 201 marker alleles differ between v2 and v4
# Converting alleles that differ to complementary base...
# 0 marker alleles differ between v2 and v4
```

Thus, I corrected the strand of those markers in the hapmap files `data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt` and `data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt` by running `scripts/correct_SNP_strands.R`. Finally, I corrected the alleles' column of these hapmap files using Tassel.

```bash
# correct markers not phased
Rscript scripts/correct_SNP_strands.R data/SNP_positions_v2-to-v4_probes-100bp.txt \
                                      data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt \
                                      data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt

# correct alleles column
run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt \
                        -export data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt \
                        -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt \
                        -export data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt \
                        -exportType HapmapDiploid
```



### Divide dataset by cross

Both SV and SNP projections will be done for each family separately. Thus, I have to create parental and RIL SNP chip files (the anchors for projections) for each specific cross with `scripts/divide_hmp_by_cross.R`.

```bash
# parents
Rscript scripts/divide_hmp_by_cross.R data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt \
                                      data/usda_biparental-crosses.txt \
                                      data/hapmap_by_cross \
                                      --parents
# rils
Rscript scripts/divide_hmp_by_cross.R data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt \
                                      data/usda_biparental-crosses.txt \
                                      data/hapmap_by_cross \
                                      --rils
```



## SV dataset

On August 9, 2019, Patrick Monnahan sent me five `.vcf` files containing structural variation calls from the software Lumpy. Each file is a SV call of 100 lines against one of the following reference genomes: B73, Mo17, PH207, PHB47, or W22. **This results are preliminary** and there are a lot of false-positives in there. However, they will be useful for me to select the 7 lines I need, project the SVs from parents to RILs, and incorporate these SV data into my simulation scripts.

The files are located at `data/SV_calls`, and I will use only the SVs called against the B73 reference genome for now.

```bash
# decompress vcf file of SV calls to B73
gunzip data/SV_calls/B73v4_2019-08-09.ls.RT.vcf.gz
```



### Hapmap format

I wrote the python script `scripts/vcf2hapmap.py` which extract information about structural variants and transform into a hapmap format file sorted by chromosome and positions. Since the `.vcf` file that Patrick sent me contains 100 inbred lines, I have to select only the 7 inbred parents used in the USDA project. Running this script will generate the file `data/usda_SVs_parents.sorted.hmp.txt`:

```bash
# for help on how to use this script
python scripts/vcf2hapmap.py --help
# run the script
python scripts/vcf2hapmap.py data/SV_calls/B73v4_2019-08-09.ls.RT.vcf \
                             data/usda_SVs_parents.sorted.hmp.txt \
                             B73,LH82,PH207,PHG35,PHG39,PHG47,PHJ40
```

The type of SV will be displayed in the marker ID (`del.[ID]` for deletions, `dup.[ID]` for duplications, etc.). Each line will have either a value of `AA` if SV is `A`bsent, or `TT` if SV is `T`here. I had to do that because I will use Tassel to project genotypes from parents to RILs, and it doesn't accept anything other than nucleotides (thus, I couldn't use numbers to indicate presence/absence of SV, as originally thought). Missing data will be coded as `NN`. Also, since SVs spam hundreds (or thousands) of bp and the exact breakpoints are hard to call, the position indicated in the hapmap file wil be the middle point of the SV.

> After projection of parental genotypes into RILs, the hapmap will be transformed into the numeric format for genomic prediction simulations. Thus, for SVs, `AA` will be transformed into `0` and `TT` to `2`.



### Divide dataset by cross

Similarly to the anchors datasets (SNP chip), I have to create parental files with SV information separately for each cross with `scripts/divide_hmp_by_cross.R`.

```bash
# parents
Rscript scripts/divide_hmp_by_cross.R data/usda_SVs_parents.sorted.hmp.txt \
                                      data/usda_biparental-crosses.txt \
                                      data/hapmap_by_cross \
                                      --parents
```



## Remove chip SNPs inside deletions

For projections, I need to create a parental and RIL files with all SNPs and SVs. However, before that, I need to address a particular issue when using these two types of variation: SNPs that are found inside deletions. Such SNPs are problematic because they will have segregation issues when you compare multiple lines that have or does not have that SV.

The `scripts/remove_SNPs_within_dels.R` will remove SNPs that are within 100kb of a deletion first. I set up 100 kb threshold because ~99% of the deletions are smaller than this and the very large deletions (>1 Mb) would make me lose a lot of SNPs. Importantly, each family will be filtered separately and a SNP will be removed only if the deletion is actually present in one of the parents (i.e. there is `TT` call in one of the parents).

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

# remove snps
for cross in $crosses; do
  echo ${cross}
  # parents
  Rscript scripts/remove_SNPs_within_dels.R data/hapmap_by_cross/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.${cross}.hmp.txt \
                                            data/hapmap_by_cross/usda_SVs_parents.sorted.${cross}.hmp.txt
  # rils
  Rscript scripts/remove_SNPs_within_dels.R data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross}.hmp.txt \
                                            data/hapmap_by_cross/usda_SVs_parents.sorted.${cross}.hmp.txt
done

# check parents and rils have same number of markers
for cross in $crosses; do
  wc -l data/hapmap_by_cross/usda_22kSNPs_parents.*.${cross}.not-in-SVs.hmp.txt
  wc -l data/hapmap_by_cross/usda_22kSNPs_rils.*.${cross}.not-in-SVs.hmp.txt
  echo ""
done

# check how many were removed per cross
for cross in $crosses; do
  before=$(wc -l data/hapmap_by_cross/usda_22kSNPs_parents.*.${cross}.hmp.txt | cut -d " " -f 1)
  after=$(wc -l data/hapmap_by_cross/usda_22kSNPs_parents.*.${cross}.not-in-SVs.hmp.txt | cut -d " " -f 1)
  echo "$((before-after)) SNPs removed for $cross"
done
```

> Only very few SNPs (1-6 SNPs) were inside the boundaries of a deletion for each family and thus removed.




## Correct SNP chip miscalls with sliding window approach

Before doing projections, it's also important to minimize the number of miscalls in the SNP chip data for RILs, otherwise the number of wrong projections may be high, which can then affect downstream predictions.

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in $crosses; do
  Rscript scripts/sliding_window_approach.R $cross \
                                            data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross}.not-in-SVs.hmp.txt \
                                            data/hapmap_by_cross/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.${cross}.not-in-SVs.hmp.txt \
                                            --window_size=15 --window_step=1 --min_snps_per_window=5
done
```

> Note: In preliminary tests, it looked like that as I increased `--window_size`, the number of hets increased as well, especially around recombination breakpoints. Increasing `--window_step` kind of controlled that a bit (at least visually), but this end up increasing the chance of having wrong recombination points.



## Resequencing SNPs

After some preliminary work on genomic prediction models with simulated data, we noticed that we end up with far more SVs (~10k) than SNPs in LD to an SV (~3k). This very likely happened because the SNP chip doesn't have the SNP density required to fully capture every SNP in LD to an SV. Fortunately, there is resequencing data available for all parents of the USDA population, and this type of data can provide the resolution needed. However, since there is no resequencing data for the RILs, I have to project these resequencing SNPs from USDA parents to respective RILs.

All parents for this population were resequenced by [Mazaheri et al (BMC Plant Biology, 2019)](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-019-1653-x), and the data can be retrieved in [this Dryad repository](https://datadryad.org/stash/landing/show?big=showme&id=doi%3A10.5061%2Fdryad.n0m260p).

> Note that the real link for for downloading this data will be sent by email.

The main file to use for this project is called `62biomAP_v_B73_SNPMatrix.txt`, which is a matrix containing all SNP calls for 56 maize inbred lines. Thus, after downloading the data, I have to select only the parents of the USDA population, transform the data to the hapmap format and then proceed to the rest of the analysis.

Unzip `62biomAP_v_B73_SNPMatrix.txt.txt.tar.gz` and it's has all the SNPs from 56 inbreds. I will just select the SNPs from the 7 parents of my population (B73, LH82, PH207, PHG35, PHG39, PHG47, and PHJ40). The dataset to be used is called `data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.txt`.

```bash
mkdir -p data/reseq_snps
cd data/reseq_snps

# note that if I download this data again, the name of the file will be different
wget http://merritt.cdlib.org/cloudcontainer/mrtstore2/21061229.tar.gzs
tar -xvf 21061229.tar.gz
gunzip 62biomAP_v_B73_SNPMatrix.txt.gz

# the only thing I care is the SNP matrix, so I'm gonna delete the other files
rm *widiv_942g_899784SNPs*

# get column number of the inbreds
head -n 1 62biomAP_v_B73_SNPMatrix.txt | tr '\t' '\n' | cat -n | grep "B73\|LH82\|PH207\|PHG35\|PHG39\|PHG47\|PHJ40"
# 2	B73v4_Ref
# 5	B73
# 21	LH82
# 32	PH207
# 36	PHG35
# 37	PHG39
# 38	PHG47
# 43	PHJ40

# get first column (position) and the ones above
cut -f 1,2,5,21,32,36,37,38,43 62biomAP_v_B73_SNPMatrix.txt > biomAP_v_B73_SNPMatrix_parents.txt

# return to project's home directory
cd ~/projects/genomic_prediction/simulation

# transform to hapmap format:
python scripts/snp_matrix2hapmap.py data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.txt data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.hmp.txt

# transform to hapmap diploid format:
run_pipeline.pl -Xmx50g -importGuess data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.hmp.txt \
                        -export data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.hmp.txt \
                        -exportType HapmapDiploid
```


### Divide dataset by cross

Similarly to the anchors datasets (SNP chip), I have to create parental files with resequencing SNPs information separately for each cross with `scripts/divide_hmp_by_cross.R`.

```bash
# parents
Rscript scripts/divide_hmp_by_cross.R data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.hmp.txt \
                                      data/usda_biparental-crosses.txt \
                                      data/reseq_snps \
                                      --parents
```




### Remove resequecing SNPs inside deletions and keep only polymorphic SNPs

Since the SNP resequencing data will be used together with SV data, I also need to make sure that such SNPs are not present within the boundaries of an SV. The `scripts/filter_reseq_snps.sh` will first remove resequencing SNPs that are within 100kb of a deletion with `scripts/remove_SNPs_within_dels.R`, and then it will keep only the polymorphic SNPs with `scripts/keep_only_poly_snps.R` to speed up projections (as it reduces a lot the total number SNPs to be projected).

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in $crosses; do
  qsub -v CROSS=${cross} scripts/filter_reseq_snps.sh
done
```

> Less than 1% of resequencing SNPs were within deletions up to 100 kb.




## Merge SNP and SV data

The markers from the SNP chip will be used as anchors for projection of both SVs and resequencing SNPs. Since these two types of markers will be projected at the same type, I need to merge the hapmap files of both marker types with the SNP chip files. First, I will combine SNP chip with SV data for parents and RILs separately with `scripts/merge_SNPs_SVs_hapmaps.R`. It's worth noting that genotypic information of SVs is only available on parental dataset. The RIL dataset have the SV names but not the genotypes (they are all `NN`).

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

# merge svs
for cross in $crosses; do
  Rscript scripts/merge_SNPs_SVs_hapmaps.R data/hapmap_by_cross/usda_SVs_parents.sorted.${cross}.hmp.txt \
                                           data/hapmap_by_cross/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.${cross}.not-in-SVs.hmp.txt \
                                           data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross}.not-in-SVs.sliding-window.hmp.txt \
                                           data/hapmap_by_cross/usda_parents_SV-SNPchip.${cross}.hmp.txt \
                                           data/hapmap_by_cross/usda_rils_SV-SNPchip.${cross}.not-projected.hmp.txt
done

# make sure the files are sorted
for cross in $crosses; do
  # parents
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile data/hapmap_by_cross/usda_parents_SV-SNPchip.${cross}.hmp.txt -outputFile data/hapmap_by_cross/usda_parents_SV-SNPchip.${cross}.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess data/hapmap_by_cross/usda_parents_SV-SNPchip.${cross}.hmp.txt -export data/hapmap_by_cross/usda_parents_SV-SNPchip.${cross}.hmp.txt -exportType HapmapDiploid
  # rils
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile data/hapmap_by_cross/usda_rils_SV-SNPchip.${cross}.not-projected.hmp.txt -outputFile data/hapmap_by_cross/usda_rils_SV-SNPchip.${cross}.not-projected.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess data/hapmap_by_cross/usda_rils_SV-SNPchip.${cross}.not-projected.hmp.txt -export data/hapmap_by_cross/usda_rils_SV-SNPchip.${cross}.not-projected.hmp.txt -exportType HapmapDiploid
done

# make sure number of rows of parental and RIL data matches in each population
for cross in $crosses; do
  wc -l data/hapmap_by_cross/usda_parents_SV-SNPchip.${cross}.hmp.txt
  wc -l data/hapmap_by_cross/usda_rils_SV-SNPchip.${cross}.not-projected.hmp.txt
  echo ""
done
```

Then, in order to project resequencing SNPs and SVs at the same time, I will merge the poymorphic resequencing SNPs with the previous files (SNP chip + SVs) for parents and RILs separately with `scripts/merge_SNPs_SVs_hapmaps.R`. Similar to previous merging, genotypic information of reseq SNPs is only available on parental dataset. The RIL dataset have the resequencing SNP names but not the genotypes (they are all `NN`). Additionally, since the same SNP information can be represented in both SNP chip and resequencing SNP data, I wrote `scripts/remove_duplicated_anchors.R` to keep only SNPs from the chip (anchors).

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

# merge poly reseq snps with snp chip markers + svs hapmaps
for cross in $crosses; do
  Rscript scripts/merge_SNPs_SVs_hapmaps.R data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.${cross}.poly.not-in-SVs.hmp.txt \
                                           data/hapmap_by_cross/usda_parents_SV-SNPchip.${cross}.hmp.txt \
                                           data/hapmap_by_cross/usda_rils_SV-SNPchip.${cross}.not-projected.hmp.txt \
                                           data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${cross}.hmp.txt \
                                           data/hapmap_by_cross/usda_rils_SV-SNPchip-polySNPreseq.${cross}.not-projected.hmp.txt
done

# remove duplicated anchors (which were created when merging GBS SNPs with )
for cross in $crosses; do
  Rscript scripts/remove_duplicated_anchors.R data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${cross}.hmp.txt
  Rscript scripts/remove_duplicated_anchors.R data/hapmap_by_cross/usda_rils_SV-SNPchip-polySNPreseq.${cross}.not-projected.hmp.txt
done

# make sure the files are sorted
for cross in $crosses; do
  # parents
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${cross}.hmp.txt -outputFile data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${cross}.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${cross}.hmp.txt -export data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${cross}.hmp.txt -exportType HapmapDiploid
  # rils
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile data/hapmap_by_cross/usda_rils_SV-SNPchip-polySNPreseq.${cross}.not-projected.hmp.txt -outputFile data/hapmap_by_cross/usda_rils_SV-SNPchip-polySNPreseq.${cross}.not-projected.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess data/hapmap_by_cross/usda_rils_SV-SNPchip-polySNPreseq.${cross}.not-projected.hmp.txt -export data/hapmap_by_cross/usda_rils_SV-SNPchip-polySNPreseq.${cross}.not-projected.hmp.txt -exportType HapmapDiploid
done

# check that parents and rils have the same number of markers
for cross in $crosses; do
  wc -l data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${cross}.hmp.txt
  wc -l data/hapmap_by_cross/usda_rils_SV-SNPchip-polySNPreseq.${cross}.not-projected.hmp.txt
  echo ""
done
```



## Project SVs and resequencing SNPs from parents to RILs

Projection of SVs and resequencing SNPs will be done on a family basis. The merged SNP chip, SV and resequencing SNP data from parents will serve as donors for projections by [FILLIN](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/FILLIN/FILLIN). The haplotype blocks generated for each parent will then be used to determine the haplotypes for projection into the RIL dataset.

> After preliminary analysis, using the options `-hapSize 200000` (haplotype block size), `-minTaxa 1` (create haplotype for each parent), and `-hybNN false` (avoid combining haplotypes) yielded the best projection results.

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_svs-snps

# project svs and snps together
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in $crosses; do
  qsub -v CROSS=$cross,HAPSIZE=200000 scripts/project_svs-snps.sh
done

# transform to diploid hapmap:
cd ~/projects/genomic_prediction/simulation/analysis/projection_svs-snps

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done

# go back to project's home folder
cd ~/projects/genomic_prediction/simulation
```



### Sliding window approach

Given the high amount of SNPs projected, I will run the sliding window approach to reduce the possible wrong projections using a bigger `--window-size`, and then plot the karyotypes after this approach to see if it was helpful or not.

```bash
# sliding window
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in $crosses; do
  qsub -v CROSS=$cross scripts/sliding_window_reseq-snps.sh
done

# make sure the number of SNPs match before and after sliding window
for cross in $crosses; do
  echo $cross
  wc -l analysis/projection_reseq-snps/*$cross*projected.hmp.txt
  wc -l analysis/projection_reseq-snps/*$cross*sliding-window.hmp.txt
  echo ""
done
```

To summarize the results of both SV and SNP projections for each family after the sliding window approach, I wrote `scripts/count_projected_SVs.R` and `scripts/count_projected_reseq-snps.R`, respectively. The average **SV projection was 93.2%** and the average **SNP projection was 94%**, with average accuracy of 94.7%% across all families.

Details about projections for each family were written into `analysis/projection_svs-snps/projection_svs_summary.txt` for SVs and `analysis/projection_svs-snps/projection_svs_summary.txt` for SNPs, but can also be visualized in different plots saved at `analysis/projection_svs-snps`.

```bash
# count svs
Rscript scripts/count_projected_SVs.R data/hapmap_by_cross \
                                      analysis/projection_svs-snps \
                                      data/usda_biparental-crosses.txt

# get average SV projection and accuracy -- use (NR - 1) to avoid counting header as a row
awk -v N=5 '{ sum += $N } END { if (NR > 1) print sum / (NR - 1) }' analysis/projection_svs-snps/projection_svs_summary.txt
awk -v N=7 '{ sum += $N } END { if (NR > 1) print sum / (NR - 1) }' analysis/projection_svs-snps/projection_svs_summary.txt

# count reseq snps
Rscript scripts/count_projected_reseq-snps.R data/reseq_snps \
                                             analysis/projection_svs-snps \
                                             data/usda_biparental-crosses.txt

# get average SNP projection and accuracy -- use (NR - 1) to avoid counting header as a row
awk -v N=5 '{ sum += $N } END { if (NR > 1) print sum / (NR - 1) }' analysis/projection_svs-snps/projection_reseq-snps_summary.txt
awk -v N=7 '{ sum += $N } END { if (NR > 1) print sum / (NR - 1) }' analysis/projection_svs-snps/projection_reseq-snps_summary.txt
```

I also ploted the karyotypes from RILs before and after SNP projection with `scripts/plot_ril_karyotypes_reseq-SNPs.R`. After running it, I don't see any major disagreements between parental and RIL data, meaning that projections seem to be very accurate indeed.

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in $crosses; do
  # svs
  Rscript scripts/plot_ril_karyotypes_reseq-SNPs.R data/B73_RefGen_V4_chrm_info.txt \
                                                   data/centromeres_Schneider-2016-pnas_v4.bed \
                                                   $cross \
                                                   analysis/projection_svs-snps \
                                                   data/hapmap_by_cross \
                                                   analysis/qc/karyotypes \
                                                   --rils=random --marker-type=sv --sliding-window
  # snps
  Rscript scripts/plot_ril_karyotypes_reseq-SNPs.R data/B73_RefGen_V4_chrm_info.txt \
                                                   data/centromeres_Schneider-2016-pnas_v4.bed \
                                                   $cross \
                                                   analysis/projection_svs-snps \
                                                   data/hapmap_by_cross \
                                                   analysis/qc/karyotypes \
                                                   --rils=random --marker-type=snp --sliding-window
done
```

> TODO: fix `scripts/plot_ril_karyotypes_reseq-SNPs.R` and also make sure to change output name based on the type of marker.




### Add monomorphic resequencing SNPs back

Since the projections were performed only with polymorphic SNPs to reduce computational time, now I have to add back the monomorphic SNPs for each family.

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in $crosses; do
  qsub -v CROSS=$cross scripts/add_mono_reseq-snps.sh
done
```



### Add back SNPs and SVs not present in a certain cross

Given that SNPs within SVs were removed for each family separately and that SVs with missing data in both parents of a family were also removed, each family will have slightly different marker composition. Thus, I wrote `scripts/add_markers_not_in_cross.R` to make sure all the families have the same markers. Markers with no information for a given cross was set as missing data (`NN`).

```bash
qsub scripts/add_markers_not_in_cross.sh

# make sure the final file is sorted and the alleles' column is correct
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in $crosses; do
  qsub -v CROSS=$cross scripts/sort_merged_projected_files.sh
done
```



### Merge projected markers from different families

Finally, I will merge marker projections from each family into a single hapmap file (`data/usda_rils_projected-SVs-SNPs.hmp.txt`), which can then be used for LD analysis and genomic predictions.

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

# exclude first columns for all crosses
cd analysis/projection_svs-snps
for cross in $crosses; do
  echo $cross
  cut -f 1-11 --complement usda_rils_projected-SVs-SNPs.$cross.all-markers.hmp.txt > tmp_hmp.$cross.txt
done
# join all rils in one file (keep entire hmp file for cross B73xLH82 though)
paste usda_rils_projected-SVs-SNPs.B73xLH82.all-markers.hmp.txt \
      tmp_hmp.B73xPH207.txt \
      tmp_hmp.B73xPHG35.txt \
      tmp_hmp.B73xPHG39.txt \
      tmp_hmp.B73xPHG47.txt \
      tmp_hmp.LH82xPH207.txt \
      tmp_hmp.LH82xPHG35.txt \
      tmp_hmp.LH82xPHG39.txt \
      tmp_hmp.LH82xPHG47.txt \
      tmp_hmp.PH207xPHG47.txt \
      tmp_hmp.PHG35xPHG39.txt \
      tmp_hmp.PHG35xPHG47.txt \
      tmp_hmp.PHG39xPHG47.txt > ~/projects/genomic_prediction/simulation/data/usda_rils_projected-SVs-SNPs.hmp.txt
# remove tmp files
rm tmp_hmp.*.txt

# go back to project's home folder
cd ~/projects/genomic_prediction/simulation

# also separate the final file into chromosomes
for chr in {1..10}; do
  echo "$chr"
  head -n 1 data/usda_rils_projected-SVs-SNPs.hmp.txt > data/usda_rils_projected-SVs-SNPs.chr$chr.hmp.txt
  awk -v chr="$chr" '$3 == chr' data/usda_rils_projected-SVs-SNPs.hmp.txt >> data/usda_rils_projected-SVs-SNPs.chr$chr.hmp.txt
done

# generate summary statistics for the final file
for chr in {1..10}; do
  qsub -v CHR=$chr scripts/tassel_qc-summary.sh
done

# plot maf and missing data distribution (just need to call one chr -- but the others should be in same folder)
Rscript scripts/plot_maf-missing_distribution.R analysis/qc/tassel_usda-SNPs-SVs_chr1_summary3.txt analysis/qc/maf_missing-data
```



## LD structure between SNPs and SVs

An important thing to do is to calculate the LD between SNPs and SVs, because if SNPs are in high LD with SVs, then genomic prediction with SVs will likely add little information to the models. Since I don't know the LD structure in the population and due to the high amount of resequencing SNPs, I'm calculating LD in a 1kb, 10kb, 100kb, and 1Mb window. Additionally, I'm filtering out variants with more than 0.25 missing data to make sure that LD is not very much influenced by missing data.

> LD decay in maize is somewhere around 5kb, so the window range tested should be enough to find out the best window

```bash
# create directory to store results
mkdir -p /scratch.global/della028/hirsch_lab/genomic_prediction/ld
mkdir -p analysis/ld

# convert hapmap to plink format
for chr in {1..10}; do
  qsub -v CHR=$chr scripts/hmp2plk.sh
done

# calculate LD
for chr in {1..10}; do
  for filter in 0.25; do
    for window in 1 10 100 1000; do
      [[ ${window} == 1000 ]] && res="walltime=24:00:00,nodes=1:ppn=1,mem=150gb" || res="walltime=8:00:00,nodes=1:ppn=1,mem=100gb"
      qsub -l ${res} -v CHR=$chr,WINDOW=$window,FILTER=$filter scripts/plink_ld_calculation.sh
    done
  done
done
```

Besides that, it's important to QC the dataset in different ways to get a broader picture of how the LD structure in my population can affect genomic predition.



### Distribution

First, I will plot the distribution of r2 values from LD between SV and its SNP with highest LD to see if all SVs are perfectly tagged by an SNP or not. Additionally, some summary statistics of LD between SNPs and SVs are also calculated.

```bash
# get only SNP-SV LD
for chr in {1..10}; do
  for filter in 0.25; do
    for window in 1 10 100 1000; do
      qsub -v CHR=$chr,WINDOW=$window,FILTER=$filter scripts/select_only_sv-snp_ld.sh
    done
  done
done

# filter plink ld files and plot distribution (only need the first chr file -- it will find the other chr)
for filter in 0.25; do
  for window in 1 10 100 1000; do
    qsub -v WINDOW=$window,FILTER=$filter scripts/distribution_ld_snps-svs.sh
  done
done

# plot extra stats for markers in high ld
for filter in 0.25; do
  for window in 1 10 100 1000; do
    echo $window
    Rscript scripts/extra_ld_stats.R analysis/ld/window-${window}kb_filter-${filter}
  done
done
```

> Most of SNPs in high LD with SVs have r2 > 0.8, and the number increases as window size increses (probably due to population structure).

Second, I will plot the distribution of each marker type (SNP or SV) along the chromosome to see if there's enrichment of a particular marker for certain parts of the chromosomes:

```bash
# create a 3 column file (id, chr, pos) to be the input for R script
for chr in {1..10}; do
  echo "$chr"
  for filter in 0.25; do
    for window in 1 10 100 1000; do
      # highest ld
      sed 1d analysis/ld/window-$window\kb_filter-$filter/ld_usda_rils_snp-sv_only.chr$chr.window-$window\kb.filter-$filter.highest-ld.ld  | awk '{ print $3 "\t" $1 "\t" $2 "\n" $7 "\t" $5 "\t" $6}' | sort -n -k 3 | uniq >> analysis/ld/window-$window\kb_filter-$filter/marker_info_highest-ld.txt
      # not in ld
      sed 1d analysis/ld/window-$window\kb_filter-$filter/ld_usda_rils_snp-sv_only.chr$chr.window-$window\kb.filter-$filter.not-in-ld.ld  | awk '{ print $3 "\t" $1 "\t" $2 "\n" $7 "\t" $5 "\t" $6}' | sort -n -k 3 | uniq >> analysis/ld/window-$window\kb_filter-$filter/marker_info_not-in-ld.txt
    done
  done
done

# plot distribution along chromosome
for filter in 0.25; do
  for window in 1 10 100 1000; do
    # highest ld
    Rscript scripts/distribution_snps-svs_chrom.R analysis/ld/window-$window\kb_filter-$filter/marker_info_highest-ld.txt analysis/ld/window-$window\kb_filter-$filter/distribution_snps-svs_chrom.png
    # not in ld
    Rscript scripts/distribution_snps-svs_chrom.R analysis/ld/window-$window\kb_filter-$filter/marker_info_not-in-ld.txt analysis/ld/window-$window\kb_filter-$filter/distribution_snps-not-ld-svs_chrom.png
  done
done

# count how many SVs have an SNP in high LD within different window sizes
for window in 1 10 100 1000; do
  count=$(grep -c -P "^del|^ins|^inv|^dup|^tra" analysis/ld/window-${window}kb_filter-0.25/marker_info_highest-ld.txt)
  echo "$count SNPs in high LD with a SV within $window kb"
done
# 5472 SNPs in high LD with a SV within 1 kb
# 7609 SNPs in high LD with a SV within 10 kb
# 7975 SNPs in high LD with a SV within 100 kb
# 7977 SNPs in high LD with a SV within 1000 kb

# also plot distribution along chromosome for anchors only
for chr in {1..10}; do
  cut -f 1,3,4 data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt | awk -v chr=$chr '$2 == chr' > analysis/projection_svs-snps/anchor_positions.chr$chr.txt
done
Rscript scripts/distribution_snps-svs_chrom.R analysis/projection_svs-snps/anchor_positions.chr$chr.txt analysis/projection_svs-snps/distribution_anchors_chrom.png

# also plot distribution snps in highest ld to svs along chromosome for by ld quarter
for chr in {1..10}; do
  for filter in 0.25; do
    for window in 1 10 100 1000; do
      # r2 < 0.25
      sed 1d analysis/ld/window-$window\kb_filter-$filter/ld_usda_rils_snp-sv_only.chr$chr.window-$window\kb.filter-$filter.highest-ld.ld | awk -v chr=$chr '$1 == chr' | awk '$9 < 0.25' | awk '{ print $3 "\t" $1 "\t" $2 "\n" $7 "\t" $5 "\t" $6}' | sort -n -k 3 | uniq >> analysis/ld/window-$window\kb_filter-$filter/marker_info_r2_0_0.25.window-${window}kb.filter-${filter}.txt
      # 0.25 < r2 < 0.5
      sed 1d analysis/ld/window-$window\kb_filter-$filter/ld_usda_rils_snp-sv_only.chr$chr.window-$window\kb.filter-$filter.highest-ld.ld | awk -v chr=$chr '$1 == chr' | awk '$9 > 0.25 && $9 < 0.5' | awk '{ print $3 "\t" $1 "\t" $2 "\n" $7 "\t" $5 "\t" $6}' | sort -n -k 3 | uniq >> analysis/ld/window-$window\kb_filter-$filter/marker_info_r2_0.25_0.5.window-${window}kb.filter-${filter}.txt
      # 0.5 < r2 < 0.75
      sed 1d analysis/ld/window-$window\kb_filter-$filter/ld_usda_rils_snp-sv_only.chr$chr.window-$window\kb.filter-$filter.highest-ld.ld | awk -v chr=$chr '$1 == chr' | awk '$9 > 0.5 && $9 < 0.75' | awk '{ print $3 "\t" $1 "\t" $2 "\n" $7 "\t" $5 "\t" $6}' | sort -n -k 3 | uniq >> analysis/ld/window-$window\kb_filter-$filter/marker_info_r2_0.5_0.75.window-${window}kb.filter-${filter}.txt
      # r2 > 0.75
      sed 1d analysis/ld/window-$window\kb_filter-$filter/ld_usda_rils_snp-sv_only.chr$chr.window-$window\kb.filter-$filter.highest-ld.ld | awk -v chr=$chr '$1 == chr' | awk '$9 > 0.75' | awk '{ print $3 "\t" $1 "\t" $2 "\n" $7 "\t" $5 "\t" $6}' | sort -n -k 3 | uniq >> analysis/ld/window-$window\kb_filter-$filter/marker_info_r2_0.75_1.window-${window}kb.filter-${filter}.txt
    done
  done
done

for filter in 0.25; do
  for window in 1 10 100 1000; do
    # r2 < 0.25
    Rscript scripts/distribution_snps-svs_chrom.R analysis/ld/window-$window\kb_filter-$filter/marker_info_r2_0_0.25.window-${window}kb.filter-${filter}.txt analysis/ld/window-$window\kb_filter-$filter/distribution_markers_chrom_r2_0_0.25.png
    # 0.25 < r2 < 0.5
    Rscript scripts/distribution_snps-svs_chrom.R analysis/ld/window-$window\kb_filter-$filter/marker_info_r2_0.25_0.5.window-${window}kb.filter-${filter}.txt analysis/ld/window-$window\kb_filter-$filter/distribution_markers_chrom_r2_0.25_0.5.png
    # 0.5 < r2 < 0.75
    Rscript scripts/distribution_snps-svs_chrom.R analysis/ld/window-$window\kb_filter-$filter/marker_info_r2_0.5_0.75.window-${window}kb.filter-${filter}.txt analysis/ld/window-$window\kb_filter-$filter/distribution_markers_chrom_r2_0.5_0.75.png
    # r2 > 0.75
    Rscript scripts/distribution_snps-svs_chrom.R analysis/ld/window-$window\kb_filter-$filter/marker_info_r2_0.75_1.window-${window}kb.filter-${filter}.txt analysis/ld/window-$window\kb_filter-$filter/distribution_markers_chrom_r2_0.75_1.png
  done
done
```

> Depletion of SNPs in high LD with SVs in centromeres because there are very few SVs called in these regions. And it's hard to visualize differences in distribution of SNPs in highest LD by different r2 values due to differences in the number of SNPs per quarter.



### Decay

Visualizing LD decay will be import to define which LD window should be used in prediction, i.e. to find SNPs that are close enough to SVs to avoid SNPs in LD due to population structure but also far enough so there's enough recombination happening between SNPs and SVs. Thus, I will plot LD decay of SNP-SV and SNP-SNP for markers up to 100kb and/or 1Mb up or downstream a SV. Additionally, I will compare these LD decay plots with LD decay patterns from SNPs of SNP chip only and from the 282 diversity panel.



#### SNP-SV

```bash
mkdir analysis/ld/decay

# merge ld files from different chr into one
for window in 100 1000; do
  cp analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-${window}kb.filter-0.25.ld analysis/ld/ld_usda_rils_snp-sv_only.window-${window}kb.filter-0.25.ld
  for chr in {2..10}; do
    sed 1d analysis/ld/ld_usda_rils_snp-sv_only.chr${chr}.window-${window}kb.filter-0.25.ld >> analysis/ld/ld_usda_rils_snp-sv_only.window-${window}kb.filter-0.25.ld
  done
  gzip analysis/ld/ld_usda_rils_snp-sv_only.window-${window}kb.filter-0.25.ld
done

# plot ld decay between all SVs and SNPs (distance from SNP to closest SV boundary)
qsub -l walltime=5:00:00,nodes=1:ppn=1,mem=120gb -v IN=analysis/ld/ld_usda_rils_snp-sv_only.window-1000kb.filter-0.25.ld.gz,OUT=analysis/ld/decay/ld_decay_usda_snps-svs_1000kb,OPT="--unequal-windows" scripts/plot_ld_decay.sh

# plot ld decay by sv type
for type in del ins dup inv tra; do
  echo $type
  zcat analysis/ld/ld_usda_rils_snp-sv_only.window-1000kb.filter-0.25.ld.gz | head -n 1 > analysis/ld/ld_usda_rils_snp-sv_only.window-1000kb.filter-0.25.${type}-only.ld
  zgrep $type analysis/ld/ld_usda_rils_snp-sv_only.window-1000kb.filter-0.25.ld.gz >> analysis/ld/ld_usda_rils_snp-sv_only.window-1000kb.filter-0.25.${type}-only.ld
  qsub -l walltime=2:00:00,nodes=1:ppn=1,mem=50gb -v IN=analysis/ld/ld_usda_rils_snp-sv_only.window-1000kb.filter-0.25.${type}-only.ld,OUT=analysis/ld/decay/ld_decay_usda_snps-svs_1000kb_${type}-only,OPT="--unequal-windows" scripts/plot_ld_decay.sh
done

# plot ld decay by chr
for chr in {1..10}; do
  echo $chr
  zcat analysis/ld/ld_usda_rils_snp-sv_only.window-1000kb.filter-0.25.ld.gz | sed 's/^ *//' | tr -s " " | tr " " "\t" | head -n 1 > analysis/ld/ld_usda_rils_snp-sv_only.chr$chr.window-1000kb.filter-0.25.ld
  zcat analysis/ld/ld_usda_rils_snp-sv_only.window-1000kb.filter-0.25.ld.gz | sed 's/^ *//' | tr -s " " | tr " " "\t" | awk -v chr=$chr '$1 == chr' >> analysis/ld/ld_usda_rils_snp-sv_only.chr$chr.window-1000kb.filter-0.25.ld
  qsub -l walltime=2:00:00,nodes=1:ppn=1,mem=50gb -v IN=analysis/ld/ld_usda_rils_snp-sv_only.chr$chr.window-1000kb.filter-0.25.ld,OUT=analysis/ld/decay/ld_decay_usda_snps-svs_1000kb_chr$chr,OPT="--unequal-windows" scripts/plot_ld_decay.sh
done
```



#### SNP-SNP

```bash
# get only SNP-SNP LD
for chr in {1..10}; do
  for filter in 0.25; do
    for window in 100 1000; do
      qsub -v CHR=$chr,WINDOW=$window,FILTER=$filter scripts/select_only_snp-snp_ld.sh
    done
  done
done

# plot ld decay between all SNPs (per chr due to size of files)
for chr in {1..10}; do
  qsub -l walltime=2:00:00,nodes=1:ppn=1,mem=120gb -v IN=/scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-snp_only.chr${chr}.window-100kb.filter-0.25.ld.gz,OUT=analysis/ld/decay/ld_decay_usda_snps-snps_chr${chr}_100kb,OPT="--unequal-windows" scripts/plot_ld_decay.sh
done
```



#### SNP chip

```bash
mkdir -p analysis/snp-chip_ld

# transform hmp to plink format
run_pipeline.pl -Xmx50g -importGuess data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt -export analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4 -exportType Plink

# calculate LD
for filter in 0.25; do
  for window in 100 1000; do
    plink --file analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window ${window}000 --ld-window-kb ${window} --geno ${filter} --out analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.window-${window}kb.filter-${filter}
  done
done

# plot decay
Rscript scripts/plot_ld_decay.R analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.window-1000kb.filter-0.25.ld.gz \
                                analysis/ld/decay/ld_decay_usda_snps-chip-only_1000kb \
                                --unequal-windows
# also plot ld decay by chromosome to compare with reseq SNPs LD decay
for chr in {1..10}; do
  zcat analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.window-1000kb.filter-0.25.ld.gz | head -n 1 > analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.chr$chr.window-1000kb.filter-0.25.ld
  zcat analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.window-1000kb.filter-0.25.ld.gz | awk -v CHR=$chr '$1 == CHR' >> analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.chr$chr.window-1000kb.filter-0.25.ld
  Rscript scripts/plot_ld_decay.R analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.chr$chr.window-1000kb.filter-0.25.ld \
                                  analysis/ld/decay/ld_decay_usda_snps-chip-only_chr${chr}_1000kb \
                                  --unequal-windows
done
```



#### 282 diversity panel

```bash
mkdir -p analysis/div-panel_ld

# download data -- notet that this data has b73v3 coordinates
cd analysis/div-panel_ld
wget https://de.cyverse.org/dl/d/E4CB9C2B-A3A0-4412-9D2B-D2DFC3A06CC9/SNP55K_maize282_AGP3_20190419.hmp.txt.gz
gunzip SNP55K_maize282_AGP3_20190419.hmp.txt.gz
# remove scaffolds
head -n 1 SNP55K_maize282_AGP3_20190419.hmp.txt > SNP55K_maize282_AGP3_20190419_no-scaff.hmp.txt
for chr in {1..10}; do
  awk -v chr="$chr" '$3 == chr' SNP55K_maize282_AGP3_20190419.hmp.txt >> SNP55K_maize282_AGP3_20190419_no-scaff.hmp.txt
done

# transform hmp into plink format
run_pipeline.pl -Xmx10g -importGuess SNP55K_maize282_AGP3_20190419_no-scaff.hmp.txt -export SNP55K_maize282_AGP3_20190419_no-scaff -exportType Plink
for filter in 0.25; do
  for window in 100 1000; do
    # calculate LD
    plink --file SNP55K_maize282_AGP3_20190419_no-scaff.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window ${window}000 --ld-window-kb ${window} --geno ${filter} --out SNP55K_maize282_AGP3_20190419_no-scaff.window-${window}kb.filter-${filter}
  done
done
cd ~/projects/genomic_prediction/simulation

# plot ld decay
Rscript scripts/plot_ld_decay.R analysis/div-panel_ld/SNP55K_maize282_AGP3_20190419_no-scaff.window-1000kb.filter-0.25.ld.gz \
                                analysis/ld/decay/ld_decay_282-div-panel_1000kb \
                                --unequal-windows
```

> LD between SNPs and SVs decays pretty fast. Within the first 100bp the median was r2 ~ 0.5, and decreased to r2 ~ 0.25 around 10kb.



### Closest SNPs to SVs

In addition to visualizing distribution of SNPs in high LD with SVs, it's also important to see how their distribution compare to the closest SNPs to SVs (as they would, in theory, be in high LD with a SV).

```bash
# get closest snps
for filter in 0.25; do
  for chr in {1..10}; do
    for type in all del ins inv dup; do
      qsub -v CHR=$chr,TYPE=$type,FILTER=$filter scripts/keep_closest_SNP_to_SV.sh
    done
  done
done

# merge chromosomes
cd ~/projects/genomic_prediction/simulation/analysis/ld/
for filter in 0.25; do
  for type in all del ins inv dup; do
    cp closest_low-missing-data-SNPs_to_SVs.chr-1.filter-0.25.all.distances.txt closest_low-missing-data-SNPs_to_SVs.filter-0.25.all.distances.txt
    cp usda_rils_projected-SVs-SNPs.chr1.closest-snps-to-svs.filter-${filter}.${type}.hmp.txt usda_rils_projected-SVs-SNPs.closest-snps-to-svs.filter-${filter}.${type}.hmp.txt
    for chr in {2..10}; do
      sed 1d closest_low-missing-data-SNPs_to_SVs.chr-${chr}.filter-0.25.all.distances.txt >> closest_low-missing-data-SNPs_to_SVs.filter-0.25.all.distances.txt
      sed 1d usda_rils_projected-SVs-SNPs.chr${chr}.closest-snps-to-svs.filter-${filter}.${type}.hmp.txt >> usda_rils_projected-SVs-SNPs.closest-snps-to-svs.filter-${filter}.${type}.hmp.txt
    done
  done
done

# plot distribution of distances of closest SNPs to SVs
module load R
cd ~/projects/genomic_prediction/simulation/
Rscript scripts/distribution_distance_closest-snps-to-svs.R analysis/ld/closest_low-missing-data-SNPs_to_SVs.filter-0.25.all.distances.txt \
                                                            analysis/ld/closest_low-missing-data-SNPs_to_SVs.filter-0.25.all.png
# 5230 SNPs were more than 0.1kb away from a SV
# 1601 SNPs were more than 0.5kb away from a SV
# 807 SNPs were more than 1kb away from a SV
# 307 SNPs were more than 5kb away from a SV
# 239 SNPs were more than 10kb away from a SV
# 92 SNPs were more than 100kb away from a SV
# 47 SNPs were more than 1000kb away from a SV

# plot LD distribution for closest snps (all SVs)
for filter in 0.25; do
  for window in 1 10 100 1000; do
    echo "filter $filter - $window kb window"
    Rscript scripts/distribution_ld_snps-svs_closest-snps.R analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-$window\kb.filter-$filter.ld data/usda_SVs_parents.sorted.hmp.txt analysis/ld/usda_rils_projected-SVs-SNPs.closest-snps-to-svs.filter-${filter}.all.hmp.txt analysis/ld/window-${window}kb_filter-${filter}_closest-snps all
  done
done
```


### LD by common parent

Just as comparison, it might be interesting to compare the LD decay for SNP-SNP (chip only) and SNP-SV, and LD distribution of SNPs and SVs for each RILs that have a common parent.



#### SNP chip

```bash
mkdir -p data/hapmap_by_common-parent

# separate hapmap by common parent
for parent in B73 PHG39 PHG47 PH207 PHG35 LH82; do
  echo $parent
  rils=$(head -n 1 data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt | tr "\t" "\n" | cat -n | tr "*" "x" | grep $parent | cut -f 1 | tr -s " " | tr "\n" "," | tr -d " " | sed '$ s/.$/ /')
  cut -f 1-11,$rils data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt > data/hapmap_by_common-parent/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${parent}-parent.hmp.txt
done

for parent in B73 PHG39 PHG47 PH207 PHG35 LH82; do
  # hmp2plk
  run_pipeline.pl -Xmx50g -importGuess data/hapmap_by_common-parent/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${parent}-parent.hmp.txt -export analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${parent}-parent -exportType Plink
  for filter in 0.25; do
    for window in 1000; do
      # calculate LD
      plink --file analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${parent}-parent.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window ${window}000 --ld-window-kb ${window} --geno ${filter} --out analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${parent}-parent.window-${window}kb.filter-${filter}
    done
  done
  # plot ld decay
  Rscript scripts/plot_ld_decay.R analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${parent}-parent.window-1000kb.filter-0.25.ld.gz \
                                  analysis/ld/decay/ld_decay_usda_snps-chip-only_1000kb_${parent}-parent \
                                  --unequal-windows
done
```



#### SNP-SV

```bash
# separate hapmap by common parent
for parent in B73 PHG39 PHG47 PH207 PHG35 LH82; do
  echo $parent
  for chr in {1..10}; do
    rils=$(head -n 1 data/usda_rils_projected-SVs-SNPs.chr$chr.hmp.txt | tr "\t" "\n" | cat -n | tr "*" "x" | grep $parent | cut -f 1 | tr -s " " | tr "\n" "," | tr -d " " | sed '$ s/.$/ /')
    cut -f 1-11,$rils data/usda_rils_projected-SVs-SNPs.chr$chr.hmp.txt > data/hapmap_by_common-parent/usda_rils_projected-SVs-SNPs.${parent}-parent.chr${chr}.hmp.txt
  done
done

# calculate ld in 1Mb windows for ld plot
for parent in B73 PHG39 PHG47 PH207 PHG35 LH82; do
  for chr in {1..10}; do
    qsub -v PARENT=$parent,CHR=$chr scripts/calculate_ld_by_parent.sh
  done
  sleep 60
done

# merge chromosomes
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld/
for parent in B73 PHG39 PHG47 PH207 PHG35 LH82; do
  echo $parent
  cp ld_usda_rils_snp-sv_only.${parent}-parent.chr1.window-1000kb.filter-0.25.ld ld_usda_rils_snp-sv_only.${parent}-parent.window-1000kb.filter-0.25.ld
  for chr in {2..10}; do
    sed 1d ld_usda_rils_snp-sv_only.${parent}-parent.chr${chr}.window-1000kb.filter-0.25.ld >> ld_usda_rils_snp-sv_only.${parent}-parent.window-1000kb.filter-0.25.ld
  done
  gzip ld_usda_rils_snp-sv_only.${parent}-parent.window-1000kb.filter-0.25.ld
done
cd ~/projects/genomic_prediction/simulation

# plot ld decay between all SVs and SNPs (distance from SNP to closest SV boundary)
for parent in B73 PHG39 PHG47 PH207 PHG35 LH82; do
  qsub -l walltime=5:00:00,nodes=1:ppn=1,mem=120gb -v IN=/scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-sv_only.${parent}-parent.window-1000kb.filter-0.25.ld.gz,OUT=analysis/ld/decay/ld_decay_usda_snps-svs_${parent}_1000kb,OPT="--unequal-windows" scripts/plot_ld_decay.sh
done


# for ld distribution plot, calculate LD in 1kb windows
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld/
for parent in B73 PHG39 PHG47 PH207 PHG35 LH82; do
  for chr in {1..10}; do
    # input 'usda_rils_projected-SVs-SNPs.${PARENT}-parent.chr${CHR}.plk*' was generated from 'scripts/calculate_ld_by_parent.sh'
    # calculate LD
    plink --file usda_rils_projected-SVs-SNPs.${parent}-parent.chr${chr}.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window 1000 --ld-window-kb 1 --geno 0.25 --out usda_rils_projected-SVs-SNPs.${parent}-parent.chr${chr}.window-1kb.filter-0.25
    # copy header for filtered file
    zcat usda_rils_projected-SVs-SNPs.${parent}-parent.chr${chr}.window-1kb.filter-0.25.ld.gz | head -n 1 > ld_usda_rils_snp-sv_only.${parent}-parent.chr${chr}.window-1kb.filter-0.25.ld
    # keep only snp and sv r2 (excluding translocations)
    zcat usda_rils_projected-SVs-SNPs.${parent}-parent.chr${chr}.window-1kb.filter-0.25.ld.gz | awk '$3 ~ /^del|^dup|^ins|^inv|^tra/ && $7 !~ /^del|^dup|^ins|^inv|^tra/ || $3 !~ /^del|^dup|^ins|^inv|^tra/ && $7 ~ /^del|^dup|^ins|^inv|^tra/' - >> ld_usda_rils_snp-sv_only.${parent}-parent.chr${chr}.window-1kb.filter-0.25.ld
  done
done
cd ~/projects/genomic_prediction/simulation

# filter plink ld files and plot distribution (only need the first chr file -- it will find the other chr)
B73 PHG39 PHG47 PH207 PHG35 LH82
for parent in B73 PHG39 PHG47 PH207 PHG35 LH82; do
  for filter in 0.25; do
    for window in 1; do
      Rscript scripts/distribution_ld_snps-svs.R /scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-sv_only.${parent}-parent.chr1.window-1kb.filter-0.25.ld data/usda_SVs_parents.sorted.hmp.txt analysis/ld/window-${window}kb_filter-${filter}_by_${parent}
    done
  done
done
```

> LD distribution from RILs with a common parent still shows few SNPs with highest LD to a SV < 0.8, but it's less than those from the LD distribution among all RILs.




### LD by family

Similarly, it might be interesting to compare the LD decay for SNP-SNP (chip only), and LD distribution of SNPs and SVs for each RIL family.



#### SNP chip

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in $crosses; do
  # hmp2plk
  run_pipeline.pl -Xmx50g -importGuess data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross}.hmp.txt -export analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross} -exportType Plink
  for filter in 0.25; do
    for window in 100 1000; do
      # calculate LD
      plink --file analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross}.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window ${window}000 --ld-window-kb ${window} --geno ${filter} --out analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross}.window-${window}kb.filter-${filter}
    done
  done
  # plot ld decay
  Rscript scripts/plot_ld_decay.R analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross}.window-1000kb.filter-0.25.ld.gz \
                                  analysis/ld/decay/ld_decay_usda_snps-chip-only_1000kb_${cross} \
                                  --unequal-windows
done
```



#### SNP-SV

```bash
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

# separate by cross
for cross in $crosses; do
  echo $cross
  for chr in {1..10}; do
    rils=$(head -n 1 data/usda_rils_projected-SVs-SNPs.chr$chr.hmp.txt | tr "\t" "\n" | cat -n | tr "*" "x" | grep $cross | cut -f 1 | tr -s " " | tr "\n" "," | tr -d " " | sed '$ s/.$/ /')
    cut -f 1-11,$rils data/usda_rils_projected-SVs-SNPs.chr$chr.hmp.txt > /scratch.global/della028/hirsch_lab/genomic_prediction/ld/usda_rils_projected-SVs-SNPs.$cross.chr$chr.hmp.txt
  done
done

# calculate ld
for cross in $crosses; do
  for chr in {1..10}; do
    qsub -v CROSS=$cross,CHR=$chr scripts/calculate_ld_by_cross.sh
  done
  sleep 15
done

# merge chromosomes
for cross in $crosses; do
  echo $cross
  cp /scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-sv_only.${cross}.chr1.window-1000kb.filter-0.25.ld /scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-sv_only.${cross}.window-1000kb.filter-0.25.ld
  for chr in {2..10}; do
    sed 1d /scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-sv_only.${cross}.chr${chr}.window-1000kb.filter-0.25.ld >> /scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-sv_only.${cross}.window-1000kb.filter-0.25.ld
  done
  gzip /scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-sv_only.${cross}.window-1000kb.filter-0.25.ld
done

# plot ld decay between all SVs and SNPs (distance from SNP to closest SV boundary)
for cross in $crosses; do
  qsub -l walltime=4:00:00,nodes=1:ppn=1,mem=120gb -v IN=/scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-sv_only.${cross}.window-1000kb.filter-0.25.ld.gz,OUT=analysis/ld/decay/ld_decay_usda_snps-svs_${cross}_1000kb,OPT="--unequal-windows" scripts/plot_ld_decay.sh
done

# for ld distribution plot, calculate LD in 1kb windows
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld/
for cross in $crosses; do
  for chr in {1..10}; do
    # input 'usda_rils_projected-SVs-SNPs.${cross}.chr${chr}.plk*' was generated from 'scripts/calculate_ld_by_parent.sh'
    # calculate LD
    plink --file usda_rils_projected-SVs-SNPs.${cross}.chr${chr}.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window 1000 --ld-window-kb 1 --geno 0.25 --out usda_rils_projected-SVs-SNPs.${cross}.chr${chr}.window-1kb.filter-0.25
    # copy header for filtered file
    zcat usda_rils_projected-SVs-SNPs.${cross}.chr${chr}.window-1kb.filter-0.25.ld.gz | head -n 1 > ld_usda_rils_snp-sv_only.${cross}.chr${chr}.window-1kb.filter-0.25.ld
    # keep only snp and sv r2 (excluding translocations)
    zcat usda_rils_projected-SVs-SNPs.${cross}.chr${chr}.window-1kb.filter-0.25.ld.gz | awk '$3 ~ /^del|^dup|^ins|^inv|^tra/ && $7 !~ /^del|^dup|^ins|^inv|^tra/ || $3 !~ /^del|^dup|^ins|^inv|^tra/ && $7 ~ /^del|^dup|^ins|^inv|^tra/' - >> ld_usda_rils_snp-sv_only.${cross}.chr${chr}.window-1kb.filter-0.25.ld
  done
done
cd ~/projects/genomic_prediction/simulation

# filter plink ld files and plot distribution (only need the first chr file -- it will find the other chr)
for cross in $crosses; do
  echo $cross
  for filter in 0.25; do
    for window in 1; do
      Rscript scripts/distribution_ld_snps-svs.R /scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-sv_only.${cross}.chr1.window-1kb.filter-0.25.ld data/usda_SVs_parents.sorted.hmp.txt analysis/ld/window-${window}kb_filter-${filter}_by_${cross}
    done
  done
done

### BUG: there is something wrong with the above script Rscript when parsing certain chromosomes from crosses
# Error in strsplit(sv, split = ".", fixed = TRUE) : non-character argument
# Calls: add.info.to.ld.df -> apply -> FUN -> unlist -> strsplit
```


### PCA

Another helpful QC is to look at PCA plots and see if SVs can also distinguish maize groups as well as SNPs. I will compare PCA plots built different marker types.

```bash
# only with SVs
qsub scripts/pca_plots_plink_svs.sh

# only with SNPs
qsub scripts/pca_plots_plink_snps.sh

# only with SNPs in high ld (1kb and <0.25 missing data)
qsub scripts/pca_plots_plink_snps-highest-ld.sh

#  only with closest snps (and by marker type)
for type in all del inv dup; do
  qsub -v TYPE=$type scripts/pca_plots_plink_closest-snps_by_type.sh
done
```

> PCA plots between SNPs and SVs are not the same, but I overall there seems to be 3 major distinct groups (and they make sense when looking at the lines within each group) and they seem to agree between SNPs and SVs. Comparing PCA from SVs and SNPs in highest LD, the pattern is the same (just with the x-axis flipped). This pattern, however, is different from the PCA from closest SNPs to SVs (althought the 3 major groups seem to be well separated). In conclusion, these plots show that there is no major problem with our datasets and also reinforce the idea that the closest SNPs to SVs are not necessarily those in highest LD to an SV.



## Trait Simulation
