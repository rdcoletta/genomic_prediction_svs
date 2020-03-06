# Structural variation and genomic prediction: a simulation approach

by Rafael Della Coletta, Alex Lipka, Martin Bohn, and Candice Hirsch (June 2019 - February 2020).

> The objective of this project is to simulate traits in one or more environments and analyze the effects of structural variants in genome prediction accuracy. The overall workflow for the simulations involves simulating traits, running genomic prediction models and validating results.




## Project folder

All data, scripts, analyses, notes and other things related to this project is located on my MSI account:

```bash
cd /home/hirschc1/della028/projects/genomic_prediction/simulation/
mkdir {analysis,data,scripts}
```




## Dataset

The USDA project contains 525 RILs generated from 7 inbred parents (B73, PHJ40, PHG39, PHG47, PH207, PHG35, LH82), which are all ex-PVPs. In addition, 400 F1 hybrids were generated from a partial diallel cross of those RILs.

There is genotypic information (22k SNP chip) for all inbred parents and their RILs. The dataset was provided by Martin Bohn from a shared folder on DropBox ("03_DispensibleGenome/EMAMP_SNP_2018"). The `.zip` files were dowloaded and decompressed on May 28, 2019. They were then transferred to `data/SNP_chip/` via FileZilla.

```
Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.csv
Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part1_122318.csv
Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part2_122318.csv
```
> There are more genotypes (parents and RILs) in those files besides the ones used in the USDA project.


### HapMap format

Before doing any analysis, it's important to transform the genotypic data that comes from the SNP chip into the HapMap format (see <https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load> for more info about the format). Here's the basic format:

| rs#   | alleles | chrom | pos | strand | assembly  | center   | protLSID | assayLSID | QCcode      | Line 1 | Line 2 | ... |
| ----- | ------- | ----- | --- | ------ | --------- | -------- | -------- | --------- | ----------- | ------ | ------ | --- |
| SNP 1 | A/C     | 1     | 20  | +      | refgen_v4 | MaizeGDB | NA       | NA        | maize_panel | CC     | AA     | ... |

The first 11 columns are required, but not all are required to have information. Usually, only the fields `rs#`, `alleles`, `chrom`, `pos` and the lines' genotypes will be used in downstream analysis (the other fields can be NA).

I used the `scripts/usda_geno2hmp.R` script to transform all parental and RIL data to hapmap format:

```bash
Rscript scripts/usda_geno2hmp.R data/SNP_chip/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.csv \
                                data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part1_122318.csv \
                                data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part2_122318.csv \
                                data
```

This script generated these files on `data/` folder:
* Genotypic data on hapmap format:
  ```
  Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.hmp.txt
  Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.hmp.txt
  ```
* Lists relating genotype name with genotype ID:
  ```
  id_table_Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.txt
  id_table_Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.txt
  ```

> RILs have their genotype ID instead of name because some RILs had multiple IDs. That's why I generated a list relating names to IDs. Also, I converted `-` to `N` to represent missing data for compatibility with TASSEL 5 software.

Since there are genotypic data for more parents and RILs than used in the USDA project, I need to keep only the genotypes described in the USDA project that will have phenotypic data collected. To do this, Candy sent the file `data/2018_field_planning.xlsx` that contains all the RILs crossesd to generate the 400 hybrids. They can be found under **USDA crosses** on `X10_Nursery_Book` and `X9_Nursery_Book` sheets. I manually copied all this information and saved on file `data/usda_RILs_2018.txt`.

> Taking a quick look at this file, I noticed that there are only 328 unique RILs and one hybrid (LH82*PHG47). I'm not sure why there is a F1 cross there, so I need to ask Candy. Also, I was expecting more RILs to generate the hybrids (328 vs 525 RILs). It will be good to talk to Candy about that.

The `scripts/remove_extra_geno-data.R` script removes extra genotypic data that will not be used in the USDA project from hapmap files:

```bash
Rscript scripts/remove_extra_geno-data.R data/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.hmp.txt \
                                         data/usda_22kSNPs_7parents.hmp.txt \
                                         data/Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.hmp.txt \
                                         data/usda_22kSNPs_325rils.hmp.txt \
                                         data/usda_RILs_2018.txt
```

This file generated the following files on `data/` folder, and will be used in further downstream analysis:

```bash
# parental data
usda_22kSNPs_7parents.hmp.txt
# RIL data (325 out of 328 RILs were genotyped)
usda_22kSNPs_325rils.hmp.txt
```


### RIL data QC

**Allele frequency of each marker**

With RILs from a biparental cross, we have the expectation that each locus will have only two alleles with frequency ~50% each (if allele is not fixed; then it will be 100% for major allele and 0% for minor). Large deviations from this expectation may indicate some errors during genotyping and the locus should be removed from analysis.

But first, I have to find out which RILs belong to the same biparental cross. I wrote the `scripts/create_list_biparental-crosses.py` script to generate a table with a biparental cross in a column and all the RIL IDs in another, and ran this in the command line to create the table:

```bash
python scripts/create_list_biparental-crosses.py data/id_table_Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.txt \
                                                 data/usda_biparental-crosses.txt
```

I will use the software [Tassel 5](https://www.maizegenetics.net/tassel) for some basic QC, because this software can generate summaries (including allele frequency) pretty quickly.

Tassel requires that the genotypic data is sorted by position and chromosome number. The genotypic data `data/usda_22kSNPs_325rils.hmp.txt`has SNPs from chromosome 10 comes after chromosome 1, so it's not sorted correctly. Thus, I used `SortGenotypeFilePlugin` from Tassel to quickly sort the genotypic data and create another hapmap file called `data/usda_22kSNPs_325rils.sorted.hmp.txt`. Then I ran Tassel again to transform the sorted data into diploid hapmap format. This file will be used for the rest of QC analysis. I also did the same procedure with the parental data.

```bash
# sort ril data
run_pipeline.pl -Xmx2g -SortGenotypeFilePlugin \
                       -inputFile data/usda_22kSNPs_325rils.hmp.txt \
                       -outputFile data/usda_22kSNPs_325rils.sorted.hmp.txt \
                       -fileType Hapmap
run_pipeline.pl -Xmx2g -importGuess data/usda_22kSNPs_325rils.sorted.hmp.txt \
                       -export data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt \
                       -exportType HapmapDiploid

# sort parents
run_pipeline.pl -Xmx2g -SortGenotypeFilePlugin \
                       -inputFile data/usda_22kSNPs_7parents.hmp.txt \
                       -outputFile data/usda_22kSNPs_7parents.sorted.hmp.txt \
                       -fileType Hapmap
run_pipeline.pl -Xmx2g -importGuess data/usda_22kSNPs_7parents.sorted.hmp.txt \
                       -export data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt \
                       -exportType HapmapDiploid
```

Since I had 59 biparental crosses in `data/usda_biparental-crosses.txt`, it would be very annoying to perform the analysis one-by-one. Thus, I decided to use the command-line version of Tassel by using the following commands:

```bash
# create directory to store QC analyses
mkdir -p analysis/qc

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
    run_pipeline.pl -Xmx6g -importGuess data/usda_22kSNPs_325rils.sorted.hmp.txt \
                    -FilterTaxaBuilderPlugin -taxaList $ril_list -endPlugin \
                    -GenotypeSummaryPlugin -endPlugin \
                    -export analysis/qc/$cross/$cross\_OverallSummary,analysis/qc/$cross/$cross\_AlleleSummary,analysis/qc/$cross/$cross\_SiteSummary,analysis/qc/$cross/$cross\_TaxaSummary
  done
} < "data/usda_biparental-crosses.txt" > "analysis/qc/tassel_log.txt"

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir analysis/qc/*
```

> **Note**: this bash script also outputs a log file at `analysis/qc/tassel_log.txt`


**Recombination frequency**

I will use the sorted diploid hapmap files from parents (`data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt`) and RILs (`data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt`) to estimate recombination frequency for each biparental population. For this purpose, I wrote the `scripts/recomb-freq_biparental-crosses.R` script. Specifically, it does three things:

1. Filter the sorted diploid hapmap files to have genotypic information only from parents and RILs that make up a specific biparental cross. It takes information from the file `data/usda_biparental-crosses.txt` to do that, and output files on `data/biparental-crosses` (two files per cross: one hapmap for the parents, and another for RILs).

  > The `alleles` column of the filtered hapmap files are not correct, because it shows the alleles present in all parents and RILs (not only for that particular biparental cross). This might be corrected in future versions of the script, if needed.

2. Combine the hapmap files for each biparental population and converts them into the correct input format required by [rqtl](http://www.rqtl.org/), the R package that estimates the allele frequency.

3. Run rqtl to estimate recombination frequencies in each biparental population, and plot the genetic distance (cM) by physical distance (Mb) for each cross as well. In fact, here is the complete list of things the function `EstimateRecombinationFreq()` does:

  * Read biparental cross information (genotypic data formated to rqtl).

  * Remove missing data (remove individual if half of the markers are missing, remove markers missing in half of individuals).

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
Rscript scripts/recomb-freq_biparental-crosses.R data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt \
                                                 data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt \
                                                 data/usda_biparental-crosses.txt \
                                                 analysis/qc \
                                                 data/biparental-crosses
```


**Plot allele frequencies of filtered markers**

I also wrote the file `scripts/usda_allele-freq_dist.R`. This script extracts the marker names from the recombination frequency tables for each cross, which were filtered to be polymorphic between parents and not show segregation distortion, to filter Tassel's `_SiteSummary` tables. This filtered table was used to make plots of the distributions of allele frequencies for each biparental population. The output was stored their respective **cross folder** in `analysis/qc`.

> Some of my first plots of the distributions had a weird shape. This was solved when I increased `binwidth` of histograms from `0.05` to `0.07`.

Then, as Candy suggested, I checked if markers that have allele frequency < 0.25 or > 0.75 were the same across populations. The majority of such markers (910) are unique of a population, while only 65 of these markers have extreme allele frequencies in more than one population.

```bash
Rscript scripts/usda_allele-freq_dist.R analysis/qc
```

**Karyotype of RILs**

To plot karyotypes I need information from chromosome length and centromere positions. The chromosome info was manually extracted from `GCA_000005005.6_B73_RefGen_v4_assembly_stats.txt` of B73_RefGen_v4 assembly at [MaizeGDB](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4) (accessed 07/19/2019) and saved in a tab-delimited text file (`data/B73_RefGen_V4_chrm_info.txt`).

Centromeres positions were extracted from Table S2 of [Schneider et al 2016](https://doi.org/10.1073/pnas.1522008113) and stored at `data/centromeres_Schneider-2016-pnas_v2-v3.bed`. However, coordinates of centromeres from chromosomes 1 to 9 are from AGPv2 assembly, and chromosome 10 is from AGPv3. Thus, I used the [Assembly Converter](http://ensembl.gramene.org/Zea_mays/Tools/AssemblyConverter?db=core) from Gramene  with the data above to convert coordinates from v2 to v4 (`data/centromeres_v2-to-v4.bed`) and v3 to v4 (`data/centromeres_v3-to-v4.bed`). Then, I wrote the R script `script/get_centromeres_v4_coord.R` and filtered the previous files to get the centromeres coordinates in the v4 assembly, which can be seen in `data/centromeres_Schneider-2016-pnas_v4.bed` and in the table below:

```bash
Rscript scripts/get_centromeres_v4_coord.R data/centromeres_v2-to-v4.bed data/centromeres_v3-to-v4.bed data/centromeres_Schneider-2016-pnas_v4.bed
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


After that, I ran `scripts/plot_ril_karyotypes.R`to see how the parental haplotypes are distributed in 5 randomnly selected RILs per population. For example, after looking at the karyotypes from RIL 2 (genotype frequency ~50/50) and RIL 3 (genotype frequency ~25/75) of cross B73 x LH82, I was not able to see any striking difference or something weird that can be causing the deviation on allele frequency. Overall, no major issues were observed in other "karyoplots".

```bash
Rscript scripts/plot_ril_karyotypes.R data/B73_RefGen_V4_chrm_info.txt \
                                      data/centromeres_Schneider-2016-pnas_v4.bed \
                                      analysis/qc
```

> It's important to check the number of recombinations per cromosome. At each generation you advance in RILs, you expect about 1-2 recombination events (although after later generations, F5-6, you should expect less because then you start recommbining blocks that have the same genotypes). For those RILs that have extreme genotypes, you should expect to see recombination more towards the end of chromosomes, that's why you would see a overrepressentation of one genotype (recombination is not always in the same place, it's a distribution).


**Percent heterozygosity per RIL population**

I wrote `scripts/markers_summary.R` to check number of missing, homozygous, and heterozygous genotypes in RILs (and parents). This script takes information from `data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt` and `data/usda_biparental-crosses.txt` to calculate these summaries. The files generated were the table `summary_markers_325rils.txt` and the boxplots `summary_markers_325rils_missing-per-cross.png` and `summary_markers_325rils_hets-per-cross.png` at the `analysis/qc` folder.

```bash
Rscript scripts/markers_summary.R data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt \
                                  analysis/qc/summary_markers_7parents.txt \
                                  data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt \
                                  analysis/qc/summary_markers_325rils.txt \
                                  data/usda_biparental-crosses.txt
```

By looking at the boxplot of missing data, what stands out is that cross `B73xPHG47` had the highest variation and amount of missing data than other populations, and few populations showed very few RILs with more extreme missing data. Although the maximum amount of missing data for a RIL was about 9%, the median amount of missing data per cross was below 2.5%.

For heterozygous markers, they were also mostly between 0 and 3% (which would be the expected for an F6 population), however population `B73xPHG47` also displayed a lot of variation and higher values of heterozygous markers. Very few RILs displayed high levels of heterozygosity, but those values were much higher than the expected (between 10 to 30%).

> I wrote this script to double check TASSEL's summary to generate plots and to do the same with parents genotypes (see below). In fact, the results are the same from the ones produced by TASSEL when using the `GenotypeSummary` plugin, but are now displayed in only one file.


**Percent heterozygosity per parent**

The script `scripts/markers_summary.R` also checks the number of missing, homozygous and heterozygous genotypes in RILs that have the same parent. This can help identify any problems with a parent during generation of RILs. The boxplots `summary_markers_325rils_missing-per-parent.png` and `summary_markers_325rils_het-per-parent.png` were generated and saved at the `analysis/qc` folder. There is nothing that stands out in these plots, so I cannot see any parent that contributes with more heterozygosity or missing data.

> Quick note that RILs for the parent PHJ40 were not genotyped (hence the 6 parents in the plots, instead of 7).


### Parental data QC

It's also important to make sure the parental lines also have good quality data, especially because parental genotypic data will be used for projection of SVs into RILs. The `scripts/markers_summary.R` also summarizes the number of missing, homozygous and heterozygous markers for the parents. It uses information from `data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt`. The files generated were the table `summary_markers_7parents.txt` and the boxplots `summary_markers_7parents_missing.png` and `summary_markers_7parents_het.png` at the `analysis/qc` folder. By looking at the bar plots, it's very clear that **PHG35 has higher missing and het markers than other parents** (~4% compared to ~0.3% of other parents).

Another way to inspect parents is by plotting the karyotypes with homo, het, and missing markers to see if there are some regions of the genome with higher proportion of hets and missing data. To do that, I wrote `scripts/plot_parents_karyotypes.R` based on my previous script to plot karyotypes for RILs. The karyotypes were saved in `analysis/qc/parents`.

```bash
Rscript scripts/plot_parents_karyotypes.R data/B73_RefGen_V4_chrm_info.txt \
                                          data/centromeres_Schneider-2016-pnas_v4.bed \
                                          data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt \
                                          analysis/qc/parents
```

Looking at the karyotypes, the distribution of both missing and het markers seem to occur more further away from the centromeres but in a bit random way. Except for PHG35, I don't see any big clusters of missing and/or het markers. Since PHG35 had much more missing and het data than other parents, it's easier to see some clusters of such markers.

After talking to Candy, it's likely that the PHG35 souce used to be genotyped had issues (either some backcross or cross-contamination). Thus, in order to correctly project the SVs from this parent to its RILs, we need to make sure the genotype from PHG35 agrees with the RILs. To do that, we will reconstruct the PHG35 genome based on the markers available from the RILs.
