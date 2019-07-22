# Simulation Genome Prediction

by Rafael Della Coletta, Alex Lipka, Martin Bohn, and Candice Hirsch (June, 2019).

> The objective of this project is to simulate traits in one or more environments and analyze the effects of structural variants in genome prediction accuracy. The overall workflow for the simulations involves simulating traits, running genomic prediction models and validating results.




## Project folder

All data, scripts, analysis, notes and other things related to this project is located on:

```
/Users/rafael/OneDrive/University of Minnesota/PhD/hirsch_lab/projects/genomic_prediction/simulation/
```




## Dataset

The USDA project contains 525 RILs generated from 7 inbred parents (B73, PHJ40, PHG39, PHG47, PH207, PHG35, LH82), which are all ex-PVPs. In addition, 400 F1 hybrids were generated from a partial diallel cross of those RILs.

There is genotypic information (22k SNP chip) for all inbred parents and their RILs. The dataset was provided by Martin Boh from a shared folder on DropBox ("03_DispensibleGenome/EMAMP_SNP_2018"). The `.zip` files were dowloaded and decompressed on May 28, 2019. They were saved on my local directory: `data/SNP_chip/`.

```
Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.csv
Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part1_122318.csv
Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part2_122318.csv
```
> **Note:** there are more genotypes (parents and RILs) in those files besides the ones used in the USDA project.


### HapMap format

Before doing any analysis, it's important to transform the genotypic data that comes from the SNP chip into the HapMap format (see <https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load> for more info about the format). Here's the basic format:

| rs#   | alleles | chrom | pos | strand | assembly  | center   | protLSID | assayLSID | QCcode      | Line 1 | Line 2 | ... |
| ----- | ------- | ----- | --- | ------ | --------- | -------- | -------- | --------- | ----------- | ------ | ------ | --- |
| SNP 1 | A/C     | 1     | 20  | +      | refgen_v4 | MaizeGDB | NA       | NA        | maize_panel | CC     | AA     | ... |

The first 11 columns are required, but not all are required to have information. Usually, only the fields `rs#`, `alleles`, `chrom`, `pos` and the lines' genotypes will be used in downstream analysis (the other fields can be NA).

I used the `scripts/usda_geno2hmp.R` script to transform all parental and RIL data to hapmap format. This script generated these files on `data/` folder:

```bash
# genotypic data on hapmap format
Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.hmp.txt
Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.hmp.txt
# lists relating genotype name with genotype ID
id_table_22kSNPs_DAS_UIUC_ParentalGenotypeData.txt
id_table_22kSNPs_DAS_UIUC_RILsGenotypeData.txt
```

> **Note:** RILs have their genotype ID instead of name because some RILs had multiple IDs. That's why I generated a list relating names to IDs.

Since there are genotypic data for more parents and RILs than used in the USDA project, I need to keep only the genotypes described in the USDA project that will have phenotypic data collected. To do this, Candy sent the file `data/2018_field_planning.xlsx` that contains all the RILs crossesd to generate the 400 hybrids. They can be found under **USDA crosses** on `X10_Nursery_Book` and `X9_Nursery_Book` sheets. I manually copied all this information and saved on file `data/usda_RILs_2018.txt`.

> **Note:** Taking a quick look at this file, I noticed that there are only 328 unique RILs and one hybrid (LH82*PHG47). I'm not sure why there is a F1 cross there, so I need to ask Candy. Also, I was expecting more RILs to generate the hybrids (328 vs 525 RILs). It will be good to talk to Candy about that.

The `scripts/remove_extra_geno-data.R` script

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
cd data/

python ../scripts/create_list_biparental-crosses.py id_table_22kSNPs_DAS_UIUC_RILsGenotypeData.txt usda_biparental-crosses.txt
```

I will use the software [Tassel 5](https://www.maizegenetics.net/tassel) for some basic QC, because this software can generate summaries (including allele frequency) pretty quickly.

Tassel requires that the genotypic data is sorted by position and chromosome number. The genotypic data `usda_22kSNPs_325rils.hmp.txt`has SNPs from chromosome 10 comes after chromosome 1, so it's not sorted correctly. Thus, I used the `Data > Sort Genotype File` plugin on Tassel (GUI version) to quickly sort the genotypic data and create another hapmap file called `usda_22kSNPs_325rils.sorted.hmp.txt`. This file will be used for the rest of QC analysis.

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


The GUI version of this script would be this:

* `File > Open` and select `usda_22kSNPs_325rils.sorted.hmp.txt`;
* `Filter > Filter Genotype Table Taxa` and add list of RILs from a biparental cross;
* Select filtered file in the left panel (under `Sequence` folder) and summarize data by clicking in `Data > Geno Summary` (the table with allele frequencies will have the suffix `_SiteSummary`)
* Export results by selecting the table of interest in the left panel (under `Result/Genotype Summary` folder) and then `Results > Table > Export (Tab)`;
* Repeat for other biparental crosses.


**Recombination frequency**

Before start working on getting the recombination frequencies for each biparental cross, I had to sort the genotypic data for parents and RILs and save them as diploid format (genotypes as "AA, AB, BB" instead of "A, H, B"). To format the parents hapmap file (`usda_22kSNPs_7parents.hmp.txt`), I used Tassel 5 GUI version by clicking on `Data > Sort Genotype File` and then `File > Save As... > Hapmap Diploid`. The genotpic file with the RILs was already sorted on previous section (`usda_22kSNPs_325rils.sorted.hmp.txt`), so I just opened this sorted version in tassel and `File > Save As... > Hapmap Diploid`.

The two files generated that will be used in this analysis are:

```
data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt
data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt
```

To estimate recombination frequency for each biparental population of all 325 RILs, I wrote the `scripts/recomb-freq_biparental-crosses.R` script. Specifically, it does three things:

1. Filter the sorted diploid hapmap files to have genotypic information only from parents and RILs that make up a specific biparental cross. It takes information from the file `data/usda_biparental-crosses.txt` to do that, and output files on `data/biparental-crosses` (two files per cross: one hapmap for the parents, and another for RILs).

  > **Note:** the `alleles` column of the filtered hapmap files are not correct, because it shows the alleles present in all parents and RILs (not only for that particular biparental cross). This might be corrected in future versions of the script, if needed.

2. Combine the hapmap files for each biparental population and converts them into the correct input format required by [rqtl](http://www.rqtl.org/), the R package that estimates the allele frequency.

3. Run rqtl to estimate recombination frequencies in each biparental population, and plot the genetic distance (cM) by physical distance (Mb) for each cross as well. In fact, here is the complete list of things the function `EstimateRecombinationFreq()` does:

  * Read biparental cross information (genotypic data formated to rqtl).

  * Remove missing data (remove individual if half of the markers are missing, remove markers missing in half of individuals).

    > **Note:** these cut-offs were arbitrary.

  * Remove duplicated individuals.

  * Remove duplicated markers.

  * Remove markers with segregation distortion (FDR adjusted p-value < 0.05).

  * Estimate recombination frequency.

  * **Plot** genotype frequencies by individual, and **write** genotypic data per RIL after all filtering.

  * **Plot** genetic and physical positions of markers per chromosome.

  * **Write summaries** of each cross before and after filtering.

  * **Write table** with genetic and physical positions of markers.

    > **Note:** all ouput goes to the the respective **cross folder** in `analysis/qc/`. It also generates a log file called `rqtl_log.txt`.


**Plot allele frequencies of filtered markers**

I also wrote the file `usda_allele-freq_dist.R`. This script extracts the marker names from the recombination frequency tables for each cross (`analysis/qc/{cross}/recomb-freq_{cross}_rils.txt`), which were filtered to be polymorphic between parents and not show segregation distortion, to filter Tassel's `_SiteSummary` tables. This filtered table was used to make plots of the distributions of allele frequencies for each biparental population. The output was stored their respective **cross folder** in `analysis/qc/`.

> **Note:** Some of my first plots of the distributions had a weird shape. This was solved when I increased `binwidth` of histograms from `0.05` to `0.07`.

Then, as Candy suggested, I checked if markers that have allele frequency < 0.25 or > 0.75 were the same across populations. The majority of such markers (1,038) are unique of a population, while only 81 of these markers have extreme allele frequencies in more than one population.


**Karyotype of RILs**

Chrm info extracted from: `GCA_000005005.6_B73_RefGen_v4_assembly_stats.txt` of B73_RefGen_v4 assembly at <https://www.maizegdb.org/genome/genome_assembly/Zm-B73-REFERENCE-GRAMENE-4.0> or <ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4> (accessed 07/19/2019).

Centromeres (extracted from Table S2 of Schneider et al 2016, <https://doi.org/10.1073/pnas.1522008113>) and stored at `data/centromeres_Schneider-2016-pnas_v2-v3.bed` (coordinates of chromosomes 1 to 9 are from AGPv2 assembly, and chromosome 10 is from AGPv3). Used Assembly Converter from Gramene <http://ensembl.gramene.org/Zea_mays/Tools/AssemblyConverter?db=core> with the data above to convert coordinates from v2 to v4 (`data/centromeres_v2-to-v4.bed`) and v3 to v4 (`data/centromeres_v3-to-v4.bed`).

As part of the R script `script/plot_ril_karyotypes.R`, I filtered the files above to get the centromeres coordinates in the v4 assembly, which can be seen in `data/centromeres_Schneider-2016-pnas_v4.bed` and in the table below:

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

After looking at the karyotypes from RIL 2 (genotype frequency ~50/50) and RIL 3 (genotype frequency ~25/75) of cross B73 x LH82, I was not able to see any striking difference or something weird that can be causing the deviation on allele frequency. I will show them to Candy, to see if she picks up anything.

**QC conclusion**

The above results indicates that there is not apparent big problem in the dataset that will be used for simulations.




## Testing simulations

I created a new folder called `tests` to store files from the tests I run while adjusting the main simulation script. I also added another markdown document called `notes/testing_simulations.md` to keep track of the changes I make.




## Projecting SVs from parents to RILs
