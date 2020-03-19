# Structural variation and genomic prediction: a simulation approach

by Rafael Della Coletta, Alex Lipka, Martin Bohn, and Candice Hirsch (June 2019 - February 2020).

> The objective of this project is to simulate traits in one or more environments and analyze the effects of structural variants in genome prediction accuracy. The overall workflow for the simulations involves simulating traits, running genomic prediction models and validating results.




## Project folder

All data, scripts, analyses, notes and other things related to this project is located on my MSI account:

```bash
cd /home/hirschc1/della028/projects/genomic_prediction/simulation/
mkdir {analysis,data,scripts}
```




## SNP dataset

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




## SV dataset

On August 9, 2019, Patrick Monnahan sent me five `.vcf` files containing structural variation calls from the software Lumpy. Each file is a SV call of 100 lines against one of the following reference genomes: B73, Mo17, PH207, PHB47, or W22. **This results are preliminary** and there are a lot of false-positives in there. However, they will be useful for me to select the 7 lines I need, project the SVs from parents to RILs, and incorporate these SV data into my simulation scripts.

The files are located at `data/SV_calls`, and I will use only the SVs called against the B73 reference genome for now.

```bash
# decompress vcf file of SV calls to B73
gunzip data/SV_calls/B73v4_2019-08-09.ls.RT.vcf.gz
```

### Hapmap format

I wrote the python script `scripts/vcf2hapmap.py` which extract information about structural variants and transform into a hapmap format file sorted by chromosome and positions. Since the `.vcf` file that Patrick sent me contains 100 inbred lines, I have to select only the 7 inbred parents used in the USDA project. Running this script will generate the file `data/usda_SVs_7parents.sorted.hmp.txt`:

```bash
# for help on how to use this script
python scripts/vcf2hapmap.py --help
# run the script
python scripts/vcf2hapmap.py data/SV_calls/B73v4_2019-08-09.ls.RT.vcf \
                             data/usda_SVs_7parents.sorted.hmp.txt \
                             B73,LH82,PH207,PHG35,PHG39,PHG47,PHJ40
```

The type of SV will be displayed in the marker ID (`del.[ID]` for deletions, `dup.[ID]` for duplications, etc.). Each line will have either a value of `AA` if SV is `A`bsent, or `TT` if SV is `T`here. I had to do that because I will use Tassel to project genotypes from parents to RILs, and it doesn't accept anything other than nucleotides (thus, I couldn't use numbers to indicate presence/absence of SV, as originally thought). Missing data will be coded as `NN`. Also, since SVs spam hundreds (or thousands) of bp and the exact breakpoints are hard to call, the position indicated in the hapmap file wil be the middle point of the SV.

> After projection of parental genotypes into RILs, the hapmap will be transformed into the numeric format for genomic prediction simulations. Thus, for SVs, `AA` will be transformed into `0` and `TT` to `2`.




## Check reference genome version of SNP chip

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
head -n 1 data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt > data/check_refgen_SNPchip/tmp.hmp.txt
# keep only markers in the first chromosome
awk '$3 == 1' data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt >> data/check_refgen_SNPchip/tmp.hmp.txt
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




## Convert SNP coordinates from B73v2 to B73v4

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
python scripts/extract_probe_seqs.py data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt data/refgen/B73v2/ 100 data/probes-100bp_22kSNPchip_B73v2.fa
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
| Processed               | 20,139 |
| Uniquely mapped         | 19,850 |
| Multimapped (supressed) | 72     |
| Unmapped                | 217    |

Now, I need to correct the positions of the SNPs in the hapmap files `data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt` and `data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt` based on new coordinates from SAM file. I wrote a python script called `scripts/convert_hmp_v2-to-v4.py` to do that and create the files `data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt` and `data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt`, and another file with SNP coordinates in both v2 and v4. This script also removes reads that uniquely mapped with no mismatch to a different chromosome in v4, since it's very unlikely that there would be such major differences between v2 and v4 assemblies (i.e., these are probably mismapped).

```bash
python scripts/convert_hmp_v2-to-v4.py data/probes-100bp_22kSNPchip_aligned-to-B73v4.sam \
                                       data/SNP_positions_v2-to-v4_probes-100bp.txt \
                                       data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt,data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt
# 23 reads discarded due to mapping into different chromosome in v4
```

I performed a quick QC (using only chromosome 1 data) to see if SNPs from v4 were actually the same as in v2 by running `scripts/check_refgen_SNPchip_v2tov4_probes.R`. The results showed that some markers that didn't match between v2 and v4, and that these mismatches were due to mapping to complementary strand:

```bash
Rscript scripts/check_refgen_SNPchip_v2tov4_probes.R data/check_refgen_SNPchip/ref-markers_chr1.txt data/SNP_positions_v2-to-v4_probes-100bp.txt
# 217 marker alleles differ between v2 and v4
# Converting alleles that differ to complementary base...
# 0 marker alleles differ between v2 and v4
```

Thus, I corrected the strand of those markers in the hapmap files `data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt` and `data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt` by running `scripts/correct_SNP_strands.R`. Finally, I corrected the alleles' column of these hapmap files using Tassel.

```bash
# correct markers not phased
Rscript scripts/correct_SNP_strands.R data/SNP_positions_v2-to-v4_probes-100bp.txt \
                                      data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt \
                                      data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt

# correct alleles column
run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt \
                        -export data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt \
                        -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt \
                        -export data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt \
                        -exportType HapmapDiploid
```




## Reconstruct PHG35 based on RIL genotypes

Based on previous QC, the PHG35 parent has much more heterozygotes than expected for a fully inbred line. This might very likely be due to polen contamination in the seeds used to do the genotyping. Now that we have SNP chip data in the refgen v4 assembly, I can use the marker genotypes in the RIL data from all PHG35 progeny and figure out the PHG35 haplotype based on the genotype of the other parent used to develop that progeny. To do that, I wrote `scripts/reconstruct_PHG35_from_RIL_data.R`.

```bash
Rscript scripts/reconstruct_PHG35_from_RIL_data.R data/usda_biparental-crosses.txt \
                                                  data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt \
                                                  data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt
```

To make sure the reconstruction worked, I ran again `scripts/markers_summary.R`. Plots at `analysis/qc` folder show that, despite higher mising data relatively to other parents, the number of heterozygotes is now zero.

```bash
Rscript scripts/markers_summary.R data/usda_22kSNPs_7parents.sorted.diploid.v4.PHG35-reconstructed.hmp.txt \
                                  analysis/qc/summary_markers_7parents_PHG35-reconstructed.txt \
                                  data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt \
                                  analysis/qc/summary_markers_325rils_v4.txt \
                                  data/usda_biparental-crosses.txt
```




## Remove SNPs inside deletions and merge SNP and SV data

For projections, I need to create parental and RIL files with all SNPs and SVs. However, before that, I need to address a particular issue when using these two types of variation: SNPs that are found inside deletions. Such SNPs are problematic because they will have segregation issues when you compare multiple lines that have or does not have that SV.

The `scripts/merge_SNPs_SVs_hapmaps.R` will remove SNPs that are within 100kb of a deletion first and then it will combine marker data for parents and RILs separately. I set up 100 kb threshold because ~99% of the deletions are smaller than this and the very large deletions (>1 Mb) would make me lose a lot of SNPs. Importantly, each family will be filtered separately.

```bash
Rscript scripts/merge_SNPs_SVs_hapmaps.R data/usda_22kSNPs_7parents.sorted.diploid.v4.PHG35-reconstructed.hmp.txt \
                                         data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt \
                                         data/usda_SVs_7parents.sorted.hmp.txt \
                                         data/usda_SNPs-SVs_7parents.not-in-PAVs.hmp.txt \
                                         data/usda_SNPs-SVs_325rils.not-in-SVs.hmp.txt \
                                         data/usda_biparental-crosses.txt \
                                         data/merged_hapmaps_by_cross
```

As shown in the table below, only very few SNPs were inside the boundaries of a deletion and thus removed:

| chr | total SNPs | SNPs removed |
| --- | ---------- | ------------ |
| 1   | 3153       | 6            |
| 2   | 2401       | 2            |
| 3   | 2560       | 2            |
| 4   | 2079       | 0            |
| 5   | 2134       | 2            |
| 6   | 1764       | 0            |
| 7   | 1602       | 0            |
| 8   | 1299       | 2            |
| 9   | 1524       | 2            |
| 10  | 1311       | 0            |

The output from this script is a merged set of SNPs and SVs not in deletions for all parents (`data/usda_SNPs-SVs_7parents.not-in-PAVs.hmp.txt`) and all RILs (`data/usda_SNPs-SVs_325rils.not-in-SVs.hmp.txt`), and also merged sets for each RIL family (files in `data/merged_hapmaps_by_cross`). The later sets will be the ones used for projection. It's worth noting that genotypic information of SVs is only available on parental dataset. The RIL dataset have the SV names but not the genotypes (they are all `NN`).




## Project SVs from parents to RILs

Projection of SVs need to be done on a family basis. The merged SNP and SV data from parents will serve as donors for projections by [FILLIN](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/FILLIN/FILLIN). The haplotype blocks generated for each parent will then be used to determine the haplotypes from the RIL dataset.

After preliminary analysis, using the options `-hapSize 2000` (haplotype block size), `-minTaxa 1` (create haplotype for each parent), and `-hybNN false` (avoid combining haplotypes) yielded the best projection results.


```bash
# create directory to store results of imputation
mkdir -p analysis/projection

# change directory
cd data/

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
    # also change "*" by "x" in the cross name to match the file names
    cross=$(echo $line | cut -f1 | tr "*" "x")
    # create haplotypes from parents
    run_pipeline.pl -FILLINFindHaplotypesPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_parents.sorted.hmp.txt \
                    -o ../analysis/projection/donors_$cross \
                    -hapSize 2000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection/donors_$cross \
                    -o ../analysis/projection/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 2000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection/*


# transform to diploid hapmap
cd ../analysis/projection
for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done

# return to project's home directory
cd ~/projects/genomic_prediction/simulation

```

To summarize the results of projections for each family, I wrote `scripts/count_projected_SVs.R`. The average **SV projection was 93.31%** with average accuracy of 92.84%. Details about projections for each family were written into `analysis/projection/projection_summary.txt`, but can also be visualized in different plots saved into `analysis/projection`.

```bash
Rscript scripts/count_projected_SVs.R data/usda_SVs_7parents.sorted.hmp.txt \
                                      data/usda_biparental-crosses.txt \
                                      analysis/projection
```

Another QC is to plot the karyotypes from RILs and compare with their parents to check whether the haplotype blocks agree. After running `scripts/plot_ril_karyotype_SVs.R`, I don't see any major disagreement between parental and RIL data, meaning that projections seem to be very accurate indeed.

```bash
Rscript scripts/plot_ril_karyotypes_SVs.R data/B73_RefGen_V4_chrm_info.txt \
                                          data/centromeres_Schneider-2016-pnas_v4.bed \
                                          data/usda_biparental-crosses.txt \
                                          analysis/projection \
                                          data/merged_hapmaps_by_cross \
                                          analysis/qc/karyotypes \
                                          --rils=random
```

In order to use these projections in genomic prediction simulations, I have to merge each projected hapmap into one file. For this purpose, I wrote `scripts/merge_projected_crosses.R` to merge projected files into the file `data/usda_SNPs-SVs_325rils.not-in-SVs.projected.hmp.txt`.

```bash
Rscript scripts/merge_projected_crosses.R data/usda_SNPs-SVs_325rils.not-in-SVs.hmp.txt \
                                          analysis/projection \
                                          data/usda_biparental-crosses.txt
```




## Resequencing SNPs

After some preliminary work on genomic prediction models with simulated data, we noticed that we end up with far more SVs (~10k) than SNPs in LD to an SV (~3k). This very likely happened because the SNP chip doesn't have the SNP density required to fully capture every SNP in LD to an SV. Fortunately, there is resequencing data available for all parents of the USDA population, and this type of data can provide the resolution needed. However, since there is no resequencing data for the RILs, I have to project these resequencing SNPs from USDA parents to respective RILs.

All parents for this population were resequenced by [Mazaheri et al (BMC Plant Biology, 2019)](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-019-1653-x), and the data can be retrieved in [this Dryad repository](https://datadryad.org/stash/landing/show?big=showme&id=doi%3A10.5061%2Fdryad.n0m260p).

> Note that the real link for for downloading this data will be sent by email.

The main file to use for this project is called `62biomAP_v_B73_SNPMatrix.txt`, which is a matrix containing all SNP calls for 56 maize inbred lines. Thus, after downloading the data, I have to select only the parents of the USDA population, transform the data to the hapmap format and then proceed to the rest of the analysis.

Unzip `62biomAP_v_B73_SNPMatrix.txt.txt.tar.gz` and it's has all the SNPs from 56 inbreds. I will just select the SNPs from the 7 parents of my population (B73, LH82, PH207, PHG35, PHG39, PHG47, and PHJ40). The dataset to be used is called `data/reseq_snps/biomAP_v_B73_SNPMatrix_7parents.txt`.

```bash
mkdir -p data/reseq_snps
cd data/reseq_snps

# note that if I download this data again, the name of the file will be different
wget http://merritt.cdlib.org/cloudcontainer/mrtstore2/21061229.tar.gz
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
cut -f 1,2,5,21,32,36,37,38,43 62biomAP_v_B73_SNPMatrix.txt > biomAP_v_B73_SNPMatrix_7parents.txt

# return to project's home directory
cd ~/projects/genomic_prediction/simulation

# transform to hapmap format:
python scripts/snp_matrix2hapmap.py data/reseq_snps/biomAP_v_B73_SNPMatrix_7parents.txt data/reseq_snps/biomAP_v_B73_SNPMatrix_7parents.hmp.txt

# transform to hapmap diploid format:
run_pipeline.pl -Xmx50g -importGuess data/reseq_snps/biomAP_v_B73_SNPMatrix_7parents.hmp.txt \
                        -export data/reseq_snps/biomAP_v_B73_SNPMatrix_7parents.hmp.txt \
                        -exportType HapmapDiploid
```


### Remove SNPs inside deletions and merge resequencing with SNP chip data

Since the SNP resequencing data will be used together with SV data, I also need to make sure that such SNPs are not present within the boundaries of an SV. The `scripts/merge_reseq_SNPs_with_SNPchip.R` will first remove resequencing SNPs that are within 100kb of a deletion, and then it will merge these SNPs with the SNPs from chip data. Since resequencing data has more than 20 million SNPs, this script will also subset resequencing SNPs by population to keep only the polymorphic ones. Using polymorphic SNPs only will speed up projections by reducing the total number of SNPs to look at.

The markers from the SNP chip will be used as anchors for projection of resequencing SNPs into RILs. The hapmap files from parents and RILs will be saved at a cross specific folder in `data/reseq_snps/`.

```bash
Rscript scripts/merge_reseq_SNPs_with_SNPchip.R data/reseq_snps/biomAP_v_B73_SNPMatrix_7parents.hmp.txt \
                                                data/usda_SVs_7parents.sorted.hmp.txt \
                                                data/usda_biparental-crosses.txt \
                                                analysis/projection \
                                                data/merged_hapmaps_by_cross \
                                                data/reseq_snps/biomAP_parents_SNPs-reseq_SNP-chip_SVs-proj.hmp.txt \
                                                data/reseq_snps/biomAP_rils_SNPs-reseq_SNP-chip_SVs-proj.hmp.txt
```

As shown in the table below, about 1% of SNPs were within deletions up to 100 kb:


| chr | total SNPs | SNPs removed |
| --- | ---------- | ------------ |
| 1   | 3368111    | 33283        |
| 2   | 2936145    | 32428        |
| 3   | 2660241    | 24650        |
| 4   | 2778443    | 20590        |
| 5   | 2349338    | 24793        |
| 6   | 1795656    | 16873        |
| 7   | 2002440    | 17866        |
| 8   | 1917030    | 21385        |
| 9   | 1850524    | 13513        |
| 10  | 1597057    | 13131        |
