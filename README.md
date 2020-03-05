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
