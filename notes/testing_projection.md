# Preliminary SV calls dataset

On August 9, 2019, Patrick sent me 5 `.vcf` files containing structural variation calls from the software Lumpy. Each file is a SV call of 100 lines against one of the following reference genomes: B73, Mo17, PH207, PHB47, or W22. **This results are preliminary** and there are a lot of false-positives in there. However, they will be useful for me to write scripts to select the 8 lines I need, project the SVs from parents to RILs, and incorporate these "real" SV data in my simulation scripts (instead of the fake toy dataset).

The files are located at `tests/data/`, and I will use only the SVs called against the B73 reference genome `B73v4_2019-08-09.ls.RT.vcf` for testing.




# VCF to Hapmap format

I wrote the python script `tests/scripts/vcf2hapmap.py` which extract information about structural variants and transform into a numeric hapmap format file sorted by chromosome and positions. Since the `.vcf` file that Patrick sent me contains 100 inbred lines, I had to select only the 7 inbred parents used in the USDA project (you provide this information as an argument in the command line).

The type of SV is displayed in the marker ID (e.g., `del.[ID]` for deletions, and `dup.[ID]` for duplications). Each line will have either a value of `AA` if SV is `A`bsent, or `TT` if SV is `T`here. I had to do that because I will use Tassel to project genotypes from parents to RILs, and it doesn't accept anything other than nucleotides (thus, I couldn't use numbers to indicate presence/absence of SV, as originally thought).

> After projection of parental genotypes into RILs, the hapmap will be transformed into the numeric format for genomic prediction simulations. Thus, for SVs, `AA` will be transformed into `0` and `TT` to `2`.

Also, since SVs spam hundreds (or thousands) of bp and the exact breakpoints are hard to call, the position indicated in the hapmap file wil be the middle point of the SV.

The conversion from VCF to Hapmap was quickly executed by the following commands to generate the file `tests/data/usda_SVs_7parents.sorted.hmp.txt`:

```bash
cd tests/scripts/
# for help on how to use this script
python vcf2hapmap.py --help
# run the script
python vcf2hapmap.py ../data/B73v4_2019-08-09.ls.RT.vcf ../data/usda_SVs_7parents.sorted.hmp.txt B73,LH82,PH207,PHG35,PHG39,PHG47,PHJ40
```




# Projection of SVs from parents to RILS

The overall idea on how to to the projection is to:

1. Merge SV hapmap with SNPs hapmap for each cross.
2. Subset merged hapmap by RIL population.
3. Project using TASSEL FILLIN for each RIL population.

Since this is the first time I'm doing it, I will run some tests first, report them to Candy, and then make changes as necessary.


## 1st test

### Merge SNPs and SVs hapmap files and subset by RIL family

I wrote `tests/scripts/merge_SNPs_SVs_hapmaps.R` to merge the SNP and SV calls of parents and RILs. The results were two hapmap files for each RIL population (one for parents and one for RILs) that can be located at `tests/data/merged_hapmaps_by_cross/`. The hapmap for parents don't have any missing data, since the SVs were called for parents. However, for RILs, all SV information (chromosome, position, etc.) was merged except for the genotypes, which were transformed into missing data (`NN`).

> Important not to have any `NA`! Otherwise, Tassel will think it's a het with alleles `N` and `A`, and will impute everything with `A`. The best to way for Tassel to read it correctly is to use the `NN` notation, so it knows that both alleles are missing.


### Projection

For projection (or imputation), I'm going to use the command-line Tassel plugins `-FILLINFindHaplotypesPlugin`, which generate haplotypes from parents, and `-FILLINImputationPlugin`, which imputes/projects the haplotypes generated from previous plugin to the missing data of the RILs (in this case, the missing data are the SVs).

Testing parameters for projection:

```bash
cd tests/data/merged_hapmap_by_cross

# 1st run:
#   - default values
run_pipeline.pl -FILLINFindHaplotypesPlugin -hmp usda_SNPs-SVs_B73xLH82_parents.sorted.hmp.txt -o donors_B73xLH82
# doesn't do anything because the haplotype block size is 8000 and I have
# between 2000-4000 sites per chromosome

# remove empty directory created before proceeding...
rmdir donors_B73xLH82

# 2nd run:
#   - preferred haplotype size = 1000
run_pipeline.pl -FILLINFindHaplotypesPlugin -hmp usda_SNPs-SVs_B73xLH82_parents.sorted.hmp.txt -o donors_B73xLH82 -hapSize 1000
# also doesn't do anything, hapsize might still be too big. But at this point I
# realized that the minimum number of taxa to generate a haplotype is set to 2.
# Since I want to generate a haplotype for each parent, I believe that this
# option should be set to 1.

# remove empty directory created before proceeding...
rmdir donors_B73xLH82

# 3rd run:
#   - preferred haplotype size = 1000
#   - minimum number of taxa to create haplotype = 1
run_pipeline.pl -FILLINFindHaplotypesPlugin -hmp usda_SNPs-SVs_B73xLH82_parents.sorted.hmp.txt -o donors_B73xLH82 -hapSize 1000 -minTaxa 1
# looks like it worked. So i will try to project this haplotypes to the RILs:
run_pipeline.pl -FILLINImputationPlugin -hmp usda_SNPs-SVs_B73xLH82_RILs.sorted.hmp.txt -d donors_B73xLH82 -o usda_SNPs-SVs_B73xLH82_RILs.projected.hmp.txt -hapSize 1000
# looks like it worked as well! Also looks like there's no more missing data (at least by quick scrolling...)

# 4th run:
#   - preferred haplotype size = 200 (perhaps 1000 was still too much)
#   - minimum number of taxa to create haplotype = 1
#   - need to change the minimum number of present sites from 500 (default) to
#     a number smaller than the haplotype size.
#   - donor name will be donors_B73xLH82_2 to compare with donors_B73xLH82
run_pipeline.pl -FILLINFindHaplotypesPlugin -hmp usda_SNPs-SVs_B73xLH82_parents.sorted.hmp.txt -o donors_B73xLH82_2 -hapSize 200 -minTaxa 1 -minPres 100
# looks like it worked. So i will try to project this haplotypes to the RILs:
run_pipeline.pl -FILLINImputationPlugin -hmp usda_SNPs-SVs_B73xLH82_RILs.sorted.hmp.txt -d donors_B73xLH82_2 -o usda_SNPs-SVs_B73xLH82_RILs.projected_2.hmp.txt -hapSize 200
# didn't work, because the missing data was still missing. Probably the very
# small haplotype block size makes it difficult to define the segments

# 5th run:
#   - preferred haplotype size = 500 (between 200 and 1000)
#   - minimum number of taxa to create haplotype = 1
#   - changed the minimum number of present sites from 500 (default) to 250
#     to allow for some missing sites.
#   - donor name will be donors_B73xLH82_3
run_pipeline.pl -FILLINFindHaplotypesPlugin -hmp usda_SNPs-SVs_B73xLH82_parents.sorted.hmp.txt -o donors_B73xLH82_3 -hapSize 500 -minTaxa 1 -minPres 250
# looks like it worked. So i will try to project this haplotypes to the RILs:
run_pipeline.pl -FILLINImputationPlugin -hmp usda_SNPs-SVs_B73xLH82_RILs.sorted.hmp.txt -d donors_B73xLH82_3 -o usda_SNPs-SVs_B73xLH82_RILs.projected_3.hmp.txt -hapSize 500
# didn't work, because the only some missing data was inputed. There is still
# lots of missing data remaning.

# therefore, the best parameters were from the 3rd run.

# to avoid confusions, I will delete all files generated and re-run the 3rd run
rm -d donors_B73xLH82*
rm donors_B73xLH82*/*
rmdir donors_B73xLH82*
rm usda_SNPs-SVs_B73xLH82_RILs.projected*
```

Based on previous tests, the imputation of parental genotypes into RILs for each population can be done using the following commands:

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation

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
                    -o merged_hapmaps_by_cross/donors_$cross \
                    -hapSize 1000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross \
                    -o imputation/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 1000 -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation/FILLIN_log.txt"

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

This time I used the option `-accuracy` from `-FILLINImputationPlugin` to check how many missing data was imputed and at what accuracy.

Although some RIL populations had all markers projected, other populations still have missing data (4 populations have accuracy < 0.5 and high proportion of unimputed data). The code above seems to be working fine, and it looks like I just need to adjust some parameters when creating or imputing haplotypes. My guess is that the problem happens when creating haplotypes. Taking a look at the `donors_B73xPHG35` folder, I noticed that some blocks of haplotypes (e.g., `donors_B73xPHG35.gc1s0.hmp.txt`) only have the haplotype for one parent. The markers within this block were not imputed. As a comparions, the markers within the block `donors_B73xPHG35.gc1s3.hmp.txt` had haplotypes for both parents and were at least partially imputed.

Since all populations have the same number of markers, it's unlikely that changing the haplotype block size would solve this problem. My guess is that missing SNPs from one parent is causing problems when generating haplotypes. A quick (and very crude) count of `--` in parental genotypes, suggested that files with genotype PHG35 (e.g., `usda_SNPs-SVs_B73xPHG35_parents.sorted.hmp.txt`) had ~3x more missing data (SVs and SNPs) than files without this parent (`usda_SNPs-SVs_B73xLH82_parents.sorted.hmp.txt`). Of course this is not very accurate, but it may be worth trying to follow this lead first. Also, I noticed that PHG35 also have some hets (perhaps more than other parents?), so changing the maximum allowed heterozygosity to be on haplotype may also help (`-mxHet`).

There are two parameters in `-FILLINFindHaplotypesPlugin` that may be worth changing to address this issue. One is called `-minPres`, which is the minimum number of present sites within input sequence to do the search, and the other is `-maxOutMiss`, which is the maximum frequency of missing data in the output haplotype (Default: 0.4). I will run some tests on a temporary folder until I find a reasonable result.

I also tested the FSFHap plugin from Tassel, but imputation was not successful:

```bash
# trying FSFHap
mkdir FSFHap

run_pipeline.pl -h merged_hapmaps_by_cross/usda_SNPs-SVs_B73xLH82_parents-and-rils_chr1_filter.sorted.hmp.txt \
                -FSFHapImputationPlugin -pedigrees usda_RILs_2018_pedigree.txt \
                -logfile FSFHap/log -cluster true -bc false -window 30 \
                -overlap 15 -endPLugin \
                -export FSFHap/usda_SNPs-SVs_B73xLH82_chr1_
```


**Partial conclusions from 1st test:**

It looks like FILLIN plugin from TASSEL should be used for projection (instead of FSFHap). Although some projections went well (i.e., all SV markers were projected from parents to RILs), others didn't work. Also, I still need to check if the version of the reference genome from the SNP chip is the same as the ones used to call SVs (v4). If they are not, I cannot proceed before adjusting the coordinates. Additionally, there is something going one with parent PHG35 that is making it difficult to do projections.


## 2nd test

### Check version of reference genome of SNP chip

One important thing that Candy remind me is to check whether the version of the reference genome used in the SNP chip is the same as used in the SV calls (i.e., refgen v4). Since I will merge the two datasets, it's important that the coordinates are aligned.

To do that, I first downloaded the `.fasta` files of the chromosome 1 for the 4 refgen versions (links below were obtained by searching the MaizeGDB website: <https://www.maizegdb.org/genome>):

```bash
# v1
http://ftp.maizesequence.org/release-4a.53/assembly/
# v2
ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/
# v3
ftp://ftp.ensemblgenomes.org/pub/plants/release-22/fasta/zea_mays/dna
# v4
ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/
```

The original file names were renamed to `v1_chr1.fasta`, `v2_chr1.fasta`, `v3_chr1.fasta`, and `v4_chr1.fasta`, and they were saved at `tests/data` folder.

Then, I extracted the genotypic data of B73 (chromosome 1) from the SNP chip dataset using this very simple R commands, to generate the file `tests/data/markers_b73_chr1.txt`:

```r
library(data.table)

setwd("/Users/rafael/OneDrive/University of Minnesota/PhD/hirsch_lab/projects/genomic_prediction/simulation")

geno.parents <- fread("data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt",
                      header = TRUE, data.table = FALSE)

geno.b73.chr1 <- geno.parents[which(geno.parents[, "chrom"] == 1),
                              c("pos", "B73")]

fwrite(geno.b73.chr1, file = "tests/data/markers_b73_chr1.txt",
       sep = "\t", quote = FALSE)
```

To find out which reference genome the SNP chip dataset was derived from, I wrote the `/tests/scripts/check_refgen_SNPchip.py`. This script extracts the positions from `tests/data/markers_b73_chr1.txt` file and writes the respective positions for each refgen version. Additionally, it counts how many matches there are between each refgen and the SNP chip.

```bash
cd tests/scripts/

python check_refgen_SNPchip.py
# Number of matches to each SNP chip position by refgen version
# B73_v1: 774 (24.2%)
# B73_v2: 2605 (81.3%)
# B73_v3: 786 (24.5%)
# B73_v4: 787 (24.6%)
```

The results suggest that refgen version 2 was used to design the SNP chip, since there is ~80% match. The rest ~20% may be due to the probes used in the chip be tagging the opposite strand of the DNA.

I asked Martin to confirm if it's indeed v2. If it is, I'm gonna need the context probe sequence (50bp probe used to hybridize with the DNA) for each position, so we can properly map them to the refgen v4, and get the right coordinates before merging with SVs.

<mark>TO DO</mark>:
* Candy told me it looks like they cannot provide the probes used, since the chip is proprietary. But they can confirm which version they used. If they say it's v2, i will have to take the positions from the markers and extract 50bp before the position to map to the refgen v4. Pay attention to the strand of the probes. Need to know which ones are positive or negative strand.


### Test accuracy of projection

I wrote `tests/scripts/projection_accuracy_FILLIN.R` to test the accuracy of projecting genotypes from parents to RILs using the FILLIN plugin from TASSEL 5. Here, I will use only SNPs from the my dataset to get a feel of how well FILLIN projection can be in my dataset. Basically, this script create haplotypes and projects SNPs for each individual RIL population, and calculates accuracy by masking random SNPs before projection and counting how many were correctly projected. You can set up how many markers are masked (default is 10%), and whether you want to keep only homozygous SNPs among all parents or use all markers available (including hets and missing markers).

The main ouputs from this script are:
* table containing details of how many genotypes were correctly projected and the source of incorrect calls
* log with accuracy summary per cross and overall accuracy

In addition, this scripts outputs for each RIL population:
* the haplotype files
* hapmap file with projected markers
* text file containing the built-in accuracy results generated by FILLIN

I tested different parameters when running this script. Chose different percentage of markers masked (1, 10 or 50%) and also tested accuracy for datasets with all markers (homo, het and missing) and with only homozygous markers. The results are saved in `tests/analysis/imputation/accuracy`. A summary figure of this test was generated at the same folder and is called `accuracy_summary_plot.png`.

Main results:
* As you **increase** the number markers masked, you **decrease** accuracy.
* Using only homozygous markers from parents to generate haplotypes **increases** projection accuracy.

  | Parameters              | Correct projection | Not imputed | Other mismatches |
  | ----------------------- | ------------------ | ----------- | ---------------- |
  | all markers, 1% masked  | 97.13%             | 0.84%       | 2.03%            |
  | all markers, 10% masked | 94%                | 3.5%        | 2.5%             |
  | all markers, 50% masked | 74.31%             | 22.44%      | 3.25%            |
  | only homo, 1% masked    | 99.12%             | 0.39%       | 0.49%            |
  | only homo, 10% masked   | 96.33%             | 1.79%       | 1.88%            |
  | only homo, 50% masked   | 81.73%             | 14.65%      | 3.62%            |

* Lowest accuracy observed was from cross PHG35xPHG39 (52.44% of genotypes were accurately projected) when using all markers and masking 50% of genotypes.
* In general, crosses with PHG35 as parent had lower projection accuracy.



<mark>TO DO</mark>:
* Map SNPs from v2 to v4 (if Martin confirms SNP chip is v2...)
* Filter out SNPs that fall inside boundaries of SV.
