# Preliminary SV calls dataset

On August 9, 2019, Patrick sent me 5 `.vcf` files containing structural variation calls from the software Lumpy. Each file is a SV call of 100 lines against one of the following reference genomes: B73, Mo17, PH207, PHB47, or W22. **This results are preliminary** and there are a lot of false-positives in there. However, they will be useful for me to write scripts to select the 8 lines I need, project the SVs from parents to RILs, and incorporate these "real" SV data in my simulation scripts (instead of the fake toy dataset).

The files are located at `tests/data/`, and I will use only the SVs called against the B73 reference genome `B73v4_2019-08-09.ls.RT.vcf` for testing.




# VCF to Hapmap format

I wrote the python script `tests/scripts/vcf2hapmap.py` which extract information about structural variants and transform into a numeric hapmap format file sorted by chromosome and positions. Since the `.vcf` file that Patrick sent me contains 100 inbred lines, I had to select only the 7 inbred parents used in the USDA project (you provide this information as an argument in the command line).

The type of SV is displayed in the marker ID (`del.[ID]` for deletions, `dup.[ID]` for duplications, etc.). Each line will have either a value of `AA` if SV is `A`bsent, or `TT` if SV is `T`here. I had to do that because I will use Tassel to project genotypes from parents to RILs, and it doesn't accept anything other than nucleotides (thus, I couldn't use numbers to indicate presence/absence of SV, as originally thought).

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

Since it was taking a while to have the confirmation, I performed additional QC (using only chromosome 1 data) to see if the mismatch observed between genotyped data and the v2 reference genome was due to complementary bases (i.e., probe used has opposite strand). The results showed that **all markers that didn't match v2 ref gen was actually the complementary nucleotide** (while for other reference genome versions it was only ~30% of markers matching with complementary base). Here are the commands I ran in R to get these results:

```r
# read data
chr1.refs <- read.table("tests/data/ref-markers_chr1.txt", header = TRUE, stringsAsFactors = FALSE)

# do the following for each version of ref gen
for (version in c("v1", "v2", "v3", "v4")) {
  # filter dataset to have only the genotyped marker and the ref gen version
  df <- chr1.refs[, c("marker", version)]
  # keep only markers that differ
  df.diff.alleles <- df[which(df[, "marker"] != df[, version]), ]

  # get proportion of markers that match between genotyped and ref gen version
  n.diff <- NROW(df.diff.alleles)
  n.total <- NROW(df)
  prop.same <- 1 - (n.diff / n.total)

  cat("Same alleles between marker and", version, "is", prop.same, "\n")

  # transform one column to its complement
  df.diff.alleles[, "marker"] <- sapply(df.diff.alleles[, "marker"], FUN = function(base) {
    if (base == "A") {
      base <- "T"
    } else if (base == "T") {
      base <- "A"
    } else if (base == "C") {
      base <- "G"
    } else if (base == "G") {
      base <- "C"
    } else {
      base <- NA
    }
  })

  # remove any missing data
  n.missing <- sum(is.na(df.diff.alleles[, "marker"]))
  cat("Missing alleles removed:", n.missing, "\n")
  df.diff.alleles <- df.diff.alleles[which(!is.na(df.diff.alleles[, "marker"])), ]

  # # get proportion of markers that match between complement genotyped and ref gen version
  n.diff.revcomp <- NROW(df[which(df.diff.alleles[, "marker"] != df.diff.alleles[, version]), ])  # 0
  n.total.revcomp <- NROW(df.diff.alleles)
  prop.same.revcomp <- 1 - (n.diff.revcomp / n.total.revcomp)

  cat("Reverse complement alleles between marker and", version, "is", prop.same.revcomp, "\n\n")
}

# Same alleles between marker and v1 is 0.2416485
# Missing alleles removed: 35
# Reverse complement alleles between marker and v1 is 0.3187135
#
# Same alleles between marker and v2 is 0.8133
# Missing alleles removed: 35
# Reverse complement alleles between marker and v2 is 1
#
# Same alleles between marker and v3 is 0.2453949
# Missing alleles removed: 35
# Reverse complement alleles between marker and v3 is 0.3375315
#
# Same alleles between marker and v4 is 0.2457071
# Missing alleles removed: 35
# Reverse complement alleles between marker and v4 is 0.3175136
```


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



### Mapping SNPs from B73v2 to B73v4 coordinates

Download refgen and extract probes sequences (100bp) at MSI:

```bash
ssh della028@login.msi.umn.edu
ssh mesabi

# home directory is: /home/hirschc1/della028
# make directory to store refgen sequence, and go to that directory
mkdir -p refgen/B73v2
cd refgen/B73v2

# download sequences for all chromosomes
for chr in $(seq 1 10); do
  wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/Zea_mays.AGPv2.dna.chromosome.$chr.fa
done;

# download README
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/README

# extract probes sequences
cd ~/projects/genomic_prediction/simulation
# 100bp probes
python scripts/extract_probe_seqs.py data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt ../../../refgen/B73v2/ 100 data/probes-100bp_22kSNPchip_B73v2.fa
# 200bp probes
python scripts/extract_probe_seqs.py data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt ../../../refgen/B73v2/ 200 data/probes-200bp_22kSNPchip_B73v2.fa
```

> Transfered script `scripts/extract_probe_seqs.py` and SNP chip data `data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt` via FileZilla from my Mac to MSI, and then transfered output `data/probes-100bp_22kSNPchip_B73v2.fa` and `data/probes-200bp_22kSNPchip_B73v2.fa` back to my Mac also with FileZilla.

Download ref gen v4 in MSI and then map probes to the v4 genome:

```bash
ssh della028@login.msi.umn.edu
ssh mesabi

# home directory is: /home/hirschc1/della028
# make directory to store refgen sequence, and go to that directory
mkdir -r refgen/B73v4
cd refgen/B73v4

# download sequences for all chromosomes
for chr in $(seq 1 10); do
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr$chr.fna.gz
done;
# decompress them
gunzip *.gz
# change extension name .fna to .fa for bowtie compatibility
for file in *.fna; do
    mv -- "$file" "${file%.fna}.fa"
done
```

Since it's only a 20k reads to map, I don't need to submit a job via `qsub`. I will just run interactively. I will build bowtie index and map probes to v4 genome directly through the command line:

```bash
ssh della028@login.msi.umn.edu
ssh mesabi

# run qsub interactively
qsub -I -l nodes=1:ppn=8,walltime=2:00:00,mem=10gb

# load bowtie1
module load bowtie/1.1.2

# build index for refgen v4 (this may take ~ an hour since maize genome is big...)
cd ~/refgen/B73v4
bowtie-build chr1.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chr10.fa B73v4

# map probes to v4 genome using bowtie1
cd ~/projects/genomic_prediction/simulation/data
# using 100bp probes
bowtie -f -v 0 -m 1 --un probes-100bp_22kSNPchip_not-aligned-to-B73v4.fa -S ~/refgen/B73v4/B73v4 probes-100bp_22kSNPchip_B73v2.fa probes-100bp_22kSNPchip_aligned-to-B73v4.sam 2> probes-100bp_22kSNPchip_aligned-to-B73v4.stats.txt
# using 200bp probes
bowtie -f -m 1 -v 0 --un probes-200bp_22kSNPchip_not-aligned-to-B73v4.fa -S ~/refgen/B73v4/B73v4 probes-200bp_22kSNPchip_B73v2.fa probes-200bp_22kSNPchip_aligned-to-B73v4.sam 2> probes-200bp_22kSNPchip_aligned-to-B73v4.stats.txt
```

> Bowtie options:
> * `-f`: input is a fasta file
> * `-v 0`: map end-to-end with 0 mismatches allowed
> * `-m 1`: don't report alignments if read map in more than 1 place
> * `--un probes-100bp_22kSNPchip_not-aligned-to-B73v4.fa`: unmapped reads are written to this file
> * `-S`: output is SAM
> * `~/refgen/B73v4/B73v4`: basename of the ref gen v4 index
> * `probes-100bp_22kSNPchip_B73v2.fa`: reads to map
> * `probes-100bp_22kSNPchip_aligned-to-B73v4.sam`: output name
> * `2> probes-100bp_22kSNPchip_aligned-to-B73v4.stats.txt`: this will print alignment statistics from std err to this file

All output generated was transfered back to my Mac via FileZilla at the folder `data/probes_v2-to-v4`. Here are the stats of the alignments:

|              | reads processed | uniquely mapped | multimapped (supressed) | unmapped |
| ------------ | --------------- | --------------- | ----------------------- | -------- |
| 100bp probes | 20,139          | 19,850          | 72                      | 217      |
| 200bp probes | 20,139          | 19,787          | 62                      | 290      |



Now, I need to correct the positions of the SNPs in the hapmap files `data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt` and `data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt` based on new coordinates from SAM file. I wrote a python script called `scripts/convert_hmp_v2-to-v4.py` to do that and create the files `data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt` and `data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt`. This script also removes reads that uniquely mapped with no mismatch to a different chromosome in v4. To do that, I ran:

```bash
# 100bp probes only, since difference in mapping was very little
python scripts/convert_hmp_v2-to-v4.py data/probes_v2-to-v4/probes-100bp_22kSNPchip_aligned-to-B73v4.sam data/probes_v2-to-v4/SNP_positions_v2-to-v4_probes-100bp.txt data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt,data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt
# 23 reads discarded due to mapping into different chromosome in v4
```

I performed a quick QC (using only chromosome 1 data) to see if SNPs from v4 were actually the same as in v2. The results showed that markers that didn't match between v2 and v4 were due to mapping to complementary strand. Here are the commands I ran in R to get these results:

```r
# read data
v4.pos <- read.table("data/probes_v2-to-v4/SNP_positions_v2-to-v4_probes-100bp.txt", header = TRUE, stringsAsFactors = FALSE)
v2.pos <- read.table("tests/data/ref-markers_chr1.txt", header = TRUE, stringsAsFactors = FALSE)

# filter v4 data to have only chromosome 1 markers
v4.pos.chr1 <- v4.pos[which(v4.pos[, "chr_v4"] == 1), ]

# make sure that the same markers are in both datasets (because after correcting SNP positions in v4,
# some SNPs from chr1 went to other chromosome; so these were excluded from this analysis)
v2.pos <- v2.pos[v2.pos[, "pos"] %in% v4.pos.chr1[, "pos_v2"], ]
v4.pos.chr1 <- v4.pos.chr1[v4.pos.chr1[, "pos_v2"] %in% v2.pos[, "pos"], ]

# create a new data frame only with the respective SNPs in v2 and v4
df <- cbind(v2.pos["v2"], v4.pos.chr1["SNP"])
colnames(df) <- c("v2", "v4")

# count how many were different in v2 and v4
df_diff_alleles <- df[which(df$v2 != df$v4), ]
NROW(df_diff_alleles)
# 217

# transform one of the columns into the complementary base
df_diff_alleles$v2 <- sapply(df_diff_alleles$v2, FUN = function(base) {
  if (base == "A") {
    base <- "T"
  } else if (base == "T") {
    base <- "A"
  } else if (base == "C") {
    base <- "G"
  } else if (base == "G") {
    base <- "C"
  }
})

# now, check again how many markers are different
NROW(df[which(df_diff_alleles$v2 != df_diff_alleles$v4), ])
# 0
```

Ran `scripts/correct_SNP_strands.R` to correct strands of `data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt` and `data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt`...


Correct the `alleles` column after correcting strands:

```bash
run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt -export data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt -export data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt -exportType HapmapDiploid
```



### Use resequencing data for PHG35

Christine sent me resequencing data contaning SNP calls for PHG35 (the bad parent). This is direct output from FreeBayes (software to call variants) with no filtering for quality, contain more maize lines other than PHG35 and also contain indels. Thus I need to filter SNPs by quality, select only PHG35 and then get only the SNPs that are in the SNP chip (~20k SNPs). This is already in v4 coordinates. The file was copied to the following path of my MSI account: `/home/hirschc1/della028/projects/genomic_prediction/simulation/data/WiDiv100_B73v4_allchr_gene1kb_skip.nbest4.k.vcf`

Before filtering the vcf file above, I need to extract the SNP positions. In my Mac, I executed the following commands on the terminal:

```bash
# get only 3rd and 4th columns, skip the header, and print to a new file
cut -f 3,4 data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt | sed 1d > tests/data/22k_SNPs_in_v4_coord.txt
```

> The file `tests/data/22k_SNPs_in_v4_coord.txt` was transfered from my Mac to MSI via FileZilla at `data/22k_SNPs_in_v4_coord.txt`

Filter VCF file Christine provided me to have only the PHG35 genotype for the positions in the SNPchip, and transform to hapmap format using TASSEL:

```bash
# filter vcf
vcftools --vcf data/WiDiv100_B73v4_allchr_gene1kb_skip.nbest4.k.vcf \
         --indv PHG35 \
         --positions data/22k_SNPs_in_v4_coord.txt \
         --out data/PHG35_reseq_B73v4_coord \
         --recode \
         --recode-INFO-all

 vcftools --vcf data/WiDiv100_B73v4_allchr_gene1kb_skip.nbest4.k.vcf \
          --indv B73 --indv LH82 --indv PH207 --indv PHG35 --indv PHG39 --indv PHG47 --indv PHJ40 \
          --positions data/22k_SNPs_in_v4_coord.txt \
          --out data/usda-parents_reseq_B73v4_coord \
          --recode \
          --recode-INFO-all

# transform to hapmap
run_pipeline.pl -Xmx10g -importGuess data/PHG35_reseq_B73v4_coord.recode.vcf -export data/PHG35_reseq_B73v4_coord.hmp.txt -exportType HapmapDiploid  # kept only 15429 markers

run_pipeline.pl -Xmx10g -importGuess data/usda-parents_reseq_B73v4_coord.recode.vcf -export data/usda-parents_reseq_B73v4_coord.hmp.txt -exportType HapmapDiploid  # kept only 15429 markers
```

> The file `data/PHG35_reseq_B73v4_coord.hmp.txt` (and parents...) was transfered from MSI to my Mac via FileZilla at `tests/data/PHG35_reseq_B73v4_coord.hmp.txt`.


<mark>TO DO:</mark>
* Filter by SNP quality as well


This kept only 15,429 markers. This means that ~4k SNPs present in the chip were not present in the resequencing PHG35. Now I have to change the genotypes calls for PHG35 in `data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt` with the data from `tests/data/PHG35_reseq_B73v4_coord.hmp.txt`, and also make sure to add `NN` to SNPs not present in the resequencing of PHG35. To do that, I ran `scripts/change_PHG35_SNPs.R`. The file generated was `tests/data/usda_22kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt`. Did the same thing for parents (resequencing)...

Correct the `alleles` column after changing the PHG35 calls:

```bash
run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt -export data/usda_22kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_7parents.sorted.diploid.v4.all-parents-corrected.hmp.txt -export data/usda_22kSNPs_7parents.sorted.diploid.v4.all-parents-corrected.hmp.txt -exportType HapmapDiploid
```



Alternatively, I excluded SNPs that were not present in PHG35 using script `scripts/keep_only_PHG35_SNPs.R`. The file generated was `tests/data/usda_15kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt`...

Correct the `alleles` column after changing the PHG35 calls:

```bash
run_pipeline.pl -Xmx10g -importGuess data/usda_15kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt -export data/usda_15kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_15kSNPs_325rils.sorted.diploid.v4.PHG35-corrected.hmp.txt -export data/usda_15kSNPs_325rils.sorted.diploid.v4.PHG35-corrected.hmp.txt -exportType HapmapDiploid
```



Alternatively, I excluded SNPs that were not present in all parents using script `scripts/keep_only_PHG35_SNPs.R`. The file generated was `tests/data/usda_15kSNPs_7parents.sorted.diploid.v4.all-parents-corrected.hmp.txt`...

```bash
run_pipeline.pl -Xmx10g -importGuess data/usda_15kSNPs_7parents.sorted.diploid.v4.all-parents-corrected.hmp.txt -export data/usda_15kSNPs_7parents.sorted.diploid.v4.all-parents-corrected.hmp.txt -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_15kSNPs_325rils.sorted.diploid.v4.all-parents-corrected.hmp.txt -export data/usda_15kSNPs_325rils.sorted.diploid.v4.all-parents-corrected.hmp.txt -exportType HapmapDiploid
```



### Filter out SNPs inside SV boundaries

Need to know SV size first. Perhaps I can add the boundaries in the name of the SV by changing line 89 of `tests/scripts/vcf2hapmap.py` to:

```python
id = sv_type + "." + chr + "." + str(sv_start) + "." + str(sv_end)
```

Then I have to filter SNPs that fall into SV boundaries by modifying the script `tests/scripts/merge_SNPs_SVs_hapmap.R`.

SNPs within any SV (dels and dups):

| chr | total SNPs | SNPs within any SV (dels and dups) | SNPs within dels only |
| --- | ---------- | ---------------------------------- | --------------------- |
| 1   | 3143       | 2594                               | 945                   |
| 2   | 2397       | 1435                               | 2                     |
| 3   | 2547       | 2301                               | 3                     |
| 4   | 2079       | 1323                               | 0                     |
| 5   | 2128       | 1979                               | 2                     |
| 6   | 1755       | 1635                               | 0                     |
| 7   | 1598       | 1521                               | 0                     |
| 8   | 1295       | 1108                               | 2                     |
| 9   | 1518       | 964                                | 2                     |
| 10  | 1306       | 218                                | 0                     |

If I use only SNPs around any SVs, I would end up with only 4,688 markers! I think it makes sense to exclude a SNP that is within a deletion and keep those that are within a duplication (since all lines would have it anyways, others just in more copies). By doing that, I can keep
18,814 markers.


Also did with `tests/data/usda_15kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt`...



### Merge SNPs and SVs hapmap files

`tests/scripts/merge_SNPs_SVs_hapmap.R`



### Project SVs from parents to RILs


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

Transform to diploid hapmap:

```bash
cd tests/data/imputation

for file in $(ls *.projected.hmp.txt); do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



Project using 15k snps present in PHG35:

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_15k

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
                    -hmp merged_hapmaps_by_cross_15k/usda_SNPs-SVs_$cross\_parents.sorted.hmp.txt \
                    -o merged_hapmaps_by_cross_15k/donors_$cross \
                    -hapSize 1000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross_15k/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross_15k/donors_$cross \
                    -o imputation_15k/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 1000 -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_15k/FILLIN_log.txt"

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross_15k/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_15k

for file in $(ls *.projected.hmp.txt); do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```

Not much difference between 20k and 15k imputation.


`scripts/count_projected_SVs.R`



Merge projected files into one with `scripts/merge_projected_crosses.R`. Created file `tests/data/usda_SNPs-SVs_325rils.not-in-SVs.projected.hmp.txt`.



#### Plot karyotypes

`tests/scripts/plot_ril_karyotype_SVs.R`




#### Reconstruct PHG35 based on RIL data

For each family with PHG35 as a parent:
- Check which allele is not from PHG35 based on resequencing data.
- Figure out which is the alternate allele based on RIL calls.
- Compare if the allele from RILs is the same from resequencing of PHG35.

`scripts/reconstruct_PHG35_from_RIL_data.R`:

Problem with recontruction:  what if I have monomorphic call for a SNP in RIL data mixed with missing data (e.g. G/N), and the parents data have the allele called in RILs (e.g G/C)? I can't know if the SNP is actually monomorphic, or if the missing data would actually be a different allele.

```
Analyzing cross LH82xPHG35 (4933 disagreements)
  145 (2.94%) of SNPs that disagree between parents and RILs are due to missing data both parents
  2002 (40.58%) of SNPs that disagree between parents and RILs are due to missing data in PHG35 parent
  165 (3.34%) of SNPs that disagree between parents and RILs are due to missing data in the non-PHG35 parent
  2486 (50.4%) of SNPs that disagree between parents and RILs are due to uncertainty in SNP being mono or polymorphic
  135 (2.74%) of SNPs that disagree between parents and RILs are due to PHG35 allele being different between resequencing and SNP data

Analyzing cross PHG35xPHG47 (2147 disagreements)
  141 (6.57%) of SNPs that disagree between parents and RILs are due to missing data both parents
  1559 (72.61%) of SNPs that disagree between parents and RILs are due to missing data in PHG35 parent
  106 (4.94%) of SNPs that disagree between parents and RILs are due to missing data in the non-PHG35 parent
  229 (10.67%) of SNPs that disagree between parents and RILs are due to uncertainty in SNP being mono or polymorphic
  112 (5.22%) of SNPs that disagree between parents and RILs are due to PHG35 allele being different between resequencing and SNP data

Analyzing cross B73xPHG35 (2993 disagreements)
  54 (1.8%) of SNPs that disagree between parents and RILs are due to missing data both parents
  2102 (70.23%) of SNPs that disagree between parents and RILs are due to missing data in PHG35 parent
  50 (1.67%) of SNPs that disagree between parents and RILs are due to missing data in the non-PHG35 parent
  642 (21.45%) of SNPs that disagree between parents and RILs are due to uncertainty in SNP being mono or polymorphic
  145 (4.84%) of SNPs that disagree between parents and RILs are due to PHG35 allele being different between resequencing and SNP data

Analyzing cross PHG35xPHG39 (4000 disagreements)
  102 (2.55%) of SNPs that disagree between parents and RILs are due to missing data both parents
  2111 (52.78%) of SNPs that disagree between parents and RILs are due to missing data in PHG35 parent
  106 (2.65%) of SNPs that disagree between parents and RILs are due to missing data in the non-PHG35 parent
  1548 (38.7%) of SNPs that disagree between parents and RILs are due to uncertainty in SNP being mono or polymorphic
  133 (3.33%) of SNPs that disagree between parents and RILs are due to PHG35 allele being different between resequencing and SNP data
```



### Adjusting projection parameters

#### `-hapSize 500`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted \
                    -hapSize 500 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted \
                    -o imputation_adjusted/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 500 -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_adjusted

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`



#### `-hapSize 500` and `-hybNN`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted2

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted \
                    -hapSize 500 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted \
                    -o imputation_adjusted2/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 500 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted2/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_adjusted2

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`




#### `-hapSize 200` and `-hybNN`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted3

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted3 \
                    -hapSize 200 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted3 \
                    -o imputation_adjusted3/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 200 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted3/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

> `-hapSize 200` was too low. Donor files were not even generated.




#### `-hapSize 300` and `-hybNN`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted4

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted4 \
                    -hapSize 300 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted4 \
                    -o imputation_adjusted4/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 300 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted4/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_adjusted4

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`
> `-hapSize 1000` > `-hapSize 500` > `-hapSize 300`




#### `-hapSize 750` and `-hybNN`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted5

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted5 \
                    -hapSize 750 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted5 \
                    -o imputation_adjusted5/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 750 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted5/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_adjusted5

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was worse than with `-hapSize 1000`
> `-hapSize 1000` > `-hapSize 750` > `-hapSize 500` > `-hapSize 300`



#### `-hapSize 1250` and `-hybNN`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted6

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted6 \
                    -hapSize 1250 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted6 \
                    -o imputation_adjusted6/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 1250 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted6/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_adjusted6

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was very similar to `-hapSize 1000`. Perhaps a bit better?


#### `-hapSize 1500` and `-hybNN`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted7

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted7 \
                    -hapSize 1500 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted7 \
                    -o imputation_adjusted7/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 1500 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted7/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_adjusted7

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was very similar to `-hapSize 1000`. Perhaps a bit better?




#### `-hapSize 2000` and `-hybNN`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted8

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted8 \
                    -hapSize 2000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted8 \
                    -o imputation_adjusted8/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 2000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted8/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_adjusted8

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was very similar to `-hapSize 1000`. Perhaps a bit better?



#### `-hapSize 2500` and `-hybNN`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted9

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted9 \
                    -hapSize 2500 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted9 \
                    -o imputation_adjusted9/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 2500 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted9/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_adjusted9

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was very similar to `-hapSize 1000`. Perhaps a bit better?



#### `-hapSize 3000` and `-hybNN`

```bash
# change directory
cd tests/data/

# create directory to store results of imputation
mkdir imputation_adjusted10

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
                    -o merged_hapmaps_by_cross/donors_$cross\_adjusted10 \
                    -hapSize 3000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross\_adjusted10 \
                    -o imputation_adjusted10/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 3000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "imputation_adjusted10/FILLIN_log.txt"
# set delimiter to whitespace
IFS=

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd tests/data/imputation_adjusted10

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Basically everything stayed the same as `-hapSize 2500`


#### Results

* `-hapSize 1000` was better than smaller sizes.
* From `-hapSize 1000` to `-hapSize 2500` some populations improved gradually, some improved and then decreased, and some stayed the same.
* So, I think I have to adjust parameters differently for each populations. Those that had projection above 75% will be run with `-hapSize 1000 -hybNN false`.
* Those below 75% will be run separatelly with specific parameters to try to increase the number of projected SVs. The populations are:
  - B73 x PHG35
  - LH82 x PHG35
  - PHG35 x PHG39
  - PHG39 x PHG47



#### Try changing `-minPres` as well

Minimum number of present sites within input sequence to do the search (Default: 500).

Scale according to


<mark>CORRECT ALLELES FIELD FROM RILS DATA AFTER PROJECTION!</mark>





## 3rd test

### VCF to Hapmap format

I wrote the python script `tests/scripts/vcf2hapmap.py` which extract information about structural variants and transform into a hapmap format file sorted by chromosome and positions. Since the `.vcf` file that Patrick sent me contains 100 inbred lines, I have to select only the 7 inbred parents used in the USDA project (you provide this information as an argument in the command line).

The type of SV is displayed in the marker ID (`del.[ID]` for deletions, `dup.[ID]` for duplications, etc.). Each line will have either a value of `AA` if SV is `A`bsent, or `TT` if SV is `T`here. I had to do that because I will use Tassel to project genotypes from parents to RILs, and it doesn't accept anything other than nucleotides (thus, I couldn't use numbers to indicate presence/absence of SV, as originally thought). Missing data was coded as `NN`.

> After projection of parental genotypes into RILs, the hapmap will be transformed into the numeric format for genomic prediction simulations. Thus, for SVs, `AA` will be transformed into `0` and `TT` to `2`.

Also, since SVs spam hundreds (or thousands) of bp and the exact breakpoints are hard to call, the position indicated in the hapmap file wil be the middle point of the SV.

The conversion from VCF to Hapmap was quickly executed by the following commands to generate the file `data/usda_SVs_7parents.sorted.hmp.txt`:

```bash
# for help on how to use this script
python vcf2hapmap.py --help
# run the script
python scripts/vcf2hapmap.py tests/data/B73v4_2019-08-09.ls.RT.vcf data/usda_SVs_7parents.sorted.hmp.txt B73,LH82,PH207,PHG35,PHG39,PHG47,PHJ40
```


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

Since it was taking a while to have the confirmation, I performed additional QC (using only chromosome 1 data) to see if the mismatch observed between genotyped data and the v2 reference genome was due to complementary bases (i.e., probe used has opposite strand). The results showed that **all markers that didn't match v2 ref gen was actually the complementary nucleotide** (while for other reference genome versions it was only ~30% of markers matching with complementary base). Here are the commands I ran in R to get these results:

```r
# read data
chr1.refs <- read.table("tests/data/ref-markers_chr1.txt", header = TRUE, stringsAsFactors = FALSE)

# do the following for each version of ref gen
for (version in c("v1", "v2", "v3", "v4")) {
  # filter dataset to have only the genotyped marker and the ref gen version
  df <- chr1.refs[, c("marker", version)]
  # keep only markers that differ
  df.diff.alleles <- df[which(df[, "marker"] != df[, version]), ]

  # get proportion of markers that match between genotyped and ref gen version
  n.diff <- NROW(df.diff.alleles)
  n.total <- NROW(df)
  prop.same <- 1 - (n.diff / n.total)

  cat("Same alleles between marker and", version, "is", prop.same, "\n")

  # transform one column to its complement
  df.diff.alleles[, "marker"] <- sapply(df.diff.alleles[, "marker"], FUN = function(base) {
    if (base == "A") {
      base <- "T"
    } else if (base == "T") {
      base <- "A"
    } else if (base == "C") {
      base <- "G"
    } else if (base == "G") {
      base <- "C"
    } else {
      base <- NA
    }
  })

  # remove any missing data
  n.missing <- sum(is.na(df.diff.alleles[, "marker"]))
  cat("Missing alleles removed:", n.missing, "\n")
  df.diff.alleles <- df.diff.alleles[which(!is.na(df.diff.alleles[, "marker"])), ]

  # # get proportion of markers that match between complement genotyped and ref gen version
  n.diff.revcomp <- NROW(df[which(df.diff.alleles[, "marker"] != df.diff.alleles[, version]), ])  # 0
  n.total.revcomp <- NROW(df.diff.alleles)
  prop.same.revcomp <- 1 - (n.diff.revcomp / n.total.revcomp)

  cat("Reverse complement alleles between marker and", version, "is", prop.same.revcomp, "\n\n")
}

# Same alleles between marker and v1 is 0.2416485
# Missing alleles removed: 35
# Reverse complement alleles between marker and v1 is 0.3187135
#
# Same alleles between marker and v2 is 0.8133
# Missing alleles removed: 35
# Reverse complement alleles between marker and v2 is 1
#
# Same alleles between marker and v3 is 0.2453949
# Missing alleles removed: 35
# Reverse complement alleles between marker and v3 is 0.3375315
#
# Same alleles between marker and v4 is 0.2457071
# Missing alleles removed: 35
# Reverse complement alleles between marker and v4 is 0.3175136
```



### Mapping SNPs from B73v2 to B73v4 coordinates

Download refgen and extract probes sequences (100bp) at MSI:

```bash
ssh della028@login.msi.umn.edu
ssh mesabi

# home directory is: /home/hirschc1/della028
# make directory to store refgen sequence, and go to that directory
mkdir -p refgen/B73v2
cd refgen/B73v2

# download sequences for all chromosomes
for chr in $(seq 1 10); do
  wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/Zea_mays.AGPv2.dna.chromosome.$chr.fa
done;

# download README
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/README

# extract probes sequences
cd ~/projects/genomic_prediction/simulation
# 100bp probes
python scripts/extract_probe_seqs.py data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt ../../../refgen/B73v2/ 100 data/probes-100bp_22kSNPchip_B73v2.fa
# 200bp probes
python scripts/extract_probe_seqs.py data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt ../../../refgen/B73v2/ 200 data/probes-200bp_22kSNPchip_B73v2.fa
```

> Transfered script `scripts/extract_probe_seqs.py` and SNP chip data `data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt` via FileZilla from my Mac to MSI, and then transfered output `data/probes-100bp_22kSNPchip_B73v2.fa` and `data/probes-200bp_22kSNPchip_B73v2.fa` back to my Mac also with FileZilla.

Download ref gen v4 in MSI and then map probes to the v4 genome:

```bash
ssh della028@login.msi.umn.edu
ssh mesabi

# home directory is: /home/hirschc1/della028
# make directory to store refgen sequence, and go to that directory
mkdir -r refgen/B73v4
cd refgen/B73v4

# download sequences for all chromosomes
for chr in $(seq 1 10); do
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr$chr.fna.gz
done;
# decompress them
gunzip *.gz
# change extension name .fna to .fa for bowtie compatibility
for file in *.fna; do
    mv -- "$file" "${file%.fna}.fa"
done
```

Since it's only a 20k reads to map, I don't need to submit a job via `qsub`. I will just run interactively. I will build bowtie index and map probes to v4 genome directly through the command line:

```bash
ssh della028@login.msi.umn.edu
ssh mesabi

# run qsub interactively
qsub -I -l nodes=1:ppn=8,walltime=2:00:00,mem=10gb

# load bowtie1
module load bowtie/1.1.2

# build index for refgen v4 (this may take ~ an hour since maize genome is big...)
cd ~/refgen/B73v4
bowtie-build chr1.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chr10.fa B73v4

# map probes to v4 genome using bowtie1
cd ~/projects/genomic_prediction/simulation/data
# using 100bp probes
bowtie -f -v 0 -m 1 --un probes-100bp_22kSNPchip_not-aligned-to-B73v4.fa -S ~/refgen/B73v4/B73v4 probes-100bp_22kSNPchip_B73v2.fa probes-100bp_22kSNPchip_aligned-to-B73v4.sam 2> probes-100bp_22kSNPchip_aligned-to-B73v4.stats.txt
# using 200bp probes
bowtie -f -m 1 -v 0 --un probes-200bp_22kSNPchip_not-aligned-to-B73v4.fa -S ~/refgen/B73v4/B73v4 probes-200bp_22kSNPchip_B73v2.fa probes-200bp_22kSNPchip_aligned-to-B73v4.sam 2> probes-200bp_22kSNPchip_aligned-to-B73v4.stats.txt
```

> Bowtie options:
> * `-f`: input is a fasta file
> * `-v 0`: map end-to-end with 0 mismatches allowed
> * `-m 1`: don't report alignments if read map in more than 1 place
> * `--un probes-100bp_22kSNPchip_not-aligned-to-B73v4.fa`: unmapped reads are written to this file
> * `-S`: output is SAM
> * `~/refgen/B73v4/B73v4`: basename of the ref gen v4 index
> * `probes-100bp_22kSNPchip_B73v2.fa`: reads to map
> * `probes-100bp_22kSNPchip_aligned-to-B73v4.sam`: output name
> * `2> probes-100bp_22kSNPchip_aligned-to-B73v4.stats.txt`: this will print alignment statistics from std err to this file

All output generated was transfered back to my Mac via FileZilla at the folder `data/probes_v2-to-v4`. Here are the stats of the alignments:

|              | reads processed | uniquely mapped | multimapped (supressed) | unmapped |
| ------------ | --------------- | --------------- | ----------------------- | -------- |
| 100bp probes | 20,139          | 19,850          | 72                      | 217      |
| 200bp probes | 20,139          | 19,787          | 62                      | 290      |



Now, I need to correct the positions of the SNPs in the hapmap files `data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt` and `data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt` based on new coordinates from SAM file. I wrote a python script called `scripts/convert_hmp_v2-to-v4.py` to do that and create the files `data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt` and `data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt`. This script also removes reads that uniquely mapped with no mismatch to a different chromosome in v4. To do that, I ran:

```bash
# 100bp probes only, since difference in mapping was very little
python scripts/convert_hmp_v2-to-v4.py data/probes_v2-to-v4/probes-100bp_22kSNPchip_aligned-to-B73v4.sam data/probes_v2-to-v4/SNP_positions_v2-to-v4_probes-100bp.txt data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt,data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt
# 23 reads discarded due to mapping into different chromosome in v4
```

I performed a quick QC (using only chromosome 1 data) to see if SNPs from v4 were actually the same as in v2. The results showed that markers that didn't match between v2 and v4 were due to mapping to complementary strand. Here are the commands I ran in R to get these results:

```r
# read data
v4.pos <- read.table("data/probes_v2-to-v4/SNP_positions_v2-to-v4_probes-100bp.txt", header = TRUE, stringsAsFactors = FALSE)
v2.pos <- read.table("tests/data/ref-markers_chr1.txt", header = TRUE, stringsAsFactors = FALSE)

# filter v4 data to have only chromosome 1 markers
v4.pos.chr1 <- v4.pos[which(v4.pos[, "chr_v4"] == 1), ]

# make sure that the same markers are in both datasets (because after correcting SNP positions in v4,
# some SNPs from chr1 went to other chromosome; so these were excluded from this analysis)
v2.pos <- v2.pos[v2.pos[, "pos"] %in% v4.pos.chr1[, "pos_v2"], ]
v4.pos.chr1 <- v4.pos.chr1[v4.pos.chr1[, "pos_v2"] %in% v2.pos[, "pos"], ]

# create a new data frame only with the respective SNPs in v2 and v4
df <- cbind(v2.pos["v2"], v4.pos.chr1["SNP"])
colnames(df) <- c("v2", "v4")

# count how many were different in v2 and v4
df_diff_alleles <- df[which(df$v2 != df$v4), ]
NROW(df_diff_alleles)
# 217

# transform one of the columns into the complementary base
df_diff_alleles$v2 <- sapply(df_diff_alleles$v2, FUN = function(base) {
  if (base == "A") {
    base <- "T"
  } else if (base == "T") {
    base <- "A"
  } else if (base == "C") {
    base <- "G"
  } else if (base == "G") {
    base <- "C"
  }
})

# now, check again how many markers are different
NROW(df[which(df_diff_alleles$v2 != df_diff_alleles$v4), ])
# 0
```

Ran `scripts/correct_SNP_strands.R` to correct strands of `data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt` and `data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt`...


Correct the `alleles` column after correcting strands:

```bash
run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt -export data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt -export data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt -exportType HapmapDiploid
```



### Reconstructed PHG35 based on RIL genotypes

Using resequencing data didn't help much because we still cannot be sure if the seed source are the same between PHG35 used in this project and the PHG35 used in resequencing.


`scripts/reconstruct_PHG35_from_RIL_data.R`


#### Summary of markers per population using reconstructed PHG35

`scripts/markers_summary.R`



### Filter out SNPs inside deletion boundaries



Then I have to filter SNPs that fall into SV boundaries by modifying the script `scripts/merge_SNPs_SVs_hapmap.R`. Created `data/usda_18kSNPs_7parents.not-in-PAVs.hmp.txt`, `data/usda_18kSNPs_325rils.not-in-PAVs.hmp.txt` and files by cross in `data/merged_hapmaps_by_cross/`.


<mark>Add function to filter SNPs by all SVs, just to compare with filtering by deletions only</mark>

If I use only SNPs around any SVs, I would end up with only ~4k markers! I think it makes sense to exclude a SNP that is within a deletion and keep those that are within a duplication (since all lines would have it anyways, others just in more copies). By doing that, I can keep ~18k markers.



### Project SVs from parents to RILs


```bash
# create directory to store results of imputation
mkdir -p analysis/reseq_projection

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
                    -o merged_hapmaps_by_cross/donors_$cross \
                    -hapSize 1000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d merged_hapmaps_by_cross/donors_$cross \
                    -o ../analysis/projection/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 1000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/reseq_projection/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir merged_hapmaps_by_cross/*
```

Transform to diploid hapmap:

```bash
cd analysis/projection

for file in $(ls *projected.hmp.txt); do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



### Adjusting projection parameters

#### Project SVs from parents to RILs (using filtered SNPs per population)


```bash
# create directory to store results of imputation
mkdir -p analysis/projection_filtered-SVs

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
                    -hmp hapmaps_not-in-dels_by_cross/usda_SNPs-SVs_$cross\_parents.sorted.hmp.txt \
                    -o hapmaps_not-in-dels_by_cross/donors_$cross \
                    -hapSize 1000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp hapmaps_not-in-dels_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d hapmaps_not-in-dels_by_cross/donors_$cross \
                    -o ../analysis/projection_filtered-SVs/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 1000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_filtered-SVs/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir hapmaps_not-in-dels_by_cross/*
```

Transform to diploid hapmap:

```bash
cd analysis/projection_filtered-SVs

for file in $(ls *projected.hmp.txt); do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```




#### `-hapSize 500`, `-minPres 50` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs500_mp50

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
                    -o ../analysis/projection_hs500_mp50/donors_$cross \
                    -hapSize 500 -minPres 50 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs500_mp50/donors_$cross \
                    -o ../analysis/projection_hs500_mp50/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 500 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs500_mp50/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs500_mp50/*


# transform to diploid hapmap:
cd ../analysis/projection_hs500_mp50

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`





#### `-hapSize 500`, `-minPres 200` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs500_mp200

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
                    -o ../analysis/projection_hs500_mp200/donors_$cross \
                    -hapSize 500 -minPres 200 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs500_mp200/donors_$cross \
                    -o ../analysis/projection_hs500_mp200/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 500 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs500_mp200/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs500_mp200/*


# transform to diploid hapmap:
cd ../analysis/projection_hs500_mp200

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`



#### `-hapSize 1000`, `-minPres 200` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs1000_mp200

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
                    -o ../analysis/projection_hs1000_mp200/donors_$cross \
                    -hapSize 1000 -minPres 200 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs1000_mp200/donors_$cross \
                    -o ../analysis/projection_hs1000_mp200/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 1000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs1000_mp200/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs1000_mp200/*


# transform to diploid hapmap:
cd ../analysis/projection_hs1000_mp200

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```



`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`






#### `-hapSize 1500`, `-minPres 500` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs1500_mp500

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
                    -o ../analysis/projection_hs1500_mp500/donors_$cross \
                    -hapSize 1500 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs1500_mp500/donors_$cross \
                    -o ../analysis/projection_hs1500_mp500/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 1500 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs1500_mp500/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs1500_mp500/*


# transform to diploid hapmap:
cd ../analysis/projection_hs1500_mp500

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```

`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`



#### `-hapSize 2000`, `-minPres 500` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs2000_mp500

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
                    -o ../analysis/projection_hs2000_mp500/donors_$cross \
                    -hapSize 2000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs2000_mp500/donors_$cross \
                    -o ../analysis/projection_hs2000_mp500/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 2000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs2000_mp500/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs2000_mp500/*


# transform to diploid hapmap:
cd ../analysis/projection_hs2000_mp500

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```

`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`





#### `-hapSize 2000`, `-minPres 1000` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs2000_mp1000

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
                    -o ../analysis/projection_hs2000_mp1000/donors_$cross \
                    -hapSize 2000 -minPres 1000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs2000_mp1000/donors_$cross \
                    -o ../analysis/projection_hs2000_mp1000/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 2000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs2000_mp1000/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs2000_mp1000/*


# transform to diploid hapmap:
cd ../analysis/projection_hs2000_mp1000

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```

`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`






#### `-hapSize 2500`, `-minPres 500` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs2500_mp500

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
                    -o ../analysis/projection_hs2500_mp500/donors_$cross \
                    -hapSize 2500 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs2500_mp500/donors_$cross \
                    -o ../analysis/projection_hs2500_mp500/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 2500 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs2500_mp500/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs2500_mp500/*


# transform to diploid hapmap:
cd ../analysis/projection_hs2500_mp500

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```

`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`




#### `-hapSize 2500`, `-minPres 1000` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs2500_mp1000

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
                    -o ../analysis/projection_hs2500_mp1000/donors_$cross \
                    -hapSize 2500 -minPres 1000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs2500_mp1000/donors_$cross \
                    -o ../analysis/projection_hs2500_mp1000/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 2500 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs2500_mp1000/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs2500_mp1000/*


# transform to diploid hapmap:
cd ../analysis/projection_hs2500_mp1000

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```

`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`




#### `-hapSize 3000`, `-minPres 500` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs3000_mp500

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
                    -o ../analysis/projection_hs3000_mp500/donors_$cross \
                    -hapSize 3000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs3000_mp500/donors_$cross \
                    -o ../analysis/projection_hs3000_mp500/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 3000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs3000_mp500/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs3000_mp500/*


# transform to diploid hapmap:
cd ../analysis/projection_hs3000_mp500

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```

`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`




#### `-hapSize 3000`, `-minPres 1000` and `-hybNN`

```bash
# create directory to store results of imputation
mkdir -p analysis/projection_hs3000_mp1000

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
                    -o ../analysis/projection_hs3000_mp1000/donors_$cross \
                    -hapSize 3000 -minPres 1000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs3000_mp1000/donors_$cross \
                    -o ../analysis/projection_hs3000_mp1000/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 3000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs3000_mp1000/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs3000_mp1000/*


# transform to diploid hapmap:
cd ../analysis/projection_hs3000_mp1000

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```

`scripts/count_projected_SVs.R`

> Percent projected was way worse than with `-hapSize 1000`



<mark>B73xPHG47 not in hapsize 2500!</mark>


<mark>hapsize 2000 seems the best!</mark>





### Summary of projected SVs


`scripts/count_projected_SVs.R`


Merge projected files into one with `scripts/merge_projected_crosses.R`. Created file `data/usda_SNPs-SVs_325rils.not-in-SVs.projected.hmp.txt`.

`scripts/plot_ril_karyotype_SVs.R`





### Test accuracy of projection

<mark> REDO THIS</mark>

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





## Resequencing SNPs

Resequencing SNPs from Mazaheri et al (2019), BMC Plant Biology
<https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-019-1653-x>
Data:
<https://datadryad.org/stash/landing/show?big=showme&id=doi%3A10.5061%2Fdryad.n0m260p>

Downloaded to `Downloads` folder and unzipped. Then unzip `62biomAP_v_B73_SNPMatrix.txt.txt.tar.gz` and it's has all the SNPs from 56 inbreds.
I will just select the SNPs from the 7 parents of my population (B73, LH82, PH207, PHG35, PHG39, PHG47, and PHJ40).

```bash
cd ~/Downloads/972727712

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
cut -f 1,2,5,21,32,36,37,38,43 62biomAP_v_B73_SNPMatrix.txt > ~/OneDrive/University\ of\ Minnesota/PhD/hirsch_lab/projects/genomic_prediction/simulation/data/resequencing_snps/biomAP_v_B73_SNPMatrix_7parents.txt
```

Resequencing SNPs from 7 parents are in the folder `tests/2020_02_20/data/resequencing_snps/` in the file called `biomAP_v_B73_SNPMatrix_7parents.txt`.

Transform it to hapmap format:

```bash
cd ~/OneDrive/University\ of\ Minnesota/PhD/hirsch_lab/projects/genomic_prediction/simulation
python scripts/snp_matrix2hapmap.py ~/Downloads/972727712/biomAP_v_B73_SNPMatrix_7parents.txt tests/2020_02_20/data/resequencing_snps/biomAP_v_B73_SNPMatrix_7parents.hmp.txt
```


Transform to hapmap diploid format:

```bash
cd ~/OneDrive/University\ of\ Minnesota/PhD/hirsch_lab/projects/genomic_prediction/simulation
run_pipeline.pl -Xmx8g -importGuess tests/2020_02_20/data/resequencing_snps/biomAP_v_B73_SNPMatrix_7parents.hmp.txt -export tests/2020_02_20/data/resequencing_snps/biomAP_v_B73_SNPMatrix_7parents.hmp.txt -exportType HapmapDiploid
```



Filter out SNPs inside deletion boundaries, keep only polymorphic:
Use donors for SVs as the anchors (SNPs only).

* Script `scripts/merge_reseq_SNPs_with_SNPs-SVs_hapmap.R`
* Created ...

Parental donor: `biomAP_parents_SNPs-reseq_and_SNPs-chip.B73xLH82.poly.hmp.txt`
Data to be projected: `biomAP_rils_SNPs-reseq_and_SNPs-chip.B73xLH82.poly.hmp.txt`





Projection of resequencing SNPs:
<mark>TO DO: find hapsize that matches number of blocks when resequencing SVs

```bash
# create directory to store results of imputation
mkdir -p tests/2020_02_20/analysis/projection_reseq_snps_

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
                    -o ../analysis/projection_hs2000_mp500/donors_$cross \
                    -hapSize 2000 -minTaxa 1
    # impute ril genotypes based on
    run_pipeline.pl -FILLINImputationPlugin \
                    -hmp merged_hapmaps_by_cross/usda_SNPs-SVs_$cross\_RILs.sorted.hmp.txt \
                    -d ../analysis/projection_hs2000_mp500/donors_$cross \
                    -o ../analysis/projection_hs2000_mp500/usda_SNPs-SVs_$cross\_RILs.projected.hmp.txt \
                    -hapSize 2000 -hybNN false -accuracy
  done
} < "usda_biparental-crosses.txt" > "../analysis/projection_hs2000_mp500/FILLIN_log.txt"
# set delimiter back to what it was
IFS=$curr_IFS

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir ../analysis/projection_hs2000_mp500/*


# transform to diploid hapmap:
cd ../analysis/projection_hs2000_mp500

for file in *"projected.hmp.txt"; do
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done
```




Count projected:

* `scripts/count_projected_SVs.R`
