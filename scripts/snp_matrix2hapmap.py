#!/usr/bin/python3

'''
created by Rafael Della Coletta
2020-02-14
'''


import argparse as ap


# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script reads a tab-delimeted file containing resequencing
             SNPs, and transform it into a hapmap file.

example: snp_matrix2hapmap.py my_file.txt my_results.hmp.txt''')
# add positional arguments
parser.add_argument("matrix_file", type=str,
                    help="SNP matrix file")
parser.add_argument("output_name", type=str,
                    help="name of the hapmap file")
# pass arguments into variables
args = parser.parse_args()
matrix_file = args.matrix_file
output_name = args.output_name


# open input file
infile = open(matrix_file, "r")

# open output file
outfile = open(output_name, "w")


# get header information
header = infile.readline()
header = header.strip()
header = header.split("\t")
# get inbred list
inbreds_list = header[1:]

# write output header
print("rs", "alleles", "chrom", "pos", "strand", "assembly",
      "center", "protLSID", "assayLSID", "panel", "QCcode",
      "\t".join(inbreds_list), sep="\t", file=outfile)


# read file line by line
for line in infile:
    line = line.strip()
    line = line.split("\t")
    # get chrom number and SNP position
    chr = line[0].split(":")[0]
    pos = line[0].split(":")[1]
    # create id based on sv type and location
    id = "snp." + chr + "." + pos
    # parse each inbred line (based on its index on header)
    inbreds_geno = []
    for index in range(1, len(header)):
        # print(line[index])
        genotype = line[index]
        inbreds_geno.append(genotype)
    # check which alleles are present for that SV
    unique_alleles = list(set(inbreds_geno))
    # remove N, if in the list
    unique_alleles = [allele for allele in unique_alleles if allele != 'N']
    unique_alleles = "/".join(unique_alleles)
    # write output
    print(id, unique_alleles, chr, pos, "NA", "NA", "NA", "NA", "NA", "NA", "NA",
          "\t".join(inbreds_geno), sep="\t", file=outfile)


# close files
infile.close()
outfile.close()
