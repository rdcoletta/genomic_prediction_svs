#!/usr/bin/python3

'''
created by Rafael Della Coletta
2019-08-15
'''


import argparse as ap
import pandas as pd


# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script reads a VCF file containing structural variant calls,
             and transform it into a numeric hapmap file.

example: vcf2hapmap.py my_file.vcf my_results.sorted.hmp.txt B73,LH82,PH207''')
# add positional arguments
parser.add_argument("vcf_file", type=str,
                    help="VCF file with SV calls")
parser.add_argument("output_name", type=str,
                    help="name of the hapmap file")
parser.add_argument("lines", type=str,
                    help="name of the lines to use (separated by comma)")
# pass arguments into variables
args = parser.parse_args()
vcf_file = args.vcf_file
output_name = args.output_name
inbreds_list = args.lines
inbreds_list = inbreds_list.split(",")


# open input file
infile = open(vcf_file, "r")

# open output file
outfile = open(output_name, "w")
# write output header
print("rs", "alleles", "chrom", "pos", "strand", "assembly", "center",
      "protLSID", "assayLSID", "panel", "QCcode", "\t".join(inbreds_list),
      sep="\t", file=outfile)

# list with important info about SVs to extract from vcf file
sv_info_list = ["CHROM", "POS", "ID", "ALT", "INFO"]

# create list to store indices of each info about sv and inbreds in the header
header_idx = []


# read file line by line
for line in infile:
    line = line.strip()
    # if line starts with ## -- skip
    if line[0:2] == "##":
        continue
    # if line starts with # -- save as header
    elif line[0:1] == "#":
        # save header
        header = line.split("\t")
        # remove # from #CHROM
        header[0] = header[0][1:]
        # get indices of columns with important information from header
        for info in sv_info_list:
            header_idx.append(header.index(info))
        for inbred in inbreds_list:
            header_idx.append(header.index(inbred))
        # print(header_idx)
    # otherwise, extract info from line
    else:
        line = line.split("\t")
        # get chrom number
        chr = line[header_idx[0]]
        # get length of sv
        sv_length = line[header_idx[4]].split("SVLEN=")
        sv_length = sv_length[1].split(";")
        sv_length = abs(int(sv_length[0]))
        # get position in the middle of the SV
        sv_start = int(line[header_idx[1]])
        sv_end = sv_start + sv_length
        pos = round((sv_start + sv_end) / 2)
        # determine type of sv
        if line[header_idx[3]].find("DEL") > -1:
            sv_type = "del"
        else:
            sv_type = "dup"
        # add sv type into id
        id = sv_type + "." + line[header_idx[2]]
        # parse each inbred line (based on its index on header)
        inbreds_geno = []
        for index in header_idx[5:]:
            # print(line[index])
            inbred_info = line[index].split(":")
            if inbred_info[0] == "1/1":
                genotype = "2"
            else:
                genotype = "0"
            inbreds_geno.append(genotype)
        # write output
        print(id, "NA", chr, pos, "NA", "NA", "NA", "NA", "NA", "NA", "NA",
              "\t".join(inbreds_geno), sep="\t", file=outfile)


# close files
infile.close()
outfile.close()


# since i was using the middle point of the SV as the its position in the
# hapmap file, it's possible that the position of the SVs are not ordered
# correctly (i.e., ascending order). Thus, i need to open previous output
# and sort SVs by chromosome and position


# open file as data frame
hapmap = pd.read_table(output_name, sep="\t", keep_default_na=False)

# sort by chromosome and then by position
hapmap_sorted = hapmap.sort_values(["chrom", "pos"])

# write sorted hapmap
hapmap_sorted.to_csv(output_name, sep="\t", index=False)
