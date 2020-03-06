#!/usr/bin/python3

'''
created by Rafael Della Coletta
2019-07-02
'''

import argparse as ap

# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script reads in a table with two columns containing the
             genotype names and IDs of RILs genotyped for the USDA project,
             and outputs a new table containing the IDs of lines that make up
             a biparental cross.''')
# add positional arguments
parser.add_argument("id_table", type=str,
                    help="table with genotype names and IDs")
parser.add_argument("output_name", type=str,
                    help="name of the output table")
# pass arguments into variables
args = parser.parse_args()
id_table = args.id_table
output_name = args.output_name


# open input file
infile = open(id_table, "r")

# skip header
infile.readline()

# create dictionary to store biparental crosses as keys and RIL IDs as values
biparental_crosses = {}

# parse file
for line in infile:
    line = line.strip()
    line = line.split("\t")
    cross = line[0].split("-")[0]
    ril_ID = line[1]
    if cross not in biparental_crosses:
        # introduce new key-value pair
        biparental_crosses[cross] = ril_ID
    else:
        # update the correspondent value by adding 1 to what's there
        biparental_crosses[cross] += ","+ril_ID

# open output file
outfile = open(output_name, "w")

# print header
print("cross", "RILs", sep="\t", file=outfile)

# write table
for key, val in biparental_crosses.items():
    print(key, val, sep="\t", file=outfile)


# close files
infile.close()
outfile.close()
