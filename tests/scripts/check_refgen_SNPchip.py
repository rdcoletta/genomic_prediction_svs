from Bio import SeqIO

# add positions to search and marker genotypes in lists
pos_to_search = []
SNPchip_markers = []

file_pos = open("../data/markers_b73_chr1.txt", "r")
file_pos.readline()

for line in file_pos:
    line = line.strip()
    line = line.split("\t")
    pos = int(line[0])
    genotype = line[1][0]  # get just the first allele (since it's homozygous)
    pos_to_search.append(pos)
    SNPchip_markers.append(genotype)

file_pos.close()


# ref genome versions
ref_versions = ["v1", "v2", "v3", "v4"]
# create empty dictionary to store results
ref_markers = {}

# get marker genotype for each reference genome
for version in ref_versions:
    # read fasta file
    seq_record = SeqIO.read("../data/" + version + "_chr1.fasta", "fasta")
    # get sequence
    my_seq = seq_record.seq
    # create list to store nucletides for each position
    pos_list = []
    # extract positions from fasta file
    for pos in pos_to_search:
        if pos <= len(my_seq):
            pos_list.append(my_seq[pos - 1]) # remember to subtract 1 because python is 0-index
        else:
            pos_list.append("N")
    ref_markers[version] = pos_list

# open file to write output
outfile = open("../data/ref-markers_chr1.txt", "w")
# write header
print("pos", "marker", "\t".join(ref_versions), sep="\t", file=outfile)
# write output
for index in range(0, len(pos_to_search)):
    print(pos_to_search[index], SNPchip_markers[index],
          ref_markers["v1"][index], ref_markers["v2"][index],
          ref_markers["v3"][index], ref_markers["v4"][index],
          sep="\t", file=outfile)

outfile.close()


# calculate how many matches each refgen has to SNP chip markers
file_results = open("../data/ref-markers_chr1.txt", "r")
# skip header
file_results.readline()

# create dictionary with refgens as keys, and values = 0
matches_dict = {"v1":0, "v2":0, "v3":0, "v4":0}
# count total number of lines
total_lines = 0

for line in file_results:
    line = line.strip()
    line = line.split("\t")
    if line[1] == line[2]:
        matches_dict["v1"] += 1
    if line[1] == line[3]:
        matches_dict["v2"] += 1
    if line[1] == line[4]:
        matches_dict["v3"] += 1
    if line[1] == line[5]:
        matches_dict["v4"] += 1
    total_lines += 1

print("Number of matches to each SNP chip position by refgen version")
print(f"B73_v1: {matches_dict['v1']} ({round(matches_dict['v1'] * 100 / total_lines, ndigits=1)}%)")
print(f"B73_v2: {matches_dict['v2']} ({round(matches_dict['v2'] * 100 / total_lines, ndigits=1)}%)")
print(f"B73_v3: {matches_dict['v3']} ({round(matches_dict['v3'] * 100 / total_lines, ndigits=1)}%)")
print(f"B73_v4: {matches_dict['v4']} ({round(matches_dict['v4'] * 100 / total_lines, ndigits=1)}%)")

file_results.close()
