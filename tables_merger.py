import re
import numpy as np
import os

lineSplit = re.compile(r'\t')

files = [
    "table_chr1",
    "table_chr2",
    "table_chr3",
    "table_chr4",
    "table_chr5",
    "table_chr6",
    "table_chr7",
    "table_chr8",
    "table_chr9",
]

merged_data = {}

if os.path.exists("table_merged"):
    print("Output table exists. Removing ...")
    os.remove("table_merged")

for file in files:
    print("Reading", file, "...")
    with open(file, 'r') as f:
        while True:
            line = f.readline()
            if line is "":
                break

            line = line.rstrip()
            line_splitted = lineSplit.split(line)
            if line_splitted[0] == "k-mer":
                merged_data["header"] = line
                continue
            
            line_splitted[1:] = list(map(int, line_splitted[1:]))
            try:
                merged_data[ line_splitted[0] ] += np.array(line_splitted[1:])
            except(KeyError):
                merged_data[ line_splitted[0] ] = np.array(line_splitted[1:])

with open("table_merged", 'a+') as f:
    f.write(merged_data.pop("header", None) + "\n")
    
    for kmer in sorted(merged_data.keys()):
        f.write(kmer + "\t" + "\t".join( list(map(str, merged_data[kmer])) ) + "\n")
