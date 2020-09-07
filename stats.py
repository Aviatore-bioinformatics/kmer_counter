import re
import os
import scipy.stats as stats


r = re.compile(r'\t')

chrom_len = {
    "chr1": 51465340,
    "chr2": 43913520,
    "chr3": 50312079,
    "chr4": 35924511,
    "chr5": 41956025,
    "chr6": 36610139,
    "chr7": 36358036,
    "chr8": 31745509,
    "chr9": 33682890,
}
total_genome_len = 361968049
normalization_size = 10000
mite_total_len = {}
with open("MITEs.bed", 'r') as f:
    while True:
        line = f.readline()
        line = line.rstrip()
        if line == "":
            break

        line_spliced = r.split(line)
        try:
            mite_total_len[line_spliced[-1]] += ( int(line_spliced[2]) - int(line_spliced[1]) + 1)
        except KeyError:
            mite_total_len[line_spliced[-1]] = ( int(line_spliced[2]) - int(line_spliced[1]) + 1)


index = []
if os.path.exists("stats_total_v2.txt"):
    os.remove("stats_total_v2.txt")
outFile = open("stats_total_v2.txt", 'a+')
mite_names = {}
mite_names_forPrint = []
with open("table_merged", 'r') as f:
    while True:
        line = f.readline()
        line = line.rstrip()
        if line == "":
            break
        kmers_in_mites_normalized = []
        kmers_out_mites_normalized = []
        kmers_in_mites_total_sum = 0
        mite_total_len_int = 0
        line_spliced = r.split(line)
        if line_spliced[0] == "k-mer":
            for i in range(2, len(line_spliced) - 2):
                if line_spliced[i][-5:] != "_edge":
                    index.append(i)
                    mite_names[i] = line_spliced[i]
                    mite_names_forPrint.append(line_spliced[i] + "_norm")
                    mite_names_forPrint.append(line_spliced[i] + "_out_norm")
            outFile.write("\t".join(["kmer", "mite_total_len", "mite", "out", "fisher_exac_p"]) + "\n")
        else:
            kmerName = line_spliced[0]
            kmerTotalOccurences = int(line_spliced[1])
            for i in index:
                if int(line_spliced[i]) > 0:
                    kmers_in_mites_total_sum += int(line_spliced[i])
                    mite_total_len_int += mite_total_len[mite_names[i]]
            if mite_total_len_int > 0:
                kmers_in_mites_total_norm = (kmers_in_mites_total_sum / mite_total_len_int) * normalization_size
                kmers_out_mites_total_norm = ((kmerTotalOccurences - kmers_in_mites_total_sum) /\
                    (total_genome_len - mite_total_len_int)) * normalization_size
                a = round( (mite_total_len_int / total_genome_len) * kmerTotalOccurences )
                b = round( ((total_genome_len - mite_total_len_int) / total_genome_len) * kmerTotalOccurences )

                kmers_in_out_total_sum = kmerTotalOccurences - kmers_in_mites_total_sum
                _, p = stats.fisher_exact([[kmers_in_mites_total_sum, kmers_in_out_total_sum], [a, b]])

                outFile.write("\t".join([kmerName, str(mite_total_len_int), str(kmers_in_mites_total_sum),\
                    str(kmers_in_out_total_sum), str(p)]) + "\n")

outFile.close()