import os
import scipy.stats as stats
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
# EXAMPLE: multipletests([0.01, 0.02, 0.03], method='bonferroni')
# RETURNS: (array([ True, False, False]), array([0.03, 0.06, 0.09]), 0.016952427508441503, 0.016666666666666666)

BONFERRONI_OUT_INDEX = 1


def progress_bar(current_value, final_value):
    offset = 100
    diff = (current_value / final_value) * 100

    if not current_value % offset:
        print(f'\r\033[0K{current_value} / {final_value} ({diff:.2f}%)', end='', flush=True)

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
with open("/media/veracrypt1/Lenovo/Dysk/2020/Projektt_Ali_MP/analizy/MITEs.bed", 'r') as f:
    for line in f:
        line = line.rstrip()

        line_spliced = line.split('\t')

        try:
            mite_total_len[line_spliced[-1]] += ( int(line_spliced[2]) - int(line_spliced[1]) + 1)
        except KeyError:
            mite_total_len[line_spliced[-1]] = ( int(line_spliced[2]) - int(line_spliced[1]) + 1)


index = []
if os.path.exists("stats_total_v2.txt"):
    os.remove("stats_total_v2.txt")

mite_names = {}

column_names = ["mite_total_len", "mite", "out", "freq", "fisher_exac_p"]
data = pd.DataFrame(columns=column_names)

print('Getting number of lines ... ', end='', flush=True)
with open("/media/veracrypt1/Lenovo/Dysk/2020/Projektt_Ali_MP/analizy/table_merged_10000_lines", 'r') as f:
    for line_index, _ in enumerate(f):
        pass
number_of_lines = line_index + 1
print(f'{number_of_lines} lines')


with open("/media/veracrypt1/Lenovo/Dysk/2020/Projektt_Ali_MP/analizy/table_merged_10000_lines", 'r') as f:
    for line_numer, line in enumerate(f):
        progress_bar(line_numer + 1, number_of_lines)

        line = line.rstrip()

        kmers_in_mites_normalized = []
        kmers_out_mites_normalized = []
        kmers_in_mites_total_sum = 0
        mite_total_len_int = 0
        line_spliced = line.split('\t')

        if line_spliced[0] == "k-mer":
            for i in range(2, len(line_spliced) - 2):
                if line_spliced[i][-5:] != "_edge":
                    index.append(i)
                    mite_names[i] = line_spliced[i]
        else:
            kmerName = line_spliced[0]
            kmerTotalOccurences = int(line_spliced[1])
            for i in index:
                if int(line_spliced[i]) > 0:
                    kmers_in_mites_total_sum += int(line_spliced[i])
                    mite_total_len_int += mite_total_len[mite_names[i]]

            if mite_total_len_int > 0:
                kmers_in_mites_total_norm = (kmers_in_mites_total_sum / mite_total_len_int) * normalization_size
                kmers_out_mites_total_norm = ((kmerTotalOccurences - kmers_in_mites_total_sum) /
                    (total_genome_len - mite_total_len_int)) * normalization_size
                a = round( (mite_total_len_int / total_genome_len) * kmerTotalOccurences )
                b = round( ((total_genome_len - mite_total_len_int) / total_genome_len) * kmerTotalOccurences )

                kmers_in_out_total_sum = kmerTotalOccurences - kmers_in_mites_total_sum
                _, p = stats.fisher_exact([[kmers_in_mites_total_sum, kmers_in_out_total_sum], [a, b]])

                freq = kmers_in_mites_total_sum / (kmers_in_mites_total_sum + kmers_in_out_total_sum)

                row = pd.Series([mite_total_len_int, kmers_in_mites_total_sum, kmers_in_out_total_sum, freq, p],
                                name=kmerName, index=column_names)
                data = data.append(row)

corrected_p = multipletests(list(data['fisher_exac_p']), method='bonferroni')[BONFERRONI_OUT_INDEX]

data['p_corrected_bon'] = corrected_p

data.to_csv('stats_total_v2.txt', sep='\t')
