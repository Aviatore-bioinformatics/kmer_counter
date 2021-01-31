import os
import scipy.stats as stats
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
from app.text_formating import warning, ok
# EXAMPLE: multipletests([0.01, 0.02, 0.03], method='bonferroni')
# RETURNS: (array([ True, False, False]), array([0.03, 0.06, 0.09]), 0.016952427508441503, 0.016666666666666666)


class Stat:
    def __init__(self, parameters):
        self.parameters = parameters
        self.BONFERRONI_OUT_INDEX = 1
        self.chrom_len = {}
        self.total_genome_len = 0
        self.mite_total_len = {}
        self.index = []
        self.mite_names = {}
        self.column_names = ["mite_total_len", "mite", "out", "freq", "fisher_exac_p"]
        self.data = pd.DataFrame(columns=self.column_names)
        self.merged_table_path = os.path.join(parameters['output_dir'], 'tables', 'table_merged.txt')
        self.output_path = os.path.join(parameters['output_dir'], 'stats', 'stats.txt')

        if not os.path.exists(os.path.join(parameters['output_dir'], 'stats')):
            os.mkdir(os.path.join(parameters['output_dir'], 'stats'))

            if not os.path.exists(os.path.join(parameters['output_dir'], 'stats')):
                print(f"{warning('Warning')} - could not create directory {os.path.join(parameters['output_dir'], 'stats')}")
                return

        if os.path.exists(os.path.join(parameters['output_dir'], 'stats', 'stats.txt')):
            os.remove(os.path.join(parameters['output_dir'], 'stats', 'stats.txt'))


    def progress_bar(self, current_value, final_value):
        offset = 100
        diff = (current_value / final_value) * 100

        if not current_value % offset:
            print(f'\r\033[0K{current_value} / {final_value} ({diff:.2f}%)', end='', flush=True)

    def run(self):
        if not os.path.exists(self.merged_table_path):
            print(
                f"{warning('Warning')} - the merged table does not exist")
            return

        self.chrom_len_calc()
        self.mite_total_len_calc()
        self.analyse()
        self.save_stats_to_file(self.stats_filtration())
        self.save_stats_to_file(self.data)

    def chrom_len_calc(self):
        for prefix in self.parameters['prefixes']:
            with open(os.path.join(self.parameters['data_dir'], f'{prefix}_oneLine.txt')) as file:
                line_len = 0

                for line in file:
                    line = line.rstrip()
                    line_len += len(line)

                self.chrom_len[prefix] = line_len
                self.total_genome_len += line_len

# chrom_len = {
#     "chr1": 51465340,
#     "chr2": 43913520,
#     "chr3": 50312079,
#     "chr4": 35924511,
#     "chr5": 41956025,
#     "chr6": 36610139,
#     "chr7": 36358036,
#     "chr8": 31745509,
#     "chr9": 33682890,
# }
# total_genome_len = 361968049
# normalization_size = 10000
# mite_total_len = {}

    def mite_total_len_calc(self):
        with open(self.parameters['bed_file'], 'r') as f:
            for line in f:
                line = line.rstrip()

                line_spliced = line.split('\t')

                try:
                    self.mite_total_len[line_spliced[-1]] += ( int(line_spliced[2]) - int(line_spliced[1]) + 1)
                except KeyError:
                    self.mite_total_len[line_spliced[-1]] = ( int(line_spliced[2]) - int(line_spliced[1]) + 1)

    def get_number_of_lines(self, file_name):
        with open(file_name, 'r') as file:
            for line_index, _ in enumerate(file):
                pass

        return line_index + 1

    def analyse(self):
        with open(self.merged_table_path, 'r') as f:
            for line_numer, line in enumerate(f):
                self.progress_bar(line_numer + 1, self.get_number_of_lines(self.merged_table_path))

                line = line.rstrip()

                kmers_in_mites_normalized = []
                kmers_out_mites_normalized = []
                kmers_in_mites_total_sum = 0
                mite_total_len_int = 0
                line_spliced = line.split('\t')

                if line_spliced[0] == "k-mer":
                    for i in range(2, len(line_spliced) - 2):
                        if line_spliced[i][-5:] != "_edge":
                            self.index.append(i)
                            self.mite_names[i] = line_spliced[i]
                else:
                    kmerName = line_spliced[0]
                    kmerTotalOccurences = int(line_spliced[1])
                    for i in self.index:
                        if int(line_spliced[i]) > 0:
                            kmers_in_mites_total_sum += int(line_spliced[i])
                            mite_total_len_int += self.mite_total_len[self.mite_names[i]]

                    if mite_total_len_int > 0:

                        # kmers_in_mites_total_norm = (kmers_in_mites_total_sum / mite_total_len_int) * self.parameters['normalization_size']
                        # kmers_out_mites_total_norm = ((kmerTotalOccurences - kmers_in_mites_total_sum) /
                        #     (self.total_genome_len - mite_total_len_int)) * self.parameters['normalization_size']
                        a = round( (mite_total_len_int / self.total_genome_len) * kmerTotalOccurences )
                        b = round( ((self.total_genome_len - mite_total_len_int) / self.total_genome_len) * kmerTotalOccurences )

                        kmers_in_out_total_sum = kmerTotalOccurences - kmers_in_mites_total_sum

                        # print(kmerTotalOccurences, kmers_in_mites_total_sum)

                        _, p = stats.fisher_exact([[kmers_in_mites_total_sum, kmers_in_out_total_sum], [a, b]])

                        freq = kmers_in_mites_total_sum / (kmers_in_mites_total_sum + kmers_in_out_total_sum)

                        row = pd.Series([mite_total_len_int, kmers_in_mites_total_sum, kmers_in_out_total_sum, freq, p],
                                        name=kmerName, index=self.column_names)
                        self.data = self.data.append(row)


        corrected_p = multipletests(list(self.data['fisher_exac_p']), method='bonferroni')[self.BONFERRONI_OUT_INDEX]

        self.data['p_corrected_bon'] = corrected_p

    def stats_filtration(self):
        # TODO First filtration step: All data with 'p_corrected_bon' > 0.05 must be filtered out
        # TODO Second filtration step: All data wich freq is out of provided threshold range must be filtered out
        if self.parameters['kmer_thresh_min'] != '':
            data_filtered_min = self.data.loc[self.data['p_corrected_bon'] > int(self.parameters['kmer_thresh_min'])]

        if self.parameters['kmer_thresh_max'] != '':
            if self.parameters['kmer_thresh_min'] != '':
                data_filtered = data_filtered_min.loc[
                    self.data['p_corrected_bon'] <= int(self.parameters['kmer_thresh_max'])]
            else:
                data_filtered = self.data.loc[
                    self.data['p_corrected_bon'] <= int(self.parameters['kmer_thresh_max'])]

        if self.parameters['kmer_thresh_min'] == '' and self.parameters['kmer_thresh_max'] == '':
            data_filtered = self.data

        return data_filtered

    def save_stats_to_file(self, data):
        data.to_csv(os.path.join(self.parameters['output_dir'], 'stats', 'stats.txt'), sep='\t')
