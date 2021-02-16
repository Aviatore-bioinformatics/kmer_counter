import re
import numpy as np
import os
from app.text_formating import red, green, print_logo, print_info, print_warning


class TableMerger:
    def __init__(self, parameters):
        self.parameters = parameters
        self.output_path = os.path.join(self.parameters['output_dir'], 'tables')
        self.merged_data = {}

    def run(self):
        print_logo("Merging tables")

        if not self.check_run():
            return True

        try:
            if self.merge_tables():
                self.write_merged_tables()

            return True
        except Exception:
            return False

    def merge_tables(self):
        lineSplit = re.compile(r'\t')

        if (os.path.exists(os.path.join(self.output_path, 'table_merged.txt'))):
            os.remove(os.path.join(self.output_path, 'table_merged.txt'))

        tables = self.get_all_table_files()
        if len(tables) == 0:
            print_warning("There is no tables to merge")
            return False

        for table in tables:
            print_info(f"Reading {table} ...")

            with open(os.path.join(self.output_path, table), 'r') as f:
                while True:
                    line = f.readline()
                    if line == "":
                        break

                    line = line.rstrip()
                    line_splitted = lineSplit.split(line)
                    if line_splitted[0] == "k-mer":
                        self.merged_data["header"] = line
                        continue

                    line_splitted[1:] = list(map(int, line_splitted[1:]))
                    try:
                        self.merged_data[line_splitted[0]] += np.array(line_splitted[1:])
                    except KeyError:
                        self.merged_data[line_splitted[0]] = np.array(line_splitted[1:])

        print_info(f"Merging completed")

        return True

    def get_all_table_files(self):
        tables = []

        for file in os.listdir(self.output_path):
            if file.split("_")[-1] in self.parameters['prefixes']:
                tables.append(file)

        return tables

    def write_merged_tables(self):
        with open(os.path.join(self.output_path, "table_merged.txt"), 'a+') as f:
            f.write(self.merged_data.pop("header", None) + "\n")

            for kmer in sorted(self.merged_data.keys()):
                f.write(kmer + "\t" + "\t".join(list(map(str, self.merged_data[kmer]))) + "\n")

    def check_run(self):
        if not os.path.exists(os.path.join(self.parameters['output_dir'], 'tables', 'table_merged.txt')):
            return True

        if self.parameters['keep_kmers_merged_table'] == 'yes':
            print_info("Keeping merged table from the previous run")
            return False

        return True