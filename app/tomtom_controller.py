import os
import subprocess
from app.text_formating import red, green, print_info, print_warning, print_logo


class Tomtom:
    def __init__(self, parameters):
        self.parameters = parameters
        self.stats_file_path = os.path.join(self.parameters['output_dir'], 'stats', 'stats.txt')

        if not os.path.exists(os.path.join(self.parameters['output_dir'], 'tomtom')):
            os.mkdir(os.path.join(self.parameters['output_dir'], 'tomtom'))

        self.output_meme_file_path = os.path.join(self.parameters['output_dir'], 'tomtom', "kmers.meme")

    def run(self):
        print_logo("K-mer comparing")

        if self.parameters['run_tomtom'] == 'no':
            print_info("Analysis canceled. To run tomtom set 'run_tomtom' parameter to 'yes'")
            return True

        self.kmers_to_meme()
        self.tomtom()

        return True

    def kmers_to_meme(self):
        """Converts kmer sequences into MEME format and saves them to the 'kmers.meme' file"""

        if os.path.exists(self.output_meme_file_path):
            os.remove(self.output_meme_file_path)

        print_info(f"Converting kmers into MEME format ... ")

        with open(self.output_meme_file_path, 'a+') as output:
            with open(self.stats_file_path, 'r') as file:
                for line in file:
                    line = line.rstrip()
                    line_splitted = line.split("\t")

                    if len(line_splitted[0]) > 0:
                        result = subprocess.run(['iupac2meme', line_splitted[0]], capture_output=True, text=True)
                        output.write(result.stdout)

    def tomtom(self):
        """Compares kmer motifs with database using tomtom"""
        print_info(f"Comparing kmer motifs with database using tomtom ... ")

        parameters = ['tomtom']

        # parameters.append('-min-overlap')
        parameters.append('-min-overlap')
        parameters.append(self.parameters['min_overlap'])

        if self.parameters['internal'] == 'yes':
            parameters.append('-internal')

        if self.parameters['threshold_type'] == 'e-value':
            parameters.append('-evalue')

        parameters.append('-thresh')
        parameters.append(self.parameters['threshold_value'])

        parameters.append('-oc')
        parameters.append(os.path.join(self.parameters['output_dir'], 'tomtom', 'tomtom_out'))

        parameters.append(self.output_meme_file_path)
        parameters.append(self.parameters['motif_database'])

        # subprocess.run(parameters, capture_output=True, text=True)
        process = subprocess.Popen(parameters, stderr=subprocess.PIPE)
        while True:
            output = process.stderr.readline().decode('utf-8').rstrip()

            if output == '' and process.poll() is not None:
                break

            if output != '':
                try:
                    if output[0] == 'P':
                        print(f"\r\033[0K{output}", end='', flush=True)
                except IndexError:
                    print(f"DEBUG: {output}")
                    raise IndexError

        print("")

        if not process.returncode:
            print_info(f"Processing completed")
            print_info(f"Report in HTML format is available at: ./output/tomtom/tomtom_out/tomtom.html")
        else:
            print_warning("Something went wrong with tomtom run. Used command:")
            print(" ".join(parameters))

        print(" ".join(parameters))



