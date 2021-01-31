def fasta_to_oneline(input_file_path, output_file_path):
    with open(input_file_path, 'r') as file_in:
        with open(output_file_path, 'w') as file_out:
            output = ""

            for line in file_in:
                line = line.rstrip()

                if line.startswith('>'):
                    continue

                output += line

            output += '\n'

            file_out.write(output)
