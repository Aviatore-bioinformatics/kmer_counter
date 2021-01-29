import os
from app.text_formating import warning, ok


def fasta_to_oneline(parameters):
    print('')
    print('Converting FASTA files to text files.')

    for file_prefix in parameters['prefixes']:
        input_file = f'{file_prefix}.fasta'
        output_file = f'{file_prefix}_oneLine.txt'

        if os.path.exists(output_file):
            print(f'The output {output_file} file already exists. Skipping ...')
            continue

        print('')
        print(f'Converting {input_file} into {output_file} ... ', end='', flush=True)

        with open(input_file, 'r') as file_in:
            try:
                with open(output_file, 'w') as file_out:
                    output = ""

                    for line in file_in:
                        line = line.rstrip()

                        if line.startswith('>'):
                            continue

                        output += line

                    output += '\n'

                    file_out.write(output)
            except Exception as e:
                print(warning('fail'))
                print(f'Something went wrong during saving to the {output_file} file.')
                print('Please, check the stderr output:\n')
                print(e)

                return False

        print(ok('ok'))

    return True
