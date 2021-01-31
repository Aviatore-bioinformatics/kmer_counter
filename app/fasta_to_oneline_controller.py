import os
from app.utils import fasta_to_oneline
from app.text_formating import warning, ok


def bulk_fasta_to_online(parameters):
    print('')
    print('Converting FASTA files to text files.')

    for file_prefix in parameters['prefixes']:
        input_file = f'{file_prefix}.fasta'
        input_file_path = os.path.join(parameters['data_dir'], input_file)
        output_file = f'{file_prefix}_oneLine.txt'
        output_file_path = os.path.join(parameters['data_dir'], f'{file_prefix}_oneLine.txt')

        if os.path.exists(output_file_path):
            print(f'The output {output_file} file already exists. Skipping ...')
            continue

        print(f'Converting {input_file} into {output_file} ... ', end='', flush=True)

        try:
            fasta_to_oneline(input_file_path, output_file_path)
        except Exception as e:
            print(warning('fail'))
            print(f'Something went wrong during saving to the {output_file} file.')
            print('Please, check the stderr output:\n')
            print(e)

            return False

    print(ok('ok'))

    return True
