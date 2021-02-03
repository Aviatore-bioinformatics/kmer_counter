import os
from app.utils import fasta_to_oneline
from app.text_formating import red, green, print_warning, print_info, print_logo


def bulk_fasta_to_oneline(parameters):
    print('')
    print_info('Converting FASTA files to text files.')

    for file_prefix in parameters['prefixes']:
        input_file = f'{file_prefix}.fasta'
        input_file_path = os.path.join(parameters['data_dir'], input_file)
        output_file = f'{file_prefix}_oneLine.txt'
        output_file_path = os.path.join(parameters['data_dir'], f'{file_prefix}_oneLine.txt')

        if os.path.exists(output_file_path):
            print_info(f'The output {output_file} file already exists. Skipping ...')
            continue

        print_info(f'Converting {input_file} into {output_file} ... ')

        try:
            fasta_to_oneline(input_file_path, output_file_path)
        except Exception as e:
            print_warning(f'Something went wrong during saving to the {output_file} file.')
            print_warning('Please, check the stderr output:\n')
            print(e)

            return False

    print_info("Conversion completed")

    return True
