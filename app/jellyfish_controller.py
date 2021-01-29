import subprocess
import os
from app.text_formating import warning, ok


def check(result):
    try:
        result.check_returncode()
        print(ok('ok'))
    except subprocess.CalledProcessError:
        print(warning('fail'))

        print('Something went wrong during k-mer counting.')
        print('Please, check the stderr output:\n')
        print(result.stderr)

        return False

    return True


def kmer_counting(fasta_file, jellyfish_file, parameters):
    print('')
    print(f'Counting k-mers in the {fasta_file} file ... ', end='', flush=True)

    result = subprocess.run(['jellyfish', 'count',
                             '-m', parameters['kmer_length'],
                             '-s', parameters['hash_size'],
                             '-t', parameters['threads_number'],
                             '-C', fasta_file,
                             '-o', jellyfish_file], capture_output=True, text=True)

    return check(result)


def dump_jf_file(output_file_full_path, jellyfish_file_full_path, jellyfish_file, output_file_name):
    print(f'Outputting counts from the {jellyfish_file} file to the {output_file_name} file ... ', end='', flush=True)

    result = subprocess.run(['jellyfish', 'dump', jellyfish_file_full_path,
                             '-o', output_file_full_path], capture_output=True, text=True)

    return check(result)


def remove_jf_file(jellyfish_file, parameters):
    if parameters['keep_intermediate_jf_files'] == 'no':
        print(f'Deleting the {jellyfish_file} file ... ', end='', flush=True)

        try:
            os.remove(jellyfish_file)
            print(ok('ok'))
        except FileNotFoundError:
            print(warning('fail'))
            print(f'The {jellyfish_file} file was not found.')


def jellyfish(parameters):
    print('Start k-mer counting using jellyfish.')

    for file_prefix in parameters['prefixes']:
        fasta_file = f'{file_prefix}.fasta'
        fasta_file_full_path = os.path.join(parameters['data_dir'], fasta_file)

        jellyfish_file = f'{file_prefix}.jf'
        jellyfish_file_full_path = os.path.join(parameters['jellyfish_out_dir'], jellyfish_file)

        output_file = f'{file_prefix}_dump.fasta'
        output_file_full_path = os.path.join(parameters['jellyfish_out_dir'], output_file)

        if os.path.exists(output_file_full_path):
            print(f'The output {output_file} file already exists. Skipping ...')
            continue

        if not kmer_counting(fasta_file_full_path, jellyfish_file_full_path, parameters):
            return False

        if not dump_jf_file(output_file_full_path, jellyfish_file_full_path, jellyfish_file, output_file):
            return False

        remove_jf_file(jellyfish_file_full_path, parameters)

    return True
