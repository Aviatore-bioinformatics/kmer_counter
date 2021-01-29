from shutil import which
import os
from app.text_formating import warning, ok

JELLYFISH_OUT_DIR = 'jellyfish'


def read_config(config_file):
    parameters = {}

    with open(config_file, 'r') as file:
        for line in file:
            line = line.rstrip()

            if line.startswith('#') or len(line) == 0:
                continue

            key, value = line.split('=')
            if key == 'prefixes':
                value = value.split(',')
            parameters[key] = value

    parameters['jellyfish_out_dir'] = os.path.join(parameters["output_dir"], JELLYFISH_OUT_DIR)

    return parameters


def is_jellyfish_installed():
    return which('jellyfish')


def test():
    print(f'Checking if jellyfish is installed ... ', end='')
    if is_jellyfish_installed():
        print(ok('ok'))
    else:
        print(warning('fail'))
        print('Install the jellyfish software.')
        return False

    return True
