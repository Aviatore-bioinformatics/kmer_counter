from shutil import which
import os
import pkg_resources
from app.text_formating import warning, ok

JELLYFISH_OUT_DIR = 'jellyfish'

PIP_REQUIREMENTS = [
    "intervaltree"
]

SOFT_REQUIREMENTS = [
    "jellyfish"
]


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

    if not os.path.exists(parameters["output_dir"]):
        os.mkdir(parameters["output_dir"])

    parameters['jellyfish_out_dir'] = os.path.join(parameters["output_dir"], JELLYFISH_OUT_DIR)

    if not os.path.exists(parameters['jellyfish_out_dir']):
        os.mkdir(parameters['jellyfish_out_dir'])

    return parameters


def soft_check():
    output = True
    print(f'Checking required software ...')

    for soft in SOFT_REQUIREMENTS:
        print(f'- {soft} ... ', end='')
        if which(soft):
            print(ok('ok'))
        else:
            print(warning('fail'))
            print(f'You need to install {soft}.')
            output = False

    return output


def pip_check():
    output = True
    installed = {pkg.key for pkg in pkg_resources.working_set}

    print("Checking required python packages ...")

    for package in PIP_REQUIREMENTS:
        print(f'- {package} ... ', end='')

        if package in installed:
            print(ok('ok'))
        else:
            print(warning('fail'))
            print(f'You need to install {package}.')
            output = False

    return output


def test():
    pip_check()
    soft_check()

    return True
