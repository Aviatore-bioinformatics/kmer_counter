from shutil import which
import os
import pkg_resources
from app.text_formating import red, green

JELLYFISH_OUT_DIR = 'jellyfish'

PIP_REQUIREMENTS = [
    "intervaltree",
    "scipy",
    "pandas",
    "statsmodels"
]

SOFT_REQUIREMENTS = [
    "jellyfish",
    "iupac2meme",
    "tomtom"
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
        print(f'- {soft} ... ', end='', flush=True)
        if which(soft):
            print(green('ok'))
        else:
            print(red('fail'), end='', flush=True)
            print(f' (You need to install {soft})')
            output = False

    return output


def pip_check():
    output = True
    installed = {pkg.key for pkg in pkg_resources.working_set}

    print("Checking required python packages ...")

    for package in PIP_REQUIREMENTS:
        print(f'- {package} ... ', end='', flush=True)

        if package in installed:
            print(green('ok'))
        else:
            print(red('fail'), end='', flush=True)
            print(f' (You need to install {package})')
            output = False

    return output


def test():
    pip_output = pip_check()
    soft_output = soft_check()

    return pip_output is True and soft_output is True
