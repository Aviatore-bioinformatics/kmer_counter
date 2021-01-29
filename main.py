from app.test_environment import read_config, test
from app.kmer_counter import jellyfish
from app.utils import fasta_to_oneline


def main():
    parameters = read_config('config.txt')

    if not test():
        return False

    # if not jellyfish(parameters):
    #     return False
    #
    # if not fasta_to_oneline(parameters):
    #     return False


if __name__ == '__main__':
    main()
