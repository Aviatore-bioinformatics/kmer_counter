from app.test_environment import read_config, test
from app.jellyfish_controller import jellyfish
from app.utils import fasta_to_oneline
import app.kmer_counter as kmer


def main():
    parameters = read_config('config.txt')

    if not test():
        return False

    if not jellyfish(parameters):
        return False

    if not fasta_to_oneline(parameters):
        return False

    kc = kmer.KmerCounter(parameters)
    # kc.run()

if __name__ == '__main__':
    main()
