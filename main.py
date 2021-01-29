from app.test_environment import read_config, test
from app.jellyfish_controller import jellyfish
from app.utils import fasta_to_oneline
import app.kmer_counter as kmer
from app.tables_merger import TableMerger
from app.stats_pandas import Stat


def main():
    parameters = read_config('config.txt')

    if not test():
        return False

    if not jellyfish(parameters):
        return False

    if not fasta_to_oneline(parameters):
        return False

    # kc = kmer.KmerCounter(parameters)
    # kc.run()
    #
    # tm = TableMerger(parameters)
    # tm.run()

    stat = Stat(parameters)
    stat.run()



if __name__ == '__main__':
    main()
