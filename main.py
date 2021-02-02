from app.test_environment import read_config, test
from app.jellyfish_controller import jellyfish
from app.utils import fasta_to_oneline
import app.kmer_counter as kmer
from app.tables_merger import TableMerger
from app.stats_pandas import Stat
from app.fasta_to_oneline_controller import bulk_fasta_to_online


def main():
    parameters = read_config('config.txt')

    if not test():
        return False

    if not jellyfish(parameters):
        return False

    if not bulk_fasta_to_online(parameters):
        return False

    kc = kmer.KmerCounter(parameters)
    if not kc.run():
        return False

    tm = TableMerger(parameters)
    if not tm.run():
        return False

    stat = Stat(parameters)
    if not stat.run():
        return False


if __name__ == '__main__':
    main()
