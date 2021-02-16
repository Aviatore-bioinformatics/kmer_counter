from app.test_environment import read_config, test

def main():
    parameters = read_config('config.txt')

    if not test():
        return False

    from app.jellyfish_controller import jellyfish
    import app.kmer_counter as kmer
    from app.tables_merger import TableMerger
    from app.stats_pandas import Stat
    from app.fasta_to_oneline_controller import bulk_fasta_to_oneline
    from app.tomtom_controller import Tomtom

    if not jellyfish(parameters):
        return False

    if not bulk_fasta_to_oneline(parameters):
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

    tomtom = Tomtom(parameters)
    if not tomtom.run():
        return False


if __name__ == '__main__':
    main()
