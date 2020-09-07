# Multitasking using Pool and map

import os
from intervaltree import Interval, IntervalTree
import time
from multiprocessing import Pool
import copy


class Timer:
    def __init__(self):
        self.start = None
        self.stop = None
    
    def startt(self):
        self.start = time.time()

    def stopp(self):
        self.stop = time.time()
        print("\n", self.stop - self.start, "sek")


def my_find(str1, str2):
    start = 0
    end = len(str1)
    all_finds = []
    while str1.find(str2, start, end) != -1:
        x = str1.find(str2, start, end)
        all_finds.append(str(x))
        start = x+1

    return all_finds


data_inputs = [
    {
        "dump_file": "chr1_dump.fasta",
        "chr_file": "chr1_oneLine_fasta.txt",
        "chr_name": "chr1",
        "output_file": "table_chr1"
    },
    {
        "dump_file": "chr2_dump.fasta",
        "chr_file": "chr2_oneLine_fasta.txt",
        "chr_name": "chr2",
        "output_file": "table_chr2"
    },
    {
        "dump_file": "chr3_dump.fasta",
        "chr_file": "chr3_oneLine_fasta.txt",
        "chr_name": "chr3",
        "output_file": "table_chr3"
    },
    {
        "dump_file": "chr4_dump.fasta",
        "chr_file": "chr4_oneLine_fasta.txt",
        "chr_name": "chr4",
        "output_file": "table_chr4"
    },
    {
        "dump_file": "chr5_dump.fasta",
        "chr_file": "chr5_oneLine_fasta.txt",
        "chr_name": "chr5",
        "output_file": "table_chr5"
    },
    {
        "dump_file": "chr6_dump.fasta",
        "chr_file": "chr6_oneLine_fasta.txt",
        "chr_name": "chr6",
        "output_file": "table_chr6"
    },
    {
        "dump_file": "chr7_dump.fasta",
        "chr_file": "chr7_oneLine_fasta.txt",
        "chr_name": "chr7",
        "output_file": "table_chr7"
    },
    {
        "dump_file": "chr8_dump.fasta",
        "chr_file": "chr8_oneLine_fasta.txt",
        "chr_name": "chr8",
        "output_file": "table_chr8"
    },
    {
        "dump_file": "chr9_dump.fasta",
        "chr_file": "chr9_oneLine_fasta.txt",
        "chr_name": "chr9",
        "output_file": "table_chr9"
    },
]
max_thread_No = 4

def worker_controller(q):
    while True:
        input = q.get()
        if input is None:
            break
        worker(input)
        q.task_done()

def worker(data_input):
    print("Loading '{}' file ...".format( data_input["dump_file"] ))
    data_kmer = {}
    name_tmp = ""
    with open(data_input["dump_file"], 'r') as file:
        cont = True
        while cont:
            line = file.readline()
            if line == '':
                cont = False
                break
            line = line.rstrip()
            if line[0] == ">":
                name_tmp = line[1:]
            else:
                data_kmer[line] = name_tmp
    
    print("Loading '{}' file ...".format( data_input["chr_file"] ))
    with open(data_input["chr_file"], 'r') as f:
        chromosome = f.read()
    
    #----------------------------------------#
    print("Loading 'MITEs.bed' file ...")
    data_mites = []
    with open("MITEs.bed", 'r') as f:
        while True:
            line = f.readline()
            if line == '':
                break
            line = line.rstrip()
            data_mites.append(line.split("\t"))
    #----------------------------------------#

    if os.path.exists( data_input["output_file"] ):
        print("The file '{}' exists. Removing ...".format(data_input["output_file"]))
        os.remove(data_input["output_file"])

    output = open(data_input["output_file"], 'a+')

    output_data_template = {}
    output_data_template["edge"] = 0
    output_data_template["genome"] = 0
    mite_names = []

    t = IntervalTree()
    for mite in data_mites:
        if mite[3] not in output_data_template.keys():
            output_data_template[mite[3]] = 0
            output_data_template[mite[3] + "_edge"] = 0
            mite_names.append(mite[3])
            mite_names.append(mite[3] + "_edge")
        if mite[0] == data_input["chr_name"]:
            t[int(mite[1]) - 0:int(mite[2])] = mite[3]
    mite_names = set(mite_names)

    output.write("\t".join(["k-mer", "total_occurences_in_{}".format(data_input["chr_name"]), "\t".join(sorted(mite_names)), "edge", "genome"]))
    output.write("\n")

    ttt = Timer()
    print("Started analysis ...")
    log = open(time.strftime('%y-%m-%d_%H-%M_') + data_input["chr_name"] + "_log.txt", 'a+')
    log.write("Analysis started at " + time.ctime() + "\n")
    log.flush()
    ttt.startt()
    kmer_No = 0
    for kmer in data_kmer.keys():

        output_data = copy.deepcopy(output_data_template)

        kmer_occurences = my_find(chromosome, kmer)

        for kmer_occurence in kmer_occurences:
            kmer_occurence = int(kmer_occurence)
            result = t[kmer_occurence + 1:kmer_occurence + 11]
            if result:
                if len( list(result) ) > 2:
                    print("\nThe interval tree length is higher than 2:", len( list(result) ), data_input["chr_name"])
                    log.write("\t".join(["intTree>2", kmer, str(kmer_occurence), mite_name]) + "\n")
                    log.flush()
                    log.close()
                    exit(1)
                elif len( list(result) ) > 1:
                    print("\nThe interval tree length is higher than 1.", data_input["chr_name"])
                    log.write("\t".join(["1<intTree<2", kmer, str(kmer_occurence), str(list(result))]) + "\n")
                    log.flush()

                    for interval in result:
                        output_data["edge"] += 1
                        output_data[interval.data + "_edge"] += 1
                else:
                    result_parsed = list(result)[0]
                    if (kmer_occurence) >= (result_parsed.begin - 1) and \
                        (kmer_occurence + 10) <= result_parsed.end:
                        output_data[result_parsed.data] += 1
                    else:
                        output_data["edge"] += 1
                        output_data[result_parsed.data + "_edge"] += 1
            else:
                output_data["genome"] += 1

        kmer_No += 1

        output.write("\t".join([kmer, str(len(kmer_occurences))]))
        for mite_name in sorted(mite_names):
            output.write("\t" + str(output_data[mite_name]))
        output.write("\t" + str(output_data["edge"]))
        output.write("\t" + str(output_data["genome"]))
        output.write("\n")
    log.close()
    output.close()
    ttt.stopp()

p = Pool(max_thread_No)
p.map(worker, data_inputs)