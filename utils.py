import shlex
import subprocess
import re


def run_command(command):
    '''
    This function runs a command using subprocess
    :param command: the command to run
    :return: the output value
    '''

    process = subprocess.Popen(shlex.split(command))
    out, err = process.communicate()

    if err is not None:
        raise Exception(err)

    return out


def miRNA_filter(filename):
    '''
     This function filters all the reads in the Unplaced.txt (the outcome of ShortStask).
     The applied filter removes all the reads that are already mapped as miRNAs.
     The output file contains all the reads that are not miRNAs.
     :param filename: the input file to filter
     :return: the filterd output file
    '''
    out_file = open("filtered_miRNA.txt", "w")
    f = open(filename) #default rt mode
    line = f.readlines()
    for x in line:
        value = x.split('\t')[4][:-1]
        if value != '>50':
            out_file.write(x)
    return out_file


def filtering(filename):
    '''
    Utilizzata per filtrare il file di Giorgio sui ncRNA. Sono stati eliminati miRNA, piRNA e lncRNA
    :param filename:
    :return:
    '''
    out_file = open("C:\\Users\Alex\Desktop\\filtered.fa", "w")
    f = open(filename)  # default rt mode
    line = f.readlines()
    stringa = ""
    for x in line:
        stringa += x
    separate = stringa.split(">")
    for x in separate:
        if re.search("gene_biotype:miRNA", x) is None:
            out_file.write(">" + x)


def split_params(string_params):
    params = string_params.split(';')
    ret = dict()
    for param in params:
        keyval = param.split('=')
        ret[keyval[0]] = keyval[1]
    return ret


def filter_dashr(in_file, out_file):
    with open(in_file) as inf:
        lines = inf.readlines()

    accepted = {'tRF5', 'tRF3', 'snoRNA', 'rRNA', 'tRNA', 'snRNA', 'scRNA'}

    with open(out_file, "w") as outf:
        outf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for line in lines:
            vals = line.split('\t')
            if vals[2] in accepted:
                g_id = split_params(vals[8])["ID"]
                chr = vals[0]
                start = vals[3]
                end = vals[4]
                strand = vals[6]
                s = '\t'.join([g_id, chr, start, end, strand]) + '\n'
                outf.write(s)


filter_dashr("assets/dashr.gff", "assets/sncRNA.saf")
