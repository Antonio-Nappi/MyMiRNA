from utils import run_command
#import matplotlib.pyplot as plt
 #import forgi.visual.mplotlib as fvm
#import forgi
import re

def fastqc(filename):
    '''
     The method checks the quality of the reads in the pre-processing step.
     The tool used to compute the QC is FastQC
     :param filename: the input files to analyse
     :return: the output of the command
    '''
    fastqc_command = "fastqc {}".format(filename)
    return run_command(fastqc_command)


def cutadapt(filename, out_file = "Norm_1_reduced_trimmed.fastq", adapter="TGGAATTCTCGGGTGCCAAGG",cores=8,quality=20):
    '''
     In the pre-processing step it is necessary to trim the adapters.
     The tool used is CutAdapt
     :param filename: the input file to analyse
     :param out_file: the output file with the results of the trimming
     :param adapter: the adapter sequence to remove from reads
     :param cores: the number of cores to perform the command
     :param quality: the minimum read quality
    '''
    in_file = filename
    command = "cutadapt -a {0} -o {1} {2} -j {3} -q {4} --discard-untrimmed -M {5} -m {6}".format(adapter,
                                                                                       out_file,
                                                                                       in_file,
                                                                                       cores, quality, 35, 10)
    out_file = run_command(command)
    print(out_file)
    fastqc(out_file)


def indexing(filename):
    '''
     This function allows to index the files to speed up the mapping process.
     The tool used is Bowtie
     :param filename: the input file to index
     :return: the output of the command
    '''
    in_file = filename
    out_file = filename+"indexed"
    command = "bowtie-build -f {0} {1}".format(in_file,out_file)
    return run_command(command)


def bam_to_fastq(filename):
    '''
     This function converts a bam file into a fastq file
     The tool used is Samtools
     :param filename: the input file to convert
     :return: the output of the command
    '''
    in_file = filename
    out_file = in_file[:len(in_file)-3]+"fastq"
    command = "samtools fastq {0} > {1}".format(in_file,out_file)
    return run_command(command)


def mapping_shortstack(in_file, ref_genome, n_core):
    '''
     This function allows to discover unique reads against the reference genome
     The tool used is ShortStack
     :param in_file: the input file that is going to be aligned against the reference genome
     :param ref_genome: the reference genome
     :param n_core: the number of cores to run for the analysis
     :return: the output of the command
    '''
    command = "ShortStack --readfile {0} --genome {1}".format(in_file, ref_genome)
    return run_command(command)

def aligning_bowtie(in_file, out_file, index_tag,noaligned_file, mismatches=1, cores=8):
    '''
     This function allows to align the trimmed reads within a reference.
     It is used to align the reads within mirBase, pirBase and Enseble (ncRNA)
    :param in_file: The input file to align
    :param out_file: the output file that contains the result
    :param mismatches: the maximum number of mismatches allowed
    :param cores: the number of cores to use for the command
    :param index_tag: the index used to build the ebwt files
    :param noaligned_file: the file that contains the unaligned reads
    :return: the output of the command
    '''
    command = "bowtie -v {0} -S -p {1} --un {2} {3} {4} {5}".format(mismatches, cores,noaligned_file,index_tag, in_file, out_file)
    return run_command(command)

def feature_counts(in_file,ref_file):
    '''
     :param in_file: the input file to analyse
     :param ref_file: the file with the reference annotation
     :return: the output of the coomand
    '''
    out_file = "out"
    command = "featureCounts -O -a {0} -o {1} {2}".format(ref_file,out_file, in_file)
    return run_command(command)

def structure():
    cg = forgi.load_rna('Cluster_229_Y.txt',allow_many = false)
    fvm.plot_rna(cg, text_kwargs={"fontweight":"black"}, lighten=0.7,
             backbone_kwargs={"linewidth":3})
    plt.show()

'''
Una volta eseguito ShortStack si indicizza nuovamente su ogni DB per avere una granularità più fine per i dati.
I passaggi successivi alla seconda indicizzazione sono comuni a tutti i percorsi.
Bisogna ottenere i FASTA da RFAM per indicizzare sui non coding e da REFSeq per i nuovi geni.
Una volta ottenuti i non mappati in nessuno dei database precedenti si effettua sempre la stessa analisi di espressione differenziale e target prediction
'''
