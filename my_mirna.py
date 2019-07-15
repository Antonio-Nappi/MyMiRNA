from utils import run_command
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi

def fastqc(filename):
    '''
     The method checks the quality of the reads in the pre-processing step.
     The tool used to compute the QC is FastQC
     :param filename: the input files to analyse
     :return: the output of the command
    '''
    fastqc_command = "fastqc {}".format(filename)
    return run_command(fastqc_command)


def cutadapt(filename, out_file = "Norm_1_reduced_trimmed.fastq", adapter="TGGAATTCTCGGGTGCCAAGG"):
    '''
     In the pre-processing step it is necessary to trim the adapters.
     The tool used is CutAdapt
     :param filename: the input file to analyse
     :param out_file: the output file with the results of the trimming
     :param adapter: the adapter sequence to remove from reads
    '''
    in_file = filename
    command = "cutadapt -a {0} -o {1} {2} -j {3} -q {4} --discard-untrimmed -M {5} -m {6}".format(adapter,
                                                                                       out_file,
                                                                                       in_file,
                                                                                       8, 20, 35, 10)
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

def mapping_tophat(in_file, ref_path, cores):
    '''
    This function allows to identify all the miRNAs by aligning the reads against mirBase.
    :param in_file: the input file with the reads that are going to be aligned
    :param ref_path: the path of the file about the reference file of mirBase
    :param cores: the number of cores to use for the alignment
    :return: the output of the command
    '''

    command = "tophat --bowtie1 -N 1 --library-type=fr-unstranded -p {0} {1} {2}".format(cores, ref_path, in_file)
    return run_command(command)

def feature_counts(in_file):
    '''
     :param in_file: the input file to analyse
     :return: the output of the coomand
    '''
    out_file = "out"
    command = "featureCounts -O -a hsa.gtf -o {0} {1}".format(out_file, in_file)
    return run_command(command)

def structure():
    cg = forgi.load_rna('Cluster_229_Y.txt',allow_many=false)
    fvm.plot_rna(cg, text_kwargs={"fontweight":"black"}, lighten=0.7,
             backbone_kwargs={"linewidth":3})
    plt.show()

structure()
