from utils import run_command


def fastqc(filename):
    '''
     The method checks the quality of the reads in the pre-processing step.
     The tool used to compute the QC is FastQC.
    '''
    fastqc_command = "fastqc {}".format(filename)
    return run_command(fastqc_command)


def cutadapt(filename, out_file = "Norm_1_reduced_trimmed.fastq", adapter="TGGAATTCTCGGGTGCCAAGG"):
    '''
     In the pre-processing step it is necessary to trim the adapters.
     The tool used is CutAdapt
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
    '''
    in_file = filename
    out_file = filename+"indexed"
    command = "bowtie-build -f {0} {1}".format(in_file,out_file)
    return run_command(command)


def bam_to_fastq(filename):
    '''
     This function converts a bam file into a fastq file
    '''
    in_file = filename
    out_file = in_file[:len(in_file)-3]+"fastq"
    command = "samtools fastq {0} > {1}".format(in_file,out_file)
    return run_command(command)


def mapping_shortstack(in_file, ref_genome, n_core):
    '''
     This function allows to discover unique reads against the reference genome
    '''
    command = "ShortStack --readfile {0} --genome {1}".format(in_file, ref_genome)
    return run_command(command)


def feature_counts(in_file):
    out_file = "out"
    command = "featureCounts -a hsa.gtf -o {0} {1}".format(out_file, in_file)
    return run_command(command)
