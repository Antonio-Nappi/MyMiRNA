from utils import run_command


def fastqc(filename):
    fastqc_command = "fastqc {}".format(filename)
    return run_command(fastqc_command)


def cutadapt(filename, adapter="TGGAATTCTCGGGTGCCAAGG"):
    in_file = filename
    out_file = "Norm_1_reduced_trimmed.fastq"
    command = "cutadapt -a {0} -o {1} {2} -j {3} -q {4} --discard-untrimmed -M {5} -m {6}".format(adapter,
                                                                                       out_file,
                                                                                       in_file,
                                                                                       8, 20, 35, 10)
    out_file = run_command(command)
    print(out_file)
    fastqc(out_file)

def indexing(filename):
    in_file = filename
    out_file = filename+"indexed"
    command = "bowtie-build -f {0} {1}".format(in_file,out_file)
    return run_command(command)
'''
This function converts a bam file into a fastq file
'''
def bam_to_fastq(filename):
    in_file = filename
    out_file = in_file[:len(in_file)-3]+"fastq"
    command = "samtools fastq {0} > {1}".format(in_file,out_file)
    return run_command(command)