from utils import run_command

'''
 The method checks the quality of the reads in the pre-processing step.
 The tool used to compute the QC is FastQC.
'''
def fastqc(filename):
    fastqc_command = "fastqc {}".format(filename)
    return run_command(fastqc_command)

'''
 In the pre-processing step it is necessary to trim the adapters.
 The tool used is CutAdapt
'''
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

'''
 This function allows to index the files to speed up the mapping process.
 The tool used is Bowtie
'''
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

'''
 This function allows to discover unique reads against the reference genome
'''
def mapping_butter(in_file, ref_genome, n_core):
    command = "butter --mismatches {0} --aln_cores {1} --max_rep {2} {3} {4}".format(1, n_core, 5000, in_file, ref_genome)
    return run_command(command)

def mapping_tophat(in_file,mode,cores):
    #bisogna indicizzare il file hairpin e mature
    #il comando cambia a seconda della modalità
    command = "tophat --bowtie1 -N 1 --library-type=fr-unstrand -p {0} /home/antonio/Scrivania/idexed/mature indicizzato {1}".format(cores,in_file)
    return run_command(command)