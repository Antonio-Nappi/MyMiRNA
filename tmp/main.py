import shlex
import subprocess


def fastqc(filename):
    fastqc_command = "fastqc {}".format(filename)
    fastqc_process = subprocess.Popen(shlex.split(fastqc_command))
    _, error = fastqc_process.communicate()
    if error is not None:
        raise Exception("FASTQC error!")

def cutadapt(filename, adapter="TGGAATTCTCGGGTGCCAAGG"):
    in_file = filename
    out_file = "Norm_1_reduced_trimmed.fastq"
    command = "cutadapt -a {0} -o {1} {2} -j {3} -q {4} --discard-untrimmed -M {5} -m {6}".format(adapter,
                                                                                       out_file,
                                                                                       in_file,
                                                                                       8, 20, 35, 10)
    p = subprocess.Popen(shlex.split(command))
    out, err = p.communicate()
    print(out, err)
    fastqc(out_file)

def bowtie():
reference = ""
