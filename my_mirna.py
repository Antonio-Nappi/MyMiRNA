from utils import run_command
import shlex
import subprocess

import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi


def fastqc(in_file, out_dir='.'):
    '''
     The method checks the quality of the reads in the pre-processing step.
     The tool used to compute the QC is FastQC
     :param in_file: the input files to analyse
     :return: the output of the command
    '''
    fastqc_command = "fastqc -o {0} {1}".format(out_dir, in_file)
    return run_command(fastqc_command)


def cutadapt(in_file, out_file="Norm_1_reduced_trimmed.fastq", adapter="TGGAATTCTCGGGTGCCAAGG", cores=8, quality=20):
    '''
     In the pre-processing step it is necessary to trim the adapters.
     The tool used is CutAdapt
     :param in_file: the input file to analyse
     :param out_file: the output file with the results of the trimming
     :param adapter: the adapter sequence to remove from reads
     :param cores: the number of cores to perform the command
     :param quality: the minimum read quality
    '''
    command = "cutadapt -a {0} -j {3} -q {4} -m {5} -o {1} {2}".format(adapter,
                                                                       out_file,
                                                                       in_file,
                                                                       cores,
                                                                       quality,
                                                                       10)
    print(command)
    run_command(command)


def multiqc(in_dir, out_dir="."):
    command = "multiqc -o {0} {1}".format(out_dir, in_dir)
    run_command(command)


def indexing(in_file):
    '''
     This function allows to index the files to speed up the mapping process.
     The tool used is Bowtie
     :param in_file: the input file to index
     :return: the output of the command
    '''
    out_file = in_file + "indexed"
    command = "bowtie-build -f {0} {1}".format(in_file, out_file)
    return run_command(command)


def bam_to_fastq(in_file):
    '''
     This function converts a bam file into a fastq file
     The tool used is Samtools
     :param in_file: the input file to convert
     :return: the output of the command
    '''
    out_file = in_file[:-3] + "fastq"
    command = "samtools fastq {0}".format(in_file)
    with open(out_file, "w") as file:
        subprocess.call(shlex.split(command), stdout=file)

    return


def mapping_shortstack(in_file, ref_genome, out_dir, multimap_number=500, n_core=6):
    '''
     This function allows to discover unique reads against the reference genome
     The tool used is ShortStack
     :param in_file: the input file that is going to be aligned against the reference genome
     :param ref_genome: the reference genome
     :param out_dir: the output directory
     :param multimap_number: the maximum number of multimapped reads allowed
     :param n_core: the number of cores to run for the analysis
     :return: the output of the command
    '''
    command = "ShortStack --outdir {4} --bowtie_m {0} --bowtie_cores {1} --readfile {2} --genome {3}".format(
        multimap_number, n_core, in_file, ref_genome, out_dir)
    return run_command(command)


def aligning_bowtie(in_file, out_file, index_tag, noaligned_file, mismatches=1, cores=8):
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
    command = "bowtie -v {0} -S -p {1} --max /dev/null --un {2} {3} {4} {5}".format(mismatches, cores, noaligned_file,
                                                                                    index_tag, in_file, out_file)
    return run_command(command)


def feature_counts(in_files, ref_file, out_file):
    '''
     :param in_files: the input file to analyse
     :param ref_file: the file with the reference annotation
     :param out_file: the output file
     :return: the output of the command
    '''
    param = ""
    for file in in_files:
        param += file + " "
    command = "featureCounts -O -a {0} -F SAF -o {1} {2}".format(ref_file, out_file, param)
    return run_command(command)


def structure():
    cg = forgi.load_rna('Cluster_229_Y.txt', allow_many=false)
    fvm.plot_rna(cg, text_kwargs={"fontweight": "black"}, lighten=0.7,
                 backbone_kwargs={"linewidth": 3})
    plt.show()


def novel_pirna(in_file, out_file, species=4):
    '''
     This function allows to discover novel piRNAs.
     The tool used is piRNN.
    :param in_file: The input file should be small RNA data in fasta format.
    :param out_file: he output file is also in fasta format.
    :param species: the code of the specie you want to analyze: 1 means C. elegans; 2 means D. melanogaster; 3 means Rat 4 means Human.
    :return: the output of the command.
    '''
    command = "python3 piRNN.py -s {0} -i {1} -o {2}".format(species, in_file, out_file)
    return run_command(command)


def differential_analysis(in_file, index, coldata, **kwargs):
    p_filter = kwargs.get("p_value", "none")
    p_adj_filter = kwargs.get("p_value_adjusted", "none")
    log_2_fold = kwargs.get("log", "none")

    command = "Rscript diff_exp.r {0} {1} {2} {3} {4} {5}".format(in_file,
                                                                  coldata,
                                                                  p_filter,
                                                                  p_adj_filter,
                                                                  log_2_fold,
                                                                  index)
    run_command(command)

    ret = [
        {
            "name": "Boxplot of relative log expression",
            "path": "assets/{}/boxplot.jpeg".format(index)
        },
        {
            "name": "Log2Fold change of all data",
            "path": "assets/{}/spread_all.jpeg".format(index)
        },
        {
            "name": "Log2Fold change of filtered data",
            "path": "assets/{}/spread_filtered.jpeg".format(index)
        },
        {
            "name": "p-value histogram of all data",
            "path": "assets/{}/hist_p_all.jpeg".format(index)
        },
        {
            "name": "p-value histogram of filtered data",
            "path": "assets/{}/hist_p_filtered.jpeg".format(index)
        },
        {
            "name": "Adjusted p-value histogram of all data",
            "path": "assets/{}/hist_padj_all.jpeg".format(index)
        },
        {
            "name": "Adjusted p-value histogram of filtered data",
            "path": "assets/{}/hist_padj_filtered.jpeg".format(index)
        },
        {
            "name": "Heatmap of the filtered data",
            "path": "assets/{}/heatmap.jpeg".format(index)
        },
        {
            "name": "Barplot of the expression of the filtered data",
            "path": "assets/{}/barplot.jpeg".format(index)
        },
        {
            "name": "Volcano plot that relates log2fold change and adjusted p-value",
            "path": "assets/{}/volcano.jpeg".format(index)
        }
    ]

    return ret

'''
Una volta eseguito ShortStack si indicizza nuovamente su ogni DB per avere una granularità più fine per i dati.
I passaggi successivi alla seconda indicizzazione sono comuni a tutti i percorsi.
Bisogna ottenere i FASTA da RFAM per indicizzare sui non coding e da REFSeq per i nuovi geni.
Una volta ottenuti i non mappati in nessuno dei database precedenti si effettua sempre la stessa analisi di espressione differenziale e target prediction
    '''
'''mirnovo è un'interfaccia web, nuovo scraping??'''