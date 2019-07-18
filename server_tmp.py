from flask import Flask, render_template
from my_mirna import aligning_bowtie, bam_to_fastq, cutadapt, fastqc, feature_counts, mapping_shortstack, multiqc
import os
import re


def quality_and_trimming():
    # TODO label samples for differential analysis
    # TODO generate index dinamically
    index = "prova"

    # TODO retrieve from request
    adapter = "TGGAATTCTCGGGTGCCAAGG"

    # TODO retrieve quality from request
    quality = 20

    # TODO retrieve cores from request
    cores = 8

    fastqc_folder = index + "/fastqc"
    trimmed_folder = index + "/trimmed"

    # TODO retrieve from request
    filenames = ["data/N1.fastq", "data/T1.fastq"]

    # Creates folders to store the fastqc output file and the trimmed fastq files
    os.mkdir(index)
    os.mkdir(fastqc_folder)
    os.mkdir(trimmed_folder)

    # First run of fastqc
    for filename in filenames:
        fastqc(filename, fastqc_folder)

    # The adapter is trimmed and the file are saved as original_nameTrimmed.fastq
    trimmed_filenames = []
    for filename in filenames:
        trimmed_filename = re.sub(r'\.fastq', '', filename.split("/")[-1]) + "Trimmed.fastq"
        trimmed_filenames.append(trimmed_filename)
        cutadapt(filename, trimmed_folder+"/"+trimmed_filename, adapter, quality, cores)

    # Second fastqc run
    for filename in trimmed_filenames:
        fastqc(trimmed_folder+"/"+filename, fastqc_folder)

    # All the fastqc results are reported in a single multiqc file
    multiqc(fastqc_folder, fastqc_folder)

    return render_template(fastqc_folder+"/multiqc_report.html")


def shortstack_mapping():
    # TODO retrieve from request
    multimap_threshold = 500

    # TODO retrieve from request
    n_cores = 8

    # TODO retrieve
    index = "prova"

    # Creates a folder for the shortstack output file
    base_dir = index+"/shortstack"
    #os.mkdir(base_dir)

    # Maps each trimmed file in the trimmed folder and converts the bam output to fastq
    shortstack_dirs = []
    for filename in os.listdir(index+"/trimmed"):
        tmp = re.sub(r'\.fastq', '', filename)
        tmp_dir = base_dir+"/"+tmp
        # os.mkdir(tmp_dir)
        shortstack_dirs.append(tmp_dir)
        #mapping_shortstack(index+"/trimmed/"+filename, "assets/human_index/genome.fa", tmp_dir, multimap_threshold, n_cores)
        bam_to_fastq(tmp_dir + "/" + filename[:-5] + "bam")

    template_parameters = {
        "unique": [],
        "mm": [],
        "mm_ign": [],
        "nm": []
    }

    # Exports data about the mapping
    for log_dir in shortstack_dirs:
        with open(log_dir+"/Log.txt") as file:
            text = file.read()

            match = re.search(r'Unique mappers: ([0-9]{5,} / [0-9]{5,}) \(([0-9]{,3}\.[0-9] %)\)', text)
            if match is None:
                template_parameters["unique"].append(("None", "None"))
            else:
                template_parameters["unique"].append((match.group(1), match.group(2)))

            match = re.search(r'Multi mappers: ([0-9]{5,} / [0-9]{5,}) \(([0-9]{,3}\.[0-9] %)\)', text)
            if match is None:
                template_parameters["mm"].append(("None", "None"))
            else:
                template_parameters["mm"].append((match.group(1), match.group(2)))

            match = re.search(r'Multi mappers ignored and marked as unmapped: ([0-9]{5,} / [0-9]{5,}) \(([0-9]{,3}\.[0-9] %)\)', text)
            if match is None:
                template_parameters["mm_ign"].append(("None", "None"))
            else:
                template_parameters["mm_ign"].append((match.group(1), match.group(2)))

            match = re.search(r'Non mappers: ([0-9]{5,} / [0-9]{5,}) \(([0-9]{,3}\.[0-9] %)\)', text)
            if match is None:
                template_parameters["nm"].append(("None", "None"))
            else:
                template_parameters["nm"].append((match.group(1), match.group(2)))

    print(template_parameters)
    return render_template("templates/shortstack.html", params=template_parameters)


def mirna():
    # TODO retrieve from request
    n_mismatch = 0

    # TODO retrieve from request
    n_cores = 8

    # TODO retrieve
    index = "prova"

    base_folder = index + "/mirna"
    os.mkdir(base_folder)

    shortstack_folder = index + "/shortstack"

    in_files = []
    for folder in os.listdir(shortstack_folder):
        tmp_dir = base_folder + "/" + folder
        os.mkdir(tmp_dir)

        out_file = tmp_dir + "/" + folder + ".sam"
        in_files.append(out_file)

        aligning_bowtie(shortstack_folder + "/" + folder + "/" + folder + ".fastq",
                        out_file,
                        "assets/mirna_index/hairpin",
                        tmp_dir + "/" + folder + "_unaligned.fastq",
                        n_mismatch,
                        n_cores)

    os.mkdir(base_folder+"/featurecounts")
    feature_counts(in_files, "assets/mirna.saf", base_folder+"/featurecounts/results")


mirna()
