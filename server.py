from flask import Flask, jsonify, render_template, request, Response, send_from_directory
from flask_cors import CORS
from mirbase_scrape import get_accession_number, get_details
from my_mirna import aligning_bowtie, bam_to_fastq, cutadapt, fastqc, feature_counts, mapping_shortstack, multiqc
import os
import re

app = Flask(__name__)
CORS(app)


@app.route('/trimming', methods=['POST'])
def quality_and_trimming():

    # TODO label samples for differential analysis
    # TODO generate index dinamically
    index = "prova1"

    params = request.get_json(silent=True, cache=False)

    adapter = params["adapter"]

    quality = params["quality"]

    # TODO retrieve cores from request
    cores = 8

    fastqc_folder = index + "/fastqc"
    trimmed_folder = index + "/trimmed"

    input_data_folder = "data"

    filenames = os.listdir(input_data_folder)

    # Creates folders to store the fastqc output file and the trimmed fastq files
    os.mkdir(index)
    os.mkdir(fastqc_folder)
    os.mkdir(trimmed_folder)

    # First fastqc run
    for filename in filenames:
        fastqc(input_data_folder + "/" + filename, fastqc_folder)

    # The adapter is trimmed and the file are saved as original_name_trimmed.fastq
    trimmed_filenames = []
    for filename in filenames:
        trimmed_filename = re.sub(r'\.fastq', '', filename) + "Trimmed.fastq"
        trimmed_filenames.append(trimmed_filename)
        cutadapt(input_data_folder + "/" + filename, trimmed_folder+"/"+trimmed_filename, adapter, quality, cores)

    # Second fastqc run
    for filename in trimmed_filenames:
        fastqc(trimmed_folder+"/"+filename, fastqc_folder)

    # All the fastqc results are reported in a single multiqc file
    multiqc(fastqc_folder, fastqc_folder)

    multiqc_path = "/assets/multiqc_report.html"

    return multiqc_path


@app.route("/shortstack", methods=['POST'])
def shortstack_mapping():
    # TODO retrieve from request
    multimap_threshold = 500

    # TODO retrieve from request
    n_cores = 8

    # TODO retrieve
    index = "prova"

    # Creates a folder for the shortstack output file
    base_dir = index+"/shortstack"
    os.mkdir(base_dir)

    # Maps each trimmed file in the trimmed folder and converts the bam output to fastq
    shortstack_dirs = []
    for filename in os.listdir(index+"/trimmed"):
        tmp = re.sub(r'\.fastq', '', filename)
        tmp_dir = base_dir+"/"+tmp
        os.mkdir(tmp_dir)
        shortstack_dirs.append(tmp_dir)
        mapping_shortstack(index+"/trimmed/"+filename, "assets/human_index/genome.fa", tmp_dir, multimap_threshold, n_cores)
        # bam_to_fastq(tmp_dir + "/" + filename[:-5] + "bam")

    # Exports data about the mapping
    template_parameters = {
        "unique": [],
        "mm": [],
        "mm_ign": [],
        "nm": []
    }
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

    return render_template("templates/shortstack.html", params=template_parameters)


@app.route("/mirna", methods=['POST'])
def mirna():
    # TODO retrieve from request
    n_mismatch = 0

    # TODO retrieve from request
    n_cores = 8

    # TODO retrieve
    index = "prova"

    base_folder = index + "/mirna"
    os.mkdir(base_folder)

    shortstack_folder = index + "/shorstack"

    in_files = []
    for folder in os.listdir(shortstack_folder):
        in_files.append(shortstack_folder + "/" + folder + "/" + folder + ".bam")

    os.mkdir(base_folder+"/featurecounts")

    with open(base_folder+"/featurecounts/coldata.tsv") as coldata_file:
        coldata_file.write("Sample\tCondition\n")
        for filename in in_files:
            type = filename.split('/')[-1][0]
            name = re.sub('/', '.', filename)
            coldata_file.write("{}\t{}\n".format(name, type))

    feature_counts(in_files, "assets/mirna.saf", base_folder+"/featurecounts/results")

    # TODO differential analysis

    # TODO return grafici (in template)



@app.route("/mirnas", methods=['GET'])
def mirnas():

    # TODO retrieve index
    index = "prova"

    with open(index+"/mirna/mirna.names") as mirna_names_file:
        mirna_names = [re.sub(r'\n', '', name) for name in mirna_names_file.readlines()]

    ret = {
        "mature": [],
        "pre": []
    }

    for name in mirna_names:
        if re.search(r'^[a-zA-Z]{2,3}-[a-zA-Z]{2,3}-[0-9]{3,5}$', name) is None:
            ret["mature"].append(name)
        else:
            ret["pre"].append(name)

    return jsonify(ret)

@app.route("/mirnas/<string: mirna>", method=['GET'])
def get_mirna(mirna):

    accession = get_accession_number(mirna)[0]

    params = get_details(mirna, accession)



if __name__ == "__main__":
    app.run('127.0.0.1', 8080)
