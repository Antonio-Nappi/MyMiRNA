from flask import render_template
from my_mirna import differential_analysis, feature_counts
import os
import re


def main_branch(rna_type, index, params):
    base_folder = index + "/" + rna_type
    os.mkdir(base_folder)

    shortstack_folder = index + "/shortstack"

    in_files = []
    for folder in os.listdir(shortstack_folder):
        in_files.append(shortstack_folder + "/" + folder + "/" + folder + ".bam")

    os.mkdir(base_folder+"/featurecounts")

    with open(base_folder+"/featurecounts/coldata.tsv", "w") as coldata_file:
        coldata_file.write("Sample\tCondition\n")
        for filename in in_files:
            sample_type = filename.split('/')[-1][0]
            name = re.sub('/', '.', filename)
            coldata_file.write("{}\t{}\n".format(name, sample_type))

    feature_counts(in_files, "assets/{}.saf".format(rna_type), base_folder+"/featurecounts/results")

    graphics = differential_analysis(base_folder+"/featurecounts/results",
                                     index,
                                     base_folder+"/featurecounts/coldata.tsv",
                                     rna_type,
                                     params)

    rendered = render_template("dearesults.html", params=graphics)

    with open("gui/src/assets/{}_diff_exp.html".format(rna_type), 'w') as file:
        file.write(rendered)

    return "assets/{}_diff_exp.html".format(rna_type)

