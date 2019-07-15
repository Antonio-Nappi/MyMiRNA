def split_params(string_params):
    params = string_params.split(';')
    ret = dict()
    for param in params:
        keyval = param.split('=')
        ret[keyval[0]] = keyval[1]
    return ret


def gff_to_saf(in_filename, out_filename, filter_primary_trainscript=True):
    with open(in_filename) as in_file:
        lines = in_file.readlines()

    i = 0
    while lines[i][0] == '#':
        i += 1

    final = []

    while i < len(lines):
        values = lines[i][:-1].split('\t')

        i += 1

        if filter_primary_trainscript and values[2] == "miRNA_primary_transcript":
            continue

        final.append({
            'Chr': values[0],
            'Start': values[3],
            'End': values[4],
            'Strand': values[6],
            'GeneID': split_params(values[8])['Name']
        })

    with open(out_filename+".saf", "w") as out_file:
        out_file.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for line in final:
            out_file.write("{}\t{}\t{}\t{}\t{}\n".format(line["GeneID"],
                                                          line["Chr"],
                                                          line["Start"],
                                                          line["End"],
                                                          line["Strand"]))


gff_to_saf("assets/hsa.gff3", "out")
