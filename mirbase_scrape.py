import requests


def get_accession_number(mirna_name):

    url = "http://www.mirbase.org/cgi-bin/textsearch_json.pl?q={}".format(mirna_name)

    r = requests.get(url)
    r.raise_for_status()

    ret = [for data in r.json()["data"]]

    return r.json()["data"][0]["accession"]


get_accession_number("hsa-miR-3605-5p")
