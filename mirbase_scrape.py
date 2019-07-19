import re
import requests
from bs4 import BeautifulSoup


def get_accession_number(mirna_name):

    url = "http://www.mirbase.org/cgi-bin/textsearch_json.pl?q={}".format(mirna_name)

    r = requests.get(url)
    r.raise_for_status()

    ret = [data["accession"] for data in r.json()["data"] if "accession" in data]

    return ret[0]


def get_details(name, accession_number):

    url = "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc={}".format(accession_number)

    r = requests.get(url)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")

    ret = {
        "name": name,
        "accession": accession_number,
        "pre": "",
        "sequence": ""
    }

    text = soup.select("tr.sectionTitle td h2")[0].text
    if text is not None:
        ret["pre"] = re.sub(r'(Stem-loop sequence )|\n|\t', '', text)
    seq_url = "http://www.mirbase.org/cgi-bin/get_seq.pl?acc={}".format(accession_number)

    seq_r = requests.get(seq_url)
    seq_r.raise_for_status()

    ret["sequence"] = re.search(r'[ATGCU]{10,}', seq_r.text).group(0)

    return ret
