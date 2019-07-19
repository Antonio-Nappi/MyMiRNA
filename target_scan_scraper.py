from bs4 import BeautifulSoup
import re
import requests

url = "http://www.targetscan.org/cgi-bin/targetscan/vert_72/targetscan.cgi"
base_url = "http://www.targetscan.org"


def target_scan_get_table(name):
    url = "http://www.targetscan.org/vert_71/temp/TargetScan7.1__{}.predicted_targets.txt".format(re.sub(r'hsa-', '', name))

    r = requests.get(url)

    if r.status_code != 200:
        return "<table></table>"

    table = "<table>"
    rows = r.text.split('\n')
    for row in rows:
        table += "<tr>"
        for col in row.split('\t'):
            table += "<td>" + col + "</td>"
        table += "</tr>"
    table += "</table>"
    return table

def target_scan_search(gid, species="Human", **kwargs):
    '''
    This function scrapes TargetScan. If more than one parameter among mir_sc, mir_c, mir_nc, mir_vcn, mirg
    is not empty, just the first entered is considered
    :param gid: the gene symbol / Ensembl gene / transcript ID
    :param species: the species (default Human)
    :param kwargs: one (and just one) additional key among mir_sc, mir_c, mir_nc, mir_vcn or mirg
    :return: A dict of URLs
    '''

    optional_params = ["mir_sc", "mir_c", "mir_nc", "mir_vcn", "mirg"]
    payload = {
        "gid": gid,
        "species": species
    }

    for key in optional_params:
        if key in kwargs:
            payload[key] = kwargs[key]
            break

    r = requests.get(url, params=payload)

    if r.status_code != 200:
        return dict()

    soup = BeautifulSoup(r.text, 'html.parser')

    check = soup.select_one("body[onload]")

    if check is None:

        links = soup.select("a[href^='/cgi-bin']")

        if len(links) > 0:
            ret = dict()
            ret["representative"] = base_url + links[0].get("href")
            ret["others"] = []

            for link in links[1:]:
                ret["others"].append(base_url+link.get("href"))

            return ret

        return dict()

    else:
        ret = {
            "representative": re.sub("'", '',re.sub(r"document\.location='", '', check.get("onload")))
        }
        return ret

print(target_scan_get_table("hsa-miR-3605-5p"))
