#! /usr/bin/env python3

import requests
import sys
import time

def get_taxa(infile):


def uniprot_search(gene, taxID):
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A{taxID}%29+AND+%28gene%3A{gene}%29%29"
    return requests.get(url).text




