#! /usr/bin/env python3

from Bio import SeqIO
import argparse

def remove_partial_seq(infile, outfasta, outremove):
    removed_headers = []
    infile = open(infile, "r")
    outfasta = open(outfasta, "w")
    kept_count = 0
    removed_count = 0

    for record in SeqIO.parse(infile, "fasta"):
        if "partial" not in record.description.lower():
            SeqIO.write(record, outfasta, "fasta")
            kept_count += 1
        else:
            removed_count += 1
            removed_headers.append(record.description)


    infile.close()
    outfasta.close()
    outremove = open(outremove, "w")

    for header in removed_headers:
        outremove.write(header + "\n")

    outremove.close()
    print(f"Complete seqeuences kept: {kept_count}")
    print(f"Partial sequences removed: {removed_count}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
            type = str,
            help = "Input FASTA file")
    parser.add_argument("-o", "--output",
            type = str,
            help = "Output FASTA file after removing partials")
    parser.add_argument("-l", "--list",
            type = str,
            help = "Output text file containing list of removed partials")
    args = parser.parse_args()
    remove_partial_seq(args.input, args.output, args.list)

if __name__=="__main__":
    main()
