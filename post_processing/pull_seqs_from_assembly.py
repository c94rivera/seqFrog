# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 14:27:49 2014

@author: RDT
adapted by: CAR
"""

from Bio import SeqIO
import csv
import sys
from tqdm import tqdm

def get_seqs(seqfile, matches, outfile):
    """Open fasta infile and return iterator of SeqRecords with protein sequences."""
    records = SeqIO.parse(seqfile, "fasta")
    seqlist=[]
    matchlist=[]
    with open(matches,"rU") as f:
        reader = csv.reader(f,delimiter="\t")
        for row in reader:
            matchlist.append(row)
        print("Blast file successfully read")
    print("Comparing contigs to Blast file")

    proteinid = []
    for rec in tqdm(records):
        proteinid.clear()
        for match in matchlist:
            if match[0] == rec.id:
                if match[1] not in proteinid:
                    proteinid.append(match[1])
                    rec2 = rec[:]

                    rec2.description += "\t" + match[1] + "\t" + match[2] + "\t" + match[3] + "\t" + match[10] + "\t" + match[13]

                    seqlist.append(rec2)
                    continue
                else:
                    continue
            else:
                continue


    print("Writing annotated contigs")
    with open(outfile, "w") as f:
        SeqIO.write(seqlist, f, "fasta")

def main():
    assert len(sys.argv) == 3, "usage: python3 pullseqs_fromtrinity.py <sequences.fasta> <blast_hits.txt>"
    infile = sys.argv[1]
    print("infile is", infile)
    blast_hits = sys.argv[2]
    print("blast file is", blast_hits)
    outfile = blast_hits[:blast_hits.index(".")]+".fasta"
    get_seqs(infile, blast_hits, outfile)
    print ("outfile is", outfile)

main()
