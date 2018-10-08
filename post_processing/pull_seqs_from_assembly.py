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

def get_seqs(seqfile, matches, outfile, comparecolumn = None, keepcolumn = None):
    """Open fasta infile and return iterator of SeqRecords with protein sequences."""
    if not comparecolumn:
        comparecolumn = 1
    if not keepcolumn:
        keepcolumn = (1, 2, 3, 10, 13)
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
                if match[comparecolumn] not in proteinid:
                    proteinid.append(match[comparecolumn])
                    rec2 = rec[:]

                    for column in keepcolumn:
                        rec2.description += "\t" + match[int(column)]

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
    #assert len(sys.argv) == 3, "usage: python3 pullseqs_fromtrinity.py <sequences.fasta> <blast_hits.txt>"
    infile = sys.argv[1]
    print("infile is", infile)
    blast_hits = sys.argv[2]
    print("blast file is", blast_hits)
    outfile = blast_hits[:blast_hits.index(".")]+".fasta"
    comparecolumn = []
    keepcolumn = []
    if sys.argv[3]:
        comparecolumn = int(sys.argv[3])
        print("Column used for comparison is", comparecolumn)
    if sys.argv[4]:
        keepcolumn = list(sys.argv[4].split(','))
        print("Columns being kept include:", keepcolumn)
    get_seqs(infile, blast_hits, outfile, comparecolumn, keepcolumn)
    print ("outfile is", outfile)

main()
