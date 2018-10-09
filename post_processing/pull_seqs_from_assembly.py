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
    #Open fasta infile and return iterator of SeqRecords with protein sequences.
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
                for x in comparecolumn:
                    if match[x] not in proteinid:
                        proteinid.append(match[x])
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
    print("infile is", infile)
    print("blast file is", blast_hits)
    outfile = infile + "_matches.fasta"#blast_hits[:blast_hits.index(".")]+".fasta"
    if comparecolumn:
        print("Column used for comparison is", comparecolumn)
    if keepcolumn:
        print("Columns being kept include:", keepcolumn)
    get_seqs(infile, blast_hits, outfile, comparecolumn, keepcolumn)
    print ("outfile is", outfile)

#grab arguments from console and pass them to python script
import argparse
parser = argparse.ArgumentParser()

#argument tags
parser.add_argument("-c", "--contig", dest = "infile", help = "Contig Input")
parser.add_argument("-b", "--blast", dest = "blast_hits", help = "Blast Input")
parser.add_argument("-comp", "--colcompare", dest = "comparecolumn", help = "Column for Comparison", default = "1")
parser.add_argument("-keep", "--colkeep", dest = "keepcolumn", help = "Columns to Keep", default = "1,2,3,4")

args = parser.parse_args()

#assign variables from command line arguments
infile = args.infile
blast_hits = args.blast_hits
comparecolumn = [int(x) for x in args.comparecolumn.split(",")]
keepcolumn = [int(x) for x in args.keepcolumn.split(",")]

#run main program
main()
