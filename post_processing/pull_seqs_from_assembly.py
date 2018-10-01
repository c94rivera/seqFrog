# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 14:27:49 2014

@author: RDT
adapted by: CAR
"""

from Bio import SeqIO
import csv
import sys

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
    for rec in records:
        for item in matchlist:
            if item[0] == rec.id:
                #######
                #find unique terms in column 1
                rec.description+="\t"
                rec.description+=item[1] + "\t" + item[2] + "\t" + item[3] + "\t" + item[10] + "\t" + item[13] + "\t"
                seqlist.append(rec)
            else:
                continue

    print("Writing annotated contigs")
    with open(outfile, "w") as f:
        SeqIO.write(seqlist, f, "fasta")

def delete_dupes(outfile):
    output_dict = SeqIO.to_dict(SeqIO.parse(outfile, "fasta"))
    final_list = []
    print("Removing Duplicates")
    for outputrec in output_dict:
        if outputrec not in final_list:
            final_list.append(outputrec)

    print("Writing output file")
    with open(outfile, "w") as o:
        SeqIO.write(final_list, o, "fasta")

def main():
    assert len(sys.argv) == 3, "usage: python3 pullseqs_fromtrinity.py <sequences.fasta> <blast_hits.txt>"
    infile = sys.argv[1]
    print("infile is", infile)
    blast_hits = sys.argv[2]
    print("blast file is", blast_hits)
    outfile = blast_hits[:blast_hits.index(".")]+".fasta"
    get_seqs(infile, blast_hits, outfile)
    delete_dupes(outfile)
    print ("outfile is", outfile)

main()
