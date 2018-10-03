# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 14:27:49 2014

@author: RDT
adapted by: CAR
"""

from Bio import SeqIO
import csv
import sys

def get_seqs(contig_file_path, blast_file_path, outfile):
    '''Open fasta infile and return iterator of SeqRecords with protein sequences.'''
    records = SeqIO.parse(contig_file_path, 'fasta')
    seqlist=[]
    matchlist=[]
    with open(blast_file_path,'rU') as f:
        reader = csv.reader(f,delimiter='\t')
        for row in reader:
            matchlist.append(row)
        print("Blast file successfully read")
    print("Comparing contigs to Blast file")
    for rec in records:
        for item in matchlist:
            if item[0] == rec.id:
                #######
                #find unique terms in column 1
                rec.description+=' '
                rec.description+=item[1] + " " + item[2] + " " + item[3] + " " + item[10] + " " + item[13] + "\t"
                seqlist.append(rec)
            else:
                continue

    print("Writing output file")
    with open(outfile, 'w') as f:
        SeqIO.write(seqlist, f, 'fasta')


def main():
    assert len(sys.argv) == 3, "usage: python3 pullseqs_fromtrinity.py <sequences.fasta> <blast_hits.txt>"
    infile = sys.argv[1]
    print("infile is", infile)
    blast_hits = sys.argv[2]
    print("blast file is", blast_hits)
    outfile = blast_hits[:blast_hits.index('.')]+'.fasta'
    get_seqs(infile, blast_hits, outfile)
    print ("outfile is", outfile)


#User input

contig_file_path = input("Please drag contig file to the terminal and press enter: ")
print('\n', "The file you selected is ", contig_file_path)

blast_file_path = input("Please drag blast file to the terminal and press enter: ")
print('\n', "The file you selected is ", blast_file_path)


main()
