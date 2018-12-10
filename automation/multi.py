#!/usr/bin/env python3
#
#Created on Sat Jun  7 14:27:49 2014
#
#@author: RDT
#adapted by: CAR
#

from Bio import SeqIO
import csv
import sys
import os
from tqdm import tqdm
from multiprocessing import Pool

def compare(rec):
    global matchlist, proteinid, seqlist, comparecolumn, keepcolumn

    proteinid.clear()
    for match in matchlist:
        if match[0] == rec.id:
            for x in comparecolumn:
                if match[x] not in proteinid:
                    proteinid.append(match[x])
                    rec2 = rec[:]

                    for column in keepcolumn:
                        rec2.description += "\t" + match[int(column)]

                    # seqlist.append(rec2)
                    return(rec2)
                continue
            else:
                continue
        else:
            continue



def get_seqs_fast(seqfile, blast_file, outfile, comparecolumn, keepcolumn):
    global proteinid, matchlist
    records = SeqIO.parse(seqfile, "fasta")
    matchlist=[]
    proteinid = []
    #return list with blast hits
    with open(blast_file,"rU") as f:
        reader = csv.reader(f,delimiter="\t")
        for row in reader:
            matchlist.append(row)
        print("Blast file successfully read")
    print("\nComparing contigs to Blast file...")

    #start multiprocessing based on number of cores and save output to results variable
    p = Pool(os.cpu_count())
    results = p.map(compare, records)

        #remove empty items in results variable
    tempResults = results
    tempResults[:] = [item for item in tempResults if item]

    #write final list of matches to file in fasta format
    print("Writing annotated contigs")
    with open(outfile, "w") as f:
        SeqIO.write(tempResults, f, "fasta")



#grab arguments from console and pass them to python script
if __name__ == '__main__':
    global records, matchlist, seqlist, comparecolumn, keepcolumn
    import argparse
    parser = argparse.ArgumentParser()

    ####flags for the required inputs
    req_grp = parser.add_argument_group(title='required arguments')
    req_grp.add_argument("-c", "--contig_file", dest = "contigs", required=True, help="Contigs input file ")
    req_grp.add_argument("-b", "--blast_hits", dest = "blast_hits", required=True, help="Blast Input File")
    ####end of flags for the required inputs

    parser.add_argument("-comp", "--colcompare", dest = "comparecolumn", help = "Column for Comparison, starts counting at '0' (default: 1)", default = "1")
    parser.add_argument("-keep", "--colkeep", dest = "keepcolumn", help = "Columns to Keep, starts counting at '0' (default: 1, 2, 3, 10,13)", default = "1,2,3,10,13")

    args = parser.parse_args()

    #assign variables from command line arguments
    contigs = args.contigs
    blast_hits = args.blast_hits
    comparecolumn = [int(x) for x in args.comparecolumn.split(",")]
    keepcolumn = [int(x) for x in args.keepcolumn.split(",")]

    print("contigs is", contigs)
    print("blast file is", blast_hits)
    outfile = contigs + "_matches.fasta"#blast_hits[:blast_hits.index(".")]+".fasta"
    print("Column used for comparison is", comparecolumn)
    print("Columns being kept include:", keepcolumn)

    #Open fasta contigs
    # records = SeqIO.parse(contigs, "fasta")
    # matchlist=[]
    # proteinid = []
    # #return list with blast hits
    # with open(blast_hits,"rU") as f:
    #     reader = csv.reader(f,delimiter="\t")
    #     for row in reader:
    #         matchlist.append(row)
    #     print("Blast file successfully read")
    # print("\nComparing contigs to Blast file...")
    #
    # #start multiprocessing based on number of cores and save output to results variable
    # p = Pool(os.cpu_count())
    # results = p.map(compare, records)
    #
    #     #remove empty items in results variable
    # tempResults = results
    # tempResults[:] = [item for item in tempResults if item]
    #
    # #write final list of matches to file in fasta format
    # print("Writing annotated contigs")
    # with open(outfile, "w") as f:
    #     SeqIO.write(tempResults, f, "fasta")
