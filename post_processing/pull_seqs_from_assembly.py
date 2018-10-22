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
from tqdm import tqdm

def get_seqs(seqfile, matches, outfile, comparecolumn, keepcolumn):
    #Open fasta contigs and return iterator of SeqRecords with protein sequences.
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


    print("Writing annotated contigss")
    with open(outfile, "w") as f:
        SeqIO.write(seqlist, f, "fasta")



def main():
    get_seqs(contigs, blast_hits, outfile, comparecolumn, keepcolumn)
    print ("outfile is", outfile)



#grab arguments from console and pass them to python script
import argparse
parser = argparse.ArgumentParser()

#argument tags
# parser.add_argument("contigs", help = "Contigs Input file")
# parser.add_argument("blast_hits", help = "Blast Input")

####flags for the required inputs
req_grp = parser.add_argument_group(title='required arguments')
req_grp.add_argument("-c", "--contig_file", dest = "contigs", required=True, help="Contigs input file ")
req_grp.add_argument("-b", "--blast_hits", dest = "blast_hits", required=True, help="Blast Input File")
####end of flags for the required inputs

parser.add_argument("-comp", "--colcompare", dest = "comparecolumn", help = "Column for Comparison, starts counting at '0' (default: 1)", default = "1")
parser.add_argument("-keep", "--colkeep", dest = "keepcolumn", help = "Columns to Keep, starts counting at '0' (default: 1, 2, 3, 4)", default = "1,2,3,4")

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

#run main program
main()
