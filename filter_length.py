#!/usr/bin/env python3

def filter_length(input_fasta, trim_length):
    sequences = {}
    good_reads = []
    bad_reads = []
    with open(input_fasta, "r") as file:
        for line in file:
            line=line.rstrip()
            if (line[0] == ">"):
                header = line
                sequences[header] = ""
            else:
                data = line
                sequences[header] += data


    # figure out which reads are good/bad
    for header in sequences.keys():
        if (len(sequences[header]) > int(trim_length)):
            good_reads.append(header)
        else:
            bad_reads.append(header)

    # write good reads
    with open("good.fasta", "w+") as good_out:
        for header in good_reads:
            good_out.write(f"{header}\n{sequences[header]}\n")

    # write bad reads
    with open("bad.fasta", "w+") as bad_out:
        for header in bad_reads:
            bad_out.write(f"{header} Excluded because too short\n{sequences[header]}")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    ####flags for the required inputs
    req_grp = parser.add_argument_group(title='required arguments')
    req_grp.add_argument("-f", "--fasta", dest = "input_fasta", required=True, help="Input FASTA file")
    req_grp.add_argument("-l", "--trim_length", dest = "trim_length", required=True, help="Minimum contig length for trim")

    args = parser.parse_args()

    #assign variables from command line arguments
    input_fasta = args.input_fasta
    trim_length = args.trim_length
    filter_length(input_fasta, trim_length)
