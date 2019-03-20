#!/usr/bin/env python3

import pandas as pd
import numpy as np


def filter_rsem2(expressdata1, expressdata2=[]):
    #open RSEM file into dataframe
    data1 = pd.read_csv(expressdata1, sep="\t")
    data2 = pd.read_csv(expressdata2, sep="\t")

    #split 'gene_id' column into contig and uniprot id
    new = data1['gene_id'].str.split('__', 5, expand=True)
    new2 = data2['gene_id'].str.split('__', 5, expand=True)

    #input the first entry from the split column into new column named gene_id
    data1['gene_id'] = new[0]
    data2['gene_id'] = new2[0]

    #input the second entry from the split column into new column named uniprot_id
    data1['uniprot_id'] = new[1]
    data2['uniprot_id'] = new2[1]

    #input the 5th entry from the split column into new column named blast match
    data1['blast_match'] = new[5]
    data2['blast_match'] = new2[5]


    #copy uniprot id to new column
    data1['transcript_id'] = data1['uniprot_id']
    data2['transcript_id'] = data2['uniprot_id']

    #sort columns
    cols = ['gene_id','uniprot_id', 'transcript_id', 'length', 'effective_length', 'expected_count', 'TPM', 'FPKM', 'blast_match']

    data1 = data1[cols]
    data2 = data2[cols]

    data1.sort_values(by=['transcript_id'])
    data2.sort_values(by=['transcript_id'])


    #add increments to duplicate uniprot id
    data1['temp_count'] = data1.groupby(['transcript_id']).cumcount()+1
    data1['transcript_id'] = data1.transcript_id.map(str) + "." + data1.temp_count.map(str)
    data1.drop('temp_count', axis=1, inplace=True)

    data2['temp_count'] = data2.groupby(['transcript_id']).cumcount()+1
    data2['transcript_id'] = data2.transcript_id.map(str) + "." + data2.temp_count.map(str)
    data2.drop('temp_count', axis=1, inplace=True)



    # #remove uniprot_ids not found in opposite dataframe
    # common = data1.merge(data2,on=['uniprot_id'])
    # unique_in_1 = data1[(~data1.uniprot_id.isin(common.uniprot_id))]
    #
    # common2 = data2.merge(data1,on=['uniprot_id'])
    # unique_in_2 = data2[(~data2.uniprot_id.isin(common.uniprot_id))]
    #
    # #remove unique uniprot ids from dataframes
    # remove_uniques1 = pd.concat([data1, unique_in_1])
    # data1 = remove_uniques1.drop_duplicates(keep=False)
    #
    # remove_uniques2 = pd.concat([data2, unique_in_2])
    # data2 = remove_uniques2.drop_duplicates(keep=False)
    #
    # unique_in_1.to_csv('~/Documents/unique_in_1.csv', sep='\t', index=False)
    # unique_in_2.to_csv('~/Documents/unique_in_2.csv', sep='\t', index=False)


    data1.to_csv('data1.genes.results', sep='\t', index=False)
    data2.to_csv('data2.genes.results', sep='\t', index=False)

def main():
    filter_rsem2(expressdata1, expressdata2)
    print ("Output is data1.genes.results/data2.genes.results")



#grab arguments from console and pass them to python script
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    ####flags for the required inputs
    req_grp = parser.add_argument_group(title='required arguments')
    req_grp.add_argument("-1", "--expression1", dest = "expression1", required=True, help="Expression data set 1")
    req_grp.add_argument("-2", "--expression2", dest = "expression2", required=True, help="Expression data set 2")
    ####end of flags for the required inputs

    args = parser.parse_args()

    #assign variables from command line arguments
    expressdata1 = args.expression1
    expressdata2 = args.expression2

    print(f"Dataset 1: {expressdata1}")
    print(f"Dataset 2: {expressdata2}")

    #run main program
    main()
