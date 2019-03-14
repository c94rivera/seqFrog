#!/usr/bin/env python3
# grab lines based on e-value, length, highest expression
# format the first column so gene id is alone and searchable (maybe put this in 2nd column and replace)
#     search for uniprot gene id
#     group all lines with same gene id and perform function to that group
# print entire line to new file
# compare both files for genes and if a gene doesn't exist in one file print it with 0's so it can still be compared
#   find differences between the dataframes, concatenate the differences into opposite dataframe?
#     use column with the blast matches?
#     use uniprot id?
# use effective length if you have duplicates in e-value column

import pandas as pd
import numpy as np

#open RSEM file into dataframe
data = pd.read_csv('/home/christopher/Documents/RSEM.genes.results', sep="\t")
data2 = pd.read_csv('/home/christopher/Documents/2RSEM.genes.results', sep="\t")

#split 'gene_id' column into contig and uniprot id
new = data['gene_id'].str.split('__', 5, expand=True)
new2 = data2['gene_id'].str.split('__', 5, expand=True)

#input the first entry from the split column into new column named gene_id
data['gene_id'] = new[0]
data2['gene_id'] = new2[0]

#input the second entry from the split column into new column named uniprot_id
data['uniprot_id'] = new[1]
data2['uniprot_id'] = new2[1]

#input the fifth entry from the split column into new column named e-value
data['e-value'] = new[4]
data2['e-value'] = new2[4]

#input the 6th entry from the split column into new column named blast_match
data['blast_match'] = new[5]
data2['blast_match'] = new2[5]

data2
#make array of all unique uniprot ids
uniques = data.uniprot_id.unique()
uniques
uniques2 = data2.uniprot_id.unique()
uniques2

#create separate dataframes for each uniprot id and place all the dataframes into a dictionary
diff_uniprot = dict(tuple(data.groupby('uniprot_id')))
diff_uniprot2 = dict(tuple(data2.groupby('uniprot_id')))

#create empty dataframe
evalue_df = pd.DataFrame()
evalue_df2 = pd.DataFrame()

#run a loop through all of the uniprot ids and take the entry with the lowest e-value(can be more than one)
#and append to the empty dataframe
for i in uniques:
    data_temp = diff_uniprot[i]
    temp = data_temp[data_temp['e-value'] == min(data_temp['e-value'])]
    evalue_df = evalue_df.append(temp, ignore_index=True)
evalue_df

for i in uniques2:
    data_temp = diff_uniprot2[i]
    temp = data_temp[data_temp['e-value'] == min(data_temp['e-value'])]
    evalue_df2 = evalue_df2.append(temp, ignore_index=True)
evalue_df2

#create empty dataframe
tpm_df = pd.DataFrame()
tpm_df2 = pd.DataFrame()

#run through data for highest TPM
for i in uniques:
    data_temp = diff_uniprot[i]
    temp = data_temp[data_temp['TPM'] == max(data_temp['TPM'])]
    tpm_df = tpm_df.append(temp, ignore_index=True)
tpm_df

for i in uniques2:
    data_temp = diff_uniprot2[i]
    temp = data_temp[data_temp['TPM'] == max(data_temp['TPM'])]
    tpm_df2 = tpm_df2.append(temp, ignore_index=True)
tpm_df2

#find uniprot ids are not present in both dataframes and create empty entries
##boolean
np.in1d(uniques, uniques2)
np.in1d(uniques2, uniques)

#actual identifiers for ones missing
s1 = np.setdiff1d(uniques, uniques2)
s2 = np.setdiff1d(uniques2, uniques)
s1

#iterate through list of missing values and add to the opposite dictionary
for i in s1:
    # {'key': 'value'}
    diff_uniprot2[i] = None
diff_uniprot2

for i in s2:
    diff_uniprot[i] = None
diff_uniprot


#find common values between files
s3 = list(set(uniques) & set(uniques2))
s3
