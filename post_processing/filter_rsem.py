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
data

#make array of all unique uniprot ids
uniques = data.uniprot_id.unique()
uniques
uniques2 = data2.uniprot_id.unique()
uniques2

#find uniprot_ids not present in both dataframes actual
s1 = np.setdiff1d(uniques, uniques2)
s2 = np.setdiff1d(uniques2, uniques)

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

#remove duplicates of same uniprot_id in DataFrame
evalue_df = evalue_df.drop_duplicates(['uniprot_id'])
evalue_df

evalue_df2 = evalue_df2.drop_duplicates(['uniprot_id'])
evalue_df2

#add missing values from between the two dataframes
# for i in s2:
#     diff_uniprot[i] = diff_uniprot.get(i, ['', '', '', '', '', '', '', i, '', ''])
#
# for i in s1:
#     diff_uniprot2[i] = diff_uniprot2.get(i, ['', '', '', '', '', '', '', i, '', ''])

#create dataframe with missing uniprot uniprot_ids
missing_in_2 = np.setdiff1d(uniques, uniques2)
missing_in_1 = np.setdiff1d(uniques2, uniques)

missing2 = pd.DataFrame(missing_in_2, columns=['uniprot_id'])
missing1 = pd.DataFrame(missing_in_1, columns=['uniprot_id'])

#add missing uniprot ids to dataframe
evalue_df = evalue_df.append((missing1), ignore_index=True)
evalue_df = evalue_df.append((missing2), ignore_index=True) #maybe dont ignore index?

missing2 = pd.concat([pd.DataFrame([i], columns=['uniprot_id']) for i in missing_in_2], ignore_index=True)
evalue_df = evalue_df.append(missing2, ignore_index=True)

#reorganize the columns
cols = ['gene_id', 'transcript_id(s)', 'length', 'effective_length', 'expected_count', 'TPM', 'FPKM', 'gene_id', 'uniprot_id', 'e-value', 'blast_match']
evalue_df = evalue_df[cols]
evalue_df
# 'gene_id', 'transcript_id(s)', 'length', 'effective_length', 'expected_count', 'TPM', 'FPKM', 'gene_id', 'uniprot_id', 'e-value', 'blast_match'


#find common values between files
s3 = list(set(uniques) & set(uniques2))
s3
#find uniprot_ids not present in both dataframes actual
merged = pd.merge(evalue_df, evalue_df2, indicator=True, how='outer')
merged[merged['_merge'] == 'right_only']

#print dataframes to file
evalue_df.to_csv('~/Documents/evalue_df.csv', sep='\t', index=False)







#trash for testing!!!
# diff_uniprot['Q9NAX4']
# list(diff_uniprot.keys())
