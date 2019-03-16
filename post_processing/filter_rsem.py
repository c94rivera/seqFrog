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
data1 = pd.read_csv('/home/christopher/Documents/1RSEM.genes.results', sep="\t")
data2 = pd.read_csv('/home/christopher/Documents/2RSEM.genes.results', sep="\t")

#split 'gene_id' column into contig and uniprot id
new = data1['gene_id'].str.split('__', 5, expand=True)
new2 = data2['gene_id'].str.split('__', 5, expand=True)

#input the first entry from the split column into new column named gene_id
data1['gene_id'] = new[0]
data2['gene_id'] = new2[0]

#input the second entry from the split column into new column named uniprot_id
data1['uniprot_id'] = new[1]
data2['uniprot_id'] = new2[1]

#input the fifth entry from the split column into new column named e-value
data1['e-value'] = new[4]
data2['e-value'] = new2[4]

#input the 6th entry from the split column into new column named blast_match
data1['blast_match'] = new[5]
data2['blast_match'] = new2[5]

#remove data with 0 effective effective_length
data1 = data1[data1.effective_length != 0]
data2 = data2[data2.effective_length != 0]

#remove uniprot_ids not found in opposite dataframe
common = data1.merge(data2,on=['uniprot_id'])
unique_in_1 = data1[(~data1.uniprot_id.isin(common.uniprot_id))]

common2 = data2.merge(data1,on=['uniprot_id'])
unique_in_2 = data2[(~data2.uniprot_id.isin(common.uniprot_id))]

#remove unique uniprot ids from dataframes
remove_uniques1 = pd.concat([data1, unique_in_1])
data1 = remove_uniques1.drop_duplicates(keep=False)

remove_uniques2 = pd.concat([data2, unique_in_2])
data2 = remove_uniques2.drop_duplicates(keep=False)


#create separate dataframes for each uniprot id and place all the dataframes into a dictionary
diff_uniprot = dict(tuple(data1.groupby('uniprot_id')))
diff_uniprot2 = dict(tuple(data2.groupby('uniprot_id')))

#make array of all unique uniprot ids
uniques = data1.uniprot_id.unique()
uniques2 = data2.uniprot_id.unique()

#create empty dataframe
evalue_df = pd.DataFrame()
evalue_df2 = pd.DataFrame()

#run a loop through all of the uniprot ids and take the entry with the lowest e-value(can be more than one)
#and append to the empty dataframe
for i in uniques:
    data_temp = diff_uniprot[i]
    temp = data_temp[data_temp['e-value'] == min(data_temp['e-value'])]
    evalue_df = evalue_df.append(temp, ignore_index=True)

for i in uniques2:
    data_temp = diff_uniprot2[i]
    temp = data_temp[data_temp['e-value'] == min(data_temp['e-value'])]
    evalue_df2 = evalue_df2.append(temp, ignore_index=True)

#remove duplicates of same uniprot_id in DataFrame
evalue_df = evalue_df.drop_duplicates(['uniprot_id'])
evalue_df2 = evalue_df2.drop_duplicates(['uniprot_id'])

#create empty dataframe
tpm_df = pd.DataFrame()
tpm_df2 = pd.DataFrame()

#run through data1 for highest TPM
for i in uniques:
    data_temp = diff_uniprot[i]
    temp = data_temp[data_temp['TPM'] == max(data_temp['TPM'])]
    tpm_df = tpm_df.append(temp, ignore_index=True)

for i in uniques2:
    data_temp = diff_uniprot2[i]
    temp = data_temp[data_temp['TPM'] == max(data_temp['TPM'])]
    tpm_df2 = tpm_df2.append(temp, ignore_index=True)

# #create dataframe with missing uniprot uniprot_ids
# item_missing_in_2 = np.setdiff1d(uniques, uniques2)
# item_missing_in_1 = np.setdiff1d(uniques2, uniques)
#
# missing2_df = pd.DataFrame(item_missing_in_2, columns=['uniprot_id'])
# missing1_df = pd.DataFrame(item_missing_in_1, columns=['uniprot_id'])


# #add missing uniprot ids to dataframe
# missing1_df = pd.concat([pd.DataFrame([i], columns=['uniprot_id']) for i in item_missing_in_1], ignore_index=True, sort=True)
# evalue_df = evalue_df.append(missing1_df, ignore_index=True)
#
# missing2_df = pd.concat([pd.DataFrame([i], columns=['uniprot_id']) for i in item_missing_in_2], ignore_index=True, sort=True)
# evalue_df2 = evalue_df2.append(missing2_df, ignore_index=True)

#reorganize the columns
# 'gene_id', 'transcript_id(s)', 'length', 'effective_length', 'expected_count', 'TPM', 'FPKM', 'gene_id', 'uniprot_id', 'e-value', 'blast_match'
cols = ['gene_id', 'transcript_id(s)', 'length', 'effective_length', 'expected_count', 'TPM', 'FPKM', 'gene_id', 'uniprot_id', 'e-value', 'blast_match']
evalue_df = evalue_df[cols]
evalue_df2 = evalue_df2[cols]

#remove duplicated columns
evalue_df = evalue_df.loc[:,~evalue_df.columns.duplicated()]
evalue_df2 = evalue_df2.loc[:,~evalue_df2.columns.duplicated()]

#sort Columns
evalue_df = evalue_df.sort_values(by=['uniprot_id'])
evalue_df2 = evalue_df2.sort_values(by=['uniprot_id'])
#find common values between files
# s3 = list(set(uniques) & set(uniques2))
# s3
# #find uniprot_ids not present in both dataframes actual
# merged = evalue_df.merge(evalue_df2, indicator=True, how='outer')
# merged[merged['_merge'] == 'right_only']

#print dataframes to file
evalue_df.to_csv('~/Documents/epithelium1.genes.results', sep='\t', index=False)
evalue_df2.to_csv('~/Documents/epithelium2.genes.results', sep='\t', index=False)

unique_in_1.to_csv('~/Documents/unique_in_1.csv', sep='\t', index=False)
unique_in_2.to_csv('~/Documents/unique_in_2.csv', sep='\t', index=False)








#trash for testing!!!
# diff_uniprot['Q9NAX4']
# list(diff_uniprot.keys())
