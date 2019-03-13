#!/usr/bin/env python3
# grab lines based on e-value, length, highest expression
# format the first column so gene id is alone and searchable (maybe put this in 2nd column and replace)
#     search for uniprot gene id
#     group all lines with same gene id and perform function to that group
# print entire line to new file

import pandas as pd

#open RSEM file into dataframe
data = pd.read_csv('/home/christopher/Documents/1RSEM.genes.results', sep="\t")

#split 'gene_id' column into contig and uniprot id
new = data['gene_id'].str.split('__', 5, expand=True)
#input the first entry from the split column into new column named gene_id
data['gene_id'] = new[0]
#input the second entry from the split column into new column named uniprot_id
data['uniprot_id'] = new[1]
#input the fifth entry from the split column into new column named e-value
data['e-value'] = new[4]
#input the 6th entry from the split column into new column named blast_match
data['blast_match'] = new[5]

#make array of all unique uniprot ids
uniques = data.uniprot_id.unique()
uniques

#create separate dataframes for each uniprot id and place all the dataframes into a dictionary
diff_uniprot = dict(tuple(data.groupby('uniprot_id')))

###not working
new_df = pd.DataFrame()
for i in uniques:
    df2 = diff_uniprot[i]
    temp = df2[df2['e-value'] == min(df2['e-value'])]
    new_df.append(temp, ignore_index=True)
new_df

df2 = diff_uniprot[uniques[1]]
df2
###


new_df = df2[df2['e-value'] == min(df2['e-value'])]
new_df
