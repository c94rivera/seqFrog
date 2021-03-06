#!/usr/bin/env python3

import pandas as pd
import numpy as np

#open RSEM file into dataframe
data1 = pd.read_csv('/home/litoria/Desktop/media/Christopher_thesis/1/quants/opumilio1/quant.sf', sep="\t")
data2 = pd.read_csv('/home/litoria/Desktop/media/Christopher_thesis/2/quants/opumilio2/quant.sf', sep="\t")

#split 'Name' column into contig and uniprot id
new = data1['Name'].str.split('__', 5, expand=True)
new2 = data2['Name'].str.split('__', 5, expand=True)

#input the first entry from the split column into new column named Name
data1['Name'] = new[0]
data2['Name'] = new2[0]

#input the second entry from the split column into new column named uniprot_id
data1['uniprot_id'] = new[1]
data2['uniprot_id'] = new2[1]

#input the 5th entry from the split column into new column named blast match
data1['blast_match'] = new[5]
data2['blast_match'] = new2[5]
data1

#copy uniprot id to new column
data1['transcript_id'] = data1['uniprot_id']
data2['transcript_id'] = data2['uniprot_id']

#sort columns
cols = ['Name','uniprot_id', 'transcript_id', 'Length', 'EffectiveLength', 'TPM', 'NumReads', 'blast_match']

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


data1.to_csv('~/Documents/epithelium1.quant.sf', sep='\t', index=False)
data2.to_csv('~/Documents/epithelium2.quant.sf', sep='\t', index=False)
