# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


# Read Text Files with Pandas using read_csv()
 
# importing pandas
import pandas as pd
import os
os.chdir('/shared/projects/braf_mutated_melanoma/mycovirus')
# read text file into pandas DataFrame
names = pd.read_table("Af293_geneName.tsv")
# display DataFrame
names

try:
    os.mkdir('deseq_v2')
except:
    pass

def get_desc(root, data_name, out):
    print('-'*50, data_name)
    df = pd.read_csv(f'{root}/{data_name}.csv', sep=';')
    df.rename(columns={'Unnamed: 0':'Symbol'}, inplace=True)
    df.Symbol = [i.split(':')[1].strip() for i in df.Symbol]

    print(f"NO data genes : {names.shape} ----> NO find genes : {df.shape}")
    final = pd.merge(names, df, how='inner',on='Symbol')
    final.to_csv(f'{out}/{data_name}.csv')
    
    
# main for deseq comparisons
root = "deseq2"
out = "deseq_v2"
datas = [i.split('.')[0].strip() for i in os.listdir(root) if i.endswith('csv')]
for data_name in datas:
    get_desc(root,data_name, out)
    


print('DONE !')