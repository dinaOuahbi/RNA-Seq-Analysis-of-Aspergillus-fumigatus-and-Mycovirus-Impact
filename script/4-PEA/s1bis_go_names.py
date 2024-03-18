#!/shared/home/dinaton/mycoenv/bin/python

# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


import pandas as pd
import os
import shutil
import matplotlib.pyplot as plt

counter=0
with open('go-basic.obo', 'r') as file:
    for line in file:
        if line.startswith("[Term]"):
            counter+=1


df = pd.DataFrame(columns=['id', 'name', 'namespace'], index=range(counter-1))


print(f"Open go-basic.obo")
with open('go-basic.obo', 'r') as file:
    current_term = {}
    index = 0
    in_term = False
    for line in file:
        if line.startswith("[Term]"):
            in_term = True
            if current_term:
                df.iloc[index, :] = list(current_term.values())[:3]
                current_term = {}
                index+=1
        elif in_term and line.strip():
            try:
                key, value = line.strip().split(': ', 1)
                current_term[key] = value
            except ValueError:
                print(line)

print('Export results')
df['namespace'].value_counts().plot.pie()
plt.savefig('plots/sets_repartition.png')

try:
    os.mkdir('PEA/')
except FileExistsError:
    pass
try:
    os.mkdir('PEA/Background_genes')
except FileExistsError:
    pass

df.to_csv('PEA/go_term_names.csv')
try:
    shutil.move(f"{os.getcwd()}/FungiDB-65_AfumigatusAf293_Curated_GO.gaf",
                "PEA/Background_genes/FungiDB-65_AfumigatusAf293_Curated_GO.gaf")
except FileNotFoundError:
    pass


print('DONE !')
