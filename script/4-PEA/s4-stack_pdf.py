#!//shared/home/dinaton/mycoenv/bin/python

# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


import PyPDF2
import os
def merge_pdfs(paths, output):
    merger = PyPDF2.PdfMerger()
    for path in paths:
        merger.append(f'{path}')
    merger.write(output)
    merger.close()
    
    
    
root = 'PEA/imgs/'
barplots = [f'{root}barplot/{i}' for i in os.listdir(f"{root}barplot")]
enrichment = [f'{root}enrichment/{i}' for i in os.listdir(f"{root}enrichment")]
GseaTable = [f'{root}GseaTable/{i}' for i in os.listdir(f"{root}GseaTable")]
# Merge the PDFs
output = "stack_plots/"
try:
    os.mkdir(output)
except:
    pass

merge_pdfs(barplots, f"{output}merged_barplots.pdf")
merge_pdfs(enrichment, f"{output}merged_enrichment.pdf")
merge_pdfs(GseaTable, f"{output}merged_GseaTable.pdf")
print("PDFs merged successfully.")
