import re
import os

with open("sdf_paper.tex", "r") as f:
    paper = f.read()

with open("empirical/empirical_PE_MKT_tables.tex", "r") as f:
    mkt = f.read()
    
with open("empirical/empirical_PE_all_factors_tables.tex", "r") as f:
    all_factors = f.read()
    
# Extract individual tables from all_factors
tables = re.findall(r'\\begin\{table\}.*?\\end\{table\}', all_factors, re.DOTALL)

# Also extract the mkt one
mkt_table = re.findall(r'\\begin\{table\}.*?\\end\{table\}', mkt, re.DOTALL)[0]

# Now let's replace them in the paper string
paper = re.sub(r'\\begin\{table\}.*?\\label\{tab:empirical_pe_mkt\}.*?\\end\{table\}', mkt_table, paper, flags=re.DOTALL)

for t in tables:
    label = re.search(r'\\label\{tab:empirical_pe_(\w+)\}', t).group(1)
    paper = re.sub(r'\\begin\{table\}.*?\\label\{tab:empirical_pe_' + label + r'\}.*?\\end\{table\}', t, paper, flags=re.DOTALL)

with open("sdf_paper.tex", "w") as f:
    f.write(paper)
    
print("Updated sdf_paper.tex")
