import re
import sys
from glob import glob

#first, parse the smk files
pattern = re.compile(r'rule (\w+):')
subgraphs = dict()

for smk in glob('../workflow/rules/*.smk'):

    smk_name = re.findall(re.compile(r'/(\w+).smk'),smk)[0]

    with open(smk,'r') as f:
        text = f.read()

    subgraphs[smk_name] = re.findall(pattern,text)


#now, need to parse the dag to get the label numbers for each rule
text_dag = sys.stdin.read()

clusnum=0
        
for smk in subgraphs.keys():

    labels = list()
    for rule in subgraphs[smk]:
        pattern_string = '(\d+)\[label = "{rule}'.format(rule=rule)
        pattern_dag = re.compile(pattern_string) 
        match = re.search(pattern_dag,text_dag)
        if match:
            labels.append(match.group(1))
       #     print(f'pattern: {pattern_string}, found match: {match.group(1)}')
       # else:
       #     print(f'NOT FOUND: pattern: {pattern_string}')


    if len(labels) >0:    
        print(f'    subgraph cluster_{clusnum}')
        print('    {')
#        print('    node [style=filled];')
        print('    '+' '.join(labels))
        print(f'    label = "{smk}";')
        print(f'    color=blue;')
        print('    }')
        clusnum = clusnum+1



