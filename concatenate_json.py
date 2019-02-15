#!/usr/bin/env python3

import os
import json
import sys
import shlex, subprocess


project_name='rna_fusco'
filename='{}.json'.format(project_name)

def get_read_pair(s):
    basename=s.split('/')[-1]
    read=basename.split('_')[3]
    if read in ['R1', 'R2']:
        return read
    else:
        sys.exit()


with open(filename, "r") as inputfile:
     data = json.loads(inputfile.read())

print(len(data['samples']))
print(len(data['units']))

reads = {}
for sample, units in data['samples'].items():
    # print(sample)
    reads[sample] = {'R1': [],
                     'R2': []}
    for unit in units:
       # print(unit)
       # print(data['units'][unit])
      
        for f in data['units'][unit]:
            reads[sample][get_read_pair(f)].append(f)
           # print(f.split('/')[-1])
           # print(get_read_pair(f))


new_samples = {}
for s, pairs in reads.items():
    cmd = ['cat']
    for p in pairs['R1']:
        cmd.append(p)
    cmd.append('> datasets/{}_R1.fastq.gz'.format(s))
    print(cmd)
    subprocess.run(' '.join(cmd), shell=True)
    cmd = ['cat']
    for p in pairs['R2']:
        cmd.append(p)
    cmd.append('> datasets/{}_R2.fastq.gz'.format(s))
    print(cmd)
    #subprocess.run(' '.join(cmd), shell=True)
    new_samples[s] = ['/ELS/els9/users/biosciences/projects/rna_fusco/datasets/{}_R1.fastq.gz'.format(s)]
    
#print(new_samples)
#cmdt = ['cat', '/tmp/RNA-Seq/RNA17-0181-R0001_S12_L001_R1_001.fastq.gz', '/tmp/RNA-Seq/RNA17-0181-R0001_S12_L004_R1_001.fastq.gz', '/tmp/RNA-Seq/RNA17-0181-R0001_S12_L005_R1_001.#fastq.gz', '> /tmp/datasets/RNA17-0196-R0001_R1.fastq.gz']
#print(cmdt)
#subprocess.run(' '.join(cmdt), shell=True)
#print(reads['RNA17-0196-R0001'])

json_template = 'config.test.data.json'
with open(json_template, "r") as inputfile:
     new_data = json.loads(inputfile.read())


new_data['samples']=new_samples
print(new_data)

json_project = 'config.project.{}.json'.format(project_name)
with open(json_project, "w") as outfile:
     json.dump(new_data, outfile, indent=4)
