#!/usr/bin/env python3

import os
import yaml
import sys
import shlex, subprocess
import argparse


class App(object):
    def __init__(self, args=None):
        self.input_file = args.input_file
        if os.path.exists(args.folder):
            raise ValueError("folder must not exist")
        else: os.mkdir(args.folder)
        self.folder = args.folder
        self.project = args.project_name if args.project_name else os.getcwd().split('/')[-1]

    def get_read_pair(self,s):
        basename = s.split('/')[-1]
        read = basename.split('_')[3]
        if read in ['R1', 'R2']:
            return read
        else:
            sys.exit()


    def run(self):
        with open(self.input_file, "r") as inputfile:
            data = yaml.load(inputfile.read())
        reads = {}
        for sample, units in data['samples'].items():
            reads[sample] = {'R1': [],
                             'R2': []}
            for unit in units:
                for f in data['units'][unit]:
                    reads[sample][self.get_read_pair(f)].append(f)
        new_samples = {}
        for s, pairs in reads.items():
            cmd = ['cat']
            for p in pairs['R1']:
                 cmd.append(p)
            cmd.append('>' + self.folder + '/{}_R1.fastq.gz'.format(s))
            print(cmd)
            subprocess.run(' '.join(cmd), shell=True)
            cmd = ['cat']
            for p in pairs['R2']:
                cmd.append(p)
            cmd.append('>' + self.folder + '/{}_R2.fastq.gz'.format(s))
            print(cmd)
        # subprocess.run(' '.join(cmd), shell=True)
            workdir = os.getcwd()
            new_samples[s] = workdir + '/' + self.folder + '/{}_R1.fastq.gz'.format(s)
        yaml_template = 'config.template.yaml'
        with open(yaml_template, "r") as inputfile:
            new_data = yaml.load(inputfile.read())
        new_data['samples'] = new_samples
        print(new_data)
        yaml_project = 'config.project.{}.yaml'.format(self.project)
        with open(yaml_project, "w") as outfile:
            yaml.dump(new_data, outfile, indent=4)


def make_parser():
    parser = argparse.ArgumentParser(description='Prepare file for pipeline')

    parser.add_argument('--input_file', '-i', type=str, required=True,
                        help='yaml input file')

    parser.add_argument('--folder', '-f', metavar="PATH", required=True,
                        help="destination folder for merged fastq files")

    parser.add_argument('--project_name', '-p', metavar="PATH",
                        help="Project name for config rename")

    return parser


def main(argv):
    parser = make_parser()
    args = parser.parse_args(argv)

    workflow = App(args=args)

    workflow.run()

#################################################à

if __name__ == '__main__':
    main(sys.argv[1:])