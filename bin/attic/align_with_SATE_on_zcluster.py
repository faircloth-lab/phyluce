#!/usr/bin/env python
# encoding: utf-8
"""
File: run_raxml_genetrees_on_zcluster.py
Author: Brant Faircloth

Created by Brant Faircloth on 07 March 2013 12:03 PST (-0800)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import glob
#import time
import argparse
import subprocess


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class CreateDir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # get the full path
        d = os.path.abspath(os.path.expanduser(values))
        # check to see if directory exists
        if os.path.exists(d):
            answer = raw_input("[WARNING] Output directory exists, REMOVE [Y/n]? ")
            if answer == "Y":
                shutil.rmtree(d)
            else:
                print "[QUIT]"
                sys.exit()
        # create the new directory
        os.makedirs(d)
        # return the full path
        setattr(namespace, self.dest, d)


def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_file_extensions(ftype):
    ext = {
        'fasta': ('.fasta', '.fsa', '.aln', '.fa'),
        'nexus': ('.nexus', '.nex'),
        'phylip': ('.phylip', '.phy'),
        'phylip-relaxed': ('.phylip', '.phy'),
        'clustal': ('.clustal', '.clw'),
        'emboss': ('.emboss',),
        'stockholm': ('.stockholm',)
    }
    return ext[ftype]


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Run SATE on zcluster (qsub) system""")
    parser.add_argument(
            "--input",
            required=True,
            action=FullPaths,
            type=is_dir,
            help="""The path holding the alignments"""
        )
    parser.add_argument(
            "--output",
            required=True,
            action=CreateDir,
            default="sate-output",
            help="""The directory to hold output files"""
        )
    parser.add_argument(
            "--control",
            default="submit-scripts",
            help="""The directory to hold control files"""
        )
    parser.add_argument(
            "--config",
            action=FullPaths,
            type=is_file,
            help="""The path to the SATE config file"""
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta'],
            default='fasta',
            help="""The input alignment format""",
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


def main():
    args = get_args()
    controldir = os.path.join(args.output, "submit-scripts")
    os.mkdir(controldir)
    # change to the controldir
    os.chdir(controldir)
    # iterate over files, get file name
    for f in get_files(args.input, args.input_format):
        name = os.path.basename(f).split(".")[0]
        output_dir = os.path.join(args.output, name)
        script = name + "_control.sh"
        cli = "#!/bin/bash\ncd {0}\ntime python2.7 /usr/local/sate/latest/run_sate.py -i {1} -j {2} -o {3} {4}\n".format(
            controldir,
            f,
            name,
            output_dir,
            args.config
            )
        control = open(script, 'w')
        control.write(cli)
        control.close()
        # run SATe
        qsub = [
                'qsub',
                '-q',
                "rcc-30d",
                script
            ]
        stderr, stdout = subprocess.Popen(
                qsub,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE
            ).communicate()
        print stdout, stderr


if __name__ == '__main__':
    main()
