#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-19 10:22:08
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import sys
import argparse
import origamiem as em
import utilities as util


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i",       "--input",       type=str, help="Particle star file")
    parser.add_argument("-pass",    "--passthrough", type=str, help="Passthrough file from cryosparc", default=None)
    parser.add_argument("-o",       "--output",      type=str, help="Output directory", default=None)
    parser.add_argument("-orig",    "--original",    type=str, help="Original star file", default=None)
    parser.add_argument("-micpath", "--micpath",     type=str, help="Micrographs path", default="Micrographs")

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':        args.input,
                 'passthrough':  args.passthrough,
                 'output':       args.output,
                 'original':     args.original,
                 'micpath':      args.micpath
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='ProjectCSparc')
    new_project.set_output_directory(args_dict['input'], args_dict['output'])

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Set cs files
    new_project.set_cs_files(blob_cs_file=args_dict['input'],
                             passthrough_cs_file=args_dict['passthrough'],
                             original_star_file=args_dict['original'])

    # Read cs files
    new_project.read_cs_files()

    # Prepare input and output files
    new_project.prepare_io_files_cs()

    # Conert cs to star
    new_project.convert_cs2star()

    # Write output files
    new_project.write_output_files(args_dict['micpath'])


if __name__ == "__main__":
    main()
