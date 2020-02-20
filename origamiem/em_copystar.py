#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-15 11:00:21
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

    parser.add_argument("-i",        "--input",       type=str, help="Particle star file")
    parser.add_argument("-o",        "--output",      type=str, help="Output directory", default=None)
    parser.add_argument("-ref",      "--reference",   type=str, help="Reference class star file")
    parser.add_argument("-cols",     "--columns",     type=str, help="Columns to copy", nargs='+')
    parser.add_argument("-comp",     "--compare",     type=str, help="Criteria to compare", choices=['img','mic'], default='img')

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'output':      args.output,
                 'reference':   args.reference,
                 'columns':     args.columns,
                 'compare':     args.compare
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Check if the reference file exists
    if args_dict['reference'] is None and not os.path.isfile(args_dict['reference']):
        parser.print_help()
        sys.exit('Reference file file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='EMCopy')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    new_project.read_class_refs(args_dict['reference'])
    print('Read class reference file {}'.format(args_dict['reference']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Add new columns
    new_project.copy_from_ref(args_dict['columns'], args_dict['compare'])

    # Write output files
    new_project.write_output_files(write_ref_class_star=False)


if __name__ == "__main__":
    main()
