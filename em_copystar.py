#!/usr/bin/env python
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
    parser.add_argument("-cols",     "--columns",     type=str, help="Columns to copy", default=None)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'output':      args.output,
                 'reference':   args.reference,
                 'columns':     args.columns
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Check if the reference file exists
    if args_dict['reference'] is None and not os.path.isfile(args_dict['reference']):
        parser.print_help()
        sys.exit('Reference file file does not exist!')

    # Get the new column parameters
    if args_dict['columns'] is not None:
        new_column_parameters = util.parse_star_parameters(args_dict['columns'])
    else:
        new_column_parameters = None

    # Create an EM project object
    new_project = em.Project(name='ProjectCopy')
    new_project.set_output_directory(args_dict['input'], args_dict['output'])

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
    new_project.copy_from_ref(new_column_parameters)

    # Write output files
    new_project.write_output_files(write_ref_class_star=False)


if __name__ == "__main__":
    main()
