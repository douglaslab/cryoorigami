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

    parser.add_argument("-i",           "--input",       type=str, help="Particle star file")
    parser.add_argument("-o",           "--output",      type=str, help="Output directory", default=None)
    parser.add_argument("-cols",        "--columns",     type=str, help="Columns to copy", nargs='*', default=None)
    parser.add_argument("-copy-priors", "--copypriors",  action='store_true', help="Copy offset and angle parameters to priors")

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'output':      args.output,
                 'columns':     args.columns,
                 'copypriors':  args.copypriors
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Get the new column parameters
    if args_dict['columns'] is not None:
        new_column_parameters = util.parse_star_parameters(args_dict['columns'])
    else:
        new_column_parameters = None

    # Copy priors overwrites
    if args_dict['copypriors']:
        new_column_parameters = {'rlnOriginX':'rlnOriginXPrior',
                                 'rlnOriginY':'rlnOriginYPrior',
                                 'rlnAnglePsi':'rlnAnglePsiPrior'}

    # Create an EM project object
    new_project = em.Project(name='ProjectCopyColumns')
    new_project.set_output_directory(args_dict['input'], args_dict['output'])

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Add new columns
    new_project.copy_columns(new_column_parameters)

    # Write output files
    new_project.write_output_files(write_ref_class_star=False)


if __name__ == "__main__":
    main()
