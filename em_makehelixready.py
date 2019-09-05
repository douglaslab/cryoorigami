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

    parser.add_argument("-i",             "--input",        type=str, help="Particle star file")
    parser.add_argument("-o",             "--output",       type=str, help="Output directory", default=None)
    parser.add_argument("-copytopriors",  "--copypriors",   action='store_true', help="Copy offset and/or angle parameters to priors", default=True)
    parser.add_argument("-copytooffsets", "--copyoffsets",  action='store_true', help="Copy priors to angle and distance offset parameters")

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':        args.input,
                 'output':       args.output,
                 'copypriors':   args.copypriors,
                 'copyoffsets':  args.copyoffsets,
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Copy offsets to priors overwrites
    if args_dict['copypriors']:
        add_column_parameters  = {'rlnAngleTilt':90,
                                  'rlnAngleTiltPrior':90,
                                  'rlnAnglePsiFlipRatio':0}

        del_column_parameters  = {'rlnAngleRotPrior':None}

        copy_column_parameters = {'rlnOriginX':  'rlnOriginXPrior',
                                  'rlnOriginY':  'rlnOriginYPrior',
                                  'rlnAnglePsi': 'rlnAnglePsiPrior',
                                  'rlnAngleTilt':'rlnAngleTiltPrior'}
        
    # Copy priors to offsets    
    if args_dict['copyoffsets']:
        add_column_parameters  = {}

        del_column_parameters  = {}

        copy_column_parameters = {'rlnOriginXPrior':  'rlnOriginX',
                                  'rlnOriginYPrior':  'rlnOriginY',
                                  'rlnAnglePsiPrior': 'rlnAnglePsi',
                                  'rlnAngleTiltPrior':'rlnAngleTilt'}

    # Create an EM project object
    new_project = em.Project(name='EMHelixReady')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Delete columns
    new_project.delete_columns(del_column_parameters)

    # Add new columns
    new_project.add_columns(add_column_parameters)

    # Add new columns
    new_project.copy_columns(copy_column_parameters)

    # Write output files
    new_project.write_output_files(write_ref_class_star=False)


if __name__ == "__main__":
    main()
