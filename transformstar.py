#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-11-15 00:07:10
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import argparse
import origamiem as em
import utilities as util


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i",       "--input",       type=str, help="Particle star file")
    parser.add_argument("-o",       "--output",      type=str, help="Output directory", default=None)
    parser.add_argument("-ref",     "--reference",   type=str, help="Reference class star file", default=None)
    parser.add_argument("-addcols", "--addcolumns",  type=str, help="Add new columns with same value for all particles", default=None)
    parser.add_argument("-recenter","--recenter",    action='store_true', help="Recenter particles")
    
    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'output':      args.output,
                 'reference':   args.reference,
                 'addcolumns':  args.addcolumns,
                 'recenter':    args.recenter
                 }


    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Get the new column parameters
    if not args_dict['addcolumns'] is None:
        new_column_parameters = parse_star_parameters(args_dict['addcolumns'])
    else:
        new_column_parameters = None

    # Create an EM project object
    new_project = em.Project()
    new_project.set_output_directory(args_dict['input'], args_dict['output'])
    
    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Check if the input file exists
    if args_dict['reference'] is not None and os.path.isfile(args_dict['reference']):
        new_project.read_class_refs(args_dict['reference'])
        print('Read class reference file {}'.format(args_dict['reference']))

    # Prepare input and output files
    new_project.prepare_io_files()

    # Perform reference based particle star transformation
    new_project.transform_particles()

    # Add new columns
    new_project.add_columns(new_column_parameters)

    # Recenter
    if args_dict['recenter']:
        new_project.recenter_particles()

    # Write output files
    new_project.write_output_files()

if __name__ == "__main__":
    main()