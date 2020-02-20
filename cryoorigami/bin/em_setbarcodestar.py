#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-15 11:00:21
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import sys
import argparse
import cryoorigami.origamiem as em
import cryoorigami.utilities as util


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i",        "--input",       type=str, help="Particle star file")
    parser.add_argument("-o",        "--output",      type=str, help="Output directory", default=None)
    parser.add_argument("-barcode",  "--barcode",     type=str, help="Barcode to set", nargs='*', default={})

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'output':      args.output,
                 'barcode':     args.barcode
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Get the new column parameters
    if args_dict['barcode'] is not None:
        barcode_parameters = util.parse_star_parameters(args_dict['barcode'])
    else:
        barcode_parameters = {}


    # Create an EM project object
    new_project = em.Project(name='EMSetBarcode')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Add new columns
    new_project.set_particle_barcode(barcode_parameters)

    # Write output files
    new_project.write_output_files(write_ref_class_star=False)


if __name__ == "__main__":
    main()
