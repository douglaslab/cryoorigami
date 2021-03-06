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

    parser.add_argument("-i",           "--input",          type=str, help="Particle star file")
    parser.add_argument("-o",           "--output",         type=str, help="Output directory", default=None)
    parser.add_argument("-selbarcode",  "--selbarcode",     type=str, help="Barcode to select", nargs='*', default={})
    parser.add_argument("-delbarcode",  "--delbarcode",     type=str, help="Barcode to delete", nargs='*', default={})

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'output':      args.output,
                 'selbarcode':  args.selbarcode,
                 'delbarcode':  args.delbarcode
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='EMSelectBarcode')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Select barcode
    new_project.select_by_barcode(args_dict['selbarcode'])

    # Delete barcode
    new_project.delete_by_barcode(args_dict['delbarcode'])

    # Write output files
    new_project.write_output_files(write_ref_class_star=False)


if __name__ == "__main__":
    main()
