#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-17 16:44:31
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

    parser.add_argument("-i",        "--input",   type=str, help="Particle star file")
    parser.add_argument("-o",        "--output",  type=str, help="Output directory", default=None)
    parser.add_argument("-orimic",   "--orimic",  type=str, help="Original micrograph star file")
    parser.add_argument("-newmic",   "--newmic",  type=str, help="New micrograph star file")

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':    args.input,
                 'output':   args.output,
                 'orimic':   args.orimic,
                 'newmic':   args.newmic
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Check if the reference file exists
    if args_dict['orimic'] is None and not os.path.isfile(args_dict['orimic']):
        parser.print_help()
        sys.exit('Original micrographs file does not exist!')

    # Check if the reference file exists
    if args_dict['newmic'] is None and not os.path.isfile(args_dict['newmic']):
        parser.print_help()
        sys.exit('New micrographs file does not exist!')

    # Create an EM project object
    new_project = em.ProjectScale(name='EMScale')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Read particle pixel size information
    new_project.read_particle_apix()

    new_project.read_micrographs(args_dict['orimic'])
    print('Read original micrographs star file {}'.format(args_dict['orimic']))

    # Read original micrographs apix
    new_project.read_mic_apix()

    new_project.read_new_micrographs(args_dict['newmic'])
    print('Read new micrographs star file {}'.format(args_dict['newmic']))

    # Read new micrographs apix
    new_project.read_new_mic_apix()

    # Perform the transformation
    new_project.run()

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
