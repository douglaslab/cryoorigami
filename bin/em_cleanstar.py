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

    parser.add_argument("-i",    "--input",       type=str, help="Particle star file")
    parser.add_argument("-o",    "--output",      type=str, help="Output directory", default=None)
    parser.add_argument("-mic",  "--micrograph",  type=str, help="Micrograph star file")

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':        args.input,
                 'output':       args.output,
                 'micrograph':   args.micrograph
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Check if the reference file exists
    if args_dict['micrograph'] is None and not os.path.isfile(args_dict['micrograph']):
        parser.print_help()
        sys.exit('Micrographs file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='EMClean')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Read particle pixel size information
    new_project.read_particle_apix()

    new_project.read_micrographs(args_dict['micrograph'])
    print('Read micrographs star file {}'.format(args_dict['micrograph']))

    # Read micrograph pixel size information
    new_project.read_mic_apix()

    # Read micrograph header
    new_project.read_mic_header()

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Check particle positions
    new_project.check_particle_pos()

    # Delete outside particles
    new_project.delete_outside_ptcls()

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
