#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2019-01-22 17:48:26
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

    parser.add_argument("-i",            "--input",         type=str,   help="Particle star file")
    parser.add_argument("-o",            "--output",        type=str,   help="Output directory", default=None)
    parser.add_argument("-batch",        "--batch",         type=int,   help="Particle batch size", default=100)
    parser.add_argument("-relion",       "--relion",        action='store_true',   help="Use relion for flipping")

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'output':        args.output,
                 'batch':         args.batch,
                 'relion':        args.relion
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.ProjectFlip()
    new_project.set_output_directory(args_dict['input'], args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Determine pixel size from particle star file
    new_project.read_particle_apix()

    # Prepare project
    new_project.prepare_project_relion()

    # Flip particles
    new_project.flip_particles_relion()


if __name__ == "__main__":
    main()
