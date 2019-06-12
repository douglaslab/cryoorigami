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

    parser.add_argument("-i",      "--input",       type=str, help="Particle star file")
    parser.add_argument("-o",      "--output",      type=str, help="Output directory", default=None)
    parser.add_argument("-Zflip",  "--Zflip",       action='store_true', help="Z-Flip particles in star file")
    parser.add_argument("-tilt90", "--tilt90",      action='store_true', help="Set tilt of particles to 90-degree side")

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'Zflip':       args.Zflip,
                 'tilt90':      args.tilt90,
                 'output':      args.output
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='EMRotate')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Make tilt angle around 90 degree
    if args_dict['tilt90']:
        new_project.tilt90()

    # Z-flip particles
    if args_dict['Zflip']:
        new_project.Zflip()

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
