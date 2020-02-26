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

    parser.add_argument("-i",              "--input",           type=str,            help="Particle star file")
    parser.add_argument("-invert-angles",  "--invertangles",    action='store_true', help="Invert eular angles")
    parser.add_argument("-invert-origin",  "--invertorigin",    action='store_true', help="Invert origins(shifts)")
    parser.add_argument("-o",              "--output",          type=str,            help="Output directory", default=None)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'invertangles':  args.invertangles,
                 'invertorigin':  args.invertorigin,
                 'output':        args.output}

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.Cistem(name='EMStar2Par')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Sort particles based on imagename
    new_project.sort_images()

    # Prepare project
    new_project.prepare_project()

    # If invert-origin
    if args_dict['invertorigin']:
        new_project.invert_origin()

    # If invert-origin
    if args_dict['invertangles']:
        new_project.invert_angles()

    # Convert to par
    new_project.convert2par()

    # Write output file
    new_project.write_output_file()

    # Write also star file
    new_project.write_star_file()


if __name__ == "__main__":
    main()
