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

    parser.add_argument("-i",     "--input",       type=str, nargs='+', help="Star file(s)")
    parser.add_argument("-o",     "--output",      type=str, help="Output directory", default=None)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'output':      args.output
                 }

    # Check if the input file exists
    if args_dict['input'] is None:
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='EMMergeStar')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particle_stars(args_dict['input'])
    print('Read particle star files {}'.format(','.join(args_dict['input'])))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
