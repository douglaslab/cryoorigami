#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2019-01-22 17:48:26
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

    parser.add_argument("-flip",     "--flip",      type=str,   help="Flip particle star file")
    parser.add_argument("-noflip",   "--noflip",    type=str,   help="No flip particle star file")
    parser.add_argument("-o",        "--output",    type=str,   help="Output directory", default=None)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'flip':     args.flip,
                 'noflip':   args.noflip,
                 'output':   args.output}

    # Check if the input file exists
    if args_dict['flip'] is None or not os.path.isfile(args_dict['flip']):
        parser.print_help()
        sys.exit('Flip file does not exist!')

    # Check if the input file exists
    if args_dict['noflip'] is None or not os.path.isfile(args_dict['noflip']):
        parser.print_help()
        sys.exit('Noflip file does not exist!')

    # Create an EM project object
    new_project = em.ProjectFlip(name='EMParticleMergeFlip')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Prepare project
    new_project.prepare_merge_project(args_dict['flip'], args_dict['noflip'])


if __name__ == "__main__":
    main()
