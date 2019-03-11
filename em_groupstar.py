#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-17 16:44:31
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

    parser.add_argument("-i",         "--input",       type=str, help="Particle star file")
    parser.add_argument("-o",         "--output",      type=str, help="Output directory", default=None)
    parser.add_argument("-maxmics",   "--maxmics",     type=int, help="Maximum number of micrographs allowed in a group", default=50)
    parser.add_argument("-method",    "--method",      type=str, help="Grouping method", choices=['defocus','intensity'], default='defocus')
    parser.add_argument("-ref",       "--reference",          type=str, help="Reference file to make the groups")
    parser.add_argument("-threshdef" ,"--thresholddefocus",   type=float, help="Defocus threshold value for grouping (Angstrom)", default=100)
    parser.add_argument("-threshint" ,"--thresholdintensity", type=float, help="Intensity scale threshold value for grouping",    default=0.1)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':      args.input,
                 'output':     args.output,
                 'maxmics':    args.maxmics,
                 'method':     args.method,
                 'reference':  args.reference,
                 'thresholddefocus':   args.thresholddefocus,
                 'thresholdintensity': args.thresholdintensity
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Check if the reference file exists
    if args_dict['reference'] is None and not os.path.isfile(args_dict['reference']):
        parser.print_help()
        sys.exit('Reference file does not exist!')

    # Create an EM project object
    new_project = em.ProjectGroup(name='EMGroup')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Read particles
    new_project.read_micrographs(args_dict['reference'])
    print('Read micrographs star file {}'.format(args_dict['reference']))

    # Set project parameters
    new_project.set_params(maxmics=args_dict['maxmics'], threshdef=args_dict['thresholddefocus'], threshint=args_dict['thresholdintensity'])

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Group particles based on a method
    new_project.group_micrographs()

    # Assign the groups
    new_project.assign_groups()

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
