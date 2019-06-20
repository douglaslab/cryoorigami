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

    parser.add_argument("-i",            "--input",         type=str,   help="Input particle star file")
    parser.add_argument("-ref",          "--reference",     type=str,   help="Reference star file")
    parser.add_argument("-addsuffix",    "--addsuffix",     type=str,   nargs='*', help="Suffix to add to the image names", default=None)
    parser.add_argument("-remsuffix",    "--removesuffix",  type=str,   help="Suffix to remove from the image names", default=None)
    parser.add_argument("-o",            "--output",        type=str,   help="Output directory", default=None)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':        args.input,
                 'reference':    args.reference,
                 'addsuffix':    args.addsuffix,
                 'removesuffix': args.removesuffix,
                 'output':       args.output}

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.ProjectImgName(name='EMReplaceImgName')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    if args_dict['addsuffix'] is not None:
        new_project.add_suffix(args_dict['addsuffix'])
    elif args_dict['removesuffix'] is not None:
        new_project.remove_suffix(args_dict['removesuffix'])
    else:
        # Check if the reference file exists
        if args_dict['reference'] is None or not os.path.isfile(args_dict['reference']):
            parser.print_help()
            sys.exit('Reference file does not exist!')

        # Read particles
        new_project.read_reference(args_dict['reference'])
        print('Read reference star file {}'.format(args_dict['reference']))

        # Prepare project
        new_project.replace_particle_img()
        
    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
