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
    parser.add_argument("-transform",    "--transform", action='store_true', help='Transform the images before writing the stacks')
    parser.add_argument("-hp",           "--highpass",      type=float, help="Highpass filter in Angstrom units", default=None)
    parser.add_argument("-lp",           "--lowpass",       type=float, help="Lowpass filter in Angstrom units", default=None)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'output':        args.output,
                 'batch':         args.batch,
                 'transform':     args.transform,
                 'highpass':      args.highpass,
                 'lowpass':       args.lowpass
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.ProjectStack(name='EMStack')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Prepare project
    new_project.prepare_project(args_dict['highpass'], args_dict['lowpass'])

    # Flip particles
    new_project.create_stack(args_dict['batch'], args_dict['transform'])


if __name__ == "__main__":
    main()