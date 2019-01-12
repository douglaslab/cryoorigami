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

    parser.add_argument("-i",        "--input",         type=str, help="Particle star file")
    parser.add_argument("-o",        "--output",        type=str, help="Output directory", default=None)
    parser.add_argument("-ref",      "--reference",     type=str, help="Reference class star file")
    parser.add_argument("-diamater", "--diameter",      type=int, help="Particle diameter in Angstroms", default=None)
    parser.add_argument("-cmask",    "--classmask",     type=str, help="Mask used for 2D classification", default=None)
    parser.add_argument("-smask",    "--submask",       type=str, help="Mask used for subtraction", default=None)
    parser.add_argument("-norm",     "--normalization", type=str, help="Normalization procedure to use for subtraction",
                        choices=['frc', 'minmax'], default='frc')

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'output':        args.output,
                 'reference':     args.reference,
                 'diameter':      args.diameter,
                 'classmask':     args.classmask,
                 'submask':       args.submask,
                 'normalization': args.normalization
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Check if the reference file exists
    # if args_dict['reference'] is None and not os.path.isfile(args_dict['reference']):
    #    parser.print_help()
    #    sys.exit('Reference file file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='ProjectSubtract2D')
    new_project.set_output_directory(args_dict['input'], args_dict['output'])

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))
    '''
    new_project.read_class_refs(args_dict['reference'])
    print('Read class reference file {}'.format(args_dict['reference']))

    # Prepare input and output files
    new_project.prepare_io_files_star()
    '''
    new_project.read_particle_mrc()
    # Add new columns
    # new_project.subtract2D(args_dict['diameter'], args_dict['classmask'], args_dict['submask'], args_dict['normalization'])

    # Write output files
    # new_project.write_output_files(write_ref_class_star=False)

if __name__ == "__main__":
    main()
