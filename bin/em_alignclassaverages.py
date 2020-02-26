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

    parser.add_argument("-i",            "--input",        type=str,  help="Class reference star file")
    parser.add_argument("-o",            "--output",       type=str,  help="Output directory", default=None)
    parser.add_argument("-use-unmasked", "--useunmasked",  action='store_true', help='Use unmasked classes for alignment')
    parser.add_argument("-apix",         "--apix",         type=float, help='Pixel size in anstroms', default=1.82)
    parser.add_argument("-diameter",     "--diameter",     type=float, help='Particle diameter in Anstroms', default=1000)
    parser.add_argument("-numiter",      "--numiter",      type=int,   help='Number of iterations', default=20)
    parser.add_argument("-sigma-psi",     "--sigmapsi",      type=float, help="Sigma-psi for alignment of classes", default=-1)
    parser.add_argument("-offset-range",  "--offsetrange",   type=int,   help="Offset range for alignment of classes", default=10)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'output':      args.output,
                 'useunmasked': args.useunmasked,
                 'apix':        args.apix,
                 'diameter':    args.diameter,
                 'numiter':     args.numiter,
                 'sigmapsi':    args.sigmapsi,
                 'offsetrange': args.offsetrange
                 }

    # Check if the input file exists
    if args_dict['input'] is None:
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.ProjectAlign2D(name='EMAlignClasses')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_class_refs(args_dict['input'])
    print('Read class reference file {}'.format(args_dict['input']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Set parameters
    new_project.set_params(args_dict['apix'], args_dict['diameter'])

    # Prepare project files
    new_project.prepare_project(use_unmasked_classes=args_dict['useunmasked'])

    # Normalize class averages
    new_project.normalize_class_refs()

    # Run relion
    new_project.set_relion_refine_args(offset_range=args_dict['offsetrange'], sigma_psi=args_dict['sigmapsi'], num_iter=args_dict['numiter'], firstiter_cc=False, T=4)
    new_project.run_refine2D()

    # Process 2D refine files
    new_project.prepare_refine2D()

    # Make transformed class stacks
    new_project.create_transformed_class_stacks()


if __name__ == "__main__":
    main()
