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

    parser.add_argument("-i",           "--input",         type=str,   help="Particle star file")
    parser.add_argument("-refclass",    "--refclass",      type=str,   help="Class star file")
    parser.add_argument("-o",           "--output",        type=str,   help="Output directory", default=None)
    parser.add_argument("-diameter",    "--diameter",      type=float, help="Particle diameter in Angstroms", default=None)
    parser.add_argument("-maskalign",   "--maskalign",     type=str,   help="Mask used for 2D classification", default=None)
    parser.add_argument("-refalign",    "--refalign",      type=str,   help="Reference mrc file used for alignment", default=None)
    parser.add_argument("-skip-rotate", "--skiprotate",    action='store_true',   help="Skip rotation in alignment of class averages to reference")
    parser.add_argument("-use-unmasked","--useunmasked",   action='store_true',   help="Use unmasked classes for alignment of classes")
    parser.add_argument("-sigma-psi",   "--sigmapsi",      type=float, help="Sigma-psi for alignment of classes", default=-1)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'refclass':      args.refclass,
                 'output':        args.output,
                 'diameter':      args.diameter,
                 'maskalign':     args.maskalign,
                 'refalign':      args.refalign,
                 'skiprotate':    args.skiprotate,
                 'useunmasked':   args.useunmasked,
                 'sigmapsi':      args.sigmapsi
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Check if the class reference file exists
    if args_dict['refclass'] is None and not os.path.isfile(args_dict['refclass']):
        parser.print_help()
        sys.exit('Class Reference file file does not exist!')

    # Create an EM project object
    new_project = em.ProjectAlign2D()
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    new_project.read_class_refs(args_dict['refclass'])
    print('Read class reference file {}'.format(args_dict['refclass']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Determine pixel size from particle star file
    new_project.read_particle_apix()

    # Set particle diameter
    new_project.set_particle_diameter(args_dict['diameter'])

    # Set alignment mask
    new_project.set_alignment_mask(args_dict['maskalign'])

    # Set alignment reference
    new_project.set_alignment_ref(args_dict['refalign'])

    # Prepare project files
    new_project.prepare_project(use_unmasked_classes=args_dict['useunmasked'])

    # Normalize class averages
    new_project.normalize_class_refs()

    # Run relion
    new_project.set_relion_refine_args(skip_rotate=args_dict['skiprotate'], sigma_psi=args_dict['sigmapsi'])
    new_project.run_refine2D()

    # Process 2D refine files
    new_project.prepare_refine2D()

    # Transform particles
    new_project.transform_particles()

    # Write output files
    new_project.write_output_files()

    # Make transformed class stacks
    new_project.create_transformed_class_stacks()


if __name__ == "__main__":
    main()
