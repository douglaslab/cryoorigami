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

    parser.add_argument("-i",             "--input",         type=str,     help="Particle star file")
    parser.add_argument("-refclass",      "--refclass",      type=str,     help="Class star file")
    parser.add_argument("-o",             "--output",        type=str,     help="Output directory", default=None)
    parser.add_argument("-diameter",      "--diameter",      type=float,   help="Particle diameter in Angstroms", default=None)
    parser.add_argument("-masksolvent",   "--masksolvent",   type=str,     help="Solvent mask used for 2D classification", default=None)
    parser.add_argument("-masksub",       "--masksub",       type=str,     help="Mask used for the subtraction", default=None)
    parser.add_argument("-maskref",       "--maskref",       type=str,     help="Reference mask that defines the boundary of the structure", default=None)
    parser.add_argument("-batch",         "--batch",         type=int,     help="Particle batch size", default=100)
    parser.add_argument("-maxptcl",       "--maxptcl",       type=int,     help="Maximum number of particles to write", default=None)
    parser.add_argument("-method",        "--method",        type=str,     help="Particle subtraction method", choices=['subctf', 'cropctf', 'crop'], default='subctf')
    parser.add_argument("-norm",          "--norm",          type=str,     help="Normalization method for subtraction", choices=['ccc', 'intensity', 'frc'], default='frc')
    parser.add_argument("-subtractbg",    "--subtractbg",    action='store_true', help="Subtract background. For crop methods only.")
    parser.add_argument("-skipfirstpeak", "--skipfirstpeak", action='store_true', help="Skip first peak in CTF.")
    parser.add_argument("-innerdiameter", "--innerdiameter", type=float,   help="Inner diameter for the region of interest", default=None)
    parser.add_argument("-clip",          "--clip",          type=int,     help="Clip box size in pixels", default=None)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'refclass':      args.refclass,
                 'output':        args.output,
                 'diameter':      args.diameter,
                 'masksolvent':   args.masksolvent,
                 'masksub':       args.masksub,
                 'maskref':       args.maskref,
                 'batch':         args.batch,
                 'maxptcl':       args.maxptcl,
                 'method':        args.method,
                 'norm':          args.norm,
                 'subtractbg':    args.subtractbg,
                 'skipfirstpeak': args.skipfirstpeak,
                 'innerdiameter': args.innerdiameter,
                 'clip':          args.clip
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
    new_project = em.ProjectSubtract2D()
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
    new_project.set_alignment_mask(args_dict['masksolvent'])

    # Set subtraction mask
    new_project.set_subtraction_mask(args_dict['masksub'])

    # Set alignment reference
    new_project.set_structure_mask(args_dict['maskref'])

    # Set inner diameter
    new_project.set_inner_diameter(args_dict['innerdiameter'])

    # Set clip box
    new_project.set_clip_box(args_dict['clip'])

    # Prepare project files
    new_project.prepare_project()

    # Subtract class mrc from particle mrc
    new_project.subtract_class_mrc(max_ptcl=args_dict['maxptcl'],
                                   batch_size=args_dict['batch'],
                                   subtract_func=args_dict['method'],
                                   subtract_bg=args_dict['subtractbg'],
                                   norm_method=args_dict['norm'],
                                   skip_to_firstpeak=args_dict['skipfirstpeak'])

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
