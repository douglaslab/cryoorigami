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

    parser.add_argument("-i",             "--input",        type=str, help="Particle star file")
    parser.add_argument("-o",             "--output",       type=str, help="Output directory", default=None)
    parser.add_argument("-cols",          "--columns",      type=str, help="Columns to copy", nargs='*', default=None)
    parser.add_argument("-copytopriors",  "--copypriors",   type=str, help="Copy offset and/or angle parameters to priors", default=None, choices=['angle','all', 'helix'])
    parser.add_argument("-copytooffsets", "--copyoffsets",  type=str, help="Copy priors to angle and distance offset parameters", default=None, choices=['angle','all', 'helix'])
    parser.add_argument("-reset-priors",  "--resetpriors",  action='store_true', help="Delete prior offset and angle columngs")
    parser.add_argument("-reset-offsets", "--resetoffsets", action='store_true', help="Assign angle and origin offsets to 0")
    parser.add_argument("-invert-psi",    "--invertpsi",    action='store_true', help="Invert psi angle")
    parser.add_argument("-invert-origin", "--invertorigin", action='store_true', help="Invert originX and originY")

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':        args.input,
                 'output':       args.output,
                 'columns':      args.columns,
                 'copypriors':   args.copypriors,
                 'copyoffsets':  args.copyoffsets,
                 'resetpriors':  args.resetpriors,
                 'resetoffsets': args.resetoffsets,
                 'invertpsi':    args.invertpsi,
                 'invertorigin': args.invertorigin
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Get the new column parameters
    if args_dict['columns'] is not None:
        new_column_parameters = util.parse_star_parameters(args_dict['columns'])
    else:
        new_column_parameters = None

    # Copy offsets to priors overwrites
    if args_dict['copypriors'] == 'all':
        new_column_parameters = {'rlnOriginX':  'rlnOriginXPrior',
                                 'rlnOriginY':  'rlnOriginYPrior',
                                 'rlnAnglePsi': 'rlnAnglePsiPrior',
                                 'rlnAngleRot': 'rlnAngleRotPrior',
                                 'rlnAngleTilt':'rlnAngleTiltPrior'}
    elif args_dict['copypriors'] == 'angle':
        new_column_parameters = {'rlnAnglePsi': 'rlnAnglePsiPrior',
                                 'rlnAngleRot': 'rlnAngleRotPrior',
                                 'rlnAngleTilt':'rlnAngleTiltPrior'}
    elif args_dict['copypriors'] == 'helix':
        new_column_parameters = {'rlnOriginX':  'rlnOriginXPrior',
                                 'rlnOriginY':  'rlnOriginYPrior',
                                 'rlnAnglePsi': 'rlnAnglePsiPrior',
                                 'rlnAngleTilt':'rlnAngleTiltPrior'}

    # Copy priors to offsets overwrites
    if args_dict['copyoffsets'] == 'all':
        new_column_parameters = {'rlnOriginXPrior':  'rlnOriginX',
                                 'rlnOriginYPrior':  'rlnOriginY',
                                 'rlnAnglePsiPrior': 'rlnAnglePsi',
                                 'rlnAngleRotPrior': 'rlnAngleRot',
                                 'rlnAngleTiltPrior':'rlnAngleTilt'}
    elif args_dict['copyoffsets'] == 'angle':
        new_column_parameters = {'rlnAnglePsiPrior': 'rlnAnglePsi',
                                 'rlnAngleRotPrior': 'rlnAngleRot',
                                 'rlnAngleTiltPrior':'rlnAngleTilt'}
    elif args_dict['copyoffsets'] == 'helix':
        new_column_parameters = {'rlnOriginXPrior':  'rlnOriginX',
                                 'rlnOriginYPrior':  'rlnOriginY',
                                 'rlnAnglePsiPrior': 'rlnAnglePsi',
                                 'rlnAngleTiltPrior':'rlnAngleTilt'}

    # Create an EM project object
    new_project = em.Project(name='EMCopyColumns')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # If reset priors option is ON
    if args_dict['resetpriors']:
        new_project.reset_priors()

    # If reset offset option is on
    if args_dict['resetoffsets']:
        new_project.reset_offsets()

    # If invert-psi option is on
    if args_dict['invertpsi']:
        new_project.invert_psi()

    # If invert-origin
    if args_dict['invertorigin']:
        new_project.invert_origin()

    # Add new columns
    new_project.copy_columns(new_column_parameters)

    # Write output files
    new_project.write_output_files(write_ref_class_star=False)


if __name__ == "__main__":
    main()
