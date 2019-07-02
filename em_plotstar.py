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

    parser.add_argument("-i",            "--input",         type=str,     help="Particle star file")
    parser.add_argument("-o",            "--output",        type=str,     help="Output directory", default=None)
    parser.add_argument("-cols",         "--columns",       type=str,     nargs='+', default=['rlnDefocusU',
                                                                                              'rlnDefocusV',
                                                                                              'rlnDefocusAngle',
                                                                                              'rlnCtfFigureOfMerit',
                                                                                              'rlnFinalResolution',
                                                                                              'rlnMaxValueProbDistribution',
                                                                                              'rlnNrOfSignificantSamples',
                                                                                              'rlnAngleRot',
                                                                                              'rlnAngleTilt',
                                                                                              'rlnAnglePsi',
                                                                                              'rlnOriginX',
                                                                                              'rlnOriginY'])

    parser.add_argument("-pairs",        "--pairs",         type=str,     nargs='+', default=['rlnCtfFigureOfMerit:rlnFinalResolution',
                                                                                              'rlnAngleTilt:rlnAngleRot',
                                                                                              'rlnMaxValueProbDistribution:rlnNrOfSignificantSamples',
                                                                                              'rlnOriginX:rlnOriginY',
                                                                                              'rlnAngleRot:rlnAngleRotPrior',
                                                                                              'rlnAngleTilt:rlnAngleTiltPrior',
                                                                                              'rlnAnglePsi:rlnAnglePsiPrior'])
    parser.add_argument("-diffs",        "--diffs",         type=str,     nargs='+', default=['rlnOriginX:rlnOriginXPrior',
                                                                                              'rlnOriginY:rlnOriginYPrior',
                                                                                              'rlnAngleRot:rlnAngleRotPrior',
                                                                                              'rlnAngleTilt:rlnAngleTiltPrior',
                                                                                              'rlnAnglePsi:rlnAnglePsiPrior'])

    parser.add_argument("-orientation",  "--orientation",   action='store_true', help="Plot orientation of the particles with respect to priors")

    parser.add_argument("-nbins",        "--nbins",         type=int,     default=20)
    parser.add_argument("-fontsize",     "--fontsize",      type=int,     default=5)
    parser.add_argument("-format",       "--format",        type=str,     default='svg', choices=['png','svg'])

    parser.add_argument("-ref",          "--reference",     type=str,   help="Reference star file", default=None)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'output':        args.output,
                 'columns':       args.columns,
                 'pairs':         args.pairs,
                 'diffs':         args.diffs,
                 'orientation':   args.orientation,
                 'fontsize':      args.fontsize,
                 'reference':     args.reference,
                 'nbins':         args.nbins,
                 'format':        args.format}

    # Create an EM project object
    new_project = em.ProjectPlot(name='EMPlot')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Read particle mrc paths
    new_project.read_ptcl_mrc_paths()

    # Prepare io files
    new_project.prepare_io_files(args_dict['format'])

    # Set tick fontsizes
    new_project.set_tick_fontsize(size=args_dict['fontsize'])

    # Read references if it is provided
    if args_dict['reference'] is not None and os.path.isfile(args_dict['reference']):
      # Read particles
      new_project.read_reference(args_dict['reference'])
      print('Read reference star file {}'.format(args_dict['reference']))

      # Run reference project
      new_project.run_ref(args_dict['pairs'], args_dict['diffs'])
    else:
      # Run ptcl project
      new_project.run_ptcl(args_dict['columns'], args_dict['pairs'], args_dict['diffs'], args_dict['orientation'], args_dict['nbins'])

    # Write output files
    new_project.write_output_files(args_dict['format'])


if __name__ == "__main__":
    main()
