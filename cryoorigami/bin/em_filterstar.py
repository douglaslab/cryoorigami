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

    parser.add_argument("-i",            "--input",         type=str,     help="Particle star file")
    parser.add_argument("-o",            "--output",        type=str,     help="Output directory", default=None)
    parser.add_argument("-maxprob",      "--maxprob",       type=float,   help="Maximum probability", default=None)
    parser.add_argument("-maxclass",     "--maxclass",      type=int,     help="Maximum number of classes assigned for a particle", default=None)
    parser.add_argument("-tilt",         "--tilt",          type=float,   nargs=2, help="Accepted range of tilt angles. [ min ,max]",       default=[0, 360])
    parser.add_argument("-dtilt",        "--dtilt",         type=float,   nargs=2, help="Accepted range of diff-tilt angles. [ min ,max]",  default=[0, 360])
    parser.add_argument("-dpsi",         "--dpsi",          type=float,   nargs=2, help="Accepted range of diff-psi angles. [ min ,max]",   default=[0, 360])
    parser.add_argument("-dalign",       "--dalign",        type=float,   nargs=2, help="Accepted range of diff-align angles. [ min ,max]", default=[-1, 1])
    
    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'output':        args.output,
                 'maxprob':       args.maxprob,
                 'maxclass':      args.maxclass,
                 'tilt':          args.tilt,
                 'dtilt':         args.dtilt,
                 'dpsi':          args.dpsi,
                 'dalign':        args.dalign}

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='EMFilter')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['input'])
    print('Read particle star file {}'.format(args_dict['input']))

    # Prepare input and output files
    new_project.prepare_io_files_star()

    # Determine pixel size from particle star file
    new_project.read_particle_apix()

    # Subtract class mrc from particle mrc
    new_project.filter_ptcls(maxprob=args_dict['maxprob'], 
                             maxclass=args_dict['maxclass'],
                             tilt_range=args_dict['tilt'],
                             dpsi_range=args_dict['dpsi'],
                             dtilt_range=args_dict['dtilt'],
                             dalign_range=args_dict['dalign'])

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
