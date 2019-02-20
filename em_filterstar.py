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
    parser.add_argument("-maxprob",      "--maxprob",       type=float,   help="Maximum probability", default=0.5)
    parser.add_argument("-maxclass",     "--maxclass",      type=int,     help="Maximum number of classes assigned for a particle", default=10)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'output':        args.output,
                 'maxprob':       args.maxprob,
                 'maxclass':      args.maxclass}

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='EMFilter2D')
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
    new_project.filter_ptcls(maxprob=args_dict['maxprob'], maxclass=args_dict['maxclass'])

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()