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

    parser.add_argument("-i",            "--input",       type=str,            help="Particle star file")
    parser.add_argument("-o",            "--output",      type=str,            help="Output directory", default=None)
    parser.add_argument("-orig",         "--original",    type=str,            help="Original star file")
    parser.add_argument("-del-classes",  "--delclasses",  type=int, nargs='*', help="Classes to be deleted", default=None)
    parser.add_argument("-sel-classes",  "--selclasses",  type=int, nargs='*', help="Selected classes", default=None)
    parser.add_argument("-score-cutoff", "--scorecutoff", type=float,          help="Frealign score lower cutoff", default=None)
    parser.add_argument("-sigma-cutoff", "--sigmacutoff", type=float,          help="Frealign sigma noise lower cutoff", default=None)
    parser.add_argument("-ml-cutoff",    "--mlcutoff",    type=float,          help="Frealign maximum-likelihood lower cutoff", default=None)
    parser.add_argument("-db",           "--db",          type=str,            help="Cistem database file", default=None)
    parser.add_argument("-dont-copy",    "--dontcopy",    action='store_true',   help="Dont copy alignment offset and angle values")


    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':       args.input,
                 'original':    args.original,
                 'delclasses':  args.delclasses,
                 'selclasses':  args.selclasses,
                 'scorecutoff': args.scorecutoff,
                 'sigmacutoff': args.sigmacutoff,
                 'mlcutoff':    args.mlcutoff,
                 'db':          args.db,
                 'dontcopy':    args.dontcopy,
                 'output':      args.output}

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Check if the input file exists
    if args_dict['original'] is None or not os.path.isfile(args_dict['original']):
        parser.print_help()
        sys.exit('Original star file does not exist!')

    # Create an EM project object
    new_project = em.Cistem(name='EMPar2Star')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_particles(args_dict['original'])
    print('Read particle star file {}'.format(args_dict['original']))

    # Read apix from star file
    new_project.read_particle_apix()

    # Prepare output star files
    new_project.prepare_io_files_star()

    # Create write formatter
    new_project.create_write_formatter()

    # Read par file
    new_project.read_par(args_dict['input'])
    print('Read particle par file {}'.format(args_dict['input']))

    # Prepare project
    new_project.set_par2star_params(args_dict['delclasses'],
                                    args_dict['selclasses'],
                                    args_dict['scorecutoff'],
                                    args_dict['sigmacutoff'],
                                    args_dict['mlcutoff'])

    # Read cistem database
    new_project.read_db(args_dict['db'])

    # Select particles that exist in the refinement package
    new_project.select_particles()

    # Convert to par
    if not args_dict['dontcopy']:
        new_project.copy2star()

    # Write output file
    new_project.write_star_file()


if __name__ == "__main__":
    main()
