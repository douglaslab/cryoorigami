#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-19 10:22:08
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

    parser.add_argument("-i",        "--input",       type=str, help="Particle cs file")
    parser.add_argument("-ref",      "--reference",   type=str, help="Template cs file", default=None)
    parser.add_argument("-pass",     "--passthrough", type=str, help="Passthrough cs file from cryosparc", default=None)
    parser.add_argument("-o",        "--output",      type=str, help="Output directory", default=None)
    parser.add_argument("-orig",     "--original",    type=str, help="Original star file", default=None)
    parser.add_argument("-micpath",  "--micpath",     type=str, help="Micrographs path", default="Micrographs")
    parser.add_argument("-projpath", "--projpath",    type=str, help="Cryosparc project path", default="")
    parser.add_argument("-delclass", "--delclass",    type=int, nargs='+', help="Classes to delete", default=[])
    parser.add_argument("-delstr",   "--delstr",      type=str, help="String to delete from Micrograph name", default='')
    parser.add_argument("-restore-offsets", "--restoreoffsets", action='store_true', help="Restore offsets from original star file")
    parser.add_argument("-merge-orig", "--mergeoriginal", action='store_true', help="Merge with original star file")

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':          args.input,
                 'reference':      args.reference,
                 'passthrough':    args.passthrough,
                 'output':         args.output,
                 'original':       args.original,
                 'micpath':        args.micpath,
                 'projpath':       args.projpath,
                 'delclass':       args.delclass,
                 'delstr':         args.delstr,
                 'restoreoffsets': args.restoreoffsets,
                 'mergeoriginal':  args.mergeoriginal
                 }

    # Check if the input file exists
    if args_dict['input'] is None or not os.path.isfile(args_dict['input']):
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Create an EM project object
    new_project = em.Project(name='EMcsparc2star')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Set cs files
    new_project.set_cs_files(blob_cs_file=args_dict['input'],
                             passthrough_cs_file=args_dict['passthrough'],
                             original_star_file=args_dict['original'],
                             ref_class_cs_file=args_dict['reference'])

    # Read cs files
    new_project.read_cs_files()

    # Prepare input and output files
    new_project.prepare_io_files_cs()

    # Conert cs to star
    new_project.convert_cs2star(args_dict['micpath'], args_dict['projpath'], args_dict['delclass'], args_dict['delstr'],
                                args_dict['restoreoffsets'], args_dict['mergeoriginal'])

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
