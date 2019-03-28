#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-15 11:00:21
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import argparse
import origamiem as em
import utilities as util


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-half1',   "--half1",   type=str, help='Half map 1')
    parser.add_argument('-half2',   "--half2",   type=str, help='Half map 2')
    parser.add_argument('-whole',   "--whole",   type=str, help='Whole map')
    parser.add_argument('-mask',    "--mask",    type=str,   default='',  help='Mask for FSC calculation')
    parser.add_argument('-apix',    "--apix",    type=float, default='',  help='Pixel size (Angstroms)')
    parser.add_argument('-bfactor', "--bfactor", type=float, default='',  help='B-factor for global sharpening (Negative)')
    parser.add_argument('-ncones',  "--ncones",  type=int,   default=500, help='Number of cones to use for dFSC calculation')
    parser.add_argument('-angle',   "--angle",   type=float, default=20,  help='Apex angle of cones for dFSC calculation')

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'half1':     args.half1,
                 'half2':     args.half2,
                 'whole':     args.whole,
                 'mask':      args.mask,
                 'apix':      args.apix,
                 'bfactor':   args.bfactor,
                 'ncones':    args.ncones,
                 'angle':     args.angle}

    # Create an EM project object
    new_project = em.ProjectFsc(name='EMFsc')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Set parameters
    new_project.set_params(args_dict['half1'], args_dict['half2'], args_dict['whole'], args_dict['mask'],
                           args_dict['apix'], args_dict['bfactor'], args_dict['ncones'], args_dict['angle'])

    # Prepare project
    new_project.prepare_project()

    # Prepare io files
    new_project.prepare_io_files()

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
