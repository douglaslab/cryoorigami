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

    parser.add_argument("-i",            "--input",         type=str,     help="Relion FSC xml file")
    parser.add_argument("-o",            "--output",        type=str,     help="Output directory", default=None)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'output':        args.output}

    # Create an EM project object
    new_project = em.ProjectPlot(name='EMFsc')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Read particles
    new_project.read_fsc(args_dict['input'])
    print('Read FSC xml file {}'.format(args_dict['input']))

    # Plot fsc
    new_project.plot_fsc()

    # Write output files
    new_project.write_fsc()


if __name__ == "__main__":
    main()
