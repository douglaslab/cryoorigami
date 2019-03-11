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

    parser.add_argument("-i",            "--input",         type=str,     nargs='+', help="Particle star file")
    parser.add_argument("-o",            "--output",        type=str,     help="Output directory", default=None)
    parser.add_argument("-cols",         "--columns",       type=str,     nargs='+', default=['rlnMaxValueProbDistribution',
                                                                                              'rlnNrOfSignificantSamples'])

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'input':         args.input,
                 'output':        args.output}

    # Create an EM project object
    new_project = em.ProjectPlot(name='EMPlott')
    new_project.set_output_directory(args_dict['output'], project_root='.')

    # Write parameters to args filename
    args_filename = new_project.output_directory+'/args.yaml'
    util.write_config_file(args_filename, args_dict)

    # Run the project
    new_project.run(args_dict['input'])

    # Write output files
    new_project.write_output_files()


if __name__ == "__main__":
    main()
