#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 13:35:55
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import yaml
import numpy as np


def write_config_file(fname, args_dict):
    '''
    Dump parameters into ymal config file
    '''
    with open(fname, 'w') as outfile:
        yaml.dump(args_dict, outfile, default_flow_style=False)


def parse_star_parameters(parameter_input):
    '''
    Parse  parameter input
    '''
    # Get the parameter groups
    groups = parameter_input.split('.')

    # Get the functions and its parameters
    parameter_lists = [group.split(':') for group in groups]

    # Create the parameters dictionary
    parameter_dict = {}
    for parameter_list in parameter_lists:
        column_name = None
        if len(parameter_list) > 0:
            column_name = parameter_list[0]
            parameter_dict[column_name] = None
        else:
            continue
        if len(parameter_list) > 1:
            parameter_dict[column_name] = parameter_list[1]

    return parameter_dict


def euler2rot2D(alpha):
    '''
    Eular angle to 2D rotation matrix
    '''
    ca = np.cos(alpha*np.pi/180)
    sa = np.sin(alpha*np.pi/180)

    r = np.array([[ca, -sa],
                  [+sa, ca]])
    return r
