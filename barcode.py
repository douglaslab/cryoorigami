#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 11:27:03
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$


import utilities as util
import numpy as np


# Barcode functions for DNA origami assisted cryo-EM
def Frame_map(barcode_num):
    # Rotation angle pitch
    angle_pitch = 34.28

    rot_angle_map = {1:  0,
                     2: -1*angle_pitch,
                     3: -2*angle_pitch,
                     4: -3*angle_pitch,
                     5:  1*angle_pitch,
                     6:  2*angle_pitch,
                     7:  3*angle_pitch}
    # Tilt angle
    tilt_angle = 90.0

    if barcode_num in rot_angle_map:
        rot_angle = rot_angle_map[barcode_num]
    else:
        rot_angle = 0.0

    return tilt_angle, rot_angle


def parse_barcode(ptcl_star):
    '''
    Parse barcode
    '''
    barcode_dict = {}
    if 'rlnParticleName' in ptcl_star:
        barcode_arr = ptcl_star['rlnParticleName'].strip().split(',')

        # Get barcode dictionary
        barcode_dict = util.parse_star_parameters(barcode_arr)

    return barcode_dict


def Frame_ptcl_angle(ptcl_star):
    '''
    Frame v3-7 barcode function
    '''
    barcode_dict = parse_barcode(ptcl_star)

    # Get bottom code
    # The code is in reverse order compared to the DNA origami design (e.g. 0 is actually 6, 6 is actually 0)

    bottom_code  = int(barcode_dict['rot'])

    tilt_angle, rot_angle = Frame_map(bottom_code)

    # Result dictionary
    result_dict = {'rlnAngleTiltPrior': tilt_angle,
                   'rlnAngleRotPrior': rot_angle}

    return result_dict


def Frame_angle(data_star):
    '''
    Frame v3-7 barcode function for data star
    '''

    # Start with empty lists
    tilt_angle_list = []
    rot_angle_list  = []

    for ptcl_index, ptcl_row in data_star.iterrows():

        barcode_dict = Frame_ptcl_angle(ptcl_row)
        tilt_angle_list.append(barcode_dict['rlnAngleTiltPrior'])
        rot_angle_list.append(barcode_dict['rlnAngleRotPrior'])

    # Assign the new values
    data_star['rlnAngleTiltPrior'] = np.array(tilt_angle_list)
    data_star['rlnAngleRotPrior']  = np.array(rot_angle_list)

    return data_star
