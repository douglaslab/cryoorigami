#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 11:27:03
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$


import cryoorigami.utilities as util
import numpy as np


# Barcode functions for DNA origami assisted cryo-EM
def Framev60(barcode_num, offsetrot=0):
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
        rot_angle = rot_angle_map[barcode_num] + offsetrot
    else:
        rot_angle = 0.0

    return tilt_angle, util.euler360to180(rot_angle)


def Framev61(barcode_num, offsetrot=0):
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
        rot_angle = rot_angle_map[barcode_num] + offsetrot
    else:
        rot_angle = 0.0

    return tilt_angle, util.euler360to180(rot_angle)


def Framev60rev(barcode_num, offsetrot=0):
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
        rot_angle = -rot_angle_map[barcode_num] +  offsetrot
    else:
        rot_angle = 0.0

    return tilt_angle, util.euler360to180(rot_angle)

# Map functions
map_funcs = {'Framev60': Framev60,
			 'Framev60rev': Framev60rev,
			 'Framev61': Framev61}

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


def Frame_ptcl_angle(ptcl_star, offsetrot):
    '''
    Frame v3-7 barcode function
    '''
    barcode_dict = parse_barcode(ptcl_star)

    # Get bottom code
    bit_code  = int(barcode_dict['bit'])

    # Frame name
    map_func  = barcode_dict['name']

    # Get the angles
    tilt_angle, rot_angle = map_funcs[map_func](bit_code, offsetrot)

    # Result dictionary
    result_dict = {'rlnAngleTiltPrior': tilt_angle,
                   'rlnAngleRotPrior': rot_angle}

    return result_dict


def Frame_angle(data_star, offsetrot):
    '''
    Frame v3-7 barcode function for data star
    '''

    # Start with empty lists
    tilt_angle_list = []
    rot_angle_list  = []

    for ptcl_index, ptcl_row in data_star.iterrows():

        barcode_dict = Frame_ptcl_angle(ptcl_row, offsetrot)
        tilt_angle_list.append(barcode_dict['rlnAngleTiltPrior'])
        rot_angle_list.append(barcode_dict['rlnAngleRotPrior'])

    # Assign the new values
    data_star['rlnAngleTiltPrior'] = np.array(tilt_angle_list)
    data_star['rlnAngleRotPrior']  = np.array(rot_angle_list)

    return data_star
