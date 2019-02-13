#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 11:27:03
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$


import utilities as util

# Barcode functions for DNA origami assisted cryo-EM

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


def Frame_


def Frame_ptcl_angle(ptcl_star):
	'''
	Frame v3-7 barcode function
	'''
	barcode_dict = parse_barcode(ptcl_star)

	# Get bottom code
	# The code is in reverse order compared to the DNA origami design

	bottom_code  = 6-int(barcode_dict['bottom'])

	# Rotation angle pitch
	angle_pitch = 34.28

	# Rotation angle
	rot_angle  = bottom_code*angle_pitch if bottom_code < 4 else -(bottom_code-3)*34.28 

	# Tilt angle
	tilt_angle = 90.0 

	# Result dictionary
	result_dict = {'rlnAngleTilt': tilt_angle,
				   'rlnAngleRot': rot_angle}

	return result_dict

def Frame_angle(data_star):
	'''
	Frame v3-7 barcode function for data star
	'''
	
	# Start with empty lists
	tilt_angle_list = []
	rot_angle_list  = []

	for ptcl_index, ptcl_row in data_star.iterrows():

		barcode_dict = Framev3_7_ptcl_angle(ptcl_row)
		tilt_angle_list.append(barcode_dict['rlnAngleTilt'])
		tilt_angle_list.append(barcode_dict['rlnAngleRot'])

	# Assign the new values
	data_star['rlnAngleTilt'] = tilt_angle_list
	data_star['rlnAngleRot']  = rot_angle_list

	return data_star

