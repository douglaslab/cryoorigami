#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2017-01-07 20:46:38
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import pyORIGAMImodel as models
import numpy as np
import argparse
import matplotlib.pyplot as py
from EMAN2 import *

#Definition of the program
usage            = ['%(prog)s --min 20 --max 50 --scale 1',
                    '------------------------------------------------------------------------',
                    'DNA Origami Hinge model generator',
                    '08/23/2017',
                    'Tural Aksel',
                    '------------------------------------------------------------------------'
                    ]

#Create the parser
parser = argparse.ArgumentParser(description='',usage='\n'.join(usage)) 
parser.add_argument('--min'          ,default = 20,   type=int,   help='Minimum hinge angle (Default:20)')
parser.add_argument('--max'          ,default = 61,   type=int,   help='Maximum hinge angle (Default:50)')
parser.add_argument('--step'         ,default = 1,    type=int,   help='Angle step (Default:1)')
parser.add_argument('--apix'         ,default = 2.448,type=float, help='Scale (Default:2.448)')
parser.add_argument('--boxsize'      ,default = 360,  type=int,   help='Box size (Default:360)')
parser.add_argument("--invert"       ,action  = 'store_true',     help='Invert the image (Default:False)',default=False)
parser.add_argument('--output'       ,default = '.',  type=str,   help='Output folder (Default:.)')
parser.add_argument('--offset'       ,default = -10,  type=int,   help='Barcode offset (Default:-10)')

logger=E2init(sys.argv, -1)

args = parser.parse_args()
args_dict = vars(args)
parser.print_help()

#Get the argument
min_angle    = args_dict['min']
max_angle    = args_dict['max']
angle_step   = args_dict['step']
apix         = args_dict['apix']
output       = args_dict['output']
boxsize      = args_dict['boxsize']
invert       = args_dict['invert']
offset       = args_dict['offset']

if not os.path.isdir(output):
	print "Output directory doesn't exist"
	sys.exit(1)

#HINGE CLASS BIT COMBINATIONS
CLASS_BITS = {0:[1,1],1:[1,2],2:[1,3],3:[2,1],4:[2,2]}

#Number of classes
num_class  = len(CLASS_BITS.keys())

#ANGLE RANGE
ANGLES     = np.arange(min_angle,max_angle,angle_step)

outhingemask_hdf     = output+'/hingemasks.hdf'
outhinge_hdf         = output+'/hinges.hdf'
	
outhingemaskbits_hdf = output+'/hingemaskbits.hdf'
outhingebits_hdf     = output+'/hingebits.hdf'

img_counter          = 0

for i in range(len(ANGLES)):
	angle                = ANGLES[i]
	
	new_hinge = models.Hinge(boxsize=boxsize,apix=apix,angle=angle,arm1_bit=1,arm2_bit=1,invert=invert,barcode_offset=offset)
	new_hinge.make_hinge()

	hinge_img = EMNumPy.numpy2em(new_hinge.hinge)
	hinge_img.write_image(outhinge_hdf,i)

	hingemask_img = EMNumPy.numpy2em(new_hinge.hingemask)
	hingemask_img.write_image(outhingemask_hdf,i)

	for key in sorted(CLASS_BITS.keys()):
		bit1 = CLASS_BITS[key][0]
		bit2 = CLASS_BITS[key][1]

		new_hinge = models.Hinge(boxsize=boxsize,apix=apix,angle=angle,arm1_bit=bit1,arm2_bit=bit2,invert=invert,barcode_offset=offset)
		new_hinge.make_hinge()
		
		hingebits_img = EMNumPy.numpy2em(new_hinge.hingebits)
		hingebits_img.write_image(outhingebits_hdf,img_counter)

		hingemaskbits_img = EMNumPy.numpy2em(new_hinge.hingebitsmask)
		hingemaskbits_img.write_image(outhingemaskbits_hdf,img_counter)
		
		#Update image counter
		img_counter += 1

		#Flip hinge around y-axis
		new_hinge.flip_y()
		new_hinge.normalize()
		new_hinge.convert_to_float()

		hingebits_img = EMNumPy.numpy2em(new_hinge.hingebits)
		hingebits_img.write_image(outhingebits_hdf,img_counter)

		hingemaskbits_img = EMNumPy.numpy2em(new_hinge.hingebitsmask)
		hingemaskbits_img.write_image(outhingemaskbits_hdf,img_counter)

		#Update image counter
		img_counter += 1

E2end(logger)

