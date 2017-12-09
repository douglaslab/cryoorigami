#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2017-10-08 19:39:29
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import pyORIGAMImodel as models
import numpy as np
import cv2
import argparse
import matplotlib.pyplot as py
from EMAN2 import *

#Definition of the program
usage            = ['%(prog)s --length 100 --boxsize 200 --anglesep 1',
                    '------------------------------------------------------------------------',
                    'Line generator',
                    '08/23/2017',
                    'Tural Aksel',
                    '------------------------------------------------------------------------'
                    ]

#Create the parser
parser = argparse.ArgumentParser(description='',usage='\n'.join(usage))
parser.add_argument('--numlines'     ,default = 1,    type=int,   help='Number of parallel lines (Default:1')
parser.add_argument('--length'       ,default = 100,  type=int,   help='Line length (Default:100)')
parser.add_argument('--width'        ,default = 10 ,  type=int,   help='Line width (Default:5)')
parser.add_argument('--boxsize'      ,default = 140,  type=int,   help='Boxsize (Default:140')
parser.add_argument('--anglesep'     ,default = 1,    type=float, help='Angle separation (Default:1)')
parser.add_argument('--output'       ,default = '.',  type=str,   help='Output folder (Default:.)')
parser.add_argument('--linesep'      ,default = 17,   type=int,   help='Seperation between lines (Default:15)')
parser.add_argument('--apix'         ,default = 2.448,type=float, help='Pixel size in Angstrom (A)')
parser.add_argument('--singleref'    ,action  = 'store_true', help='Write only single model',default=False)
parser.add_argument("--ppid"         , type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

logger=E2init(sys.argv, -1)

args = parser.parse_args()
args_dict = vars(args)
parser.print_help()

#Get the argument
num_lines    = args_dict['numlines']
line_length  = args_dict['length']
line_width   = args_dict['width']
line_sep     = args_dict['linesep']
boxsize      = args_dict['boxsize']
anglesep     = args_dict['anglesep']
output       = args_dict['output']
apix         = args_dict['apix']
singleref    = args_dict['singleref']

if not os.path.isdir(output):
     print "Output directory doesn't exist"
     sys.exit(1)

#Create bitted hinges
lines_out_hdf         = output+"/lines-n%d-w%d-l%d-b%d-s%d.hdf"%(num_lines,line_width,line_length,boxsize,line_sep)
masks_out_hdf         = output+"/masks-n%d-w%d-l%d-b%d-s%d.hdf"%(num_lines,line_width,line_length,boxsize,line_sep)

#Create line object
new_line = models.Line(numlines=num_lines, length=line_length, width = line_width, boxsize = boxsize, angle_sep = anglesep, line_sep=line_sep)
new_line.create_line_series()
num_models = len(new_line.rotlines_list)

#Check if there is only single model to write
if singleref:
     num_models = 1

for i in range(num_models):
     angle,rotline = new_line.rotlines_list[i]
     rotline_em = EMNumPy.numpy2em(rotline)
     
     #Assign pixel size
     rotline_em['apix_x'] = apix
     rotline_em['apix_y'] = apix
     rotline_em['apix_z'] = apix
     rotline_em.write_image(lines_out_hdf,i)

     angle,rotmask = new_line.rotmask_list[i]
     rotmask_em = EMNumPy.numpy2em(rotmask)
     
     #Assign pixel size
     rotmask_em['apix_x'] = apix
     rotmask_em['apix_y'] = apix
     rotmask_em['apix_z'] = apix
     rotmask_em.write_image(masks_out_hdf,i)

E2end(logger)


