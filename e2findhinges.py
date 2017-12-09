#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2017-08-24 12:54:50
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

from EMAN2 import *


import pyORIGAMIem
import numpy as np
import scipy.ndimage
import scipy.signal
import matplotlib.pyplot as py
import os
import sys

progname = os.path.basename(sys.argv[0])
usage = """prog <input> [options]

Detect the hinges from arm particle data 
"""
    
parser = EMArgumentParser(usage=usage,version=EMANVERSION)

parser.add_argument("--particles"       ,type=str  , help="Arm particle file. Required", default=None)
parser.add_argument("--armref"          ,type=str  , help="Reference file for arm center detection. Required", default=None)
parser.add_argument("--simmx"           ,type=str  , help="Simmx file. Required", default=None)
parser.add_argument("--micrograph"      ,type=str  , help="Micrograph file for which the arms are detected. Required", default=None)
parser.add_argument("--linerefs"        ,type=str  , help="Directory for the line references. Required", default=None)
parser.add_argument("--modelrefs"       ,type=str  , help="Dricetory for the model references. Required", default=None)
parser.add_argument("--apix"            ,type=float, help="Pixel size of the particles and micrographs", default=2.44)
parser.add_argument("--radius"          ,type=float, help="Particle radius in Angstrom", default=380)
parser.add_argument("--refangle"        ,type=float, help="Reference angle", default=-26.08)
parser.add_argument("--plot"            ,action  = 'store_true',help='Plot the results on the micrograph',default=False)
parser.add_argument("--saveplot"        ,action  = 'store_true',help='Save the micrograph plots',default=False)
parser.add_argument("--singleref"       ,action  = 'store_true',help='There is only a single class to compare with.',default=False)
parser.add_argument("--keepedges"       ,action  = 'store_true',help='Keep the particles close to edges.',default=False)
parser.add_argument("--writeptcls"      ,action  = 'store_true',help='Write particles.',default=False)
parser.add_argument("--output"          ,type=str  , help="Particles output directory. Required", default="particles")

(options, args) = parser.parse_args()

#Read the parameters

#Similarity file
simmx_file      = options.simmx

#Micrograph file
micrograph_file = options.micrograph

#Particle file
particle_file   = options.particles

#Reference angle for simmx file
refangle        = options.refangle

#Plot
plot_results    = options.plot

#Pixel size
apix            = options.apix

#Keep edges
keepedges       = options.keepedges

#Single ref
singleref       = options.singleref

#Line references directory
line_refs_dir   = options.linerefs

#Model references directory
model_refs_dir  = options.modelrefs

#Write particles
write_particles = options.writeptcls

#Save micrograph plots
save_plot       = options.saveplot

#Particles output
particles_output= options.output


#Check if output directory exists
if not os.path.isdir(particles_output):
	os.mkdir(particles_output)

#Check the parameters
if simmx_file == None or not os.path.isfile(simmx_file):
    sys.exit('Simmilarity comparison file is missing.')

if micrograph_file == None or not os.path.isfile(micrograph_file):
    sys.exit('Micrograph is missing.')

if particle_file == None or not os.path.isfile(particle_file):
    sys.exit('Particles file is missing.')

#Lines references
linerefs_file          = line_refs_dir+'/lines-n1-w10-l60-b140-s17.hdf'

#Small arms references
smallarmrefs_file      = line_refs_dir+'/lines-n4-w10-l90-b280-s17.hdf'

#Long arms references
longarmrefs_file       = line_refs_dir+'/lines-n4-w10-l186-b280-s17.hdf'

#Short barcode files
shortbarcode_ref_file  = line_refs_dir+'/masks-n2-w10-l66-b140-s17.hdf'

#Long barcode file
longbarcode_ref_file   = line_refs_dir+'/masks-n2-w10-l102-b140-s17.hdf'

#Hinge angle models
hinge_refs_file        = model_refs_dir+'/hinges.hdf'

#Hinge angle models
hingebits_refs_file    = model_refs_dir+'/hingebits.hdf'

#Read the micrograph
new_micrograph = pyORIGAMIem.EMMicrograph(particle_output_dir = particles_output)

#Read line refs
new_micrograph.read_line_refs(linerefs_file)

#Read small arm refs
new_micrograph.read_small_arm_refs(smallarmrefs_file)

#Read long arm refs
new_micrograph.read_long_arm_refs(longarmrefs_file)

#Read short barcode references
new_micrograph.read_short_barcode_refs(shortbarcode_ref_file)

#Read long barcode references
new_micrograph.read_long_barcode_refs(longbarcode_ref_file)

#Read hinge references
new_micrograph.read_hinge_refs(hinge_refs_file)

#Read hinge bits references
new_micrograph.read_hingebits_refs(hingebits_refs_file)

#Create micrograph
new_micrograph.read_micrograph(micrograph_file)

#Normalize micrograph
new_micrograph.normalize_micrograph()

#Make a flat version of the micrograph
new_micrograph.make_flat()

#Calculate micrograph mean intensity
new_micrograph.calculate_mean()

#Read simmx output
new_micrograph.read_e2simmx_output(simmx_file)

#Classify particles based on e2simmx score
new_micrograph.classify_e2simmx(skip_zero=singleref)

#Set reference angles
new_micrograph.set_reference_angles({0:refangle,1:refangle})

#Read e2simmx particles
new_micrograph.read_e2simmx_particles(particle_file)

#Make arm particles
new_micrograph.make_arm_particles()

#Plot results
if plot_results:
	new_micrograph.plot_micrograph()
	
	#Save raw micrograph
	if save_plot: new_micrograph.save_plot(fname='micrograph_raw.png')

#Detect hinges
new_micrograph.detect_hinges()

#Remove hinge conflicts
new_micrograph.remove_hinge_conflicts()

#Count the number of valid hinges
new_micrograph.count_hinges()

#Write particles
if write_particles:
	
	new_micrograph.extract_hinges()

#Plot results
if plot_results:
	#Plot text
	new_micrograph.PLOT_TEXT=True

	#Plot particle locations
	new_micrograph.plot_e2simmx_particles()

	#Save the arms only figure
	if save_plot: new_micrograph.save_plot(fname='micrograph_arms.png')

	#Plot hinges
	new_micrograph.plot_hinge_particles()

	#Plot mask boundaries
	if write_particles:
		
		new_micrograph.plot_hinge_masks()

	#Save the hinge figure
	if save_plot: new_micrograph.save_plot(fname='micrograph_hinges.png')

	#Show plot
	new_micrograph.show_plot()

