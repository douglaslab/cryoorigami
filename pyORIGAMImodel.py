#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2017-01-04 12:21:30
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import numpy as np
import matplotlib.pyplot as py
import scipy.ndimage

from EMAN2 import *

#DEFAULT HINGE PARAMETERS
ARM_LENGTH  = 186
ARM_STRUTS  = 4

#BIT PARAMETERS
BIT_STRUTS = 2


#BARCODE PARAMETERS

AVERAGE_BIT    = {'originx': 0       , 'length': 186 ,'width':   BIT_STRUTS}

LONG_BIT1      = {'originx':93+42-51 , 'length': 102 ,'width':   BIT_STRUTS}
LONG_BIT2      = {'originx':93-42-51 , 'length': 102 ,'width':   BIT_STRUTS}

SHORT_BIT1     = {'originx':93+60-33 , 'length': 66  ,'width':   BIT_STRUTS}
SHORT_BIT2     = {'originx':93   -33 , 'length': 66  ,'width':   BIT_STRUTS}
SHORT_BIT3     = {'originx':93-60-33 , 'length': 66  ,'width':   BIT_STRUTS}

def normalize_intensity(img):
    '''
    Normalize img intensity
    '''
    return (img - np.min(img))/(np.max(img)-np.min(img))

class Line:
    def __init__(self, numlines=1, length=100, width = 10, boxsize = 140, angle_sep = 1, line_sep=17):
        #Line parameters
        self.length   = length
        self.width    = width
        self.boxsize  = boxsize
        self.sigma    = width/2.0
        self.numlines = numlines
        self.linesep  = line_sep

        #The box centers
        self.boxcenterx = int(0.5*self.boxsize)
        self.boxcentery = int(0.5*self.boxsize)
        self.boxcenter  = np.array([self.boxcenterx,self.boxcentery])
        
        #The line centers
        self.linecenterx  = int(0.5*self.boxsize)
        self.linecentery  = int(0.5*self.boxsize)

        #Line left origin
        self.originx  = self.linecenterx - int(0.5*self.length)
        self.originy  = self.linecentery - int(0.5*self.width)
        self.origin   = np.array([self.originy,self.originx])

        #Right upper corner
        self.rcornery = self.linecentery + int(0.5*self.width)
        self.rcornerx = self.linecenterx + int(0.5*self.length)
        
        #Angle series
        self.angles   = np.arange(0,360,step=angle_sep)
        self.lines    = []

    def shift_linecenterx(self,shift_value=0):
        '''
        Shift linecenterx
        '''
        self.linecenterx += int(shift_value)

    def set_linecentery(self,linecentery):
        '''
        Set line centery
        '''
        self.linecentery = linecentery

    def generate_line(self):
        '''
        Generate single horizontal line
        '''
        #Make the meshgrids
        x = np.arange(self.boxsize)
        y = np.arange(self.boxsize)
        xv, yv = np.meshgrid(x, y)
        
        #Make a mask for the line
        self.linemask = np.zeros((self.boxsize,self.boxsize),dtype=np.float32)
        self.linemask[self.linecentery-int(0.5*self.width):self.linecentery+int(0.5*self.width),self.linecenterx-int(0.5*self.length):self.linecenterx+int(0.5*self.length)] = 1.0
        
        self.line  = np.exp(-(yv-self.linecentery)**2/self.sigma**2)
        self.line *= self.linemask

        #Make a mask for the n-lines
        self.mask  = np.zeros((self.boxsize,self.boxsize),dtype=np.float32)
        self.mask[self.linecentery-int(0.5*self.linesep):self.linecentery+int(0.5*self.width),self.linecenterx-int(0.5*self.length):self.linecenterx+int(0.5*self.length)] = 1.0
    
    def flip_y(self):
        '''
        Flip lines around y axis
        '''
        self.lines     = self.lines[::-1,:]
        self.linemasks = self.linemasks[::-1,:]
        self.mask      = self.mask[::-1,:]

    def flip_x(self):
        '''
        Flip lines around x axis
        '''
        self.lines     = self.lines[:,::-1]
        self.linemasks = self.linemasks[:,::-1]
        self.mask      = self.mask[:,::-1]

    def shift_lines(self, shift_array):
        '''
        Shift lines by the shift array
        '''
        self.lines      = scipy.ndimage.interpolation.shift(self.lines,     shift_array)
        self.linemasks  = scipy.ndimage.interpolation.shift(self.linemasks, shift_array)
        self.mask       = scipy.ndimage.interpolation.shift(self.mask     , shift_array)

    def shift_lines_to_left_origin(self):
        '''
        Shift lines to left origin
        '''

        #Determine the shift array
        self.shift_array = self.boxcenter - self.origin
        
        #Shift lines by the shift array
        self.shift_lines(self.shift_array)

    def generate_nlines(self):
        '''
        Generate n-parallel lines
        '''
        self.lines     = np.zeros((self.boxsize,self.boxsize),dtype=np.float32)
        self.linemasks = np.zeros((self.boxsize,self.boxsize),dtype=np.float32)
        
        for i in range(self.numlines):
            self.generate_line()
            self.lines       += self.line
            self.linemasks   += self.linemask
            self.linecentery += self.linesep

        #Update the mask parameters
        self.rcornery += (self.numlines-1)*self.linesep
        
        #Make the mask
        self.mask = np.zeros((self.boxsize,self.boxsize),dtype=np.float32)
        self.mask[self.originy:self.rcornery,self.originx:self.rcornerx] = 1

    def rotate_lines(self,angle):
        '''
        Rotate the line
        '''
        self.rotlines     = scipy.ndimage.rotate(self.lines,-angle,reshape=False)
        self.rotlines     = self.rotlines/np.max(self.rotlines)

        self.rotlinemasks = scipy.ndimage.rotate(self.linemasks,-angle,reshape=False)
        self.rotlinemasks = self.rotlinemasks/np.max(self.rotlinemasks)

        self.rotmask      = scipy.ndimage.rotate(self.mask,-angle,reshape=False)
        self.rotmask      = self.rotmask/np.max(self.rotmask)

    def normalize_lines(self):
        '''
        Normalize line intensity
        '''
        
        self.lines     = normalize_intensity(self.lines)
        self.linemasks = normalize_intensity(self.linemasks)
        self.mask      = normalize_intensity(self.mask)

    def normalize_rots(self):
        '''
        Normalize the intensity
        '''

        self.rotlines     = normalize_intensity(self.rotlines)
        self.rotlinemasks = normalize_intensity(self.rotlinemasks)
        self.rotmask      = normalize_intensity(self.rotmask)

    def create_line_series(self,shift_value=0):
        '''
        Create line series 
        '''
        #Shift linecenterx
        self.shift_linecenterx(shift_value)
        
        #Generate the horizontal line
        self.generate_nlines()

        #Prepare all the rotation lines
        self.rotlines_list     = []
        self.rotlinemasks_list = []
        self.rotmask_list      = []

        for angle in self.angles:
            #Rotate lines
            self.rotate_lines(angle)
            
            self.rotlines_list.append([angle,self.rotlines])
            self.rotlinemasks_list.append([angle,self.rotlinemasks])
            self.rotmask_list.append([angle,self.rotmask])

class Origami:
    def __init__(self, cadnano_file = None):
        self.xyz = []   #3D binary representation of an origami structure
        self.xy  = []   #2D binary representation of an origami structure

class Hinge(Origami):
    def __init__(self, boxsize = 360, apix = 2.448, angle = 37, arm1_length = ARM_LENGTH, arm2_length = ARM_LENGTH, arm1_struts = ARM_STRUTS, arm2_struts = ARM_STRUTS, arm1_bit = 1, arm2_bit = 1, flip=False, invert=False, barcode_offset=0):
        
        #Scale and canvas parameters
        self.apix        = apix
        self.scale       = 1.0/apix
        self.canvas      = int(boxsize)
        
        self.angle       = angle
        self.arm1_length = arm1_length
        self.arm2_length = arm2_length

        self.arm1_width  = arm1_struts
        self.arm2_width  = arm2_struts
        self.offset      = 0

        self.arm1        = []
        self.arm2        = []
        self.mask        = []
        self.hinge       = []

        #Process images
        self.invert        = invert

        #Bits parameters
        self.arm1_bit_num  = arm1_bit
        self.arm2_bit_num  = arm2_bit

        self.arm1_bits     = {1:LONG_BIT1 ,2:LONG_BIT2} 
        self.arm2_bits     = {1:SHORT_BIT1,2:SHORT_BIT2,3:SHORT_BIT3}

        #Arm bits
        self.arm1_bit      = self.arm1_bits[arm1_bit]
        self.arm2_bit      = self.arm2_bits[arm2_bit]

        #Tether location
        self.tether_pos    = []
        self.tether_origin = []

        #Canvas parameters
        self.double_canvas = int(2*self.canvas)
        self.half_canvas   = int(0.5*self.canvas)

        self.cx            = self.canvas
        self.cy            = self.canvas

        #Barcode offset
        self.barcode_offset= barcode_offset
        self.length_offset = barcode_offset

        #Flip parameter
        self.flip        = flip
        self.flipparam   = 1
        if self.flip : self.flipparam   = -1 

    def make_hinge(self):
        #Create 2xcanvas-size
        self.double_canvas = int(2*self.canvas)
        self.half_canvas   = int(0.5*self.canvas)
        
        #Canvas centers
        self.cx            = self.canvas
        self.cy            = self.canvas

        #Create the line objects - Adjust the arm lengths based on the arm length offset
        self.arm1_lines     = Line(numlines=self.arm1_width, length=self.arm1_length+self.length_offset, boxsize = self.double_canvas, width = 10, line_sep=17)
        self.arm2_lines     = Line(numlines=self.arm2_width, length=self.arm2_length+self.length_offset, boxsize = self.double_canvas, width = 10, line_sep=17)

        self.arm1bits_lines = Line(numlines=self.arm1_bit['width'], length=self.arm1_bit['length'], boxsize = self.double_canvas, width = 10, line_sep=17)
        self.arm2bits_lines = Line(numlines=self.arm2_bit['width'], length=self.arm2_bit['length'], boxsize = self.double_canvas, width = 10, line_sep=17)

        #Draw the line objects
        self.arm1_lines.generate_nlines()
        self.arm2_lines.generate_nlines()

        self.arm1bits_lines.generate_nlines()
        self.arm2bits_lines.generate_nlines()
        
        #Shift the lines to origins
        self.arm1_lines.shift_lines_to_left_origin()
        self.arm2_lines.shift_lines_to_left_origin()

        self.arm1bits_lines.shift_lines_to_left_origin()
        self.arm2bits_lines.shift_lines_to_left_origin()

        #Shift the bits to register with the arms - adjust te barcode origin based on the barcode offset
        arm1bit_shift_array = np.array([self.arm1_lines.linecentery-self.arm1_lines.boxcentery,self.arm1_bit['originx']+self.barcode_offset])
        arm2bit_shift_array = np.array([self.arm2_lines.linecentery-self.arm2_lines.boxcentery,self.arm2_bit['originx']+self.barcode_offset])

        self.arm1bits_lines.shift_lines(arm1bit_shift_array)
        self.arm2bits_lines.shift_lines(arm2bit_shift_array)

        #Flip arm 2 components
        self.arm2_lines.flip_y()
        self.arm2bits_lines.flip_y()

        #Rotate arm 1 and arm2 in opposite angles
        arm1_angle =  0.5*self.angle
        arm2_angle = -0.5*self.angle

        self.arm1_lines.rotate_lines(arm1_angle)
        self.arm1bits_lines.rotate_lines(arm1_angle)

        self.arm2_lines.rotate_lines(arm2_angle)
        self.arm2bits_lines.rotate_lines(arm2_angle)

        #Create two canvas for the model
        self.arm1          = self.arm1_lines.rotlines
        self.arm2          = self.arm2_lines.rotlines
        self.hinge         = self.arm1 + self.arm2

        #Create the canvas for mask
        self.arm1mask      = self.arm1_lines.rotmask
        self.arm2mask      = self.arm2_lines.rotmask
        self.hingemask     = self.arm1mask + self.arm2mask

        #Create two canvas for the model
        self.arm1bits      = self.arm1 + self.arm1bits_lines.rotlines
        self.arm2bits      = self.arm2 + self.arm2bits_lines.rotlines
        self.hingebits     = self.arm1bits + self.arm2bits

        #Create the canvas for mask
        self.arm1bitsmask  = self.arm1mask + self.arm1bits_lines.rotmask
        self.arm2bitsmask  = self.arm2mask + self.arm2bits_lines.rotmask
        self.hingebitsmask = self.arm1bitsmask + self.arm2bitsmask

        #Shift the models
        self.shift_to_cmass()

        #Flip around x-axis
        self.flip_x()

        #Extract hinge
        self.extract_hinge()

        #Convert to 16 bit
        self.convert_to_float()

        #Normalize
        self.normalize()

        #Invert
        if self.invert: self.invert_images()

    def normalize(self):
        '''
        Normalize the data
        '''
        self.hinge     = normalize_intensity(self.hinge)
        self.hingemask = normalize_intensity(self.hingemask)

        self.hingebits     = normalize_intensity(self.hingebits)
        self.hingebitsmask = normalize_intensity(self.hingebitsmask)

    def shift_to_cmass(self):
        '''
        Shift to cmass
        '''

        #Determine the center of mass in x and y dimesion
        self.cmassy_1 = 0.5*self.arm1_length*np.sin( 0.5*self.angle*np.pi/180.0)
        self.cmassy_2 = 0.5*self.arm2_length*np.sin(-0.5*self.angle*np.pi/180.0)
        
        self.cmassx_1 = 0.5*self.arm1_length*np.cos( 0.5*self.angle*np.pi/180.0)
        self.cmassx_2 = 0.5*self.arm2_length*np.cos(-0.5*self.angle*np.pi/180.0)

        self.cmassy   = int(0.5*(self.cmassy_1+self.cmassy_2))
        self.cmassx   = int(0.5*(self.cmassx_1+self.cmassx_2))

        #Determine the shift array
        shift_array   = np.array([-self.cmassy,-self.cmassx])

        #Shift the images
        self.hinge        = scipy.ndimage.interpolation.shift(self.hinge    , shift_array)
        self.hingemask    = scipy.ndimage.interpolation.shift(self.hingemask, shift_array)

        self.hingebits    = scipy.ndimage.interpolation.shift(self.hingebits    , shift_array)
        self.hingebitsmask= scipy.ndimage.interpolation.shift(self.hingebitsmask, shift_array)

    def invert_images(self):
        '''
        Invert images
        '''
        self.hinge         = 1.0 - self.hinge
        self.hingemask     = 1.0 - self.hingemask

        self.hingebits     = 1.0 - self.hingebits
        self.hingebitsmask = 1.0 - self.hingebitsmask 

    def flip_x(self):
        '''
        Flip around x-axis
        '''
        self.hinge         = self.hinge[:,::-1] 
        self.hingemask     = self.hingemask[:,::-1]

        self.hingebits     = self.hingebits[:,::-1]
        self.hingebitsmask = self.hingebitsmask[:,::-1]

    def flip_y(self):
        '''
        Flip around x-axis
        '''
        self.hinge         = self.hinge[::-1,:] 
        self.hingemask     = self.hingemask[::-1,:]

        self.hingebits     = self.hingebits[::-1,:]
        self.hingebitsmask = self.hingebitsmask[::-1,:]

    def extract_hinge(self):
        '''
        Extract the hinge
        '''
        self.hinge         = self.hinge[self.cy-self.half_canvas:self.cy+self.half_canvas,self.cx-self.half_canvas:self.cx+self.half_canvas]
        self.hingemask     = self.hingemask[self.cy-self.half_canvas:self.cy+self.half_canvas,self.cx-self.half_canvas:self.cx+self.half_canvas]

        self.hingebits     = self.hingebits[self.cy-self.half_canvas:self.cy+self.half_canvas,self.cx-self.half_canvas:self.cx+self.half_canvas]
        self.hingebitsmask = self.hingebitsmask[self.cy-self.half_canvas:self.cy+self.half_canvas,self.cx-self.half_canvas:self.cx+self.half_canvas]

    def convert_to_float(self):
        '''
        Convert 16 bit
        '''
        self.hinge     = np.array(self.hinge,dtype=np.float32)
        self.hingemask = np.array(self.hingemask,dtype=np.float32)

        self.hingebits     = np.array(self.hingebits,dtype=np.float32)
        self.hingebitsmask = np.array(self.hingebitsmask,dtype=np.float32)

if __name__ == '__main__':

    #Example run
    new_hinge = Hinge(flip=False)
    new_hinge.make_hinge()
    py.figure(0)
    py.subplot(2,2,1)
    py.imshow(new_hinge.hinge,cmap='gray')
    py.subplot(2,2,2)
    py.imshow(new_hinge.hingemask,cmap='gray')
    py.subplot(2,2,3)
    py.imshow(new_hinge.hingebits,cmap='gray')
    py.subplot(2,2,4)
    py.imshow(new_hinge.hingebitsmask,cmap='gray')

    py.show()