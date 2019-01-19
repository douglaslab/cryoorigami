#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 13:35:55
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import yaml
import numpy as np
import matplotlib.pyplot as py
import scipy.ndimage

from pyfftw.builders import fft2
from pyfftw.builders import ifft2


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

    # Get the functions and its parameters
    parameter_lists = [group.split(':') for group in parameter_input]

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


def normalize_minmax(data):

    return (data - np.min(data)) / (np.max(data) - np.min(data))


def rad2deg(rad):
    '''
    Convert radian to eular degree
    '''
    return 180.0*rad/np.pi


def euler2rot2D(alpha):
    '''
    Eular angle to 2D rotation matrix
    '''
    ca = np.cos(alpha*np.pi/180)
    sa = np.sin(alpha*np.pi/180)

    r = np.array([[ca, -sa],
                  [+sa, ca]])
    return r


def get_distances(single_coordinate, list_coordinate):
    '''
    Get distances between a single coordinate and a coordinate list
    '''
    distances = list_coordinate - single_coordinate

    return np.sqrt(np.sum(distances**2, axis=1))

### FFT utilities ###


def eval_fft2(data, num_threads=12):
    '''
    Evaluate rfft2
    '''
    ft =  fft2(data,
               threads=num_threads,
               planner_effort="FFTW_ESTIMATE",
               overwrite_input=False,
               auto_align_input=True,
               auto_contiguous=True)

    # Return rfft2 result
    return ft(data.copy(), np.zeros(ft.output_shape, dtype=ft.output_dtype)).copy()


def eval_ifft2(data, num_threads=12):
    '''
    Evaluate rfft2
    '''
    ift =  ifft2(data,
                 threads=num_threads,
                 planner_effort="FFTW_ESTIMATE",
                 auto_align_input=True,
                 auto_contiguous=True)

    # Return rfft2 result
    return ift(data.copy(), np.zeros(ift.output_shape, dtype=ift.output_dtype)).copy()


### PLOT Utilities ###

def fft_plot(fft_output):
    '''
    Plot fft output
    '''
    py.imshow(np.fft.fftshift(fft_output), cmap='gray')
    py.show()


def close_plot():
    py.close()


### MASK utilities ###

def circular_mask(shape, center=None, radius=None):

    # Determine the radius
    radius = int(radius)

    # use the middle of the image
    if center is None:
        center = [shape[1]//2, shape[0]//2]

    # use the smallest distance between the center and image walls
    if radius is None:
        radius = min(center[0], center[1], shape[1]-center[0], shape[0]-center[1])

    Y, X = np.ogrid[:shape[0], :shape[1]]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return np.array(mask, dtype='float32')


def threshold_mask(data, threshold_val=0.05, erode=2, dilate=5):
    '''
    Create a threshold based mask from a 2D img data
    '''
    mask = np.array(data > threshold_val, dtype=int)

    # Perform erosion to remove small disconnected areas
    if erode > 0:
        mask = scipy.ndimage.binary_erosion(mask, structure=np.ones((erode, erode))).astype(mask.dtype)

    # Perform erosion to remove small disconnected areas
    if dilate > 0:
        mask = scipy.ndimage.binary_dilation(mask, structure=np.ones((dilate, dilate))).astype(mask.dtype)

    return np.array(mask, dtype='float32')
