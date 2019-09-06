#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 13:35:55
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as py
import scipy.ndimage

from pyfftw.builders import fft2
from pyfftw.builders import ifft2


def expmap(e):
    '''
    Convert axis-angle vector into 3D rotation matrix
    '''
    theta = np.linalg.norm(e)
    if theta < 1e-16:
        return np.identity(3, dtype=e.dtype)
    w = e / theta
    k = np.array([[0, w[2], -w[1]],
                  [-w[2], 0, w[0]],
                  [w[1], -w[0], 0]], dtype=e.dtype)
    r = np.identity(3, dtype=e.dtype) + np.sin(theta) * k + (1 - np.cos(theta)) * np.dot(k, k)

    return r

def rodriguez2euler(e):
    '''
    Convert rodriguez rotation convention to euler
    '''
    R = expmap(e)
    return rot2euler(R)

def rot2euler(r):
    '''
    Decompose rotation matrix into Euler angles
    '''
    # assert(isrotation(r))
    # Shoemake rotation matrix decomposition algorithm with same conventions as Relion.
    epsilon = np.finfo(np.double).eps
    abs_sb = np.sqrt(r[0, 2] ** 2 + r[1, 2] ** 2)
    if abs_sb > 16 * epsilon:
        gamma = np.arctan2(r[1, 2], -r[0, 2])
        alpha = np.arctan2(r[2, 1], r[2, 0])
        if np.abs(np.sin(gamma)) < epsilon:
            sign_sb = np.sign(-r[0, 2]) / np.cos(gamma)
        else:
            sign_sb = np.sign(r[1, 2]) if np.sin(gamma) > 0 else -np.sign(r[1, 2])
        beta = np.arctan2(sign_sb * abs_sb, r[2, 2])
    else:
        if np.sign(r[2, 2]) > 0:
            alpha = 0
            beta = 0
            gamma = np.arctan2(-r[1, 0], r[0, 0])
        else:
            alpha = 0
            beta = np.pi
            gamma = np.arctan2(r[1, 0], -r[0, 0])
    return alpha, beta, gamma

def get_min_angle_diff(angle_diff):
    '''
    Get minimum of the angle difference between two vectors
    '''
    data = angle_diff%360
    dataComp = 360 - data
    dataMin  = pd.DataFrame([data, dataComp]).min()

    return dataMin

def estimate_cistem_params(fsc50, prev_percent, prev_limit, num_particles, mask_radius, num_classes):
    '''
    Estimate cistem auto-refine parameters based on Cistem paper
    '''
    if fsc50 is not None and mask_radius is not None:
        new_limit = 1.0/(1.0/fsc50 - 1.0/mask_radius)
    else:
        new_limit = None

    if prev_limit is not None and num_particles is not None and prev_percent is not None:
        # Calculate percent-particles
        calc_percent  = 8000.0*num_classes*np.exp(75/prev_limit**2)/num_particles

        # Better-resolution percent
        new_percent_betterRes = max([prev_percent, calc_percent])

        # Worse-resolution percent
        new_percent_worseRes  = 1.5*prev_percent 
        if new_percent_worseRes > 1:
            new_percent_worseRes = 1
    else:
        new_percent_betterRes = 1.0
        new_percent_worseRes  = 1.0

    return new_limit, new_percent_betterRes, new_percent_worseRes


def write_config_file(fname, args_dict):
    '''
    Dump parameters into ymal config file
    '''
    with open(fname, 'w') as outfile:
        yaml.dump(args_dict, outfile, default_flow_style=False)


def convert_dict2str(dictionary={}):
    '''
    Convert dictionary to string
    '''
    conversion_str = ''
    if dictionary:
        conversion_str = ','.join(sorted([str(k)+':'+str(v) for k, v in dictionary.items()]))

    return conversion_str


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

def circular_mask(shape, center=None, radius=None, soft_edge=None):

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

    mask = np.array(dist_from_center <= radius, dtype=np.float32)

    # Check for the soft edge width
    if soft_edge is not None:
        mask = scipy.ndimage.filters.gaussian_filter(mask, soft_edge)

    return mask


def threshold_mask(data, threshold_high=0.05, threshold_low=None, erode=2, dilate=5):
    '''
    Create a threshold based mask from a 2D img data
    '''
    mask = np.array(data > threshold_high, dtype=int)

    # Apply lower threshold
    if threshold_low is not None:
        mask += np.array(data < threshold_low, dtype=int)

    # Perform erosion to remove small disconnected areas
    if erode > 0:
        mask = scipy.ndimage.binary_erosion(mask, structure=np.ones((erode, erode))).astype(mask.dtype)

    # Perform erosion to remove small disconnected areas
    if dilate > 0:
        mask = scipy.ndimage.binary_dilation(mask, structure=np.ones((dilate, dilate))).astype(mask.dtype)

    return np.array(mask, dtype='float32')
