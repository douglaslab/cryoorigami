#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 13:35:55
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import yaml
import numpy as np
import matplotlib.pyplot as py

from pyfftw.builders import rfft2
from pyfftw.builders import irfft2
from numba import jit, prange


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


def ctf_freqs(shape, apix=1.0):
    """
    :param shape: Shape tuple.
    :param apix: Pixel size in Angstrom.
    """
    assert shape[0] == shape[1]
    xfreq   = np.fft.rfftfreq(shape[1], apix)
    yfreq   = np.fft.rfftfreq(shape[0], apix)

    sx, sy  = np.meshgrid(xfreq, yfreq)
    s       = np.sqrt(sx**2 + sy**2)

    # Determine angle grids
    a = np.arctan2(sy, sx)

    # Determine r
    spacing = 1.0/(shape[0]*apix)
    r  = np.round(s/spacing)*spacing

    return s, sx, sy, a, r


@jit(nopython=True, parallel=True)
def frc(ft1, ft2, r, rbins):
    '''
    Compute frc between two ft
    '''
    out = np.zeros(ft1.shape, dtype=np.float32)
    for i in prange(len(rbins)):
        mask = r == rbins[i]
        corr  = np.sum(ft1[mask]*np.conj(ft2[mask]))
        norm1 = np.sqrt(np.sum(ft1[mask]**2))
        norm2 = np.sqrt(np.sum(ft2[mask]**2))
        out[mask] = corr/(norm1*norm2)

    return out


def eval_pshift(sx, sy, originx, originy):
    """
    Evaluate phase shift

    :param sx: Precomputed frequency-x grid for CTF evaluation.
    :param sy: Precomputed frequency-y grid for CTF evaluation.
    """

    return np.exp(-2 * np.pi * 1j * (-originx * sx + -originy * sy))


def eval_ctf(s, a, def1, def2, angast=0, phase=0, kv=300, ac=0.1, cs=2.0, bf=0, lp=0):
    """
    :param s: Precomputed frequency grid for CTF evaluation.
    :param a: Precomputed frequency grid angles.
    :param def1: 1st prinicipal underfocus distance (Å).
    :param def2: 2nd principal underfocus distance (Å).
    :param angast: Angle of astigmatism (deg) from x-axis to azimuth.
    :param phase: Phase shift (deg).
    :param kv:  Microscope acceleration potential (kV).
    :param ac:  Amplitude contrast in [0, 1.0].
    :param cs:  Spherical aberration (mm).
    :param bf:  B-factor, divided by 4 in exponential, lowpass positive.
    :param lp:  Hard low-pass filter (Å), should usually be Nyquist.
    """
    angast = np.deg2rad(angast)
    kv = kv * 1e3
    cs = cs * 1e7
    lamb = 12.2643247 / np.sqrt(kv * (1. + kv * 0.978466e-6))
    def_avg = -(def1 + def2) * 0.5
    def_dev = -(def1 - def2) * 0.5

    # k paramaters
    k1 = np.pi / 2. * 2 * lamb
    k2 = np.pi / 2. * cs * lamb**3
    k3 = np.sqrt(1 - ac**2)
    k4 = bf / 4.            # B-factor, follows RELION convention.
    k5 = np.deg2rad(phase)  # Phase shift.

    # Hard low-pass filter
    if lp != 0:
        s *= s <= (1. / lp)
    s2 = s**2
    s4 = s2**2
    dZ = def_avg + def_dev * (np.cos(2 * (a - angast)))
    gamma = (k1 * dZ * s2) + (k2 * s4) - k5
    ctf = -(k3 * np.sin(gamma) - ac*np.cos(gamma))

    # Enforce envelope
    if bf != 0:
        ctf *= np.exp(-k4 * s2)
    return ctf


def eval_rfft2(data, num_threads=12):
    '''
    Evaluate rfft2
    '''
    ft = rfft2(data,
               threads=num_threads,
               planner_effort="FFTW_ESTIMATE",
               overwrite_input=False,
               auto_align_input=True,
               auto_contiguous=True)

    # Return rfft2 result
    return ft(data.copy(), np.zeros(ft.output_shape, dtype=ft.output_dtype)).copy()


def eval_irfft2(data, num_threads=12):
    '''
    Evaluate rfft2
    '''
    ift = irfft2(data,
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

def circular_mask(shape, center=None, diameter=None):

    # Determine the radius
    radius = int(diameter//2)

    # use the middle of the image
    if center is None:
        center = [shape[1]//2, shape[0]//2]

    # use the smallest distance between the center and image walls
    if radius is None:
        radius = min(center[0], center[1], shape[1]-center[0], shape[0]-center[1])

    Y, X = np.ogrid[:shape[0], :shape[1]]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask.view('float32')


def binary_mask(data, threshold_val=0.05):
    '''
    Create a binary mask from a 2D img data
    '''
    mask = data > threshold_val
    return mask.view('float32')
