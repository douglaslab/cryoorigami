#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2019-01-21 08:17:52
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import mrcfile
import numpy as np
import math
import scipy.ndimage
import random
from numba import jit


def fibonacci_sphere(num_points=1, randomize=True):
    '''
    Get points from a fibonacci sphere
    '''

    random_num = 1.0
    if randomize:
        random_num = random.random()*num_points

    # The offset for sampling
    offset = 2./num_points

    # The increment
    increment = math.pi*(3.-math.sqrt(5.))

    # Points list
    points_list = np.arange(num_points)

    # Y coordinate
    points_y = points_list*offset - 1 + 0.5*offset

    # Get distance from center
    points_r = np.sqrt(1.0-points_y**2)

    # Get the phi angle
    points_phi = ((points_list+random_num) % num_points)*increment

    # X-Z coordinates
    points_x = np.cos(points_phi)*points_r
    points_z = np.sin(points_phi)*points_r

    # Stack the x,y,z coordinates
    points_coor = np.hstack((np.vstack(points_x), np.vstack(points_y), np.vstack(points_z)))

    return points_coor


def calc_ccc(current_img2D, other_img2D, mask=None):
    '''
    Measure ccc with other img2D
    '''
    current_mean, current_std = calc_mean_std_intensity(current_img2D, mask)
    other_mean, other_std     = calc_mean_std_intensity(other_img2D, mask)

    if mask is not None and np.sum(mask) > 0:
        cross_correlation = np.average((current_img2D-current_mean)*(other_img2D-other_mean), weights=mask)
    else:
        cross_correlation = np.mean((current_img2D-current_mean)*(other_img2D-other_mean))

    return cross_correlation/(current_std*other_std)


def flipX_ptcl(ptcl_star):
    '''
    Flip particle
    '''
    img2D = read_ptcl_mrc(ptcl_star)
    img2D = flipX_img2D(img2D, ptcl_star)

    return img2D


def is_flipX(ptcl_star):
    '''
    Check if particle needs to be flipped
    '''
    if 'rlnIsFlip' in ptcl_star and ptcl_star['rlnIsFlip'] == 1:
        return True
    else:
        return False


def flipX_img2D(img2D, ptcl_star):
    '''
    Flip image around X-axis
    '''
    new_img2D = img2D
    if is_flipX(ptcl_star):
        new_img2D = img2D[:, ::-1]

    return new_img2D


def generate_noise(img2D, mask=None, noise_mean=0.0, noise_std=1.0, sigma=0):
    '''
    Make Noise
    '''
    np.random.seed()
    noise     = np.random.normal(noise_mean, noise_std, img2D.shape)
    new_img2D = img2D.copy()

    try:
        if mask is not None:
            new_img2D = noise*mask + new_img2D*(1.0-mask)
        else:
            new_img2D = noise

        return np.array(new_img2D, dtype='float32')

    except IndexError:
        return None


def clip_img2D(img2D, clipbox=None, origin=[0, 0]):
    '''
    Clip img2D
    '''
    # Check for clip
    if clipbox is not None and clipbox < img2D.shape[0] and clipbox < img2D.shape[1]:

        # Get the origin
        originX, originY = origin

        # Get the img center
        centerY = img2D.shape[0] // 2 + int(originY)
        centerX = img2D.shape[1] // 2 + int(originX)

        # Determine halfbox size
        halfbox = clipbox // 2

        # Clip image
        img2D = img2D[centerY-halfbox:centerY+halfbox, centerX-halfbox:centerX+halfbox]

    return img2D


def read_ptcl_mrc(ptcl_star, transform=False, fft_mask=None, clipbox=None, background_mask=None, recenter=False):
    '''
    Particle mrc data
    '''

    # Read particle data
    particle_image_num, particle_image_name = ptcl_star['rlnImageName'].split('@')

    with mrcfile.mmap(particle_image_name, mode='r') as mrc_data:
        if len(mrc_data.data.shape) == 3:
            img2D = mrc_data.data[int(particle_image_num)-1]
        else:
            img2D = mrc_data.data

        # If transform is ON, transform the ptcl image and rest the offset distance and angles
        if transform:
            img2D = transform_ptcl_img2D(img2D, ptcl_star)

        # Check pass filter options
        if fft_mask is not None:
            img2D = fft_filter(img2D, fft_mask)

        # Get img center
        img_origin = [0, 0]

        # If recenter option is on
        if recenter:
            fracX, intX = math.modf(ptcl_star['rlnOriginX'])
            fracY, intY = math.modf(ptcl_star['rlnOriginY'])
            img_origin  = [intX, intY]

        # Check for clip
        img2D = clip_img2D(img2D, clipbox, origin=img_origin)

        # Normalize using the background mask
        if background_mask is not None:
            img2D = norm_intensity(img2D, background_mask, new_mean=0.0, new_std=1.0)

    return np.copy(img2D)


def fft_filter(img2D, fft_mask):
    '''
    Perform highpass and/or lowpass filter on the image
    '''
    fft_img2D       = np.fft.fft2(img2D)
    fft_img2D_shift = np.fft.fftshift(fft_img2D)
    fft_img2D       = np.fft.ifftshift(fft_img2D_shift*fft_mask)

    return np.real(np.fft.ifft2(fft_img2D))


def write_img2D(img2D, ptcl_mrc, ptcl_index):
    '''
    Write img2D
    '''
    ptcl_mrc.data[ptcl_index] = img2D


def fft_img2D(img2D, mask=None):
    '''
    Take fft of image
    '''
    if mask is not None and mask.shape == img2D.shape:
            img2D_fft = np.fft.rfft2(mask*img2D)
    else:
            img2D_fft = np.fft.rfft2(img2D)

    return img2D_fft


def ifft_img2D(fft2D):
    '''
    Take ifft of a fft
    '''
    return np.real(np.fft.irfft2(fft2D))


def rotate_img2D(img2D, angle):
    '''
    Rotate img2D by an eularangle
    '''
    return scipy.ndimage.rotate(img2D, angle=-angle, axes=(0, 1), reshape=False)


def shift_img2D(img2D, shiftX, shiftY):
    '''
    Shift img2D by a vector
    '''
    return scipy.ndimage.shift(img2D, shift=[shiftY, shiftX])


def inv_rotate_ptcl_img2D(img2D, ptcl_star):
    '''
    Inverse rotation based on ptcl data
    '''
    psi = ptcl_star['rlnAnglePsi']
    return rotate_img2D(img2D, -psi)


def inv_shift_ptcl_img2D(img2D, ptcl_star):
    '''
    Inverse shift based on ptcl data
    '''
    originX = ptcl_star['rlnOriginX']
    originY = ptcl_star['rlnOriginY']

    return shift_img2D(img2D, -originX, -originY)


def inv_transform_ptcl_img2D(img2D, ptcl_star):
    '''
    Inverse ptcl transform
    '''

    img2D_tr = inv_rotate_ptcl_img2D(img2D, ptcl_star)
    return inv_shift_ptcl_img2D(img2D_tr, ptcl_star)


def transform_ptcl_img2D(img2D, ptcl_star):
    img2D_tr = shift_ptcl_img2D(img2D, ptcl_star)
    return rotate_ptcl_img2D(img2D_tr, ptcl_star)


def rotate_ptcl_img2D(img2D, ptcl_star):
    '''
    Rotation based on ptcl data
    '''
    psi = ptcl_star['rlnAnglePsi']
    return rotate_img2D(img2D, psi)


def shift_ptcl_img2D(img2D, ptcl_star):
    '''
    Shift based on ptcl data
    '''
    originX = ptcl_star['rlnOriginX']
    originY = ptcl_star['rlnOriginY']

    return shift_img2D(img2D, originX, originY)


def eval_ptcl_ctf(ctf_s, ctf_a, ptcl_star, bf=0, do_intact_until_first_peak=True):
    '''
    Determine ctf from particle data
    '''
    defU       = ptcl_star['rlnDefocusU']
    defV       = ptcl_star['rlnDefocusV']
    defA       = ptcl_star['rlnDefocusAngle']
    phaseShift = ptcl_star['rlnPhaseShift']
    kV         = ptcl_star['rlnVoltage']
    ac         = ptcl_star['rlnAmplitudeContrast']
    cs         = ptcl_star['rlnSphericalAberration']

    return eval_ctf(ctf_s, ctf_a, defU, defV, defA, phaseShift, kV, ac, cs, bf, do_intact_until_first_peak)


def subtract_ctf(img2D, ctf):
    '''
    Subtract ctf from img2D
    '''
    fft2D  = fft_img2D(img2D)
    fft2D  = fft2D/ctf
    ifft2D = ifft_img2D(fft2D)
    return ifft2D


def correct_fft_ctf(fft2D, ctf):
    '''
    Correct ctf
    '''
    return fft2D*ctf


@jit(nopython=True)
def eval_ctf(ctf_s, ctf_a, defU, defV, defA=0, phaseShift=0, kv=300, ac=0.1, cs=2.0, bf=0, do_intact_until_first_peak=False):
    '''
    :param defU: 1st prinicipal underfocus distance (Å).
    :param defV: 2nd principal underfocus distance (Å).
    :param defA: Angle of astigmatism (deg) from x-axis to azimuth.
    :param phaseShift: Phase shift (deg).
    :param kv:  Microscope acceleration potential (kV).
    :param ac:  Amplitude contrast in [0, 1.0].
    :param cs:  Spherical aberration (mm).
    :param bf:  B-factor, divided by 4 in exponential, lowpass positive.
    :param lp:  Hard low-pass filter (Å), should usually be Nyquist.
    :param hp:  High-pass filter (Å)
    '''

    # parameter unit conversions
    defA = np.deg2rad(defA)
    kv = kv * 1e3
    cs = cs * 1e7
    lamb = 12.2643247 / np.sqrt(kv * (1. + kv * 0.978466e-6))
    def_avg = -(defU + defV) * 0.5
    def_dev = -(defU - defV) * 0.5

    # k paramaters
    k1 = np.pi / 2. * 2 * lamb
    k2 = np.pi / 2. * cs * lamb**3
    k3 = np.arctan(ac/np.sqrt(1 - ac**2))
    k4 = bf / 4.                # B-factor, follows RELION convention.
    k5 = np.deg2rad(phaseShift)  # Phase shift.

    # Assign s grid
    s = ctf_s

    s2 = s**2
    s4 = s2**2
    dZ = def_avg + def_dev * (np.cos(2 * (ctf_a - defA)))
    gamma = (k1 * dZ * s2) + (k2 * s4) - k5 - k3

    # Determine ctf
    img2D_ctf = -np.sin(gamma)

    # Do intact until first peak
    if do_intact_until_first_peak:

        # Loop for the numba jit
        for i in range(gamma.shape[0]):
            for j in range(gamma.shape[1]):
                if np.abs(gamma[i, j]) < np.pi/2:
                    img2D_ctf[i, j] = 1.0

    # Enforce envelope
    if bf != 0:
        img2D_ctf *= np.exp(-k4 * s2)

    return img2D_ctf


def calc_fcc(current_fft, ref_fft, fft_r, fft_s, fft_x, fft_y, fft_z, cone_point, angle):
    '''
    Calculate fcc for a cone point
    '''

    # Determine the cone mask
    cone_mask = create_fft_cone_mask(fft_s, fft_x, fft_y, fft_z, cone_point, angle)

    # Calc fsc
    fsc3D = calc_fsc(current_fft, ref_fft, fft_r, cone_mask)

    return fsc3D, cone_mask


def create_fft_cone_mask(fft_s, fft_x, fft_y, fft_z, cone_point, angle=10):
    '''
    Create cone mask
    '''
    # Initialize the mask
    cone_mask = np.zeros(fft_s.shape)
    cone_mask = fft_x*cone_point[0] + fft_y*cone_point[1] + fft_z*cone_point[2]

    # Get nonzero r-values
    nonzero_r = fft_s > 0
    zero_r    = fft_s == 0

    # Get cone length
    cone_len = np.sqrt(np.sum(cone_point**2))

    # Assign the normalized values
    cone_mask[nonzero_r] = cone_mask[nonzero_r]/(fft_s[nonzero_r]*cone_len)

    # Determine the cosangle threshold
    cosangle_thresh = np.cos(angle/180.0*np.pi)

    # Make values higher than the threshold  1
    valid = cone_mask >= cosangle_thresh

    cone_mask = np.zeros(fft_s.shape)

    # Assign 1 to valid entries
    cone_mask[valid]  = 1
    cone_mask[zero_r] = 1

    return cone_mask


def calc_fsc(current_fft, ref_fft, fft_r, fft_mask=None):
    '''
    Compute frc between two ft
    '''
    fsc3D = np.zeros(current_fft.shape, dtype=np.float32)
    rbins = np.sort(np.unique(fft_r))
    nbins = len(rbins)

    # Get cross correlations
    cross_cc = current_fft*np.conj(ref_fft)
    ref_cc   = np.abs(ref_fft)**2

    # Calculate for half the FFT
    for i in range(nbins):
        mask  = fft_r == rbins[i]
        # Incorporate fft_mask
        if fft_mask is not None:
            mask *= (fft_mask > 0)

        corr  = np.sum(cross_cc[mask])
        norm  = np.sum(ref_cc[mask])

        # Make the assignment only if mask has positive elements
        if np.sum(mask) > 0:
            fsc3D[mask] = np.abs(corr/norm)

    return fsc3D


@jit(nopython=True)
def calc_fsc_numba(current_fft, ref_fft, fft_r, fft_mask=None):
    '''
    Compute frc between two ft
    '''
    fsc3D = np.zeros(current_fft.shape)
    max_r = int(np.max(fft_r))

    # Initialize arrays
    corr_sum = np.zeros(max_r+1)
    self_sum = np.zeros(max_r+1)

    # Get cross correlations
    cross_cc = np.abs(current_fft*np.conj(ref_fft))
    ref_cc   = np.abs(ref_fft)**2

    # Calculate for half the FFT
    for i in range(fsc3D.shape[0]):
        for j in range(fsc3D.shape[1]):
            for k in range(fsc3D.shape[2]):
                r = int(fft_r[i, j, k])

                # Determine mask coef
                if fft_mask is not None:
                    mask_coef = fft_mask[i, j, k]
                else:
                    mask_coef = 1.0

                corr_sum[r] += mask_coef*cross_cc[i, j, k]
                self_sum[r] += mask_coef*ref_cc[i, j, k]

    # Assign to final bin
    corr_sum[max_r] = 1.0
    self_sum[max_r] = 1.0

    # Determine norm sums
    norm_sum = corr_sum/self_sum

    # Assign values
    for i in range(fsc3D.shape[0]):
        for j in range(fsc3D.shape[1]):
            for k in range(fsc3D.shape[2]):
                r = fft_r[i, j, k]

                # Assign the new value
                fsc3D[i, j, k] = norm_sum[r]

    return fsc3D


def calc_frc(current_fft, ref_fft, ctf_r):
    '''
    Compute frc between two ft
    '''
    frc2D = np.ones(current_fft.shape, dtype=np.float32)
    rbins = np.sort(np.unique(ctf_r))
    nbins = len(rbins)
    max_r = np.max(rbins)

    # Get cross correlations
    cross_cc = current_fft*np.conj(ref_fft)
    ref_cc   = np.abs(ref_fft)**2

    # Calculate for half the FFT
    for i in range(nbins-1):
        mask  = ctf_r == rbins[i]
        corr  = np.sum(cross_cc[mask])
        norm  = np.sum(ref_cc[mask])

        frc2D[mask] = np.abs(corr/norm)

    # For the rest assign it to 1
    mask = ctf_r == max_r
    frc2D[mask] = 1.0

    return frc2D


@jit(nopython=True)
def calc_frc_numba(current_fft, ref_fft, ctf_r):
    '''
    Compute frc between two ft
    '''
    frc2D = np.ones(current_fft.shape)
    max_r = np.max(ctf_r)

    # Initialize arrays
    corr_sum = np.zeros(max_r+1)
    self_sum = np.zeros(max_r+1)

    # Get cross correlations
    cross_cc = np.abs(current_fft*np.conj(ref_fft))
    ref_cc   = np.abs(ref_fft)**2

    # Calculate for half the FFT
    for i in range(frc2D.shape[0]):
        for j in range(frc2D.shape[1]):
            r = ctf_r[i, j]
            corr_sum[r] += cross_cc[i, j]
            self_sum[r] += ref_cc[i, j]

    # Assign to final bin
    corr_sum[max_r] = 1.0
    self_sum[max_r] = 1.0

    # Determine norm sums
    norm_sum = corr_sum/self_sum

    # Assign values
    for i in range(frc2D.shape[0]):
        for j in range(frc2D.shape[1]):
            r = ctf_r[i, j]

            # Assign the new value
            frc2D[i, j] = norm_sum[r]

    return frc2D


def normalize_frc(img2D_fft, frc2D):
    '''
    Normalize frc
    '''

    return img2D_fft*frc2D


def calc_dot(current_img2D, ref_img2D, mask=None):
    '''
    Calc img dot
    '''
    if mask is not None and mask.shape == current_img2D.shape:
        cross_cc = np.average(current_img2D*ref_img2D, weights=mask)
        self_cc  = np.average(ref_img2D*ref_img2D, weights=mask)
    else:
        cross_cc = np.mean(current_img2D*ref_img2D)
        self_cc  = np.mean(ref_img2D*ref_img2D)

    return cross_cc/self_cc


def calc_intensity_ratio(current_img2D, ref_img2D, mask=None):
    '''
    Calc the mean intensity ratio between the images
    '''
    if mask is not None and mask.shape == current_img2D.shape:
        current_mean = np.average(current_img2D, weights=mask)
        ref_mean     = np.average(ref_img2D, weights=mask)
    else:
        current_mean = np.mean(current_img2D)
        ref_mean     = np.mean(ref_img2D)

    return 1.0*current_mean/ref_mean


def normalize_dot(img2D, scale):
    '''
    Normalize by dot
    '''
    return scale*img2D


def calc_mean_std_intensity(img2D, mask):
    '''
    Calculate mean and std intensity
    '''
    if mask is not None and mask.shape == img2D.shape:
        mean_intensity = np.average(img2D, weights=mask)
        std_intensity  = np.sqrt(np.average((img2D-mean_intensity)**2, weights=mask))
    else:
        mean_intensity = np.mean(img2D)
        std_intensity  = np.std(img2D)

    return mean_intensity, std_intensity


def threshold_above(img2D, threshold=1.0):
    '''
    Threshold function
    '''
    img2D[img2D > 1.0] = 1.0

    return img2D


def threshold_below(img2D, threshold=0.0):
    '''
    Threshold function
    '''
    img2D[img2D < 0.0] = 0.0

    return img2D


def create_bg_mask(mask_align_img2D, mask_subtract_img2D):
    '''
    Create bg mask for subtraction from align and subtract masks
    '''
    mask_bg_img2D = 1 - mask_align_img2D + mask_subtract_img2D
    mask_bg_img2D = threshold_above(mask_bg_img2D)
    mask_bg_img2D = threshold_below(mask_bg_img2D)

    return mask_bg_img2D


def subtract_class_ctf(class_img2D, ctf_a, ctf_s, ctf_r, ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D, subtract_bg, norm_method, do_intact_until_first_peak):
    '''
    Subtract class
    '''
    ptcl_img2D  = read_ptcl_mrc(ptcl_star)
    class_ctf   = eval_ptcl_ctf(ctf_s, ctf_a, ptcl_star, 0, do_intact_until_first_peak)
    class_img2D = inv_transform_imgs(class_img2D, ptcl_star)
    mask_align_img2D, mask_structure_img2D, mask_subtract_img2D = inv_transform_masks(ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D)

    # Take class FFT
    class_fft2D = fft_img2D(class_img2D, mask_structure_img2D)

    # CTF correct
    class_fft2D = correct_fft_ctf(class_fft2D, class_ctf)

    # Create frc coefficients
    fft_coeff   = np.ones(class_fft2D.shape)
    real_coeff  = 1.0

    if norm_method == 'frc':
        ptcl_fft2D = fft_img2D(ptcl_img2D)
        fft_coeff  = calc_frc(ptcl_fft2D, class_fft2D,  ctf_r)
    elif norm_method == 'ccc':
        masked_class_img2D = ifft_img2D(class_fft2D)
        real_coeff = calc_dot(ptcl_img2D, masked_class_img2D)
    elif norm_method == 'intensity':
        masked_class_img2D = ifft_img2D(class_fft2D)
        real_coeff = calc_intensity_ratio(ptcl_img2D, masked_class_img2D)

    ptcl_img2D  = subtract_class_from_ptcl(class_img2D, class_ctf, ptcl_img2D, mask_subtract_img2D, fft_coeff, real_coeff)

    return ptcl_img2D


def crop_class_ctf(class_img2D, ctf_a, ctf_s, ctf_r, ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D, subtract_bg, norm_method, do_intact_until_first_peak):
    '''
    Subtract class
    '''

    ptcl_img2D  = read_ptcl_mrc(ptcl_star)
    class_ctf   = eval_ptcl_ctf(ctf_s, ctf_a, ptcl_star, 0, do_intact_until_first_peak)
    mask_align_img2D, mask_structure_img2D, mask_subtract_img2D = inv_transform_masks(ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D)

    class_img2D = generate_noise(class_img2D)
    class_img2D = ctf_correct_class_img(class_img2D, class_ctf)
    class_img2D = norm_intensity_noise_img(class_img2D, ptcl_img2D, mask_align_img2D, mask_structure_img2D)

    # Determine the mask for bg subtraction
    if subtract_bg:
        mask_bg_img2D = create_bg_mask(mask_align_img2D, mask_subtract_img2D)
    else:
        mask_bg_img2D = mask_subtract_img2D

    ptcl_img2D  = crop_class_from_ptcl(class_img2D, ptcl_img2D, mask_bg_img2D)

    return ptcl_img2D


def crop_class(class_img2D, ctf_a, ctf_s, ctf_r, ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D, subtract_bg, norm_method, do_intact_until_first_peak):
    '''
    Crop class with no ctf correction
    '''
    ptcl_img2D  = read_ptcl_mrc(ptcl_star)
    class_img2D = inv_transform_imgs(class_img2D, ptcl_star)
    mask_align_img2D, mask_structure_img2D, mask_subtract_img2D = inv_transform_masks(ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D)

    # Get the mask for noise stats
    mask_background_img2D = mask_align_img2D - mask_structure_img2D
    background_mean, background_std = calc_mean_std_intensity(ptcl_img2D, mask_background_img2D)

    # Determine the mask for bg subtraction
    if subtract_bg:
        mask_bg_img2D = create_bg_mask(mask_align_img2D, mask_subtract_img2D)
    else:
        mask_bg_img2D = mask_subtract_img2D

    ptcl_img2D = generate_noise(ptcl_img2D, mask_bg_img2D, background_mean, background_std)

    return ptcl_img2D


def inv_transform_imgs(class_img2D, ptcl_star):
    '''
    Inverse transform imgs
    '''

    # Inverse transform class img2D
    class_img2D = inv_transform_ptcl_img2D(class_img2D, ptcl_star)

    return class_img2D


def inv_transform_masks(ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D):
    '''
    Inverse transform only masks
    '''
    mask_align_img2D     = inv_transform_ptcl_img2D(mask_align_img2D, ptcl_star)
    mask_structure_img2D = inv_transform_ptcl_img2D(mask_structure_img2D, ptcl_star)
    mask_subtract_img2D  = inv_transform_ptcl_img2D(mask_subtract_img2D, ptcl_star)

    return mask_align_img2D, mask_structure_img2D, mask_subtract_img2D


def ctf_correct_class_img(class_img2D, class_ctf):
    '''
    Compare class to ptcl and adjust intensities
    '''
    # Take class FFT
    class_fft2D = fft_img2D(class_img2D)

    # CTF correct
    class_fft2D = correct_fft_ctf(class_fft2D, class_ctf)

    # Take inverse FFT
    class_img2D = ifft_img2D(class_fft2D)

    return class_img2D


def gaussian_filter(img2D, sigma):
    '''
    Gaussian filter
    '''
    return scipy.ndimage.filters.gaussian_filter(img2D, sigma)


def norm_intensity(img2D, mask=None, new_mean=0, new_std=None):
    '''
    Normalize intensity to match a new gauss-distribution
    '''
    mean_intensity, std_intensity = calc_mean_std_intensity(img2D, mask)

    # Adjust the mean and std-values value
    # Zero mean value
    img2D -= mean_intensity

    # Adjust stdev
    if new_std is not None:
        img2D *= new_std/std_intensity

    # Bring mean value to new mean-value
    img2D += new_mean

    return img2D


def norm_intensity_bg_and_area(img2D, mask_bg, new_val_bg, mask_area, new_val_area):
    '''
    Normalize the intensity - background and an area of interest
    '''
    # Get the background intensity
    background_intensity = np.average(img2D, weights=mask_bg)

    # Subtract background intensity
    img2D -= background_intensity

    # Get the area intensity
    area_intensity = np.average(img2D, weights=mask_area)

    # Normalize the area intensity
    img2D    *= (new_val_area-new_val_bg)/area_intensity

    # Finally add the new background intenstiy
    img2D    += new_val_bg

    return img2D


def norm_intensity_class_img(class_img2D, class_ctf, ptcl_img2D, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D):
    '''
    Intensity normalize class img2D
    '''

    # Set background mask
    mask_background_img2D = mask_align_img2D - mask_structure_img2D

    # Calculate background mean, std intensity
    background_mean, background_std = calc_mean_std_intensity(ptcl_img2D, mask_background_img2D)

    # Calculate structure mean, std intensity
    structure_mean, structure_std   = calc_mean_std_intensity(ptcl_img2D, mask_structure_img2D)

    # Calculate subtract mean, std intensity
    subtract_mean, subtract_std     = calc_mean_std_intensity(ptcl_img2D, mask_subtract_img2D)

    # Set background intensity to particle background mean intensity
    class_img2D = norm_intensity_bg_and_area(class_img2D, mask_background_img2D, background_mean, mask_structure_img2D, structure_mean)

    # Store the original class mrc
    class_img2D = subtract_ctf(class_img2D, class_ctf)

    return class_img2D


def norm_intensity_noise_img(class_img2D, ptcl_img2D, mask_align_img2D, mask_structure_img2D):
    '''
    Intensity normalize class img2D
    '''

    # Set background mask
    mask_background_img2D = mask_align_img2D - mask_structure_img2D

    # Calculate background mean, std intensity
    background_mean, background_std = calc_mean_std_intensity(ptcl_img2D, mask_background_img2D)

    # Set background intensity to particle background mean intensity
    class_img2D = norm_intensity(class_img2D, mask=mask_background_img2D, new_mean=background_mean, new_std=background_std)

    return class_img2D


def subtract_class_from_ptcl(class_img2D, class_ctf, ptcl_img2D, mask_subtract_img2D, fft_coeff, real_coeff):
    '''
    Subtract class from ptcl
    '''

    # Take class FFT
    class_fft2D = fft_img2D(class_img2D, mask_subtract_img2D)

    # CTF correct class image
    class_fft2D = fft_coeff*correct_fft_ctf(class_fft2D, class_ctf)

    # Take inverse FFT
    class_img2D = real_coeff*ifft_img2D(class_fft2D)

    # Subtract the class image from ptcl image
    return ptcl_img2D - class_img2D


def crop_class_from_ptcl(class_img2D, ptcl_img2D, mask_subtract_img2D):
    '''
    Subtract class from ptcl
    '''

    # Subtract the class image from ptcl image
    mask_background_img2D = 1.0-mask_subtract_img2D
    ptcl_img2D = ptcl_img2D*mask_background_img2D + class_img2D*mask_subtract_img2D

    return ptcl_img2D
