#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2019-01-21 08:17:52
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import mrcfile
import numpy as np
import scipy.ndimage


def ccc(current_img2D, other_img2D, mask=None):
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
    img2D = flipX(img2D, ptcl_star)

    return img2D


def is_flipX(ptcl_star):
    '''
    Check if particle needs to be flipped
    '''
    if 'rlnIsFlip' in ptcl_star and ptcl_star['rlnIsFlip'] == 1:
        return True
    else:
        return False


def flipX(img2D, ptcl_star):
    '''
    Flip image around X-axis
    '''
    new_img2D = img2D
    if is_flipX(ptcl_star):
        new_img2D = img2D[:, ::-1]

    return new_img2D


def create_noise(img2D, mask=None, noise_mean=0.0, noise_std=1.0, sigma=0):
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


def read_ptcl_mrc(ptcl_star, transform=False, fft_mask=None):
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
            img2D = pass_filter(img2D, fft_mask)


    return np.copy(img2D)

def pass_filter(img2D, fft_mask):
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


def eval_ptcl_ctf(ctf_s, ctf_a, ptcl_star, bf=0, lp=None, hp=None):
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

    return eval_ctf(ctf_s, ctf_a, defU, defV, defA, phaseShift, kV, ac, cs, bf, lp, hp)


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


def eval_ctf(ctf_s, ctf_a, defU, defV, defA=0, phaseShift=0, kv=300, ac=0.1, cs=2.0, bf=0, lp=None, hp=None):
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
    k3 = np.sqrt(1 - ac**2)
    k4 = bf / 4.                # B-factor, follows RELION convention.
    k5 = np.deg2rad(phaseShift)  # Phase shift.

    # Assign s grid
    s = ctf_s

    # Hard low-pass filter
    if lp is not None:
        s = ctf_s*(ctf_s <= (1. / lp))

    # Hard high-pass filter
    if hp is not None:
        s = ctf_s*(ctf_s >= (1. / hp))

    s2 = s**2
    s4 = s2**2
    dZ = def_avg + def_dev * (np.cos(2 * (ctf_a - defA)))
    gamma = (k1 * dZ * s2) + (k2 * s4) - k5

    # Determine ctf
    img2D_ctf = -(k3 * np.sin(gamma) - ac*np.cos(gamma))

    # Enforce envelope
    if bf != 0:
        img2D_ctf *= np.exp(-k4 * s2)

    return img2D_ctf


def calc_frc(current_fft, other_fft, ctf_r):
    '''
    Compute frc between two ft
    '''
    frc2D = np.ones(current_fft.shape, dtype=np.float32)
    rbins = np.sort(np.unique(ctf_r))

    for i in range(len(rbins)):
        mask = ctf_r == rbins[i]
        corr  = np.sum(current_fft[mask]*np.conj(other_fft[mask]))
        norm1 = np.sqrt(np.sum(np.abs(current_fft[mask])**2))
        norm2 = np.sqrt(np.sum(np.abs(other_fft[mask])**2))

        frc2D[mask] = np.real(corr)/(norm1*norm2)

    return frc2D


def normalize_frc(img2D_fft, frc2D):
    '''
    Normalize frc
    '''
    
    return img2D_fft*frc2D


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


def normalize_intensity(img2D, mask=None, new_mean=0, new_std=None):
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


def normalize_bg_area_intensity(img2D, mask_bg, new_val_bg, mask_area, new_val_area):
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


def subtract_class_ctf(class_img2D, ctf_a, ctf_s, ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D, lowpass=None, highpass=None, subtract_bg=True):
    '''
    Subtract class
    '''
    ptcl_img2D  = read_ptcl_mrc(ptcl_star)
    class_ctf   = eval_ptcl_ctf(ctf_s, ctf_a, ptcl_star, bf=0, lp=lowpass, hp=None)
    class_img2D, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D = inv_transform_imgs(class_img2D, ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D)

    class_img2D = ctf_correct_class_img(class_img2D, class_ctf)

    class_img2D = intensity_norm_class_img(class_img2D, class_ctf, ptcl_img2D, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D*mask_structure_img2D)
    ptcl_img2D  = subtract_class_from_ptcl(class_img2D, class_ctf, ptcl_img2D, mask_subtract_img2D)

    return ptcl_img2D


def crop_class_ctf(class_img2D, ctf_a, ctf_s, ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D, lowpass=None, highpass=None, subtract_bg=True):
    '''
    Subtract class
    '''

    ptcl_img2D  = read_ptcl_mrc(ptcl_star)
    class_ctf   = eval_ptcl_ctf(ctf_s, ctf_a, ptcl_star, bf=0, lp=lowpass, hp=None)
    mask_align_img2D, mask_structure_img2D, mask_subtract_img2D = inv_transform_masks(ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D)
    class_img2D = create_noise(class_img2D)
    class_img2D = ctf_correct_class_img(class_img2D, class_ctf)
    class_img2D = intensity_norm_noise_img(class_img2D, ptcl_img2D, mask_align_img2D, mask_structure_img2D)

    # Determine the mask for bg subtraction
    if subtract_bg:
        mask_bg_img2D = 1 - mask_align_img2D + mask_subtract_img2D
        mask_bg_igm2D = threshold_above(mask_bg_img2D)
        mask_bg_igm2D = threshold_below(mask_bg_img2D)
    else:
        mask_bg_img2D = mask_subtract_img2D

    ptcl_img2D  = crop_class_from_ptcl(class_img2D, ptcl_img2D, mask_bg_img2D)

    return ptcl_img2D


def crop_class(class_img2D, ctf_a, ctf_s, ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D, lowpass=None, highpass=None, subtract_bg=True):
    '''
    Crop class with no ctf correction
    '''
    ptcl_img2D  = read_ptcl_mrc(ptcl_star)
    class_img2D, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D = inv_transform_imgs(class_img2D, ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D)

    # Get the mask for noise stats
    mask_noise_img2D = mask_align_img2D-mask_structure_img2D
    background_mean, background_std = calc_mean_std_intensity(ptcl_img2D, mask_noise_img2D)

    # Determine the mask for bg subtraction
    if subtract_bg:
        mask_bg_img2D = 1 - mask_align_img2D + mask_subtract_img2D
        mask_bg_igm2D = threshold_above(mask_bg_img2D)
        mask_bg_igm2D = threshold_below(mask_bg_img2D)
    else:
        mask_bg_img2D = mask_subtract_img2D

    ptcl_img2D = create_noise(ptcl_img2D, mask_bg_img2D, background_mean, background_std)

    return ptcl_img2D


def inv_transform_imgs(class_img2D, ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D):
    '''
    Inverse transform imgs and masks
    '''
    # Inverse transform the masks
    mask_align_img2D = inv_transform_ptcl_img2D(mask_align_img2D, ptcl_star)
    mask_structure_img2D = inv_transform_ptcl_img2D(mask_structure_img2D, ptcl_star)
    mask_subtract_img2D = inv_transform_ptcl_img2D(mask_subtract_img2D, ptcl_star)

    # Inverse transform class img2D
    class_img2D = inv_transform_ptcl_img2D(class_img2D, ptcl_star)

    return class_img2D, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D

def inv_transform_masks(ptcl_star, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D):
    '''
    Inverse transform only masks
    '''
    mask_align_img2D = inv_transform_ptcl_img2D(mask_align_img2D, ptcl_star)
    mask_structure_img2D = inv_transform_ptcl_img2D(mask_structure_img2D, ptcl_star)
    mask_subtract_img2D = inv_transform_ptcl_img2D(mask_subtract_img2D, ptcl_star)

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


def intensity_norm_class_img(class_img2D, class_ctf, ptcl_img2D, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D):
    '''
    Intensity normalize class img2D
    '''

    # Set background mask
    mask_background_img2D = mask_align_img2D - mask_structure_img2D

    # Calculate background mean, std intensity
    background_mean, background_std = calc_mean_std_intensity(ptcl_img2D, mask_background_img2D)

    # Calculate structure mean, std intensity
    structure_mean, structure_std = calc_mean_std_intensity(ptcl_img2D, mask_structure_img2D)

    # Calculate subtract mean, std intensity
    subtract_mean, subtract_std = calc_mean_std_intensity(ptcl_img2D, mask_subtract_img2D)

    # Set background intensity to particle background mean intensity
    class_img2D = normalize_bg_area_intensity(class_img2D, mask_background_img2D, background_mean, mask_subtract_img2D, subtract_mean)

    # Store the original class mrc
    class_img2D = subtract_ctf(class_img2D, class_ctf)

    return class_img2D


def intensity_norm_noise_img(class_img2D, ptcl_img2D, mask_align_img2D, mask_structure_img2D):
    '''
    Intensity normalize class img2D
    '''

    # Set background mask
    mask_background_img2D = mask_align_img2D - mask_structure_img2D

    # Calculate background mean, std intensity
    background_mean, background_std = calc_mean_std_intensity(ptcl_img2D, mask_background_img2D)

    # Set background intensity to particle background mean intensity
    class_img2D = normalize_intensity(class_img2D, mask=mask_background_img2D, new_mean=background_mean, new_std=background_std)

    return class_img2D


def subtract_class_from_ptcl(class_img2D, class_ctf, ptcl_img2D, mask_subtract_img2D):
    '''
    Subtract class from ptcl
    '''

    # Take class FFT
    class_fft2D = fft_img2D(class_img2D, mask_subtract_img2D)

    # CTF correct class image
    class_fft2D = correct_fft_ctf(class_fft2D, class_ctf)

    # Take inverse FFT
    class_img2D = ifft_img2D(class_fft2D)

    # Subtract the class image from ptcl image
    return ptcl_img2D - class_img2D


def crop_class_from_ptcl(class_img2D, ptcl_img2D, mask_subtract_img2D):
    '''
    Subtract class from ptcl
    '''

    # Subtract the class image from ptcl image
    ptcl_img2D = ptcl_img2D*(1.0-mask_subtract_img2D)+class_img2D*mask_subtract_img2D

    return ptcl_img2D
