#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 11:27:03
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import re
import glob
import numpy as np
import pandas as pd
import utilities as util
import mrcfile
import subprocess
import sys
import scipy.ndimage
import parallelem
import multiprocessing
import shutil
import sqlite3
import matplotlib.pyplot as py
import barcode

from collections import Counter
from mpl_toolkits.mplot3d import Axes3D
from shutil import copyfile
from matplotlib.colors import LogNorm


class Relion:
    def __init__(self):
            self.name = None


class Project:
    '''
        Project File
    '''
    def __init__(self, name='EMProject'):
        self.name            = name
        self.particle_star   = None
        self.ref_class_star  = None
        self.micrograph_star = None
        self.class_ids       = None

        # Particle props
        self.particle_diameter_A  = None
        self.particle_radius_pix  = None

        # Input files
        self.particle_star_file   = None
        self.ref_class_star_file  = None
        self.micrograph_star_file = None

        # Output files
        self.particle_out_file  = None
        self.ref_class_out_file = None

        # Micrograph pixel size
        self.mic_apix           = None
        self.particle_apix      = None
        self.ref_apix           = None

        # Micrograph files
        self.first_mic_file  = None
        self.first_mic_mrc   = None

        # MRC files and objects
        self.particle_mrc        = None
        self.ref_class_mrc       = None
        self.consensus_class_mrc = None

        self.particle_mrc_file        = None
        self.ref_class_mrc_file       = None
        self.consensus_class_mrc_file = None

        # Metadata file
        self.metadata      = None
        self.metadata_file = None

        # Micrograph dimensions
        self.mic_NX = 0
        self.mic_NY = 0
        self.mic_NZ = 0

        # Cryosparc objects
        self.particle_cs    = None
        self.ref_class_cs   = None

        # Cryosparc files
        self.blob_cs_file        = None
        self.ref_class_cs_file   = None
        self.passthrough_cs_file = None
        self.original_star_file  = None

        # Additional data frames
        self.particle_data_props = pd.DataFrame(columns=['insideFrame'])

        # Alignment references
        self.ref_align_star_file  = None
        self.ref_align_mrc_file   = None

        # Alignment mask file
        self.mask_align_mrc_file    = None
        self.mask_subtract_mrc_file = None

        # First particle and class mrc files
        self.first_particle_mrc       = None
        self.first_particle_mrc_file  = None

        self.first_ref_class_mrc      = None
        self.first_ref_class_mrc_file = None

        # Low and highpass filters
        self.highpass_angstrom        = None
        self.lowpass_angstrom         = None

        # Relion code
        self.relion_refine_exe    = 'relion_refine'
        self.relion_refine_args   = []

        self.relion_norm_exe      = 'relion_preprocess'
        self.relion_norm_args     = []

        self.relion_image_handler_exe = 'relion_image_handler'
        self.relion_flip_args         = []
        self.relion_noflip_args       = []

        # Particles and class2D models
        self.particles = []
        self.class2Ds  = []
        self.class3Ds  = []

        # Star files to merge
        self.particle_star_files = []
        self.other_star          = None

        # Metadata file
        self.metadata      = None
        self.metadata_file = None

        # Geometrical parameters
        self.direction = None
        self.diff_tilt = None
        self.diff_psi  = None
        self.diff_rot  = None

    def write_metadata(self):
        '''
        Write metadata
        '''
        if self.metadata is not None and self.metadata_file is not None:
            self.metadata.to_csv(self.metadata_file, header=True, index=False, sep='\t')

    def set_particle_num(self, num=None):
        '''
        Set particle number
        '''
        if num is not None and self.particle_star is not None:
            self.particle_star.data_block = self.particle_star.data_block.loc[:num, :]

    def read_ptcl_mrc_paths(self):
        '''
        Read ptcl mrc paths
        '''
        self.particle_star.read_mrc_paths()

    def group_ptcls(method='defocus'):
        '''
        Grouping method
        '''

        # Sort particles using defocus or intensity scale

    def pick_random_set(self, num_ptcls, rand_seed=1):
        '''
        Pick random subset of the particles
        '''
        self.particle_star.pick_random_set(num_ptcls, rand_seed)

    def invert_psi(self):
        '''
        Invert psi angle
        '''
        if self.particle_star.has_label('rlnAnglePsi'):
            self.particle_star.data_block['rlnAnglePsi'] *= -1.0

    def invert_angles(self):
        '''
        Invert psi/rot/tilt angles
        '''
        if self.particle_star.has_label('rlnAnglePsi'):
            self.particle_star.data_block['rlnAnglePsi'] *= -1.0

        if self.particle_star.has_label('rlnAngleRot'):
            self.particle_star.data_block['rlnAngleRot'] *= -1.0

        if self.particle_star.has_label('rlnAngleTilt'):
            self.particle_star.data_block['rlnAngleTilt'] *= -1.0

    def invert_origin(self):
        '''
        Invert origin
        '''
        if self.particle_star.has_label('rlnOriginX'):
            self.particle_star.data_block['rlnOriginX'] *= -1.0

        if self.particle_star.has_label('rlnOriginY'):
            self.particle_star.data_block['rlnOriginY'] *= -1.0

    def set_highpass_filter(self, hp=None):
        '''
        Set highpass filter
        '''
        if hp is not None:
            self.highpass_angstrom = hp

    def set_lowpass_filter(self, lp=None):
        '''
        Set lowpass filter
        '''
        if lp is not None:
            self.lowpass_angstrom = lp

    def set_particle_radius(self):
        '''
        Set particle radius from particle diameter in Angstrom
        '''
        if self.particle_diameter_A is not None and self.particle_apix is not None:
            self.particle_radius_pix = int(self.particle_diameter_A//(2*self.particle_apix))

    def select_by_barcode(self, barcode_list):
        '''
        Select by barcode
        '''
        self.particle_star.select_by_barcode(barcode_list)

    def delete_by_barcode(self, barcode_list):
        '''
        Select by barcode
        '''
        self.particle_star.delete_by_barcode(barcode_list)

    def apply_barcode(self):
        '''
        Apply barcode on particle star
        '''
        new_data_block = barcode.Frame_angle(self.particle_star.get_data_block())
        self.particle_star.set_data_block(new_data_block)

    def append_particle_barcode(self, barcode={}):
        '''
        Append particle barcode
        '''
        if self.particle_star is not None:
            self.particle_star.append_barcode(barcode)

    def set_particle_barcode(self, barcode={}):
        '''
        Append particle barcode
        '''
        if self.particle_star is not None:
            self.particle_star.set_barcode(barcode)

    def rename_columns(self, column_params):
        '''
        Rename columns
        '''
        self.particle_star.rename_columns(column_params)

    def flipX_particles(self):
        '''
        Flip particles in star file
        '''
        self.particle_star.flipX()

    def check_particle_pos(self):
        '''
        Check location of all particles
        '''
        num_ptcls = self.particle_star.get_numRows()
        ptcl_list = np.arange(num_ptcls)
        ptcl_pos_list = []

        for ptcl in ptcl_list:
            isInside = self.particle_star.is_particle_inside(ptcl, self.mic_apix, self.mic_NX, self.mic_NY)
            ptcl_pos_list.append(isInside)

        ptcl_pos_list = np.array(ptcl_pos_list)
        self.particle_data_props['insideFrame'] = ptcl_pos_list

    def delete_outside_ptcls(self):
        '''
        Delete particles outside the frame
        '''
        delete_ptcls = np.where(self.particle_data_props['insideFrame'] == False)[0]
        self.particle_star.delete_ptcls(delete_ptcls)

    def read_mic_header(self):
        '''
        Read the header from first micrograph
        '''
        if self.micrograph_star is not None:
            if self.micrograph_star.has_label('rlnMicrographName'):
                self.first_mic_file = self.micrograph_star.data_block.loc[0, 'rlnMicrographName']
                self.first_mic_mrc  = MRC(self.first_mic_file)
                self.set_mic_dimensions()

    def set_mic_dimensions(self):
        '''
        Set micrograph dimensions
        '''
        if 'NX' in self.first_mic_mrc.header.dtype.names:
            self.mic_NX = self.first_mic_mrc.header['NX']
            self.mic_NY = self.first_mic_mrc.header['NY']
            self.mic_NZ = self.first_mic_mrc.header['NZ']
        elif len(self.first_mic_mrc.img3D.shape) == 3:
            self.mic_NZ, self.mic_NY, self.mic_NX = self.first_mic_mrc.img3D.shape
        else:
            self.mic_NY, self.mic_NX = self.first_mic_mrc.img3D.shape

    def set_mic_apix(self, apix=1.82):
        '''
        Set micrograph apix
        '''
        self.mic_apix = apix

    def set_particle_diameter(self, diameter):
        '''
        Set particle diameter in Angstroms
        '''
        self.particle_diameter_A = diameter

    def set_particle_apix(self, apix=1.82):
        '''
        Set particle apix
        '''
        self.particle_apix = apix

    def set_ref_class_apix(self, apix=1.82):
        '''
        Set particle apix
        '''
        self.ref_class_apix = apix

    def read_mic_apix(self):
        '''
        Read and set micrograph apix
        '''
        self.micrograph_star.determine_star_apix()
        self.set_mic_apix(self.micrograph_star.get_star_apix())

    def read_particle_apix(self):
        '''
        Read and set micrograph apix
        '''
        self.particle_star.determine_star_apix()
        self.set_particle_apix(self.particle_star.get_star_apix())

    def read_ref_apix(self):
        '''
        Read ref apix
        '''
        self.ref_class_star.determine_star_apix()
        self.set_ref_class_apix(self.ref_class_star.get_star_apix())

    def read_class_refs(self, file, new_classids=False):
        '''
        Read class refs
        '''
        self.ref_class_star_file = os.path.abspath(file)
        self.ref_class_star = Star(file)

        # Generate new classids from particle numbers in image names
        if new_classids:
            self.ref_class_star.num2className()

    def read_particle_stars(self, files):
        '''
        Read a batch of star files to process
        '''
        # Count number of files
        file_counter = 0
        for file in files:
            glob_files = glob.glob(file)

            # Iterate over glob files
            for glob_file in glob_files:

                # Get the file extension
                head, ext = os.path.splitext(glob_file)
                if os.path.isfile(glob_file) and ext == '.star':
                    # Read the first file
                    if file_counter == 0:
                        self.read_particles(glob_file)
                    else:
                        self.read_other_particles(glob_file)

                    # Update file counter
                    file_counter += 1

                    # Merge the star files
                    self.particle_star.merge_star(self.other_star)

    def read_particles(self, file):
        '''
        Read particle star
        '''
        self.particle_star_file = os.path.abspath(file)
        self.particle_star = Star(file)

    def read_other_particles(self, file):
        '''
        Read other star
        '''
        self.other_star_file = os.path.abspath(file)
        self.other_star      = Star(file)

    def read_micrographs(self, file):
        '''
        Read micrograph star file
        '''
        self.micrograph_star_file = os.path.abspath(file)
        self.micrograph_star      = Star(file)

    def sort_micrographs(self, column='rlnDefocusU'):
        '''
        Sort micrographs
        '''
        self.micrograph_star.sort(column=column)

    def get_class_ids(self):
        '''
        Get class names
        '''
        self.class_ids =  self.particle_star.get_class_ids()

    def recenter_particles(self):
        '''
        Recenter particles
        '''
        self.particle_star.determine_star_apix()
        self.particle_star.recenter2D(mic_apix=self.mic_apix)

    def copy_from_ref(self, columns=[], compare='img'):
        '''
        Copy columns from reference star
        '''

        if compare == 'img':
            cmp_columns = ['shortImageName']

            # Create short Image name
            self.particle_star.create_shortImageName()
            self.ref_class_star.create_shortImageName()

        else:
            cmp_columns = ['rlnMicrographName', 'rlnCoordinateX', 'rlnCoordinateY']

        # Delete copy columns
        self.particle_star.delete_columns(columns)

        # Shrink reference star data block
        self.ref_class_star.keep_columns(columns+cmp_columns)

        # Merge the two data sets
        self.particle_star.merge(self.ref_class_star, cmp_columns)
        
        # Delete accesory columns
        if compare == 'img':
            self.particle_star.delete_shortImageName()

    def copy_columns(self, column_params):
        '''
        Copy from one column to another new column in particle star file
        '''
        if column_params is not None:
            for from_column, to_column in column_params.items():
                self.particle_star.copy_column2column(from_column, to_column)

    def reset_offsets(self):
        '''
        Make all offsets 0
        '''
        offset_columns = ['rlnOriginX',
                          'rlnOriginY',
                          'rlnAnglePsi',
                          'rlnAngleTilt',
                          'rlnAngleRot']
        for column in offset_columns:
            self.particle_star.set_column(column, 0)

    def add_columns(self, column_params=None):
        '''
        Add new columns
        '''
        if column_params is not None:
            for label, value in column_params.items():
                self.particle_star.set_column(label, value)

    def delete_columns(self, column_params=None):
        '''
        Delete columns
        '''
        if column_params is not None:
            for label, value in column_params.items():
                self.particle_star.delete_column(label)

    def reset_priors(self):
        '''
        Delete prior columns
        '''
        prior_columns = ['rlnOriginXPrior',
                         'rlnOriginYPrior',
                         'rlnAnglePsiPrior',
                         'rlnAngleRotPrior',
                         'rlnAngleTiltPrior']

        for label in prior_columns:
            self.particle_star.delete_column(label)

    def toggle_flip(self):
        '''
        Set flip on for particles
        '''
        if self.particle_star:
            if self.particle_star.has_label('rlnIsFlip'):
                flipON  = self.particle_star.data_block['rlnIsFlip'] == 1
                flipOFF = self.particle_star.data_block['rlnIsFlip'] == 0

                self.particle_star.data_block.loc[flipON, 'rlnIsFlip'] = 0
                self.particle_star.data_block.loc[flipOFF,'rlnIsFlip'] = 1
            else:
                self.particle_star.set_column('rlnIsFlip', 1)

    def transform_particles(self, final_offset=[0, 0], com_offset=False, rotate_psi=0):
        '''
        Transform particle star file based on the class star file
        '''
        if self.particle_star is None:
            print('No transformation due to missing particle data')
            return 0

        if self.consensus_class_mrc is not None:
            final_offset = self.consensus_class_mrc.determine_com(img_num=0)

        if self.ref_class_star is not None:
            # Ref data block
            ref_data_block = self.ref_class_star.get_data_block()

            # Iterate through every class
            for i in range(ref_data_block.shape[0]):
                # Get class id
                class_id = ref_data_block['rlnClassNumber'][i]

                # Get rotangle
                if self.ref_class_star.has_label('rlnAnglePsi'):
                    rot_angle = ref_data_block['rlnAnglePsi'][i]
                else:
                    rot_angle = 0.0

                # Get offsetx, offsety
                if self.ref_class_star.has_label('rlnOriginX') and self.ref_class_star.has_label('rlnOriginY'):
                    offset_x = ref_data_block['rlnOriginX'][i]
                    offset_y = ref_data_block['rlnOriginY'][i]
                else:
                    offset_x = 0.0
                    offset_y = 0.0

                # Get class rows
                class_ptcls = self.particle_star.get_class_rows(class_id)

                print("Processing class {:d}. Number of particles {:d}".format(class_id, len(class_ptcls)))

                # Make the transformation
                self.particle_star.rotate2D(rotangle=rot_angle,
                                            offset=[offset_x, offset_y],
                                            final_offset=final_offset,
                                            ptcls=class_ptcls)
        else:
            # Determine the particles
            ptcls = np.arange(self.particle_star.data_block.shape[0])
            # Make the transformation
            self.particle_star.rotate2D(rotangle=0.0,
                                        offset=[0.0, 0.0],
                                        final_offset=final_offset,
                                        ptcls=ptcls)

        # Rotate psi
        self.particle_star.rotate_psi(rotangle=rotate_psi)

        return 1

    def write_output_files(self, write_particle_star=True, write_ref_class_star=True, write_cs_star=True):
        '''
        Write output files
        '''
        if self.particle_star is not None and write_particle_star:
            self.particle_star.write(self.particle_out_file)

        if self.ref_class_star is not None and write_ref_class_star:
            self.ref_class_star.write(self.ref_class_out_file)

        if self.particle_cs is not None and write_cs_star:
            self.particle_cs.star.write(self.particle_out_file)

        if self.ref_class_cs is not None and write_cs_star:
            self.ref_class_cs.star.write(self.ref_class_out_file)

        if self.metadata is not None:
            self.write_metadata()

    def set_output_directory(self, output_directory=None, project_root='.'):
        '''
        Set output directory
        '''

        if output_directory is not None:
            self.output_directory = output_directory
        else:
            # Get project root
            head = project_root

            # List existing output directories
            potential_directories = list(filter(lambda x: os.path.isdir(x),
                                         glob.glob(head+'/'+self.name+'_em_[0-9][0-9][0-9]')))

            # Get the number extensions
            number_extensions = [int(x[-3:]) for x in potential_directories]

            # Get the counter
            output_counter = 1
            if len(number_extensions) > 0:
                output_counter = max(number_extensions)+1

            self.output_directory = head+'/'+self.name+"_em_%03d" % (output_counter)

        # Make directory
        os.makedirs(self.output_directory, exist_ok=True)

    def write_mirror_files(self):
        '''
        Write output files
        '''
        if self.left_star is not None:
            self.left_star.write(self.mirror_left_out_file)

        if self.right_star is not None:
            self.right_star.write(self.mirror_right_out_file)

    def split_mirrors(self):
        '''
        Split mirrors
        '''

        # If the particle star has the flip variable
        if self.particle_star.has_label('rlnIsFlip'):
            # Create left and right stars
            self.left_star  = Star()
            self.right_star = Star()

            # Create masks
            left_mask  = self.particle_star.data_block['rlnIsFlip'] == 0
            right_mask = self.particle_star.data_block['rlnIsFlip'] == 1

            self.left_star.data_block  = self.particle_star.data_block.loc[left_mask, :]
            self.right_star.data_block = self.particle_star.data_block.loc[right_mask, :]

    def prepare_metadata_file(self):
        '''
        Prepare metadata file
        '''
        if self.metadata_file is None:
            self.metadata_file = self.output_directory+'/metadata.txt'

    def prepare_io_files_star(self):
        # Prepare metadata file
        self.prepare_metadata_file()

        # Copy input file to output directory
        if self.particle_star_file is not None:
            head, tail = os.path.split(self.particle_star_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.particle_star_file, self.output_directory+'/particle_input'+ext)
            self.particle_out_file = self.output_directory+'/particle_output'+ext

            # Make symlink
            self.make_symlink2parent(self.particle_star_file)

        if self.ref_class_star_file is not None:
            head, tail = os.path.split(self.ref_class_star_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.ref_class_star_file, self.output_directory+'/class2D_input'+ext)
            self.ref_class_out_file = self.output_directory+'/class2D_output'+ext

    def prepare_mirror_files_star(self):
        # Copy input file to output directory
        if self.particle_star_file is not None:
            head, tail = os.path.split(self.particle_star_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.particle_star_file, self.output_directory+'/'+root+'_particle_input'+ext)
            self.mirror_left_out_file  = self.output_directory+'/'+root+'_particle_left'+ext
            self.mirror_right_out_file = self.output_directory+'/'+root+'_particle_right'+ext

    def prepare_io_files_cs(self):
        # Copy input files to output directory
        if self.blob_cs_file is not None:
            head, tail = os.path.split(self.blob_cs_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.blob_cs_file, self.output_directory+'/blob_input'+ext)
            self.particle_out_file = self.output_directory+'/particle_output.star'

            # Make symlink
            self.make_symlink2parent(self.blob_cs_file)

        if self.passthrough_cs_file is not None:
            head, tail = os.path.split(self.passthrough_cs_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.passthrough_cs_file, self.output_directory+'/passthrough_input'+ext)

        if self.original_star_file is not None:
            head, tail = os.path.split(self.original_star_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.original_star_file, self.output_directory+'/original_particle'+ext)

        if self.ref_class_cs_file is not None:
            head, tail = os.path.split(self.ref_class_cs_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.ref_class_cs_file, self.output_directory+'/class2D_input'+ext)
            self.ref_class_out_file = self.output_directory+'/class2D_output.star'

    def set_cs_files(self, blob_cs_file=None, passthrough_cs_file=None, original_star_file=None, ref_class_cs_file=None):
        '''
        Set input cs files
        '''
        self.blob_cs_file        = blob_cs_file
        self.passthrough_cs_file = passthrough_cs_file
        self.original_star_file  = original_star_file
        self.ref_class_cs_file   = ref_class_cs_file

    def read_cs_files(self):
        '''
        Read cs file
        '''
        self.particle_cs =  CryoSparc()

        if self.blob_cs_file is not None:
            self.particle_cs.read_blob(self.blob_cs_file)

        if self.passthrough_cs_file is not None:
            self.particle_cs.read_passthrough(self.passthrough_cs_file)

        if self.original_star_file is not None:
            self.particle_cs.read_original_star(self.original_star_file)

        if self.ref_class_cs_file is not None:
            self.ref_class_cs = CryoSparc()
            self.ref_class_cs.read_blob(self.ref_class_cs_file)

    def convert_cs2star(self, mic_path='Micrographs', proj_path='', del_classes=[], del_str='', restore_offsets=False, merge_original=False, reset_classes=True):
        '''
        Convert to cs to star file
        '''

        self.particle_cs.set_project_path(proj_path)
        self.particle_cs.convert2star()
        self.particle_cs.copy_from_original(mic_path)

        # Delete classes
        self.particle_cs.delete_classes(del_classes)

        # Reset class names
        if reset_classes:
            self.particle_cs.reset_class_numbers()

        # Remove a string from micrographs name
        self.particle_cs.remove_str_from_micrograph_names(del_str)

        if self.ref_class_cs is not None:
            self.ref_class_cs.set_project_path(proj_path)
            self.ref_class_cs.convert2star()
            self.ref_class_cs.convert_idx_to_classnumber()
            self.ref_class_cs.rename_star_columns(columns={'rlnImageName': 'rlnReferenceImage'})

            # Convert template mrc file to mrcs
            self.ref_class_cs.get_ref_mrc_file()
            self.ref_class_cs.convert_ref_mrc_to_mrcs()
            self.ref_class_cs.rename_ref_star_to_mrcs()

            # Delete unwanted classes
            self.ref_class_cs.delete_classes(del_classes)

            # Reset class names
            if reset_classes:
                self.ref_class_cs.reset_class_numbers()

        # Merge the data from original star file
        if self.particle_cs.original_star is not None and merge_original:
            self.particle_cs.merge_with_original_star(restore_offsets)

    def read_particle_mrc(self, particle_id=0):
        '''
        Read particle mrc
        '''
        if self.particle_star is not None:
            particle_data_block = self.particle_star.get_data_block()

            # Get only single particle data
            if particle_id < particle_data_block.shape[0]:

                # Get particle image num and filename
                image_name = particle_data_block.loc[particle_id, "rlnImageName"]
                image_num, image_file = image_name.split('@')

                # Corrected image number
                cor_image_num             = int(image_num) - 1
                self.particle_mrc         = MRC(image_file, cor_image_num)
                self.particle_mrc.project = self

                # Set star information for the particle
                if self.particle_apix is None:
                    self.particle_apix = self.particle_star.determine_star_apix()

                # Set particle apix
                self.particle_mrc.set_apix(self.particle_apix)

            else:
                self.particle_mrc = None

    def read_ref_class_mrc(self, class_number=1):
        '''
        Read class mrc
        '''
        if self.ref_class_mrc_file is not None:
            # Get corrected image number
            cor_image_num = class_number - 1

            # Get Mrc
            self.ref_class_mrc = MRC(self.ref_class_mrc_file, cor_image_num)

        elif self.ref_class_star is not None:
            ref_class_data_block = self.ref_class_star.get_data_block()

            # Class mask
            class_mask = ref_class_data_block["rlnClassNumber"] == class_number

            if np.sum(class_mask) > 0:

                # Get class image num and filename
                image_name = ref_class_data_block.loc[class_mask, "rlnReferenceImage"]
                image_num, image_file = image_name.split('@')

                # Corrected image number
                cor_image_num      = int(image_num) - 1
                self.ref_class_mrc = MRC(image_file, cor_image_num)
            else:
                self.ref_class_mrc = None

        else:
            self.particle_mrc = None

    def set_relion_output_str(self, name='run'):
        '''
        Set relion output str
        '''
        self.relion_output_str = self.output_directory+'/'+name

    def set_relion_refine_exe(self):
        '''
        Set relion exe
        '''
        relion_process = subprocess.run(['which', 'relion_refine'], stdout=subprocess.PIPE, universal_newlines=True)
        self.relion_refine_exe = relion_process.stdout.strip()

    def set_relion_norm_exe(self):
        '''
        Set relion exe
        '''
        relion_process = subprocess.run(['which', 'relion_preprocess'], stdout=subprocess.PIPE, universal_newlines=True)
        self.relion_norm_exe = relion_process.stdout.strip()

    def set_relion_image_handler_exe(self):
        '''
        Set relion image handler exe
        '''
        relion_process = subprocess.run(['which', 'relion_image_handler'], stdout=subprocess.PIPE, universal_newlines=True)
        self.relion_image_handler_exe = relion_process.stdout.strip()

    def set_relion_stack_create_exe(self):
        relion_process = subprocess.run(['which', 'relion_stack_create'], stdout=subprocess.PIPE, universal_newlines=True)
        self.relion_stack_create_exe = relion_process.stdout.strip()

    def set_structure_mask(self, mask_file):
        '''
        Set structure mask
        '''
        self.mask_structure_mrc_file = mask_file

    def set_ref_class_num(self, num=1):
        '''
        Set ref class number
        '''
        self.ref_class_number = num

    def set_alignment_mask(self, mask_file):
        '''
        Set alignment mask
        '''
        self.mask_align_mrc_file    = mask_file

    def set_alignment_ref(self, ref_file):
        '''
        Set alignment reference
        '''
        self.ref_align_mrc_file     = ref_file

    def set_subtraction_mask(self, mask_file):
        '''
        Set subtraction mask
        '''
        self.mask_subtract_mrc_file = mask_file

    def read_first_particle_mrc(self):
        '''
        Read first particle mrc
        '''
        if self.particle_star is not None:
            # Get first particle image to determine shape parameters
            image_num, image_name = self.particle_star.get_image_num_name(0)

            # Read image
            self.first_particle_mrc = MRC(image_name, int(image_num)-1)

    def read_first_ref_class_mrc(self, ref_num=1, write_ref=True):
        '''
        Read first particle mrc
        '''
        if self.ref_class_star is not None:
            # Get first particle image to determine shape parameters
            image_num, image_name = self.ref_class_star.get_image_num_name(ref_num-1)

            # Read image
            self.first_ref_class_mrc = MRC(image_name, int(image_num)-1)

            if write_ref:
                # Write first ref class image
                self.write_first_ref_class_mrc(image_name, image_num)

    def write_first_ref_class_mrc(self, ref_image_name, ref_image_num):
        '''
        Write first ref class mrc
        '''
        if self.first_ref_class_mrc is not None:
            head, ext = os.path.splitext(ref_image_name)

            # Set 2D data to 3D
            self.first_ref_class_mrc.set_img3D(self.first_ref_class_mrc.get_img2D())

            # Define first class mrc filename
            self.first_ref_class_mrc_file = head+'_%03d.mrcs' % (ref_image_num)

            # Write image
            self.first_ref_class_mrc.write_img(fname=self.first_ref_class_mrc_file, apix=self.particle_apix)

    def read_ptcl_mrc(self, ptcl_star):

        # Read particle data
        particle_image_num, particle_image_name = ptcl_star['rlnImageName'].split('@')

        # Get particle image
        particle_mrc = MRC(particle_image_name, int(particle_image_num)-1)

        return particle_mrc

    def filter_ptcls(self, 
                     maxprob=0.5, 
                     maxclass=10,
                     tilt_range=[0, 360],
                     dpsi_range=[0, 360],
                     dtilt_range=[0, 360],
                     dalign_range=[0, 360]):
        '''
        Filter ptcls
        '''
        if self.particle_star is not None:
            # Filter based on classification parameters
            self.particle_star.filter_classification(maxprob, maxclass)

            # Filter based on orientation parameters
            self.particle_star.filter_orientation(tilt_range,
                                                  dpsi_range,
                                                  dtilt_range,
                                                  dalign_range)


    def make_symlink2parent(self, input_file, out_path='particle_input'):
        '''
        Make symlink to input file folder
        '''
        # Split input file
        head, tail = os.path.split(input_file)

        # Get relative directory to input file folder
        relative_input_dir = os.path.abspath(head)

        # Destination directory
        relative_output_dir = os.path.relpath(os.path.abspath(self.output_directory+'/'+out_path))

        # Create symlink
        os.symlink(relative_input_dir, relative_output_dir)

    def tilt90_particles(self):
        '''
        Make all tilts at 90
        '''
        self.particle_star.tilt90()

    def Zflip_particles(self):
        '''
        Z-flip particles
        '''
        self.particle_star.Zflip()

    def align_to_priors(self, columns):
        '''
        Align particle column values to their priors
        '''
        for column in columns:
            self.particle_star.align_to_prior(column)



class ProjectScale(Project):
    '''
    Star scale project with new set of micrographs
    '''
    def __init__(self, name='EMScale'):
        super().__init__(name)

        self.new_micrograph_star      = None
        self.new_micrograph_star_file = None

    def read_new_micrographs(self, file):
        '''
        Read micrograph star file
        '''
        self.new_micrograph_star_file = os.path.abspath(file)
        self.new_micrograph_star      = Star(file)

    def read_new_mic_apix(self):
        '''
        Read and set micrograph apix
        '''
        self.new_micrograph_star.determine_star_apix()
        self.set_new_mic_apix(self.new_micrograph_star.get_star_apix())

    def set_new_mic_apix(self, apix=1.82):
        '''
        Set micrograph apix
        '''
        self.new_mic_apix = apix

    def scale_particle_star(self):
        '''
        Scale particle star
        '''
        self.coor_scale_ratio   = 1.0*self.mic_apix/self.nex_mic_apix

        # Scale the coordinates first
        if self.particle_star.has_label('rlnCoordinateX'):
            self.particle_star.data_block['rlnCoordinateX'] =  np.round(self.particle_star.data_block['rlnCoordinateX']*self.coor_scale_ratio)

        if self.particle_star.has_label('rlnCoordinateY'):
            self.particle_star.data_block['rlnCoordinateY'] =  np.round(self.particle_star.data_block['rlnCoordinateY']*self.coor_scale_ratio)

    def rename_micrographs(self):
        '''
        Rename micrograph names in particle star
        '''
        if self.particle_star.has_label('rlnMicrographName'):
            # Split head and tail for micrograph name
            ori_head, ori_tail = os.path.split(self.particle_star.data_block['rlnMicrographName'][0])

            # New head and new tail
            new_head, new_tail = os.path.split(self.new_micrograph_star.data_block['rlnMicrographName'][0])

            # Determine deletion and replace strings
            del_str = ori_head+'/'
            new_str = new_head+'/'

            # Replace the Micrographs folder
            self.particle_star.data_block['rlnMicrographName'] = self.particle_star.data_block.rlnMicrographName.replace({del_str: new_str}, regex=False)

    def run(self):
        '''
        Run the project
        '''
        self.scale_particle_star()
        self.rename_micrographs()


class ProjectFsc(Project):
    '''
    Particle flip project
    '''
    def __init__(self, name='EMFsc'):
        super().__init__(name)

        self.half_map1_file = None
        self.half_map2_file = None
        self.whole_map_file = None
        self.mask_file      = None

        # Map objects
        self.half_map1 = None
        self.half_map2 = None
        self.whole_map = None
        self.mask      = None

        self.apix      = None
        self.bfactor   = None
        self.ncones    = None
        self.angle     = None
        self.batch     = None

        # Fibonacci points
        self.fib_points = None

        # Final map
        self.fcc_img3D = None
        self.fcc_map   = None

    def set_params(self, half_map1, half_map2, whole_map, mask, apix, bfactor, ncones, angle, batch):
        '''
        Set project parameters
        '''
        self.half_map1_file = half_map1
        self.half_map2_file = half_map2
        self.whole_map_file = whole_map
        self.mask_file      = mask

        self.apix    = apix
        self.bfactor = bfactor
        self.ncones  = ncones
        self.angle   = angle
        self.batch   = batch

    def prepare_project(self):
        '''
        Read the maps and take ffts
        '''
        self.half_map1 = Map(self.half_map1_file)
        self.half_map2 = Map(self.half_map2_file)
        if self.whole_map_file is not None:
            self.whole_map = Map(self.whole_map_file)
        if self.mask_file is not None:
            self.mask      = Map(self.mask_file)

        # Get fibonacci points
        self.fib_points = parallelem.fibonacci_sphere(self.ncones)

        # Apply mask
        self.apply_mask()

        # Prepare ffts
        self.prepare_ffts()

        # Prepare fft grid
        self.prepare_fft_grid()

    def prepare_io_files(self):
        '''
        Prepare io files
        '''
        # Copy input files
        root, file = os.path.split(self.half_map1_file)
        copyfile(self.half_map1_file, self.output_directory+'/'+file)

        root, file = os.path.split(self.half_map2_file)
        copyfile(self.half_map2_file, self.output_directory+'/'+file)

        # Copy whole map
        if self.whole_map_file is not None:
            root, file = os.path.split(self.whole_map_file)
            copyfile(self.whole_map_file, self.output_directory+'/'+file)

        # Copy mask
        if self.mask_file is not None:
            root, file = os.path.split(self.mask_file)
            copyfile(self.mask_file, self.output_directory+'/'+file)

        # Fcc output file
        self.fcc_out_file = self.output_directory+'/dfsc.mrc'

    def prepare_ffts(self):
        '''
        Prepare ffts
        '''
        self.half_map1.take_fft()
        self.half_map1.center_fft()

        self.half_map2.take_fft()
        self.half_map2.center_fft()

    def prepare_fft_grid(self):
        '''
        Prepare fft grid
        '''
        # Get map dimensions
        ydim, xdim, zdim = self.half_map1.get_map_dimensions()

        xfreq   = np.fft.fftshift(np.fft.fftfreq(ydim, self.apix))
        yfreq   = np.fft.fftshift(np.fft.fftfreq(xdim, self.apix))
        zfreq   = np.fft.fftshift(np.fft.fftfreq(zdim, self.apix))

        self.fft_sx, self.fft_sy, self.fft_sz = np.meshgrid(xfreq, yfreq, zfreq)
        self.fft_s                = np.sqrt(self.fft_sx**2 + self.fft_sy**2 + self.fft_sz**2)

        # Determine r
        spacing = 1.0/(xdim*self.apix)
        self.fft_r = np.round(self.fft_s/spacing)

        self.fft_r = self.fft_r.astype(int)

    def apply_mask(self):
        '''
        Apply mask on the maps
        '''
        if self.mask is not None:
            if self.half_map1 is not None:
                self.half_map1.apply_map_mask(self.mask)
            if self.half_map2 is not None:
                self.half_map2.apply_map_mask(self.mask)
            if self.whole_map is not None:
                self.whole_map.apply_map_mask(self.mask)

    def init_results(self):
        '''
        Initialize the results
        '''
        self.fcc_img3D   = np.zeros(self.half_map1.get_fft().shape)
        self.count_img3D = np.zeros(self.half_map1.get_fft().shape)

    def run(self):
        '''
        Run the computation
        '''
        # Create a pool
        mp_pool = multiprocessing.Pool(multiprocessing.cpu_count())

        # Fcc results
        self.pool_results = []

        # Initialize the results
        self.init_results()

        # Get ffts
        half_map1_fft = self.half_map1.get_fft()
        half_map2_fft = self.half_map2.get_fft()

        for i in range(self.fib_points.shape[0]):

            # Determine cone point
            cone_point = self.fib_points[i,:]

            print('Calculating FSC for cone %d/%d %d/100' % (i+1,
                                                      self.ncones,
                                                      100.0*(i+1)/self.ncones))
            # Create a new process
            worker_result = mp_pool.apply_async(parallelem.calc_fcc, args=(half_map1_fft,
                                                                           half_map2_fft,
                                                                           self.fft_r,
                                                                           self.fft_s,
                                                                           self.fft_sx,
                                                                           self.fft_sy,
                                                                           self.fft_sz,
                                                                           cone_point,
                                                                           self.angle))

            self.pool_results.append(worker_result)

            if len(self.pool_results) == self.batch:
                self.process_results()

        # Process results
        self.process_results()

        # Normalize results
        self.normalize_results()

    def process_results(self):
        '''
        Process results
        '''

        # Iterate over the results
        for result in self.pool_results:
            fcc_img3D, count_img3D = result.get()

            self.fcc_img3D   += fcc_img3D
            self.count_img3D += count_img3D.astype(int)

        self.pool_results = []

    def center_results(self):
        '''
        Center results
        '''
        self.ffc_img3D = np.fft.fftshift(self.fcc_img3D)

    def normalize_results(self):
        # Normalize the results by number of cones
        # Get only the non-zero values
        valid = self.count_img3D > 0
        self.fcc_img3D[valid] /= self.count_img3D[valid]

    def write_output_files(self):
        '''
        Write output files
        '''
        self.fcc_map = Map()
        self.fcc_map.set_img3D(self.fcc_img3D)
        self.fcc_map.write_img3D(self.fcc_out_file, self.apix)


class Map:
    '''
    3D Map object
    '''
    def __init__(self, file=None):
        self.map_file = file
        self.map_mrc  = None
        self.img3D    = None
        self.fft3D    = None

        # Rad the file
        if self.map_file is not None:
            self.read_map(self.map_file)

    def write_img3D(self, fname, pixel_size):
        '''
        Write img3D
        '''
        mrc = mrcfile.new(fname, overwrite=True)
        mrc.set_data(self.img3D)
        mrc.voxel_size = pixel_size
        mrc.close()

    def take_fft(self):
        '''
        Take fft
        '''
        self.fft3D = np.fft.fftn(self.img3D)

    def center_fft(self):
        '''
        Center fft
        '''
        self.fft3D = np.fft.fftshift(self.fft3D)

    def get_fft(self):
        '''
        Get fft
        '''
        return self.fft3D

    def set_img3D(self, img3D):
        '''
        Set img3D
        '''
        self.img3D = np.array(img3D.copy(), dtype=np.float32)

    def read_map(self, file):
        '''
        Read map
        '''
        self.map_mrc = MRC(file)

        # Assign Img3D
        if self.map_mrc.img3D is not None:
            self.img3D = self.map_mrc.img3D

    def get_map_dimensions(self):
        '''
        Get 3D map dimensions
        '''
        if self.img3D is not None:
            return self.img3D.shape
        else:
            return None

    def apply_map_mask(self, mask):
        '''
        Apply mask
        '''
        if mask is not None and self.img3D.shape == mask.img3D.shape:
            self.img3D = self.img3D*mask.img3D


class ProjectFlip(Project):
    '''
    Particle flip project
    '''
    def __init__(self, name='EMParticleFlip'):
        super().__init__(name)

        self.flipped_mrc_file = None
        self.flipped_mrc      = None

        self.flipped_star_file = None
        self.flipped_star      = None

        self.flipped_results   = []

        # For relion
        self.tmp1_star_file = None
        self.tmp2_star_file = None

        self.tmp_flip_star_file   = None
        self.tmp_noflip_star_file = None

        self.tmp_flip_star   = None
        self.tmp_noflip_star = None

        self.combined_flip_star = None

        # Relion arguments
        self.relion_flip_args   = []

        # Flip and noflip stars for merge project
        self.flip_star   = None
        self.noflip_star = None

    def read_flip_star(self, flip_star_file):
        '''
        Read flip star file
        '''
        if os.path.isfile(flip_star_file):
            self.flip_star =  Star(flip_star_file)

    def read_noflip_star(self, noflip_star_file):
        '''
        Read noflip star file
        '''
        if os.path.isfile(noflip_star_file):
            self.noflip_star = Star(noflip_star_file)

    def prepare_merge_project(self, flip_star_file, noflip_star_file):
        '''
        prepare merge project
        '''
        combined_flip_star_file = os.path.relpath(os.path.abspath(self.output_directory+'/particle_combined.star'))
        self.read_flip_star(flip_star_file)
        self.read_noflip_star(noflip_star_file)

        # Make the symlinks to input folders
        self.make_symlink2parent(flip_star_file, 'flip_input')
        self.make_symlink2parent(noflip_star_file, 'noflip_input')

        # Create and set rlnIsFlip column
        self.flip_star.set_column(label='rlnIsFlip', value=1)
        self.noflip_star.set_column(label='rlnIsFlip', value=0)

        # Get first data block
        flip_data_block   = self.flip_star.get_data_block()
        noflip_data_block = self.noflip_star.get_data_block()

        # Create merge project
        combined_star = Star()
        combined_star.set_data_block(pd.concat([flip_data_block, noflip_data_block]))

        # Write the combined file
        combined_star.write(combined_flip_star_file)

    def prepare_output_files(self):
        # Copy input file to output directory
        if self.particle_star_file is not None:
            head, tail = os.path.split(self.particle_star_file)
            root, ext  = os.path.splitext(tail)

            # Prepare symlink to input folder
            self.make_symlink2parent(self.particle_star_file)

            self.flipped_mrc_file  = os.path.relpath(os.path.abspath(self.output_directory+'/flipped.mrcs'))
            self.flipped_star_file = os.path.relpath(os.path.abspath(self.output_directory+'/flipped.star'))

    def create_output_mrc(self):
        if self.flipped_mrc is None:
            # Determine shape parameters
            num_particles = self.particle_star.data_block.shape[0]
            NY, NX        = self.first_particle_mrc.img2D.shape

            # Create output MRC file
            self.flipped_mrc = MRC(file=self.flipped_mrc_file, shape=(num_particles, NY, NX))

    def create_output_star(self):
        if self.flipped_star is None:
            self.flipped_star = Star()
            self.flipped_star.copy(self.particle_star)

    def prepare_project(self):
        '''
        Prepare project
        '''
        self.read_first_particle_mrc()
        self.prepare_output_files()
        self.create_output_mrc()
        self.create_output_star()

    def set_relion_args(self):
        '''
        Set relion arguments
        '''
        self.relion_flip_args = [self.relion_image_handler_exe,
                                 '--i', self.tmp1_star_file,
                                 '--o', 'flip',
                                 '--flipX']

        self.relion_noflip_args = [self.relion_image_handler_exe,
                                   '--i', self.tmp2_star_file,
                                   '--o', 'noflip']

    def prepare_project_relion(self):
        '''
        Prepare relion project
        '''
        self.create_files_relion()
        self.set_relion_image_handler_exe()
        self.set_relion_args()
        self.split_star_relion()

    def write_results(self):
        '''
        Write results
        '''
        # Get number of particles
        num_ptcls = len(self.flipped_results)

        # If there is not ptcl to write, return only
        if num_ptcls == 0:
            return

        # Show status
        print('Writing  %d particles' % (num_ptcls))

        # Get all the data
        ptcl_list = [ptcl_index for ptcl_result, ptcl_index in self.flipped_results]
        ptcl_data = [ptcl_result.get() for ptcl_result, ptcl_index in self.flipped_results]

        # Write mrc file
        self.flipped_mrc.mrc_data.data[ptcl_list] = ptcl_data

        # Write star file
        new_image_names = []
        for ptcl_index in ptcl_list:
            new_image_names.append('%07d@%s' % (ptcl_index+1, self.flipped_mrc_file))

        self.flipped_star.data_block.loc[ptcl_list, 'rlnImageName'] = new_image_names

        # Reset the containers
        self.flipped_results = []

    def create_files_relion(self):
        '''
        Create tmp files
        '''
        if self.particle_star_file is not None:
            head, tail = os.path.split(self.particle_star_file)
            root, ext  = os.path.splitext(tail)

            self.tmp1_star_file = os.path.relpath(os.path.abspath(self.output_directory+'/particle_tmp1.star'))
            self.tmp2_star_file = os.path.relpath(os.path.abspath(self.output_directory+'/particle_tmp2.star'))

            self.tmp_flip_star_file   = os.path.relpath(os.path.abspath(self.output_directory+'/particle_tmp1_flip.star'))
            self.tmp_noflip_star_file = os.path.relpath(os.path.abspath(self.output_directory+'/particle_tmp2_noflip.star'))

            self.combined_flip_star_file = os.path.relpath(os.path.abspath(self.output_directory+'/particle_combined.star'))

    def merge_star_relion(self):
        '''
        Merge the tmp star files
        '''
        # Make star files
        star1 = Star(self.tmp_flip_star_file)
        star2 = Star(self.tmp_noflip_star_file)

        # For star1 update the geometry parameters to reflect the flip operation
        star1.flipX()

        # Get first data block
        star1_data_block = star1.get_data_block()
        star2_data_block = star2.get_data_block()

        # Data block list
        data_block_list = []
        if star1_data_block is not None:
            data_block_list.append(star1_data_block)
        if star2_data_block is not None:
            data_block_list.append(star2_data_block)

        combined_data_block = pd.concat(data_block_list)

        # Make a combined star file
        combined_star = Star()
        combined_star.set_data_block(combined_data_block)
        combined_star.write(self.combined_flip_star_file)

    def split_star_relion(self):
        '''
        Split star file
        '''
        # Get particle data block
        particle_data_block = self.particle_star.get_data_block()

        # Make the masks
        if 'rlnIsFlip' in particle_data_block:
            flip_mask   = particle_data_block['rlnIsFlip'] == 1
            noflip_mask = particle_data_block['rlnIsFlip'] == 0

            # Make star files
            star1 = Star()
            star1.set_data_block(particle_data_block.loc[flip_mask, :])
            star1.write(self.tmp1_star_file)

            star2 = Star()
            star2.set_data_block(particle_data_block.loc[noflip_mask, :])
            star2.write(self.tmp2_star_file)
        else:
            self.particle_star.add_column('rlnIsFlip')
            self.particle_star.write(self.tmp2_star_file)

    def flip_particles_relion(self):
        '''
        Flip particles using relion
        '''
        if len(self.relion_image_handler_exe) > 0:

            # Execute only if the file exists
            if os.path.isfile(self.tmp1_star_file):
                self.relion_flip_subprocess = subprocess.run(self.relion_flip_args,
                                                             universal_newlines=True)
            if os.path.isfile(self.tmp2_star_file):
                self.relion_noflip_subprocess = subprocess.run(self.relion_noflip_args,
                                                               universal_newlines=True)
        else:
            sys.exit('Relion image handler doesnt exist')

        # Merge the star files
        self.merge_star_relion()

    def flip_particles(self, batch_size=100):
        '''
        Flip particles
        '''
        particle_data_block = self.particle_star.get_data_block()

        # Create a pool
        mp_pool = multiprocessing.Pool(multiprocessing.cpu_count())

        # Initialize results list
        self.flipped_results = []

        # Get number of particles
        num_ptcls = particle_data_block.shape[0]

        # Iterate over all the particles
        for ptcl_index, ptcl_row in particle_data_block.iterrows():

            print('Flipping Particle %d/%d %d/100' % (ptcl_index+1,
                                                      num_ptcls,
                                                      100.0*(ptcl_index+1)/num_ptcls))

            # Create a new process
            worker_result = mp_pool.apply_async(parallelem.flipX_ptcl, args=(ptcl_row,))

            self.flipped_results.append([worker_result, ptcl_index])

            # Write results
            if len(self.flipped_results) == batch_size:
                self.write_results()

        # Complete writing the remainings
        self.write_results()

        # Flip star file parameters
        self.flipped_star.flipX()

        # Close the output files and write
        self.flipped_star.write(self.flipped_star_file)
        self.flipped_mrc.close()


class ProjectStack(Project):
    '''
    Create particle stack
    '''
    def __init__(self, name='EMStack'):
        super().__init__(name)
        self.stack_mrc_file  = None
        self.stack_star_file = None

        self.stack_mrc       = None
        self.stack_star      = None

        # Fft grid parameters
        self.fft_r           = None

        # For normalization
        self.background_mask = None

        # Original micrograph apix
        self.orig_apix       = None

    def prepare_star(self):
        '''
        Prepare star file
        '''
        if self.particle_star.has_label('rlnReferenceImage') and not self.particle_star.has_label('rlnImageName'):
            self.particle_star.rename_columns({'rlnReferenceImage':'rlnImageName'})

    def prepare_background_mask(self, clipbox=None):
        '''
        Prepare background maks
        '''
        if self.particle_diameter_A is not None:
            self.background_mask = util.circular_mask(self.first_particle_mrc.get_img2D().shape, radius=self.particle_radius_pix)

        if clipbox is not None and self.background_mask is not None:
            self.background_mask = parallelem.clip_img2D(self.background_mask, clipbox)

    def prepare_output_files(self):
        '''
        Create output files
        '''

        # Make symlink to input folder
        self.make_symlink2parent(self.particle_star_file)

        self.stack_star_file = os.path.relpath(os.path.abspath(self.output_directory+'/stack.star'))
        self.stack_mrc_file  = os.path.relpath(os.path.abspath(self.output_directory+'/stack.mrcs'))

    def prepare_project(self, highpass=None, lowpass=None, clipbox=None):
        '''
        Prepare meta objects using reference class avarages
        '''
        self.prepare_star()
        self.sort_star()
        self.read_particle_apix()
        self.set_particle_radius()
        self.read_first_particle_mrc()
        self.eval_fft_grid()
        self.make_filter_mask(highpass, lowpass)
        self.prepare_output_files()
        self.create_output_stack_mrc(clipbox)
        self.create_output_stack_star()
        self.prepare_background_mask(clipbox)

    def eval_fft_grid(self):
        '''
        Eval fft grid
        '''
        self.first_particle_mrc.eval_fft_grid(apix=self.particle_apix)
        self.fft_r = self.first_particle_mrc.get_fft_r()

    def make_filter_mask(self, highpass=None, lowpass=None, sigma=2):
        '''
        Make filter pass - centered at 0 frequency
        '''
        # Shift fft_r grid
        self.fft_r_shift = np.fft.fftshift(self.fft_r)

        lowpass_mask  = np.ones(self.fft_r_shift.shape, dtype=bool)
        highpass_mask = np.ones(self.fft_r_shift.shape, dtype=bool)

        # Highpass and lowpass filter masks
        if highpass:
            highpass_mask = self.fft_r > 1.0/highpass
        if lowpass:
            lowpass_mask  = self.fft_r < 1.0/lowpass

        # Get the combination of two masks
        self.fft_mask = np.logical_or(highpass_mask, lowpass_mask)

        # Apply gaussian filter to smooth the mask
        self.fft_mask = scipy.ndimage.filters.gaussian_filter(self.fft_mask, sigma)

    def create_output_stack_star(self):
        '''
        Create output star file
        '''
        self.stack_star = Star()
        self.stack_star.copy(self.particle_star)

    def create_output_stack_mrc(self, clipbox=None):
        '''
        Create output subtract mrc object and file
        '''
        # Create MRCS output file
        if self.stack_mrc is None:
            # Determine shape parameters
            num_particles = self.particle_star.data_block.shape[0]
            NY, NX        = self.first_particle_mrc.img2D.shape

            # Check the clip box size
            if clipbox is not None and clipbox < NY and clipbox < NX:
                NX = clipbox
                NY = clipbox

            # Create output MRC file
            self.stack_mrc = MRC(file=self.stack_mrc_file, shape=(num_particles, NY, NX))

    def write_results(self):
        '''
        Write results
        '''
        # Get number of particles
        num_ptcls = len(self.stack_results)

        # If number of particles is 0, then quit early
        if num_ptcls == 0:
            return

        # Show status
        print('Writing  %d particles' % (num_ptcls))

        # Get all the data
        ptcl_list = [ptcl_index for ptcl_result, ptcl_index in self.stack_results]
        ptcl_data = [ptcl_result.get() for ptcl_result, ptcl_index in self.stack_results]

        # Write mrc file
        self.stack_mrc.mrc_data.data[ptcl_list] = ptcl_data

        # Write star file
        new_image_names = []
        for ptcl_index in ptcl_list:
            new_image_names.append('%07d@%s' % (ptcl_index+1, self.stack_mrc_file))

        self.stack_star.data_block.loc[ptcl_list, 'rlnImageName'] = new_image_names

        # Reset the containers
        self.stack_results = []

    def set_orig_apix(self, apix):
        '''
        Set original apix
        '''
        self.orig_apix = apix

    def sort_star(self, colname='rlnClassDistribution', ascending=False):
        '''
        Sort based on class distribution
        '''
        if self.particle_star.has_label(colname):
            self.particle_star.sort(column=colname, ascending=ascending)


    def center_stack_star(self):
        '''
        Center stack star
        '''
        data_block = self.stack_star.get_data_block()

        # Get the origin
        fracOX, intOX = np.modf(data_block['rlnOriginX'].values)
        fracOY, intOY = np.modf(data_block['rlnOriginY'].values)

        # Leav only the fractions
        self.stack_star.data_block['rlnOriginX'] =  fracOX
        self.stack_star.data_block['rlnOriginY'] =  fracOY

        # Adjust coordinates
        fracCX, intCX = np.modf(intOX*self.particle_apix/self.orig_apix)
        fracCY, intCY = np.modf(intOY*self.particle_apix/self.orig_apix)

        self.stack_star.data_block['rlnCoordinateX'] += intCX
        self.stack_star.data_block['rlnCoordinateY'] += intCY

        # Add the leftovers to origins
        self.stack_star.data_block['rlnOriginX'] += fracCX*self.orig_apix/self.particle_apix
        self.stack_star.data_block['rlnOriginY'] += fracCY*self.orig_apix/self.particle_apix

    def create_stack(self, batch_size=100, transform=False, clipbox=None, recenter=False):
        '''
        Flip particles
        '''
        particle_data_block = self.stack_star.get_data_block()

        # Create a pool
        mp_pool = multiprocessing.Pool(multiprocessing.cpu_count())

        # Initialize results list
        self.stack_results = []

        # Get number of particles
        num_ptcls = particle_data_block.shape[0]

        # Iterate over all the particles
        for ptcl_index, ptcl_row in particle_data_block.iterrows():

            print('Writing Particle %d/%d %d/100' % (ptcl_index+1,
                                                     num_ptcls,
                                                     int(100.0*(ptcl_index+1)/num_ptcls)))

            # Create a new process
            worker_result = mp_pool.apply_async(parallelem.read_ptcl_mrc, args=(ptcl_row, transform, self.fft_mask, clipbox, self.background_mask, recenter))

            self.stack_results.append([worker_result, ptcl_index])

            # Write results
            if len(self.stack_results) == batch_size:
                self.write_results()

        # Complete writing the remainings
        self.write_results()

        # If transform is ON, reset the offsets
        if transform:
            self.stack_star.reset_offsets()

        # If recentter only fix the coordinates
        if recenter:
            self.center_stack_star()

        # Close the output files and write
        self.stack_star.write(self.stack_star_file)
        self.stack_mrc.close()


class ProjectSubtract2D(Project):
    '''
    Particle subtraction project
    '''
    def __init__(self, name='EMParticleSubtract2D'):
        super().__init__(name)

        # Instantenous class mrc
        self.class_mrc = None

        # Mask files and objects
        self.mask_align_mrc_file     = None
        self.mask_structure_mrc_file = None
        self.mask_subtract_mrc_file  = None

        self.mask_align_mrc     = None
        self.mask_structure_mrc = None
        self.mask_subtract_mrc  = None

        # Alignment references
        self.ref_align_star_file  = None
        self.ref_align_mrc_file   = None

        # Output files and objects
        self.subtracted_star      = None
        self.subtracted_star_file = None

        self.subtracted_mrc       = None
        self.subtracted_mrc_file  = None

        # Particle props
        self.particle_diameter_A  = None
        self.particle_radius_pix  = None

        # Circular and threshold masks
        self.circular_mask        = None
        self.threshold_mask       = None

        # Intensity statistics and masks
        self.background_mask     = None
        self.structure_mask      = None

        self.background_mean = None
        self.background_std  = None

        self.structure_mean  = None
        self.structure_std   = None

        # Subtraction results
        self.subtraction_results = []

        # Subtraction functions
        self.sub_funcs = {'subctf':  parallelem.subtract_class_ctf,
                          'cropctf': parallelem.crop_class_ctf,
                          'crop':    parallelem.crop_class}

        # Other mask parameters
        self.inner_diameter   = None
        self.inner_radius_pix = None
        self.inner_mask       = None

        # Clip size in pixels
        self.clip_box         = None

    def write_results(self):
        '''
        Write results
        '''
        # Get number of particles
        num_ptcls = len(self.subtraction_results)

        # If number of particles is 0, then quit early
        if num_ptcls == 0:
            return

        # Show status
        print('Writing  %d particles' % (num_ptcls))

        # Get all the data
        ptcl_list = [ptcl_index for ptcl_result, ptcl_index in self.subtraction_results]
        ptcl_data = [ptcl_result.get() for ptcl_result, ptcl_index in self.subtraction_results]

        # Write mrc file
        self.subtracted_mrc.mrc_data.data[ptcl_list] = ptcl_data

        # Write star file
        new_image_names = []
        for ptcl_index in ptcl_list:
            new_image_names.append('%07d@%s' % (ptcl_index+1, self.subtracted_mrc_file))

        self.subtracted_star.data_block.loc[ptcl_list, 'rlnImageName'] = new_image_names

        # Reset the containers
        self.subtraction_results = []

    def duplicate_imgs(self):
        '''
        Duplicate images for parallel processing
        '''
        class_img2D          = self.class_mrc.get_img2D().copy()
        mask_align_img2D     = self.mask_align_mrc.get_img2D().copy()
        mask_structure_img2D = self.mask_structure_mrc.get_img2D().copy()
        mask_subtract_img2D  = self.mask_subtract_mrc.get_img2D().copy()

        return class_img2D, mask_align_img2D, mask_structure_img2D, mask_subtract_img2D

    def create_threshold_mask(self, class_mrc, threshold_high=0.05, threshold_low=None):
        '''
        Create threshold mask
        '''
        if self.threshold_mask is None:
            self.threshold_mask = class_mrc.make_threshold_mask(threshold_high=threshold_high, threshold_low=threshold_low)

    def subtract_class_mrc(self, threshold_val=0.05, max_ptcl=None, batch_size=100, subtract_func='subctf', subtract_bg=False, norm_method='ccc', skip_to_firstpeak=True):
        '''
        Subtract class mrc file
        '''

        # Replace ReferenceImage with ImageName
        self.ref_class_star.rename_column('rlnReferenceImage', 'rlnImageName')

        # Get class data
        class_data_block = self.ref_class_star.get_data_block()

        # Get particle data
        particle_data_block = self.particle_star.get_data_block()

        # Particle counter
        particle_counter = 0

        # Create a pool
        mp_pool = multiprocessing.Pool(multiprocessing.cpu_count())

        # Results container
        self.subtraction_results = []

        # Iterate over each class
        for class_index, class_row in class_data_block.iterrows():

            # Get class info
            class_image_num, class_image_name = self.ref_class_star.get_image_num_name(class_index)

            # Read image
            self.class_mrc = MRC(class_image_name, int(class_image_num)-1)

            # Transform the class mrc with
            self.class_mrc.transform_ptcl_img2D(class_row)

            # Store to original
            self.class_mrc.store_to_original()

            # Prepare for CTF
            self.class_mrc.eval_ctf_grid(self.particle_apix)

            # Get class numbers for the current class
            current_class_number = class_row['rlnClassNumber']
            class_mask           = particle_data_block['rlnClassNumber'] == current_class_number
            particle_data        = particle_data_block.loc[class_mask, :]
            num_ptcls            = particle_data_block.shape[0]

            # Make threshold mask
            self.threshold_mask = self.class_mrc.make_threshold_mask(threshold_high=threshold_val, threshold_low=-threshold_val)

            # If stucture mask is not provided, assign to circular mask
            if self.mask_structure_mrc_file is None:
                self.mask_structure_mrc = MRC()
                self.mask_structure_mrc.set_img2D(self.circular_mask)
                self.mask_structure_mrc.store_to_original()

            # If subtraction mask is not defined, assign structure mask to it
            if self.mask_subtract_mrc_file is None:
                self.mask_subtract_mrc = MRC()
                self.mask_subtract_mrc.set_img2D(self.mask_structure_mrc.get_img2D())

                # If inner mask exists, use that to create subtract mask
                if self.inner_mask is not None:
                    self.mask_subtract_mrc.apply_mask(1-self.inner_mask)

                self.mask_subtract_mrc.store_to_original()

            for ptcl_index, ptcl_row in particle_data.iterrows():
                # Update particle counter
                particle_counter += 1

                # If exceed max number of particles, quit
                if max_ptcl is not None and particle_counter > max_ptcl:
                    break

                print('Subtracting Particle %d/%d %d/100' % (particle_counter,
                                                             num_ptcls,
                                                             100.0*particle_counter/num_ptcls))

                # Copy masks and images
                (class_img2D,
                 mask_align_img2D,
                 mask_structure_img2D,
                 mask_subtract_img2D) = self.duplicate_imgs()

                # Organize parallel processing input parameters
                pl_class_img2D        = class_img2D
                pl_ctf_grid_angle     = self.class_mrc.ctf_a
                pl_ctf_grid_s         = self.class_mrc.ctf_s
                pl_ctf_grid_r         = self.class_mrc.ctf_r
                pl_ptcl_star          = ptcl_row
                pl_subtract_bg        = subtract_bg

                # Create a new process
                worker_result = mp_pool.apply_async(self.sub_funcs[subtract_func], args=(pl_class_img2D,
                                                                                         pl_ctf_grid_angle,
                                                                                         pl_ctf_grid_s,
                                                                                         pl_ctf_grid_r,
                                                                                         pl_ptcl_star,
                                                                                         mask_align_img2D,
                                                                                         mask_structure_img2D,
                                                                                         mask_subtract_img2D,
                                                                                         pl_subtract_bg,
                                                                                         norm_method,
                                                                                         skip_to_firstpeak,
                                                                                         self.clip_box))

                self.subtraction_results.append([worker_result, ptcl_index])

                # Write results
                if len(self.subtraction_results) == batch_size:
                    self.write_results()

        # Complete writing the remainings
        self.write_results()

        # Close the output files and write
        self.subtracted_star.write(self.subtracted_star_file)
        self.subtracted_mrc.close()

    def prepare_output_files(self):
        '''
        Create output files
        '''

        self.subtracted_star_file = os.path.relpath(os.path.abspath(self.output_directory+'/subtracted.star'))
        self.subtracted_mrc_file  = os.path.relpath(os.path.abspath(self.output_directory+'/subtracted.mrcs'))

    def read_masks(self):
        '''
        Read masks
        '''
        # 1. Alignment mask - ideally a circular mask
        if self.mask_align_mrc_file is not None:
            self.mask_align_mrc = MRC(self.mask_align_mrc_file)
        else:
            self.mask_align_mrc = MRC()
            self.mask_align_mrc.set_img2D(self.circular_mask)

        # 2. Structure mask - mask that defines the boundaries of structure
        if self.mask_structure_mrc_file is not None:
            self.mask_structure_mrc = MRC(self.mask_structure_mrc_file)

        # 3. Subtract mask - mask used for subtraction
        if self.mask_subtract_mrc_file is not None:
            self.mask_subtract_mrc = MRC(self.mask_subtract_mrc_file)

        # Store the originals of masks
        self.mask_align_mrc.store_to_original()

    def set_inner_diameter(self, inner_diameter):
        '''
        Set inner diameter
        '''
        self.inner_diameter = inner_diameter

    def prepare_project(self):
        '''
        Prepare project
        '''
        self.set_particle_radius()
        self.prepare_output_files()
        self.prepare_meta_objects()
        self.read_masks()

    def create_output_subtract_star(self):
        '''
        Create output star file
        '''
        self.subtracted_star = Star()
        self.subtracted_star.copy(self.particle_star)

    def create_output_subtract_mrc(self):
        '''
        Create output subtract mrc object and file
        '''
        # Create MRCS output file
        if self.subtracted_mrc is None:
            # Determine shape parameters
            num_particles = self.particle_star.data_block.shape[0]
            NY, NX        = self.first_particle_mrc.img2D.shape

            if self.clip_box is not None and self.clip_box < NY and self.clip_box < NX:
                NX = self.clip_box
                NY = self.clip_box

            # Create output MRC file
            self.subtracted_mrc = MRC(file=self.subtracted_mrc_file, shape=(num_particles, NY, NX))

    def create_circular_mask(self, sigma=4):
        '''
        Create circular mask
        '''
        if self.particle_radius_pix is None:
            self.particle_radius_pix = int(np.min(self.first_particle_mrc.img2D.shape)//2)

        # Create circular mask
        if self.circular_mask is None:
            self.circular_mask = util.circular_mask(self.first_particle_mrc.img2D.shape,
                                                    center=None,
                                                    radius=self.particle_radius_pix)
            # Make the mask soft
            self.circular_mask = scipy.ndimage.filters.gaussian_filter(self.circular_mask, sigma)

    def create_inner_mask(self):
        '''
        Use inner diameter to create an outer mask
        '''
        if self.inner_diameter is not None:
            self.inner_radius_pix = int(0.5*self.inner_diameter/self.particle_apix)

            # Create inner mask
            self.inner_mask = util.circular_mask(self.first_particle_mrc.img2D.shape,
                                                 center=None,
                                                 radius=self.inner_radius_pix)

    def set_clip_box(self, clipbox):
        '''
        Set clip box
        '''
        self.clip_box = clipbox

    def prepare_meta_objects(self):
        '''
        Prepare meta objects using reference class avarages
        '''
        self.read_first_particle_mrc()
        self.create_output_subtract_mrc()
        self.create_output_subtract_star()
        self.create_circular_mask()
        self.create_inner_mask()


class ProjectImgName(Project):
    '''
    Intersection project
    '''
    def __init__(self, name='EMReplaceImgName'):
        super().__init__(name)

        self.particle_star        = None
        self.particle_star_file   = None

        self.ref_star             = None
        self.ref_star_file        = None

        self.replace_str          = None
        self.ref_img_name         = None

    def read_reference(self, file):
        '''
        Read ref star file
        '''
        self.ref_star_file = os.path.abspath(file)
        self.ref_star = Star(file)

    def add_suffix(self, suffix_list):
        '''
        Add suffix to image name
        '''
        if len(suffix_list) == 1:
            self.particle_star.data_block['rlnImageName'] = self.particle_star.data_block['rlnImageName'].str.replace('.mrcs',suffix_list[0]+'.mrcs')
        elif len(suffix_list) == 2:
            self.particle_star.data_block['rlnImageName'] = self.particle_star.data_block['rlnImageName'].str.replace(suffix_list[0]+'.mrcs',suffix_list[1]+'.mrcs')

    def remove_suffix(self, remove_suffix):
        '''
        Add suffix to image name
        '''
        if remove_suffix is not None:
            self.particle_star.data_block['rlnImageName'] = self.particle_star.data_block['rlnImageName'].str.replace(remove_suffix+'.mrcs','.mrcs')

    def get_reference_img(self):
        '''
        Get reference img path
        '''
        if not self.ref_star.has_label('rlnImageName'):
            sys.exit('Reference star doesnt have rlnImageName column')

        # Get image name
        img_num, img_name = self.ref_star.data_block['rlnImageName'][0].split('@')

        self.ref_img_name = img_name

        return img_name

    def replace_particle_img(self):
        '''
        Replace particle img name
        '''
        
        # Get reference image name
        self.get_reference_img()

        orig_img_names = self.particle_star.data_block['rlnImageName'].tolist()

        # Convert image names
        new_img_names = list(map(lambda x: x.split('@')[0]+'@'+self.ref_img_name, orig_img_names))

        # Assign new image names
        self.particle_star.data_block['rlnImageName'] = new_img_names



class ProjectGroup(Project):
    '''
    Intersection project
    '''
    def __init__(self, name='EMGroup'):
        super().__init__(name)

        self.particle_star        = None
        self.particle_star_file   = None

        self.micrograph_star      = None
        self.micrograph_star_file = None

        self.maxmics              = None
        self.thresholddef         = None
        self.thresholdint         = None

        # Column value to compare
        self.cmp_column           = 'rlnDefocusU'

    def set_params(self, maxmics=50, threshdef=100, threshint=0.1):
        '''
        Set grouping parameters
        '''
        self.maxmics      = maxmics
        self.thresholddef = threshdef
        self.thresholdint = threshint

    def group_micrographs(self):
        '''
        Make the groups
        '''

        # If the micrograph star is not compatible, skip
        if self.micrograph_star is None or not self.micrograph_star.has_label(self.cmp_column) or self.micrograph_star.has_label('rlnGroupName'):
            return

        # Sort groups based on defocus groups
        self.sort_micrographs(self.cmp_column)

        # Group names
        group_name_list    = []
        current_group_name = 1

        # Get the first defocus value
        previous_defocus = self.micrograph_star.data_block[self.cmp_column].values.tolist()[0]

        # Iterate through the list
        for mic_index, mic_row in self.micrograph_star.data_block.iterrows():
            if mic_row[self.cmp_column] - previous_defocus > self.thresholddef:
                # Set the new group name
                current_group_name += 1

                # Set the previous defocus
                previous_defocus = mic_row[self.cmp_column]

            # Add the group name
            group_name_list.append(str(current_group_name))

        # Set the group name in reference star
        self.micrograph_star.data_block['rlnGroupName'] = np.array(group_name_list)

    def group_particles(self):
        '''
        Group particles
        '''
        # If the micrograph star is not compatible, skip
        if self.particle_star is None or not self.particle_star.has_label(self.cmp_column):
            return

        # Sort groups based on defocus groups
        self.particle_star.sort(self.cmp_column)

        # Group names
        group_name_list    = []
        current_group_name = 1

        # Get the first defocus value
        previous_defocus = self.particle_star.data_block[self.cmp_column].values.tolist()[0]

        # Iterate through the list
        for ptcl_index, ptcl_row in self.particle_star.data_block.iterrows():
            if ptcl_row[self.cmp_column] - previous_defocus > self.thresholddef:
                # Set the new group name
                current_group_name += 1

                # Set the previous defocus
                previous_defocus = ptcl_row[self.cmp_column]

            # Add the group name
            group_name_list.append(str(current_group_name))

        # Set the group name in reference star
        self.particle_star.data_block['rlnGroupName'] = np.array(group_name_list)

        # Finally sort particles based on rlnImageName
        self.particle_star.sort('rlnImageName')

    def assign_groups(self):
        '''
        Assign groups to particles
        '''
        # If the star data doesnt have the right columns, skip
        if not self.micrograph_star.has_label('rlnGroupName') or self.particle_star is None:
            return

        # Set a new column for the particle star
        self.particle_star.set_column('rlnGroupName', 0)

        # Iterate over the class names
        for mic_index, mic_row in self.micrograph_star.data_block.iterrows():
            # Get micrograph name
            micrograph_name = mic_row['rlnMicrographName']

            # Micrograph maske
            mask = self.particle_star.data_block['rlnMicrographName'] == micrograph_name 

            # Set the new group name
            self.particle_star.data_block.loc[mask, 'rlnGroupName'] = mic_row['rlnGroupName']


class ProjectIntersect(Project):
    '''
    Intersection project
    '''
    def __init__(self, name='EMIntersect'):
        super().__init__(name)

        self.particle_star      = None
        self.particle_star_file = None

        self.stars = []
        self.files = []

    def set_star_files(self, star_files):
        '''
        Set star files
        '''
        self.files = star_files

    def read_particle_star_file(self):
        '''
        Read first star file
        '''
        if len(self.files) > 0 and os.path.isfile(self.files[0]):

            print('Reading first star file %s' % (self.files[0]))
            self.particle_star_file = self.files[0]
            self.particle_star      = Star(self.particle_star_file)

    def intersect_stars(self):
        '''
        Intersect star files
        '''
        # Comparison columns 
        cmp_columns = ['rlnMicrographName', 'rlnCoordinateX', 'rlnCoordinateY']

        if self.particle_star is not None:
            for i in range(1, len(self.files)):

                print('Reading star file:%d %s' % (i, self.files[i]))

                # Read new star file
                current_star = Star(self.files[i])

                # Intersect with the first star
                self.particle_star.intersect(current_star, cmp_columns)

    def run(self, star_files):
        '''
        Run project
        '''
        self.set_star_files(star_files)
        self.read_particle_star_file()
        self.intersect_stars()
        self.prepare_io_files_star()


class ProjectPlot(Project):
    '''
    Plotting project
    '''
    def __init__(self, name='EMPlot'):
        super().__init__(name)

        self.column_names = []
        self.column_pairs = []
        self.column_diffs = []

        self.ref_star_file = None
        self.ref_star      = None

    def read_reference(self, file):
        '''
        Read ref star file
        '''
        self.ref_star_file = os.path.abspath(file)
        self.ref_star = Star(file)

    def plot_diff(self, diff_pair, nbins, ref_star=None):
        '''
        Plot difference
        '''
        if len(diff_pair) != 2:
            return

        # Get column names
        column1, column2 = diff_pair
        if ref_star is not None:
            diff = self.particle_star.ref_norm_diff(ref_star, column1, column2)
        else:
            diff = self.particle_star.get_norm_diff(column1, column2)
        if diff is not None:
            py.hist(diff, density=True, bins=nbins)
            py.xlabel(column1+'-'+column2)

    def calc_orientation(self):
        '''
        Calculate orientation information
        '''
        if (self.particle_star.has_label('rlnAnglePsi') and
            self.particle_star.has_label('rlnAnglePsiPrior') and
            self.particle_star.has_label('rlnAngleTilt') and
            self.particle_star.has_label('rlnAngleTiltPrior')):

            self.direction, self.diff_psi, self.diff_tilt = self.particle_star.get_align_diff()

            self.metadata = pd.DataFrame(data={'direction': self.direction,
                                               'deltaPsi':  self.diff_psi,
                                               'deltaTilt': self.diff_tilt})
    def plot_orientation(self, nbins):
        '''
        Plot orientation
        '''
        if self.direction is not None:
            py.hist(self.direction, density=True, bins=[-1.25, -0.75, 0.75, 1.25])
            py.xlim(-1.5, 1.5)
            py.xlabel('Alignment')

    def plot_diff_angles(self, nbins):
        '''
        Plot diff_psi and diff_tilt
        '''
        if self.diff_psi is not None and self.diff_tilt is not None:
            py.hist2d(self.diff_psi, self.diff_tilt, bins=nbins, cmap='Blues', norm=LogNorm(), range=[[0, 180], [0, 180]])
            py.xlabel(r'$\Delta\psi$')
            py.ylabel(r'$\Delta\theta$')

    def plot_mic_paths(self):
        '''
        Plot mic path counter
        '''
        # Read the dates
        regex = r"\d+-\d+-\d+"
        dates = []
        for path in self.particle_star.mic_counter.keys():
            match = re.search(regex, path)
            dates.append(match.group(0))
        py.bar(dates, self.particle_star.mic_counter.values(), color='orange')

    def plot_hist(self, column_name, nbins):
        '''
        Plot histogram
        '''
        if self.particle_star.has_label(column_name):
            # Get data
            py.hist(self.particle_star.get_norm_data(column_name), density=True, bins=nbins)
            py.xlabel(column_name)

    def plot_scatter(self, column_pair):
        '''
        Plot histogram
        '''
        if len(column_pair) != 2:
            return

        # Get column names
        column1, column2 = column_pair
        if self.particle_star.has_label(column1) and self.particle_star.has_label(column2):
            py.plot(self.particle_star.get_norm_data(column1), self.particle_star.get_norm_data(column2),'*')
            py.xlabel(column1)
            py.ylabel(column2)

    def plot_hist2D(self, column_pair, nbins=20):
        '''
        Plot histogram
        '''
        if len(column_pair) != 2:
            return

        # Get column names
        column1, column2 = column_pair
        if self.particle_star.has_label(column1) and self.particle_star.has_label(column2):
            py.hist2d(self.particle_star.get_norm_data(column1), self.particle_star.get_norm_data(column2),bins=nbins**2, cmap=py.cm.Blues)
            py.xlabel(column1)
            py.ylabel(column2)

    def run_ref(self, column_pairs, column_diffs, nbins=20):
        '''
        Plot utility comparing ref and ptcl star
        '''

        # Get column diffs
        self.column_diffs = [pair.split(':') for pair in column_diffs]

        # Get number of diffs
        num_diffs = len(self.column_diffs)

        # Determine number of rows and columns
        num_rows = int(np.sqrt(num_diffs)) + 1

        # Create figure
        py.figure(figsize=(20,20))

        # Plot difference plots
        for i in range(num_diffs):
            py.subplot(num_rows, num_rows, i+1)
            self.plot_diff(self.column_diffs[i], nbins, self.ref_star)

        # Tight layout
        py.tight_layout()

    def set_tick_fontsize(self, size=8):
        '''
        Set tick fontsizes
        '''
        py.rc('xtick',labelsize=size)
        py.rc('ytick',labelsize=size)

    def run_ptcl(self, column_names, column_pairs, column_diffs, orientation, nbins=20):
        '''
        Plot utility for ptcl star
        '''
        self.column_names = column_names

        # Get column pairs
        self.column_pairs = [pair.split(':') for pair in column_pairs]

        # Get column diffs
        self.column_diffs = [pair.split(':') for pair in column_diffs]

        # Get singles number
        num_singles = len(self.column_names)

        # Get number of pairs
        num_pairs = len(self.column_pairs)

        # Get number of diffs
        num_diffs = len(self.column_diffs)

        # Determine number of rows and columns
        num_rows = int(np.sqrt(num_singles+num_pairs+num_diffs+3)) + 1

        # Create figure
        py.figure(figsize=(20,20))

        # Plot histograms for each column
        for i in range(num_singles):
            py.subplot(num_rows, num_rows, i+1)
            self.plot_hist(self.column_names[i], nbins)

        # Update num plots
        num_plots = num_singles

        # Plot scatter plots
        for i in range(num_pairs):
            py.subplot(num_rows, num_rows, num_plots+i+1)
            self.plot_hist2D(self.column_pairs[i], nbins)

        # Update num plots
        num_plots += num_pairs

        # Plot difference plots
        for i in range(num_diffs):
            py.subplot(num_rows, num_rows, num_plots+i+1)
            self.plot_diff(self.column_diffs[i], nbins)

        # Update num plots
        num_plots += num_diffs
        
        # Plot orientation angle
        if orientation:
            self.calc_orientation()
            py.subplot(num_rows, num_rows, num_plots+1)
            self.plot_orientation(nbins)

            py.subplot(num_rows, num_rows, num_plots+2)
            self.plot_diff_angles(nbins)

        # Plot micrographs paths
        if self.particle_star.mic_counter is not None:
            py.subplot(num_rows, num_rows, num_plots+3)
            self.plot_mic_paths()

        # Tight layout
        py.tight_layout()

    def prepare_io_files(self, output_format='svg'):
        '''
        Prepare output files
        '''
        
        head, tail = os.path.split(self.particle_star_file)
        root, ext  = os.path.splitext(tail)
        copyfile(self.particle_star_file, self.output_directory+'/particle_input'+ext)
        self.particle_plot_file = self.output_directory+'/particle_output.'+output_format

        # Make symlink
        self.make_symlink2parent(self.particle_star_file)

        # Prepare metadata file
        self.prepare_metadata_file()

    def write_output_files(self,output_format='svg'):
        '''
        Write output files
        '''
        # Write metadata
        self.write_metadata()

        # Save plot
        py.savefig(self.particle_plot_file, dpi=100,transparent=False, format=output_format)


class ProjectAlign2D(Project):
    '''
    Class Align2D Project Class
    '''
    def __init__(self, name='EMClassAlign2D'):
        super().__init__(name)

        # Temporary class star and mrc files
        self.ref_class_tmp_star_file      = None
        self.ref_class_tmp_star_norm_out  = None
        self.ref_class_tmp_star_norm_file = None

        # Alignment outfile
        self.ref_class_tmp_star_align_file = None
        self.ref_class_tmp_mrc_file        = None

        # Final aligned outfiles
        self.ref_class_final_out           = None

        # Relion output string
        self.relion_output_str    = None
        self.ref_class_number     = 1

        # Trasnformed mrc file
        self.ref_class_transformed_star_file = None
        self.ref_class_transformed_mrc_file  = None
        self.ref_class_transformed_star      = None
        self.ref_class_transformed_mrc       = None

        # Relion refine2D parameters
        self.relion_num_iter                 = 1

    def prepare_project(self, use_unmasked_classes=False):
        '''
        Prepare project
        '''
        self.set_particle_radius()
        self.create_tmp_files_star()
        self.prepare_tmp_input(use_unmasked_classes)
        self.write_tmp_files()
        self.set_relion_output_str()
        self.set_relion_refine_exe()
        self.set_relion_norm_exe()
        self.set_relion_stack_create_exe()
        self.read_first_ref_class_mrc(self.ref_class_number, write_ref=True)

    def set_params(self, apix, diameter):
        '''
        Set particle apix and diamter
        '''
        self.particle_apix       = apix
        self.particle_diameter_A = diameter

    def create_output_transformed_mrc(self):
        '''
        Create output subtract mrc object and file
        '''

        # Read first ref class mrc file
        self.read_first_ref_class_mrc(write_ref=False)

        # Determine the output transformed mrcs and star files
        self.ref_class_transformed_mrc_file  = os.path.relpath(os.path.abspath(self.output_directory+'/Class2D_output_transformed.mrcs'))
        self.ref_class_transformed_star_file = os.path.relpath(os.path.abspath(self.output_directory+'/Class2D_output_transformed.star'))

        # Create MRCS output file
        if self.ref_class_transformed_mrc is None:
            # Determine shape parameters
            num_particles = self.ref_class_star.data_block.shape[0]
            NY, NX        = self.first_ref_class_mrc.img2D.shape

            # Create output MRC file
            self.ref_class_transformed_mrc = MRC(file=self.ref_class_transformed_mrc_file, shape=(num_particles, NY, NX))

    def write_transformed_stack(self):
        '''
        Write transformed stack
        '''
        # Get class data
        class_data_block   = self.ref_class_star.get_data_block()

        # Keep transformed imgs
        transformed_img2Ds = []
        ptcl_list          = []

        # Iterate over each class
        for class_index, class_row in class_data_block.iterrows():

            # Get class info
            class_image_num, class_image_name = self.ref_class_star.get_image_num_name(class_index)

            # Read image
            class_mrc = MRC(class_image_name, int(class_image_num)-1)

            # Transform image
            class_mrc.transform_ptcl_img2D(class_row)

            # Add img to list
            transformed_img2Ds.append(class_mrc.img2D.copy())

            # ptcl list
            ptcl_list.append(class_index)

        # Write the transfomed img2D
        self.ref_class_transformed_mrc.mrc_data.data[ptcl_list] = transformed_img2Ds
        self.ref_class_transformed_mrc.close()

    def write_transformed_star(self):
        '''
        Write transformed star
        '''
        if self.ref_class_star.has_label('rlnImageName'):
            self.ref_class_star.data_block['rlnImageName'] = self.ref_class_star.data_block.rlnImageName.replace({r'@.*':'@'+self.ref_class_transformed_mrc_file},regex=True)

        # Reset the offsets and angles
        self.ref_class_star.data_block['rlnOriginX']  = 0.0
        self.ref_class_star.data_block['rlnOriginY']  = 0.0
        self.ref_class_star.data_block['rlnAnglePsi'] = 0.0

        self.ref_class_star.write(self.ref_class_transformed_star_file)

    def create_transformed_class_stacks(self):
        '''
        Create transformed class2D stacks
        '''
        self.create_output_transformed_mrc()
        self.write_transformed_stack()
        self.write_transformed_star()

    def normalize_class_refs(self):
        '''
        Normalize class references
        '''
        self.relion_norm_args = [self.relion_norm_exe,
                                 '--norm',
                                 '--operate_on',  self.ref_class_tmp_star_file,
                                 '--operate_out', self.ref_class_tmp_star_norm_out,
                                 '--bg_radius',   str(self.particle_radius_pix),
                                 ]
        self.relion_norm_subprocess = subprocess.run(self.relion_norm_args,
                                                     stdout=subprocess.PIPE,
                                                     universal_newlines=True)

    def prepare_tmp_input(self, use_unmasked_classes=False):
        '''
        Prepare the class reference star file for alignment of class averages
        '''
        if self.ref_class_star is not None:
            self.ref_class_star.change_label('rlnReferenceImage', 'rlnImageName')

            # Add psi prior column
            self.ref_class_star.set_column('rlnAnglePsiPrior',0)

            # Add offset prior columns
            self.ref_class_star.set_column('rlnOriginXPrior',0)

            # Add offset prior columns
            self.ref_class_star.set_column('rlnOriginYPrior',0)

            # If unmasked class option is on, use unmasked classes
            if use_unmasked_classes:
                self.ref_class_star.replace_with_unmasked_classes()

    def create_tmp_files_star(self):
        # Copy input file to output directory
        if self.ref_class_star_file is not None:
            head, tail = os.path.split(self.ref_class_star_file)
            root, ext  = os.path.splitext(tail)
            self.ref_class_tmp_star_file      = self.output_directory+'/'+root+'_class_tmp'+ext

            # Set normalization files
            self.ref_class_tmp_star_norm_out  = os.path.relpath(os.path.abspath(self.output_directory+'/'+root+'_class_tmp_norm'))
            self.ref_class_tmp_star_norm_file = os.path.relpath(os.path.abspath(self.output_directory+'/'+root+'_class_tmp_norm'+ext))

    def set_refine2D_files(self):
        '''
        Set Refine 2D files
        '''
        self.refine2D_it0_star_file = self.output_directory + '/' + 'run_it000_data.star'
        self.refine2D_it1_star_file = self.output_directory + '/' + 'run_it%03d_data.star' %(self.relion_num_iter)

    def read_refine2D_files(self):
        '''
        Read Refine 2D files
        '''
        self.refine2D_it0_star = Star(self.refine2D_it0_star_file)
        self.refine2D_it1_star = Star(self.refine2D_it1_star_file)

    def prepare_refine2D(self):
        # Set refine 2D files
        self.set_refine2D_files()

        # Read the files
        self.read_refine2D_files()

        # Copy rlnClassNumber from start file to final file
        self.refine2D_it1_star.copy_columns(self.refine2D_it0_star, {'rlnClassNumber': None})

        # Assign to class reference star
        self.ref_class_star = self.refine2D_it1_star

    def write_tmp_files(self):
        '''
        Write tmp files
        '''
        self.ref_class_star.write(self.ref_class_tmp_star_file)

    def set_relion_refine_args(self, skip_rotate=False, sigma_psi=-1, offset_range=10, offset_step=1, psi_step=1, gpu=0, num_refs=1, num_iter=1, firstiter_cc=True, T=2):

        # Check if reference file exists
        if self.ref_align_mrc_file is None:
            self.ref_align_mrc_file = self.first_ref_class_mrc_file

        # Set number of iterations
        self.relion_num_iter = num_iter

        # Get the maximum offset range possible
        if offset_range is None:
            offset_range_max = int(self.first_ref_class_mrc.img2D.shape[0]//2)
        else:
            offset_range_max = offset_range

        self.relion_refine_args = [self.relion_refine_exe,
                                   '--i', self.ref_class_tmp_star_norm_file,
                                   '--o', self.relion_output_str,
                                   '--dont_combine_weights_via_disc',
                                   '--flatten_solvent',
                                   '--zero_mask',
                                   '--oversampling', '1',
                                   '--norm',
                                   '--scale',
                                   '--offset_step',  str(offset_step),
                                   '--offset_range', str(offset_range_max),
                                   '--psi_step',     str(psi_step),
                                   '--j', '3',
                                   '--pool', '50',
                                   '--pad', '2',
                                   '--iter', str(num_iter),
                                   '--tau2_fudge', str(T),
                                   '--particle_diameter', str(self.particle_diameter_A),
                                   '--angpix', str(self.particle_apix),
                                   '--K', str(num_refs)]

        # First itercc
        if firstiter_cc:
            self.relion_refine_args.append('--firstiter_cc')

        # GPU option
        if gpu is not None:
            self.relion_refine_args+=['--gpu', str(gpu)]

        # Check rotation option for 2D class alignment
        if skip_rotate:
            self.relion_refine_args.append('--skip_rotate')

        # Sigma psi option to limit psi search
        if sigma_psi > 0:
            self.relion_refine_args += ['--sigma_psi', str(sigma_psi)]

        if self.ref_align_mrc_file is not None:
            self.relion_refine_args += ['--ref', self.ref_align_mrc_file]

        if self.mask_align_mrc_file is not None:
            self.relion_refine_args += ['--solvent_mask', self.mask_align_mrc_file]

    def run_refine2D(self):
        '''
        Run relion_refine to align classes
        '''

        if len(self.relion_refine_exe) > 0:
            self.relion_refine_subprocess = subprocess.run(self.relion_refine_args,
                                                           stdout=subprocess.PIPE,
                                                           universal_newlines=True)
        else:
            sys.exit('Relion refine doesnt exist')


class Micrograph:
    '''
        Micrograph class
    '''
    def __init__(self):
        self.name = None


class Class2D:
    '''
    2D class model
    '''
    def __init__(self):
        self.name = None
        self.star = None
        self.idx  = None


class Particle:
    '''
        Particle Class
    '''
    def __init__(self):
        self.name = None
        self.star = None
        self.idx  = None


class MRC:
    '''
        MRC  class
    '''
    def __init__(self, file=None, image_num=None, shape=None):
        self.name             = None
        self.mrc_data         = None
        self.project          = None
        self.header           = None
        self.img3D            = None
        self.img2D_original   = None

        self.img2D            = None
        self.img2D_fft        = None
        self.img2D_ifft       = None

        self.img2D_pshift     = None
        self.img2D_ctf        = None

        # CTF grid parameters
        self.ctf_s            = None
        self.ctf_sx           = None
        self.ctf_sy           = None
        self.ctf_a            = None
        self.ctf_r            = None
        self.star_data        = None
        self.star             = None
        self.apix             = None

        # FRC parameters
        self.frc2D            = None

        # Masks
        self.mask_align      = None
        self.mask_subtract   = None
        self.mask_circular   = None
        self.mask_threshold  = None
        self.mode2type = {0: np.dtype(np.int8),
                          1: np.dtype(np.int16),
                          2: np.dtype(np.float32),
                          6: np.dtype(np.uint16)}

        self.type2mode = {np.dtype(np.int8):    0,
                          np.dtype(np.int16):   1,
                          np.dtype(np.float32): 2,
                          np.dtype(np.uint16):  6}

        self.HEADER_LEN = int(1024)                     # In bytes.

        if file is not None:
            if os.path.isfile(file):
                self.read(file, image_num)
            else:
                self.create(file, shape)

    def gaussian_filter(self, sigma=2):
        '''
        Apply gaussian filter to img2D
        '''

    def ccc(self, other, mask=None):
        '''
        Measure ccc with other img2D
        '''
        current_mean, current_std = self.calc_mean_std_intensity(mask)
        other_mean, other_std     = other.calc_mean_std_intensity(mask)

        if mask is not None:
            cross_correlation = np.average((self.img2D-current_mean)*(other.img2D-other_mean), weights=mask)
        else:
            cross_correlation = np.mean((self.img2D-current_mean)(other.img2D-other_mean))

        return cross_correlation/(current_std*other_std)

    def flipX(self):
        '''
        Flip on X-axis
        '''
        self.img2D = self.img2D[:, ::-1]

    def intersect(self, other):
        '''
        Take the intersection with other
        '''
        if self.img2D.shape == other.img2D.shape:
            self.img2D = self.img2D*other.img2D

    def copy(self):
        '''
        Copy contents of other mrc to current one
        '''
        other = MRC()
        other.img2D_original = np.copy(self.img2D_original)
        other.img2D          = np.copy(self.img2D_original)
        other.ctf_s          = np.copy(self.ctf_s)
        other.ctf_sx         = np.copy(self.ctf_sx)
        other.ctf_sy         = np.copy(self.ctf_sy)
        other.ctf_a          = np.copy(self.ctf_a)
        other.ctf_r          = np.copy(self.ctf_r)

        return other

    def create_noise(self, mask=None, noise_mean=0.0, noise_std=1.0):
        '''
        Make Noise
        '''
        noise = np.random.normal(noise_mean, noise_std, self.img2D.shape)
        if mask is not None:
            self.img2D[mask > 0] = noise[mask > 0]
        else:
            self.img2D = noise

    def normalize_bg_area_intensity(self, mask_bg, new_val_bg, mask_area, new_val_area):
        '''
        Normalize the intensity - background and an area of interest
        '''
        # Get the background intensity
        background_intensity = np.mean(self.img2D[mask_bg > 0])

        # Subtract background intensity
        self.img2D -= background_intensity

        # Get the area intensity
        area_intensity = np.mean(self.img2D[mask_area > 0])

        # Normalize the area intensity
        self.img2D    *= (new_val_area-new_val_bg)/area_intensity

        # Finally add the new background intenstiy
        self.img2D    += new_val_bg

    def set_background_intensity(self, mask, new_val=0):
        '''
        Set background inetensity
        '''

        background_intensity = np.mean(self.img2D[mask > 0])

        # Subtract background intensity
        self.img2D -= background_intensity

        # Set new intensity
        self.img2D += new_val

    def set_area_intensity(self, mask, new_val=1.0):
        '''
        Set area intensity
        '''
        area_intensity = np.mean(self.img2D[mask > 0])
        self.img2D    *= new_val/area_intensity

    def convert_to_binary_mask(self):
        '''
        Convert img2D to binary mask
        '''
        mask = np.zeros(self.img2D.shape, dtype='float32')
        mask[self.img2D > 0] = 1
        self.img2D = np.copy(mask)

    def set_norm_intensity_params(self, mean, std):
        '''
        Set norm intensity params
        '''
        self.norm_mean = mean
        self.norm_std  = std

    def apply_mask(self, mask):
        '''
        Apply mask
        '''
        if mask.shape == self.img2D.shape:
            self.img2D = self.img2D*mask

    def make_circular_mask(self):
        '''
        Make circular Mask
        '''
        self.mask_circular = util.circular_mask(self.img2D.shape)
        return self.mask_circular

    def make_threshold_mask(self, threshold_high=0.05, threshold_low=None):
        ''''
        Make threshold mask
        '''
        self.mask_threshold = util.threshold_mask(self.img2D, threshold_high, threshold_low)
        return self.mask_threshold

    def eval_ctf_grid(self, apix=1.0):
        '''
        Create ctf freq grids
        '''
        if self.img2D is not None:
            assert self.img2D.shape[0] == self.img2D.shape[1]

            # Get the pixel size information
            if self.apix is not None:
                apix = self.apix

            xfreq   = np.fft.rfftfreq(self.img2D.shape[1], apix)
            yfreq   = np.fft.fftfreq(self.img2D.shape[0], apix)

            self.ctf_sx, self.ctf_sy  = np.meshgrid(xfreq, yfreq)
            self.ctf_s                = np.sqrt(self.ctf_sx**2 + self.ctf_sy**2)

            # Determine angle grids
            self.ctf_a = np.arctan2(self.ctf_sy, self.ctf_sx)

            # Determine r
            spacing = 1.0/(self.img2D.shape[0]*apix)
            self.ctf_r = np.round(self.ctf_s/spacing)

            # Determine r-cutoff
            r_cutoff = self.img2D.shape[0] // 2
            mask     = self.ctf_r >= r_cutoff

            self.ctf_r[mask] = r_cutoff 

            # Change ctf_r data dtype
            self.ctf_r = np.array(self.ctf_r, dtype=int)

    def get_ctf_r(self):
        '''
        Get ctf_r
        '''
        return self.ctf_r.copy()

    def eval_fft_grid(self, apix=1.0):
        '''
        Create fft freq grids
        '''
        if self.img2D is not None:
            assert self.img2D.shape[0] == self.img2D.shape[1]

            # Get the pixel size information
            if self.apix is not None:
                apix = self.apix

            xfreq   = np.fft.fftfreq(self.img2D.shape[1], apix)
            yfreq   = np.fft.fftfreq(self.img2D.shape[0], apix)

            self.fft_sx, self.fft_sy  = np.meshgrid(xfreq, yfreq)
            self.fft_s                = np.sqrt(self.fft_sx**2 + self.fft_sy**2)

            # Determine r
            spacing = 1.0/(self.img2D.shape[0]*apix)
            self.fft_r  = np.round(self.fft_s/spacing)*spacing

    def get_fft_r(self):
        '''
        Get fft_r
        '''
        return self.fft_r.copy()

    def rotate_img2D(self, angle):
        '''
        Rotate img2D by an eularangle
        '''
        self.img2D = scipy.ndimage.rotate(self.img2D, angle=-angle, axes=(0, 1), reshape=False)

    def shift_img2D(self, shiftX, shiftY):
        '''
        Shift img2D by a vector
        '''
        self.img2D = scipy.ndimage.shift(self.img2D, shift=[shiftY, shiftX])

    def inv_rotate_ptcl_img2D(self, ptcl_star):
        '''
        Inverse rotation based on ptcl data
        '''
        psi = ptcl_star['rlnAnglePsi']
        self.rotate_img2D(-psi)

    def inv_shift_ptcl_img2D(self, ptcl_star):
        '''
        Inverse shift based on ptcl data
        '''
        originX = ptcl_star['rlnOriginX']
        originY = ptcl_star['rlnOriginY']

        self.shift_img2D(-originX, -originY)

    def inv_transform_ptcl_img2D(self, ptcl_star):
        '''
        Inverse ptcl transform
        '''

        self.inv_rotate_ptcl_img2D(ptcl_star)
        self.inv_shift_ptcl_img2D(ptcl_star)

    def transform_ptcl_img2D(self, ptcl_star):
        self.shift_ptcl_img2D(ptcl_star)
        self.rotate_ptcl_img2D(ptcl_star)

    def rotate_ptcl_img2D(self, ptcl_star):
        '''
        Rotation based on ptcl data
        '''
        if 'rlnAnglePsi' not in ptcl_star:
            return

        psi = ptcl_star['rlnAnglePsi']
        self.rotate_img2D(psi)

    def shift_ptcl_img2D(self, ptcl_star):
        '''
        Shift based on ptcl data
        '''
        if 'rlnOriginX' not in ptcl_star or 'rlnOriginY' not in ptcl_star:
            return

        originX = ptcl_star['rlnOriginX']
        originY = ptcl_star['rlnOriginY']

        self.shift_img2D(originX, originY)

    def eval_ptcl_ctf(self, ptcl_star, bf=0):
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

        return self._eval_ctf(defU, defV, defA, phaseShift, kV, ac, cs, bf)

    def subtract_ctf(self, ctf):
        '''
        Subtract ctf from img2D
        '''
        self.fft_img2D()
        self.img2D_fft = self.img2D_fft/ctf
        self.ifft_img2D()
        self.copy_to_img2D(self.img2D_ifft)

    def correct_fft_ctf(self):
        '''
        Correct ctf
        '''
        if self.img2D_ctf is not None:
            self.img2D_fft *= self.img2D_ctf

    def correct_fft_pshift(self):
        '''
        Correct fft with pshift
        '''
        if self.img2D_pshift is not None:
            self.img2D_fft *= self.img2D_pshift

    def eval_ptcl_fft_pshift(self, ptcl_star):
        '''
        Determine fft-pshift
        '''
        originX = ptcl_star['rlnOriginX']
        originY = ptcl_star['rlnOriginY']

        return self._eval_fft_pshift(originX, originY)

    def _eval_fft_pshift(self, originx, originy):
        '''
        Evaluate pshift
        '''

        self.img2D_pshift = np.exp(-2 * np.pi * 1j * (-originx * self.ctf_sx + -originy * self.ctf_sy))
        return self.img2D_pshift

    def _eval_ctf(self, defU, defV, defA=0, phaseShift=0, kv=300, ac=0.1, cs=2.0, bf=0, lp=None, hp=None, do_intact_until_first_peak=False):
        '''
        :param defU: 1st prinicipal underfocus distance (Å).
        :param defV: 2nd principal underfocus distance (Å).
        :param defA: Angle of astigmatism (deg) from x-axis to azimuth.
        :param phaseShift: Phase shift (deg).
        :param kv:  Microscope acceleration potential (kV).
        :param ac:  Amplitude contrast in [0, 1.0].
        :param cs:  Spherical aberration (mm).
        :param bf:  B-factor, divided by 4 in exponential, lowpass positive.
        :param lp:  Hard low-pass filter (Å), should usually be Nyquist
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
        s = self.ctf_s

        s2 = s**2
        s4 = s2**2
        dZ = def_avg + def_dev * (np.cos(2 * (self.ctf_a - defA)))
        gamma = (k1 * dZ * s2) + (k2 * s4) - k5 - k3

        # Determine ctf
        img2D_ctf = -np.sin(gamma)

        # Do intact until first peak
        if do_intact_until_first_peak:

            # Mask for low angles
            low_mask = np.abs(gamma) < np.pi/2

            # Assign 1
            img2D_ctf[low_mask] = 1.0

        # Enforce envelope
        if bf != 0:
            img2D_ctf *= np.exp(-k4 * s2)

        return img2D_ctf

    def calc_frc(self, other):
        '''
        Compute frc between two ft
        '''
        self.frc2D = np.zeros(self.img2D_fft.shape, dtype=np.float32)
        rbins = np.sort(np.unique(self.ctf_r))
        self.frc1D = np.zeros(len(rbins), dtype=np.float32)
        for i in range(len(rbins)):
            mask = self.ctf_r == rbins[i]
            corr  = np.sum(self.img2D_fft[mask]*np.conj(other.img2D_fft[mask]))
            norm1 = np.sqrt(np.sum(np.abs(self.img2D_fft[mask])**2))
            norm2 = np.sqrt(np.sum(np.abs(other.img2D_fft[mask])**2))

            self.frc1D[i]    = np.real(corr)/(norm1*norm2)
            self.frc2D[mask] = self.frc1D[i]

    def normalize_frc(self):
        '''
        Normalize frc
        '''

        if self.frc2D is not None:
            self.img2D_fft *= self.frc2D

    def calc_mean_std_intensity(self, mask=None):
        '''
        Calculate mean and std intensity
        '''
        if mask is not None and mask.shape == self.img2D.shape:
            self.mean_intensity = np.mean(self.img2D[mask > 0])
            self.std_intensity  = np.std(self.img2D[mask > 0])
        else:
            self.mean_intensity = np.mean(self.img2D)
            self.std_intensity  = np.std(self.img2D)

        return self.mean_intensity, self.std_intensity

    def normalize_intensity(self, mask=None, new_mean=0, new_std=None):
        '''
        Normalize intensity to match a new gauss-distribution
        '''
        self.calc_mean_std_intensity(mask)

        # Adjust the mean and std-values value
        # Zero mean value
        self.img2D -= self.mean_intensity

        # Adjust stdev
        if new_std is not None:
            self.img2D *= new_std/self.std_intensity

        # Bring mean value to new mean-value
        self.img2D += new_mean

    def get_img2D(self):
        '''
        Get img2D
        '''
        return self.img2D

    def fft_img2D(self, mask=None):
        '''
        FFT img2D
        '''
        if mask is not None and mask.shape == self.img2D.shape:
            self.img2D_fft = np.fft.rfft2(mask*self.img2D)
        else:
            self.img2D_fft = np.fft.rfft2(self.img2D)

    def ifft_img2D(self):
        '''
        Inverse FFT 2D
        '''
        self.img2D_ifft = np.real(np.fft.irfft2(self.img2D_fft))
        return self.img2D_ifft

    def set_img2D_fft(self, data):
        '''
        Set img2D_fft
        '''
        self.img2D_fft = np.copy(data)

    def set_img2D(self, data):
        '''
        Set img2D
        '''
        self.img2D = np.copy(data)

    def set_img3D(self, data):
        '''
        Set img3D
        '''
        self.img3D = np.copy(data)

    def set_apix(self, apix):
        '''
        Set pixel size
        '''
        self.apix = apix

    def set_star(self, star):
        '''
        Set star for the mrc file
        '''
        self.star = star

    def set_star_data(self, star_data):
        '''
        Set star data
        '''
        self.star_data = star_data

    def create(self, file, shape=None):
        '''
        Create file
        '''
        if shape is not None:
            self.mrc_data = mrcfile.new_mmap(file, shape, mrc_mode=2)

    def append_img(self, data, img_counter):
        '''
        Append img
        '''
        if self.mrc_data is not None:
            self.mrc_data.data[img_counter] = data

    def close(self):
        '''
        Close file
        '''
        self.mrc_data.close()

    def flush(self):
        '''
        Flush data to file
        '''
        self.mrc_data.flush()

    def read(self, file, image_num=None):
        '''
        Read MRC file
        '''

        with mrcfile.mmap(file, permissive=True, mode='r') as self.mrc_data:
            self.header   = self.mrc_data.header
            self.img3D    = self.mrc_data.data

            # If image is 2D
            if len(self.img3D.shape) == 2:
                self.img2D = np.copy(self.img3D)

            # If image is 3D with a single image
            if len(self.img3D.shape) == 3 and self.img3D.shape[0] == 1:
                self.img2D = np.copy(self.img3D[0])

            # If image num is defined and it is 3D data
            if image_num is not None and image_num < self.img3D.shape[0] and len(self.img3D.shape) == 3:
                self.img2D = np.copy(self.img3D[image_num])

            # Store an original copy of img2D
            if self.img2D is not None:
                self.img2D_original = np.copy(self.img2D)

    def store_from_original(self):
        '''
        Store from original
        '''
        self.img2D     = np.copy(self.img2D_original)

    def store_to_original(self):
        '''
        Store to original
        '''
        self.img2D_original = np.copy(self.img2D)

    def copy_to_original(self, data):
        '''
        Copy to original
        '''
        self.img2D_original = np.copy(data)

    def copy_to_img2D(self, data):
        '''
        Copy to img2D
        '''
        if data.shape == self.img2D.shape:
            self.img2D = np.copy(data)

    def subtract_from_img2D(self, data, mask=None):
        '''
        Subtract from img2D
        '''
        if mask is not None and data.shape == self.img2D.shape and mask.shape == self.img2D.shape:
            self.img2D -= data*mask
        else:
            self.img2D -= data

    def copy_img(self, other):
        '''
        Deep copy of the other image
        '''

        # Set mrc data
        self.mrc_data.set_data(other.img)

        # Set the img
        self.img3D = self.mrc_data.data

    def set_img(self, img):
        '''
        Set image data
        '''
        self.mrc_data.set_data(img)

    def write_img(self, fname, apix=1, origin=None, fast=False):
        """
        Write a MRC file. Fortran axes order is assumed.
        :param fname: Destination path.
        :param apix: Pixel size in Å for MRC header.
        :param origin: Coordinate of origin voxel.
        :param fast: Skip computing density statistics in header. Default is False.
        """

        # Create an mrc data file
        self.mrc_data = mrcfile.new(fname, overwrite=True)

        # Set mrc data
        self.mrc_data.set_data(self.img3D)

        # Update header stats and header
        self.mrc_data.update_header_from_data()

        # Set origin
        if origin is not None:
            self.mrc_data.header['origin'] = origin

        # Update stats in no fast mode
        if not fast:
            self.mrc_data.update_header_stats()

        # Set pixel size
        self.mrc_data.voxel_size = apix

        # Close the file
        self.mrc_data.close()

    def write_imgs(self, fname, idx):
        with mrcfile.mmap(fname, mode='w+') as mrc:
            # Mrc data shape
            mrc_nz, mrc_ny, mrc_nx = mrc.data.shape

            # Img shape
            img_nz, img_ny, img_nx = self.img3D.shape

            # Check the two data
            if img_ny == mrc_ny and img_nx == mrc_nx and mrc_nz >= idx+img_nz:
                mrc.data[idx:idx+img_nz] = self.img3D

    def read_imgs(self, fname, idx, num=1):
        with mrcfile.mmap(fname, mode='w+') as mrc:
            # Mrc data shape
            mrc_nz, mrc_ny, mrc_nx = mrc.data.shape

            # Check the two data
            if mrc_nz >= idx+num:
                self.img3D = mrc.data[idx:idx+num]

    def get_img(self, img_num=0):
        '''
        Get a single image
        '''
        if self.img3D is not None and img_num < self.img3D.shape[0]:
            return self.img3D[img_num]
        else:
            return None

    def determine_com(self, img_num=0, threshold_val=0):
        '''
        Determine center-of-mass at a given image number
        '''
        if self.img3D is not None and img_num < self.img3D.shape[0]:
            # Get normalized image
            self.img2D = self.img3D[img_num]

        if self.img2D is not None:
            # Create a mask
            mask = np.array(self.img2D > threshold_val, dtype='float32')
            origin_x = int(0.5*self.img2D.shape[1])
            origin_y = int(0.5*self.img2D.shape[0])
            x, y = np.meshgrid(np.arange(self.img2D.shape[1], dtype=np.double),
                               np.arange(self.img2D.shape[0], dtype=np.double))

            com_x = np.sum(x*mask*self.img2D)/np.sum(mask*self.img2D)
            com_y = np.sum(y*mask*self.img2D)/np.sum(mask*self.img2D)

            return [com_x-origin_x, com_y-origin_y]
        else:
            return None

    def read_header(self, file):
        '''
         HEADER FORMAT
        # 0      (0,4)      NX  number of columns (fastest changing in map)
        # 1      (4,8)      NY  number of rows
        # 2      (8,12)     NZ  number of sections (slowest changing in map)
        # 3      (12,16)    MODE  data type:
        #                       0   image: signed 8-bit bytes range -128 to 127
        #                       1   image: 16-bit halfwords
        #                       2   image: 32-bit reals
        #                       3   transform: complex 16-bit integers
        #                       4   transform: complex 32-bit reals
        # 4      (16,20)    NXSTART number of first column in map
        # 5      (20,24)    NYSTART number of first row in map
        # 6      (24,28)    NZSTART number of first section in map
        # 7      (28,32)    MX      number of intervals along X
        # 8      (32,36)    MY      number of intervals along Y
        # 9      (36,40)    MZ      number of intervals along Z
        # 10-13  (40,52)    CELLA   cell dimensions in angstroms
        # 13-16  (52,64)    CELLB   cell angles in degrees
        # 16     (64,68)    MAPC    axis corresp to cols (1,2,3 for X,Y,Z)
        # 17     (68,72)    MAPR    axis corresp to rows (1,2,3 for X,Y,Z)
        # 18     (72,76)    MAPS    axis corresp to sections (1,2,3 for X,Y,Z)
        # 19     (76,80)    DMIN    minimum density value
        # 20     (80,84)    DMAX    maximum density value
        # 21     (84,88)    DMEAN   mean density value
        # 22     (88,92)    ISPG    space group number, 0 for images or 1 for volumes
        # 23     (92,96)    NSYMBT  number of bytes in extended header
        # 24-49  (96,196)   EXTRA   extra space used for anything
        #           26  (104)   EXTTYP      extended header type("MRCO" for MRC)
        #           27  (108)   NVERSION    MRC format version (20140)
        # 49-52  (196,208)  ORIGIN  origin in X,Y,Z used for transforms
        # 52     (208,212)  MAP     character string 'MAP ' to identify file type
        # 53     (212,216)  MACHST  machine stamp
        # 54     (216,220)  RMS     rms deviation of map from mean density
        # 55     (220,224)  NLABL   number of labels being used
        # 56-256 (224,1024) LABEL(80,10)    10 80-character text labels
        '''

        with mrcfile.open(file, header_only=True) as self.mrc_data:
            self.header   = self.mrc_data.header
        return self.header


class EMfile:
    '''
        EM file class
    '''
    def __init__(self):
        self.name = None


class Star(EMfile):
    '''
        Star class
    '''
    def __init__(self, file=None):
        self._file_path     = None                        # Path to directory where this script stays
        self.star_file      = None
        self.name           = None
        self.data_block     = None
        self.data_props     = None
        self.data_name      = None
        self.data_labels    = None
        self.data_formats   = None
        self.data_dtypes    = None
        self.data_skip_rows = None
        self.metadata_file  = 'relion_metadata_labels.dat'
        self.PARAMETERS     = {}
        self.str2type       = {'double': float,
                               'string': str,
                               'int':    int,
                               'bool':   lambda x: bool(int(x))}
        self.str2nptype     = {'double': 'f4',
                               'string': 'U100',
                               'int': 'i4',
                               'bool': 'b'}

        self.type2format    = {'double': "%13.6f",
                               'string': "%s",
                               'int':    "%12d",
                               'bool':   "%2d"}

        # Read metadata labels
        self._read_metadata_labels()

        # Star pixel size
        self.star_apix     = None

        # Micrograph pixel size
        self.mic_apix      = 1.82

        # Micrographs
        self.micrographs   = None
        self.mic_folders   = None
        self.mic_counter   = None

        # Read file
        if file is not None and os.path.isfile(file):
            self.read(file)

    def read_mrc_paths(self):
        '''
        Read mrc paths
        '''
        if self.has_label('rlnMicrographName'):
            self.micrographs = self.data_block['rlnMicrographName'].tolist()
            self.mic_folders = [os.path.split(os.path.realpath(mic))[0] for mic in self.micrographs]
            self.mic_counter = Counter(self.mic_folders)

    def get_align_diff(self):
        '''
        Get alignment difference in angles
        '''
        # Get psi difference
        diffPsi   = self.get_norm_diff('rlnAnglePsi', 'rlnAnglePsiPrior')

        # Get tilt difference
        diffTilt  = self.get_norm_diff('rlnAngleTilt', 'rlnAngleTiltPrior')

        # Cumulative difference
        directionPsi  = 2*(diffPsi < 90).astype(int)-1
        directionTilt = 2*(diffTilt < 90).astype(int)-1 
        
        return directionPsi*directionTilt, diffPsi, diffTilt

    def get_norm_diff(self, column1, column2):
        '''
        Get angle difference
        '''
        diff = None
        if self.has_label(column1) and self.has_label(column2):
            diff = self.get_norm_data(column1) - self.get_norm_data(column2)
            
            # If the columns are angles, perform a normalization procedure
            if self.is_angle_column(column1) and self.is_angle_column(column2):
                # If the columns are angles perform a different analysis
                diff = util.get_min_angle_diff(diff)

        return diff

    def ref_norm_diff(self, other, column1, column2):
        '''
        Get angle differences between the columns of self and other star
        '''
        diff = None
        if self.has_label(column1) and other.has_label(column2):
            diff = self.get_norm_data(column1) - other.get_norm_data(column2)
            
            # If the columns are angles, perform a normalization procedure
            if self.is_angle_column(column1) and other.is_angle_column(column2):
                # If the columns are angles perform a different analysis
                diff = util.get_min_angle_diff(diff)

        return diff

    def is_angle_column(self, column):
        '''
        Check if it is an angle column
        '''
        return ('Angle' in column)


    def align_to_prior(self, column):
        '''
        Align angle values to its prior values
        '''

        # Get prior name for the angle
        column_prior = column+'Prior'
        if self.has_label(column) and self.has_label(column_prior) and self.is_angle_column(column):
            # Angle difference
            diff = (self.data_block[column]%360 - self.data_block[column_prior]%360)%360
            
            # Complement angle difference
            diffComp = 360 - diff
            
            # Get minimum difference
            diffMin = pd.DataFrame([diff, diffComp]).min()
            
            # If minimum difference is over 90, flip it
            self.data_block.loc[diffMin > 90, column] += 180

            # Normalize angles
            self.data_block[column] %= 360

    def pick_random_set(self, num_ptcls, rand_seed=1):
        '''
        Pick random set of ptcls
        '''
        total_ptcls = len(self.data_block.index)

        if num_ptcls < total_ptcls:
            self.data_block = self.data_block.sample(n=num_ptcls, random_state=rand_seed)

    def get_norm_data(self, label):
        '''
        Get normalized data
        '''
        data = None
        if self.has_label(label):
            data = self.data_block[label]

        if data is not None and 'Angle' in label:
            data %= 360

        return data

    def create_shortImageName(self):
        '''
        Create short image name column from imagename
        '''
        new_column = []
        for ptcl_index, ptcl_row in self.data_block.iterrows():
            # Parse imagename
            image_id, image_name = ptcl_row['rlnImageName'].split('@')
            head, tail = os.path.split(image_name)

            # New image name
            new_image_name = str(int(image_id))+'@'+tail

            # Append image to list
            new_column.append(new_image_name)

        # Assign new column to short imagename
        self.data_block['shortImageName'] = new_column

        # Add new parameter
        self.PARAMETERS['shortImageName'] = self.PARAMETERS['rlnImageName']

    def delete_shortImageName(self):
        '''
        Delete short image name
        '''
        self.delete_column('shortImageName')

    def set_barcode(self, barcode={}):
        '''
        Set the particle barcode
        '''
        barcode_str = util.convert_dict2str(barcode)

        # Set barcode
        self.set_column('rlnParticleName', barcode_str)

    def set_ptcl_barcode(self, ptcl_index, barcode={}):
        '''
        Set ptcl barcode
        '''

        if ptcl_index < self.data_block.shape[0] and self.has_label('rlnParticleName'):
            current_barcode = self.read_ptcl_barcode(ptcl_index)
            new_barcode     = {**current_barcode, **barcode}

            self.data_block.loc[ptcl_index,'rlnParticleName'] = util.convert_dict2str(new_barcode)

    def append_barcode(self, barcode={}):
        '''
        Append the particle barcode
        '''

        if not self.has_label('rlnParticleName'):
            self.set_barcode(barcode)
        else:
            barcode_str_list = []
            for ptcl_index, ptcl_row in self.data_block.iterrows():
                # Get new barcode
                new_barcode = self.read_ptcl_barcode(ptcl_index, barcode)

                # Get barcode information for the particle
                barcode_str_list.append(util.convert_dict2str(new_barcode))

            self.data_block['rlnParticleName'] = barcode_str_list

    def read_ptcl_barcode(self, ptcl_index, barcode={}):
        '''
        Read barcode
        '''
        current_barcode = {}
        new_barcode     = {}

        if self.has_label('rlnParticleName'):
            barcode_str = self.data_block.loc[ptcl_index, 'rlnParticleName']
            current_barcode = util.parse_star_parameters(barcode_str.split(','))

        # Update new barcode
        new_barcode = {**current_barcode, **barcode}

        return new_barcode

    def delete_by_barcode(self, select_barcode_list):
        '''
        Select ptcls by barcode
        '''

        if not self.has_label('rlnParticleName'):
            return

        selected_ptcls = []

        for ptcl_index, ptcl_row in self.data_block.iterrows():
            ptcl_barcode_list = ptcl_row['rlnParticleName'].strip().split(',')

            # Get intersection elements
            intersect = [barcode for barcode in select_barcode_list if barcode in ptcl_barcode_list]

            if len(intersect) != len(select_barcode_list):
                selected_ptcls.append(ptcl_index)

        # Make list numpy array
        selected_ptcls = np.array(selected_ptcls)

        self.data_block = self.data_block.loc[selected_ptcls, :]

    def select_by_barcode(self, select_barcode_list):
        '''
        Select ptcls by barcode
        '''

        if not self.has_label('rlnParticleName'):
            return

        selected_ptcls = []

        for ptcl_index, ptcl_row in self.data_block.iterrows():
            ptcl_barcode_list = ptcl_row['rlnParticleName'].strip().split(',')

            # Get intersection elements
            intersect = [barcode for barcode in select_barcode_list if barcode in ptcl_barcode_list]

            if len(intersect) == len(select_barcode_list):
                selected_ptcls.append(ptcl_index)

        # Make list numpy array
        selected_ptcls = np.array(selected_ptcls)

        self.data_block = self.data_block.loc[selected_ptcls, :]

    def intersect(self, other, cmp_columns):
        '''
        Intersect this star with other star object
        '''

        intersect_data_block = pd.merge(self.data_block, other.data_block[cmp_columns], how='inner')
        self.set_data_block(intersect_data_block)

    def merge(self, other, cmp_columns):
        '''
        Merge the two data set
        '''
        intersect_data_block = pd.merge(self.data_block, other.data_block,on=cmp_columns, how='inner')
        self.set_data_block(intersect_data_block)

    def sort(self, column='rlnDefocusU', ascending=False):
        '''
        Sort the star object based on a column
        '''
        if self.has_label(column):
            self.data_block.sort_values(column, inplace=True, ascending=ascending)
            self.data_block.reset_index(drop=True, inplace=True)

    def filter_classification(self, maxprob=0.5, maxclass=10):
        '''
        Filter the data using max probability in star data
        '''
        if self.has_label('rlnMaxValueProbDistribution'):
            if maxprob is not None:
                prob_mask = self.data_block['rlnMaxValueProbDistribution'] >= maxprob
                self.data_block = self.data_block.loc[prob_mask, :]

        if self.has_label('rlnNrOfSignificantSamples'):
            if maxclass is not None:
                class_mask = self.data_block['rlnNrOfSignificantSamples'] <= maxclass
                self.data_block = self.data_block.loc[class_mask, :]

    def filter_orientation(self, tilt_range=[0, 360], dpsi_range=[0, 360], dtilt_range=[0, 360], dalign_range=[0, 360]):
        
        # Tilt Angle
        if self.has_label('rlnAngleTilt'):
            tilt = self.get_norm_data('rlnAngleTilt')
            
            # Make selection
            self.data_block = self.data_block.loc[tilt>=tilt_range[0],:]
            self.data_block = self.data_block.loc[tilt<=tilt_range[1],:]

        # diffPsi
        if self.has_label('rlnAnglePsi') and self.has_label('rlnAnglePsiPrior'):
            diffPsi = self.get_norm_diff('rlnAnglePsi', 'rlnAnglePsiPrior')
            
            # Make selection
            self.data_block = self.data_block.loc[diffPsi>=dpsi_range[0],:]
            self.data_block = self.data_block.loc[diffPsi<=dpsi_range[1],:]

        # diffTilt
        if self.has_label('rlnAngleTilt') and self.has_label('rlnAngleTiltPrior'):
            diffTilt = self.get_norm_diff('rlnAngleTilt', 'rlnAngleTiltPrior')

            # Make selection
            self.data_block = self.data_block.loc[diffTilt>=dtilt_range[0],:]
            self.data_block = self.data_block.loc[diffTilt<=dtilt_range[1],:]

        # diffAlign
        if(self.has_label('rlnAnglePsi') and self.has_label('rlnAnglePsiPrior') and
           self.has_label('rlnAngleTilt') and self.has_label('rlnAngleTiltPrior')):
            
            diffAlign  = self.get_align_diff()

            # Make selection
            self.data_block = self.data_block.loc[diffAlign>=dalign_range[0],:]
            self.data_block = self.data_block.loc[diffAlign<=dalign_range[1],:]

    def set_data_block(self, data):
        '''
        Set data block
        '''
        self.data_block = data.copy()

    def replace_with_unmasked_classes(self):
        '''
        Replace with unmasked class averages
        '''
        for ptcl_index, ptcl_row in self.data_block.iterrows():
            image_num, image_name = ptcl_row['rlnImageName'].split('@')
            head, tail = os.path.split(image_name)

            # Construct new image name
            new_image_name = head + '/run_unmasked_classes.mrcs'

            # Replace the image name
            self.data_block.loc[ptcl_index, 'rlnImageName'] = image_num+'@'+new_image_name

    def get_image_num_name(self, ptcl):
        '''
        Get ptcl image num and name
        '''
        if ptcl < self.data_block.shape[0] and self.has_label('rlnImageName'):
            image_num, image_name = self.data_block.loc[ptcl, 'rlnImageName'].split('@')

        return int(image_num), image_name

    def rename_columns(self, column_params):
        if column_params is not None:
            for old_column, new_column in column_params.items():

                # Check that the column name exists and the new name is a proper star varaible
                if self.has_label(old_column) and new_column in self.PARAMETERS:
                    self.data_block = self.data_block.rename(index=int, columns={old_column:new_column})

    def reset_offsets(self):
        '''
        Reset all the offsets and the classification angles
        '''
        offset_params = ['rlnOriginX',
                         'rlnOriginY',
                         'rlnAnglePsi',
                         'rlnAngleRot',
                         'rlnAngleTilt']

        prior_params = ['rlnOriginXPrior',
                        'rlnOriginYPrior',
                        'rlnAnglePsiPrior',
                        'rlnAngleRotPrior',
                        'rlnAngleTiltPrior']

        # Set offsets to 0
        for param in offset_params:
            self.set_column(param, 0)

        # Delete the prior columns
        for param in prior_params:
            self.delete_column(param)

    def flipX(self):
        '''
        Modify the geometrical values for flip around X
        '''

        if self.has_label('rlnIsFlip'):
            valid_rows = self.data_block['rlnIsFlip'] == 1
        elif self.data_block is not None:
            valid_rows = np.arange(self.data_block.shape[0])
            self.add_column('rlnIsFlip')
        else:
            valid_rows = None

        # Invert X
        if self.has_label('rlnOriginX'):
            self.data_block.loc[valid_rows,'rlnOriginX'] = -self.data_block.loc[valid_rows,'rlnOriginX']

        # Invert X-Prior
        if self.has_label('rlnOriginXPrior'):
            self.data_block.loc[valid_rows,'rlnOriginXPrior'] = -self.data_block.loc[valid_rows,'rlnOriginXPrior']

        # Update Psi
        if self.has_label('rlnAnglePsi'):
            self.data_block.loc[valid_rows,'rlnAnglePsi'] = 180 - self.data_block.loc[valid_rows,'rlnAnglePsi']

        # Update Psi-Prior
        if self.has_label('rlnAnglePsiPrior'):
            self.data_block.loc[valid_rows,'rlnAnglePsiPrior'] = 180 - self.data_block.loc[valid_rows,'rlnAnglePsiPrior']

    def create_micname_from_imagename(self, mic_path='Micrographs'):
        '''
        Create micrographname from imagename
        '''

        if self.has_label('rlnImageName'):
            # Add micropraph name
            self.add_column('rlnMicrographName')

            # Create micrograph names
            new_mic_name_list = []
            for i in range(self.data_block.shape[0]):
                # Parse imagename
                image_head, image_tail = os.path.split(self.data_block['rlnImageName'][i])
                file_head, file_ext    = os.path.splitext(image_tail)

                # New micgropraph name
                new_mic_name = mic_path+'/'+file_head+'.mrc'

                # Add micropgraph names
                new_mic_name_list.append(new_mic_name)

            self.data_block.loc[:, 'rlnMicrographName'] = new_mic_name_list

    def build_proximity_map(self, pix_range=50, pix_step=50):
        '''
        Build ptcl map based on location on micrographs
        '''
        self.ptcl_map = {}

        for index, ptcl in self.data_block.iterrows():

            coord_x  = ptcl['rlnCoordinateX']
            coord_y  = ptcl['rlnCoordinateY']
            mic_name = ptcl['rlnMicrographName']

            x_floor = np.floor(1.0*(coord_x-pix_range)/pix_step)*pix_step
            x_ceil  = np.ceil(1.0*(coord_x+pix_range)/pix_step)*pix_step

            y_floor = np.floor(1.0*(coord_y-pix_range)/pix_step)*pix_step
            y_ceil  = np.ceil(1.0*(coord_y+pix_range)/pix_step)*pix_step

            x_list   = np.arange(x_floor, x_ceil+1, pix_step, dtype=int)
            y_list   = np.arange(y_floor, y_ceil+1, pix_step, dtype=int)

            xv, yv   = np.meshgrid(x_list, y_list)
            mic_list = [mic_name]*xv.size

            xym = list(zip(mic_list, xv.flatten(), yv.flatten()))
            ptcl_list = [ptcl]*len(xym)

            new_map = dict(zip(xym, ptcl_list))

            self.ptcl_map.update(new_map)

    def build_ptcl_map(self):
        '''
        Build ptcl map based on img id
        '''
        self.ptcl_map = {}

        for index, ptcl in self.data_block.iterrows():

            # Get ptcl key
           ptcl_key = self.get_ptcl_key(ptcl)

           self.ptcl_map[ptcl_key] = ptcl


    def set_micrograph_apix(self, apix=1.82):
        '''
        Set micrograph pixel size
        '''
        self.mic_apix = apix

    def _read_metadata_labels(self):
        '''
        Read the metadata labels
        '''
        self._file_path  = os.path.dirname(os.path.abspath(__file__))

        f = open(self._file_path + '/' + self.metadata_file, 'r')
        lines = f.readlines()
        f.close()

        # Parse through each line
        for line in lines:
            m = re.search('(\w+) \((\w+)\)\s+: (.*)', line)
            if m:
                param_name  = m.group(1)
                type_name   = m.group(2)
                description = m.group(3)
                self.PARAMETERS[param_name] = {'typename': type_name,
                                               'nptype': self.str2nptype[type_name],
                                               'type': self.str2type[type_name],
                                               'description': description}

    def _read_header(self, file):
        '''
        Read star header file
        '''
        f = open(file, 'r')

        # Assign star file
        self.star_file = file

        # Data block found
        data_block_found    = False
        data_labels_found   = False
        self.data_skip_rows = 0

        for line in f:
            if not data_block_found:
                m = re.search('^data_(\w*)', line)
                if m:
                    self.data_name = m.group(1)
                    data_block_found = True
            elif not data_labels_found:
                m = re.search('^loop_(\w*)', line)
                if m:
                    self.data_labels  = []
                    data_labels_found = True
            else:
                m = re.search('^_(\w+)', line)
                if m:
                    self.data_labels.append(m.group(1))
                else:
                    break

            # Update the number of rows that need to be skipped
            self.data_skip_rows += 1

        # Close the file
        f.close()

    def _prepare_data_types(self):
        '''
        Prepare the data types from header
        '''
        if self.data_block is not None:
            self.data_labels   = list(self.data_block.columns.values)
        self.data_formats  = [self.PARAMETERS[label]['nptype'] for label in self.data_labels]
        self.data_dtypes   = {'names': self.data_labels,
                              'formats': self.data_formats}

    def get_numRows(self):
        '''
        Get number of raws
        '''
        return self.data_block.shape[0]

    def delete_ptcls(self, ptcls):
        '''
        Delete particles
        '''
        prev_numRows = self.data_block.shape[0]
        self.data_block = self.data_block.drop(ptcls)
        new_numRows  = self.data_block.shape[0]

        print('Old row number: %d - New row number: %d - Deleted rows: %d' % (prev_numRows, new_numRows, len(ptcls)))

    def read(self, file):
        '''
        Read Star file and create the data blocks
        '''
        self._read_header(file)
        self._prepare_data_types()
        self.data_block = np.loadtxt(file,
                                     skiprows=self.data_skip_rows,
                                     dtype=self.data_dtypes)

        # Convert to data frame
        self.data_block = pd.DataFrame(self.data_block)

        # Assign number of data points
        self.num_data_points = self.data_block.shape[0]

    def delete_column(self, label):
        if self.has_label(label):
            self.data_block = self.data_block.drop(columns=label)

    def delete_columns(self, columns):
        '''
        Delete columns
        '''
        for label in columns:
            self.delete_column(label)

    def keep_columns(self, columns):
        '''
        Keep columns
        '''
        column_list = list(filter(self.has_label, columns))

        self.data_block = self.data_block[column_list]

    def add_column(self, label=None):
        '''
        Add new column to data block
        '''
        if label not in self.PARAMETERS:
            print('%s is not a valid Star label.' % (label))
            return None
        elif label in self.data_block.columns.values:
            print('%s exists. Not creating a new column.' % (label))
            return None

        # Create new column
        new_data_column = np.empty([self.data_block.shape[0], 1], dtype=self.PARAMETERS[label]['nptype'])

        # Append the data column
        self.data_block[label] = new_data_column

        # Initialize the column
        if self.PARAMETERS[label]['nptype'][0] == 'U':
            self.data_block.loc[:, label] = ''
        else:
            self.data_block.loc[:, label] = 0

        # Recreate data types
        self._prepare_data_types()

        return 1

    def set_column(self, label, value=None):
        '''
        Set a column value
        '''

        if self.has_label(label) and value is not None:
            self.data_block.loc[:, label] = self.PARAMETERS[label]['type'](value)
        else:
            success = self.add_column(label)
            if success and value is not None:
                self.data_block.loc[:, label] = self.PARAMETERS[label]['type'](value)

    def copy(self, other=None):
        '''
        Deep copy other star to self
        '''
        self.data_block   = other.data_block.copy()
        self.data_labels  = other.data_labels.copy()
        self.data_formats = other.data_formats.copy()
        self.data_dtypes  = other.data_dtypes.copy()

    def copy_columns(self, other, columns={}, ptcl_copy_list=None):
        '''
        Copy columns from another star object
        '''

        # If the sizes don't match don't perform copy
        if ptcl_copy_list is None and other.data_block.shape[0] != self.data_block.shape[0]:
            return 0
        elif ptcl_copy_list is None and other.data_block.shape[0] == self.data_block.shape[0]:
            ptcl_list = np.arange(self.data_block.shape[0])
            ptcl_copy_list = pd.DataFrame({'self': ptcl_list, 'other': ptcl_list})

        # Iterate over all columns
        for label, value in columns.items():
            if other.has_label(label):
                if self.has_label(label):
                    self.data_block.loc[ptcl_copy_list['self'].tolist(), label] = other.data_block.loc[ptcl_copy_list['other'].tolist(), label].tolist()
                else:
                    self.add_column(label)
                    self.set_column(label, value)

        # Assign only the portion of the data frame
        self.data_block = self.data_block.loc[ptcl_copy_list['self'].tolist(), :]

    def copy_column2column(self, from_column, to_column):
        '''
        Copy from one column to another
        '''
        if not self.has_label(from_column):
            return 0

        if not self.has_label(to_column):
            self.add_column(to_column)

        # Copy from-column to to-column
        self.data_block[to_column] = self.data_block[from_column]

    def _has_coordinate(self):
        '''
        Has coordinate
        '''
        if(self.has_label('rlnMicrographName') and
           self.has_label('rlnCoordinateX') and
           self.has_label('rlnCoordinateY')):

            return True
        else:
            return False

    def _get_coordinates(self, ptcls):
        '''
        Get coordinate information for a particle list
        '''
        if self._has_coordinate():

            # Get the coordinates
            coordinate_x = self.data_block['rlnCoordinateX'][ptcls]*self.mic_apix
            coordinate_y = self.data_block['rlnCoordinateY'][ptcls]*self.mic_apix

            if self.has_label('rlnOriginX') and self.has_label('rlnOriginY'):
                coordinate_x = coordinate_x - self.data_block['rlnOriginX'][ptcls]*self.star_apix
                coordinate_y = coordinate_y - self.data_block['rlnOriginY'][ptcls]*self.star_apix

            return np.hstack((np.vstack(coordinate_x.tolist()), np.vstack(coordinate_y.tolist())))
        else:
            return None

    def set_comment(self, ptcl, comment=''):
        '''
        Set comment
        '''
        if self.has_label('rlnParticleName'):
            self.data_block.loc[ptcl, 'rlnParticleName'] = comment

    def get_comment(self, ptcl):
        '''
        Read the comment
        '''
        if self.has_label('rlnParticleName'):
            return self.data_block['rlnParticleName'][ptcl]
        else:
            return None

    def _get_micrograph_name(self, ptcl):
        '''
        Get the micrograph name
        '''
        if self.has_label('rlnMicrographName'):
            return self.data_block['rlnMicrographName'][ptcl]
        else:
            return None

    def _get_ptcls_on_micrograph(self, micrograph_name=None):
        '''
        Get ptcl ids with a micrograph name
        '''
        ptcl_list = None
        if self.has_label('rlnMicrographName') and micrograph_name is not None:
            ptcl_list = self.data_block.index[self.data_block['rlnMicrographName'] == micrograph_name].tolist()

            return ptcl_list

    def get_proximity_key(self, ptcl, pixel_step=10):
        '''
        Get proximity key
        '''
        coordinate_x = ptcl["rlnCoordinateX"]
        coordinate_y = ptcl["rlnCoordinateY"]

        mic_name = ptcl["rlnMicrographName"]

        int_x = int(np.floor(1.0*coordinate_x/pixel_step)*pixel_step)
        int_y = int(np.floor(1.0*coordinate_y/pixel_step)*pixel_step)

        return (mic_name, int_x, int_y)

    def expand_img_column(self):
        '''
        Expand img column to its id-number and stack-name
        '''
        if self.has_label('rlnImageName'):
            def parse_img_name(img):
                img_id, img_name   = img.split('@')
                img_head, img_tail = os.path.split(img_name)

                return int(img_id), img_tail

            img_num_list  = []
            img_head_list = []

            img_data = np.array(list(map(parse_img_name, self.data_block['rlnImageName'].tolist())))

            # Expand the columns
            self.data_block['idx'] = img_data[:,0]
            self.data_block['img'] = img_data[:,1]

    def delete_img_columns(self):
        '''
        Delete expanded image columns
        '''
        if self.has_label('idx') and self.has_label('img'):
            self.data_block = self.data_block.drop(columns=['idx', 'img'])

    def get_ptcl_key(self, ptcl):
        '''
        Get ptcl key
        '''
        # Get img number and full image name
        img_num, img_name = ptcl['rlnImageName'].split('@')
        
        # Get img tail
        img_head, img_tail = os.path.split(img_name)

        # Construct ptcl key
        ptcl_key = (int(img_num), img_tail)

        return ptcl_key

    def has_label(self, label):
        '''
        Check if the label exists in data frame
        '''
        if self.data_block is not None and label in self.data_block.columns:
            return True
        else:
            return False

    def get_data_block(self):
        '''
        Get data block
        '''
        return self.data_block

    def is_particle_inside(self, ptcl, mic_apix, NX, NY):
        '''
        Is particle inside
        '''
        # Relative scale of pixel sizes
        apix_scale  = 1.0*self.star_apix/mic_apix
        cor_offsetx = self.data_block['rlnOriginX'][ptcl]*apix_scale
        int_offsetx = np.round(cor_offsetx)

        cor_offsety = self.data_block['rlnOriginY'][ptcl]*apix_scale
        int_offsety = np.round(cor_offsety)

        new_coordx = self.data_block.loc[ptcl, 'rlnCoordinateX'] - int_offsetx
        new_coordy = self.data_block.loc[ptcl, 'rlnCoordinateY'] - int_offsety

        if(new_coordx < NX and new_coordx > 0 and
           new_coordy < NY and new_coordy > 0):
            return True
        else:
            return False

    def addOffset2D(self, t=[0, 0], ptcls=None):
        '''
        Translate
        '''
        if len(t) == 2:
            dx = float(t[0])
            dy = float(t[1])

            # If necessary create the new data columns
            if not self.has_label('rlnOriginX'):
                self.add_column('rlnOriginX')
            if not self.has_label('rlnOriginY'):
                self.add_column('rlnOriginY')

            if ptcls is None:
                ptcls = np.arange(self.num_data_points)

            self.data_block.loc[ptcls, 'rlnOriginX'] += dx
            self.data_block.loc[ptcls, 'rlnOriginY'] += dy

    def set_star_apix(self, apix=None):
        '''
        Set star apix
        '''
        if type(apix) == float:
            self.star_apix = apix
        else:
            self.star_apix = 1.0

    def get_star_apix(self):
        '''
        Get star apix
        '''
        return self.star_apix

    def determine_star_apix(self):
        '''
        Determine star apix
        '''
        if self.has_label('rlnDetectorPixelSize') and self.has_label('rlnMagnification'):
            self.star_apix = 10000*self.data_block.loc[0, 'rlnDetectorPixelSize']/self.data_block.loc[0, 'rlnMagnification']
        else:
            print('Warning: No pixel size information in star file %s' % (self.star_file))
            self.star_apix = 1.0

        return self.star_apix

    def recenter2D(self, mic_apix=1.82):
        '''
        Recenter particles
        '''

        # Relative scale of pixel sizes
        self.apix_scale = 1.0*self.star_apix/mic_apix

        if(self.has_label('rlnOriginX') and
           self.has_label('rlnOriginY')):

            # Center x-coordinate
            cor_offsetx = self.data_block['rlnOriginX']*self.apix_scale
            int_offsetx = np.round(cor_offsetx)
            dif_offsetx = cor_offsetx - int_offsetx
            self.data_block.loc[:, 'rlnOriginX']      = dif_offsetx/self.apix_scale
            self.data_block.loc[:, 'rlnCoordinateX'] -= int_offsetx

            # Center y-coordinate
            cor_offsety = self.data_block['rlnOriginY']*self.apix_scale
            int_offsety = np.round(cor_offsety)
            dif_offsety = cor_offsety - int_offsety
            self.data_block.loc[:, 'rlnOriginY']      = dif_offsety/self.apix_scale
            self.data_block.loc[:, 'rlnCoordinateY'] -= int_offsety

    def change_label(self, old_label, new_label):
        '''
        Change label name
        '''
        if self.has_label(old_label) and new_label in self.PARAMETERS and not self.has_label(new_label):
            self.data_block.rename(columns={old_label: new_label},
                                   inplace=True)

    def rename_column(self, old_label, new_label):
        '''
        Rename column
        '''
        self.change_label(old_label, new_label)

    def dublicate_column(self, label, new_label):
        '''
        Duplicate a column with a label
        '''
        if self.has_label(label) and new_label in self.PARAMETERS:
            self.data_block.loc[:, new_label] = self.data_block[label]

    def rotate_psi(self, rotangle=0):
        '''
        Rotate psi angle
        '''
        self.data_block.loc[:, 'rlnAnglePsi'] += rotangle

        # Normalize psi
        self.normalize_psi()

    def merge_star(self, other_star):
        '''
        Merge with the current star
        '''
        if other_star is not None and self.data_block is not None:
            self.data_block = self.data_block.append(other_star.data_block, ignore_index=True)

    def normalize_psi(self):
        '''
        Normalize psi angle
        '''
        self.data_block.loc[:, 'rlnAnglePsi'] %= 360

        # Find angles higher than 180
        mask = self.data_block['rlnAnglePsi'] > 180

        # Subtract 180 from angles higher than 180
        self.data_block.loc[mask, 'rlnAnglePsi'] -= 360

    def rotate2D(self, rotangle=0, offset=[0, 0], final_offset=[0, 0], ptcls=None):
        '''
        Rotate particles
        '''

        # Check if the offset columns exist
        if(not self.has_label('rlnOriginX') or
           not self.has_label('rlnOriginY')):
            self.add_column('rlnOriginX')
            self.add_column('rlnOriginY')

        # Check if the Psi (rot for 2D transformation)
        if not self.has_label('rlnAnglePsi'):
            self.add_column('rlnAnglePsi')

        # Update Offsets
        if ptcls is None:
            ptcls = np.arange(self.num_data_points)

        # Check if there is any particle to transform
        if len(ptcls) == 0:
            return

        # Iterate through each particle to get the corrected offset
        new_offsets = []

        for ptcl in ptcls:
            oldangle = self.data_block.loc[ptcl, 'rlnAnglePsi']
            rotM = util.euler2rot2D(float(oldangle))

            # Get the transpose
            rotMT = rotM.T

            # Get the corrected offset
            corrected_offset = rotMT.dot(np.array(offset))

            # Final offset
            final_rotM  = util.euler2rot2D(float(oldangle+rotangle))
            final_rotMT = final_rotM.T

            final_corrected_offset = final_rotMT.dot(np.array(final_offset))

            new_offsets.append(corrected_offset+final_corrected_offset)

        # Update offsets (Needs to be investigated)
        new_offsets = np.array(new_offsets)
        self.data_block.loc[ptcls, 'rlnOriginX'] += new_offsets[:, 0]
        self.data_block.loc[ptcls, 'rlnOriginY'] += new_offsets[:, 1]

        # Update psi angles
        self.data_block.loc[ptcls, 'rlnAnglePsi'] += rotangle

        # Normalize psi angle
        self.normalize_psi()

    def num2className(self, ptcls=None):
        '''
        Assign class names from particle numbers in ImageName
        '''
        if not self.has_label('rlnClassNumber'):
            self.add_column('rlnClassNumber')

        # Get the particle ids
        particle_nums = [int(image_name.split('@')[0]) for image_name in self.data_block['rlnImageName']]
        particle_nums = np.array(particle_nums)

        if ptcls is None:
            ptcls = np.arange(self.num_data_points)

        self.data_block.loc[ptcls, 'rlnClassNumber'] = particle_nums[ptcls]

    def get_class_rows(self, class_id=1):
        '''
        Get rows with the defined class id
        '''
        if self.has_label('rlnClassNumber'):
            return np.nonzero(self.data_block['rlnClassNumber'] == class_id)[0]
        else:
            return None

    def get_class_ids(self):
        '''
        Return class ids
        '''
        return np.unique(self.data_block['rlnClassNumber'])

    def get_column(self, label=None):
        '''
        Get data column
        '''
        if self.has_label(label):
            return self.data_block[label]

    def tilt90(self):
        '''
        Tilt 90 star
        '''
        if self.has_label('rlnAngleTilt') and self.has_label('rlnAngleRot'):
            valid_rows = (self.data_block['rlnAngleTilt'] % 360) > 180

            # Update tilt angle
            self.data_block.loc[valid_rows, 'rlnAngleTilt'] -= 180
            self.data_block.loc[valid_rows, 'rlnAngleTilt'] %= 360

            # Update rot angle
            self.data_block.loc[valid_rows, 'rlnAngleRot'] += 180
            self.data_block.loc[valid_rows, 'rlnAngleRot'] %= 360

    def Zflip(self):
        '''
        Z-flip star
        '''
        if self.has_label('rlnAngleRot'):

            # Update tilt angle
            self.data_block['rlnAngleRot'] *= -1
            self.data_block['rlnAngleRot'] %= 360

    def create_write_formatter(self):
        '''
        Create write formatter
        '''
        formatter = []

        for label in self.data_block.columns:
            type_name = self.PARAMETERS[label]["typename"]
            formatter.append(self.type2format[type_name])

        # Create write formatter
        self.write_formatter = '  '.join(formatter)

    def write(self, out_fname, verbose=True):
        '''
        Write star file
        '''

        # Create the formatter
        self.create_write_formatter()

        # Create header
        header = []

        # Write data block name
        header.append("data_%s" % self.data_name)
        header.append("")
        header.append("loop_")

        # Write the data labels
        for label in self.data_block.columns:
            header.append("_%s" % label)

        # Make header string
        header = '\n'.join(header)

        # Print particle number info
        if verbose:
            print('Writing %d particles in %s' % (self.data_block.shape[0], out_fname))

        # Save file
        np.savetxt(out_fname, self.data_block.values, fmt=self.write_formatter, header=header, comments='')


class Cistem(Project):
    '''
    Cistem class
    '''
    def __init__(self, name='EMStar2Par'):
        super().__init__(name)

        self.par_file       = None
        self.par_data_block = None 

        self.original_star      = None
        self.original_star_file = None

        # Par to star parameters
        self.delclasses  = []
        self.selclasses  = []
        self.scorecutoff = None
        self.sigmacutoff = None
        self.mlcutoff    = None

        # Database variables
        self.db_conn        = None
        self.ref_packages   = []
        self.ref_nums       = []
        self.ref_num        = None
        self.orig_positions = None

    def set_par2star_params(self, delclasses=None, selclasses=None, scorecutoff=None, sigmacutoff=None, mlcutoff=None):
        '''
        Set par2star params
        '''
        self.delclasses  = delclasses
        self.selclasses  = selclasses
        self.scorecutoff = scorecutoff
        self.sigmacutoff = sigmacutoff
        self.mlcutoff    = mlcutoff

    def sort_images(self):
        '''
        Sort based on the image name
        '''
        self.particle_star.sort('rlnImageName')

    def read_par(self, fname):
        '''
        Read par file from Cistem
        '''
        if os.path.isfile(fname):
            self.par_data_block = np.genfromtxt(fname,
                                                skip_header=self.par_skip_rows,
                                                dtype=self.par_dtypes,
                                                comments='C')

    def read_db(self, db_file):
        '''
        Read Cistem database
        '''
        if db_file is not None and os.path.isfile(db_file):
            self.db_conn = sqlite3.connect(db_file)
        else:
            return

        # Get refinement package info
        self.ref_packages   = []
        self.ref_nums       = []
        self.orig_positions = []

        # Set a db cursor
        c = self.db_conn.cursor()

        # Iterate over the refinement packages
        for row in c.execute("select * from refinement_package_assets"):
            print('Refinement ID: %d - Refinement PackageName: %s' % (row[0], row[1]))
            self.ref_packages.append(row)
            self.ref_nums.append(int(row[0]))

        # Get reference info
        while self.ref_num not in self.ref_nums:
            self.ref_num = int(input('Enter the refinement package ID: '))

        # Get the original position ids
        for row in c.execute("select original_particle_position_asset_id from refinement_package_contained_particles_%d" % (self.ref_num)):
            self.orig_positions.append(int(row[0])-1)

    def select_particles(self):
        '''
        Select from particles based on the cistem database information
        '''
        if self.orig_positions is not None and len(self.orig_positions) > 0:
            self.particle_star.data_block = self.particle_star.data_block.iloc[self.orig_positions]

    def copy2star(self):
        '''
        Convert par data to star object
        '''
        if self.particle_star is not None and self.par_data_block is not None:
            if self.particle_star.data_block.shape[0] != self.par_data_block.shape[0]:
                sys.exit('Particle star and par file rows dont match. Exiting!')

            # Copy the data columns from par to star
            self.particle_star.data_block['rlnOriginX'] = -self.par_data_block['SHX']/self.particle_apix
            self.particle_star.data_block['rlnOriginY'] = -self.par_data_block['SHY']/self.particle_apix

            self.particle_star.data_block['rlnAnglePsi']  = self.par_data_block['PSI']
            self.particle_star.data_block['rlnAngleTilt'] = self.par_data_block['THETA']
            self.particle_star.data_block['rlnAngleRot']  = self.par_data_block['PHI']

    def create_write_formatter(self):
        '''
        Create write formatter
        '''
        self.par_skip_rows   = 1
        self.header_list     = ('C', 'PSI', 'THETA', 'PHI', 'SHX', 'SHY', 'MAG', 'INCLUDE', 'DF1', 'DF2', 'ANGAST', 'PSHIFT','OCC', 'LogP', 'SIGMA', 'SCORE', 'CHANGE')
        self.data_formats    = ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4', 'f4', 'f4', 'f4')
        self.write_header    = "%-7s%8s%8s%8s%10s%10s%8s%9s%6s%9s%8s%8s%8s%10s%11s%8s%8s" % self.header_list
        self.write_formatter = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%8.2f%10d%11.4f%8.2f%8.2f"

        # Par data types
        self.par_dtypes   = {'names': self.header_list,
                             'formats': self.data_formats}

    def write_output_file(self, verbose=True):
        # Save file
        # Print particle number info
        if verbose:
            print('Writing %d particles in %s' % (self.par_data.shape[0], self.particle_par_file))
        np.savetxt(self.particle_par_file, self.par_data.values, fmt=self.write_formatter, header=self.write_header, comments='')

    def write_star_file(self):
        if self.particle_star is not None:
            self.particle_star.write(self.particle_out_file)

    def prepare_project(self):
        '''
        Prepare project
        '''
        self.prepare_io_files_par()
        self.read_particle_apix()
        self.create_write_formatter()

    def prepare_io_files_par(self):
        # Copy input file to output directory
        if self.particle_star_file is not None:
            head, tail = os.path.split(self.particle_star_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.particle_star_file, self.output_directory+'/particle_input'+ext)
            self.particle_par_file = self.output_directory+'/particle_output.par'
            self.particle_out_file = self.output_directory+'/particle_output.star'

            # Make symlink
            self.make_symlink2parent(self.particle_star_file)

    def convert2par(self):
        '''
        Convert star dataframe to par data frame
        '''
        # 1. Initialize data frame with C column
        self.par_data      = pd.DataFrame(columns=self.header_list)
        self.par_data['C'] = np.arange(self.particle_star.data_block.shape[0], dtype=int) + 1 

        # 2. PSI angle
        if self.particle_star.has_label('rlnAnglePsi'):
            self.par_data['PSI'] = self.particle_star.data_block['rlnAnglePsi']
        else:
            self.par_data['PSI'] = 0

        # 3. THETA ANGLE
        if self.particle_star.has_label('rlnAngleTilt'):
            self.par_data['THETA'] = self.particle_star.data_block['rlnAngleTilt']
        else:
            self.par_data['THETA'] = 0

        # 4. PHI ANGLE
        if self.particle_star.has_label('rlnAngleRot'):
            self.par_data['PHI'] = self.particle_star.data_block['rlnAngleRot']
        else:
            self.par_data['PHI'] = 0

        # 5. SHX
        if self.particle_star.has_label('rlnOriginX'):
            self.par_data['SHX'] = -self.particle_star.data_block['rlnOriginX']*self.particle_apix
        else:
            self.par_data['SHX'] = 0

        # 6. SHY
        if self.particle_star.has_label('rlnOriginY'):
            self.par_data['SHY'] = -self.particle_star.data_block['rlnOriginY']*self.particle_apix
        else:
            self.par_data['SHY'] = 0

        # 7. MAG
        if self.particle_star.has_label('rlnMagnification'):
            self.par_data['MAG'] = self.particle_star.data_block['rlnMagnification']
        else:
            self.par_data['MAG'] = 0

        # 8. INCLUDE
        if self.particle_star.has_label('rlnGroupNumber'):
            self.par_data['INCLUDE'] = self.particle_star.data_block['rlnGroupNumber']
        else:
            self.par_data['INCLUDE'] = 0

        # 9. DF1
        if self.particle_star.has_label('rlnDefocusU'):
            self.par_data['DF1'] = self.particle_star.data_block['rlnDefocusU']
        else:
            self.par_data['DF1'] = 0

        # 10. DF2
        if self.particle_star.has_label('rlnDefocusV'):
            self.par_data['DF2'] = self.particle_star.data_block['rlnDefocusV']
        else:
            self.par_data['DF2'] = 0

        # 11. ANGAST
        if self.particle_star.has_label('rlnDefocusAngle'):
            self.par_data['ANGAST'] = self.particle_star.data_block['rlnDefocusAngle']
        else:
            self.par_data['ANGAST'] = 0

        # REST
        self.par_data['PSHIFT'] = 0.0
        self.par_data['OCC']    = 100
        self.par_data['LogP']   = 0.0
        self.par_data['SIGMA']  = 0.5
        self.par_data['SCORE']  = 0.0
        self.par_data['CHANGE'] = 0.0


class CryoSparc(EMfile):
    '''
        Cryosparc class
    '''
    def __init__(self):
        self.name = None
        self.data_block_dict        = {}
        self.data_block_blob        = None
        self.data_block_passthrough = None
        self.cs2star_blob           = {'blob/res_A': 'rlnEstimatedResolution',
                                       'ctf/amp_contrast': 'rlnAmplitudeContrast',
                                       'ctf/accel_kv': 'rlnVoltage',
                                       'ctf/cs_mm': 'rlnSphericalAberration',
                                       'ctf/df1_A': 'rlnDefocusU',
                                       'ctf/df2_A': 'rlnDefocusV',
                                       'ctf/ctf_fit_to_A': "rlnCtfMaxResolution",
                                       'alignments2D/class_posterior': "rlnMaxValueProbDistribution",
                                       'alignments2D/class': "rlnClassNumber"}

        self.cs2star_passthrough    = {'location/micrograph_path': 'rlnMicrographName',
                                       'ctf/amp_contrast': 'rlnAmplitudeContrast',
                                       'ctf/accel_kv': 'rlnVoltage',
                                       'ctf/cs_mm': 'rlnSphericalAberration',
                                       'ctf/df1_A': 'rlnDefocusU',
                                       'ctf/df2_A': 'rlnDefocusV',
                                       'ctf/ctf_fit_to_A': "rlnCtfMaxResolution",
                                       'alignments2D/class_posterior': "rlnMaxValueProbDistribution",
                                       'alignments2D/class': "rlnClassNumber"
                                       }

        self.star                   = None
        self.original_star          = None
        self.original_path          = ''

        self.project_path           = ''
        self.ref_mrc_file           = None
        self.ref_mrcs_file          = None

    def set_relion_image_handler_exe(self):
        '''
        Set relion image handler exe
        '''
        relion_process = subprocess.run(['which', 'relion_image_handler'], stdout=subprocess.PIPE, universal_newlines=True)
        self.relion_image_handler_exe = relion_process.stdout.strip()

    def remove_str_from_micrograph_names(self, del_str=''):
        '''
        Remove a string from micrograph names
        '''
        if self.star.has_label('rlnMicrographName'):
            self.star.data_block['rlnMicrographName'] = self.star.data_block.rlnMicrographName.replace({del_str: ""},regex=True)

    def adjust_class_number(self):
        '''
        Adjust class number so that it is in relion format starting from 1 to N
        '''
        if self.star.has_label('rlnClassNumber'):
            self.star.data_block['rlnClassNumber'] += 1

    def reset_class_numbers(self):
        '''
        Reset class numbers
        '''
        if self.star.has_label('rlnClassNumber'):
            # Get class number list
            prev_class_number_list = self.star.data_block['rlnClassNumber'].tolist()

            # Get previous numbers
            prev_class_numbers = np.sort(np.unique(prev_class_number_list))

            # New class numbers
            new_class_numbers = np.arange(1, len(prev_class_numbers)+1)

            # Blank class list
            blank_class_list  = np.zeros(len(prev_class_number_list))

            # Iterate over all the classes
            for i in range(len(prev_class_numbers)):
                prev_class_number = prev_class_numbers[i]

                # Get valid class numbers
                valid = prev_class_number_list == prev_class_number

                # Assign new class number
                blank_class_list[valid] = new_class_numbers[i]

            # Assign the list
            self.star.data_block.loc[:, 'rlnClassNumber'] = blank_class_list

    def delete_classes(self, del_classes=[]):
        '''
        Delete classes
        '''
        if self.star.has_label('rlnClassNumber'):
            keep_classes = ~self.star.data_block['rlnClassNumber'].isin(del_classes)
            self.star.data_block = self.star.data_block.loc[keep_classes, :]

    def get_ref_mrc_file(self):
        '''
        Get first ref mrc file
        '''
        if self.star is not None and self.star.has_label('rlnReferenceImage'):
            ref_index, self.ref_mrc_file = self.star.data_block['rlnReferenceImage'][0].split('@')

    def merge_with_original_star(self, restore_offsets=False):
        '''
        Merge with original star
        '''
        # Create shortImageName and delete rlnImageName for ptcl star
        self.star.create_shortImageName()

        # Create shortImage nam
        self.original_star.create_shortImageName()

        # Delete micrograph name in the new cs-star file
        self.star.delete_column('rlnMicrographName')

        # Comparison list
        cmp_list = ['shortImageName']

        # Delete the first set of candidate list in star file
        candidate_list = ['rlnCoordinateX','rlnCoordinateY', 'rlnImageName', 'rlnMicrographName']
        self.star.delete_columns(candidate_list)

        # Candidate column list
        candidate_list += ['shortImageName', 'rlnIsFlip', 'rlnParticleName']

        # Add the priors
        candidate_list +=  ['rlnOriginXPrior',
                            'rlnOriginYPrior',
                            'rlnAnglePsiPrior',
                            'rlnAngleRotPrior',
                            'rlnAngleTiltPrior']

        # If to restore the coordinate and angle offsets
        if restore_offsets:
            offset_params = ['rlnOriginX',
                             'rlnOriginY',
                             'rlnAnglePsi',
                             'rlnAngleRot',
                             'rlnAngleTilt']

            self.star.delete_columns(offset_params)

            # Add offset params to the candidate list
            candidate_list += offset_params

        # Final list
        final_list = [label for label in candidate_list if self.original_star.has_label(label)]

        # Get the original block
        original_block = self.original_star.data_block[final_list]

        # Get particle data block
        ptcl_block = self.star.get_data_block()

        # Merge star with original
        intersected_block = pd.merge(ptcl_block, original_block, on=cmp_list, how='inner')

        # Set data block for star
        self.star.set_data_block(intersected_block)

        # Finally remove shortImageName
        self.star.delete_column('shortImageName')

    def convert_ref_mrc_to_mrcs(self):
        '''
        Rename img file to mrcs
        '''
        # Get the path for relion image handler
        self.set_relion_image_handler_exe()

        if self.ref_mrc_file is not None:
            mrc_file, ext = os.path.splitext(self.ref_mrc_file)
            if os.path.isfile(self.ref_mrc_file) and ext == '.mrc':

                # Define the new files
                self.ref_mrcs_file        = mrc_file + '.mrcs'
                self.ref_mrcs_flipXY_file = mrc_file + '_flipXY.mrcs'

                # Make an mrcs copy of the mrc file
                shutil.copy(self.ref_mrc_file, self.ref_mrcs_file)

                # RUn relion image handler to flipXY ref mrc
                relion_args = [self.relion_image_handler_exe,
                               '--flipXY',
                               '--i', self.ref_mrcs_file,
                               '--o', self.ref_mrcs_flipXY_file]
                relion_subprocess = subprocess.run(relion_args,
                                                   stdout=subprocess.PIPE,
                                                   universal_newlines=True)

    def rename_ref_star_to_mrcs(self):
        '''
        Rename rlnImagename to *.mrcs
        '''
        if self.ref_mrc_file is not None:
            mrc_file, ext = os.path.splitext(self.ref_mrc_file)
            if os.path.isfile(self.ref_mrc_file) and ext == '.mrc':
                self.star.data_block.loc[:,'rlnReferenceImage'] = self.star.data_block.rlnReferenceImage.replace({".mrc": ".mrcs"},regex=True)

    def read_blob(self, file):
        '''
        Read cryosparc blob file
        '''
        data = np.load(file)
        self.data_block_blob = data

    def read_passthrough(self, file):
        '''
        Read passthrough file
        '''
        data = np.load(file)
        self.data_block_passthrough = data

    def has_label_blob(self, label):
        '''
        Check if label exists in blob data block
        '''
        if self.data_block_blob is not None and label in self.data_block_blob.dtype.names:
            return True
        else:
            return False

    def has_label_passthrough(self, label):
        '''
        Check if label exists in passthrough data block
        '''
        if self.data_block_passthrough is not None and label in self.data_block_passthrough.dtype.names:
            return True
        else:
            return False

    def read_original_star(self, fname):
        '''
        Read original star file
        '''
        # Read original path
        head, tail = os.path.split(fname)
        self.original_path = head+'/Micrographs/'

        self.original_star = Star(fname)

    def set_project_path(self, project_path=''):
        '''
        Set project path
        '''
        if len(project_path) > 0:
            self.project_path = project_path + '/'
        else:
            self.project_path = ''

    def convert2star(self):
        '''
        Convert to star format
        '''
        self.star            = Star()
        self.data_block_dict = {}

        # Get the parameters from blob file
        if self.data_block_blob is not None:

            # Do the  direct copies
            for cs_label, rln_label in self.cs2star_passthrough.items():
                if self.has_label_passthrough(cs_label):
                    self.data_block_dict[rln_label] = np.array(self.data_block_passthrough[cs_label],
                                                               dtype=self.star.PARAMETERS[rln_label]['nptype'])

            # Do the  direct copies
            for cs_label, rln_label in self.cs2star_blob.items():
                if self.has_label_blob(cs_label):
                    self.data_block_dict[rln_label] = np.array(self.data_block_blob[cs_label],
                                                               dtype=self.star.PARAMETERS[rln_label]['nptype'])

            # rlnImageName
            new_data_column = []
            if self.has_label_blob('blob/path') and self.has_label_blob('blob/idx'):
                for i in range(self.data_block_blob.shape[0]):

                    # Read the root and the file
                    image_name = "%010d@%s" % (self.data_block_blob['blob/idx'][i]+1,
                                               self.project_path+self.data_block_blob['blob/path'][i].decode("utf-8"))
                    new_data_column.append(image_name)

                self.data_block_dict['rlnImageName'] = np.array(new_data_column,
                                                                dtype=self.star.PARAMETERS['rlnImageName']['nptype'])

            # rlnClassNumber
            if self.has_label_blob('alignments2D/class'):
                self.data_block_dict['rlnClassNumber'] = np.array(self.data_block_blob['alignments2D/class'],
                                                                  dtype=self.star.PARAMETERS['rlnClassNumber']['nptype'])
            if self.has_label_passthrough('alignments2D/class'):
                self.data_block_dict['rlnClassNumber'] = np.array(self.data_block_passthrough['alignments2D/class'],
                                                                  dtype=self.star.PARAMETERS['rlnClassNumber']['nptype'])
            # rlnOriginX/Y
            if self.has_label_passthrough('alignments2D/shift'):
                self.data_block_dict['rlnOriginX'] = np.array(self.data_block_passthrough['alignments2D/shift'][:, 0],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])

                self.data_block_dict['rlnOriginY'] = np.array(self.data_block_passthrough['alignments2D/shift'][:, 1],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])
            if self.has_label_blob('alignments2D/shift'):
                self.data_block_dict['rlnOriginX'] = np.array(self.data_block_blob['alignments2D/shift'][:, 0],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])

                self.data_block_dict['rlnOriginY'] = np.array(self.data_block_blob['alignments2D/shift'][:, 1],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])
            if self.has_label_passthrough('alignments3D/shift'):
                self.data_block_dict['rlnOriginX'] = np.array(self.data_block_passthrough['alignments3D/shift'][:, 0],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])

                self.data_block_dict['rlnOriginY'] = np.array(self.data_block_passthrough['alignments3D/shift'][:, 1],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])
            if self.has_label_blob('alignments3D/shift'):
                self.data_block_dict['rlnOriginX'] = np.array(self.data_block_blob['alignments3D/shift'][:, 0],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])

                self.data_block_dict['rlnOriginY'] = np.array(self.data_block_blob['alignments3D/shift'][:, 1],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])

            # rlnAngleRot, rlnAngleTilt, rlnAnglePsi
            if self.has_label_passthrough('alignments2D/pose'):
                self.data_block_dict['rlnAnglePsi'] = np.array(util.rad2deg(self.data_block_passthrough['alignments2D/pose']),
                                                               dtype=self.star.PARAMETERS['rlnAnglePsi']['nptype'])
            if self.has_label_blob('alignments2D/pose'):
                self.data_block_dict['rlnAnglePsi'] = np.array(util.rad2deg(self.data_block_blob['alignments2D/pose']),
                                                               dtype=self.star.PARAMETERS['rlnAnglePsi']['nptype'])
            if self.has_label_passthrough('alignments3D/pose'):
                
                # Convert rotation matrix to eular angles
                radAngles   = np.array([util.rodriguez2euler(x) for x in self.data_block_passthrough['alignments3D/pose']])
                eulerAngles = util.rad2deg(radAngles) 
                
                self.data_block_dict['rlnAngleRot']  = np.array(eulerAngles[:, 0],
                                                               dtype=self.star.PARAMETERS['rlnAngleRot']['nptype'])
                self.data_block_dict['rlnAngleTilt'] = np.array(eulerAngles[:, 1],
                                                                dtype=self.star.PARAMETERS['rlnAngleTilt']['nptype'])
                self.data_block_dict['rlnAnglePsi']  = np.array(eulerAngles[:, 2],
                                                               dtype=self.star.PARAMETERS['rlnAnglePsi']['nptype'])
            if self.has_label_blob('alignments3D/pose'):
                
                # Convert rotation matrix to eular angles
                radAngles = np.array([util.rodriguez2euler(x) for x in self.data_block_blob['alignments3D/pose']])
                eulerAngles = util.rad2deg(radAngles)

                self.data_block_dict['rlnAngleRot']  = np.array(eulerAngles[:, 0],
                                                               dtype=self.star.PARAMETERS['rlnAngleRot']['nptype'])
                self.data_block_dict['rlnAngleTilt'] = np.array(eulerAngles[:, 1],
                                                                dtype=self.star.PARAMETERS['rlnAngleTilt']['nptype'])
                self.data_block_dict['rlnAnglePsi']  = np.array(eulerAngles[:, 2],
                                                               dtype=self.star.PARAMETERS['rlnAnglePsi']['nptype'])

            # rlnCoordianteX/Y
            if(self.has_label_passthrough('location/center_x_frac') and
               self.has_label_passthrough('location/center_y_frac')):
                coordinate_x = np.round(self.data_block_passthrough['location/micrograph_shape'][:, 1]*self.data_block_passthrough['location/center_x_frac'])
                self.data_block_dict['rlnCoordinateX'] = np.array(coordinate_x, dtype=self.star.PARAMETERS['rlnCoordinateX']['nptype'])

                coordinate_y = np.round(self.data_block_passthrough['location/micrograph_shape'][:, 0]*self.data_block_passthrough['location/center_y_frac'])
                self.data_block_dict['rlnCoordinateY'] = np.array(coordinate_y, dtype=self.star.PARAMETERS['rlnCoordinateY']['nptype'])

            # Defocus and phase-shift angles
            if self.has_label_passthrough('ctf/df_angle_rad'):
                self.data_block_dict['rlnDefocusAngle'] = np.array(util.rad2deg(self.data_block_passthrough['ctf/df_angle_rad']),
                                                                   dtype=self.star.PARAMETERS['rlnDefocusAngle']['nptype'])

            if self.has_label_passthrough('ctf/phase_shift_rad'):
                self.data_block_dict['rlnPhaseShift'] = np.array(util.rad2deg(self.data_block_passthrough['ctf/phase_shift_rad']),
                                                                 dtype=self.star.PARAMETERS['rlnPhaseShift']['nptype'])

            # Defocus and phase-shift angles
            if self.has_label_blob('ctf/df_angle_rad'):
                self.data_block_dict['rlnDefocusAngle'] = np.array(util.rad2deg(self.data_block_blob['ctf/df_angle_rad']),
                                                                   dtype=self.star.PARAMETERS['rlnDefocusAngle']['nptype'])

            if self.has_label_blob('ctf/phase_shift_rad'):
                self.data_block_dict['rlnPhaseShift'] = np.array(util.rad2deg(self.data_block_blob['ctf/phase_shift_rad']),
                                                                 dtype=self.star.PARAMETERS['rlnPhaseShift']['nptype'])

            # Create the data block for star
            self.star.data_block = pd.DataFrame.from_dict(self.data_block_dict)

    def copy_from_original(self, mic_path='Micrographs'):
        '''
        Copy from original star
        '''
        if self.original_star is not None:
            if self.original_star.has_label('rlnDetectorPixelSize'):
                self.star.data_block['rlnDetectorPixelSize'] = self.original_star.data_block['rlnDetectorPixelSize'][0]

            if self.original_star.has_label('rlnMagnification'):
                self.star.data_block['rlnMagnification'] = self.original_star.data_block['rlnMagnification'][0]

        if mic_path is not None:
            self.star.create_micname_from_imagename(mic_path)

    def rename_star_columns(self, columns={}):
        '''
        Replace column names
        '''
        self.star.rename_columns(columns)

    def convert_idx_to_classnumber(self):
        '''
        Convert idx to classnumber
        '''
        if self.star is not None and self.star.data_block is not None:
            self.star.data_block['rlnClassNumber'] = np.array(self.data_block_blob['blob/idx'],
                                                              dtype=self.star.PARAMETERS['rlnClassNumber']['nptype'])