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

from shutil import copyfile


class Relion:
    def __init__(self):
            self.name = None


class Project:
    '''
        Project File
    '''
    def __init__(self, name='ProjectOrigami'):
        self.name            = name
        self.particle_star   = None
        self.ref_class_star  = None
        self.micrograph_star = None
        self.class_ids       = None

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

        # Micrograph files
        self.first_mic_file  = None
        self.first_mic_mrc   = None

        # Micrograph dimensions
        self.mic_NX = 0
        self.mic_NY = 0
        self.mic_NZ = 0

        # Cryosparc objects
        self.particle_cs    = None

        # Cryosparc files
        self.blob_cs_file        = None
        self.passthrough_cs_file = None
        self.original_star_file  = None

        # Additional data frames
        self.particle_data_props = pd.DataFrame(columns=['insideFrame'])

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
        self.mic_NX = self.first_mic_mrc.header['NX']
        self.mic_NY = self.first_mic_mrc.header['NY']
        self.mic_NZ = self.first_mic_mrc.header['NZ']

    def set_mic_apix(self, apix=1.82):
        '''
        Set micrograph apix
        '''
        self.mic_apix = apix

    def set_particle_apix(self, apix=1.82):
        '''
        Set particle apix
        '''
        self.particle_apix = apix

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

    def read_class_refs(self, file, new_classids=False):
        '''
        Read class refs
        '''
        self.ref_class_star_file = os.path.abspath(file)
        self.ref_class_star = Star(file)

        # Generate new classids from particle numbers in image names
        if new_classids:
            self.ref_class_star.num2className()

    def read_particles(self, file):
        '''
        Read particle star
        '''
        self.particle_star_file = os.path.abspath(file)
        self.particle_star = Star(file)

    def read_micrographs(self, file):
        '''
        Read micrograph star file
        '''
        self.micrograph_star_file = os.path.abspath(file)
        self.micrograph_star      = Star(file)

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

    def copy_from_ref(self, columns={}):
        '''
        Copy columns from reference star
        '''
        if self.particle_star.get_numRows() != self.ref_class_star.get_numRows():
            print('Number of rows in the particle and reference star files dont match')
            return 0

        # Iterate over all columns
        self.particle_star.copy_columns(self.ref_class_star, columns)

    def add_columns(self, column_params=None):
        '''
        Add new columns
        '''
        if column_params is not None:
            for label, value in column_params.items():
                self.particle_star.set_column(label, value)

    def transform_particles(self, final_offset=[0, 0]):
        '''
        Transform particle star file based on the class star file
        '''
        if self.particle_star is None:
            print('No transformation due to missing particle data')
            return 0

        if self.ref_class_star is not None:
            # Ref data block
            ref_data_block = self.ref_class_star.get_data_block()

            # Iterate through every class
            for i in range(ref_data_block.shape[0]):
                # Get class id
                class_id = ref_data_block['rlnClassNumber'][i]

                # Get rotangle
                rot_angle = ref_data_block['rlnAnglePsi'][i]

                # Get offsetx, offsety
                offset_x = ref_data_block['rlnOriginX'][i]
                offset_y = ref_data_block['rlnOriginY'][i]

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

        return 1

    def write_output_files(self, write_particle_star=True, write_ref_class_star=True, write_cs_star=True):
        '''
        Write output files
        '''
        if self.particle_star is not None and write_particle_star:
            self.particle_star.write(self.particle_out_file)

        if self.ref_class_star is not None and write_ref_class_star:
            self.ref_class_star.write(self.ref_class_out_file)

        print(self.particle_out_file)
        if self.particle_cs is not None and write_cs_star:
            self.particle_cs.star.write(self.particle_out_file)

    def set_output_directory(self, input_filename, output_directory=None):
        '''
        Set output directory
        '''

        # Parse input file
        head, tail       = os.path.split(os.path.abspath(input_filename))
        root, ext        = os.path.splitext(tail)

        if output_directory:
            self.output_directory = output_directory
        else:
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

    def prepare_io_files_star(self):
        # Copy input file to output directory
        if self.particle_star_file is not None:
            head, tail = os.path.split(self.particle_star_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.particle_star_file, self.output_directory+'/'+root+'_particle_input'+ext)
            self.particle_out_file = self.output_directory+'/'+root+'_particle_output'+ext

        if self.ref_class_star_file is not None:
            head, tail = os.path.split(self.ref_class_star_file)
            copyfile(self.ref_class_star_file, self.output_directory+'/'+root+'_class_input'+ext)
            self.ref_class_out_file = self.output_directory+'/'+root+'_class_output'+ext

    def prepare_io_files_cs(self):
        # Copy input files to output directory
        if self.blob_cs_file is not None:
            head, tail = os.path.split(self.blob_cs_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.blob_cs_file, self.output_directory+'/'+root+'_blob_input'+ext)
            self.particle_out_file = self.output_directory+'/'+root+'_particle_output.star'

        if self.passthrough_cs_file is not None:
            head, tail = os.path.split(self.passthrough_cs_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.passthrough_cs_file, self.output_directory+'/'+root+'_passthrough_input'+ext)

        if self.original_star_file is not None:
            head, tail = os.path.split(self.original_star_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.original_star_file, self.output_directory+'/'+root+'_original_particle'+ext)

    def set_cs_files(self, blob_cs_file=None, passthrough_cs_file=None, original_star_file=None):
        '''
        Set input cs files
        '''
        self.blob_cs_file        = blob_cs_file
        self.passthrough_cs_file = passthrough_cs_file
        self.original_star_file  = original_star_file

    def read_cs_files(self):
        '''
        Read cs file
        '''
        self.particle_cs = CryoSparc()

        if self.blob_cs_file is not None:
            self.particle_cs.read_blob(self.blob_cs_file)

        if self.passthrough_cs_file is not None:
            self.particle_cs.read_passthrough(self.passthrough_cs_file)

        if self.original_star_file is not None:
            self.particle_cs.read_original_star(self.original_star_file)

    def convert_cs2star(self):
        '''
        Convert to cs to star file
        '''
        self.particle_cs.convert2star()
        self.particle_cs.copy_from_original()


class Micrograph:
    '''
        Micrograph class
    '''
    def __init__(self):
        self.name = None


class Particle:
    '''
        Particle Class
    '''
    def __init__(self):
        self.name = None


class MRC:
    '''
        MRC  class
    '''
    def __init__(self, file):
        self.name = None
        self.header =  {'NX': 0,
                        'NY': 0,
                        'NZ': 0,
                        'CELLAX': 0,
                        'CELLAY': 0,
                        'CELLAZ': 0}

        if file is not None:
            self.read_header(file)

    def read(self):
        '''
        Read MRC file
        '''

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

        with open(file) as f:
            # Go to 0 pointer
            f.seek(0)

            # Read header
            header_int   = np.fromfile(f, dtype=np.int32, count=256)
            header_float = header_int.view(np.float32)

            # Read the micrograph and cell dimensions
            [self.header['NX'], self.header['NY'], self.header['NZ'], self.header['MODE']] = header_int[:4]
            [self.header['CELLAX'], self.header['CELLAY'], self.header['CELLAZ']] = header_float[10:13]
            if self.header['CELLAX'] == self.header['CELLAY'] == self.header['CELLAZ'] == 0:
                self.header['CELLAX'] = self.header['NX']
                self.header['CELLAY'] = self.header['NY']
                self.header['CELLAZ'] = self.header['NZ']

            # Close file
            f.close()
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
                               'bool':   bool}
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

        # Read file
        if file is not None:
            self.read(file)

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

    def add_column(self, label=None):
        '''
        Add new column to data block
        '''
        if label not in self.PARAMETERS:
            print('%s is not a valid Star label.' % (label))
            return None
        elif label in self.data_labels:
            print('%s exists. Not creating a new column.' % (label))
            return None

        # Create new column
        new_data_column = np.empty([self.num_data_points, 1], dtype=self.PARAMETERS[label]['nptype'])

        # Append the data column
        self.data_block[label] = new_data_column
        self.data_labels.append(label)

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
            self.data_block.loc[:, label] = value
        else:
            success = self.add_column(label)
            if success and value is not None:
                self.data_block.loc[:, label] = self.PARAMETERS[value]['type'](value)

    def copy(self, other=None):
        '''
        Deep copy other star to self
        '''
        self.data_block   = other.data_block.copy()
        self.data_labels  = other.data_labels.copy()
        self.data_formats = other.data_formats.copy()
        self.data_dtypes  = other.data_dtypes.copy()

    def copy_columns(self, other, columns={}):
        '''
        Copy columns from another star object
        '''
        # If the sizes don't match don't perform copy
        if other.data_block.shape[0] != self.data_block.shape[0]:
            return 0

        # Iterate over all columns
        for label, value in columns.items():
            if other.has_label(label):
                if self.has_label(label):
                    self.data_block.loc[:, label] = other.data_block[label]
                else:
                    self.add_column(label)
                    self.set_column(label, value)

    def has_label(self, label):
        '''
        Check if the label exists in data frame
        '''
        if label in self.data_block.columns:
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
        apix_scale = 1.0*self.star_apix/mic_apix
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
        if self.has_label(old_label) and new_label in self.PARAMETERS:
            self.data_block.rename(columns={old_label: new_label},
                                   inplace=True)

    def dublicate_column(self, label, new_label):
        '''
        Duplicate a column with a label
        '''
        if self.has_label(label) and new_label in self.PARAMETERS:
            self.data_block.loc[:, new_label] = self.data_block[label]

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
        self.data_block.loc[ptcls, 'rlnAnglePsi'] += rotangle % 360

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

    def write(self, out_fname):
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

        # Save file
        np.savetxt(out_fname, self.data_block.values, fmt=self.write_formatter, header=header, comments='')


class CryoSparc(EMfile):
    '''
        Cryosparc class
    '''
    def __init__(self):
        self.name = None
        self.data_block_dict        = {}
        self.data_block_blob        = None
        self.data_block_passthrough = None
        self.cs2star_blob           = {}
        self.cs2star_passthrough    = {'location/micrograph_path': 'rlnMicrographName',
                                       'ctf/amp_contrast': 'rlnAmplitudeContrast',
                                       'ctf/accel_kv': 'rlnVoltage',
                                       'ctf/cs_mm': 'rlnSphericalAberration',
                                       'ctf/df1_A': 'rlnDefocusU',
                                       'ctf/df2_A': 'rlnDefocusV',
                                       'ctf/ctf_fit_to_A': "rlnCtfMaxResolution"
                                       }

        self.star                   = None
        self.original_star          = None

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
        if label in self.data_block_blob.dtype.names:
            return True
        else:
            return False

    def has_label_passthrough(self, label):
        '''
        Check if label exists in passthrough data block
        '''
        if label in self.data_block_passthrough.dtype.names:
            return True
        else:
            return False

    def read_original_star(self, fname):
        '''
        Read original star file
        '''
        self.original_star = Star(fname)

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

            # rlnImageName
            new_data_column = []
            if self.has_label_blob('blob/path') and self.has_label_blob('blob/idx'):
                for i in range(self.data_block_blob.shape[0]):
                    image_name = "%010d@%s" % (self.data_block_blob['blob/idx'][i], self.data_block_blob['blob/path'][i])
                    new_data_column.append(image_name)

                self.data_block_dict['rlnImageName'] = np.array(new_data_column,
                                                                dtype=self.star.PARAMETERS['rlnImageName']['nptype'])

            # rlnOriginX/Y
            if self.has_label_passthrough('alignments2D/shift'):
                self.data_block_dict['rlnOriginX'] = np.array(self.data_block_passthrough['alignments2D/shift'][:, 0],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])

                self.data_block_dict['rlnOriginY'] = np.array(self.data_block_passthrough['alignments2D/shift'][:, 1],
                                                              dtype=self.star.PARAMETERS['rlnOriginX']['nptype'])

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

            # Create the data block for star
            self.star.data_block = pd.DataFrame.from_dict(self.data_block_dict)

    def copy_from_original(self):
        '''
        Copy from original star
        '''

        if self.original_star is not None:
            if self.original_star.has_label('rlnDetectorPixelSize'):
                self.star.data_block['rlnDetectorPixelSize'] = self.original_star.data_block['rlnDetectorPixelSize'][0]

            if self.original_star.has_label('rlnMagnification'):
                self.star.data_block['rlnMagnification'] = self.original_star.data_block['rlnMagnification'][0]

            # Get micrographs path
            if self.original_star.has_label('rlnMicrographName'):
                original_mic_head, original_mic_tail = os.path.split(self.original_star.data_block['rlnMicrographName'][0])

                new_mic_names = []
                for i in range(self.star.data_block.shape[0]):
                    mic_head, mic_tail = os.path.split(self.star.data_block['rlnMicrographName'][i])

                    # Create the micrograph name
                    new_mic_name = original_mic_head+'/'+mic_tail
                    new_mic_names.append(new_mic_name)
                self.star.data_block['rlnMicrographName'] = new_mic_names
