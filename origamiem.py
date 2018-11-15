#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 11:27:03
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import re
import glob
import sys
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
        self.name           = name
        self.particle_star  = None
        self.ref_class_star = None
        self.class_ids      = None

        # Input files
        self.particle_star_file  = None
        self.ref_class_star_file = None

        # Output files
        self.particle_out_file  = None
        self.ref_class_out_file = None

    def read_class_refs(self, file):
        '''
        Read class refs
        '''
        self.ref_class_star_file = os.path.abspath(file)
        self.ref_class_star = Star(file)
        self.ref_class_star.num2className()

    def read_particles(self, file):
        '''
        Read particle star
        '''
        self.particle_star_file = os.path.abspath(file)
        self.particle_star = Star(file)

    def get_class_ids(self):
        '''
        Get class names
        '''
        self.class_ids =  self.particle_star.get_class_ids()

    def recenter_particles(self):
        '''
        Recenter particles
        '''
        self.particle_star.recenter2D()

    def add_columns(self, column_params=None):
        '''
        Add new columns
        '''
        if not column_params is None:
            for label, value in column_params.items():
                self.particle_star.set_column(label, value)

    def transform_particles(self):
        '''
        Transform particle star file based on the class star file
        '''
        if self.particle_star is None or self.ref_class_star is None:
            print('No transformation due to missing particle or reference data')
            return 0

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
                                        ptcls=class_ptcls)
        return 1

    def write_output_files(self):
        '''
        Write output files
        '''
        if not self.particle_star is None:
            self.particle_star.write(self.particle_out_file)

        if not self.ref_class_star is None:
            self.ref_class_star.write(self.ref_class_out_file)

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

    def prepare_io_files(self):
        # Copy input file to output directory
        if not self.particle_star_file is None:
            head, tail = os.path.split(self.particle_star_file)
            root, ext  = os.path.splitext(tail)
            copyfile(self.particle_star_file, self.output_directory+'/'+root+'_particle_input'+ext)
            self.particle_out_file = self.output_directory+'/'+root+'_particle_output'+ext

        if not self.ref_class_star_file is None:
            head, tail = os.path.split(self.ref_class_star_file)
            copyfile(self.ref_class_star_file, self.output_directory+'/'+root+'_class_input'+ext)
            self.ref_class_out_file = self.output_directory+'/'+root+'_class_output'+ext

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


class Mrc:
    '''
        MRC  class
    '''
    def __init__(self):
        self.name = None


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
        self._file_path     = None                         #Path to directory where this script stays 
        self.star_file      = None
        self.name           = None
        self.data_block     = None
        self.data_name      = None
        self.data_labels    = None
        self.data_formats   = None
        self.data_dtypes    = None
        self.data_skip_rows = None
        self.metadata_file  = 'relion_metadata_labels.dat' 
        self.PARAMETERS     = {}
        self.str2type       = {'double': 'f4',
                               'string':'U100',
                               'int':'i4',
                               'bool':'b'}

        self.type2format    = {'double': "%13.6f",
                               'string': "%s",
                               'int':    "%12d",
                               'bool':   "%2d"}


        # Read metadata labels
        self._read_metadata_labels()

        # Read file
        self.read(file)

    def _read_metadata_labels(self):
        '''
        Read the metadata labels
        '''
        self._file_path  = os.path.dirname(os.path.abspath(__file__))

        f = open(self._file_path + '/'+ self.metadata_file,'r')
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
                                               'type': self.str2type[type_name], 
                                               'description': description}

    def _read_header(self, file):
        '''
        Read star header file
        '''
        f = open(file,'r')

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
        self.data_formats  = [self.PARAMETERS[label]['type'] for label in self.data_labels]
        self.data_dtypes   = { 'names': self.data_labels,
                               'formats': self.data_formats }

    def read(self, file):
        '''
        Read Star file and create the data blocks
        '''
        self._read_header(file)
        self._prepare_data_types()
        self.data_block = np.loadtxt(file,
                                    skiprows = self.data_skip_rows,
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
            print('%s is not a valid Star label.'%(label))
            return None
        elif label in self.data_labels:
            print('%s exists. Not creating a new column.'%(label))
            return None

        # Create new column
        new_data_column = np.empty([self.num_data_points,1], dtype=self.PARAMETERS[label]['type'])

        # Append the data column
        self.data_block[label] = new_data_column
        self.data_labels.append(label)

        # Initialize the column
        if self.PARAMETERS[label]['type'][0] == 'U':
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
        if self.has_label(label) and not value is None:
            self.data_block.loc[:, label] = value
        else:
            success = self.add_column(label)
            if success and not value is None:
                self.data_block.loc[:, label] = value

    def copy(self, other=None):
        '''
        Deep copy other star to self
        '''
        self.data_block   = other.data_block.copy()
        self.data_labels  = other.data_labels.copy()
        self.data_formats = other.data_formats.copy()
        self.data_dtypes  = other.data_dtypes.copy()

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

    def recenter2D(self):
        '''
        Recenter particles
        '''
        if(self.has_label('rlnOriginX') and
           self.has_label('rlnOriginY')):

            # Center x-coordinate
            remx, intoffsetx   = np.modf(self.data_block['rlnOriginX'])
            self.data_block.loc[:, 'rlnOriginX']     = remx
            self.data_block.loc[:,'rlnCoordinateX'] -= intoffsetx

            # Center y-coordinate
            remy, intoffsety   = np.modf(self.data_block['rlnOriginY'])
            self.data_block.loc[:, 'rlnOriginY']      = remy
            self.data_block.loc[:, 'rlnCoordinateY'] -= intoffsety

    def rotate2D(self, rotangle=0, offset=[0, 0], ptcls=None):
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
        for ptcl in ptcls:
            oldangle = self.data_block['rlnAnglePsi'][ptcl]
            rotM = util.euler2rot2D(float(-oldangle))

            # Get the corrected offset
            corrected_offset = np.array(offset).dot(rotM)
            
            # Update the offsets
            self.data_block.loc[ptcl, 'rlnOriginX'] += corrected_offset[0]
            self.data_block.loc[ptcl, 'rlnOriginY'] += corrected_offset[1]


        # Update psi angles
        self.data_block.loc[ptcls, 'rlnAnglePsi'] += rotangle%360

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
        self.write_formatter = '  '.join(formatter) + '\n'

    def write(self, out_fname):
        '''
        Write star file
        '''

        # Create the formatter
        self.create_write_formatter()

        # Create header
        header = []

        # Write data block name
        header.append("data_%s"%self.data_name)

        #Write the data labels
        for label in self.data_block.columns:
            header.append("_%s"%label)

        # Make header string
        header='\n'.join(header)

        # Save file
        np.savetxt(out_fname, self.data_block.values,fmt=self.write_formatter,header=header,comments='')


class CryoSparc(EMfile):
    '''
        Cryosparc class
    '''
    def __init__(self):
        self.name = None
        self.blob_data        = None
        self.passthrough_data = None
        self.cs2starmap       = {}

    def read_blob(self, file):
        '''
        Read cryosparc blob file
        '''

    def read_passthrough(self, file):
        '''
        Read passthrough file
        '''

    def convert2star(self):
        '''
        Convert to star format
        '''
