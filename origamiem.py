#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 11:27:03
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import numpy as np
import scipy


class Relion:
    def __init__(self):
            self.name = None


class CS:
    def __init__(self):
        self.name = None


class Project:
    '''
        Project File
    '''
    def __init__(self):
        self.name = None


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
        self.cs2starmap = {}


class Star(EMfile):
    '''
        Star class
    '''
    def __init__(self):
        self._file_path    = None                         #Path to directory where this script stays 

        self.name          = None
        self.data_blocks   = []
        self.metadata_file = 'relion_metadata_labels.dat' 

    def _read_metadata_labels(self):
        '''
        Read the metadata labels
        '''
        self._file_path  = os.path.dirname(os.path.abspath(__file__))

        f = read(self._file_path + '/'+ self.metadata_file)
        

    def read(file):
        '''
        Read Star file and create the data blocks
        '''




class CryoSparc(EMfile):
    '''
        Cryosparc class
    '''
    def __init__(self):
        self.name = None


def main():
    print('Hello World')
