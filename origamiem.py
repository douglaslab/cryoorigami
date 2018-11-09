#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-11-09 11:27:03
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import numpy as np
import scipy


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
        self.name = None


class CryoSparc(EMfile):
    '''
        Cryosparc class
    '''
    def __init__(self):
        self.name = None


def main():
    print('Hello World')
