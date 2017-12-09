#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2017-10-18 10:40:08
# @Author  : Tural Aksel (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

from EMAN2 import *

import EMtools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

import scipy.ndimage
import scipy.signal
import os
import sys
import copy

#GLOBAL variables
line_edge_effect = [pe.Stroke(linewidth=3, foreground='k'), pe.Normal()]

class EMProject:
    def __init__(self):
        self.micrograph_set     = None
        self.initial_3D_model   = None 
        self.refined_3D_model   = None
        self.class_averages_2D  = None  
        self.class_averages_3D  = None
        self.class_references   = None

        self.LINE_REFS          = None
        self.SMALL_ARM_REFS     = None
        self.LONG_ARM_REFS      = None
        self.BARCODE_REFS       = None
        self.HINGE_REFS         = None

    def read_line_refs(self):
        '''
        Read line references 
        '''

    def read_small_arm_refs(self):
        '''
        Read small arm refs
        '''

    def read_long_arm_refs(self):
        '''
        Read long arm references
        '''
    def read_hinge_refs(self):
        '''
        Read hinge references
        '''

    def read_hinge_barcode_refs(self):
        '''
        Read hinge barcode refs
        '''

class EM3DModel:
    def __init__(self):
        self.model_map         = None

class EMClassAverage:
    def __init__(self):
        self.class_average     = None

class EMMicrographSet:
    def __init__(self):
        self.particles           = None
        self.micrographs         = None
        self.em_models           = None
        self.project_name        = None
        self.folder_name         = None

        self.voltage             = None
        self.amplitude_constrast = None
        self.pixel_size          = None
        self.cs                  = None

        self.particle_radius     = None
        self.arm_radius          = None
        self.hinge_radius        = None 

class EMMicrograph:
    def __init__(self, particle_output_dir = "particles"):
        self.micrograph_file        = None
        self.micrograph_eman        = None
        self.micrograph_img         = None
        self.micrograph_nx          = None
        self.micrograph_ny          = None
        self.micrograph_mean        = None
        self.micrograph_flat        = None
        
        self.e2simmx_file           = None
        self.e2simmx_img            = None
        self.e2simmx_scores         = None
        self.e2simmx_dx             = None
        self.e2simmx_dy             = None
        self.e2simmx_drot           = None
        self.e2simmx_class_ids      = None
        self.e2simmx_class_dict     = None
        self.e2simmx_particles      = None
        self.e2simmx_particle_file  = None
        self.e2simmx_particle_coords= None
        self.e2simmx_particle_rots  = None

        self.arm_particles          = None
        self.arm_coords             = None
        self.arm_rots               = None
        self.arm_available          = None

        self.hinge_particles_1D     = None
        self.hinge_particles_2D     = None
        self.hinge_quality_scores   = None
        self.hinge_cmp_scores       = None
        
        self.LINE_REFS              = None
        self.SMALL_ARM_REFS         = None
        self.LONG_ARM_REFS          = None
        
        self.SHORT_BARCODE_REFS     = None
        self.LONG_BARCODE_REFS      = None

        self.ARM_WIDTH_SPAN         = 120
        self.ARM_LENGTH_SPAN        = 186
        self.HALFARM_WIDTH          = 30
        self.HALFARM_LENGTH         = 93

        self.PLOT_TEXT              = True
        self.PLOT_BARCODE           = True

        self.HINGE_BOXSIZE          = 360

        self.DISTANCE_CUTOFF        = {"low":1*self.HALFARM_WIDTH, "high":2.2*self.HALFARM_LENGTH}

        self.ctf_parameters         = None

        self.particle_output_dir    = particle_output_dir

    def _read_refs(self,ref_file, dtype='numpy'):
        '''
        Read references
        '''
        ref_imgs = []
        num_refs = EMUtil.get_image_count(ref_file)

        for i in range(num_refs):
            line_ref_data = EMData()
            line_ref_data.read_image(ref_file,i)
            
            #If the desired data type is numpy, transform data to numpy array
            if dtype == 'numpy': line_ref_data = np.copy(EMNumPy.em2numpy(line_ref_data))
            
            ref_imgs.append(line_ref_data)

        return ref_imgs

    def normalize_micrograph(self):
        '''
        Normalize micrograph
        '''
        self.micrograph_max = np.max(self.micrograph_img)
        self.micrograph_min = np.min(self.micrograph_img)

        self.micrograph_img = (self.micrograph_img-self.micrograph_min)/(self.micrograph_max-self.micrograph_min)

    def make_flat(self):
        '''
        Make a flat copy of the micrograph
        '''
        self.micrograph_flat = self.micrograph_img.flatten()


    def calculate_mean(self):
        '''
        Determine mean intensity
        '''
        self.micrograph_mean = np.mean(self.micrograph_img)

    def read_line_refs(self,line_ref_file):
        '''
        Read line refs
        '''
        self.LINE_REFS = self._read_refs(line_ref_file)

    def read_small_arm_refs(self,small_arm_ref_file):
        '''
        Read small arm refs
        '''
        self.SMALL_ARM_REFS = self._read_refs(small_arm_ref_file)

    def read_long_arm_refs(self,long_arm_ref_file):
        '''
        Read long arm refs
        '''
        self.LONG_ARM_REFS = self._read_refs(long_arm_ref_file)

    def read_short_barcode_refs(self, short_barcode_ref_file):
        '''
        Read short barcode refs
        '''
        self.SHORT_BARCODE_REFS = self._read_refs(short_barcode_ref_file)


    def read_long_barcode_refs(self, long_barcode_ref_file):
        '''
        Read short barcode refs
        '''
        self.LONG_BARCODE_REFS  = self._read_refs(long_barcode_ref_file)

    def read_hinge_refs(self,hinge_refs_file):
        '''
        Read no bit hinge ref files
        '''
        self.HINGE_REFS     = self._read_refs(hinge_refs_file, dtype='eman')

    def read_hingebits_refs(self,hingebits_refs_file):
        '''
        Read no bit hinge ref files
        '''
        self.HINGEBITS_REFS = self._read_refs(hingebits_refs_file ,dtype='eman')

    def read_micrograph(self,micrograph_file):
        '''
        Read micrograph file
        '''
        self.micrograph_file = micrograph_file
        self.micrograph_eman = EMData(self.micrograph_file)
        self.micrograph_img  = np.copy(EMNumPy.em2numpy(self.micrograph_eman))
        self.micrograph_nx,self.micrograph_ny = self.micrograph_img.shape

    def read_e2simmx_output(self,e2simmx_file):
        '''
        Read e2simmx output
        '''
        self.e2simmx_file = e2simmx_file
        
        #1. Get the scores
        simmx_scores_data = EMData()
        simmx_scores_data.read_image(self.e2simmx_file,0)
        self.e2simmx_scores = np.copy(EMNumPy.em2numpy(simmx_scores_data))

        #2. Get dx
        simmx_dx_data = EMData()
        simmx_dx_data.read_image(self.e2simmx_file,1)
        self.e2simmx_dx = np.copy(EMNumPy.em2numpy(simmx_dx_data))

        #3. Get dy
        simmx_dy_data = EMData()
        simmx_dy_data.read_image(self.e2simmx_file,2)
        self.e2simmx_dy = np.copy(EMNumPy.em2numpy(simmx_dy_data))

        #3. Get rot
        simmx_rot_data = EMData()
        simmx_rot_data.read_image(self.e2simmx_file,3)
        self.e2simmx_drot = np.copy(EMNumPy.em2numpy(simmx_rot_data))

        #Determine the number of e2simmx classses
        self.e2simmx_num_classes = self.e2simmx_scores.shape[1]

    def classify_e2simmx(self, skip_zero=True, score_cutoff={'low':-0.1,'high':-1.0}):
        '''
        Determine classes based on e2simmx output scores
        If skip zero, consider the rest of the scores
        '''

        #e2simmx class dictionary
        self.e2simmx_class_dict = {}

        if skip_zero:
            self.e2simmx_class_ids = np.argmin(self.e2simmx_scores[:,1:], axis=1) + 1
        else:
            self.e2simmx_class_ids = np.argmin(self.e2simmx_scores, axis=1)

        #Determine class members
        for class_id in range(self.e2simmx_num_classes):
            self.e2simmx_class_dict[class_id] = np.nonzero((self.e2simmx_class_ids == class_id)*
                                                           (self.e2simmx_scores[:,class_id]< score_cutoff['low'])*
                                                           (self.e2simmx_scores[:,class_id]>score_cutoff['high']))[0]

    def plot_micrograph(self):
        '''
        Plot micrograph in gray scale
        '''
        #Create a micrograph plot
        plt.figure(0,figsize=(self.micrograph_nx/100.0,self.micrograph_ny/100.0),dpi=100)
        ax = plt.gca()
        plt.imshow(self.micrograph_img,cmap='gray')
        plt.plot(0,0,'.')

    def save_plot(self,fname='micrograph.png',dpi=100):
        '''
        Save the micrograph plot
        '''
        #Get the current axis
        ax = plt.gca()

        #Make axis invisible
        ax.axis('off')
        
        #Save figure
        plt.savefig(fname,dpi=dpi,transparent=True)

    def show_plot(self):
        '''
        Show plot
        '''
        plt.xlim((0,self.micrograph_nx))
        plt.ylim((0,self.micrograph_ny))
        plt.show()

    def set_reference_angles(self,refangle_dict):
        '''
        Assign reference angles
        '''
        self.e2simmx_refangles = refangle_dict.copy()

    def read_e2simmx_particles(self,e2simmx_particle_file):
        '''
        Read e2simmx particles
        '''
        self.e2simmx_particle_file   = e2simmx_particle_file
        self.e2simmx_particles       = {}
        self.e2simmx_particle_coords = {}
        self.e2simmx_particle_rots   = {}

        #Read particles
        for class_id in range(self.e2simmx_num_classes):
            self.e2simmx_particles[class_id]       = []
            self.e2simmx_particle_coords[class_id] = []
            self.e2simmx_particle_rots[class_id]   = []

            for particle_id in self.e2simmx_class_dict[class_id]:
                #Read particle file
                new_particle_data          = EMData(self.e2simmx_particle_file,particle_id)
                
                #Create new particle
                new_particle                = EMParticle()
                new_particle.particle_file  = self.e2simmx_particle_file
                new_particle.particle_id    = particle_id
                new_particle.particle_eman  = new_particle_data
                new_particle.particle_coord = np.array(new_particle_data['ptcl_source_coord'])
                new_particle.particle_x     = new_particle.particle_coord[0]
                new_particle.particle_y     = new_particle.particle_coord[1]
                new_particle.particle_rot   = self.e2simmx_drot[particle_id,class_id] + self.e2simmx_refangles[class_id]

                #Add particle to list
                self.e2simmx_particles[class_id].append(new_particle)

                #Add coordinates and rotations to particle lists
                self.e2simmx_particle_coords[class_id].append(new_particle.particle_coord)
                self.e2simmx_particle_rots[class_id].append(new_particle.particle_rot)

            #Make the lists numpy array
            self.e2simmx_particle_coords[class_id] = np.array(self.e2simmx_particle_coords[class_id])
            self.e2simmx_particle_rots[class_id]   = np.array(self.e2simmx_particle_rots[class_id])


    def make_arm_particles(self,class_id = 1):
        '''
        Make arm particles from simmx particle list
        '''
        #Initialize the container
        self.arm_particles = []
        self.arm_coords    = []
        self.arm_rots      = []

        for particle in self.e2simmx_particles[class_id]:

            new_origami_arm = EMOrigamiArm()
            
            #Copy particle
            new_origami_arm.copy_particle(particle)

            #Make the arm available
            new_origami_arm.arm_available = True

            #Make the arm valid
            new_origami_arm.arm_valid     = True 

            #Assign the coordinate
            new_origami_arm.arm_center    = particle.particle_coord

            #Assign the angle
            new_origami_arm.arm_angle     = particle.particle_rot

            #Add arm particle to list
            self.arm_particles.append(new_origami_arm)
            self.arm_coords.append(particle.particle_coord)
            self.arm_rots.append(particle.particle_rot)

        #Make the coords and rots array to numpy array
        self.arm_coords = np.array(self.arm_coords)
        self.arm_rots   = np.array(self.arm_rots)

        #Initialize arm available
        self.arm_available = np.ones(len(self.arm_particles)) 

    def detect_hinges(self):
        '''
        Detect origami hinges from the arm particles
        '''

        #Initialize the 2D hinge container
        self.hinge_particles_1D   = []
        self.hinge_particles_2D   = [[None]*len(self.arm_particles) for i in range(len(self.arm_particles))]
        self.hinge_quality_scores = [[None]*len(self.arm_particles) for i in range(len(self.arm_particles))]
        self.hinge_cmp_scores     = [[None]*len(self.arm_particles) for i in range(len(self.arm_particles))]

        #Number of valid particles
        for i in range(len(self.arm_particles)-1):
            
            #Get the current particle
            current_arm_particle = self.arm_particles[i]
            
            #If current particle is not available, move with the next particle
            if current_arm_particle.arm_available == False:
                continue

            #Arm distances
            neighbor_distances = EMtools.eular_distance(current_arm_particle.arm_center, self.arm_coords[i+1:,:])
            neighbor_valid     = np.nonzero((neighbor_distances > self.DISTANCE_CUTOFF["low"])*
                                            (neighbor_distances < self.DISTANCE_CUTOFF["high"])*
                                            (self.arm_available[i+1:] == 1))[0]

            #If the size of valid array is 0, continue with next 
            if len(neighbor_valid) == 0:
                #Make the particle unavailable for the next iterations
                current_arm_particle.arm_valid     = False
                current_arm_particle.arm_available = False
                continue
            
            #Adjust the numbering in neighbor valid
            neighbor_valid +=  (i+1)
            
            #Go through each neighbor, keep the measures in an array
            local_results = []
            for neighbor in neighbor_valid:
                
                #Pick up the neighbour arm particle
                neighbor_arm_particle = self.arm_particles[neighbor]

                #Pick up the coordinates
                current_coord  = current_arm_particle.arm_center
                neighbor_coord = self.arm_coords[neighbor,:]

                #Pick up the angles
                current_angle  = current_arm_particle.particle_rot
                neighbor_angle = self.arm_rots[neighbor]

                #If the angles are equal skip it
                if current_angle == neighbor_angle:
                    continue

                #Find the intersection point and triangle parameters for each particle
                hinge_tip_point = EMtools.intersection_point(current_coord,current_angle,neighbor_coord,neighbor_angle)
                
                #If there is a problem with finding the hinge tip point, skip it
                if len(hinge_tip_point) == 1:
                    continue

                #Determine the which arm is left or right
                vector_1 = current_coord  - hinge_tip_point
                vector_2 = neighbor_coord - hinge_tip_point

                #Angle between vectors
                cross_product  = np.cross(vector_1,vector_2)
                 
                #Make a new hinge
                new_origami_hinge = EMOrigamiHinge()
                
                #Assign eman data
                new_origami_hinge.micrograph_eman = self.micrograph_eman

                #Assign micrograph image
                new_origami_hinge.micrograph_img  = self.micrograph_img

                #Assign micrograph mean intensity
                new_origami_hinge.micrograph_mean = self.micrograph_mean

                #Assign micrograph file name
                new_origami_hinge.micrograph_file = self.micrograph_file

                #Assign micrograph flat file
                new_origami_hinge.micrograph_flat = self.micrograph_flat

                #Assign the references
                new_origami_hinge.LINE_REFS           = self.LINE_REFS
                new_origami_hinge.SMALL_ARM_REFS      = self.SMALL_ARM_REFS
                new_origami_hinge.LONG_ARM_REFS       = self.LONG_ARM_REFS 
                
                new_origami_hinge.SHORT_BARCODE_REFS  = self.SHORT_BARCODE_REFS
                new_origami_hinge.LONG_BARCODE_REFS   = self.LONG_BARCODE_REFS

                new_origami_hinge.HINGE_REFS          = self.HINGE_REFS
                new_origami_hinge.HINGEBITS_REFS      = self.HINGEBITS_REFS

                if cross_product > 0:
                    new_origami_hinge.left_arm = EMOrigamiArm()
                    new_origami_hinge.left_arm.copy_arm(current_arm_particle)

                    new_origami_hinge.right_arm = EMOrigamiArm()
                    new_origami_hinge.right_arm.copy_arm(neighbor_arm_particle)

                    new_origami_hinge.left_arm_center = current_coord
                    new_origami_hinge.left_arm_angle  = current_angle

                    new_origami_hinge.right_arm_center = neighbor_coord
                    new_origami_hinge.right_arm_angle  = neighbor_angle
                else:
                    new_origami_hinge.left_arm = EMOrigamiArm()
                    new_origami_hinge.left_arm.copy_arm(neighbor_arm_particle)

                    new_origami_hinge.right_arm = EMOrigamiArm()
                    new_origami_hinge.right_arm.copy_arm(current_arm_particle)

                    new_origami_hinge.left_arm_center = neighbor_coord
                    new_origami_hinge.left_arm_angle  = neighbor_angle

                    new_origami_hinge.right_arm_center = current_coord
                    new_origami_hinge.right_arm_angle  = current_angle   

                #Perform the full optimization
                optimization_result = new_origami_hinge.optimization_protocol()

                #Check the result
                if optimization_result == False: new_origami_hinge.valid = False

                #Keep hinge if only it is valid
                if new_origami_hinge.valid:
                    
                    #Add hinge to 2D hinge container
                    self.hinge_particles_2D[i][neighbor]   = new_origami_hinge
                    self.hinge_particles_2D[neighbor][i]   = new_origami_hinge
                    
                    #Add the scores
                    self.hinge_quality_scores[i][neighbor] = new_origami_hinge.hinge_quality_score
                    self.hinge_cmp_scores[i][neighbor]     = new_origami_hinge.hinge_cmp_score
                    self.hinge_quality_scores[neighbor][i] = new_origami_hinge.hinge_quality_score
                    self.hinge_cmp_scores[neighbor][i]     = new_origami_hinge.hinge_cmp_score
                    
                    #Add hinge to hinge list
                    self.hinge_particles_1D.append(new_origami_hinge)
        
        #Convert arrays to numpy array
        self.hinge_quality_scores = np.array(self.hinge_quality_scores)
        self.hinge_cmp_scores     = np.array(self.hinge_cmp_scores)

    def remove_hinge_conflicts(self):
        '''
        Remove hinge conflicts. Some arms may be part of a several hinges.
        Make sure that each arm is part of only single hinge.
        '''
        
        for i in range(len(self.hinge_quality_scores)):
            
            #Determine the valid index
            valid_index = np.nonzero(self.hinge_quality_scores[i,:] != None)[0]

            #If there is no element continue
            if len(valid_index) == 0:
                continue
            
            #Determine the min quality score value
            best_index  = np.argmin(self.hinge_quality_scores[i,valid_index])
            best_index  = valid_index[best_index]

            #Make invalid all hinges except the best one 
            for index in list(set(valid_index)-set([best_index])):
                self.hinge_particles_2D[i][index].valid = False

    def count_hinges(self):
        '''
        Count the valid hinges
        '''
        self.num_valid_hinges = 0
        for hinge_particle in self.hinge_particles_1D:
            if hinge_particle.valid:
                self.num_valid_hinges += 1

        print('Number of valid hinges:%d'%(self.num_valid_hinges))
    
    def make_hinge_masks(self):
        '''
        For the valid hinges make the triangular and random masks
        '''
        for hinge_particle in self.hinge_particles_1D:
            if hinge_particle.valid:
                
                #Make triangular and random masks
                hinge_particle.make_masks()

    def extract_hinges(self):
        '''
        Extract hinge particles
        '''
        
        #Get the micrograph basename 
        self.micrograph_base = os.path.splitext(os.path.basename(self.micrograph_file))[0]

        #Hinge particle filenames for the micrograph
        out_hinge_fname            = self.particle_output_dir+'/'+self.micrograph_base+'-hinges.hdf'
        out_hingemodel_fname       = self.particle_output_dir+'/'+self.micrograph_base+'-hinges-model.hdf'

        out_maskedhinge_fname      = self.particle_output_dir+'/'+self.micrograph_base+'-masked-hinges.hdf'
        out_maskedhingemodel_fname = self.particle_output_dir+'/'+self.micrograph_base+'-masked-hinges-model.hdf'
        
        for hinge_particle in self.hinge_particles_1D:
            if hinge_particle.valid:

                #1. Make hinge eman 
                hinge_particle.make_hinge_eman()

                #2. Align to hinge model
                hinge_particle.align_to_hinge_angle_model()

                #3. Compare to barcoded models
                hinge_particle.compare_to_hinge_barcode_models()

                #4. Flip hinge based on the long barcode side
                if hinge_particle.hinge_barcode_flip: 
                    EMtools.eman_flip_y(hinge_particle.hinge_eman)
                    EMtools.eman_flip_y(hinge_particle.masked_hinge_eman)

                #5. Flip hinge based on the long barcode side
                if hinge_particle.hinge_model_flip: 
                    EMtools.eman_flip_y(hinge_particle.hinge_eman_model)
                    EMtools.eman_flip_y(hinge_particle.masked_hinge_eman_model)

                #6. Clip masked hinges
                hinge_particle.clip_masked_hinge()

                #7. Write hinges in a micrograph file
                hinge_particle.hinge_eman.append_image(out_hinge_fname)
                hinge_particle.hinge_eman_model.append_image(out_hingemodel_fname)

                #8. Write masked hinges in a micrograph file
                hinge_particle.masked_hinge_eman.append_image(out_maskedhinge_fname)
                hinge_particle.masked_hinge_eman_model.append_image(out_maskedhingemodel_fname)

                #9. Write hinges in different angle classes
                out_hinge_angle_fname            = self.particle_output_dir+'/hinges-%2d.hdf'%(int(hinge_particle.hinge_angle))
                out_hingebits_angle_fname        = self.particle_output_dir+'/hinges-%2d-%d.hdf'%(int(hinge_particle.hinge_angle),hinge_particle.hinge_barcode_class)

                out_hingemodel_angle_fname       = self.particle_output_dir+'/hinges-model-%2d.hdf'%(int(hinge_particle.hinge_angle))
                out_hingebitsmodel_angle_fname   = self.particle_output_dir+'/hinges-model-%2d-%d.hdf'%(int(hinge_particle.hinge_angle),hinge_particle.hinge_model_class)

                #Hinge based on barcode detection
                hinge_particle.hinge_eman.append_image(out_hinge_angle_fname)
                hinge_particle.hinge_eman.append_image(out_hingebits_angle_fname)

                #Hinge based on model comparison
                hinge_particle.hinge_eman_model.append_image(out_hingemodel_angle_fname)
                hinge_particle.hinge_eman_model.append_image(out_hingebitsmodel_angle_fname)

                #10. Write masked hinges in different angle classes
                out_maskedhinge_angle_fname            = self.particle_output_dir+'/masked-hinges-%2d.hdf'%(int(hinge_particle.hinge_angle))
                out_maskedhingebits_angle_fname        = self.particle_output_dir+'/masked-hinges-%2d-%d.hdf'%(int(hinge_particle.hinge_angle),hinge_particle.hinge_barcode_class)

                out_maskedhingemodel_angle_fname       = self.particle_output_dir+'/masked-hinges-model-%2d.hdf'%(int(hinge_particle.hinge_angle))
                out_maskedhingebitsmodel_angle_fname   = self.particle_output_dir+'/masked-hinges-model-%2d-%d.hdf'%(int(hinge_particle.hinge_angle),hinge_particle.hinge_model_class)

                #Hinge based on barcode detection
                hinge_particle.masked_hinge_eman.append_image(out_maskedhinge_angle_fname)
                hinge_particle.masked_hinge_eman.append_image(out_maskedhingebits_angle_fname)

                #Hinge based on model comparison
                hinge_particle.masked_hinge_eman_model.append_image(out_maskedhingemodel_angle_fname)
                hinge_particle.masked_hinge_eman_model.append_image(out_maskedhingebitsmodel_angle_fname)


    def plot_hinge_particles(self):
        '''
        Plot hinge particles
        '''
        LINEWIDTH=3
        for hinge_particle in self.hinge_particles_1D:
            color_spec = '#00FF00'
            if not hinge_particle.valid:
                color_spec = '#FF0000'
            
            #Make the tether for ideal coordinates
            plt.plot([hinge_particle.ideal_left_arm_inner_side[0],hinge_particle.ideal_right_arm_inner_side[0]],
                     [hinge_particle.ideal_left_arm_inner_side[1],hinge_particle.ideal_right_arm_inner_side[1]], color=color_spec, linestyle='-',linewidth=5)

            #Make the tether for original coordinates
            plt.plot([hinge_particle.left_arm_inner_side[0],hinge_particle.right_arm_inner_side[0]],
                     [hinge_particle.left_arm_inner_side[1],hinge_particle.right_arm_inner_side[1]],color=color_spec,linestyle='--',linewidth=5)

            #Make the left side
            plt.plot([hinge_particle.ideal_left_arm_inner_side[0],hinge_particle.hinge_tip[0]],
                     [hinge_particle.ideal_left_arm_inner_side[1],hinge_particle.hinge_tip[1]]     ,color='yellow'   ,linestyle='-',linewidth=LINEWIDTH)

            #Make the right side
            plt.plot([hinge_particle.ideal_right_arm_inner_side[0],hinge_particle.hinge_tip[0]],
                     [hinge_particle.ideal_right_arm_inner_side[1],hinge_particle.hinge_tip[1]]   ,color='#00FFFF'   ,linestyle='-',linewidth=LINEWIDTH)

            #Plot outer edge centers
            plt.plot([hinge_particle.ideal_left_arm_outer_edge[0]] ,[hinge_particle.ideal_left_arm_outer_edge[1]] ,color='yellow' ,marker='*',markersize=7)
            plt.plot([hinge_particle.ideal_right_arm_outer_edge[0]],[hinge_particle.ideal_right_arm_outer_edge[1]],color='#00FFFF',marker='*',markersize=7)

            #Make edge lines
            left_outer_edge_line  = EMtools.make_line(hinge_particle.ideal_left_arm_outer_edge ,hinge_particle.left_arm.four_line_angle ,2*self.HALFARM_LENGTH)
            right_outer_edge_line = EMtools.make_line(hinge_particle.ideal_right_arm_outer_edge,hinge_particle.right_arm.four_line_angle,2*self.HALFARM_LENGTH)

            plt.plot(left_outer_edge_line[:,0] ,left_outer_edge_line[:,1] ,color='yellow' ,linestyle='-',linewidth=LINEWIDTH)
            plt.plot(right_outer_edge_line[:,0],right_outer_edge_line[:,1],color='#00FFFF',linestyle='-',linewidth=LINEWIDTH)

            #Plot the hinge specs
            if self.PLOT_TEXT:
                plt.text(hinge_particle.hinge_tip[0] ,hinge_particle.hinge_tip[1] ,"%d-%.2f-%.2f"%(hinge_particle.hinge_angle,hinge_particle.hinge_cmp_score,hinge_particle.hinge_quality_score),fontsize=10)

            #Plot vertical vectors
            left_vertical_point  = hinge_particle.left_arm.ideal_arm_inner_edge  + self.HALFARM_LENGTH*hinge_particle.left_arm.vertical_vec
            right_vertical_point = hinge_particle.right_arm.ideal_arm_inner_edge + self.HALFARM_LENGTH*hinge_particle.right_arm.vertical_vec

            plt.plot([hinge_particle.left_arm.ideal_arm_inner_edge[0],left_vertical_point[0]],[hinge_particle.left_arm.ideal_arm_inner_edge[1],left_vertical_point[1]]    ,color='yellow' ,linestyle='-',linewidth=LINEWIDTH)
            plt.plot([hinge_particle.right_arm.ideal_arm_inner_edge[0],right_vertical_point[0]],[hinge_particle.right_arm.ideal_arm_inner_edge[1],right_vertical_point[1]],color='#00FFFF',linestyle='-',linewidth=LINEWIDTH)

            #Plot horizontal vectors
            left_horizontal_point  = hinge_particle.left_arm.ideal_arm_inner_edge  + self.HALFARM_WIDTH*hinge_particle.left_arm.horizontal_vec
            right_horizontal_point = hinge_particle.right_arm.ideal_arm_inner_edge + self.HALFARM_WIDTH*hinge_particle.right_arm.horizontal_vec

            plt.plot([hinge_particle.left_arm.ideal_arm_inner_edge[0],left_horizontal_point[0]],[hinge_particle.left_arm.ideal_arm_inner_edge[1],left_horizontal_point[1]]    ,color='yellow',linestyle='-' ,linewidth=LINEWIDTH)
            plt.plot([hinge_particle.right_arm.ideal_arm_inner_edge[0],right_horizontal_point[0]],[hinge_particle.right_arm.ideal_arm_inner_edge[1],right_horizontal_point[1]],color='#00FFFF',linestyle='-',linewidth=LINEWIDTH)

            #Plot barcode locations
            if self.PLOT_BARCODE: hinge_particle.plot_barcode_locations(self.PLOT_TEXT)

    def plot_hinge_masks(self):
        '''
        Plot hinge masks
        '''
        LINEWIDTH  = 2
        color_spec = 'black'
        for hinge_particle in self.hinge_particles_1D:

            if hinge_particle.valid:

                #Correct the mask points based on Hinge origin 
                micrograph_trapezoid_points       = hinge_particle.trapezoid_points       + hinge_particle.HINGE_BOXORIGIN
                micrograph_higher_triangle_points = hinge_particle.higher_triangle_points + hinge_particle.HINGE_BOXORIGIN
                micrograph_protein_points         = hinge_particle.protein_points         + hinge_particle.HINGE_BOXORIGIN

                #Find the center of mass for the masks
                com_trapezoid_points       = np.mean(micrograph_trapezoid_points,axis=0)
                com_higher_triangle_points = np.mean(micrograph_higher_triangle_points,axis=0)
                com_protein_points         = np.mean(micrograph_protein_points,axis=0)

                #Point lists
                micrograph_trapezoid_point_list_x       = micrograph_trapezoid_points[:,0].tolist() + [micrograph_trapezoid_points[0,0]]
                micrograph_trapezoid_point_list_y       = micrograph_trapezoid_points[:,1].tolist() + [micrograph_trapezoid_points[0,1]]

                micrograph_higher_triangle_point_list_x = micrograph_higher_triangle_points[:,0].tolist() + [micrograph_higher_triangle_points[0,0]]
                micrograph_higher_triangle_point_list_y = micrograph_higher_triangle_points[:,1].tolist() + [micrograph_higher_triangle_points[0,1]]

                #Plot the mask edges
                plt.plot(micrograph_trapezoid_point_list_x,       micrograph_trapezoid_point_list_y,       color=color_spec ,linestyle='-',linewidth=LINEWIDTH)
                plt.plot(micrograph_higher_triangle_point_list_x, micrograph_higher_triangle_point_list_y, color=color_spec ,linestyle='-',linewidth=LINEWIDTH)

                #Write mean intensities in each mask area
                if self.PLOT_TEXT:
                    plt.text(com_trapezoid_points[0],      com_trapezoid_points[1]      ,'%.5f'%(hinge_particle.hinge_lower_random_mean_intensity), fontsize=10 ,color=color_spec)
                    plt.text(com_higher_triangle_points[0],com_higher_triangle_points[1],'%.5f'%(hinge_particle.hinge_higher_random_mean_intensity),fontsize=10,color=color_spec)
                    plt.text(com_protein_points[0],        com_protein_points[1]        ,'%.5f'%(hinge_particle.hinge_protein_mean_intensity),fontsize=10,color=color_spec)


    def plot_optimization_paths(self, point_spacing = 10):
        '''
        Plot the optimization paths taken by the hinges
        '''
        for hinge_particle in self.hinge_particles_1D:
            
            #Pick the left and right arms
            left_arm_particle  = hinge_particle.left_arm
            right_arm_particle = hinge_particle.right_arm

            #1. Plot the first optimization path
            left_cmp_scores  = left_arm_particle.arm_optimization[2]['cmp_scores']
            right_cmp_scores = right_arm_particle.arm_optimization[2]['cmp_scores']

            for cmp_score in left_cmp_scores:
                plt.plot([cmp_score[0]],[cmp_score[1]],'yo',markersize=5)

            for cmp_score in right_cmp_scores:
                plt.plot([cmp_score[0]],[cmp_score[1]],'bo',markersize=5)

            #2. Plot the second optimization path
            left_cmp_scores  = left_arm_particle.arm_optimization[3]['cmp_scores']
            right_cmp_scores = right_arm_particle.arm_optimization[3]['cmp_scores']  

            for cmp_score in left_cmp_scores[::point_spacing]:
                plt.plot([cmp_score[0]],[cmp_score[1]],'yv',markersize=5)

            for cmp_score in right_cmp_scores[::point_spacing]:
                plt.plot([cmp_score[0]],[cmp_score[1]],'bv',markersize=5)
            
            #3. Plot the third optimization path
            left_cmp_scores  = left_arm_particle.arm_optimization[4]['cmp_scores']
            right_cmp_scores = right_arm_particle.arm_optimization[4]['cmp_scores']

            for cmp_score in left_cmp_scores:
                plt.plot([cmp_score[0]],[cmp_score[1]],'ys',markersize=5)

            for cmp_score in right_cmp_scores:
                plt.plot([cmp_score[0]],[cmp_score[1]],'bs',markersize=5)

            #4. Plot the fifth optimization path
            left_cmp_scores  = left_arm_particle.arm_optimization[5]['cmp_scores']
            right_cmp_scores = right_arm_particle.arm_optimization[5]['cmp_scores']

            for cmp_score in left_cmp_scores[::point_spacing]:
                plt.plot([cmp_score[0]],[cmp_score[1]],'yp',markersize=5)

            for cmp_score in right_cmp_scores[::point_spacing]:
                plt.plot([cmp_score[0]],[cmp_score[1]],'bp',markersize=5)

            #5. Plot the sixth optimization path
            left_cmp_scores  = left_arm_particle.arm_optimization[6]['cmp_scores']
            right_cmp_scores = right_arm_particle.arm_optimization[6]['cmp_scores']

            for cmp_score in left_cmp_scores[::point_spacing]:
                plt.plot([cmp_score[0]],[cmp_score[1]],'y*',markersize=5)

            for cmp_score in right_cmp_scores[::point_spacing]:
                plt.plot([cmp_score[0]],[cmp_score[1]],'b*',markersize=5)

    def plot_e2simmx_particles(self):
        '''
        Plot e2simmx particles
        '''
        #Tag the bad particles in red - class 0 is bad

        self.bad_class_id = 0 
        
        if self.PLOT_TEXT:
            for particle in self.e2simmx_particles[self.bad_class_id]:
                plt.text(particle.particle_x,particle.particle_y,str(particle.particle_id),color='red') 

        #Tag the good particles in blue
        for class_id in range(1,self.e2simmx_num_classes):
            for particle in self.e2simmx_particles[class_id]:
                if self.PLOT_TEXT:
                    plt.text(particle.particle_x,particle.particle_y,str(particle.particle_id),color='blue')
               
                #Make the lines
                line_1 = EMtools.make_line(particle.particle_coord,   particle.particle_rot,2*self.HALFARM_LENGTH)
                line_2 = EMtools.make_line(particle.particle_coord,90+particle.particle_rot,2*self.HALFARM_WIDTH)

                #Plot the lines centered around particle
                plt.plot([line_1[0,0],line_1[1,0]],[line_1[0,1],line_1[1,1]],color='#FF00FF',linestyle='-',linewidth=3)
                plt.plot([line_2[0,0],line_2[1,0]],[line_2[0,1],line_2[1,1]],color='#FF00FF',linestyle='-',linewidth=3)

class EMParticle:
    def __init__(self):
        self.particle_file         = None
        self.particle_id           = None
        self.particle_eman         = None
        self.particle_source_coord = None
        self.particle_x            = None
        self.particle_y            = None
        self.particle_rot          = None

        self.micrograph_eman       = None
        self.micrograph_img        = None

        self.PARTICLE_RADIUS       = None
        self.PARTICLE_BOXSIZE      = None


    def copy_particle(self,other_particle):
        '''
        Elementwise copy of the particle object
        '''
        self.particle_file         = other_particle.particle_file
        self.particle_id           = other_particle.particle_id
        self.particle_eman         = other_particle.particle_eman
        self.particle_source_coord = other_particle.particle_source_coord
        self.particle_x            = other_particle.particle_x
        self.particle_y            = other_particle.particle_y
        self.particle_rot          = other_particle.particle_rot
        
        self.PARTICLE_RADIUS       = other_particle.PARTICLE_RADIUS
        self.PARTICLE_BOXSIZE      = other_particle.PARTICLE_BOXSIZE

    def extract_particle(self,boxsize):
        '''
        Function: Extracts particle image
        
        Parameters:

            boxsize: Particle boxsize 

        Returns: 2D numpy array
        ''' 


class EMParticleSet:
    def __init__(self):
        self.particles      = None

class EMOrigami(EMParticle):
    def __init__(self):
        self.origami_center = None

class EMOrigamiArm(EMOrigami):
    def __init__(self):
        self.arm_id           = None
        self.arm_center       = None
        self.arm_angle        = None
        self.arm_inner_edge   = None
        self.arm_outer_edge   = None
        self.arm_cmp_score    = None
        self.arm_available    = None
        self.arm_valid        = None
        self.arm_optimization = {} 

    def copy_arm(self,other_arm):
        #1. Copy the particle attributes
        self.copy_particle(other_arm)

        #2. Copy arm attributes
        self.arm_id           = other_arm.arm_id
        self.arm_available    = other_arm.arm_available
        self.arm_valid        = other_arm.arm_valid
        self.arm_angle        = other_arm.arm_angle

        self.arm_center       = copy.deepcopy(other_arm.arm_center)
        self.arm_inner_edge   = copy.deepcopy(other_arm.arm_inner_edge)
        self.arm_outer_edge   = copy.deepcopy(other_arm.arm_outer_edge)
        self.arm_optimization = other_arm.arm_optimization.copy()

class EMOrigamiHinge(EMOrigami):
    def __init__(self):
        self.left_arm                   = None
        self.left_arm_angle             = None
        self.left_arm_center            = None
        self.left_arm_inner_edge        = None
        self.left_arm_inner_side        = None
        self.left_arm_outer_edge        = None
        self.left_arm_outer_side        = None
        self.left_arm_mask_side         = None

        self.right_arm                  = None
        self.right_arm_angle            = None
        self.right_arm_center           = None
        self.right_arm_inner_edge       = None
        self.right_arm_inner_side       = None
        self.right_arm_outer_edge       = None
        self.right_arm_outer_side       = None
        self.right_arm_mask_side        = None

        self.ideal_left_arm_center      = None
        self.ideal_left_arm_inner_edge  = None
        self.ideal_left_arm_outer_edge  = None

        self.ideal_right_arm_center     = None 
        self.ideal_right_arm_inner_edge = None 
        self.ideal_right_arm_outer_edge = None

        self.ideal_center               = None
        self.ideal_cmp_score            = None

        self.hinge_angle                = None
        self.hinge_tip                  = None
        self.hinge_tip_angle            = None
        self.hinge_center               = None
        self.hinge_cmp_score            = None
        self.hinge_quality_score        = None
        
        self.hinge_barcode_class        = 0
        self.hinge_barcode_flip         = 0
        self.hinge_model_class          = 0
        self.hinge_model_flip           = 0
        self.hinge_model_cmp_score      = 0
        self.hingebits_model_cmp_score  = 0
        
        self.hinge_triangle_mask        = None
        self.hinge_random_mask          = None
        
        self.hinge_optimization         = []

        self.model_transform            = None
        self.horizontal_rotate          = None

        self.micrograph_mean            = None
        self.micrograph_eman            = None
        self.micrograph_img             = None
        self.micrograph_file            = None

        self.valid                      = True
        self.out_of_bounds              = False
        self.bad_cmp_score              = False
        self.bad_quality_score          = False
        self.bad_angle                  = False
        self.bad_displacement           = False

        self.HALFARM_WIDTH              = 30
        self.HALFARM_LENGTH             = 93
        self.ARM_WIDTH                  = 60
        self.ARM_LENGTH                 = 186

        self.BARCODE_FULL_LENGTH        = 186
        self.BARCODE_LONG_LENGTH        = 102
        self.BARCODE_SHORT_LENGTH       =  66
        self.BARCODE_OFFSET             = -10
        self.CENTER_OFFSET              = -5
        self.PROTEIN_OFFSET             = 25

        self.ARM_WIDTH_SPAN             = 120
        self.ARM_LENGTH_SPAN            = 180
        self.ARM_BOXSIZE                = 140
        self.ARM_HALF_BOXSIZE           = 70
        self.ARM_ANGLE_CUTOFF           = 10

        self.DIST_TO_EDGE               = 5
        self.DIST_TO_MASK               = 9
        self.LINE_TO_LINE_DIST          = 17
        self.MAXIMA_WINDOW              = 7
        
        self.HINGE_ANGLE_CUTOFF         = {'low':20,'high':60}
        self.HINGE_QUAL_CUTOFF          = 0.5   #Old values: 0.3 (2017-11-29)
        self.HINGE_CMP_CUTOFF           = 0.9   #Old values:1.01 (2017-11-29)
        self.HINGE_BOXSIZE              = 360
        self.HINGE_BOXORIGIN            = np.array([0,0])
        self.HINGE_HALFBOXSIZE          = 180
        self.HINGE_NUM_CLASSES          = 2*5

        self.PARTICLE_BOXSIZE           = 100

        self.LINE_REFS                  = None
        self.SMALL_ARM_REFS             = None
        self.LONG_ARM_REFS              = None

        self.SHORT_BARCODE_REFS         = None
        self.LONG_BARCODE_REFS          = None
        self.FULL_BARCODE_REFS          = None

        self.HINGE_REFS                 = None
        self.HINGEBITS_REFS             = None

        self.BARCODE_BOXSIZE            = 140
        self.FULLBARCODE_BOXSIZE        = 280

        self.PLOT_TEXT                  = True
        self.PLOT_BARCODE               = True

        self.long_barcode_index         = None
        self.short_barcode_index        = None

    def clip_masked_hinge(self):
        '''
        Clip masked eman images
        '''
        self.masked_hinge_eman       = self.masked_hinge_eman.get_clip(Region(self.HINGE_HALFBOXSIZE-0.5*self.PARTICLE_BOXSIZE,self.HINGE_HALFBOXSIZE-0.5*self.PARTICLE_BOXSIZE,self.PARTICLE_BOXSIZE,self.PARTICLE_BOXSIZE))
        self.masked_hinge_eman_model = self.masked_hinge_eman_model.get_clip(Region(self.HINGE_HALFBOXSIZE-0.5*self.PARTICLE_BOXSIZE,self.HINGE_HALFBOXSIZE-0.5*self.PARTICLE_BOXSIZE,self.PARTICLE_BOXSIZE,self.PARTICLE_BOXSIZE))

    def make_triangle_mask(self):
        '''
        Make triangular mask for hinge inside
        '''

        #Determine inner side base coordinates
        self.ideal_left_arm_mask_base  = self.ideal_left_arm_mask_side  + self.left_arm.vertical_vec*self.HALFARM_LENGTH
        self.ideal_right_arm_mask_base = self.ideal_right_arm_mask_side + self.right_arm.vertical_vec*self.HALFARM_LENGTH

        #Determine the box origin
        self.HINGE_BOXORIGIN = self.ideal_center - 0.5*self.HINGE_BOXSIZE

        #Triangle points
        self.triangle_points = np.vstack((self.ideal_mask_tip,self.ideal_left_arm_mask_base,self.ideal_right_arm_mask_base))

        #make the triangular mask
        self.hinge_triangle_mask = EMtools.create_polygon_mask(self.triangle_points,self.HINGE_BOXORIGIN,self.HINGE_BOXSIZE)

    def make_protein_mask(self):
        '''
        Make triangular mask for hinge inside
        '''

        #In order to determine mask side locations, we need to move along DNA tether
        #Horizontal vector needs to be corrected by the 1.0/cos(0.5*hinge_angle)

        #Determine lower base coordinates
        self.ideal_left_arm_protein_higher_base  = self.ideal_left_arm_mask_side   - self.left_arm.vertical_vec*self.PROTEIN_OFFSET/self.cos_half_hinge_angle
        self.ideal_right_arm_protein_higher_base = self.ideal_right_arm_mask_side  - self.right_arm.vertical_vec*self.PROTEIN_OFFSET/self.cos_half_hinge_angle

        #Determine lower base coordinates
        self.ideal_left_arm_protein_lower_base  = self.ideal_left_arm_mask_side  + self.left_arm.vertical_vec*self.PROTEIN_OFFSET/self.cos_half_hinge_angle
        self.ideal_right_arm_protein_lower_base = self.ideal_right_arm_mask_side + self.right_arm.vertical_vec*self.PROTEIN_OFFSET/self.cos_half_hinge_angle

        #Triangle points
        self.protein_points      = np.vstack((self.ideal_left_arm_protein_lower_base, self.ideal_right_arm_protein_lower_base, self.ideal_right_arm_protein_higher_base, self.ideal_left_arm_protein_higher_base))

        #make the triangular mask
        self.hinge_protein_mask  = EMtools.create_polygon_mask(self.protein_points,self.HINGE_BOXORIGIN,self.HINGE_BOXSIZE)

        #Determine mean random intensity
        self.hinge_protein_mean_intensity = np.sum(self.hinge_protein_mask*self.hinge_img)/np.sum(self.hinge_protein_mask)

        #Make small triangle mask
        self.higher_triangle_points = np.vstack((self.ideal_mask_tip, self.ideal_left_arm_protein_higher_base,self.ideal_right_arm_protein_higher_base))

        #Make higher triangle mask
        self.hinge_higher_triangle_mask = EMtools.create_polygon_mask(self.higher_triangle_points,self.HINGE_BOXORIGIN,self.HINGE_BOXSIZE)

        #Determine mean random intensity
        self.hinge_higher_random_mean_intensity = np.sum(self.hinge_higher_triangle_mask*self.hinge_img)/np.sum(self.hinge_higher_triangle_mask)

    def make_random_mask(self):
        '''
        Make mask with random intensities from micrograph
        '''

        #Trapezoid points
        self.trapezoid_points     = np.vstack((self.ideal_left_arm_protein_lower_base,self.ideal_left_arm_mask_base,self.ideal_right_arm_mask_base,self.ideal_right_arm_protein_lower_base))

        #Make a trapezoid mask
        self.hinge_trapezoid_mask = EMtools.create_polygon_mask(self.trapezoid_points, self.HINGE_BOXORIGIN, self.HINGE_BOXSIZE)
        
        #Retrieve the random points
        self.hinge_random         = self.hinge_trapezoid_mask*self.hinge_img
        self.hinge_random_flat    = self.hinge_random[np.nonzero(self.hinge_random)]

        #Determine mean random intensity
        self.hinge_lower_random_mean_intensity = np.sum(self.hinge_random)/np.sum(self.hinge_trapezoid_mask)

        #make the random mask
        self.hinge_random_mask   = EMtools.random_mask(self.hinge_random_flat,self.HINGE_BOXSIZE)


    def make_masks(self):
        '''
        Make triangle and the random mask
        '''
        #1. Make triangle mask
        self.make_triangle_mask()

        #2. Make protein mask
        self.make_protein_mask()

        #3. Make random mask
        self.make_random_mask()

    def make_hinge_eman(self):
        '''
        Extract hinge particle in the eman format
        '''

        #Clip for hinge_eman
        self.hinge_eman = self.micrograph_eman.get_clip(Region(self.ideal_center[0]-self.HINGE_HALFBOXSIZE,self.ideal_center[1]-self.HINGE_HALFBOXSIZE, self.HINGE_BOXSIZE, self.HINGE_BOXSIZE))

        #Save a numpy copy
        self.hinge_img  = EMtools.get_region(self.micrograph_img, self.ideal_center[0]-self.HINGE_HALFBOXSIZE, self.ideal_center[1]-self.HINGE_HALFBOXSIZE, self.HINGE_BOXSIZE, self.HINGE_BOXSIZE)

        #Make the hinge masks
        self.make_masks()

        #Store the hinge parameters

        #Center coordinate
        self.hinge_eman["ptcl_source_coord"]        = (self.ideal_center[0],self.ideal_center[1])
        self.hinge_eman["ptcl_source_image"]        = self.micrograph_file
        
        self.hinge_eman["hinge_angle"]              = self.hinge_angle
        self.hinge_eman["hinge_tip_angle"]          = self.hinge_tip_angle
        self.hinge_eman["hinge_tip"]                = (self.hinge_tip[0],self.hinge_tip[1])
        self.hinge_eman["hinge_valid"]              = int(self.valid)
        self.hinge_eman["hinge_cmp_score"]          = self.hinge_cmp_score
        self.hinge_eman["hinge_quality_score"]      = self.hinge_quality_score
        self.hinge_eman["hinge_barcode_class"]      = self.hinge_barcode_class
        self.hinge_eman["hinge_barcode_flip"]       = self.hinge_barcode_flip
        self.hinge_eman["hinge_model_class"]        = self.hinge_model_class
        self.hinge_eman["hinge_model_flip"]         = self.hinge_model_flip

        self.hinge_eman["left_arm_four_line_angle"] =  self.left_arm.four_line_angle
        self.hinge_eman["left_arm_vertical_vec"]    = (self.left_arm.vertical_vec[0]    ,self.left_arm.vertical_vec[1])
        self.hinge_eman["left_arm_horizontal_vec"]  = (self.left_arm.horizontal_vec[0]  ,self.left_arm.horizontal_vec[1])
        self.hinge_eman["left_arm_inner_side"]      = (self.ideal_left_arm_inner_side[0],self.ideal_left_arm_inner_side[1])

        self.hinge_eman["right_arm_four_line_angle"]=  self.right_arm.four_line_angle
        self.hinge_eman["right_arm_vertical_vec"]   = (self.right_arm.vertical_vec[0]    ,self.right_arm.vertical_vec[1])
        self.hinge_eman["right_arm_horizontal_vec"] = (self.right_arm.horizontal_vec[0]  ,self.right_arm.horizontal_vec[1])
        self.hinge_eman["right_arm_inner_side"]     = (self.ideal_right_arm_inner_side[0],self.ideal_right_arm_inner_side[1])

        self.hinge_eman["mask_tip"]                 = (self.ideal_mask_tip[0],self.ideal_mask_tip[1])
        self.hinge_eman["left_arm_mask_side"]       = (self.ideal_left_arm_mask_side[0] ,self.ideal_left_arm_mask_side[1])
        self.hinge_eman["right_arm_mask_side"]      = (self.ideal_right_arm_mask_side[0],self.ideal_right_arm_mask_side[1])
        self.hinge_eman["triangle_points"]          = (self.triangle_points[0,0],  self.triangle_points[0,1],  self.triangle_points[1,0],  self.triangle_points[1,1],  self.triangle_points[2,0],  self.triangle_points[2,1])
        self.hinge_eman["trapezoid_points"]         = (self.trapezoid_points[0,0], self.trapezoid_points[0,1], self.trapezoid_points[1,0] ,self.trapezoid_points[1,1], self.trapezoid_points[2,0], self.trapezoid_points[2,1], self.trapezoid_points[3,0], self.trapezoid_points[3,1])
        self.hinge_eman["protein_points"]           = (self.protein_points[0,0],   self.protein_points[0,1],   self.protein_points[1,0],   self.protein_points[1,1],   self.protein_points[2,0],   self.protein_points[2,1],   self.protein_points[3,0],   self.protein_points[3,1])

        self.hinge_eman["hinge_higher_random_mean_intensity"]  = self.hinge_higher_random_mean_intensity
        self.hinge_eman["hinge_lower_random_mean_intensity"]   = self.hinge_lower_random_mean_intensity
        self.hinge_eman["hinge_protein_mean_intensity"]        = self.hinge_protein_mean_intensity

        #Apply the masks onto masked hinge eman image
        self.apply_masks()

        #Rotate hinge so that tip points to right
        self.horizontal_rotate = Transform()
        self.horizontal_rotate.set_params({"type":"eman","phi":self.hinge_tip_angle})

        self.hinge_eman.transform(self.horizontal_rotate)
        self.masked_hinge_eman.transform(self.horizontal_rotate)

        #Store transformation information
        self.hinge_eman['xform.hinge']        = self.horizontal_rotate
        self.masked_hinge_eman['xform.hinge'] = self.horizontal_rotate

    def align_to_hinge_angle_model(self):
        '''
        Align the hinge eman to the hinge angle model
        '''
        hinge_angle_index     = int(self.hinge_angle) - int(self.HINGE_ANGLE_CUTOFF["low"])
        hinge_angle_model     = self.HINGE_REFS[hinge_angle_index] 
        self.hinge_eman_model = self.hinge_eman.align('translational',hinge_angle_model,{'maxshift':50},'dot',{})
        self.model_transform  = self.hinge_eman_model['xform.align2d']

        #Make a deep copy for the aligned hinge eman
        self.masked_hinge_eman_model = EMData(self.masked_hinge_eman)

        #Transform the image with model alignment parameters
        self.masked_hinge_eman_model.transform(self.model_transform)

    def apply_masks(self):
        '''
        Apply the masks onto masked hinge eman images
        '''

        #Sum the background and the masked image
        self.masked_hinge_img   = (1-self.hinge_triangle_mask)*self.hinge_random_mask + self.hinge_triangle_mask*self.hinge_img
        self.masked_hinge_eman  = EMData(EMNumPy.numpy2em(self.masked_hinge_img))

        #Make the masked protein image
        self.masked_protein_img = self.hinge_protein_mask*self.hinge_img
        self.masked_protein_eman= EMData(EMNumPy.numpy2em(self.masked_protein_img))

    def compare_to_hinge_barcode_models(self):
        '''
        Compare hinge to barcode models of a particular angle
        '''
        
        hinge_angle_index   = int(self.hinge_angle) - int(self.HINGE_ANGLE_CUTOFF["low"])
        hinge_barcode_index = hinge_angle_index*self.HINGE_NUM_CLASSES

        cmp_scores = []
        for i in range(self.HINGE_NUM_CLASSES):
            barcode_model = self.HINGEBITS_REFS[i+hinge_barcode_index]
            qual_match    = self.hinge_eman_model.cmp('dot',barcode_model,{})

            cmp_scores.append([i,qual_match])

        #Sort scores and assign the model state to hinge
        cmp_scores     = np.array(cmp_scores)
        sorted_results = cmp_scores[cmp_scores[:,1].argsort()]
        best_result    = sorted_results[0,:]

        #Assign the scores and model flip parameter
        self.hinge_model_flip          = best_result[0]%2
        self.hinge_model_class         = int(best_result[0]/2)
        self.hingebits_model_cmp_score = best_result[1] 

        #Assign the new values to hinge eman model
        self.hinge_eman_model["hinge_model_flip"]                 = self.hinge_model_flip
        self.hinge_eman_model["hinge_model_class"]                = self.hinge_model_class
        self.hinge_eman_model["hingebits_model_cmp_score"]        = self.hingebits_model_cmp_score

        #Assign the new values to masked hinge eman model
        self.masked_hinge_eman_model["hinge_model_flip"]          = self.hinge_model_flip
        self.masked_hinge_eman_model["hinge_model_class"]         = self.hinge_model_class
        self.masked_hinge_eman_model["hingebits_model_cmp_score"] = self.hingebits_model_cmp_score

    def determine_ideal_hinge_parameters(self, dist_to_edge = 5):
        '''
        Determine and assign ideal hinge parameters
        '''
        if not self.valid:
            return False

        (hinge_ideal_left_arm_edge, 
         hinge_ideal_right_arm_edge, 
         hinge_ideal_origin, 
         hinge_angle, 
         hinge_inter_angle, 
         hinge_inter_point, 
         hinge_ideal_inter_dist) = EMtools.determine_hinge_parameters(self.left_arm_inner_edge, 
                                                                      self.left_arm.arm_angle, 
                                                                      self.right_arm_inner_edge, 
                                                                      self.right_arm.arm_angle, 
                                                                      self.HALFARM_LENGTH, 
                                                                      dist_to_edge)
        #Check if the results are valid
        if hinge_angle == None:
            self.valid = False
            return False

        #Assign the ideal values
        self.hinge_angle                = hinge_angle
        self.hinge_tip_angle            = hinge_inter_angle
        self.hinge_tip                  = hinge_inter_point
        self.ideal_center               = hinge_ideal_origin
        self.ideal_left_arm_inner_edge  = hinge_ideal_left_arm_edge
        self.ideal_right_arm_inner_edge = hinge_ideal_right_arm_edge
        self.ideal_inter_dist           = hinge_ideal_inter_dist

        #Assign the ideal coordinates to local arm parameters
        self.left_arm.ideal_arm_inner_edge  = hinge_ideal_left_arm_edge
        self.right_arm.ideal_arm_inner_edge = hinge_ideal_right_arm_edge

        #Check if the solutions are on the right side of the tip
        left_arm_vertical_vec = hinge_ideal_left_arm_edge - hinge_inter_point
        if np.sum(self.left_arm.vertical_vec*left_arm_vertical_vec) < 0:
            return False 

        return True

    def determine_inner_side_parameters(self):
        '''
        Determine the inner side parameters
        '''
            
        if not self.valid:
            return False 

        left_arm_inner_edge_coord  = self.left_arm_inner_edge
        left_arm_horizontal_vec    = self.left_arm.horizontal_vec

        right_arm_inner_edge_coord = self.right_arm_inner_edge
        right_arm_horizontal_vec   = self.right_arm.horizontal_vec

        self.left_arm_inner_side   = left_arm_inner_edge_coord  + self.DIST_TO_EDGE*left_arm_horizontal_vec
        self.right_arm_inner_side  = right_arm_inner_edge_coord + self.DIST_TO_EDGE*right_arm_horizontal_vec

        (hinge_ideal_left_arm_side, 
         hinge_ideal_right_arm_side, 
         hinge_ideal_origin, 
         hinge_angle, 
         hinge_inter_angle, 
         hinge_inter_point, 
         hinge_ideal_inter_dist) = EMtools.determine_hinge_parameters(self.left_arm_inner_side, 
                                                                      self.left_arm.arm_angle, 
                                                                      self.right_arm_inner_side, 
                                                                      self.right_arm.arm_angle, 
                                                                      self.HALFARM_LENGTH, 
                                                                      0)
        #Check if the results are valid
        if hinge_angle == None:
            self.valid = False
            return False

        #Assign the ideal values
        self.hinge_angle                = hinge_angle
        self.hinge_tip_angle            = hinge_inter_angle
        self.hinge_tip                  = hinge_inter_point
        self.ideal_center               = hinge_ideal_origin
        self.ideal_left_arm_inner_side  = hinge_ideal_left_arm_side
        self.ideal_right_arm_inner_side = hinge_ideal_right_arm_side
        self.ideal_inter_dist           = hinge_ideal_inter_dist

        #Assign the ideal coordinates to local arm parameters
        self.left_arm.ideal_arm_inner_side  = hinge_ideal_left_arm_side
        self.right_arm.ideal_arm_inner_side = hinge_ideal_right_arm_side

        #Check if the solutions are on the right side of the tip
        left_arm_vertical_vec = hinge_ideal_left_arm_side - hinge_inter_point
        if np.sum(self.left_arm.vertical_vec*left_arm_vertical_vec) < 0:
            return False

        return True

    def determine_mask_side_parameters(self):
        '''
        Determine the mask side parameters
        '''
            
        if not self.valid:
            return False 

        left_arm_inner_edge_coord  = self.ideal_left_arm_inner_side
        left_arm_horizontal_vec    = self.left_arm.horizontal_vec

        right_arm_inner_edge_coord = self.ideal_right_arm_inner_side
        right_arm_horizontal_vec   = self.right_arm.horizontal_vec

        #In order to determine mask side locations, we need to move along DNA tether
        #Horizontal vector needs to be corrected by the 1.0/cos(0.5*hinge_angle)
        
        #Determine cos(0.5*hinge_angle)
        self.cos_half_hinge_angle = np.cos(0.5*self.hinge_angle*np.pi/180.0)

        self.left_arm_mask_side   = left_arm_inner_edge_coord  + self.DIST_TO_MASK*left_arm_horizontal_vec/self.cos_half_hinge_angle
        self.right_arm_mask_side  = right_arm_inner_edge_coord + self.DIST_TO_MASK*right_arm_horizontal_vec/self.cos_half_hinge_angle

        (ideal_left_arm_mask_side, 
         ideal_right_arm_mask_side, 
         ideal_mask_origin, 
         hinge_angle, 
         hinge_inter_angle, 
         ideal_mask_inter_point, 
         ideal_mask_inter_dist) = EMtools.determine_hinge_parameters(self.left_arm_mask_side, 
                                                                      self.left_arm.arm_angle, 
                                                                      self.right_arm_mask_side, 
                                                                      self.right_arm.arm_angle, 
                                                                      self.HALFARM_LENGTH, 
                                                                      0)
        #Check if the results are valid
        if hinge_angle == None:
            self.valid = False
            return False

        #Assign the ideal values
        self.ideal_mask_tip             = ideal_mask_inter_point
        self.ideal_mask_center          = 0.5*(self.left_arm_mask_side+self.right_arm_mask_side)
        self.ideal_left_arm_mask_side   = self.left_arm_mask_side
        self.ideal_right_arm_mask_side  = self.right_arm_mask_side
        self.ideal_mask_inter_dist      = EMtools.eular_distance(self.left_arm_mask_side,self.right_arm_mask_side)

        #Assign the ideal coordinates to local arm parameters
        self.left_arm.ideal_arm_mask_side  = self.left_arm_mask_side
        self.right_arm.ideal_arm_mask_side = self.right_arm_mask_side

        return True

    def determine_outer_edge_parameters(self):
        '''
        Determine the inner side parameters
        '''
        if not self.valid:
            return False 

        #1. Outer edge parameters for the optimized coordinates
        self.left_arm_outer_edge        = self.left_arm_inner_edge        - 3.0*self.LINE_TO_LINE_DIST*self.left_arm.horizontal_vec
        self.right_arm_outer_edge       = self.right_arm_inner_edge       - 3.0*self.LINE_TO_LINE_DIST*self.right_arm.horizontal_vec

        self.left_arm.arm_outer_edge    = self.left_arm.arm_inner_edge    - 3.0*self.LINE_TO_LINE_DIST*self.left_arm.horizontal_vec
        self.right_arm.arm_outer_edge   = self.right_arm.arm_inner_edge   - 3.0*self.LINE_TO_LINE_DIST*self.right_arm.horizontal_vec

        #2. Outer edge parameters for the ideal coordinates
        self.ideal_left_arm_outer_edge  = self.ideal_left_arm_inner_edge  - 3.0*self.LINE_TO_LINE_DIST*self.left_arm.horizontal_vec
        self.ideal_right_arm_outer_edge = self.ideal_right_arm_inner_edge - 3.0*self.LINE_TO_LINE_DIST*self.right_arm.horizontal_vec

        self.left_arm.ideal_arm_outer_edge  = self.left_arm.ideal_arm_inner_edge  - 3.0*self.LINE_TO_LINE_DIST*self.left_arm.horizontal_vec
        self.right_arm.ideal_arm_outer_edge = self.right_arm.ideal_arm_inner_edge - 3.0*self.LINE_TO_LINE_DIST*self.right_arm.horizontal_vec

        #3. Optimize ideal outer edge parameter
        result_left  = self.outer_edge_optimization(self.left_arm,   self.LONG_ARM_REFS)
        result_right = self.outer_edge_optimization(self.right_arm,  self.LONG_ARM_REFS)

        #Check the results first
        if len(result_left) == 0 or len(result_right) == 0: return False

        #4. Assign the local arm values to hinge parameters
        self.ideal_left_arm_outer_edge  = self.left_arm.ideal_arm_outer_edge
        self.ideal_right_arm_outer_edge = self.right_arm.ideal_arm_outer_edge

        return True

    def determine_barcode_edge_parameters(self):
        '''
        Determine barcode edge parameters
        '''
        if not self.valid:
            return False 
        
        #1. Assign first to local arm values 
        self.left_arm.arm_barcode_edge  = self.left_arm.ideal_arm_outer_edge  - 1.0*self.LINE_TO_LINE_DIST*self.left_arm.horizontal_vec
        self.right_arm.arm_barcode_edge = self.right_arm.ideal_arm_outer_edge - 1.0*self.LINE_TO_LINE_DIST*self.right_arm.horizontal_vec

        #2. Assign local arm values to hinge parameters
        self.left_arm_barcode_edge  = self.ideal_left_arm_outer_edge  - 1.0*self.LINE_TO_LINE_DIST*self.left_arm.horizontal_vec
        self.right_arm_barcode_edge = self.ideal_right_arm_outer_edge - 1.0*self.LINE_TO_LINE_DIST*self.right_arm.horizontal_vec

        return True

    def outer_edge_optimization(self, left_arm_particle, arm_refs):
        '''
        Optimize the outer edge center
        '''

        left_arm_coord  = left_arm_particle.ideal_arm_outer_edge
        
        left_arm_angle  = left_arm_particle.arm_angle

        horizontal_vec  = left_arm_particle.horizontal_vec

        four_line_angle = left_arm_particle.four_line_angle

        #Get the reference model

        ref_numpy = arm_refs[int((four_line_angle+180)%360)]
        
        #Determine the 4 line points
        distance    = int(0.6*self.LINE_TO_LINE_DIST)
        start_point = left_arm_coord - 0.3*self.LINE_TO_LINE_DIST*horizontal_vec
        line_points = [start_point + i*horizontal_vec for i in range(distance)]

        #Find the best point
        cmp_scores = []

        for new_point in line_points:

            #Check if image is in image bounds
            if EMtools.img_out_of_bounds(self.micrograph_eman, new_point, self.ARM_BOXSIZE):
                continue

            boxim_numpy = EMtools.get_region(self.micrograph_img, new_point[0]-self.ARM_BOXSIZE, new_point[1]-self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE)

            #Mean score
            mean_score  = EMtools.dot_cmp(boxim_numpy,ref_numpy,self.micrograph_mean)

            #Add to list
            cmp_scores.append([new_point[0],new_point[1],mean_score])

        #If there are no results return 0 size array
        if len(cmp_scores) == 0: return []

        cmp_scores     = np.array(cmp_scores)
        sorted_results = cmp_scores[cmp_scores[:,2].argsort()[::-1]]
        best_result    = sorted_results[0,:]

        left_arm_particle.arm_optimization.append({'cmp_scores':cmp_scores, 
                                                   'sorted_results':sorted_results,
                                                   'best_result':best_result,
                                                   'arm_outer_edge':np.array([best_result[0],best_result[1]]),
                                                   'arm_angle':left_arm_angle, 
                                                   'four_line_angle':four_line_angle})

        #Assign the new edge value
        left_arm_particle.ideal_arm_outer_edge = np.array([best_result[0],best_result[1]])

        return best_result


    def determine_barcode_parameters(self):
        '''
        Determine barcode parameters
        '''

        if not self.valid:
            return False 

        #1. Full length barcodes
        self.left_full_barcode          = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.left_arm.four_line_angle),  barcode_relative_coord= 0+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_FULL_LENGTH,  barcode_boxsize=self.FULLBARCODE_BOXSIZE)
        self.right_full_barcode         = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.right_arm.four_line_angle), barcode_relative_coord= 0+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_FULL_LENGTH,  barcode_boxsize=self.FULLBARCODE_BOXSIZE)

        self.left_full_barcode.set_actual_coord(self.left_arm_barcode_edge   ,self.left_arm.vertical_vec)
        self.right_full_barcode.set_actual_coord(self.right_arm_barcode_edge ,self.right_arm.vertical_vec)

        #2. Left short barcodes
        self.left_short_barcode_bottom  = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.left_arm.four_line_angle),  barcode_relative_coord= 60+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_SHORT_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)
        self.left_short_barcode_middle  = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.left_arm.four_line_angle),  barcode_relative_coord=  0+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_SHORT_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)
        self.left_short_barcode_top     = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.left_arm.four_line_angle),  barcode_relative_coord=-60+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_SHORT_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)

        self.left_short_barcode_bottom.set_actual_coord(self.left_arm_barcode_edge ,self.left_arm.vertical_vec)
        self.left_short_barcode_middle.set_actual_coord(self.left_arm_barcode_edge ,self.left_arm.vertical_vec)
        self.left_short_barcode_top.set_actual_coord(self.left_arm_barcode_edge    ,self.left_arm.vertical_vec)

        #3. Left long barcodes
        self.left_long_barcode_bottom   = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.left_arm.four_line_angle),  barcode_relative_coord= 42+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_LONG_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)
        self.left_long_barcode_top      = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.left_arm.four_line_angle),  barcode_relative_coord=-42+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_LONG_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)

        self.left_long_barcode_bottom.set_actual_coord(self.left_arm_barcode_edge ,self.left_arm.vertical_vec)
        self.left_long_barcode_top.set_actual_coord(self.left_arm_barcode_edge    ,self.left_arm.vertical_vec)

        #4. Right short barcodes
        self.right_short_barcode_bottom = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.right_arm.four_line_angle), barcode_relative_coord= 60+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_SHORT_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)
        self.right_short_barcode_middle = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.right_arm.four_line_angle), barcode_relative_coord=  0+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_SHORT_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)
        self.right_short_barcode_top    = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.right_arm.four_line_angle), barcode_relative_coord=-60+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_SHORT_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)

        self.right_short_barcode_bottom.set_actual_coord(self.right_arm_barcode_edge ,self.right_arm.vertical_vec)
        self.right_short_barcode_middle.set_actual_coord(self.right_arm_barcode_edge ,self.right_arm.vertical_vec)
        self.right_short_barcode_top.set_actual_coord(self.right_arm_barcode_edge    ,self.right_arm.vertical_vec)

        #5. Right long barcodes
        self.right_long_barcode_bottom  = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.right_arm.four_line_angle), barcode_relative_coord= 42+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_LONG_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)
        self.right_long_barcode_top     = EMOrigamiBarcode(barcode_num_lines = 2, barcode_angle= int(self.right_arm.four_line_angle), barcode_relative_coord=-42+self.BARCODE_OFFSET  ,barcode_length=self.BARCODE_LONG_LENGTH, barcode_boxsize=self.BARCODE_BOXSIZE)

        self.right_long_barcode_bottom.set_actual_coord(self.right_arm_barcode_edge ,self.right_arm.vertical_vec)
        self.right_long_barcode_top.set_actual_coord(self.right_arm_barcode_edge    ,self.right_arm.vertical_vec)

        #Keep the barcode objects in an array
        self.full_barcodes  = [self.left_full_barcode,self.right_full_barcode]

        self.short_barcodes = [[self.left_short_barcode_bottom,self.right_short_barcode_bottom],
                               [self.left_short_barcode_middle,self.right_short_barcode_middle],
                               [self.left_short_barcode_top   ,self.right_short_barcode_top   ]]

        self.long_barcodes  = [[self.left_long_barcode_bottom,self.right_long_barcode_bottom ],
                               [self.left_long_barcode_top   ,self.right_long_barcode_top    ]]

        return True

    def plot_barcode_locations(self,plot_text=True):
        '''
        Plot barcode locations
        '''
        plt.plot([self.left_short_barcode_bottom.barcode_actual_coord[0]],[self.left_short_barcode_bottom.barcode_actual_coord[1]],'yo',markersize=6)
        plt.plot([self.left_short_barcode_middle.barcode_actual_coord[0]],[self.left_short_barcode_middle.barcode_actual_coord[1]],'yo',markersize=7)
        plt.plot([self.left_short_barcode_top.barcode_actual_coord[0]],   [self.left_short_barcode_top.barcode_actual_coord[1]],'yo',markersize=8)

        plt.plot([self.left_long_barcode_bottom.barcode_actual_coord[0]],[self.left_long_barcode_bottom.barcode_actual_coord[1]],'ys',markersize=7)
        plt.plot([self.left_long_barcode_top.barcode_actual_coord[0]],   [self.left_long_barcode_top.barcode_actual_coord[1]],'ys',markersize=8)

        plt.plot([self.right_short_barcode_bottom.barcode_actual_coord[0]],[self.right_short_barcode_bottom.barcode_actual_coord[1]],'bo',markersize=6)
        plt.plot([self.right_short_barcode_middle.barcode_actual_coord[0]],[self.right_short_barcode_middle.barcode_actual_coord[1]],'bo',markersize=7)
        plt.plot([self.right_short_barcode_top.barcode_actual_coord[0]],   [self.right_short_barcode_top.barcode_actual_coord[1]],'bo',markersize=8)

        plt.plot([self.right_long_barcode_bottom.barcode_actual_coord[0]],[self.right_long_barcode_bottom.barcode_actual_coord[1]],'bs',markersize=7)
        plt.plot([self.right_long_barcode_top.barcode_actual_coord[0]],   [self.right_long_barcode_top.barcode_actual_coord[1]],'bs',markersize=8)

        #Print the long barcode barcode scores
        if plot_text:
            plt.text(self.left_long_barcode_bottom.barcode_actual_coord[0],self.left_long_barcode_bottom.barcode_actual_coord[1],'%.2f'%(self.left_long_barcode_bottom.barcode_score),color='k',fontsize=12)
            plt.text(self.left_long_barcode_top.barcode_actual_coord[0],self.left_long_barcode_top.barcode_actual_coord[1],'%.2f'%(self.left_long_barcode_top.barcode_score),color='k',fontsize=12)

            plt.text(self.right_long_barcode_bottom.barcode_actual_coord[0],self.right_long_barcode_bottom.barcode_actual_coord[1],'%.2f'%(self.right_long_barcode_bottom.barcode_score),color='k',fontsize=12)
            plt.text(self.right_long_barcode_top.barcode_actual_coord[0],self.right_long_barcode_top.barcode_actual_coord[1],'%.2f'%(self.right_long_barcode_top.barcode_score),color='k',fontsize=12)

        #Plot the detected barcodes
        detected_long_barcode  = self.long_barcodes[self.long_barcode_index][self.long_barcode_side]
        detected_short_barcode = self.short_barcodes[self.short_barcode_index][1-self.long_barcode_side]

        plt.plot([detected_long_barcode.barcode_actual_coord[0]],[detected_long_barcode.barcode_actual_coord[1]],'ks',markersize=8)
        plt.plot([detected_short_barcode.barcode_actual_coord[0]],[detected_short_barcode.barcode_actual_coord[1]],'ko',markersize=8)

    def check_hinge_quality(self):
        '''
        Check hinge quality based on the optimization results
        '''

        #Check if the hinge is valid
        if not self.valid:
            return False

        #Check the displacements of left and right arm origins
        self.left_arm_displacement  = EMtools.eular_distance(self.ideal_left_arm_inner_edge  , self.left_arm.arm_center)
        self.right_arm_displacement = EMtools.eular_distance(self.ideal_right_arm_inner_edge , self.right_arm.arm_center)

        #Determine hinge quality score
        left_diff   = EMtools.eular_distance(self.left_arm.arm_inner_edge , self.ideal_left_arm_inner_edge)
        right_diff  = EMtools.eular_distance(self.right_arm.arm_inner_edge, self.ideal_right_arm_inner_edge)
        center_diff = np.fabs(EMtools.eular_distance(self.left_arm.arm_inner_edge , self.right_arm.arm_inner_edge) - self.ideal_inter_dist)
        
        ideal_half_arm_length = EMtools.eular_distance(self.ideal_left_arm_inner_edge,self.hinge_tip)

        #Calculate ideal comparison score
        self.hinge_quality_score = (np.sqrt(left_diff**2 + right_diff**2 + center_diff**2)/
                                    np.sqrt(2*ideal_half_arm_length**2 + self.ideal_inter_dist**2))

        #Check if the new ideal points are within image bounds
        if (EMtools.img_out_of_bounds(self.micrograph_eman, self.ideal_left_arm_inner_edge, self.HALFARM_LENGTH) or
            EMtools.img_out_of_bounds(self.micrograph_eman, self.ideal_right_arm_inner_edge, self.HALFARM_LENGTH) or
            EMtools.img_out_of_bounds(self.micrograph_eman, self.hinge_tip, self.HALFARM_LENGTH)):
            self.img_out_of_bounds = True

        #Displacement comparisons
        if self.left_arm_displacement > self.ARM_LENGTH or self.right_arm_displacement > self.ARM_LENGTH:
            self.bad_displacement = True

        #Quality score
        if np.isnan(self.hinge_quality_score) or self.hinge_quality_score > self.HINGE_QUAL_CUTOFF:
            self.bad_quality_score = True
        
        #Cmp score
        if np.isnan(self.hinge_cmp_score) or self.hinge_cmp_score < self.HINGE_CMP_CUTOFF:
            self.bad_cmp_score = True

        #Hinge angle
        if self.hinge_angle < self.HINGE_ANGLE_CUTOFF['low'] or self.hinge_angle > self.HINGE_ANGLE_CUTOFF['high']:
            self.bad_angle = True

        #Final verdict
        if self.bad_displacement or self.bad_angle or self.bad_quality_score or self.bad_cmp_score:
            self.valid = False

        return True

    def copy_hinge_to_arm(self):
        '''
        Copy hinge coordinates to arm coordinates
        '''
        self.left_arm.arm_inner_edge  = self.left_arm_inner_edge.copy()
        self.right_arm.arm_inner_edge = self.right_arm_inner_edge.copy()

    def arm_optimization_first_step(self, left_arm_particle, right_arm_particle, arm_width_span, arm_length_span, line_refs, small_arm_refs):
        '''
        Optimize detection for arm features
        '''

        #Initialize optimization vector
        left_arm_particle.arm_optimization = []

        #Initialize optimization dictionary
        left_arm_particle.arm_valid = True

        #Retreive the paramaters from arm particles
        left_arm_center = left_arm_particle.arm_center
        left_arm_angle  = left_arm_particle.arm_angle

        right_arm_center= right_arm_particle.arm_center
        right_arm_angle = right_arm_particle.arm_angle

        #Left arm line and vector
        left_arm_line  = EMtools.make_line(left_arm_center, left_arm_angle+90, self.ARM_WIDTH_SPAN)
        left_arm_vec   = left_arm_line[1,:] - left_arm_line[0,:]

        left_to_right  = right_arm_center - left_arm_center
        
        #Determine left vector direction
        dot_angle = EMtools.ray_angle(left_arm_vec, left_to_right)
        
        #Check if the angle determination worked well
        if dot_angle == None:
            return []
        
        if dot_angle > 90:
            left_arm_vec  = -left_arm_vec
            left_arm_line =  left_arm_line[::-1,:]

        #Normalize vector
        left_arm_vec = EMtools.normalize_vector(left_arm_vec)

        #Make the right and left line points
        left_x  = np.around(np.linspace(left_arm_line[0,0], left_arm_line[0,0]+1.0*left_arm_vec[0]*self.ARM_WIDTH_SPAN, self.ARM_WIDTH_SPAN))
        left_y  = np.around(np.linspace(left_arm_line[0,1], left_arm_line[0,1]+1.0*left_arm_vec[1]*self.ARM_WIDTH_SPAN, self.ARM_WIDTH_SPAN))

        #make the points
        left_points  = np.unique(zip(left_x,left_y),axis=0)

        #Score and coordinates container
        left_scores  = []
        left_centers = []

        #1. Make an intensity profile perpendicular to arm
        for i in range(len(left_points)):
            left_center = left_points[i,:]
            
            #Check if image is in image bounds
            if EMtools.img_out_of_bounds(self.micrograph_eman, left_center, self.ARM_HALF_BOXSIZE):
                continue
            
            #Crop particle from image
            boxim_numpy = EMtools.get_region(self.micrograph_img, left_center[0]-self.ARM_HALF_BOXSIZE, left_center[1]-self.ARM_HALF_BOXSIZE, self.ARM_BOXSIZE, self.ARM_BOXSIZE)

            #Get reference
            ref_numpy   = line_refs[int(left_arm_angle%180)]

            #Mean score
            mean_score  = EMtools.dot_cmp(boxim_numpy,ref_numpy,self.micrograph_mean)

            left_scores.append(mean_score)
            left_centers.append([left_center[0],left_center[1]])

        #If there are no results return 0 size array
        if len(left_centers) == 0: return []

        max_left_centers = scipy.signal.argrelmax(np.array(left_scores),order=self.MAXIMA_WINDOW)[0]
        left_scores      = np.array(left_scores)
        left_centers     = np.array(left_centers)

        #Store optimization results
        left_arm_particle.arm_optimization.append({'scores':left_scores[max_left_centers],
                                                   'centers':left_centers[max_left_centers]})
        
        #2. Determine the best matching lines along the profile
        nRef    = len(line_refs)
        results = []

        for i in range(len(max_left_centers)):
            maxim_index  = max_left_centers[i]
            maxim_score  = left_scores[maxim_index]
            maxim_center = left_centers[maxim_index]
            
            #Measure distance with respect to arm center (Outer direction is 1)
            distance    = EMtools.eular_distance(np.array(maxim_center),np.array(left_arm_center))
            direction   = 1

            #Maximum center difference
            maxim_center_diff = maxim_center-left_arm_center

            #Determine direction angle
            direction_angle = EMtools.ray_angle(maxim_center_diff, left_arm_vec)

            #Check if direction angle worked fine
            if direction_angle == None:
                continue

            if maxim_center_diff.dot(maxim_center_diff) > 0 and direction_angle > 90:
                direction = -1

            #Correct the distance
            distance   = distance*direction 

            #Crop the particle centered around the maxim center
            boxim_numpy = EMtools.get_region(self.micrograph_img, maxim_center[0]-self.ARM_HALF_BOXSIZE, maxim_center[1]-self.ARM_HALF_BOXSIZE, self.ARM_BOXSIZE, self.ARM_BOXSIZE)

            #Store the comparison socres
            cmp_scores  = []

            #Compare with 180 angle models
            for j in range(180):
                ref_numpy = line_refs[j]

                #Mean score
                mean_score  = EMtools.dot_cmp(boxim_numpy,ref_numpy,self.micrograph_mean)
                
                cmp_scores.append(mean_score)

            #Convert the comparison scores to numpy
            cmp_scores = np.array(cmp_scores)
            max_index  = np.argmax(cmp_scores)

            #Determine the angle difference with respect to arm angle
            angle_diff = np.fabs(left_arm_angle%180 - max_index)
            if angle_diff > 90: angle_diff = np.fabs(angle_diff - 180)

            #Add the best results to the results list
            if angle_diff < self.ARM_ANGLE_CUTOFF and cmp_scores[max_index] > 0 and np.sum(cmp_scores) > 0:
                results.append([maxim_center, maxim_score, max_index, angle_diff, distance, cmp_scores[max_index]])
        
        #Convert results to numpy array
        results = np.array(results)
        
        #If there are no results return 0 size array
        if len(results) == 0: return []

        #Sort results based on positive (towards outer) and negative (towards inner) distance
        pos_valid           = np.nonzero(results[:,4] > 0)[0]
        pos_results         = results[pos_valid,:]
        sorted_pos_results  = pos_results[pos_results[:,4].argsort()]  

        minus_valid          = np.nonzero(results[:,4] < 0)[0]
        minus_results        = results[minus_valid,:]
        sorted_minus_results = minus_results[minus_results[:,4].argsort()[::-1]]

        #Ideally pick a close positive result, if not pick the negative result
        if len(sorted_pos_results) > 0 and sorted_pos_results[0,4] < self.LINE_TO_LINE_DIST:
            best_result  = sorted_pos_results[0,:]
        elif len(sorted_minus_results) > 0:
            best_result  = sorted_minus_results[0,:]
        else:
            left_arm_particle.arm_valid = False
            self.valid                  = False
            return []

        #Determine the maximum score and normalize the scores
        best_score     = np.max(results[:,5])
        best_i         = np.argmax(results[:,5])
        best_angle     = results[best_i,2]

        #Normalize the scores
        results[:,5] = 100.0*results[:,5]/best_score
        results[:,1] = 100.0*results[:,1]/best_score

        #Sort results based on the comparison score
        sorted_score_results = results[results[:,5].argsort()[::-1]]

        #Get the best result location
        best_result_coord  = best_result[0]

        #Adjust the angles
        #1. First subtract the best angle
        sorted_score_results[:,2]     -= best_angle
        
        #2. If the difference is higher than 90, add 180
        valid = np.nonzero(np.abs(sorted_score_results[:,2]) > 90)[0]
        sorted_score_results[valid,2] += 180

        #3. Add back the best angle
        sorted_score_results[:,2]     += best_angle
        
        #4. If the angle is higher than 180, subtract 180
        valid = np.nonzero(sorted_score_results[:,2] > 180)[0]
        sorted_score_results[:,2]     -= 180

        best_avg_angle     = int(round(np.mean(sorted_score_results[:2,2])%180))

        #Update left arm vec
        left_arm_line  = EMtools.make_line(left_arm_center, best_avg_angle+90, self.ARM_WIDTH_SPAN)
        left_arm_vec   = left_arm_line[1,:] - left_arm_line[0,:]

        #Determine left-vector directions
        dot_angle = EMtools.ray_angle(left_arm_vec, left_to_right)

        #Check if the angle determination worked fine
        if dot_angle == None:
            return []

        if dot_angle > 90:
            left_arm_vec  = -left_arm_vec
            left_arm_line =  left_arm_line[::-1,:]

        #Normalize vector
        left_arm_vec = EMtools.normalize_vector(left_arm_vec)

        #Store the optimization results
        left_arm_particle.arm_optimization.append({'arm_angle':best_avg_angle, 
                                                   'positive_results':sorted_pos_results, 
                                                   'minus_results':sorted_minus_results, 
                                                   'score_results':sorted_score_results})

        #3. Determine the inner arm edge
        previous_coord     = best_result_coord
        four_line_distance = EMtools.eular_distance(previous_coord, left_arm_center)
        four_line_angle    = best_avg_angle

        #Make a unit vector oriented at best angle
        four_line_vec      = EMtools.unit_vector(best_avg_angle)

        #Determine the orientation of the 4-line aligner
        cross_product = np.cross(left_arm_vec, four_line_vec)
        if cross_product < 0: four_line_angle = best_avg_angle + 180

        #Store the scores
        cmp_scores = []

        #Finally, Move in the left arm vec direction to find a better solution than the existing one
        line_counter = 0 
        
        #Get the reference image
        ref_numpy = small_arm_refs[four_line_angle]

        while four_line_distance < 1.0*self.HALFARM_WIDTH:
            
            #Update the coordinate
            new_coord  = previous_coord+line_counter*self.LINE_TO_LINE_DIST*left_arm_vec

            four_line_distance = EMtools.eular_distance(new_coord,left_arm_center) 

            #Update line counter
            line_counter += 1

            #Check if image is in image bounds
            if EMtools.img_out_of_bounds(self.micrograph_eman, new_coord, self.ARM_BOXSIZE):
                continue
            
            #If the distance is higher than the limit exit the loop
            if four_line_distance > 2*self.HALFARM_WIDTH:
                break

            #Crop the image
            boxim_numpy = EMtools.get_region(self.micrograph_img, new_coord[0]-self.ARM_BOXSIZE, new_coord[1]-self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE)

            #Mean score
            mean_score  = EMtools.dot_cmp(boxim_numpy,ref_numpy,self.micrograph_mean)
            
            #Add the scores to comparison list
            cmp_scores.append([new_coord[0], new_coord[1], mean_score])


        #If there are no results return 0 size array
        if len(cmp_scores) == 0: return []

        #Sort results and get the best one
        cmp_scores     = np.array(cmp_scores)
        sorted_results = cmp_scores[cmp_scores[:,2].argsort()[::-1]]
        best_result    = np.array([sorted_results[0][0],sorted_results[0][1], best_avg_angle, sorted_results[0][2]])
        
        #Add the results to arm optimization dictionary 
        left_arm_particle.arm_optimization.append({'cmp_scores':cmp_scores,
                                                   'sorted_results':sorted_results,
                                                   'best_result':best_result,
                                                   'arm_inner_edge':np.array([best_result[0],best_result[1]]),
                                                   'arm_angle':best_avg_angle,
                                                   'four_line_angle':four_line_angle})
        
        #Determine vertical vector
        hinge_inter_point = EMtools.intersection_point(np.array([best_result[0],best_result[1]]),best_result[2],right_arm_center,right_arm_angle)
        
        #Check if intersection returned a valid answer
        if hinge_inter_point[0] == None: return []

        vertical_vec      = np.array([best_result[0],best_result[1]]) - hinge_inter_point
        vertical_vec      = EMtools.normalize_vector(vertical_vec)

        #Assign the final results
        left_arm_particle.arm_inner_edge  = np.array([best_result[0],best_result[1]])
        left_arm_particle.arm_angle       = best_avg_angle
        left_arm_particle.four_line_angle = four_line_angle
        left_arm_particle.horizontal_vec  = left_arm_vec
        left_arm_particle.vertical_vec    = vertical_vec

        #Assign the new results
        return best_result

    def arm_optimization_angle(self, left_arm_particle, arm_refs, angle_range=20):
        '''
        Final angle optimization
        '''
        #Get the coordinates
        left_arm_coord  = left_arm_particle.arm_inner_edge

        four_line_angle = left_arm_particle.four_line_angle

        if EMtools.img_out_of_bounds(self.micrograph_eman, left_arm_coord, self.ARM_BOXSIZE):
            self.valid = False
            return []

        #Crop the box
        boxim_numpy = EMtools.get_region(self.micrograph_img, left_arm_coord[0]-self.ARM_BOXSIZE, left_arm_coord[1]-self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE)

        #Possible angles for optimization
        possible_angles = np.arange(four_line_angle-0.5*angle_range,four_line_angle+0.5*angle_range)

        #Store the results
        cmp_scores      = []
        
        for possible_angle in possible_angles:
            #Reference image
            ref_numpy = arm_refs[int(possible_angle%360)]

            #Mean score
            mean_score  = EMtools.dot_cmp(boxim_numpy,ref_numpy,self.micrograph_mean)

            #Add to list
            cmp_scores.append([possible_angle,mean_score])

        #Comparison results
        cmp_scores     = np.array(cmp_scores)
        sorted_results = cmp_scores[cmp_scores[:,1].argsort()[::-1]]
        best_result    = sorted_results[0,:]
        best_angle     = int(best_result[0]%360)

        left_arm_particle.arm_optimization.append({'cmp_scores':cmp_scores, 
                                                   'sorted_results':sorted_results,
                                                   'best_result':best_result,
                                                   'arm_angle':best_angle,
                                                   'four_line_angle':best_angle})

        #Assign the final results
        left_arm_particle.arm_angle        = best_angle
        left_arm_particle.four_line_angle  = best_angle

        return best_result

    def determine_the_vectors(self):
        '''
        Determine the horizontal and vertical vectors from angle values
        '''
        (hinge_left_arm_coord, 
         hinge_right_arm_coord, 
         hinge_origin, 
         hinge_angle, 
         hinge_inter_angle, 
         hinge_inter_point, 
         hinge_ideal_inter_dist) = EMtools.determine_hinge_parameters(self.left_arm.arm_inner_edge, 
                                                                      self.left_arm.arm_angle, 
                                                                      self.right_arm.arm_inner_edge, 
                                                                      self.right_arm.arm_angle, 
                                                                      self.HALFARM_LENGTH, 
                                                                      self.DIST_TO_EDGE)
         #Check if hinge parameters are determined properly
        if hinge_angle == None:
            self.valid = False
            return False
        
        #Determine vertical vectors
        self.left_arm.vertical_vec   = EMtools.normalize_vector(self.left_arm.arm_inner_edge  - hinge_inter_point)
        self.right_arm.vertical_vec  = EMtools.normalize_vector(self.right_arm.arm_inner_edge - hinge_inter_point)

        #1. Determine left arm horizontal vector
        left_arm_horizontal_vec = EMtools.unit_vector(self.left_arm.arm_angle+90)
        cross_product           = np.cross(self.left_arm.vertical_vec,left_arm_horizontal_vec)

        if cross_product < 0: left_arm_horizontal_vec = - left_arm_horizontal_vec

        #2. Determine right arm horizontal vector

        right_arm_horizontal_vec = EMtools.unit_vector(self.right_arm.arm_angle+90)
        cross_product            = np.cross(self.right_arm.vertical_vec,right_arm_horizontal_vec)

        if cross_product > 0: right_arm_horizontal_vec = - right_arm_horizontal_vec
         
        #Assign the vectors
        self.left_arm.horizontal_vec  = left_arm_horizontal_vec
        self.right_arm.horizontal_vec = right_arm_horizontal_vec
        
        return True
    def arm_optimization_vertical(self, left_arm_particle, right_arm_particle, arm_refs, point_separation=5):

        #Retrieve previous optimization results
        left_arm_coord  = left_arm_particle.arm_inner_edge
        right_arm_coord = right_arm_particle.arm_inner_edge
        
        left_arm_angle  = left_arm_particle.arm_angle
        right_arm_angle = right_arm_particle.arm_angle

        horizontal_vec  = left_arm_particle.horizontal_vec
        vertical_vec    = left_arm_particle.vertical_vec

        four_line_angle = left_arm_particle.four_line_angle

        #Determine hinge parameters
        (hinge_left_arm_coord, 
         hinge_right_arm_coord, 
         hinge_origin, 
         hinge_angle, 
         hinge_inter_angle, 
         hinge_inter_point, 
         hinge_ideal_inter_dist) = EMtools.determine_hinge_parameters(left_arm_coord, 
                                                                      left_arm_angle, 
                                                                      right_arm_coord, 
                                                                      right_arm_angle, 
                                                                      self.HALFARM_LENGTH, 
                                                                      self.DIST_TO_EDGE)
        #Check if hinge parameters are determined properly
        if hinge_angle == None:
            self.valid = False
            return []

        #Start point
        start_point     = hinge_inter_point
        
        #Possible points
        possible_points = EMtools.make_vector_points(start_point,vertical_vec,int(3.0*self.HALFARM_LENGTH))

        #Get the reference model
        ref_numpy = arm_refs[int(four_line_angle)]

        #Keep the results
        cmp_scores = []

        for new_coord in possible_points[::point_separation,:]:

            #Check if image is in image bounds
            if EMtools.img_out_of_bounds(self.micrograph_eman, new_coord, self.ARM_BOXSIZE):
                cmp_scores.append([new_coord[0],new_coord[1],0.0])
                continue

            boxim_numpy = EMtools.get_region(self.micrograph_img, new_coord[0]-self.ARM_BOXSIZE, new_coord[1]-self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE)
           
            #Mean score
            mean_score  = EMtools.dot_cmp(boxim_numpy,ref_numpy,self.micrograph_mean)

            #Add to list
            cmp_scores.append([new_coord[0],new_coord[1],mean_score])

        #If there are no results return 0 size array
        if len(cmp_scores) == 0: return []

        cmp_scores     = np.array(cmp_scores)
        sorted_results = cmp_scores[cmp_scores[:,2].argsort()[::-1]]
        best_result    = sorted_results[0,:]
        
        left_arm_particle.arm_optimization.append({'cmp_scores':cmp_scores, 
                                                   'sorted_results':sorted_results,
                                                   'best_result':best_result,
                                                   'arm_inner_edge':np.array([best_result[0],best_result[1]]),
                                                   'arm_angle':left_arm_angle,
                                                   'four_line_angle':four_line_angle,
                                                   'hinge_tip':hinge_inter_point})

        #Assign the final results
        left_arm_particle.arm_inner_edge  = np.array([best_result[0],best_result[1]])

        return best_result

    def arm_optimization_horizontal_coarse(self,left_arm_particle, arm_refs):
        '''
        Perform post vertical arm optimization
        '''
        left_arm_coord  = left_arm_particle.arm_inner_edge
        
        left_arm_angle  = left_arm_particle.arm_angle

        horizontal_vec  = left_arm_particle.horizontal_vec

        four_line_angle = left_arm_particle.four_line_angle

        #Get the reference model
        ref_numpy = arm_refs[int(four_line_angle)]
        
        #Determine the 4 line points
        line_points = [left_arm_coord -   self.LINE_TO_LINE_DIST*horizontal_vec, 
                       left_arm_coord, 
                       left_arm_coord +   self.LINE_TO_LINE_DIST*horizontal_vec]

        #Find the best point
        cmp_scores = []

        for new_point in line_points:

            #Check if image is in image bounds
            if EMtools.img_out_of_bounds(self.micrograph_eman, new_point, self.ARM_BOXSIZE):
                continue

            boxim_numpy = EMtools.get_region(self.micrograph_img, new_point[0]-self.ARM_BOXSIZE, new_point[1]-self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE)

            #Mean score
            mean_score  = EMtools.dot_cmp(boxim_numpy,ref_numpy,self.micrograph_mean)

            #Add to list
            cmp_scores.append([new_point[0],new_point[1],mean_score])

        #If there are no results return 0 size array
        if len(cmp_scores) == 0: return []

        cmp_scores     = np.array(cmp_scores)
        sorted_results = cmp_scores[cmp_scores[:,2].argsort()[::-1]]
        best_result    = sorted_results[0,:]

        left_arm_particle.arm_optimization.append({'cmp_scores':cmp_scores, 
                                                   'sorted_results':sorted_results,
                                                   'best_result':best_result,
                                                   'arm_inner_edge':np.array([best_result[0],best_result[1]]),
                                                   'arm_angle':left_arm_angle, 
                                                   'four_line_angle':four_line_angle})

        #Assign the new edge value
        left_arm_particle.arm_inner_edge = np.array([best_result[0],best_result[1]])
        
        return best_result

    def arm_optimization_horizontal_fine(self, left_arm_particle, arm_refs):
        '''
        Perform fine optimization for the edge detection
        '''
        left_arm_coord  = left_arm_particle.arm_inner_edge
        
        left_arm_angle  = left_arm_particle.arm_angle

        horizontal_vec  = left_arm_particle.horizontal_vec

        four_line_angle = left_arm_particle.four_line_angle

        #Get the reference model
        ref_numpy = arm_refs[int(four_line_angle)]
        
        #Determine the 4 line points
        distance    = int(2.2*self.LINE_TO_LINE_DIST)
        start_point = left_arm_coord - 1.1*self.LINE_TO_LINE_DIST*horizontal_vec
        line_points = [start_point + i*horizontal_vec for i in range(distance)]

        #Find the best point
        cmp_scores = []

        for new_point in line_points:

            #Check if image is in image bounds
            if EMtools.img_out_of_bounds(self.micrograph_eman, new_point, self.ARM_BOXSIZE):
                continue

            boxim_numpy = EMtools.get_region(self.micrograph_img, new_point[0]-self.ARM_BOXSIZE, new_point[1]-self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE)

            #Mean score
            mean_score  = EMtools.dot_cmp(boxim_numpy,ref_numpy,self.micrograph_mean)

            #Add to list
            cmp_scores.append([new_point[0],new_point[1],mean_score])

        #If there are no results return 0 size array
        if len(cmp_scores) == 0: return []

        cmp_scores     = np.array(cmp_scores)
        sorted_results = cmp_scores[cmp_scores[:,2].argsort()[::-1]]
        best_result    = sorted_results[0,:]

        left_arm_particle.arm_optimization.append({'cmp_scores':cmp_scores, 
                                                   'sorted_results':sorted_results,
                                                   'best_result':best_result,
                                                   'arm_inner_edge':np.array([best_result[0],best_result[1]]),
                                                   'arm_angle':left_arm_angle, 
                                                   'four_line_angle':four_line_angle})

        #Assign the new edge value
        left_arm_particle.arm_inner_edge = np.array([best_result[0],best_result[1]])

        return best_result

    def optimize_hinge_center(self):
        '''
        Final hinge optimization. Comes after a vertical optimization
        '''
        #Retrieve the results
        left_scores  = self.left_arm.arm_optimization[-1]['cmp_scores']
        right_scores = self.right_arm.arm_optimization[-1]['cmp_scores']

        left_angle   = self.left_arm.arm_optimization[-1]['arm_angle']
        right_angle  = self.right_arm.arm_optimization[-1]['arm_angle']

        #Combine the scores and find the maximum
        total_scores = 0.5*(left_scores[:,2] + right_scores[:,2])

        max_i = np.argmax(total_scores)

        #Determine the final results
        hinge_left_arm_coord  = np.array([left_scores[max_i,0]  ,left_scores[max_i,1]])
        hinge_right_arm_coord = np.array([right_scores[max_i,0] ,right_scores[max_i,1]])

        #Store final results
        self.hinge_optimization.append({'left_arm_inner_edge': hinge_left_arm_coord,
                                        'left_arm_angle':left_angle,
                                        'right_arm_inner_edge': hinge_right_arm_coord,
                                        'right_arm_angle': right_angle,
                                        'cmp_score':total_scores[max_i],
                                        'total_scores':total_scores,
                                        'left_scores':left_scores,
                                        'right_scores':right_scores})

        #Update edge coordinates
        self.left_arm_inner_edge  = hinge_left_arm_coord
        self.right_arm_inner_edge = hinge_right_arm_coord

        #Hinge comparison score
        self.hinge_cmp_score      = total_scores[max_i]

        return hinge_left_arm_coord, hinge_right_arm_coord, total_scores[max_i]

    def calculate_barcode_scores(self):
        '''
        Calculate barcode scores
        '''
        if not self.valid:
            return False 

        self.long_barcode_scores  = np.zeros((2,2))
        self.short_barcode_scores = np.zeros((3,2))

        #Calculate short barcode scores
        self.short_barcode_scores[0,0] = self.left_short_barcode_bottom.determine_barcode_score(self.micrograph_img, self.SHORT_BARCODE_REFS, self.micrograph_mean)
        self.short_barcode_scores[1,0] = self.left_short_barcode_middle.determine_barcode_score(self.micrograph_img, self.SHORT_BARCODE_REFS, self.micrograph_mean)
        self.short_barcode_scores[2,0] = self.left_short_barcode_top.determine_barcode_score(self.micrograph_img,    self.SHORT_BARCODE_REFS, self.micrograph_mean)
        
        self.short_barcode_scores[0,1] = self.right_short_barcode_bottom.determine_barcode_score(self.micrograph_img, self.SHORT_BARCODE_REFS, self.micrograph_mean)
        self.short_barcode_scores[1,1] = self.right_short_barcode_middle.determine_barcode_score(self.micrograph_img, self.SHORT_BARCODE_REFS, self.micrograph_mean)
        self.short_barcode_scores[2,1] = self.right_short_barcode_top.determine_barcode_score(self.micrograph_img   , self.SHORT_BARCODE_REFS, self.micrograph_mean)
        
        #Calculate long barcode scores
        self.long_barcode_scores[0,0] = self.left_long_barcode_bottom.determine_barcode_score(self.micrograph_img, self.LONG_BARCODE_REFS, self.micrograph_mean)
        self.long_barcode_scores[1,0] = self.left_long_barcode_top.determine_barcode_score(self.micrograph_img,    self.LONG_BARCODE_REFS, self.micrograph_mean)
        
        self.long_barcode_scores[0,1] = self.right_long_barcode_bottom.determine_barcode_score(self.micrograph_img, self.LONG_BARCODE_REFS, self.micrograph_mean)
        self.long_barcode_scores[1,1] = self.right_long_barcode_top.determine_barcode_score(self.micrograph_img, self.LONG_BARCODE_REFS   , self.micrograph_mean)
     
        return True

    def optimize_barcode_location(self,left_arm_particle, barcode_refs):
        '''
        Optimize barcode location
        '''

        #Retrieve relevant parameters
        left_arm_coord  = left_arm_particle.arm_barcode_edge

        four_line_angle = left_arm_particle.four_line_angle

        vertical_vec    = left_arm_particle.vertical_vec

        #Start point
        start_point     = left_arm_coord - 1.1*self.HALFARM_LENGTH*vertical_vec 
        
        #Possible points
        possible_points = EMtools.make_vector_points(start_point,vertical_vec,int(2.2*self.HALFARM_LENGTH))

        #Get the reference model
        ref_numpy = barcode_refs[int(four_line_angle)]

        #Keep the results
        cmp_scores = []

        for new_coord in possible_points[::point_separation,:]:

            #Check if image is in image bounds
            if EMtools.img_out_of_bounds(self.micrograph_eman, new_coord, self.ARM_BOXSIZE):
                cmp_scores.append([new_coord[0],new_coord[1],0.0])
                continue

            boxim_numpy = EMtools.get_region(self.micrograph_img, new_coord[0]-self.ARM_BOXSIZE, new_coord[1]-self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE, 2*self.ARM_BOXSIZE)
           
            #Mean score
            mean_score  = EMtools.dot_cmp(boxim_numpy,ref_numpy,self.micrograph_mean)

            #Add to list
            cmp_scores.append([new_coord[0],new_coord[1],mean_score])

        #If there are no results return 0 size array
        if len(cmp_scores) == 0: return []

        cmp_scores     = np.array(cmp_scores)
        sorted_results = cmp_scores[cmp_scores[:,2].argsort()[::-1]]
        best_result    = sorted_results[0,:]
        
        left_arm_particle.arm_optimization.append({'cmp_scores':cmp_scores, 
                                                   'sorted_results':sorted_results,
                                                   'best_result':best_result,
                                                   'barcode_location':np.array([best_result[0],best_result[1]]),
                                                   'four_line_angle':four_line_angle})

        return best_result

    def assign_barcode_state(self):
        '''
        Assign the barcode state 
        '''

        #1. Determine the long barcode side
        self.normalized_long_barcode_scores = self.long_barcode_scores/np.min(self.long_barcode_scores,axis=0) 
        self.long_barcode_side              = np.argmax(np.sum(self.normalized_long_barcode_scores,axis=0))

        #2. Determine the best location for long barcode
        self.long_barcode_index = np.argmax(self.normalized_long_barcode_scores[:,self.long_barcode_side])
        
        #3. If the long barcode is top, conside only the bottom and middle short bits on the opposite side
        if self.long_barcode_index == 1:
            self.short_barcode_index = np.argmax(self.short_barcode_scores[:2,1-self.long_barcode_side])
        else:
            self.short_barcode_index = np.argmax(self.short_barcode_scores[:, 1-self.long_barcode_side])

        #4. Assign the raw scores
        self.long_barcode_cmp_score  = self.long_barcode_scores[self.long_barcode_index  ,  self.long_barcode_side]
        self.short_barcode_cmp_score = self.short_barcode_scores[self.short_barcode_index,1-self.long_barcode_side]

        #5. Assign hinge barcode class and flip parameter
        self.hinge_barcode_class     = self.long_barcode_index*len(self.short_barcode_scores[:,0])+self.short_barcode_index
        self.hinge_barcode_flip      = self.long_barcode_side

        print(self.left_arm.particle_id, self.right_arm.particle_id, self.long_barcode_side, self.long_barcode_index, self.short_barcode_index, self.hinge_barcode_class, self.hinge_barcode_flip)
        

    def optimization_protocol(self):
        '''
        Complete optimization protocol for left and right arm detection
        '''
        #1. First step in arm detection (horizontal)
        result_left = self.arm_optimization_first_step(self.left_arm,  
                                         self.right_arm, 
                                         self.ARM_WIDTH_SPAN, 
                                         self.ARM_LENGTH_SPAN, 
                                         self.LINE_REFS, 
                                         self.SMALL_ARM_REFS)
        
        result_right = self.arm_optimization_first_step(self.right_arm, 
                                         self.left_arm,  
                                         self.ARM_WIDTH_SPAN, 
                                         self.ARM_LENGTH_SPAN, 
                                         self.LINE_REFS, 
                                         self.SMALL_ARM_REFS)
        
        #Check results
        if len(result_left) == 0 or len(result_right) == 0: return False

        #Continue if both arms are valid
        if self.left_arm.arm_valid and self.right_arm.arm_valid:
            
            #21. Perform an angle optimization
            result_left  = self.arm_optimization_angle(self.left_arm,  self.SMALL_ARM_REFS, angle_range=20)
            result_right = self.arm_optimization_angle(self.right_arm, self.SMALL_ARM_REFS, angle_range=20)

            #Check results
            if len(result_left) == 0 or len(result_right) == 0: return False

            #22. Determine the horizontal and vertical vectors
            vector_result = self.determine_the_vectors()
            if not vector_result: return False

            #23. Second perform vertical optimization
            result_left = self.arm_optimization_vertical(self.left_arm,  
                                           self.right_arm, 
                                           self.LONG_ARM_REFS,
                                           point_separation=2)
            
            result_right = self.arm_optimization_vertical(self.right_arm, 
                                           self.left_arm,  
                                           self.LONG_ARM_REFS,
                                           point_separation=2)
            
            #Check results
            if len(result_left) == 0 or len(result_right) == 0: return False

            #24. Merge the independent left and right arm vertical optimizations 
            self.optimize_hinge_center()

            #25. Pass hinge coordinates to arm coordinates
            self.copy_hinge_to_arm()

            #31. Second round of horizontal optimization
            result_left  = self.arm_optimization_horizontal_coarse(self.left_arm,  self.SMALL_ARM_REFS)
            result_right = self.arm_optimization_horizontal_coarse(self.right_arm, self.SMALL_ARM_REFS)

            #Check results
            if len(result_left) == 0 or len(result_right) == 0: return False

            #32. Second round of vertical optimization
            result_left  = self.arm_optimization_vertical(self.left_arm,  self.right_arm, self.LONG_ARM_REFS, point_separation=2)
            result_right = self.arm_optimization_vertical(self.right_arm, self.left_arm,  self.LONG_ARM_REFS, point_separation=2)

            #Check results
            if len(result_left) == 0 or len(result_right) == 0: return False

            #33. Merge the independent left and right arm vertical optimizations
            self.optimize_hinge_center()

            #34. Pass hinge coordinates to arm coordinates
            self.copy_hinge_to_arm()

            #41. Perform fine horizontal optimization
            result_left  = self.arm_optimization_horizontal_fine(self.left_arm,  self.SMALL_ARM_REFS)
            result_right = self.arm_optimization_horizontal_fine(self.right_arm, self.SMALL_ARM_REFS)

            #Check results
            if len(result_left) == 0 or len(result_right) == 0: return False

            #42. Third round of vertical optimization
            result_left  = self.arm_optimization_vertical(self.left_arm,  self.right_arm, self.LONG_ARM_REFS, point_separation=1)
            result_right = self.arm_optimization_vertical(self.right_arm, self.left_arm,  self.LONG_ARM_REFS, point_separation=1)

            #Check results
            if len(result_left) == 0 or len(result_right) == 0: return False

            #43. Merge the independent left and right arm vertical optimizations
            self.optimize_hinge_center()

            #51. Determine and assign the ideal hinge parameters
            ideal_hinge_result = self.determine_ideal_hinge_parameters()
            if not ideal_hinge_result: return False

            #52. Determine inner side parameters
            inner_side_result = self.determine_inner_side_parameters()
            if not inner_side_result: return False

            #53. Check hinge quality
            self.check_hinge_quality()

            #54. Determine inner side parameters
            mask_side_result = self.determine_mask_side_parameters()
            if not mask_side_result: return False

            #55. Determine outer edge parameters
            self.determine_outer_edge_parameters()

            #56. Determine barcode edge parameters
            self.determine_barcode_edge_parameters()

            #57. Determine barcode parameters
            self.determine_barcode_parameters()

            #58. Calculate barcode scores
            barcode_score_result = self.calculate_barcode_scores()
            if not barcode_score_result: return False

            #59. Assign barcode state
            self.assign_barcode_state()

            return True


class EMParticleFinder:
    def __init__(self):
        self.particle_distance_cutoff  = None
        self.particle_score_cutoff     = None

class EMHingeFinder(EMParticleFinder):
    def __init__(self):
        self.hinge_angle_cutoff           = None
        self.inter_arm_distance_tolerance = None
        self.arm_to_inter_point_tolerance = None

class EMModel:
    def __init__(self):
        self.boxsize     = None
        self.background  = 0
        self.orientation = 0
        self.mask        = None
        self.img         = None

class EMLineModel(EMModel):
    def __init__(self):
        self.boxsize     = None
        self.background  = 0
        self.orientation = 0 
        self.linewidth   = 0
        self.gauss_sigma = None  

class EMOrigamiModel(EMModel):
    def __init__(self):
        self.modelname   = None

class EMOrigamiArmModel(EMOrigamiModel):
    def __init__(self):
        self.modelname   = None
        self.numlines    = None
        self.linesep     = None

class EMOrigamiHingeModel(EMOrigamiModel):
    def __init__(self):
        self.left_arm_length  = None
        self.left_arm_width   = None

        self.right_arm_length = None
        self.right_arm_width  = None

class EMOrigamiBarcode:
    def __init__(self, barcode_center=None, barcode_num_lines=None, barcode_length=None, barcode_angle=None, barcode_score=None, barcode_relative_coord=None, barcode_actual_coord=None, barcode_boxsize = None):
        self.barcode_relative_coord = barcode_relative_coord
        self.barcode_actual_coord   = barcode_actual_coord 
        self.barcode_center         = barcode_center
        self.barcode_num_lines      = barcode_num_lines
        self.barcode_length         = barcode_length
        self.barcode_angle          = barcode_angle
        self.barcode_score          = barcode_score
        self.barcode_boxsize        = barcode_boxsize
    
    def get_barcode_hash(self):
        '''
        Set barcode hash
        '''
        self.barcode_hash = "n%d-w%d-l%d-b%d-s%d"%(self.barcode_num_lines,10,self.barcode_length,self.barcode_boxsize,17)

    def set_actual_coord(self,start_coord, vertical_vec):
        '''
        Set barcode actual coordinate
        '''
        self.vertical_vec           = vertical_vec
        self.barcode_actual_coord   = start_coord + self.barcode_relative_coord*vertical_vec

    def barcode_distance(self,other):
        '''
        Measure the distance between this and the other barcode
        '''
        return EMtools.eular_distance(self.barcode_actual_coord,other.barcode_actual_coord)

    def determine_barcode_score(self, micrograph_img, barcode_refs, average_value):
        '''
        Determine the barcode score
        '''
        #Get the center coordinate
        barcode_coord = self.barcode_actual_coord
        
        #Crop the barcode img
        boxim_img     = EMtools.get_region(micrograph_img, barcode_coord[0]-0.5*self.barcode_boxsize, barcode_coord[1]-0.5*self.barcode_boxsize, self.barcode_boxsize, self.barcode_boxsize)
        
        #Check if the cropped image has the right dimensions
        if not boxim_img.shape[0] == self.barcode_boxsize or not boxim_img.shape[1] == self.barcode_boxsize: 
            return None

        #Get the reference img
        ref_img       = barcode_refs[self.barcode_angle]

        #Assign the comparison result
        self.barcode_score = EMtools.dot_cmp(boxim_img,ref_img, average_value)
        
        #Return barcode score
        return self.barcode_score