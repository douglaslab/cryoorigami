#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2017-10-19 15:02:13
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import numpy as np

#Constants
CLASS_COLORS = {0:'r',1:'b'}


def random_mask(micrograph_flat,boxsize):
    '''
    Make a random mask at boxsize*boxsize from a flat micrograph image
    '''
    return np.random.choice(micrograph_flat, size=(boxsize,boxsize))

def vector_length(vector):
    '''
    Calculate vector length
    '''
    return np.sqrt(np.sum(vector**2))

def create_polygon_mask(polygon_corners,origin=np.array([0,0]),boxsize=400):
    '''
    Create concave polygon mask 
    '''
    
    #Number of points
    num_corners = len(polygon_corners)

    #Put the corners into correct reference frame
    polygon_corners -= origin

    #Determine center of mass
    cmass = np.mean(polygon_corners,axis=0)

    #Create numpy array mask
    mask = np.ones((boxsize,boxsize))

    #Create meshgrid 
    X,Y = np.meshgrid(np.arange(boxsize),np.arange(boxsize))
    XY  = np.vstack(np.dstack((X,Y)))

    #Go over all the corners
    for i in range(num_corners):
        corner_1 = polygon_corners[i,:]
        corner_2 = polygon_corners[(i+1)%num_corners,:]

        #Create edge vector
        edge_vec    = normalize_vector(corner_2 - corner_1)

        #Cmass vector
        cmass_vec   = cmass - corner_1

        #Length projection of cmass vector over edge vector
        proj_length = np.sum(cmass_vec*edge_vec)

        #Determine projection point
        proj_point  = corner_1 + proj_length*edge_vec

        #Determine edge normal
        edge_normal = normalize_vector(cmass - proj_point)

        #Determine the cosa for all the elements fo XY
        XY_diff      = XY - corner_1
        XY_dot_edge  = np.sum(XY_diff*edge_normal,axis=1)

        #Return the indexes
        inside_1D = np.nonzero(XY_dot_edge >= 0)[0]
        inside_2D = np.unravel_index(inside_1D,(boxsize,boxsize))

        #Create new mask
        new_mask            = np.zeros((boxsize,boxsize))
        new_mask[inside_2D] = 1

        #Multiply the mask
        mask *= new_mask
            
    return mask


def eman_flip_y(eman_img):
    '''
    Flip eman image along y-axis
    '''
    eman_img.process_inplace("xform.flip",{"axis":"y"})

def eman_flip_x(eman_img):
    '''
    Flip eman image along x-axis
    '''
    eman_img.process_inplace("xform.flip",{"axis":"x"})

def get_region(micrograph_numpy, left_x, left_y, length_x, length_y):
    '''
    Return a region from a numpy 2D array
    '''

    round_left_x  = int(round(left_x))
    round_left_y  = int(round(left_y))

    round_right_x = int(round_left_x+length_x)
    round_right_y = int(round_left_y+length_y)

    #Store the difference values
    left_x_diff  = 0
    left_y_diff  = 0

    right_x_diff = length_x
    right_y_diff = length_y

    #Adjust the boundaries
    if round_left_x < 0: left_x_diff = -round_left_x; round_left_x = 0  
    if round_left_y < 0: left_y_diff = -round_left_y; round_left_y = 0

    if round_right_x >= micrograph_numpy.shape[1]: right_x_diff -= (round_right_x-micrograph_numpy.shape[1]);round_right_x = micrograph_numpy.shape[1]
    if round_right_y >= micrograph_numpy.shape[0]: right_y_diff -= (round_right_y-micrograph_numpy.shape[0]);round_right_y = micrograph_numpy.shape[0]

    #Create zero array
    background_img= np.zeros((length_x,length_y))

    #Assign the values from micrograph
    background_img[left_y_diff:right_y_diff,left_x_diff:right_x_diff] = micrograph_numpy[round_left_y:round_right_y, round_left_x:round_right_x]

    return background_img

def dot_cmp(boxim_img,ref_img,average_value=None):
    '''
    Comparison function
    '''
    dot_product = 1.0*np.sum(boxim_img*ref_img)/np.sum(ref_img)

    if average_value == None:
        cmp_score = dot_product/np.mean(boxim_img)
    else:
        cmp_score = dot_product/average_value

    return cmp_score

def ccc_cmp(boxim_img,ref_img,mask_img):
    '''
    Comparison function
    '''

    #Calculate the average intensity within the reference box
    average_value  = 1.0*np.sum(boxim_img*mask_img)/np.sum(mask_img)
    
    #The dot product
    cmp_score      = 1.0*np.sum((boxim_img-average_value)*ref_img)/np.sum(ref_img)

    return cmp_score

def determine_hinge_parameters(left_arm_coord, left_arm_angle, right_arm_coord, right_arm_angle, halfarm_length, dist_to_inner_edge):
    '''
    Find hinge parameters
    '''
    #Bad result vector
    bad_result = [None]*7

    inter_point    = intersection_point(left_arm_coord, left_arm_angle, right_arm_coord ,right_arm_angle)
    
    #Check if the intersection point is valid
    if inter_point[0] == None:
        return bad_result

    left_direction = left_arm_coord - inter_point 
    left_length    = np.sqrt(left_direction.dot(left_direction))
    
    #Check if the length is 0
    if left_length == 0:
        return bad_result
    
    left_direction = left_direction/left_length

    right_direction = right_arm_coord - inter_point
    right_length    = np.sqrt(right_direction.dot(right_direction))
    
    #Check if the length is 0
    if right_length == 0:
        return bad_result

    right_direction = right_direction/right_length

    #Determine the hinge angle
    hinge_angle  = ray_angle(left_direction,right_direction)

    #Check hinge angle result
    if hinge_angle == None or hinge_angle == 0:
        return bad_result

    x = dist_to_inner_edge/np.tan(0.5*hinge_angle/180*np.pi)+halfarm_length

    #Find the new arm centers
    new_left_arm_coord  = inter_point + x*left_direction
    new_right_arm_coord = inter_point + x*right_direction

    #Determine new origin
    new_origin  = 0.5*(new_left_arm_coord+new_right_arm_coord)

    #Ideal inter arm distance
    ideal_inter_dist = 2.0*x*np.sin(0.5*hinge_angle*np.pi/180)

    #Get the intersection angle
    dirx,diry   = inter_point - new_origin
    hipo        = np.sqrt(dirx**2+diry**2)

    #Check the hipotenus
    if not hipo > 0:
        return bad_result

    inter_angle = np.arccos(1.0*dirx/hipo)*180.0/np.pi
    
    if diry < 0:
        inter_angle = 360 - inter_angle
    
    return new_left_arm_coord, new_right_arm_coord, new_origin, hinge_angle, inter_angle, inter_point, ideal_inter_dist

def eular_distance(particle_coord, neighbors_coord):
    '''
    Calculate the distance between a particle and potential particle neigbors
    ''' 

    dist = neighbors_coord - particle_coord
    if len(neighbors_coord.shape) == 1:
        return np.sqrt(np.sum(dist**2))

    return np.sqrt(np.sum(dist**2,axis=1))

def unit_vector(angle):
    '''
    Make a unit vector oriented at an angle
    '''
    x = np.cos(np.pi*angle/180.0)
    y = np.sin(np.pi*angle/180.0)

    return np.array([x,y]) 

def intersection_point(particle1_coord,particle1_angle,particle2_coord,particle2_angle):
    '''
    Returns the intersection point for two lines 
    passing through two points at an angle
    '''

    rad1 = particle1_angle*np.pi/180.0
    rad2 = particle2_angle*np.pi/180.0
    
    #If the slopes of the lines are the same return None
    if rad1 == rad2:
        return np.array([None,None])

    sin1 = np.sin(rad1)
    cos1 = np.cos(rad1)

    sin2 = np.sin(rad2)
    cos2 = np.cos(rad2)

    x1,y1 = particle1_coord
    x2,y2 = particle2_coord

    c1 = -sin1*x1+cos1*y1
    c2 = -sin2*x2+cos2*y2

    M   = np.array([[-sin1,cos1],[-sin2,cos2]])
    c   = np.array([c1,c2])

    #Find the inteesection point
    try:
        sol = np.linalg.solve(M, c)
    except np.linalg.linalg.LinAlgError as err:
        if 'Singular matrix' in err.message:
            return np.array([None,None])
        else:
            raise

    return sol

def ray_angle(vector_1,vector_2):
    '''
    Returns the angle between two vectors
    '''
    vector_length_product = np.linalg.norm(vector_1)*np.linalg.norm(vector_2)
    
    #If either of the vector has 0 length, return None
    if vector_length_product == 0:
        return None
    
    #Cos value
    cos_value = np.dot(vector_1,vector_2)/vector_length_product

    #Check cos value
    if cos_value > 1 or cos_value < -1:
        return None
    
    return np.arccos(cos_value)*180.0/np.pi

def normalize_vector(vector):
    '''
    Normalize vector by its length
    '''
    return vector/np.sqrt(vector.dot(vector)) 

def make_vector_points(start_point, vector, length):
    '''
    Make points along a vector
    '''
    int_points    = np.arange(round(length))
    vector_points = start_point+np.array([x*vector for x in int_points])

    return vector_points

def img_out_of_bounds(EM_img,coord,particle_radius=0):
    '''
    Determine if coord is image out of bounds
    '''
    nx,ny = EM_img["nx"], EM_img["ny"]

    #Check if the coordinate is within image bounds
    if coord[0]-particle_radius < 0 or coord[0]+particle_radius >=nx or coord[1]-particle_radius < 0 or coord[1]+particle_radius >= ny:
        return True
    else:
        return False 

def make_line(center_coord, line_angle, line_length):
    '''
    Determine the start and end points
    for a line center around center coordinate
    '''
    rad1 = line_angle*np.pi/180.0
    
    sin1 = np.sin(rad1)
    cos1 = np.cos(rad1)

    x0,y0 = center_coord
    half_length = 0.5*round(line_length)

    dy = half_length*sin1
    dx = half_length*cos1

    #Point 1
    x1 = x0 - dx
    y1 = y0 - dy

    #Point 2
    x2 = x0 + dx
    y2 = y0 + dy

    return np.array([[x1,y1],[x2,y2]])

def make_line_points(center_coord, angle, length):
    '''
    Make points along a line
    '''

    start_point,end_point = make_line(center_coord,angle,length)

    x_points  = np.around(np.linspace(start_point[0],end_point[0],round(length)),dtype=int)
    y_points  = np.around(np.linspace(start_point[1],end_point[1],round(length)),dtype=int)

    return x_points,y_points

def rotation_matrix(angle):
    '''
    Make a rotation matrix
    for rotation angle
    '''

    theta = np.radians(angle)
    c, s  = np.cos(theta), np.sin(theta)
    R     = np.matrix([[c, -s], [s, c]])

    return R