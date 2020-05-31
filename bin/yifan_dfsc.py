#!/usr/bin/env python

# Adapted from Jun Rui from Yifan Cheng Lab
# 03/28/19

import sys,math,random,numpy,struct,argparse
import matplotlib
import os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob

import time
import timeit
import cryoorigami.utilities as util

# mrc image class ------------------------------------------------------------
class mrc_image:
    def __init__(self,filename):
        # Header items -------------------------------------------------------
        # 0-11    = nx,ny,nz (i)
        # 12-15   = mode (i)
        # 16-27   = nxstart,nystart,nzstart (i)
        # 28-39   = mx,my,mz (i)
        # 40-51   = cell size (x,y,z) (f)
        # 52-63   = cell angles (alpha,beta,gamma) (f)
        # 64-75   = mapc,mapr,maps (1,2,3) (i)
        # 76-87   = amin,amax,amean (f)
        # 88-95   = ispg,nsymbt (i)
        # 96-207  = 0 (i)
        # 208-215 = cmap,stamp (c)
        # 216-219 = rms (f)
        # 220-223 = nlabels (i)
        # 224-1023 = 10 lines of 80 char titles (c)
        # variables ----------------------------------------------------------
        self.filename=filename
        self.byte_pattern1='=' + '10i' + '6f' + '3i' + '3f' + '30i' + '4s'*2 + 'fi'
        #self.byte_pattern1='=' + 'i'*10 + 'f'*6 + 'i'*3 + 'f'*3 + 'i'*30 + '4s'*2 + 'fi'
        self.byte_pattern2='=' + '800s'
        print("byte_pattern1", self.byte_pattern1)

    # read header only -------------------------------------------------------
    def read_head(self):
        input_image=open(self.filename,'rb')
        self.header1=input_image.read(224)
        self.header2=input_image.read(800)
        self.dim=struct.unpack(self.byte_pattern1,self.header1)[:3]   #(dimx,dimy,dimz)
        input_image.close()

    # read entire image/stack ------------------------------------------------
    def read(self):
        input_image=open(self.filename,'rb')
        self.header1=input_image.read(224)
        self.header2=input_image.read(800)
        header = struct.unpack(self.byte_pattern1,self.header1)
        self.dim=header[:3]   #(dimx,dimy,dimz)
        self.imagetype=header[3]
        print("header", header)
        self.cellsize=header[10:13] # cell size in angstroms 
        #0: 8-bit signed, 1:16-bit signed, 2: 32-bit float, 6: unsigned 16-bit (non-std)
        if (self.imagetype == 0):imtype=numpy.uint8
        elif (self.imagetype ==1):imtype='h'
        elif (self.imagetype ==2):imtype='f4'
        elif (self.imagetype ==6):imtype='H'
        else:imtype='unknown'   #should put a fail here
        input_image_dimension=(self.dim[2],self.dim[1],self.dim[0]) 
        self.image_data=numpy.fromfile(file=input_image,dtype=imtype,count=self.dim[0]*\
                                           self.dim[1]*self.dim[2]).reshape(input_image_dimension)
        self.image_data=self.image_data.astype(numpy.float32)
        input_image.close()

    # read image slice -------------------------------------------------------
    def read_slice(self,nslice):
        input_image=open(self.filename,'rb')
        self.header1=input_image.read(224)
        self.header2=input_image.read(800)
        header=struct.unpack(self.byte_pattern1,self.header1)
        self.dim=header[:3]   #(dimx,dimy,dimz)
        self.imagetype=header[3]  
        #0: 8-bit signed, 1:16-bit signed, 2: 32-bit float, 6: unsigned 16-bit (non-std)
        if (self.imagetype == 0):
            imtype=numpy.uint8
            nbytes=1
        elif (self.imagetype ==1):imtype='h'
        elif (self.imagetype ==2):
            imtype='f4'
            nbytes=4
        elif (self.imagetype ==6):imtype='H'
        else:imtype='unknown'   #should put a fail here
        input_image_dimension=(self.dim[1],self.dim[0]) 
        input_image.seek(nbytes*nslice*self.dim[0]*self.dim[1],1) #find slice
        self.image_data=numpy.fromfile(input_image,dtype=imtype,count=self.dim[0]*\
                                           self.dim[1]).reshape(input_image_dimension)
        self.image_data=self.image_data.astype(numpy.float32)
        input_image.close()

    # write mrc --------------------------------------------------------------
    def write(self,image,apix,output=numpy.ones(1)):
        if len(image.shape)==2:image=image.reshape((1,image.shape[0],image.shape[1])) #3D image
        output_image=open(self.filename,'w')
        # output all images --------------------------------------------------
        if output.shape[0]==1:
            dim=[image.shape[2],image.shape[1],image.shape[0]]
            amin=image.min()
            amax=image.max()
            amean=numpy.average(image)
            dtype=image.dtype
            if dtype==numpy.uint8:dtype=0
            elif dtype==numpy.float32:dtype=2
            else:print('WARNING: Unknown data type')
            header1=struct.pack(self.byte_pattern1,dim[0],dim[1],dim[2],dtype,0,0,0,dim[0],dim[1],dim[2],
                                dim[0]*apix,dim[1]*apix,dim[2]*apix,90,90,90,1,2,3,amin,amax,amean,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'MAP','DA\x00\x00',0.0,1)
            header2=struct.pack(self.byte_pattern2,' ') #comments
            print("header1", dim[0]*apix)
            output_image.write(header1)
            output_image.write(header2)
            image.tofile(output_image)
            output_image.close()
        # output only select images ------------------------------------------
        else:
            nparts=output.shape[0]
            loc=numpy.where(output==1)[0]
            nout=loc.shape[0]
            if ntrue==0:print('ERROR: no labeled images')
            else:
                dum=image.shape # Get dimensions of image/stk
                dim=[image.shape[2],image.shape[1],nout]
                image2=numpy.zeros((dim[2],dim[1],dim[0]),dtype='f4')
                for i in range(nout):image2[i,:,:]=image[int(loc[i]),:,:] # Take true particles only
                amin=image2.min()
                amax=image2.max()
                amean=numpy.average(image2)
                header1=struct.pack(self.byte_pattern1,dim[0],dim[1],dim[2],2,0,0,0,dim[0],dim[1],dim[2],
                                dim[0]*apix,dim[1]*apix,dim[2]*apix,90,90,90,1,2,3,amin,amax,amean,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'MAP','DA\x00\x00',0.0,1)
                header2=struct.pack(self.byte_pattern2,' ')
                output_image.write(header1)
                output_image.write(header2)
                image2.tofile(output_image)
                output_image.close()

# arguments from user -------------------------------------------------------
class argument_parser():
    def __init__(self):
        parser=argparse.ArgumentParser(description='Calculate directional Fourier shell correlation (dFSC). dFSC curves are output to dfsc_table.txt. Each row corresponds to the dFSC in the direction specified by the vector in Fibonacci_points.txt. Plot using GNUplot command: plot "dfsc_1d.txt" matrix w l')
        parser.add_argument('half1',help='Half-map 1')
        parser.add_argument('half2',help='Half-map 2')
        parser.add_argument('whole',help='Whole-map')
        parser.add_argument('--o',      default=None, help='Output directory')
        parser.add_argument('--mask',   default='',  help='Mask for FSC calculation')
        parser.add_argument('--apix',   default='',  help='Pixel size (Angstroms)')
        parser.add_argument('--bfactor',default='',  help='B-factor for global sharpening (Negative)')
        parser.add_argument('--ncones', default=500, help='Number of cones to use for dFSC calculation')
        parser.add_argument('--angle',  default=20., help='Apex angle of cones for dFSC calculation')
        self.args=parser.parse_args()

# main code -----------------------------------------------------------------
class Main():
    def __init__(self):
        # get user input ----------------------------------------------------
        parser=argument_parser()
        args=parser.args
        file_half1=args.half1; print('Half-map 1:',file_half1)
        file_half2=args.half2; print('Half-map 2:',file_half2)
        file_whole=args.whole; print('Whole-map:',file_whole)
        file_mask=args.mask; print('Mask:',file_mask)
        apix=args.apix; print('Pixel size:',apix)
        adhoc_bfactor=args.bfactor; print('B-factor:',adhoc_bfactor)
        if apix=='':
            apix=0
        else:
            apix=float(apix)
        ncones=int(args.ncones); print('Number of cones:',ncones)
        halfangle=float(args.angle); print('Angular spread:',halfangle)
        # Set output directory
        self.set_output_directory(args.o)

        # Write the input parameters
        # Prepare args dict
        args_dict = {'half1':     args.half1,
                     'half2':     args.half2,
                     'whole':     args.whole,
                     'mask':      args.mask,
                     'apix':      args.apix,
                     'bfactor':   args.bfactor,
                     'ncones':    args.ncones,
                     'angle':     args.angle,
                     'output':    args.o}

        args_filename = self.output_directory+'/args.yaml'
        util.write_config_file(args_filename, args_dict)

        # load half maps ----------------------------------------------------
        print('Reading maps')

        dummy1=mrc_image(file_half1)
        dummy1.read()
        dummy2=mrc_image(file_half2)
        dummy2.read()
        dummy4=mrc_image(file_whole)
        dummy4.read()
        if file_mask!='':
            dummy3=mrc_image(file_mask)
            dummy3.read()
            maskmap=dummy3.image_data
        else:
            imdim=dummy1.image_data.shape
            maskmap=numpy.ones((imdim[0],imdim[1],imdim[2]),dtype=numpy.float32)
        if apix == 0:
            apix = dummy4.cellsize[0] / dummy4.image_data.shape[0] 
            print("cellsize", apix)
        map_half1=dummy1.image_data*maskmap
        map_half2=dummy2.image_data*maskmap
        map_whole=dummy4.image_data*maskmap

        # Fourier transforms ------------------------------------------------
        print('Computing FFTs')

        print("map_half1.shape", map_half1.shape)
        half1_ft=numpy.fft.rfftn(map_half1)
        print("half1_ft.shape", half1_ft.shape)
        half1_ft=self.fft_wrap(half1_ft)
        print("half1_ft", half1_ft.shape)
        half2_ft=numpy.fft.rfftn(map_half2)
        half2_ft=self.fft_wrap(half2_ft)
        whole_ft=numpy.fft.rfftn(map_whole)
        print("whole_ft", whole_ft.shape)
        whole_ft=self.fft_wrap(whole_ft)

        # print "map_half1_shape", map_half1.shape
        # print "half1_ft_shape", half1_ft.shape


        # compute points on a sphere ----------------------------------------
        print('Computing Fibonacci sphere')

        arsize=map_half1.shape[0]
        arsizeby2=arsize/2
        points=((self.fibonacci_sphere(ncones)+1)*arsizeby2).astype(int)

        # print "points_shape", points.shape
        # print "points", points

        fccmap=numpy.zeros((arsize,arsize,arsize),dtype=numpy.float64)
        for i in range(points.shape[0]):
            fccmap[points[i,0],points[i,1],points[i,2]]=1.
        # create mask -------------------------------------------------------
        print('Computing mask')

        xyzgrid=numpy.mgrid[-arsizeby2:arsizeby2,-arsizeby2:arsizeby2,-arsizeby2:arsizeby2]
        dist2=xyzgrid[0]**2+xyzgrid[1]**2+xyzgrid[2]**2
        mask=numpy.zeros((arsize,arsize,arsize),dtype=numpy.uint8)
        coord=numpy.where(dist2<arsizeby2**2)
        mask[coord[0],coord[1],coord[2]]=1


        # computer upper half fsc
        # mask_half = mask
        # mask_half[:,:,arsizeby2:]=0
        # coord_upper = numpy.where( mask_half > 0)
        # x = coord_upper[0]
        # y = coord_upper[1]
        # z = coord_upper[2]
        # xyz = numpy.array([x,y,z])
        # fcc_upper_half = self.fcc_compute(half1_ft, half2_ft, xyz, dist2)
        # mask2 = numpy.zeros((arsize, arsize, arsize), dtype=numpy.float32)
        # cref = self.fcc_Cref_mask(mask2, xyz, fcc_upper_half)
        # weighted_ft = whole_ft * cref
        # fit_minres = 10
        # fit_maxres = 0
        # figure = "upper_half_GuinierPlot.png"
        # bfactor_global = self.bfactor_compute_save(weighted_ft, xyz, fcc_upper_half, apix, dist2, fit_minres, fit_maxres, figure)



        # compute global fsc and apply global bfactor sharpening
        fcc_global = numpy.zeros(arsizeby2,dtype=numpy.float64)
        x = coord[0]
        y = coord[1]
        z = coord[2]
        xyz = numpy.array([x,y,z])
        fcc_global = self.fcc_compute(half1_ft, half2_ft, xyz, dist2)
        print("fcc_global.shape[0]", fcc_global.size)
        mask2 = numpy.zeros((arsize, arsize, arsize), dtype=numpy.float32)
        cref = self.fcc_Cref_mask(mask2, xyz, fcc_global)
        outmrc=mrc_image(self.output_directory+'/global_Cref.mrc')
        outmrc.write(cref.astype(numpy.float32),apix)
        weighted_ft = whole_ft * cref
        weighted = numpy.fft.irfftn(self.fft_unwrap(weighted_ft))
        outmrc=mrc_image(self.output_directory+'/gWeighted_uSharpened.mrc')
        outmrc.write(weighted.astype(numpy.float32),apix)
        fit_minres = 10
        fit_maxres = 0
        figure = self.output_directory+"/global_GuinierPlot.png"
        if adhoc_bfactor == '':
            bfactor_global = self.bfactor_compute_save(weighted_ft, xyz, fcc_global, apix, dist2, fit_minres, fit_maxres, figure)
        else:
            bfactor_global = float(adhoc_bfactor)
        #bfactor_global = -200
        weighted_ft *= numpy.exp(-(bfactor_global / 4.) * dist2*mask/(arsize * apix)**2)
        whole_sharpened = numpy.fft.irfftn(self.fft_unwrap(weighted_ft))
        filename = self.output_directory+'/gWeighted_gSharpened' + str(bfactor_global) + '.mrc'
        outmrc=mrc_image(filename)
        outmrc.write(whole_sharpened.astype(numpy.float32),apix) 





        # compute cone points -----------------------------------------------
        points2=points.astype(numpy.float64)/numpy.float64(arsizeby2)-1.

        # print "points2", points2

        dispoints=numpy.sqrt(points2[:,0]**2+points2[:,1]**2+points2[:,2]**2)

        # print "dispoints", dispoints.shape, dispoints

        points2[:,0]=points2[:,0]/dispoints
        points2[:,1]=points2[:,1]/dispoints
        points2[:,2]=points2[:,2]/dispoints

        # print "points2", points2

        fcc = numpy.zeros((points.shape[0],arsizeby2),dtype=numpy.float64)
        bfactor = numpy.zeros(points.shape[0], dtype=numpy.float64)
        #angthres=numpy.pi*numpy.sqrt(180*360/points.shape[0])/180./2.
        angthres=numpy.pi*halfangle/180.

        # print "angthres", angthres

        cosangthres=numpy.cos(angthres)
        coord=numpy.where(dist2==0)
        dist3=dist2[:,:,:]
        dist3[coord[0],coord[1],coord[2]]=1
        vecarray=xyzgrid/numpy.sqrt(dist3)

        # print "vecarray", vecarray

        # calculate first cone ----------------------------------------------
        print('Computing first cone')

        # print time.strftime("%Y-%m-%d %H:%M:%S")

        cosangarray0=numpy.dot(vecarray[0],points2[0,0])+\
                         numpy.dot(vecarray[1],points2[0,1])+\
                         numpy.dot(vecarray[2],points2[0,2])
        cosangarray0=cosangarray0*mask
        conepoints0=(numpy.where((cosangarray0>cosangthres)))

        # print "conepoints0", numpy.ndim(conepoints0), conepoints0

        x0=conepoints0[0]
        y0=conepoints0[1]
        z0=conepoints0[2]
        xyz0=numpy.array([x0,y0,z0])

        # print "xyz0", xyz0.shape, xyz0

        fcc[0,:]=self.fcc_compute(half1_ft,half2_ft,xyz0,dist2)
        mask2=numpy.zeros((arsize,arsize,arsize),dtype=numpy.uint8)
        mask2[x0,y0,z0]=1                                                  # mask2 is a mask for a cone
        dist4=dist2*mask2
        mask2=numpy.zeros((arsize,arsize,arsize),dtype=numpy.float32)      # mask2 also works as an initialization
        mask3=self.fcc_mask(mask2,xyz0,fcc[0])                             # mask3 is fcc for a cone
        mask4=numpy.zeros((arsize,arsize,arsize),dtype=numpy.float32)
        mask4+=mask3                                                       # mask4 is fcc for the whole sphere 
        mask5=numpy.zeros((arsize,arsize,arsize),dtype=numpy.float32)      
        mask5[x0,y0,z0]=1                                                  # mask5 records overlaps between cones

        cref_cone = self.fcc_Cref_mask(mask2, xyz0, fcc[0])
        mask_Cref = numpy.zeros((arsize,arsize,arsize),dtype=numpy.float32)
        mask_Cref += cref_cone                                             
        cone_ft = whole_ft * cref_cone
        bfactor[0] = self.bfactor_compute(cone_ft, xyz0, fcc[0], apix, dist2, fit_minres, fit_maxres)
        bfactor_cone = numpy.zeros((arsize,arsize,arsize),dtype=numpy.uint8)
        bfactor_cone[x0, y0, z0] = bfactor[0]
        mask_bfactor = numpy.zeros((arsize,arsize,arsize),dtype=numpy.float32)
        mask_bfactor += bfactor_cone

        
        # print time.strftime("%Y-%m-%d %H:%M:%S")
        totaltime = 0
        # iterate over all points -------------------------------------------
        for i in range(1,points2.shape[0]):
            #sys.stdout.flush()
            #sys.stdout.write("\rComputing dFSC...")
            #sys.stdout.flush()
            print("cone", i)

            start_i = timeit.default_timer()

            cosangarray1=numpy.dot(vecarray[0],points2[i,0])+\
                             numpy.dot(vecarray[1],points2[i,1])+\
                             numpy.dot(vecarray[2],points2[i,2])
            cosangarray1=cosangarray1*mask
            conepoints1=(numpy.where((cosangarray1>cosangthres)))
            print("fibnaci points i", points2[i])
            # print "conepoints0", numpy.ndim(conepoints0), conepoints0

            x1=conepoints1[0]
            y1=conepoints1[1]
            z1=conepoints1[2]
            xyz1=numpy.array([x1,y1,z1])
            print("xyz1.size", xyz1.size)
 
            # rotaxis=points2[i]+points2[0]
            # rotaxismag=numpy.linalg.norm(rotaxis)

            # # print "rotaxis", rotaxis.shape, rotaxis
            # # print "rotaxismag", rotaxismag.shape, rotaxismag


            # if rotaxismag==0.:
            #     xyz1=numpy.array([-x0,-y0,-z0])
            # else:
            #     rotaxis=rotaxis/rotaxismag
            #     rotmatrix=self.rotation_matrix_180(rotaxis)
            #     xyz1=numpy.round(numpy.dot(rotmatrix,xyz0-arsizeby2)).astype(int)+arsizeby2
            # 
            #     # print "xyz1", xyz1.shape, xyz1

            #     xyz1[xyz1>arsize-1] = arsize-1

                # if numpy.max(xyz1)>arsize-1:
                #     coord=numpy.where(xyz1>arsize-1)
                #     xyz1[coord[0],coord[1]]=arsize-1
                # elif numpy.min(xyz1)<-arsize+1:
                #     coord=numpy.where(xyz1<-arsize+1)
                #     xyz1[coord[0],coord[1]]=-arsize+1

            # print "coord_shape", numpy.ndim(coord), coord[0].shape

            #  before = timeit.default_timer()
            # print "before", before - start_i 

            # fcc computation time
            # start_fcc = timeit.default_timer()

            fcc[i,:]=self.fcc_compute(half1_ft,half2_ft,xyz1,dist2)

            stop_fcc = timeit.default_timer()
            # print "fcc", i, stop_fcc - start_fcc
            totaltime += stop_fcc - start_i
            
            mask2[:,:,:]=0
            mask2[xyz1[0],xyz1[1],xyz1[2]]=1
           
            # print "mask2", mask2.shape, mask2 
            
            dist4=dist2*mask2
        
            start_mask2 = timeit.default_timer()
            
            mask2[:,:,:]=0
            mask3=self.fcc_mask(mask2,xyz1,fcc[i])
            
            stop_mask2 = timeit.default_timer()

            # print "mask2", stop_mask2 - start_mask2        
            
            #mask4[xyz1[0],xyz1[1],xyz1[2]]=0
            mask4+=mask3
        
            mask2[:,:,:]=0
            mask2[xyz1[0],xyz1[1],xyz1[2]]=1
        
            mask5+=mask2
            
            mask2[:,:,:]=0
            cref_cone = self.fcc_Cref_mask(mask2, xyz1, fcc[i])
            mask_Cref += cref_cone
            cone_ft = whole_ft * cref_cone
            bfactor[i] = self.bfactor_compute(cone_ft, xyz1, fcc[i], apix, dist2, fit_minres, fit_maxres)
            bfactor_cone = numpy.zeros((arsize,arsize,arsize),dtype=numpy.uint8)
            bfactor_cone[x1, y1, z1] = bfactor[i]
            mask_bfactor += bfactor_cone
             
            stop_i = timeit.default_timer()
            # print "after" , stop_i - stop_fcc
            # print "iteration", i, stop_i - start_i

            #sys.stdout.flush()
            #sys.stdout.write("\rComputing dFSC...%s"%(int(100.*float(i)/float(points.shape[0]-1))))
            #sys.stdout.write('%')
            #sys.stdout.flush()
        xyz2=numpy.where(mask5>0)
        mask4[xyz2[0],xyz2[1],xyz2[2]]=mask4[xyz2[0],xyz2[1],xyz2[2]]\
                                        /mask5[xyz2[0],xyz2[1],xyz2[2]]
        mask_Cref[xyz2[0],xyz2[1],xyz2[2]] = mask_Cref[xyz2[0],xyz2[1],xyz2[2]]\
                                             /mask5[xyz2[0],xyz2[1],xyz2[2]]
        mask_bfactor[xyz2[0],xyz2[1],xyz2[2]] = mask_bfactor[xyz2[0],xyz2[1],xyz2[2]]\
                                                /mask5[xyz2[0],xyz2[1],xyz2[2]]
        gWeighted_dSharpened_ft = weighted_ft * numpy.exp(-(mask_bfactor / 4.) * dist2*mask/(arsize*apix)**2)
        gWeighted_dSharpened = numpy.fft.irfftn(self.fft_unwrap(gWeighted_dSharpened_ft))
        outmrc=mrc_image(self.output_directory+'/gWeighted_dSharpened.mrc')
        outmrc.write(gWeighted_dSharpened.astype(numpy.float32),apix)

        directional_weighted_ft = whole_ft * mask_Cref
        directional_weighted = numpy.fft.irfftn(self.fft_unwrap(directional_weighted_ft))
        outmrc=mrc_image(self.output_directory+'/dWeighted_uSharpened.mrc')
        outmrc.write(directional_weighted.astype(numpy.float32),apix)

        directional_sharpened_ft = whole_ft * numpy.exp(-(mask_bfactor / 4.) * dist2*mask/(arsize*apix)**2)
        directional_sharpened = numpy.fft.irfftn(self.fft_unwrap(directional_sharpened_ft))
        outmrc=mrc_image(self.output_directory+'/uWeighted_dSharpened.mrc')
        outmrc.write(directional_sharpened.astype(numpy.float32),apix) 

        directional_weighted_sharpened_ft = whole_ft * mask_Cref * numpy.exp(-(mask_bfactor / 4.) * dist2*mask/(arsize*apix)**2)
        directional_weighted_sharpened = numpy.fft.irfftn(self.fft_unwrap(directional_weighted_sharpened_ft))
        outmrc=mrc_image(self.output_directory+'/dWeighted_dSharpened.mrc')
        outmrc.write(directional_weighted_sharpened.astype(numpy.float32),apix)

        outmrc=mrc_image(self.output_directory+'/mask_Cref.mrc')
        outmrc.write(mask_Cref.astype(numpy.float32),apix)

        mask_Cref_fsc = numpy.sqrt(2*mask4/(1+mask4))
        outmrc=mrc_image(self.output_directory+'/directional_Cref.mrc')
        outmrc.write(mask_Cref_fsc.astype(numpy.float32),apix)

        weight_sharpened_ft = whole_ft * mask_Cref * numpy.exp(-(bfactor_global / 4.) * dist2*mask/(arsize*apix)**2)
        weight_sharpened = numpy.fft.irfftn(self.fft_unwrap(weight_sharpened_ft))
        #weight_sharpened = numpy.fft.irfftn(weight_sharpened_ft)
        outmrc=mrc_image(self.output_directory+'/dWeighted_gSharpened.mrc')
        outmrc.write(weight_sharpened.astype(numpy.float32),apix)

        print('totaltime', totaltime)

        #print time.strftime("%Y-%m-%d %H:%M:%S")

        # plot 1d dFSCs -----------------------------------------------------
        x=numpy.arange(float(fcc.shape[1]))+1
        if apix==0:
            plt.plot(x,fcc[0,:],'purple',linewidth=1.0,alpha=0.05,label='dFSC')
        else:
            x2=x/(2.*fcc.shape[1]*apix)
            plt.plot(x2,fcc[0,:],'purple',linewidth=1.0,alpha=0.05,label='dFSC')
        for i in range(1,fcc.shape[0]):

            # print "fcc", i
            # print fcc[i,0:9]

            if apix==0:
                plt.plot(x,fcc[i,:],'purple',linewidth=1.0,alpha=0.05)
            else:
                plt.plot(x2,fcc[i,:],'purple',linewidth=1.0,alpha=0.05)
        fcc_avg=numpy.average(fcc,axis=0)
        if apix==0:
            plt.plot(x,fcc_avg,'green',linewidth=2.0,alpha=1.0,label='average dFSC')
        else:
            plt.plot(x2,fcc_avg,'green',linewidth=2.0,alpha=1.0,label='average dFSC')
        plt.legend(loc='best')
        if apix==0:
            plt.axis([0,fcc.shape[1],0,1])
            plt.xlabel('Spatial frequency (pixels)')
        else:
            plt.axis([0.,1./(2.*apix),0,1])
            #plt.gca().invert_xaxis()
            plt.xlabel('Spatial frequency (1/Angstroms)')
        plt.ylabel('dFSC')
        plt.savefig(self.output_directory+"/dfsc_1d_plot.svg",format="svg")
        # write out fcc table and points ------------------------------------

        # print "write", time.strftime("%Y-%m-%d %H:%M:%S")

        numpy.savetxt(self.output_directory+'/fibonacci_points.txt',points2,fmt="%.3f")
        numpy.savetxt(self.output_directory+'/dfsc_1d.txt',fcc,fmt="%.3f")
        with open(self.output_directory+'/dfsc_axes.txt','w') as f:
            if apix==0:
                f.write('# Pixels\n')
                numpy.savetxt(f,x,newline=" ",fmt="%.3f")
            else:
                f.write('# Pixels\n')
                f.write('# Spatial frequency (1/Angstroms)\n')
                f.write('# 1/Spatial frequency (Angstroms)\n')
                numpy.savetxt(f,(x,x2,1./x2),fmt="%.3f")
        outmrc=mrc_image(self.output_directory+'/fibonacci_sphere.mrc')
        outmrc.write(fccmap.real.astype(numpy.float32),apix)
        # write out dFSC map ------------------------------------------------
        outmrc=mrc_image(self.output_directory+'/dfsc_3d.mrc')
        outmrc.write(mask4.astype(numpy.float32),apix)

        # print "write", time.strftime("%Y-%m-%d %H:%M:%S")


    # Fibonacci sphere ------------------------------------------------------
    def fibonacci_sphere(self,samples=1,randomize=True):
        rnd=1.
        if randomize:
            rnd=random.random()*samples
        points=numpy.zeros((samples,3),dtype=numpy.float64)
        offset=2./samples
        increment=math.pi*(3.-math.sqrt(5.))
        for i in range(samples):
            y=((i*offset)-1)+(offset/2)
            r=math.sqrt(1-pow(y,2))
            phi=((i+rnd)%samples)*increment
            x=math.cos(phi)*r
            z=math.sin(phi)*r
            points[i,0]=x
            points[i,1]=y
            points[i,2]=z
        return points

    # rotation matrix -------------------------------------------------------
    def rotation_matrix_180(self,axis):
        x=axis[0]
        y=axis[1]
        z=axis[2]
        R=numpy.zeros((3,3),dtype=numpy.float64)
        R[0,0]=-1+2*x**2
        R[0,1]=2*x*y
        R[0,2]=2*x*z
        R[1,0]=2*x*y
        R[1,1]=-1+2*y**2
        R[1,2]=2*y*z
        R[2,0]=2*x*z
        R[2,1]=2*y*z
        R[2,2]=-1+2*z**2
        return R

    # compute Fourier cone correlation --------------------------------------
    def fcc_compute(self,half1,half2,points,dist):
        # before_bin = timeit.default_timer()
        x0=points[0]
        y0=points[1]
        z0=points[2]
        ftsize=half1.shape[0]/2
        fcc=numpy.zeros(ftsize,dtype=numpy.float64)
        fc = (half1[x0,y0,z0] * half2[x0,y0,z0].conjugate()).real
        # mid_bin = timeit.default_timer()
        dist10 = numpy.round(numpy.linalg.norm(points - dist.shape[0]/2, axis=0)).astype(int)

        # dist10 = numpy.round(numpy.sqrt(x2**2+y2**2+z2**2)).astype(int)
 
        # start_bin = timeit.default_timer()       
        sum1_bin = numpy.bincount(dist10, fc, minlength=0)
        sum2_bin = numpy.sqrt(numpy.bincount(dist10, (half1[x0,y0,z0]*half1[x0,y0,z0].conjugate()).real, minlength=0) *\
                              numpy.bincount(dist10, (half2[x0,y0,z0]*half2[x0,y0,z0].conjugate()).real, minlength=0))
        # stop_bin = timeit.default_timer()
        # print "before_bin", start_bin - before_bin
        # print "mid_bin", mid_bin - before_bin
        # print "mid2_bin", start_bin - mid_bin
        # print "bin_time", stop_bin - start_bin
        # print "sum1_bin", sum1_bin[:ftsize].size, sum1_bin[0:5]
        # print "sum2_bin", sum2_bin[:ftsize].size, sum2_bin[0:5]
        sum1_bin[0:3][sum1_bin[0:3] == 0.] = 1          # there may be few points in this small shell with radius less than 3
        sum1_bin[0:3][sum1_bin[0:3] <= 0.] = numpy.abs(sum1_bin[0:3][sum1_bin[0:3] <= 0.])
        sum2_bin[0:3][sum2_bin[0:3] == 0.] = 1
        fcc = sum1_bin[:ftsize] / sum2_bin[:ftsize]
        # print "sum1_bin", sum1_bin, "\n sum2_bin", sum2_bin

        # mask=numpy.zeros((dist.shape[0],dist.shape[1],dist.shape[2]),dtype=numpy.float64)
        # mask[x0,y0,z0]=1.
        # dist2=dist*mask

        # for i in range(ftsize):

        #     # index computation time
        #     start_index = timeit.default_timer()

        #     coords=numpy.where((dist2>(i)**2) & (dist2<(i+2)**2))

        #     stop_index = timeit.default_timer()
        #     print "index", i, stop_index - start_index

        #     x1=coords[0]
        #     y1=coords[1]
        #     z1=coords[2]
        #     if x1!=[]:
        #         fc=half1[x1,y1,z1]*half2[x1,y1,z1].conjugate()
        #         sum1=numpy.sum(fc).real
        #         print "half1*half2", x1.size, type(fc), fc.size
        #         sum2=numpy.sqrt(numpy.sum(half1[x1,y1,z1]*half1[x1,y1,z1].conjugate())*\
        #                         numpy.sum(half2[x1,y1,z1]*half2[x1,y1,z1].conjugate())).real
        #         fcc[i]=sum1/sum2

        return fcc


    #def fcc_wight(self,half1,half2,whole,points,dist):
    #    # before_bin = timeit.default_timer()
    #    x0=points[0]
    #    y0=points[1]
    #    z0=points[2]
    #    ftsize=half1.shape[0]/2
    #    fcc=numpy.zeros(ftsize,dtype=numpy.float64)
    #    fc = (half1[x0,y0,z0] * half2[x0,y0,z0].conjugate()).real
    #    # mid_bin = timeit.default_timer()
    #    dist10 = numpy.round(numpy.linalg.norm(points - dist.shape[0]/2, axis=0)).astype(int)

    #    # dist10 = numpy.round(numpy.sqrt(x2**2+y2**2+z2**2)).astype(int)

    #    # start_bin = timeit.default_timer()       
    #    sum1_bin = numpy.bincount(dist10, fc, minlength=0)
    #    sum2_bin = numpy.sqrt(numpy.bincount(dist10, (half1[x0,y0,z0]*half1[x0,y0,z0].conjugate()).real, minlength=0) *\
    #                          numpy.bincount(dist10, (half2[x0,y0,z0]*half2[x0,y0,z0].conjugate()).real, minlength=0))
    #    # stop_bin = timeit.default_timer()
    #    # print "before_bin", start_bin - before_bin
    #    # print "mid_bin", mid_bin - before_bin
    #    # print "mid2_bin", start_bin - mid_bin
    #    # print "bin_time", stop_bin - start_bin
    #    # print "sum1_bin", sum1_bin[:ftsize].size, sum1_bin[0:5]
    #    # print "sum2_bin", sum2_bin[:ftsize].size, sum2_bin[0:5]
    #    sum1_bin[0:3][sum1_bin[0:3] == 0.] = 1          # there may be few points in this small shell with radius less than 3
    #    sum1_bin[0:3][sum1_bin[0:3] <= 0.] = numpy.abs(sum1_bin[0:3][sum1_bin[0:3] <= 0.])
    #    sum2_bin[0:3][sum2_bin[0:3] == 0.] = 1
    #    fcc = sum1_bin[:ftsize] / sum2_bin[:ftsize]

    #    return fcc


    def set_output_directory(self, output_directory=None, project_root='.'):
        '''
        Set output directory
        '''

        if output_directory is not None:
            self.output_directory = output_directory
        else:
            # Get project root
            head = project_root

            # Directory head
            dir_head = 'YifanFsc'

            # List existing output directories
            potential_directories = list(filter(lambda x: os.path.isdir(x),
                                         glob.glob(head+'/'+dir_head+'_em_[0-9][0-9][0-9]')))

            # Get the number extensions
            number_extensions = [int(x[-3:]) for x in potential_directories]

            # Get the counter
            output_counter = 1
            if len(number_extensions) > 0:
                output_counter = max(number_extensions)+1

            self.output_directory = head+'/'+dir_head+"_em_%03d" % (output_counter)


        # Make directory
        os.mkdir(self.output_directory)

    # wrapping FFT ----------------------------------------------------------
    def fft_wrap(self,array):
        x=array.shape[0]
        y=array.shape[1]
        z=array.shape[2]
        z=2*(z-1)
        x2=x/2
        y2=y/2
        z2=z/2

        dummy=numpy.zeros((x,y,z),dtype=numpy.complex128)
        dummy[:x2,:y2,z2+1:]=array[x2:,y2:,1:z2]
        dummy[:x2,y2:,z2+1:]=array[x2:,:y2,1:z2]
        dummy[x2:,:y2,z2+1:]=array[:x2,y2:,1:z2]
        dummy[x2:,y2:,z2+1:]=array[:x2,:y2,1:z2]
        arrayflip=array[::-1,::-1,::-1]
        dummy[x2:,y2:,:z2+1]=arrayflip[:x2,:y2,:]
        dummy[:x2,y2:,:z2+1]=arrayflip[x2:,:y2,:]
        dummy[x2:,:y2,:z2+1]=arrayflip[:x2,y2:,:]
        dummy[:x2,:y2,:z2+1]=arrayflip[x2:,y2:,:]
        return dummy

    # unwrapping FFT --------------------------------------------------------
    def fft_unwrap(self,array):
        x=array.shape[0]
        y=array.shape[1]
        z=array.shape[2]
        x2=x/2
        y2=y/2
        z2=z/2
        array=array[:,:,:z2+1]
        arrayflip=array[::-1,::-1,::-1]
        dummy=numpy.zeros((x,y,z2+1),dtype=numpy.complex128)
        dummy[x2:,y2:,:]=arrayflip[:x2,:y2,:]
        dummy[:x2,y2:,:]=arrayflip[x2:,:y2,:]
        dummy[x2:,:y2,:]=arrayflip[:x2,y2:,:]
        dummy[:x2,:y2,:]=arrayflip[x2:,y2:,:]
        return dummy

    # mask for Fourier cone correlation -------------------------------------
    def fcc_mask(self,mask,points,fcc):

        # print "mask_fcc_before", fcc

        # coords=numpy.where(fcc<0)
        # fcc[:coords[0][0]]=1

        # print "mask_coords", type(coords), coords[0][0], coords
        # print "mask_fcc_after", fcc

        x0=points[0]
        y0=points[1]
        z0=points[2]

        dist10 = numpy.round(numpy.linalg.norm(points - mask.shape[0]/2, axis=0)).astype(int)
        # print "dist10", dist10.size, dist10

        for i in range(mask.shape[0]/2):
            # if fcc[i] < 0 :
            #     fcc[i] = -fcc[i]
            if  fcc[i]>=0.143:
                coords=numpy.where(dist10==i)
                x1=x0[coords[0]]
                y1=y0[coords[0]]
                z1=z0[coords[0]]
                if x1 != []:
                    # print "coords_max", numpy.amax(y1), coords[0]
                    mask[x1, y1, z1] = fcc[i]    
        return mask

    def fcc_Cref_mask(self,mask,points,fcc):

        x0=points[0]
        y0=points[1]
        z0=points[2]

        dist10 = numpy.round(numpy.linalg.norm(points - mask.shape[0]/2, axis=0)).astype(int)

        for i in range(mask.shape[0]/2):
            # if fcc[i] < 0 :
            #     fcc[i] = -fcc[i]
            if  fcc[i]>=0.143:
                coords=numpy.where(dist10==i)
                x1=x0[coords[0]]
                y1=y0[coords[0]]
                z1=z0[coords[0]]
                if x1 != []:
                    # print "coords_max", numpy.amax(y1), coords[0]
                    mask[x1, y1, z1] = numpy.sqrt(2*fcc[i]/(1+fcc[i]))
        return mask

    def bfactor_compute(self, map_ft, points, fcc, apix, dist, fit_minres, fit_maxres):
        x=points[0]
        y=points[1]
        z=points[2]
        radius = numpy.round(numpy.sqrt(dist[x,y,z])).astype(int)
        count = numpy.bincount(radius, minlength=0)
        structureF = numpy.bincount(radius, numpy.absolute(map_ft[x,y,z]), minlength=0)[:fcc.size]
        res = numpy.arange(float(fcc.size))
        res[1:] = (map_ft.shape[0] * apix) / res[1:]
        index_fit = numpy.where((res <= fit_minres) & (res >= fit_maxres) & (res >= apix*2) & (structureF > 0))
        lnF_fit = numpy.log(structureF[index_fit] / count[index_fit])
        slope, intercept = numpy.polyfit(1/res[index_fit]**2, lnF_fit, 1)
        bfactor = 4. * slope
        print("slope, intercept, bfactor", slope, intercept, bfactor)
        return bfactor

    def bfactor_compute_save(self, map_ft, points, fcc, apix, dist, fit_minres, fit_maxres, file_name):
        x=points[0]
        y=points[1]
        z=points[2]
        radius = numpy.round(numpy.sqrt(dist[x,y,z])).astype(int)
        count = numpy.bincount(radius, minlength=0)
        structureF = numpy.bincount(radius, numpy.absolute(map_ft[x,y,z]), minlength=0)[:fcc.size]
        res = numpy.arange(float(fcc.size))
        res[1:] = (map_ft.shape[0] * apix) / res[1:]
        index_fit = numpy.where((res <= fit_minres) & (res >= fit_maxres) & (res >= apix*2) & (structureF > 0))
        res_fit = 1/res[index_fit]**2
        lnF_fit = numpy.log(structureF[index_fit] / count[index_fit])
        print("structureF", structureF)
        print("res", res)
        print("index_fit", index_fit)
        print("res_fit", res_fit)
        print("lnF_fit", lnF_fit)
        slope, intercept = numpy.polyfit(res_fit, lnF_fit, 1)
        # plt.plot(res_fit, lnF_fit, 'm.-', label='lnF')
        # plt.plot(res_fit, res_fit * slope + intercept, 'b')
        if apix==0:
            plt.axis([0, fcc.shape[1], -10, 10])
            plt.xlabel('Spatial frequency (pixels)')
        else:
            plt.axis([0., 1./(2.*apix)**2, 0, 5])
            #plt.gca().invert_xaxis()
            plt.xlabel('Spatial frequency (1/Angstroms^2)')
        plt.savefig(file_name)
        bfactor = 4. * slope
        print("slope, intercept, bfactor", slope, intercept, bfactor)
        return bfactor

    # def bfactor_mask(self, mask, points, bfactor, apix):

    #     x0=points[0]
    #     y0=points[1]
    #     z0=points[2]

    #     dist10 = numpy.round(numpy.linalg.norm(points - mask.shape[0]/2, axis=0)).astype(int)

    #     for i in range(mask.shape[0]/2):
    #         # if fcc[i] < 0 :
    #         #     fcc[i] = -fcc[i]
    #         coords=numpy.where(dist10==i)
    #         x1=x0[coords[0]]
    #         y1=y0[coords[0]]
    #         z1=z0[coords[0]]
    #         if x1 != []:
    #             # print "coords_max", numpy.amax(y1), coords[0]
    #             mask[x1, y1, z1] = numpy.sqrt(2*fcc[i]/(1+fcc[i]))
    #     return mask

# run program ---------------------------------------------------------------
Main()
print('dFSC calculation complete!')
