#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import tifffile
from copy import deepcopy


import skimage
from skimage.segmentation import find_boundaries
from scipy import ndimage as ndi

from . import balloon
# import sys
# sys.path.append('./src/PombeTrack/library')
# #  sys.path.append('../PombeTrack/library')
# import balloon

# #  Load image
# def load_brightfield(image_path):
#     im_frames=tifffile.TiffFile(image_path)
#     z_stack=np.zeros(im_frames.pages[0].asarray().shape +\
#             (len(im_frames.pages[::3]),), dtype=int)
#     for ii in range(0, z_stack.shape[2]):
#         z_stack[:,:,ii] = im_frames.pages[3*ii].asarray()
#     z_stack=skimage.img_as_uint(z_stack)

#     #  im_mid=z_stack[:,:,int(np.floor(z_stack.shape[2]/2))]
#     #  im_up=z_stack[:,:,int(np.floor(z_stack.shape[2]/2)-1)]
#     #  im=np.maximum(im_mid,im_up)

#     return z_stack



#  Preprocessing
def preprocessing(im, eq_hist=True):

    #histogram equalisation
    # im_eq=exposure.ualize_hist(im)
    if eq_hist:
        im_eq = skimage.exposure.equalize_adapthist(skimage.filters.gaussian(im,2), clip_limit=0.1)
    else:
        im_eq=skimage.filters.gaussian(im,1)/skimage.filters.gaussian(im,1).max()
    # im_eq=im/im.max()

    # mean filter to remove shading
    im_mf=skimage.filters.rank.mean(im_eq,skimage.morphology.disk(40))
    im_pp=im_eq-im_mf/255
    return im_pp

#  four filters for finding edges
def find_edges(im):
    f1=np.array([[-1,-1,-1],[2,2,2],[-1,-1,-1]])
    f2=np.array([[-1,2,-1],[-1,2,-1],[-1,2,-1]])
    f3=np.array([[2,-1,-1],[-1,2,-1],[-1,-1,2]])
    f4=np.array([[-1,-1,2],[-1,2,-1],[2,-1,-1]])


    im_f1=ndi.convolve(im, f1, mode='constant', cval=0.0)
    im_f2=ndi.convolve(im, f2, mode='constant', cval=0.0)
    im_f3=ndi.convolve(im, f3, mode='constant', cval=0.0)
    im_f4=ndi.convolve(im, f4, mode='constant', cval=0.0)

    im_f=(im_f1**2+im_f2**2+im_f3**2+im_f4**2)**(1/2)
    return im_f/im_f.max()

#  Otsu thresholding to find cell interior
def find_cellinterior(im_pp):
    # Gradient
    im_gr=find_edges(im_pp)

    # Otsu thresholding
    ot=skimage.filters.threshold_otsu(im_gr)
    miu_b=np.mean(im_pp[im_gr<ot])

    # Otsu for membrane
#     ot_m=threshold_otsu(im_pp[im_pp<miu_b])

    # Otsu for interior
    ot_i=skimage.filters.threshold_otsu(im_pp[im_pp>miu_b])


    # Regions
#     im_m=im_pp<ot_m
    im_i=im_pp>ot_i
#     im_b=(im_pp>=ot_m)*(im_pp<miu_b)

    # Close opening
    im_i=skimage.morphology.binary_closing(im_i)
    im_i=skimage.morphology.binary_opening(im_i)
    im_i=skimage.morphology.remove_small_holes(im_i,1000)
    im_i=skimage.morphology.remove_small_objects(im_i,500,connectivity = 1)
    return im_i


#  Use watershed algorithm to find boundaries between connected cells
def find_watershed(im_i):
    # Watershed
    im_dt=ndi.distance_transform_edt(im_i)

    local_maxi = im_dt>6
    local_maxi=skimage.morphology.remove_small_objects(local_maxi,100,connectivity = 1)

    markers = ndi.label(local_maxi)[0]
    im_wat=skimage.morphology.watershed(-im_dt, markers, mask=im_i)
    return im_wat


# Find boundaries
def sort_in_radial(im):
    bd_ii=np.transpose(np.asarray(im.nonzero()))
    bd_ii_mean=[np.mean(bd_ii[:,0]).astype(int), np.mean(bd_ii[:,1]).astype(int)]
    bd_ii_recentered=bd_ii-bd_ii_mean
    bd_z=bd_ii_recentered[:,0]+1j*bd_ii_recentered[:,1]
    bd_ang=np.angle(bd_z)
    bd_ii_ang=np.concatenate((bd_ii,bd_ang.reshape((bd_ang.shape[0],1))),axis=1)
    bd_ii_sorted=bd_ii_ang[bd_ii_ang[:,2].argsort()[::-1]]
    return bd_ii_sorted


def sort_in_order(im):
    bd_ii=np.transpose(np.asarray(im.nonzero()))
    bd_ii_sorted=np.zeros(bd_ii.shape)
    bd_ii_sorted[0,:]=bd_ii[0,:]
    bd_ii=bd_ii[~np.all(bd_ii==bd_ii[0,:],1),:]
    for i in range(1,len(bd_ii_sorted)):
        dist_list=np.linalg.norm(bd_ii-bd_ii_sorted[i-1,:],axis=1)
        bd_ii_sorted[i,:]=bd_ii[((dist_list==np.min(dist_list)).nonzero()[0])[0],:]
        # remove selected item from bd_ii
        bd_ii=bd_ii[~np.all(bd_ii==bd_ii_sorted[i,:],1),:]
    return bd_ii_sorted


def find_bd(im_wat):
    bd = []
    for ii in range(1,im_wat.max()+1):
        im_ii=im_wat==ii
        if not (np.any(np.asarray(im_ii.nonzero())==0) or
                np.any(np.asarray(im_ii.nonzero())==2047)):
            im_ii_bd=find_boundaries(im_ii,mode='inner')
            # Sort in radial
            bd_ii_sorted=sort_in_order(im_ii_bd)
            bd.append(bd_ii_sorted[:,0:2].astype(int))
    return bd





# # Define Classes for cells and images
# class Cell(object):
#     def __init__(self, balloon_obj, origin_y, origin_x, halfwidth):
#         self.coordinates = balloon_obj.get_coordinates(accept=True) + [origin_y, origin_x]
#         self.origin_y=origin_y
#         self.origin_x=origin_x
#         self.halfwidth=halfwidth
#         self.area = balloon_obj.get_area()
#         self.length=0
#         self.width=0
#         self.GFP_sig=0
#         self.RFP_sig=0
#         self.nuclei_mask=np.array([])
#         self.nuclei_num=0
#         self.nuclei_sig=0

# class Image(object):
#     def __init__(self, impath, base_image, CellList):
#         self.experiment_ID = impath[0]
#         self.movie_ID=impath[1]
#         self.position_ID=impath[2]
#         self.base_image = base_image
#         self.CellList = CellList
#     def plotall(self):
#         f = plt.figure(figsize=(20,20))
#         ax = f.add_subplot(111)
#         ax.imshow(self.base_image, cmap="gray")
#         ax.set_aspect("equal")
#         for ii in range(0, len(self.CellList)):
#             coor = (self.CellList[ii]).coordinates
#             ax.plot(coor[:, 1], coor[:, 0], marker="None", color="r", linewidth=0.5)
#             ax.text(np.mean(coor[: ,1]), np.mean(coor[:, 0]), str(ii), fontsize = 14, 
#                     horizontalalignment='center', verticalalignment='center', color = 'k')
#         # ax.set_title("Auto Segmentation")
#         return f
#     def rm(self, index_CellList):
#         self.CellList=np.delete(self.CellList, index_CellList)
#         return self.plotall()
#     def ad(self, cell_obj):
#         self.CellList=np.append(self.CellList, cell_obj)
#         return self.plotall()














# Find balloon object in a small window
def find_balloon_obj(edges, base_image):
   halfwidth=int(np.amax(np.amax(edges,0)-np.amin(edges,0))/2)+50
   # Centre of the cell
   cy=np.mean(edges[:,0]).astype(int)
   cx=np.mean(edges[:,1]).astype(int)
#
   # Set initial origin and Find offset
   oy=cy-halfwidth
   ox=cx-halfwidth
   if oy<0:
       offset_y=oy
   elif oy>2047:
       offset_y=oy-2047
   else:
       offset_y=0
   if ox<0:
       offset_x=ox
   elif ox>2047:
       offset_x=ox-2047
   else:
       offset_x=0
#
   # Find origin(top left corner) for window
   origin_y=oy-offset_y
   origin_x=ox-offset_x
#
   # Set initial nodes
   init_nodes=deepcopy(edges)
   init_nodes[:,0]=init_nodes[:,0]-origin_y
   init_nodes[:,1]=init_nodes[:,1]-origin_x
#
   # Set im_window and balloon object
   im_window=base_image[origin_y:origin_y+2*halfwidth,origin_x:origin_x+2*halfwidth]
#
   balloon_obj=balloon.Balloon(init_nodes,im_window)
   return balloon_obj, origin_y, origin_x, halfwidth
