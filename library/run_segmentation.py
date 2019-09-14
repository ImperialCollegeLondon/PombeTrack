#!/usr/bin/env python3
import sys
sys.path.append('./src/ImageAnalysis')
from segmentation import *
from analysis import *

import numpy as np
import matplotlib.pyplot as plt

import time

import pickle

from copy import deepcopy

import skimage


###########################################################################
# Load
experiment_ID=sys.argv[1]
movie_ID=sys.argv[2]
position_ID=sys.argv[3]

#  experiment_ID='2019-02-11-Ibidi-03'
#  movie_ID=3
#  position_ID=18

image_path='/data/alvin/imagedata/'+experiment_ID+'/movie_'\
        +str(movie_ID)+'/movie_'+str(movie_ID)+'_MMStack_Pos'+str(position_ID)+'.ome.tif'


# Find rough boundaries
z_stack=load_brightfield(image_path)

im_mid=z_stack[:,:,int(np.floor(z_stack.shape[2]/2))]
im_up=z_stack[:,:,int(np.floor(z_stack.shape[2]/2)-1)]
im=np.maximum(im_mid,im_up)

im_pp=preprocessing(im)
im_i=find_cellinterior(im_pp)
im_wat=find_watershed(im_i)
bd=find_bd(im_wat)


# Create CellList
CellList = np.array([])

for index in range(0, len(bd)):
    balloon_obj, origin_y, origin_x, halfwidth=find_balloon_obj(bd[index][::5], im)

    # Evolve the contour
    try:
#         sensitivity=0.5
#         for i in range(100):
#             balloon_obj.evolve(display=False,image_percentile=sensitivity)
        sensitivity=0.4
        area_init=balloon_obj.get_area()
        for i in range(200):
            balloon_obj.evolve(display=False,image_percentile=sensitivity)
            if balloon_obj.get_area()>1.5*area_init or balloon_obj.get_area()<0.5*area_init:
                raise ValueError()
        CellList=np.append(CellList, Cell(balloon_obj, origin_y, origin_x, halfwidth))
    except:
        pass

# Create Image object
Im_obj=Image([experiment_ID, movie_ID, position_ID], im_mid, CellList)
fig=Im_obj.plotall()




#########################################################################
# Analysis
#  GFP_im,RFP_im = load_fluorescence(image_path)
#  RFP_pp=preprocessing(RFP_im,False)
#  for idd in range(len(Im_obj.CellList)):
    #  coords=Im_obj.CellList[idd].coordinates
    #  Signal
    #  Im_obj.CellList[idd].GFP_sig,_=get_signals(coords,GFP_im)
    #  Im_obj.CellList[idd].RFP_sig,_=get_signals(coords,RFP_im)
    #  Nulcei
    #  nu,_=find_nuclei(coords,RFP_pp)
    #  Im_obj.CellList[idd].nuclei_mask=nu
    #  Im_obj.CellList[idd].nuclei_num=nu.max()
    #  Im_obj.CellList[idd].nuclei_sig,_=get_nuclei_signals(nu, RFP_im)
    #  Shape:
    #  Im_obj.CellList[idd].length,Im_obj.CellList[idd].width = find_shape(coords, RFP_pp)


#############################################################################
# Save outputs
# Saving images:
output_path='./data/segmentation_outputs_2/'
fig.savefig(output_path+'experiment_'+experiment_ID+'_movie_'+str(movie_ID)+'_Pos'+str(position_ID)+'.png')
# Saving the objects:
with open(output_path+'experiment_'+experiment_ID+'_movie_'+str(movie_ID)+'_Pos'+str(position_ID)+'.pkl', 'wb') as file:  # Python 3: open(..., 'wb')
    pickle.dump(Im_obj, file)

