#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy import units as u
import csv
from sedfitter.filter import Filter
from sedfitter.convolve import convolve_model_dir

"""
Created on Tue May  7 18:42:40 2019

@author: mspovich

Convolve a set of model SEDs for the Robitaille sedfitter with 1 or more new 
filter profiles downloaded from 
http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=UKIRT/UKIDSS.J&&mode=browse&gname=UKIRT&gname2=UKIDSS#filter
Central wavelength is also listed on SVO Filter Profiles as lambda_cen

Here's an example for the Johnson UBVRI filters:
fprofile = ['/Users/mspovich/research/ysomodels/Generic_Bessell.U.dat.txt',
            '/Users/mspovich/research/ysomodels/Generic_Bessell.B.dat.txt',
            '/Users/mspovich/research/ysomodels/Generic_Bessell.V.dat.txt',
            '/Users/mspovich/research/ysomodels/Generic_Bessell.R.dat.txt',
            '/Users/mspovich/research/ysomodels/Generic_Bessell.I.dat.txt'] 

flambda = [0.360452,0.435833,0.544579,0.643647,0.805453]

fname = ['BU','BB','BV','BR','BI']
"""

#Model sets that I typically use
#model_dir = '01_sp--s-i'
#model_dir = '02_sp--h-i'
#model_dir = '14_spu-smi'
#model_dir = '15_spu-hmi'
#model_dir = '16_spubsmi'
#model_dir = '17_spubhmi'
#model_dir = 'models_pms'

#

def convolve_new_filt(model_dir, fprofile, flambda, fname='NewFilt'):

    c = 2.99792e8

    for i in range(0,len(fprofile)):
        with open(fprofile[i]) as pf:
            fprofreader = csv.reader(pf, delimiter=' ', quotechar='|', quoting=csv.QUOTE_NONNUMERIC)
            fprofdata = [r for r in fprofreader]

        fu = Filter()

        fu.name = fname[i] 
        fu.central_wavelength = flambda[i] * u.micron

        n = len(fprofdata)
        nu = [c / (fprofdata[j][0] * 1.e-10) for j in range(n)]
        nu.reverse()
        fu.nu = nu * u.Hz
        response = [fprofdata[j][1] for j in range(n)]
        response.reverse()
        fu.response = response

        fu.normalize()
    
        convolve_model_dir(model_dir, [fu])

