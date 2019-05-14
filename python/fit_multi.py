#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy import units as u
from sedfitter.extinction import Extinction
from sedfitter import fit
from sedfitter import filter_output
from sedfitter import write_parameters

def fit_multiyso(data, filter_names, apertures, models_topdir, n_data_min=3,
        extinction_file=None, av_range=None, distance_range=None,
        output_format=('F', 6.), output_convolved=False,
        remove_resolved=False, cpd=4.):
    """
    Created on Wed May  8 11:17:09 2019

    @author: mspovich

    Call the sedfitter fit function multiple times to fit various R17 YSO 
    model sets to a given source sample.
    
    All arguments in call are passed along to the fit module.
    """

    # Extinction law
    extinction = Extinction.from_file(extinction_file, 
                                      columns=[0, 3], wav_unit=u.micron, 
                                      chi_unit=u.cm**2 / u.g)
    
    # DISK-only models
    model_dir_spsi = models_topdir + '/01_sp--s-i' 
    model_dir_sphi = models_topdir + '/02_sp--h-i' 

    # Disks with inner holes at Rsublimation..
    fit(data, filter_names, apertures, model_dir_spsi, '01_sp--s-i.fitinfo', 
        n_data_min=3, extinction_law=extinction, 
        distance_range=distance_range, av_range=av_range, 
        output_format=output_format)

    # Disks with varying inner hole size.. 
    fit(data, filter_names, apertures, model_dir_sphi, '02_sp--h-i.fitinfo', 
        n_data_min=3, extinction_law=extinction,
        distance_range=distance_range, av_range=av_range, 
        output_format=output_format)
    
    # DISK+ENVELOPE models
    model_dir_spusmi = models_topdir + '/14_spu-smi' 
    model_dir_spuhmi = models_topdir + '/15_spu-hmi'  

    # Disks+Ulrich envelopes with inner holes 
    fit(data, filter_names, apertures, model_dir_spusmi, '14_spu-smi.fitinfo', 
        n_data_min=3, extinction_law=extinction, 
        distance_range=distance_range, av_range=av_range, 
        output_format=output_format)
        
    fit(data, filter_names, apertures, model_dir_spuhmi, '15_spu-hmi.fitinfo', 
        n_data_min=3, extinction_law=extinction, 
        distance_range=distance_range, av_range=av_range, 
        output_format=output_format)
    
    # DISK+ENVELOPE+BipolarCavity models
    model_dir_spubsmi = models_topdir + '/16_spubsmi' 
    model_dir_spubhmi = models_topdir + '/17_spubhmi' 

    # Disks+Ulrich envelopes with inner holes and bipolar cavities
    fit(data, filter_names, apertures, model_dir_spubsmi, '16_spubsmi.fitinfo', 
        n_data_min=3, extinction_law=extinction, 
        distance_range=distance_range, av_range=av_range, 
        output_format=output_format)
     
    fit(data, filter_names, apertures, model_dir_spubhmi, '17_spubhmi.fitinfo', 
        n_data_min=3, extinction_law=extinction, 
        distance_range=distance_range, av_range=av_range, 
        output_format=output_format)
         
    # Splitting output to identify well-fit sources.
    print("Splitting output into well-fit vs. poorly-fit using cpd = %4.1f" % cpd)
    filter_output('01_sp--s-i.fitinfo', cpd=cpd) 
    filter_output('02_sp--h-i.fitinfo', cpd=cpd)
    filter_output('14_spu-smi.fitinfo', cpd=cpd) 
    filter_output('15_spu-hmi.fitinfo', cpd=cpd)    
    filter_output('16_spubsmi.fitinfo', cpd=cpd) 
    filter_output('17_spubhmi.fitinfo', cpd=cpd) 

    #Write out parameters of WELL-FIT models from each set 
    #for subsequent cross-set analysis
    cpdstr = '{:02.0f}'.format(cpd)
    print("Writing ASCII parameter files pars_*g4.txt")
    write_parameters('01_sp--s-i.fitinfo_good','pars_01g'+cpdstr+'.txt',
                     select_format=('A',-1)) 
    write_parameters('02_sp--h-i.fitinfo_good','pars_02g'+cpdstr+'.txt',
                     select_format=('A',-1)) 
    write_parameters('14_spu-smi.fitinfo_good','pars_14g'+cpdstr+'txt',
                     select_format=('A',-1)) 
    write_parameters('15_spu-hmi.fitinfo_good','pars_15g'+cpdstr+'.txt',
                     select_format=('A',-1)) 
    write_parameters('16_spubsmi.fitinfo_good','pars_16g'+cpdstr+'.txt',
                     select_format=('A',-1)) 
    write_parameters('17_spubhmi.fitinfo_good','pars_17g'+cpdstr+'.txt',
                     select_format=('A',-1)) 

    