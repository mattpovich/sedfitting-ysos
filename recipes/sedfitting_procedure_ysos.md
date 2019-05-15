# `sedfitting_procedure_ysos.txt`

**Authors: M. S. Povich & T. P. Robitaille**

## Description

This recipe guides you through the steps of fitting multiple SED models from [Robitaille (2017)](https://doi.org/10.1051/0004-6361/201425486) to combined *Spitzer*, 2MASS, and UKIDSS (or VVV) candidate YSOs in a specified target field.

This pipeline creates two final science products:
1. A summary table of the YSOs and SED fit classifications called `<target>.ised.fits`.
1. A directory called `results_ysoc` containing IDL save files of the SED model parameters for each candidate YSO. 

## Version History

* Original sedfitting_procedure_ukidss-mir — 2011-08-24 by M. S. Povich 
* Updated to use the [`python sedfitter`](https://github.com/astrofrog/sedfitter) — 2018-07-23 by M. S. Povich 
* CURRENT (v1.0) — 2019-05-15 by M. S. Povich

I recommend using two different terminal windows simultaneously, one in `bash` and the other in `tcsh`. The correct shell in which to execute each block of command is indicated as follows:
 
  * `python` command blocks are preceded by **>>>**
  * `tcsh/IDL` command blocks are preceded by **%/IDL>**


## INITIAL SETUP

Get Tom Robitaille's sedfitter software and install it following the instructions
at http://sedfitter.readthedocs.io/en/stable/installation.html

Download the following 6 SED model sets from https://doi.org/10.5281/zenodo.166732 (there are 17 model sets total; these are the only ones used by this pipeline), expane the `*tar.gz` archives and —> prepend the model set number to each directory

* `sp--s-i.tar.gz —> 01_sp--s-i.tar` 
* `sp--h-i.tar.gz`
* `spu-smi.tar.gz`
* `spu-hmi.tar.gz`
* `spubsmi.tar.gz`
* `spubhim.tar.gz`

Make sure you have the following libraries on your IDL path:

* [IDL astronomy users library](https://github.com/wlandsman/IDLAstro)
* [IDL coyote library](https://github.com/idl-coyote/coyote)
* [IDL library bundled with this recipe](https://github.com/mattpovich/sedfitting-phrds)


## TARGET-SPECIFIC SETUP

**>>>**

    TARGET='m17swex'  # python EXAMPLE. **Choose your own target name.**


**%**

    setenv TARGET m17swex  # tcsh EXAMPLE. **Choose your own target name (must match above)**


It is a good idea to name your working directory `$TARGET/sedfitter`, e.g. (in `tcsh`):

**%**

    mkdir -p $TARGET/sedfitter  #SKIP if directory exists from previous recipe!
    cd $TARGET/sedfitter
    
[Download](https://irsa.ipac.caltech.edu/data/SPITZER/MIPSGAL/) a **MIPSGAL 24 µm mosaic** covering your target field-of-view and name it `../$TARGET.mm24.fits`.

**Optional.** If you are interested in flagging YSOs found within a specfied "cropped" sub-region of your target catalog, make an `SAOImage ds9` regionfile containing one or more `polygon` regions and save it as `../$TARGET.xfov.reg` (required coordinate format is fk5 DECIMAL DEGREES). I personally use this feature to identify YSOs falling within the *Chandra*/ACIS field-of-view within a larger *Spitzer* mosaic.

####SED FITTING STARTS HERE####

mkdir sedfitter
cd sedfitter





	 




On MacOS, I find I can just page through the large set of *pdf files made in a finder window with the preview pane made large.  

Reset 4.5 um fluxes to upper limits for suspected excess sources (outflow candidates; Povich et al. 2011)

  %
  IDL>
        ugoflag,'data_glimpse+ysoc','data_glimpse+ysoc_ugos',nwav=nwav,/regions		
---------------------
Get MIPSGAL 24 um photometry for strong IRE sources. Using MUCH LESS aggressive error resetting than previously to avoid overly downweighting 24 um compared to NIR+IRAC datapoints.
NOTE: It is a good idea to review the header of the MIPS mosaic to check if information is available for mean exposure time and zodaical background subtraction, otherwise just use the default keyword values.

   IDL>
	.r mipsfitapphot
	 spawn,'echo $TARGET', target
 	 mipsmosaic='../'+target+'.mm24.fits'  ;This is your MIPS 24 um mosaic covering the field of view.
  	 data_in='data_glimpse+ysoc_ugos'
  	 data_out='data_glimpse+ysoc_ugos_24um'
  	 photlist='mipsmerger.txt'
  	 regfile='mipsmerger_apertures.reg'
	 snlim = 5. ;Other values are possible, experiment if desired
	 mfe = 0.10 ;Recommend NOT going lower than this.
 	 mipsfitapphot,mipsmosaic,data_in,data_out,photlist,apreg=regfile,aprad=3.5,anrad=7.,snlim=snlim,min_frac_errs=mfe,nwav=nwav

OPTIONAL: Review the extraction results, iterate the above with modified values of aprad, anrad, or snlim:

  ds9 ../$TARGET.mm24.fits -region mipsmerger_ysoc_apertures.reg &

#OPTIONAL: If another iteration is desired, archive previous as, e.g.:
   mkdir mipsphot_sn7
   cp mipsmerger* mipsphot_sn7

I'm sticking with min S/N = 5 (it's never gonna be perfect).

---------------------------------
Fit the YSO candidates list with R17 YSO model sets.

## NOTE: The R17 models do not come pre-convolved with UKIDSS and VVV filter
   profiles. You can download the SVO filter profiles from
   http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php
   and colve multiple filters with R17 models in one shot
   using new_filt.py.

   I recommend placing new_filt.py in your top models directory
   and fit_multi.py
   in your sedfitter working directory.

>>>

from fit_multi import fit_multiyso
extinction24_file = '/Users/mspovich/research/ysomodels/ex_law_mc.par'
models_topdir = '/Users/mspovich/research/ysomodels'

data = 'data_glimpse+ysoc_ugos_24um'

#UKIDSS example
filters = ['UDSSJ', 'UDSSH', 'UDSSK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 'M1'] 
apertures = [2., 2., 2., 3., 3., 3., 3., 3., 3., 3., 7.] * u.arcsec

#VVV example
filters = ['VVVZ', 'VVVY', 'VVVJ', 'VVVH', 'VVVK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 'M1'] 
apertures = [2., 2., 2., 2., 2., 3., 3., 3., 3., 3., 3., 3., 7.] * u.arcsec

fit_multiyso(data, filters, apertures, models_topdir,
		   extinction_file=extinction24_file,
		   n_data_min=3,  output_format=('F',6),
		   distance_range=[2.9, 3.3] * u.kpc,
		   av_range=[0.,100.], cpd=9.)

#OPTIONAL: Plot up the good (and bad) SED fits (FOR EXAMPLE):
from sedfitter import plot
plot('01_sp--s-i.fitinfo_good', 'plots_01g4') 
plot('01_sp--s-i.fitinfo_bad', 'plots_01b4') 
plot('02_sp--h-i.fitinfo_good', 'plots_02g4') 
plot('02_sp--h-i.fitinfo_bad', 'plots_02b4') 


##OPTIONAL: Parameter plots
from sedfitter import plot_params_2d

##EXAMPLE of a 2-d parameter plot, inclination vs. disk mass to suss out physically improbable models (very massive disks, inclinations near face- or edge-on). 
plot_params_2d('02_sp--h-i.fitinfo_good', 'disk.mass', 'inclination', 'plots_02g4_mdisk_incl', log_x=True, log_y=False, select_format=('A', -1)) 

##EXAMPLE of a 2-d parameter plot, key STELLAR parameters. 
plot_params_2d('02_sp--h-i.fitinfo_good', 'star.temperature', 'star.radius', 'plots_02g4_reff_rstar', log_x=True, log_y=True, select_format=('A', -1)) 


#Send all well-fit models into IDL save files, in the results_ysoc
 directory.
 Combines all model sets fit into SINGLE parameter tables for each source.

%
   IDL>
   data = 'data_glimpse+ysoc_ugos_24um'
   cpdlab = '04'  ;String label for cpd cut given to fit_multiyso in Python
   spawn,'ls pars_*g'+cpdlab+'.txt',infiles 
   r17_pars2idl, infiles, 'results_ysoc', data_parent=data,filter='F6'

%
 # Once you've successfully made the IDL save files, free up disk
 space! EXAMPLE:
 
tar cvzf pars_r17g04.tgz pars_*g04.txt  #if you want to keep these at all!
'rm' pars_*g04.txt  

=============================================================================
CREATE FITS TABLE SUMMARIZING RESULTS  ($TARGET.ised.fits)

Included steps:
FLAG CANDIDATE GALAXIES AND AGN (galflag.pro)
FLAG CANDIDATE AGB STARS
ASSIGN YSO STAGES (stageflag_r17.pro)

NOTE: The regionfile *xfov.reg must contain one or more non-overlapping polygons in the fk5 decimal degree coordinate system.

  idl
     xfov = '../'+target+'.xfov.reg'	
     create_ised_r17, xfov=xfov
  exit

  ln ised_r17.fits ../$TARGET.ised.fits
  ln ised_stageflag.reg ../$TARGET.stageflag.reg

END OF RECIPE



   
