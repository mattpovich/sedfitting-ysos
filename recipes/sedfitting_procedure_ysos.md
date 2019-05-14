sedfitting_procedure_ukidss-mir_PYTHON.txt
Original sedfitting_procedure_ukidss-mir 2011-08-24 by M. S. Povich
UPDATED to use PYTHON sedfitter software 2018-07-23 by M. S. Povich
UPDATED to take Pat's new GLIMPSE+ matched tables as input 2019-05-03
by M. S. Povich

This recipe guides you through the steps of performing a basic SED
fitting classification for the combined Spitzer, 2MASS, and UKIDSS (or
VVV) catalog sources in a specified target field. UNLIKE the original MYStIX version, this assumes you have a MIPS 24 um mosaic available as well!

The final science product is an IDL save file called <target>.ised.sav. 

INITIAL SETUP

Get Tom's sedfitter software and install it following the instructions
at http://sedfitter.readthedocs.io/en/stable/installation.html

#I am using this version of python, which I need to
 run in a bash shell:

 Python 3.6.5 :: Anaconda, Inc.
 
 IDL commands will need to be executed in tcsh,
 so we require two different terminal windows open in different shells.
 
  ### Python command blocks are preceded by ">>>"
  ### tcsh/IDL command blocks are preceded by "%/IDL>"

>>>
TARGET='m17swex'  #FOR EXAMPLE

%
setenv TARGET m17swex

TARGET-SPECIFIC SETUP

#Link to the IRAC and MIPS mosaics for the target, wherever you've put them...
ln -s /Users/mspovich/research/M17SWex/mosaics/GLM_01424-058_mosaic_I1.fits $TARGET.mch1.fits
ln -s /Users/mspovich/research/M17SWex/mosaics/GLM_01424-058_mosaic_I2.fits $TARGET.mch2.fits
ln -s /Users/mspovich/research/M17SWex/mosaics/GLM_01424-058_mosaic_I3.fits $TARGET.mch3.fits
ln -s /Users/mspovich/research/M17SWex/mosaics/GLM_01424-058_mosaic_I4.fits $TARGET.mch4.fits
ln -s /Users/mspovich/research/M17SWex/mosaics/GLM_01424-058_mosaic_24p.fits $TARGET.mm24.fits
#OR make new ones from the original GLIMPSE/MIPSGAL mosaics, for
example (nessie):

% IDL>
  m = readfits('mosaic_enhanced_mips24_339.0.fits',hd)         
  xc = 7738 & yc = 1691 & dx = 2352 & dy = 2348                 
  x0 = xc-dx/2 & x1 = xc+dx/2
  y0 = yc-dy/2 & y1 = yc+dy/2
  hextract,im,hd,x0,x1,y0,y1
  writefits,'nessie_mm24.fits',im,hd   

#Get regionfile for ACIS field-of-view (fk5 DECIMAL DEGREES) and name it $TARGET.xfov.reg

####SED FITTING STARTS HERE####

mkdir sedfitter
cd sedfitter


=============================================================================
CREATE A TEXT FILE TO CONTAIN NOTES FROM YOUR SED FITTING TO THIS SPECIFIC TARGET. 
Save as sedfitting_notes.txt


=============================================================================
FIT ALL SOURCES WITH REDDENED STELLAR PHOTOSPHERES


## Prepare the fitter datafile from Pat's GLIMPSE+ data structure

  %/IDL>
  glimpseplusfile = '../glimpse+/Catalog_plus.fits.gz'
  datafile = 'data_glimpse+'
  vvv = 0 ;Set to 1 if using VVV catalog instead of UKIDSS
  select = 1  ;Set to 1 if coordinate box is desired, 0 for all sources set to 0
  subset = 0
  .run
    if select eq 1 then begin  ;eg M17 SWex field - COPY to sedfitting_notes.txt
      lmin = 13.694
      lmax = 14.759
      bmin = -1.0
      bmax = -0.262
      subset = [0.5*(lmin+lmax),0.5*(bmin+bmax),lmax-lmin,bmax-bmin]
    endif
  end
  glimpseplus2data,glimpseplusfile,select,datafile,vvv=vvv,subset=subset
;Copy printed output to sedfitting_notes.txt

  nwav = 10   ;Set to the number of filter bandpasses in DATAFILE
  if nwav gt 10 then bs = 2 else bs = 0
  twomass = 0   ; Set to 1 if your JHKs colors are on the 2MASS system
  if not twomass then mk = 1
  magfromdata,datafile,0+bs,j,nwav=nwav,mk=mk
  magfromdata,datafile,1+bs,h,nwav=nwav,mk=mk
  magfromdata,datafile,2+bs,k,nwav=nwav,mk=mk
  plot_nircc_rv,j,h,k,twomass=twomass

;Estimate the maximum reddening in Av by comparing the reddening locus of stars to the reddening ;vectors (marked at intervals of Av=5 mag). Note the value that you used.

; Filter out sources with colors consistent with stars, malmquist biases, and other general photometry problems

  maxav = 30. ;DEFAULT for IRDCs. Others should be LESS.
  data_good = datafile + 'keep'
  data_bad = datafile + 'toss'
  malmcullav,datafile,data_good,data_bad,maxav=maxav,nwav=nwav

###Fit all sources with stellar atmospheres.

  >>>
  python
from astropy import units as u
from sedfitter import fit
from sedfitter.extinction import Extinction

	#python seems to dislike the ~ for home directory character, using ABSOLUTE paths.
	# Define path to models
model_dir_kurucz = '/Users/mspovich/research/ysomodels/models_kurucz'

	# Read in extinction law
extinction = Extinction.from_file('/Users/mspovich/research/ysomodels/ex_law_gl.par', columns=[0, 3], wav_unit=u.micron, chi_unit=u.cm**2 / u.g)

	# Define filters and apertures (this example is for UKIDSS+2MASS+IRAC)
filters = ['UDSSJ', 'UDSSH', 'UDSSK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']
apertures = [2.0, 2.0, 2.0, 3., 3., 3., 3., 3., 3., 3.] * u.arcsec

        #VVV+2MASS+IRAC
filters = ['VVVZ','VVVY','VVVJ', 'VVVH', 'VVVK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']
apertures = [2.0,2.0,2.0, 2.0, 2.0, 3., 3., 3., 3., 3., 3., 3.] * u.arcsec

	# Fit! Note distance_range=[] is required but does nothing
	with these models! Use the Max Av you determined from the CCD;
	this can (maybe should) exceed the MAXAV used for MALMCULL AV
fit('data_glimpse+keep', filters, apertures, model_dir_kurucz, 'stellar_keep.fitinfo', n_data_min=4, extinction_law=extinction, distance_range=[1.5, 2.] * u.kpc, av_range=[0., 50.], output_format=('N',1))

##IMPORTANT! The maxav above is the one from the JHK CCD, and it exceeds the value used for MALMCULLAV. This is OK!

###Split the output into well-fit versus poorly-fit.
#After some testing, settled on cpd=9. to filter output of stellar
photosphere fits

>>>
from sedfitter import filter_output
filter_output('stellar_keep.fitinfo', cpd=9.) 

from sedfitter import write_parameters
write_parameters('stellar_keep.fitinfo_bad','pars_stellar_bad9.txt')

	 
% IDL>
	 fitinfobad9 = 'pars_stellar_bad9.txt'
	 data_parent='data_glimpse+keep'
	 fitinfo2data,fitinfobad9,data_parent,'data_glimpse+sb9'
	 ds9regfromdata,'data_glimpse+sb9','data_glimpse+sb9.reg',color='red'
	 data_activate,'data_glimpse+sb9','data_glimpse+ysoc',nwav=nwav

;Copy the output of the fitinfo2data command to sedfitting_notes.txt


###Let's check the spatial distributions of these sources...

  ds9 ../$TARGET.mch1.fits -region data_glimpse+sb9.reg &


###OPTIONAL Plot up some SEDs for sanity checks.

  >>>
from sedfitter import plot
plot('stellar_keep.fitinfo_bad', 'plots_sb9')

#It's OK to kill this with Control-C if you don't want to make thousands of plots!

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



   
