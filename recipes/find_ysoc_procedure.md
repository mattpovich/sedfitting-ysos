# `find_ysoc_procedure`

**Authors: M. S. Povich & T. P. Robitaille**

## Description
This recipe guides you through a series of color selections and SED fitting of reddened stellar photospheres to identify a sample of candidate young stellar objects (YSOs) with mid-IR excess emission in a specified region of the GLIMPSE or similar survey with *Spitzer* (or *WISE*). The final science product is an ASCII file called `data_glimpse+ysoc` that can be used for subsequent SED fitting classification analysis of YSOs.

## Version History

* Original sedfitting_procedure_ukidss-mir — 2011-08-24 by M. S. Povich 
* Updated to use  [`python sedfitter`](https://github.com/astrofrog/sedfitter) — 2018-07-23 by M. S. Povich 
* CURRENT (v1.0) — 2019-05-15 by M. S. Povich

 I recommend using two different terminal windows simultaneously, one in `bash` and the other in `tcsh`. The correct shell in which to execute each block of command is indicated as follows:
 
  * `python` command blocks are preceded by **>>>**
  * `tcsh/IDL` command blocks are preceded by **%/IDL>**


## INITIAL SETUP 

Get Tom Robitaille's sedfitter software and install it following the instructions
at http://sedfitter.readthedocs.io/en/stable/installation.html

Download the [`models_kurucz` SED model set](ftp://ftp.astro.wisc.edu/outgoing/tom/model_packages/models_kurucz_05sep11.tgz). 

Make sure you have the following libraries on your IDL path:

* [IDL astronomy users library](https://github.com/wlandsman/IDLAstro)
* [IDL coyote library](https://github.com/idl-coyote/coyote)
* [IDL library bundled with this recipe](https://github.com/mattpovich/sedfitting-phrds)

**>>>**

    TARGET='m17swex'  # python EXAMPLE. **Choose your own target name.**


**%**

    setenv TARGET m17swex  # tcsh EXAMPLE. **Choose your own target name (must match above)**


It is a good idea to name your working directory `$TARGET/sedfitter`, e.g. (in `tcsh`):

**%**

    mkdir -p $TARGET/sedfitter
    cd $TARGET/sedfitter
    
## PREPARE SOURCE PHOTOMETRY FOR SED FITTING

You will need to assemble a catalog of point sources with available 1-8 µm broadband photometry. This pipeline assumes you are using the [GLIMPSE point-source catalog](http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=II/293), which includes 2MASS *JHKs* photometry. In addition, we recommend supplementing 2MASS/GLIMPSE by positionally matching to either the [UKIDSS](http://wsa.roe.ac.uk/index.html) (Northern hemisphere targets) or [VVV](https://vvvsurvey.org/data-releases/) (Southern hemisphere) point-source catalogs.

Prepare the input fitter data file (**`data_glimpse+`**) following the guidelines at https://sedfitter.readthedocs.io/en/stable/data.html. This pipeline assumes you have, *at minimum*, a total of `n=nwav=10` wavelengths/filters, in this order:

* UKIDSS: *JHK* (filters 1-3)
* 2MASS:  *JHKs* (filters 4-6)
* *Spitzer*/IRAC: [3.6], [4.5], [5.8], [8.0] (filters 7-10)

Other surveys could be substituted (e.g. *WISE* for *Spitzer*/IRAC is an obvious choice), but the ORDER that photometry points appear in the data file should not be changed. This pipeline also supports adding VVV/UKIDSS *YZ*, if available, in the filters 1 & 2 positions, in which case `n=12`.

### An automated option for producing the fitter data file 
I have included with this distribution an `IDL` procedure, `glimpseplus2data.pro`, that automatically creates the fitter data file from a `.fits` table containing the requisite photometry columns. **WARNING: Currently this works only for GLIMPSE, not *WISE* or other MIR catalog formats.** 

Once you have downloaded both the GLIMPSE and UKIDSS/VVV catalog tables for your desired target field, perform a positional cross-match (the "matchstick" function in `TOPCAT` works well for this) between the catalogs. The new, matched catalog, `GLIMPSEc_plus.fits.gz`, **must** contain the following columns: `GLIMPSE` (source ID), `Glat` and `Glon` (source positions), 7 flux (`F*` and `e_F*`) columns, and 3 or 5 the UKIDSS/VVV magnitudes with uncertainties and error bit flags (`*AperMag3`, `*AperMag3err`, and `*ppeErrBits` columns).

**IDL>**

    glimpseplusfile = 'GLIMPSEc_plus.fits.gz'
    datafile = 'data_glimpse+'
    vvv = 0 ;Set to 1 if using VVV catalog instead of UKIDSS
    select = 0  ;Set to 1 if CROPPING to coordinate box is desired, and specify below
    subset = 0
    .run
      if select eq 1 then begin  ;eg M17 SWex field
        lmin = 13.694
        lmax = 14.759
        bmin = -1.0
        bmax = -0.262
        subset = [0.5*(lmin+lmax),0.5*(bmin+bmax),lmax-lmin,bmax-bmin]
      endif
    end
    glimpseplus2data,glimpseplusfile,select,datafile,vvv=vvv,subset=subset

## Estimate the maximum reddening to stars in the sample 

Plot and visually inspect the *J – H* vs. *H - K* color-color diagram of the source population:

**IDL>**

   	data = 'data_xir'
	   vvv = 0   ; Set to 1 if if using VVV catalog instead of UKIDSS
	   if not vvv then mk = 1
	   magfromdata,data,0,j,nwav=10,mk=mk
	   magfromdata,data,1,h,nwav=10,mk=mk
	   magfromdata,data,2,k,nwav=10,mk=mk
	   plot_nircc_rv,j,h,k,twomass=vvv
	
Estimate the maximum reddening in *Av* magnitudes  by comparing the locus
  of stars to the reddening vectors (marked at intervals of *Av* = 5
  mag). Record the maximum value you estimated. *OPTIONAL: estimate and record a minimum Av if this is not zero.*
  
## Filter out sources with colors consistent with stars, malmquist biases, and other photometry issues generally affecting longer wavelengths

This set of color cuts, described by [Povich et al. (2011)](https://ui.adsabs.harvard.edu/abs/2011ApJS..194...14P/abstract), gets rid of the vast majority of non-excess and "marginal"-excess sources, which greatly speeds up the next step and saves disk space.

**IDL>**

    datafile = 'data_glimpse+'
    maxav = 30. ;This is the MAXIMUM advisable value. Choose your own based on the maximum Av you determined above.
    data_good = datafile + 'keep'
    data_bad = datafile + 'toss'
    malmcullav,datafile,data_good,data_bad,maxav=maxav,nwav=nwav


## FIT REMAINING SOURCES WITH REDDENED STELLAR PHOTOSPHERES

**Preliminary step:** The two extinction laws that I use are packaged in the `ex_law` subdirectory of this repository. Place a copy of `ex_law_gl.par` into your `$TARGET/sedfitter` directory.

**>>>

    from astropy import units as u
    from sedfitter import fit
    from sedfitter.extinction import Extinction
  
    # Read in extinction law
    extinction = Extinction.from_file('/Users/mspovich/research/ysomodels/ex_law_gl.par', columns=[0, 3], 
        wav_unit=u.micron,    chi_unit=u.cm**2 / u.g)

  	 # Define path to models
    model_dir_kurucz = '<your path here> + /models_kurucz'

	   # Define filters and apertures
    filters = ['UDSSJ', 'UDSSH', 'UDSSK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']
    apertures = [2.0, 2.0, 2.0, 3., 3., 3., 3., 3., 3., 3.] * u.arcsec

The SED `models_pms` set is aperture-independent, so the `apertures` variable doesn't matter other than requiring the correct length. Similarly, `distance_range=[min,max]` is required but does nothing since these models are scale-free.
Note these filter names are *specific* to the pre-convolved model SEDs and must match the filenames in the `models_pms/convolved` directory. You can [convolve the SED models with new broadband filters](https://sedfitter.readthedocs.io/en/stable/convolution.html) (I cannot guarantee that all UKIDSS and VVV filters are included in the `models_kurucz` distribution) using the `new_filt.py` module included with this package.
 
	# Fit!  Use the Max Av you determined from the CCD;
	this can (maybe should) exceed the MAXAV used for MALMCULL AV
fit('data_glimpse+keep', filters, apertures, model_dir_kurucz, 'stellar_keep.fitinfo', n_data_min=4, extinction_law=extinction, distance_range=[1.5, 2.] * u.kpc, av_range=[0., 50.], output_format=('N',1))

   
 
