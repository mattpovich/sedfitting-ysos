# `find_ysoc_procedure`

**Authors: M. S. Povich & T. P. Robitaille**

## Description
This recipe guides you through a series of color selections and SED fitting of reddened stellar photospheres to identify a sample of candidate young stellar objects (YSOs) with mid-IR excess emission in a specified region of GLIMPSE or a similar survey with *Spitzer* (or *WISE*). This general approach has been applied extensively to catalog YSOs in Galactic star-forming regions that suffer a high degree of contamination from field stars, in particular reddened background giants ([Povich et al. 2009](https://doi.org/10.1088/0004-637X/696/2/1278), [Povich & Whitney 2010](https://doi.org/10.1088/2041-8205/714/2/L285), Povich et al. [2011](https://doi.org/10.1088/0067-0049/194/1/14), [2013](https://doi.org/10.1088/0067-0049/209/2/31), [2016](https://doi.org/10.3847/0004-637X/825/2/125)). The end product is an ASCII file called `data_glimpse+ysoc` that can be used for subsequent SED fitting classification analysis of YSOs.

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

Download the [`models_kurucz` SED model set](ftp://ftp.astro.wisc.edu/outgoing/tom/model_packages/models_kurucz_05sep11.tgz). 

Make sure you have the following libraries on your IDL path:

* [IDL astronomy users library](https://github.com/wlandsman/IDLAstro)
* [IDL coyote library](https://github.com/idl-coyote/coyote)
* [IDL library bundled with this recipe](https://github.com/mattpovich/sedfitting-phrds)


## TARGET-SPECIFIC SETUP

**>>>** `TARGET='m17swex'  # python EXAMPLE. **Choose your own target name.**`

**%** `setenv TARGET m17swex  # tcsh EXAMPLE. **Choose your own target name (must match above)**`


It is a good idea to name your working directory `$TARGET/sedfitter`, e.g. (in `tcsh`):

**%**

    mkdir -p $TARGET/sedfitter
    cd $TARGET/sedfitter
    
**Recommended** Download *Spitzer*/IRAC mosaics of your target field and place them in the $TARGET directory, with filenames `../$TARGET.mch[1234].fits` for bands `[3.6], [4.5], [5.8], [8.0]`. GLIMPSE mosaic images can be obtained from the NASA/IPAC Infrared Science Archive at https://irsa.ipac.caltech.edu/data/SPITZER/GLIMPSE/.
    
## PREPARE SOURCE PHOTOMETRY FOR SED FITTING

You will need to assemble a catalog of point sources with available 1-8 µm broadband photometry. This pipeline assumes you are using the [GLIMPSE point-source catalog](http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=II/293), which includes 2MASS *JHKs* photometry. In addition, we recommend supplementing 2MASS/GLIMPSE by positionally matching to either the [UKIDSS](http://wsa.roe.ac.uk/index.html) (Northern hemisphere targets) or [VVV](https://vvvsurvey.org/data-releases/) (Southern hemisphere) point-source catalogs.

Prepare the input fitter data file (**`data_glimpse+`**) following the guidelines at https://sedfitter.readthedocs.io/en/stable/data.html. This pipeline assumes you have, *at minimum*, a total of `n=nwav=10` wavelengths/filters, in this order:

* UKIDSS: *JHK* (filters 1-3)
* 2MASS:  *JHKs* (filters 4-6)
* *Spitzer*/IRAC: [3.6], [4.5], [5.8], [8.0] (filters 7-10)

Other surveys could be substituted (e.g. *WISE* for *Spitzer*/IRAC is an obvious choice), but the ORDER that photometry points appear in the data file should not be changed. This pipeline also supports adding VVV/UKIDSS *YZ*, if available, in the filters 1 & 2 positions, in which case `n=12`.

I **strongly recommend** selecting the better of the duplicated 2MASS and UKIDSS/VVV *JHK(s)* fluxes and setting the flags for the *unused* fluxes to 0 initially. This will facilitate the filtering out of non-excess sources even in cases where there is disagreement in the near-IR fluxes due to variability or blending of sources. My own procedure for creating `data_glimpse+` (below) does exactly this, and the unused photometry can be re-activated later when fitting with YSO models.

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

This set of color cuts, described by [Povich et al. (2011)](https://doi.org/10.1088/0067-0049/194/1/14), gets rid of the vast majority of non-excess and "marginal"-excess sources, which greatly speeds up the next step and saves disk space.

**IDL>**

    datafile = 'data_glimpse+'
    maxav = 30. ;This is the MAXIMUM advisable value. Choose your own based on the maximum Av you determined above.
    data_good = datafile + 'keep'
    data_bad = datafile + 'toss'
    malmcullav,datafile,data_good,data_bad,maxav=maxav,nwav=nwav


## FIT REMAINING SOURCES WITH REDDENED STELLAR PHOTOSPHERES

**Preliminary step:** The two extinction laws that I use are packaged in the `ex_law` subdirectory of this repository. Place a copy of `ex_law_gl.par` into your `$TARGET/sedfitter` directory.

**>>>**

    from astropy import units as u
    from sedfitter import fit
    from sedfitter.extinction import Extinction
  
    # Read in extinction law
    extinction = Extinction.from_file('ex_law_gl.par', columns=[0, 3], 
        wav_unit=u.micron,    chi_unit=u.cm**2 / u.g)

  	 # Define path to models
    model_dir_kurucz = '<your path here> + /models_kurucz'

	   # Define filters and apertures
    filters = ['UDSSJ', 'UDSSH', 'UDSSK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']
    apertures = [2.0, 2.0, 2.0, 3., 3., 3., 3., 3., 3., 3.] * u.arcsec

The `models_kurucz` stellar atmosphere SEDs are aperture-independent, so the `apertures` variable doesn't matter as long as it has the correct length. Similarly, `distance_range=[min,max]` is required but does nothing since these models are scale-free.
Note these filter names are *specific* to the pre-convolved model SEDs and must match the filenames in the `models_pms/convolved` directory. You can [convolve the SED models with new broadband filters](https://sedfitter.readthedocs.io/en/stable/convolution.html) (I cannot guarantee that all UKIDSS and VVV filters are included in the `models_kurucz` distribution) using the `new_filt.py` module included with this package.

 **Important:** In the `fit()` call below make sure the `av_range=[]` reflects the maximum (and minimum if nonzero) extinction estimated from the *JHK* color-color diagram above. (It is fine—perhaps even preferable—if the maximum extinction exceeds the value of `maxav` used in the IDL> `malmcullav` call above.)
 
**>>>** 

    fit('data_glimpse+keep', filters, apertures, model_dir_kurucz, 'stellar_keep.fitinfo', n_data_min=4, 
    	extinction_law=extinction, distance_range=[1.5, 2.] * u.kpc, av_range=[0., 50.], output_format=('N',1))
		
Split the output into well-fit versus poorly-fit. After much testing, I recommend using `cpd=9.` as the cutoff for well-fit models from the `models_kurucz` set applied to this set of broadband photometry used in `glimpse_data+keep`.
		
**>>>**

    from sedfitter import filter_output
    filter_output('stellar_keep.fitinfo', cpd=9.) 

Write SED fit parameters to a text file (*badly-fit sources only*) and create `data_glimpse+ysoc` for subsequent fitting with YSO models.

**>>>**

    from sedfitter import write_parameters
    write_parameters('stellar_keep.fitinfo_bad','pars_stellar_bad9.txt')
    
**IDL>**
	
    fitinfobad9 = 'pars_stellar_bad9.txt'
    data_parent='data_glimpse+keep'
    fitinfo2data,fitinfobad9,data_parent,'data_glimpse+sb9'
    ds9regfromdata,'data_glimpse+sb9','data_glimpse+sb9.reg',color='red'
    data_activate,'data_glimpse+sb9','data_glimpse+ysoc',nwav=nwav
    

## RECOMMENDED QUALITY-CONTROL CHECKS

### Examine the spatial distributions of the candidate YSOs on an image of the target field.

**%** `ds9 ../$TARGET.mch1.fits -region data_glimpse+sb9.reg &`

Candidate YSOs should appear clustered in or near molecular clouds, IR dark clouds, on the rims of bubbles, and/or with known young stellar clusters in the field. If you judge that your YSO candidates dominated by a distributed, likely contaminating population, consider setting `cpd>9.` and rerunning `filter_output`.

### Plot up some badly-fit SEDs to assess why they went wrong.

**>>>**

    from sedfitter import plot
    plot('stellar_keep.fitinfo_bad', 'plots_sb9')

(Some target regions could contain thousands of YSO candidates. It's fine to kill the above command after awhile with Control-C if you don't want to make thousands of plots!)

If you judge that a substantial number of YSO candidates in these plots *should* have been well-fit by stellar photospheres, you can rerun `filter_output` above with `cpd>9.`
 

## NEXT STEPS:

Fit the sources in `data_glimpse+ysoc` with YSO models from Robitaille (2017) following the recipe in `sedfitting_procedure_ysos`.



   
 
