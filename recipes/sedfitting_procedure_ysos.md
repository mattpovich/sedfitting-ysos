# `sedfitting_procedure_ysos.txt`

**Authors: M. S. Povich & T. P. Robitaille**

## Description

This recipe guides you through the steps of fitting multiple SED models from [Robitaille (2017)](https://doi.org/10.1051/0004-6361/201425486) to combined *Spitzer*, 2MASS, and UKIDSS (or VVV) candidate young stellar objects (YSOs) with mid-infrared excess emission in a specified target field. The data analysis pipeline described here is inteded to be run *after* you have identified a set of candidate YSOs via the `find_ysoc_procedure` (included in this repo) or [`sedfitting_procedure_xpms`](https://github.com/mattpovich/sedfitting-phrds/blob/master/recipes/sedfitting_procedure_xpms.md) pipelines.

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

Create a directory (I call mine `ysomodels`) to hold all the model set subdirectories (note: if you've already downloaded model sets previously, for example in conjunction with one of my other recipes, it's a good idea to keep all the models in this same directory). Download the following 6 SED model sets from https://doi.org/10.5281/zenodo.166732 (there are 18 model sets total; these are the only ones used by this pipeline), expand each `*tar.gz` archive, and —> prepend the model set index number to each directory created:

* `sp--s-i.tar.gz —> 01_sp--s-i` 
* `sp--h-i.tar.gz —> 02_sp--h-i`
* `spu-smi.tar.gz —> 14_spu-smi`
* `spu-hmi.tar.gz —> 15_spu-hmi`
* `spubsmi.tar.gz —> 16_spubsmi`
* `spubhim.tar.gz —> 17_spubhim`

Copy the module `new_filt.py` from the `python` subdirectory of this library into the directory containing all of the above yso model subdirectories. The R17 models do not come pre-convolved with UKIDSS and VVV filter profiles. You can download the SVO filter profiles from http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php and colve multiple filters with R17 models in one shot using `new_filt.py,` which is a simple wrapper for `convolve_model_dir`. See https://sedfitter.readthedocs.io/en/stable/convolution.html for details.

Make sure you have the following libraries on your IDL path:

* [IDL astronomy users library](https://github.com/wlandsman/IDLAstro)
* [IDL coyote library](https://github.com/idl-coyote/coyote)
* [IDL library bundled with this recipe](https://github.com/mattpovich/sedfitting-phrds)


## TARGET-SPECIFIC SETUP

**>>>** `TARGET='m17swex'  # python EXAMPLE. **Choose your own target name.**`


**%** `setenv TARGET m17swex  # tcsh EXAMPLE. **Choose your own target name (must match above)**`

It is a good idea to name your working directory `$TARGET/sedfitter`, e.g. (in `tcsh`):

**%**

    mkdir -p $TARGET/sedfitter  #SKIP if directory exists from previous recipe!
    cd $TARGET/sedfitter
    
**Critical!** Copy the module `fit_multi.py` from the `python` subdirectory of this library into your `sedfitter` working directory.    
    
[Download](https://irsa.ipac.caltech.edu/data/SPITZER/MIPSGAL/) a **MIPSGAL 24 µm mosaic** covering your target field-of-view and name it `../$TARGET.mm24.fits`.

**Optional.** If you are interested in flagging YSOs found within a specfied "cropped" sub-region of your target catalog, make an `SAOImage ds9` regionfile containing one or more `polygon` regions and save it as `../$TARGET.xfov.reg` (required coordinate format is fk5 DECIMAL DEGREES). I personally use this feature to identify YSOs falling within the *Chandra*/ACIS field-of-view within a larger *Spitzer* mosaic.

## PREPARE SOURCE PHOTOMETRY FOR SED FITTING

There are three main ways you may have carved out a set of YSO candidates, so select the option below that matches your situation:

I. *Proceeding from `find_ysoc_procedure`.* YSO candidates were identified from a combined NIR/MIR catalog using my "blind" IR-excess selection criteria for highly-contaminated fields. Congratulations, you already made `data_glimpse+ysoc`, so proceed to **Flag UGOs** below.

II. *Proceeding from [`sedfitting_procedure_xpms`](https://github.com/mattpovich/sedfitting-phrds/blob/master/recipes/sedfitting_procedure_xpms.md).* YSO candidates were identified from a sample of IR counterparts to X-ray sources (or some other indicator of youth, e.g. H-alpha sources and/or a parallax/proper-motion-selected members of a young stellar cluster). In this case, your `$TARGET/sedfitter` directory should already contain two files, `xpms.fitinfo_bad` and `data_xir`. Use the following IDL command to create the new fitter data file:

**IDL>**

    parfile_bad = 'pars_xpms_bad.txt' 
    data_parent='data_xir'
    fitinfo2data,parfile_bad,data_parent,'data_glimpse+ysoc'

III. *Using your own custom set of YSO candidates.* If you already have a list of candidate YSOs that you know are mid-IR excess sources (from a previous color selection, for example) you can create your own `data_glimpse+ysoc` file conforming to the criteria described under **PREPARE SOURCE PHOTOMETRY FOR SED FITTING** in `find_ysoc_procedure` and skip the subsequent filtering steps in that recipe.
	 
## Flag "unresolved green objects" (UGOs)
Reset 4.5 µm fluxes to upper limits for suspected excess sources based on MIR colors (unresolved molecular outflow candidates; see [Povich et al. 2013](https://doi.org/10.1088/0067-0049/209/2/31)).

**IDL>**
        
	nwav = 10    ;set to 12 if ZY photometry is included in data_glimpse+ysoc
	ugoflag,'data_glimpse+ysoc','data_glimpse+ysoc_ugos',nwav=nwav,/regions		

## Perform aperture photometry on the MIPS 24 µm mosaic at each YSO position

I developed this pipeline back when no point-source catalog existed from MIPSGAL. Users of this pipeline may prefer to use the 24 µm catalog published by [Gutermuth & Heyer (2015)](https://doi.org/10.1088/0004-6256/149/2/64), which includes matches to GLIMPSE sources. *This is fine, but if you choose to go this route you are responsible for adding the 24 µm data to glimpse+ysoc_ugos_24um.*
In my experience, the majority of YSO candidates identifiable with IRAC lack firm detections at 24 µm, and the aperture photometry procedure below estimates useful upper limits for these sources.

**IDL>**

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

OPTIONAL. Visually review the extraction results on the 24 µm image. In the regionfile created, the photometry apertures are colored according to detections (creen circles), bad-pixels/saturation (black), or upper limits (red). Background annuli are colored white.

**%** `ds9 ../$TARGET.mm24.fits -region mipsmerger_ysoc_apertures.reg &`

If desired, iterate the above with modified values of aprad (aperture radius in arcsecond), anrad (annular background inner radius, in arcsec), or snlim (minimum S/N for a detection).


## Fit the 1-24 µm SEDs of YSO candidates with R17 YSO model sets.

**Preliminary step:** The two extinction laws that I use are packaged in the `ex_law` subdirectory of this repository. Place a copy of `ex_law_mc.par` into your `$TARGET/sedfitter` directory.

**Important!** After some experimentation on two test fields, I selected `cpd=4.` as the default value to separate well-fit from poorly-fit sources. But I recommend you experiment to settle on the optimal value, as this is likely to depend on the data set and selection criteria used.

*Note that the final command in the block below can take some time to run. On my MacBook Pro (2.7 GHz Intel Core i5, 16 GB RAM) the SED fitting for all included model sets takes ~1 hour per 500 sources fit.*

**>>>**

    from astropy import units as u
    from fit_multi import fit_multiyso
    extinction24_file = 'ex_law_mc.par'
    models_topdir = <path to your directory containing all the YSO model subdirectories> # Mine is called 'ysomodels'

    data = 'data_glimpse+ysoc_ugos_24um'

    #UKIDSS example. NOTE filter names MUST match filenames in the `ysomodels/*/convolved` subdirectories
    filters = ['UDSSJ', 'UDSSH', 'UDSSK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 'M1'] 
    apertures = [2., 2., 2., 3., 3., 3., 3., 3., 3., 3., 7.] * u.arcsec

    #VVV example
    filters = ['VVVZ', 'VVVY', 'VVVJ', 'VVVH', 'VVVK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 'M1'] 
    apertures = [2., 2., 2., 2., 2., 3., 3., 3., 3., 3., 3., 3., 7.] * u.arcsec

    fit_multiyso(data, filters, apertures, models_topdir,
		   extinction_file=extinction24_file,
		   n_data_min=3,  output_format=('F',6),
		   distance_range=[2.9, 3.3] * u.kpc,
		   av_range=[0.,100.], cpd=4.)


OPTIONAL: Plot up some good (and bad) SED fits. The example below compares one disk-only with one disk+envelope model set:

**>>>**
    
    from sedfitter import plot
    plot('01_sp--s-i.fitinfo_good', 'plots_01g4') 
    plot('01_sp--s-i.fitinfo_bad', 'plots_01b4') 
    plot('02_sp--h-i.fitinfo_good', 'plots_16g4') 
    plot('02_sp--h-i.fitinfo_bad', 'plots_16b4') 

Based on reviewing plots like these, you can revise the goodness-of-fit criterion `cpd` and iterate the above.

Once you are satisfied with your goodness-of-fit criterion, send all well-fit models into IDL save files in the `results_ysoc` directory.
This combines the best-fit models from all sets into SINGLE parameter tables for each source, choosing between model sets using the "Occam's razor" rule described by Robitaille (2017).

**IDL>**

    data = 'data_glimpse+ysoc_ugos_24um'
    spawn,'ls pars_*g04.txt',infiles 
    r17_pars2idl, infiles, 'results_ysoc', data_parent=data,filter='F6'
    
**Optional.** Free up disk space, for example (cpd=4):

**%**
    tar cvzf fitinfo_r17_cpd04.tgz *sp*.fitinfo*
    'rm' pars_*g04.txt  

## CREATE THE `$TARGET.ised.fits` SUMMARIZING THE SED FITTING RESULTS  

Included steps:
* Flag candidate starburst galaxies and AGN background contaminants (`galflag_phot.pro`)
* Flag candidate dusty AGB star contaminants
* Assign evolutionary stages to YSOs (`stageflag_r17.pro`)

NOTE: The (OPTIONAL) regionfile *xfov.reg must contain one (or more) non-overlapping polygons in the fk5 decimal degree coordinate system.

**IDL>**

     xfov = '../'+target+'.xfov.reg'	;OPTIONAL!
     create_ised_r17, xfov=xfov

**%**
    
    ln ised_r17.fits ../$TARGET.ised.fits
    ln ised_stageflag.reg ../$TARGET.stageflag.reg





   
