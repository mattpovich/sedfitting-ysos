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

Prepare the input fitter data file following the guidelines at https://sedfitter.readthedocs.io/en/stable/data.html. This pipeline assumes you have a total of `n=nwav=10` wavelengths/filters, in this order:

* UKIDSS: *JHK* (filters 1-3)
* 2MASS:  *JHKs* (filters 4-6)
* *Spitzer*/IRAC: [3.6], [4.5], [5.8], [8.0] (filters 7-10)

Other surveys could be substituted (e.g. *WISE* for *Spitzer*/IRAC is an obvious choice), but the ORDER that photometry points appear in the data file should not be changed. This pipeline also supports adding VVV/UKIDSS *YZ*, if available, in the filters 1 & 2 positions, in which case `n=12`.

### An automated option for producing the fitter data file 
I have included with this distribution an `IDL` procedure, `glimpseplus2data.pro`, that automatically creates the fitter data file from a `.fits` table containing the requisite photometry columns. Once you have downloaded both the GLIMPSE and UKIDSS/VVV catalogs for your desired target field, 
