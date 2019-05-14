pro mipsfitapphot, mipsmosaic, data_in, data_out, photlist, $
                   zody=zody, aprad=aprad, anrad=anrad, apreg=apreg, $
                   satlev=satlev, meanexp=meanexp, snlim=snlim, $
                   min_frac_errs=min_frac_errs, fk5=fk5, nwav=nwav

;+
;NAME:
;     MIPSFITAPPHOT
;PURPOSE:
;Master routine for custom aperture photometry of a MIPS 24 um
;mosaic. 
;EXPLANATION:
;This routine is intended for a very specific purpose: extract
;fluxes or place upper limits for 24 um counterparts to candidate YSOs
;identified by IR excess in NIR+IRAC surveys. Note that the
;photometric accuracy is not expected to be better than 5%, and in
;some cases may be much, much worse.
;
;CALLING SEQUENCE:
;     mipsfitapphot, mipsmosaic, data_in, data_out, photlist,
;                   [ zody=zody, aprad=aprad, anrad=anrad, apreg=apreg,
;                   satlev=satlev, meanexp=meanexp, snlim=snlim, /fk5 ]
;
;INPUTS:
;     MIPSMOSAIC - (string) Path to .fits file containing MIPS 24 um
;                  mosaic.
;     DATA_IN    - (string) Path to fitter datafile (ASCII text)
;                  containing sources found in the field of
;                  MIPSMOSAIC. File must be in standard fitter data
;                  format (see manual by T. P. Robitaille) with n=7
;                  datapoints (2MASS+IRAC) OR n=11 (+MSX).  
;
;OUTPUTS:
;     DATA_OUT   - (string) Path to new fitter datafile (ASCII text)
;                  to write. This will be DATA_IN with 24 um
;                  photometry information added.
;     PHOTLIST   - (string) Path to file containing ASCII table with
;                  24 um photometry results.
;
;KEYWORD PARAMETERS:
;     ZODY       - (float) Value of zodaical background intensity
;                  subtracted from MIPSMOSAIC, in MJy/sr. DEFAULT =
;                  60 MJy/sr, which is near the max Zody at 24 um;
;                  Kelsall et al. 1998, ApJ, 508, 44.
;     APRAD      - (float) Radius of extraction aperture
;                  (arcsec). DEFAULT = 7". NOTE: MIPSFITAPPHOT
;                  will use APRAD to choose the nearest of the following
;                  values: 3.5", 7", 13", or 20".
;     ANRAD      - (float) Inner radius of background annulus
;                  (arcsec). DEFAULT = 7". NOTE: ANRAD must be
;                  >= APRAD. MIPSFITAPPHOT will use ANRAD to choose
;                  the nearest of the following annulus sizes:
;                  6-8", 7-13", 20-32", or 40-50".
;     APREG      - (string) Path to ds9 regionfile containing
;                  apertures and background annuli in image
;                  coordinates.
;     SATLEV     - (float) Intensity value for which MIPSMOSAIC enters
;                  non-linear saturation regime. DEFAULT = 1700 MJy/sr.
;     MEANEXP    - (float) Mean exposure time per pixel in MIPSMOSAIC
;                  (s). DEFAULT = 34.5 s, appropriate for M17 in MIPSGAL.
;     SNLIM      - (float) S/N limit for detection, below which
;                  extracted flux is used as an upper limit only. 
;                  DEFAULT = 5.
;     MIN_FRAC_ERRS - Minimum fractional error to set DATA_OUT. DEFAULT=0.05
;     COORDSEL   - (integer) Set if DATA_IN and MIPSMOSAIC use different coordinate systems.
;                  Follows EULER convention: 1 = FK5 (DATA) to Galactic (MOSAIC),
;                  2 = Galactic (DATA) to FK5 (MOSAIC)
;     NWAV=      - (integer) Number of photometry bandpasses in
;                  data_in. DEFAULT = 10
;
;PROCEDURES CALLED:
;     READFITS, ADXY, EXTAST, APER_SAT, APREGIONS, MIPS24APCORR
;
;REVISION HISTORY:
;     Written January 2010       M. S. Povich
;        Modified saturation flag to suppress extraction if average of
;aperture values exceeds SATLEV.  MSP 2010/2/23
;
;     PRODUCTION VERSION (v1.0)      6 May 2019   MSP
;        - Added NWAV keyword and updated code to read column-separated
;          DATA_IN file.
;        - Updated replaced FK5 keyword with COORDSEL, expanded functionality
;-

if n_params() lt 3 then begin
  print, 'syntax -  mipsfitapphot, mipsmosaic, data_in, data_out, photlist,'
  print, '              zody=60., aprad=7., anrad=7., apreg=,'
  print, '               satlev=1700., meanexp=34.5, snlim=5., '
  print, '               min_frac_errs=0.05, coordsel='
  return
endif

;Set Keyword defaults
  if not keyword_set(nwav) then nwav = 10
  if not keyword_set(snlim) then snlim = 5.
  if not keyword_set(min_frac_errs) then min_frac_errs = 0.15
  if not keyword_set(zody) then zody = 60.
  if not keyword_set(aprad) then aprad = 7.
  if not keyword_set(anrad) then anrad = 7. else begin
     if anrad lt aprad then anrad = aprad
  endelse 
  if not keyword_set(satlev) then satlev = 1700.
  maxint = satlev + zody ;Set saturation value to give APER_SAT
  if not keyword_set(meanexp) then meanexp = 34.5 ;This could be determined from FITS header except there's no guarantee that MIPSMOSAIC had a full header..
  fconv = 0.0454   ; (DN/s)/(MJy/sr) Typical for MIPS24 um, +/- 4%.
  gain = 5.0       ; e-/DN
  rn = 40          ; Readnoise e- (Rieke et al. 2004, SPIE)
  

;Establish extraction aperture and background annulus. Choose from a
;predtermined set so that aperture corrections from the MIPS Data
;Handbook v. 3.3.1 can be used (see MIPS24APCORR).

  aprs = [ 3.5, 7.0, 13.0, 20.0 ]
  anrs = [ [6.,8.], [7.,13.], [20.,32.], [40.,50.] ]

  case 1 of
     aprad lt 6. : aprsel = 0
     aprad ge 6. and aprad lt 10. : aprsel = 1
     aprad ge 10. and aprad lt 17. : aprsel = 2
     aprad ge 17. : aprsel = 3
  endcase
  case 1 of
     anrad lt 6. : anrsel = 0
     anrad ge 6. and anrad lt 10. : anrsel = 1
     anrad ge 10. and anrad lt 35. : anrsel = 2
     anrad ge 35. : anrsel = 3
  endcase

;Read fitter datafile
  nlines = file_lines(data_in)
  lines = strarr(nlines)
  openr,u,data_in,/get_lun
  line = ''
  for i=0L,nlines-1 do begin
     readf,u,line
     lines[i] = line
  endfor 
  close,u
  free_lun,u
  
 ;Find the start, end position for each COLUMN
  ncol = 3*(nwav+1)
  colstart = intarr(ncol)
  colend = intarr(ncol)
  pos = 0
  for i=0,ncol-1 do begin
     while strpos(line,' ',pos) eq pos do pos++
     colstart[i] = pos
     while strpos(line,' ',pos) ne pos do pos++
     if pos ne 0 then colend[i] = pos else $
        colend[i] = strlen(line)
  endfor
  
;These could be RA, Dec, but I always use Galactic.
  l = double(strmid(lines,colstart[1],colend[1]-colstart[1]))
  b = double(strmid(lines,colstart[2],colend[2]-colstart[2]))
    
;Read in image mosaic
  image = readfits(mipsmosaic, head)

  image = image + zody  ;Add Zodaical background to mosaic (for Poisson noise calculations.)
  
;Find locations of sources from DATA_IN in x,y coordinates for IMAGE
  if keyword_set(coordsel) then begin ;Convert Galactic to FK5 coords if necessary
     if coordsel ne 1 then coordsel = 2 ; Idiot-proofing!
     euler, l, b, xcoord, ycoord, coordsel
  endif else begin      
     xcoord = l
     ycoord = b
  endelse 
  adxy, head, xcoord, ycoord, x, y

;Get mosaic plate scale from FITS header
  extast, head, astr, nopars

  ;Plate scale in ARCSEC/pix
  case 1 of 
     nopars eq 1 : begin
        xdelt = astr.cdelt[0]*3600.
        ydelt = astr.cdelt[1]*3600.
     end
     nopars eq 2 : begin
        xdelt = astr.cd[0,0]*3600.
        ydelt = astr.cd[1,1]*3600.
     end 
     else : begin
        print, 'FITS header for '+mipsmosaic+' does not contain standard astrometry keywords! RETURNING!'
        return
     endelse 
  endcase 

;Convert aperture radii into image (pixel) coordinates
  
  apradpix = aprs[aprsel]/abs(ydelt)     ;Extraction radius
  anradpix = anrs[*,anrsel]/abs(ydelt)   ;Background annulus radii [inner,outer]

;Conversion from MJy/sr to mJy/pix
  k = 0.135*(ydelt/2.4)^2

;Compute phpadu values for APER_SAT.
  int2counts = fconv*gain*meanexp

;Perform aperture photometry
  aper_sat, image, x, y, f, df, sky, dsky, int2counts, apradpix, anradpix $
            , [0,0], satflags, satval=maxint, print=photlist, /silent, /flux, $
            readnoise = rn

  f = f*k & df = df*k & sky = sky*k & dsky = dsky*k
  
;Flag each source as unextracted (0), detected (1), lower limit (2) or
;upper limit (3)
  flags = intarr(nlines)
  ind_det = where(finite(f),n_ex,complement=ind_noex,ncomplement=n_noex)
  if n_ex eq 0 then message,'No sources extracted!!'
  if n_noex ne 0 then begin
     f[ind_noex] = 0
     df[ind_noex] = 0
  endif 
  flags[ind_det] = 1
  ind_sat = where(satflags eq 1,n_sat)
  if n_sat ne 0 then begin 
     flags[ind_sat] = 2
     df[ind_sat] = 1
  endif 
  ind_ul = where(f/df le snlim,n_ul)
  if n_ul ne 0 then begin
     flags[ind_ul] = 3
     ;Treat positive and negative extractions differently
     ind_pos = where(f[ind_ul] gt 0,n_pos,complement=ind_neg,ncomplement=n_neg)  
     if n_pos ne 0 then f[ind_ul[ind_pos]] = f[ind_ul[ind_pos]] + snlim*df[ind_ul[ind_pos]]
     if n_neg ne 0 then f[ind_ul[ind_neg]] = snlim*df[ind_ul[ind_neg]]
     df[ind_ul] = 1.
  endif
  ind_allsat = where(satflags eq 2,n_asat)
  if n_asat ne 0 then begin
     flags[ind_allsat] = 0
     f[ind_allsat] = 0
     df[ind_allsat] = 0
  endif 

;Make ds9 regionfile of extraction apertures, if desired.

  if keyword_set(apreg) then begin
     xreg = x + 1
     yreg = y + 1
     apregions, apreg, xreg, yreg, apradpix, anradpix[0], anradpix[1], colors=flags
  endif 

  ind_fg = where(flags eq 1,nfg)
  if nfg eq 0 then message,'No good detections!'

;Perform aperture correction
  mips24apcorr, aprsel, anrsel, apcorr
  f = f*apcorr 
  df[ind_fg] = df[ind_fg]*apcorr

;Reset errors to minimum
  ind_reset = where(df[ind_fg]/f[ind_fg] lt min_frac_errs,nreset)
  if nreset ne 0 then df[ind_fg[ind_reset]] = min_frac_errs*f[ind_fg[ind_reset]]
  print,'Reset errors to ',strtrim(100*min_frac_errs,2),'% for ',strtrim(nreset,2),'/',strtrim(nfg,2),' detected sources.'

;Print some stats
  print,'Total sources:     ',nlines
  print,'Extracted sources: ',n_elements(where(flags ne 0))
  print,'Good fluxes:       ',n_elements(where(flags eq 1))
  print,'Upper limits:      ',n_elements(where(flags eq 3))

;Add 24 um fluxes to lines from DATA_IN
;  linelen = strlen(lines[0])
  ;; case 1 of
  ;;    linelen eq 232 : cut = 64   ;2MASS + IRAC
  ;;    linelen eq 310 : cut = 70   ;UKIDSS + 2MASS + IRAC
  ;;    linelen ge 334 : cut = 72   ;2MASS + IRAC + MSX
  ;;    else : begin
  ;;       print,'DATA_IN is not in correct format. Cannot print DATA_OUT. Returning.'
  ;;       return
  ;;    end 
  ;; endcase 

  cut = colend[(nwav + 3) - 1]
  
  lines1 = strmid(lines,0,cut)
  lines2 = strmid(lines,cut)
  openw,v,data_out,/get_lun
  for i=0L,nlines-1 do printf,v,format='(A,I2,A,2(2X,E10.4))', $
                              lines1[i],flags[i],lines2[i],f[i],df[i]
  close,v
  free_lun,v


end


pro mips24apcorr, aprsel, anrsel, apcorr

;+
;NAME:
;    MIPS24APCORR
;PURPOSE:
;Lookup aperture correction for MIPS 24 um photometry given a
;combination of extraction aperture and background annulus. Called by
;MIPSFITAPPHOT.
;
;
;CALLING SEQUENCE:
;     mips24apcorr, aprsel, anrsel, apcorr
;
;INPUTS:
;     APRSEL -- (integer) Selects the aperture radius used as follows:
;        0, 1, 2, 3 ==> 3.5", 7.0", 13.0", 20.0"
;     ANRSEL -- (integer) Selects the background annulus used as follows:
;        0, 1, 2, 3 ==> 6-8", 7-13", 20-32", 40-50"
; 
;OUTPUT:
;     APCORR -- (float) Aperture correction from the lookup table.
;
;REVISION HISTORY
;     Written January 2010        M. S. Povich
;
;NOTE: Lookup table values are taken from the MIPS Data Handbook v.3.3.1
;(http://ssc.spitzer.caltech.edu/mips/dh/mipsdatahandbook3.3.1.pdf)
;
;-

if n_params() lt 2 then begin
   print, 'syntax - mips24apcorr, aprsel, anrsel, apcorr'
   return
endif

;Define aperture correction lookup table
  apcorrtable = $
     [ [2.78, 2.80, 2.57, 2.56], $
       [-99., 2.05, 1.61, 1.61], $
       [-99., -99., 1.17, 1.16], $
       [-99., -99., -99., 1.08]]

  apcorr = apcorrtable[anrsel,aprsel]

  print,'MIPS24APCORR: Aperture correction = ',strtrim(apcorr,2),'.'

end


pro apregions, regfile, x, y, rap, ranin, ranout, colors=colors

;+
;NAME:
;    APREGIONS
;PURPOSE:
;Generate a ds9 region file in IMAGE coordinates showing photometric
;extraction apertures. Main purpose is for checking aperture
;placement. Called by MIPSFITAPPHOT.
;
;CALLING SEQUENCE:
;     apregions, regfile, x, y, rap, ranin, ranout, colors=colors
;
;INPUTS:
;     REGFILE    -- (string) Path to ds9 regionfile to write.
;     X,Y          -- (float[nsrc]) x,y coordinates of aperture centers.
;     RAP          -- (float[nsrc]) extraction aperture radii.
;     RANIN,RANOUT -- (float[nsrc]) inner, outer radii for background
;                     annuli.
;
;KEYWORD:
;     COLORS=      -- (integer[nsrc]) Vector containing keys to
;                     color-coding for extraction apertures. DEFAULT:
;                     All apertures colored green.
;
;REVISION HISTORY
;     Written January 2010       M. S. Povich
;
;-

if n_params() lt 6 then begin
   print,'syntax - apregions, regfile, x, y, rap, ranin, ranout, [colors=]'
   return
endif

  nsrc = n_elements(x)

  if not keyword_set(colors) then colors = replicate(3,nsrc)
  colorkey = [ 'black', 'green', 'blue', 'red', 'white', 'yellow' ] 

  openw,10,regfile
  for i=0L,nsrc-1 do begin
     printf,10,'image; circle(' + strtrim(x[i],2) + ',' + strtrim(y[i],2) $
            + ',' + strtrim(rap,2) + ') # color=' + colorkey[colors[i]]
     printf,10,'image; annulus(' + strtrim(x[i],2) + ',' + strtrim(y[i],2) $
            + ',' + strtrim(ranin,2) + ',' + strtrim(ranout,2) $
            + ') # color=white'
  endfor 
  close,10

  print,'APREGIONS: Wrote ds9 regionfile ',regfile

end

pro aper_sat,image,xc,yc,mags,errap,sky,skyerr,phpadu,apr,skyrad, $
             badpix, satsrc, $
       SETSKYVAL = setskyval,PRINT = print, SILENT = silent, FLUX=flux, $
       EXACT = exact, Nan = nan, READNOISE = readnoise, MEANBACK = meanback, $
       CLIPSIG=clipsig, MAXITER=maxiter,CONVERGE_NUM=converge_num $
             , satval=satval
;+
; NAME:
;      APER_SAT
; PURPOSE:
;      Compute concentric aperture photometry (adapted from DAOPHOT) 
; EXPLANATION:
;     APER can compute photometry in several user-specified aperture radii.  
;     A separate sky value is computed for each source using specified inner 
;     and outer sky radii.   
;
; CALLING SEQUENCE:
;     APER, image, xc, yc, [ mags, errap, sky, skyerr, phpadu, apr, skyrad, 
;                       badpix, /NAN, /EXACT, /FLUX, PRINT = , /SILENT, 
;                       /MEANBACK, SETSKYVAL = ]
; INPUTS:
;     IMAGE -  input image array
;     XC     - vector of x coordinates. 
;     YC     - vector of y coordinates
;
; OPTIONAL INPUTS:
;     PHPADU - Photons per Analog Digital Units, numeric scalar.  Converts
;               the data numbers in IMAGE to photon units.  (APER assumes
;               Poisson statistics.)  
;     APR    - Vector of up to 12 REAL photometry aperture radii.
;     SKYRAD - Two element vector giving the inner and outer radii
;               to be used for the sky annulus.   Ignored if the SETSKYVAL
;              keyword is set.
;     BADPIX - Two element vector giving the minimum and maximum value
;               of a good pixel.   If badpix is not supplied or if BADPIX[0] is
;               equal to BADPIX[1] then it is assumed that there are no bad
;               pixels.     Note that fluxes will not be computed for any star
;               with a bad pixel within the aperture area, but that bad pixels
;               will be simply ignored for the sky computation.    The BADPIX
;               parameter is ignored if the /NAN keyword is set.
;
; OPTIONAL KEYWORD INPUTS:
;     CLIPSIG - if /MEANBACK is set, then this is the number of sigma at which 
;             to clip the background.  Default=3
;     CONVERGE_NUM:  if /MEANBACK is set then if the proportion of 
;           rejected pixels is less than this fraction, the iterations stop.  
;           Default=0.02, i.e., iteration stops if fewer than 2% of pixels 
;           excluded.
;     /EXACT -  By default, APER counts subpixels, but uses a polygon 
;             approximation for the intersection of a circular aperture with
;             a square pixel (and normalizes the total area of the sum of the
;             pixels to exactly match the circular area).   If the /EXACT 
;             keyword, then the intersection of the circular aperture with a
;             square pixel is computed exactly.    The /EXACT keyword is much
;             slower and is only needed when small (~2 pixels) apertures are
;             used with very undersampled data.    
;     /FLUX - By default, APER uses a magnitude system where a magnitude of
;               25 corresponds to 1 flux unit.   If set, then APER will keep
;              results in flux units instead of magnitudes.    
;     MAXITER if /MEANBACK is set then this is the ceiling on number of 
;             clipping iterations of the background.  Default=5
;     /MEANBACK - if set, then the background is computed using the 3 sigma 
;             clipped mean (using meanclip.pro) rather than using the mode 
;             computed with mmm.pro.    This keyword is useful for the Poisson 
;             count regime or where contamination is known  to be minimal.
;     /NAN  - If set then APER will check for NAN values in the image.   /NAN
;             takes precedence over the BADPIX parameter.   Note that fluxes 
;             will not be computed for any star with a NAN pixel within the 
;             aperture area, but that NAN pixels will be simply ignored for 
;             the sky computation.
;     PRINT - if set and non-zero then APER will also write its results to
;               a file aper.prt.   One can specify the output file name by
;               setting PRINT = 'filename'.
;     READNOISE - Scalar giving the read noise (or minimum noise for any
;              pixel.   This value is passed to the procedure mmm.pro when
;              computing the sky, and is only need for images where
;              the noise is low, and pixel values are quantized.   
;     /SILENT -  If supplied and non-zero then no output is displayed to the
;               terminal.
;     SETSKYVAL - Use this keyword to force the sky to a specified value 
;               rather than have APER compute a sky value.    SETSKYVAL 
;               can either be a scalar specifying the sky value to use for 
;               all sources, or a 3 element vector specifying the sky value, 
;               the sigma of the sky value, and the number of elements used 
;               to compute a sky value.   The 3 element form of SETSKYVAL
;               is needed for accurate error budgeting.
;     SATVAL  - Specify maximum IMAGE value above which a source is
;               saturated and the extracted FLUX will be considered a
;               lower limit.
;
; OUTPUTS:
;     MAGS   -  NAPER by NSTAR array giving the magnitude for each star in
;               each aperture.  (NAPER is the number of apertures, and NSTAR
;               is the number of stars).   If the /FLUX keyword is not set, then
;               a flux of 1 digital unit is assigned a zero point magnitude of 
;               25.
;     ERRAP  -  NAPER by NSTAR array giving error for each star.  If a 
;               magnitude could not be determined then  ERRAP = 9.99 (if in 
;                magnitudes) or ERRAP = !VALUES.F_NAN (if /FLUX is set).
;     SKY  -    NSTAR element vector giving sky value for each star in 
;               flux units
;     SKYERR -  NSTAR element vector giving error in sky values
;     SATSRC -  NSTAR element vector flagging saturated sources (used
;              only if SATVAL keyword is defined).
;
; EXAMPLE:
;       Determine the flux and error for photometry radii of 3 and 5 pixels
;       surrounding the position 234.2,344.3 on an image array, im.   Compute
;       the partial pixel area exactly.    Assume that the flux units are in
;       Poisson counts, so that PHPADU = 1, and the sky value is already known
;       to be 1.3, and that the range [-32767,80000] for bad low and bad high
;       pixels
;      
;
;       IDL> aper, im, 234.2, 344.3, flux, eflux, sky,skyerr, 1, [3,5], -1, $
;            [-32767,80000],/exact, /flux, setsky = 1.3
;       
; PROCEDURES USED:
;       GETOPT, MMM, PIXWT(), STRN(), STRNUMBER()
; NOTES:
;       Reasons that a valid magnitude cannot be computed include the following:
;      (1) Star position is too close (within 0.5 pixels) to edge of the frame
;      (2) Less than 20 valid pixels available for computing sky
;      (3) Modal value of sky could not be computed by the procedure MMM
;      (4) *Any* pixel within the aperture radius is a "bad" pixel
;      (5) The total computed flux is negative
;
;
;       For the case where the source is fainter than the background, APER will
;       return negative fluxes if /FLUX is set, but will otherwise give 
;       invalid data (since negative fluxes can't be converted to magnitudes) 
; 
;       APER was modified in June 2000 in two ways: (1) the /EXACT keyword was
;       added (2) the approximation of the intersection of a circular aperture
;       with square pixels was improved (i.e. when /EXACT is not used) 
; REVISON HISTORY:
;       Adapted to IDL from DAOPHOT June, 1989   B. Pfarr, STX
;       FLUX keyword added                       J. E. Hollis, February, 1996
;       SETSKYVAL keyword, increase maxsky       W. Landsman, May 1997
;       Work for more than 32767 stars           W. Landsman, August 1997
;       Don't abort for insufficient sky pixels  W. Landsman  May 2000
;       Added /EXACT keyword                     W. Landsman  June 2000 
;       Allow SETSKYVAL = 0                      W. Landsman  December 2000 
;       Set BADPIX[0] = BADPIX[1] to ignore bad pixels W. L.  January 2001     
;       Fix chk_badpixel problem introduced Jan 01 C. Ishida/W.L. February 2001
;       Set bad fluxes and error to NAN if /FLUX is set  W. Landsman Oct. 2001 
;       Remove restrictions on maximum sky radius W. Landsman  July 2003
;       Added /NAN keyword  W. Landsman November 2004
;       Set badflux=0 if neither /NAN nor badpix is set  M. Perrin December 2004
;       Added READNOISE keyword   W. Landsman January 2005
;       Added MEANBACK keyword   W. Landsman October 2005
;       Correct typo when /EXACT and multiple apertures used.  W.L. Dec 2005
;       Remove VMS-specific code W.L. Sep 2006
;       Add additional keywords if /MEANBACK is set W.L  Nov 2006
;       Allow negative fluxes if /FLUX is set  W.L.  Mar 2008
;       Previous update would crash if first star was out of range  W.L. Mar 2008
;
;       Added SATLEV keyword, SATSRC output, and changed name to APER_SAT
;       M. S. Povich Jan 2010
;       Added SATLEV=2 option to flag sources where the mean value in
;       the aperture is saturated.
;
;-
 COMPILE_OPT IDL2
 On_error,2
;             Set parameter limits
 minsky = 20   ;Smallest number of pixels from which the sky may be determined
 maxsky = 10000         ;Maximum number of pixels allowed in the sky annulus.
;                                
if N_params() LT 3 then begin    ;Enough parameters supplied?
  print, $
  'Syntax - APER_SAT, image, xc, yc, [ mags, errap, sky, skyerr, phpadu, apr, '
  print,'             skyrad, badpix, /EXACT, /FLUX, SETSKYVAL = ,PRINT=, '
  print,'             /SILENT, /NAN, SATVAL= ]'
  return
endif 

 s = size(image)
 if ( s[0] NE 2 ) then message, $
       'ERROR - Image array (first parameter) must be 2 dimensional'
 ncol = s[1] & nrow = s[2]           ;Number of columns and rows in image array

  silent = keyword_set(SILENT)

 if not keyword_set(nan) then begin
 if (N_elements(badpix) NE 2) then begin ;Bad pixel values supplied
GET_BADPIX:  
   ans = ''
   print,'Enter low and high bad pixel values, [RETURN] for defaults'
   read,'Low and high bad pixel values [none]: ',ans
   if ans EQ  '' then badpix = [0,0] else begin
   badpix = getopt(ans,'F')
   if ( N_elements(badpix) NE 2 ) then begin
        message,'Expecting 2 scalar values',/continue
        goto,GET_BADPIX
   endif
   endelse
 endif 

 chk_badpix = badpix[0] LT badpix[1]     ;Ignore bad pixel checks?
 endif

 if ( N_elements(apr) LT 1 ) then begin              ;Read in aperture sizes?
   apr = fltarr(10)
   read, 'Enter first aperture radius: ',ap
   apr[0] = ap
   ap = 'aper'
   for i = 1,9 do begin                                                   
GETAP: 
      read,'Enter another aperture radius, [RETURN to terminate]: ',ap
      if ap EQ '' then goto,DONE  
      result = strnumber(ap,val)
      if result EQ 1 then apr[i] = val else goto, GETAP   
   endfor
DONE: 
   apr = apr[0:i-1]
 endif


 if N_elements(SETSKYVAL) GT 0 then begin
     if N_elements( SETSKYVAL ) EQ 1 then setskyval = [setskyval,0.,1.]
     if N_elements( SETSKYVAL ) NE 3 then message, $
        'ERROR - Keyword SETSKYVAL must contain 1 or 3 elements'
     skyrad = [ 0., max(apr) + 1]
 endif

;Get radii of sky annulii

 if N_elements(skyrad) NE 2 then begin
   skyrad = fltarr(2)
   read,'Enter inner and outer sky radius (pixel units): ',skyrad
 endif else skyrad = float(skyrad)

 if ( N_elements(phpadu) LT 1 ) then $ 
   read,'Enter scale factor in Photons per Analog per Digital Unit: ',phpadu

 Naper = N_elements( apr )                        ;Number of apertures
 Nstars = min([ N_elements(xc), N_elements(yc) ])  ;Number of stars to measure

 ms = strarr( Naper )       ;String array to display mag for each aperture
 if keyword_set(flux) then $
          fmt = '(F8.1,1x,A,F7.1)' else $           ;Flux format
          fmt = '(F9.3,A,F5.3)'                  ;Magnitude format
 fmt2 = '(I5,2F8.2,F7.2,3A,3(/,28x,4A,:))'       ;Screen format
 fmt3 = '(I4,5F8.2,6A,2(/,44x,9A,:))'            ;Print format

 mags = fltarr( Naper, Nstars) & errap = mags           ;Declare arrays
 sky = fltarr( Nstars )        & skyerr = sky     
 area = !PI*apr*apr                 ;Area of each aperture

 if keyword_set(EXACT) then begin
      bigrad = apr + 0.5
      smallrad = apr/sqrt(2) - 0.5 
 endif
     

 if N_elements(SETSKYVAL) EQ 0 then begin

     rinsq =  (skyrad[0]> 0.)^2 
     routsq = skyrad[1]^2
 endif 

 if keyword_set(PRINT) then begin      ;Open output file and write header info?
   if size(PRINT,/TNAME) NE 'STRING'  then file = 'aper.prt' $
                                   else file = print
   message,'Results will be written to a file ' + file,/INF
   openw,lun,file,/GET_LUN
   printf,lun,'Program: APER: '+ systime(), '   User: ', $
      getenv('USER'),'  Host: ',getenv('HOST')
   for j = 0, Naper-1 do printf,lun, $
               format='(a,i2,a,f4.1)','Radius of aperture ',j,' = ',apr[j]
   if N_elements(SETSKYVAL) EQ 0  then begin
   printf,lun,f='(/a,f4.1)','Inner radius for sky annulus = ',skyrad[0]
   printf,lun,f='(a,f4.1)', 'Outer radius for sky annulus = ',skyrad[1]
   endif else printf,lun,'Sky values fixed at ', strtrim(setskyval[0],2)
   if keyword_set(FLUX) then begin
       printf,lun,f='(/a)', $
           'Star   X       Y        Sky   SkySig    SkySkw   Fluxes'
      endif else printf,lun,f='(/a)', $
           'Star   X       Y        Sky   SkySig    SkySkw   Magnitudes'
 endif
 print = keyword_set(PRINT)

;         Print header
 if not SILENT then begin
    if (KEYWORD_SET(FLUX)) then begin
       print, format="(/1X,'Star',5X,'X',7X,'Y',6X,'Sky',8X,'Fluxes')"
    endif else print, $ 
       format="(/1X,'Star',5X,'X',7X,'Y',6X,'Sky',8X,'Magnitudes')" 
 endif

;  Compute the limits of the submatrix.   Do all stars in vector notation.

 lx = fix(xc-skyrad[1]) > 0           ;Lower limit X direction
 ux = fix(xc+skyrad[1]) < (ncol-1)    ;Upper limit X direction
 nx = ux-lx+1                         ;Number of pixels X direction
 ly = fix(yc-skyrad[1]) > 0           ;Lower limit Y direction
 uy = fix(yc+skyrad[1]) < (nrow-1);   ;Upper limit Y direction
 ny = uy-ly +1                        ;Number of pixels Y direction
 dx = xc-lx                         ;X coordinate of star's centroid in subarray
 dy = yc-ly                         ;Y coordinate of star's centroid in subarray

 edge = (dx-0.5) < (nx+0.5-dx) < (dy-0.5) < (ny+0.5-dy) ;Closest edge to array
 badstar = ((xc LT 0.5) or (xc GT ncol-1.5) $  ;Stars too close to the edge
        or (yc LT 0.5) or (yc GT nrow-1.5))
;
 badindex = where( badstar, Nbad)              ;Any stars outside image
 if ( Nbad GT 0 ) then message, /INF, $
      'WARNING - ' + strn(nbad) + ' star positions outside image'
      if keyword_set(flux) then begin 
          badval = !VALUES.F_NAN
	  baderr = badval
      endif else begin 
          badval = 99.999
	  baderr = 9.999
      endelse	  	  
 
;Set up SATSRC array
      satsrc = intarr(nstars)

 for i = 0L, Nstars-1 do begin           ;Compute magnitudes for each star
   apmag = replicate(badval, Naper)   & magerr = replicate(baderr, Naper) 
   skymod = 0. & skysig = 0. &  skyskw = 0.  ;Sky mode sigma and skew
   if badstar[i] then goto, BADSTAR         
   error1=apmag   & error2 = apmag   & error3 = apmag

   rotbuf = image[ lx[i]:ux[i], ly[i]:uy[i] ] ;Extract subarray from image
;  RSQ will be an array, the same size as ROTBUF containing the square of
;      the distance of each pixel to the center pixel.

 
    dxsq = ( findgen( nx[i] ) - dx[i] )^2
    rsq = fltarr( nx[i], ny[i], /NOZERO )
   for ii = 0, ny[i]-1 do rsq[0,ii] = dxsq + (ii-dy[i])^2


 if keyword_set(exact) then begin 
       nbox = lindgen(nx[i]*ny[i])
       xx = reform( (nbox mod nx[i]), nx[i], ny[i])
       yy = reform( (nbox/nx[i]),nx[i],ny[i])
       x1 = abs(xx-dx[i]) 
       y1 = abs(yy-dy[i])
  endif else begin 
   r = sqrt(rsq) - 0.5    ;2-d array of the radius of each pixel in the subarray
 endelse

;  Select pixels within sky annulus, and eliminate pixels falling
;       below BADLO threshold.  SKYBUF will be 1-d array of sky pixels
 if N_elements(SETSKYVAL) EQ 0 then begin

 skypix = ( rsq GE rinsq ) and ( rsq LE routsq )
 if keyword_set(nan) then skypix = skypix and finite(rotbuf) $
 else if chk_badpix then skypix = skypix and ( rotbuf GT badpix[0] ) and $
        (rotbuf LT badpix[1] )
 sindex =  where(skypix, Nsky) 
 Nsky =   Nsky < maxsky   ;Must be less than MAXSKY pixels
 if ( nsky LT minsky ) then begin                       ;Sufficient sky pixels?
    if not silent then $
        message,'There aren''t enough valid pixels in the sky annulus.',/con
    goto, BADSTAR
 endif
  skybuf = rotbuf[ sindex[0:nsky-1] ]     

  if keyword_set(meanback) then $
   meanclip,skybuf,skymod,skysig, $ 
         CLIPSIG=clipsig, MAXITER=maxiter, CONVERGE_NUM=converge_num  else $
     mmm, skybuf, skymod, skysig, skyskw, readnoise=readnoise
           
 

;  Obtain the mode, standard deviation, and skewness of the peak in the
;      sky histogram, by calling MMM.

 skyvar = skysig^2    ;Variance of the sky brightness
 sigsq = skyvar/nsky  ;Square of standard error of mean sky brightness

;If the modal sky value could not be determined, then all apertures for this
; star are bad

 if ( skysig LT 0.0 ) then goto, BADSTAR 

 skysig = skysig < 999.99      ;Don't overload output formats
 skyskw = skyskw >(-99)<999.9
 endif else begin
    skymod = setskyval[0]
    skysig = setskyval[1]
    nsky = setskyval[2]
    skyvar = skysig^2
    sigsq = skyvar/nsky
    skyskw = 0
endelse



 for k = 0,Naper-1 do begin      ;Find pixels within each aperture

    if ( edge[i] GE apr[k] ) then begin ;Does aperture extend outside the image?
       if keyword_set(EXACT) then begin
          mask = fltarr(nx[i],ny[i])
          good = where( ( x1 LT smallrad[k] ) and (y1 LT smallrad[k] ), Ngood)
          if Ngood GT 0 then mask[good] = 1.0
          bad = where(  (x1 GT bigrad[k]) or (y1 GT bigrad[k] )) ;Fix 05-Dec-05
          mask[bad] = -1

          gfract = where(mask EQ 0.0, Nfract) 
          if Nfract GT 0 then mask[gfract] = $
             PIXWT(dx[i],dy[i],apr[k],xx[gfract],yy[gfract]) > 0.0
          thisap = where(mask GT 0.0)
          thisapd = rotbuf[thisap]
          fractn = mask[thisap]
       endif else begin
;
          thisap = where( r LT apr[k] ) ;Select pixels within radius
          thisapd = rotbuf[thisap]
          thisapr = r[thisap]
          fractn = (apr[k]-thisapr < 1.0 >0.0 ) ;Fraction of pixels to count
          full = fractn EQ 1.0
          gfull = where(full, Nfull)
          gfract = where(1 - full)
          factor = (area[k] - Nfull ) / total(fractn[gfract])
          fractn[gfract] = fractn[gfract]*factor
       endelse

;Check aperture for pixels hitting the saturation level. If so, flag
;the aperture to place a lower limit on the flux (upper lim on mag)
;and set all NaN pixels to SATLEV.
       if keyword_set(satval) then begin
          if max(thisapd) ge satval and min(thisapd) lt satval then begin
             satsrc[i] = 1
             undef = where(~finite(thisapd,/nan),nundef)
             if nundef ne 0 then thisapd[undef] = satval
          endif 
          if min(thisapd) ge satval then satsrc[i] = 2
       endif 
;if i eq 580 then stop ;debug       

;     If the pixel is bad, set the total counts in this aperture to a large
;        negative number
;
       if keyword_set(NaN) then $
          badflux =  min(finite(thisapd)) EQ 0   $
       else if chk_badpix then begin
          minthisapd = min(thisapd, max = maxthisapd)
          badflux = (minthisapd LE badpix[0] ) or ( maxthisapd GE badpix[1])
       endif else badflux = 0
  
       if not badflux then $ 
          apmag[k] = total(thisapd*fractn) ;Total over irregular aperture
    endif
 endfor                         ;k
 if keyword_set(flux) then g = where(finite(apmag), Ng)  else $
    g = where(apmag NE badval, Ng)
 if Ng GT 0 then begin 
  apmag[g] = apmag[g] - skymod*area[g]  ;Subtract sky from the integrated brightnesses

   error1[g] = area[g]*skyvar   ;Scatter in sky values
   error2[g] = (apmag[g] > 0)/phpadu  ;Random photon noise 
   error3[g] = sigsq*area[g]^2  ;Uncertainty in mean sky brightness
   magerr[g] = sqrt(error1[g] + error2[g] + error3[g])

  if not keyword_set(FLUX) then begin
    good = where (apmag GT 0.0, Ngood)     ;Are there any valid integrated fluxes?
    if ( Ngood GT 0 ) then begin               ;If YES then compute errors
      magerr[good] = 1.0857*magerr[good]/apmag[good]   ;1.0857 = log(10)/2.5
      apmag[good] =  25.-2.5*alog10(apmag[good])  
   endif
 endif  
 endif

 BADSTAR:   
 
;Print out magnitudes for this star

 for ii = 0,Naper-1 do $              ;Concatenate mags into a string

    ms[ii] = string( apmag[ii],'+-',magerr[ii], FORM = fmt)
   if PRINT then  printf,lun, $      ;Write results to file?
      form = fmt3,  i, xc[i], yc[i], skymod, skysig, skyskw, ms
   if not SILENT then print,form = fmt2, $       ;Write results to terminal?
          i,xc[i],yc[i],skymod,ms

   sky[i] = skymod    &  skyerr[i] = skysig  ;Store in output variable
   mags[0,i] = apmag  &  errap[0,i]= magerr
 endfor                                              ;i

 if PRINT then free_lun, lun             ;Close output file

 return
 end
