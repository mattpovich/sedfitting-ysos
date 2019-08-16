;Create a standard data table <target>.ised.fits as the
;output of the recipe sedfitting_procedure_ukidss-mir_PYTHON.txt used
;in conjunction with the Robitaille (2017) YSO models. This
;program is NOT general in the sense that it assumes typical combinations of filters in the order that I generally use, assumes Galactic coordinates in the input, and a bunch of other things like that. If anyone asks me to work on generalizing this later, I'll consider it!

;KEYWORD PARAMETERS
;TARGET_DIR            - 'string'  Path to the directory containing the SED fitting results as IDL save files. DEFAULT = 'target_ysoc'
;SOURCELIST_DIR        - 'string'  Name of the list of IDL save files to process. DEFAULT = 'sourcelist.txt'
;XFOV                  - 'string'  Path to an optional ds9 region file containing the polygon outlining an ACIS field-of-view.
;
;The columns of the saved table ised are as follows:
;(1)  MIR_NAME           = (A30) Source designation from Spitzer catalog.
;(2)  RA              = (D9.5) Right ascension of Spitzer source (J2000 deg)
;(3)  DEC              = (D9.5) Declination of Spitzer source (J2000 deg)
;(4)  IRMAG              = NBAND(F5.2) Vector of NIR+MIR mags used in the SED fitting (NaN=non-detection)
;(5)  IRMAG_ERR          = NBAND(F6.2) Uncertainties on IRMAG (NaN=non-detection, -99.99=upper limit used in fitting). 
;(X)  NIRPHOT_CAT        = 3(I2) Provenance of each JHK mag (0 = 2MASS, 1 = UKIDSS, -1 = none) DEPRECATED!
;(X)  NIRPHOT_NUM_SM     = (I1) Number of UKIDSS secondary matches to Spitzer catalog source. DEPRECATED!
;(6)  SED_FLAG            = (I1)  [-1] No models
;, 0 = likely YSO/diskless PMS, 1 = candidate galaxy/PAH, 2 = candidate AGN, (NEW) 3 = candidate AGB star
;(7)  SED_CHISQ_NORM     = (F6.3) Normalized chi-squared for best-fit model 
;(8) SED_AV             = (F7.2) Visual extinction (mag) to source from SED fits.
;(9) SED_AV_ERR          = (F7.2) Uncertainty on SED_AV (mag).
;(10) SED_STAGE          = (I2) 1 = Stage I (envelopes), 2 = Stage II
;(disks), -1 = ambiguous (I/II), [3] = Stage III (no MIR excess, PMS models)
;(XX) PROB_DENS          = (F5.3) Probability that source is a member
;based on local source density. DEPRECATED!

pro create_ised_r17,target_ysoc=target_ysoc,sourcelist_ysoc=sourcelist_ysoc,xfov=xfov, target_xpms=target_xpms, sourcelist_xpms=sourcelist_xpms,fk5=fk5

;KEYWORD PARAMETERS
;    TARGET_YSOC =     - 'string': Path to directory containing IDL save files with fitting results from R17 YSO models. DEFAULT = 'results_ysoc'
;    TARGET_XPMS =     - 'string': Optional path to directory (e.g. 'results_xpms' containing IDL save files with fitting results from PMS models.
;    SOURCELIST_YSOC = - 'string': Path to file within TARGET_YSOC that contains list of save file names to process. DEFAULT = 'sourcelist.txt'
;    SOURCELIST_XPMS = - 'string': Path to file within TARGET_XPMS that contains list of save file names to process. DEFAULT = 'sourcelist.txt' (unused if TARGET_XPMS is not specified).


;CALLS
;    GALFLAG, STAGEFLAG_R17, DS9_POLYGONVERTICES, MAKEREG_XY
;
;FILES CREATED
;    ised.sav, su_cat_yso_final_galflag.reg, su_cat_yso_final_stageflag.reg
;
;VERSION HISTORY
;  Adapted from CREATE_ISED
;  PRODUCTION (v1.0) - 6 May 2019   M. S. Povich
;               v1.1 -  16 May 2019 M. S. Povich
;    Added target_xpms=, sourcelist_xmps keyword functionality to
;    combine SCIM-X and MIRES catalog information into a single table.
;    Defined SED_FLAG < 0 and SED_STAGE = 3 for sources fith with PMS models.
;    Added /FK5 switch in case input coordinates are celesial, not
;    Galactic.
;               v1.2 - 15 August 2019 M. S. Povich
;    Patch to allow inclusion in table of photometry NOT used in SED fitting
;    (e.g. when H-band is ommitted).
  
  if not keyword_set(target_ysoc) then target_ysoc = 'results_ysoc'

  if not keyword_set(sourcelist_ysoc) then sourcelist_ysoc = 'sourcelist.txt'
  readcol,target_ysoc + '/' + sourcelist_ysoc,format='(A)',names
  files = target_ysoc + '/' + names
  n_src = n_elements(files)

  if keyword_set(target_xpms) then begin
     if not keyword_set(sourcelist_xpms) then sourcelist_xpms = 'sourcelist.txt'
     readcol,target_xpms + '/' + sourcelist_xpms,format='(A)',names2
     files = [files,target_xpms + '/' + names2]
     n_yso = n_src
     n_src = n_elements(files) 
  endif 
  
;Initialization for photometry variables
  restore,files[0]  ;JUST to get the number of filters right!
  nband = n_elements(s.valid)
  if nband le 11 then $  ;This is UKIDSS + 2MASS + IRAC + MIPS
     zerof = 1.d3*[1531,1006.,663.1,1577.,1050.,674.9,280.9,179.7,115.0,64.13,7.17] $
  else $   ; VVV + 2MASS + IRAC + MIPS
     zerof = 1.d3*[2128.,2072.,1577.,1050.,674.9,1577.,1050.,674.9,280.9,179.7,115.0,64.13,7.17]

  ;Below assumes a standard combination of UKIDSS/VVV + 2MASS + IRAC + MIPS24 photometry
  case 1 of
     nband eq 7: zerof = 1.e3*[1577.,1050.,674.9,280.9,179.7,115.0,64.13] ;2MASS+IRAC
     nband eq 8: zerof = 1.e3*[1577.,1050.,674.9,280.9,179.7,115.0,64.13,7.17] ;+MIPS24
     nband eq 10: zerof = 1.e3*[1531,1006.,663.1, $ ;+UKIDSS JHK
                                1577.,1050.,674.9,280.9,179.7,115.0,64.13]
     nband eq 11: zerof = 1.d3*[1531,1006.,663.1, $ ;+UKIDSS JHK
                                1577.,1050.,674.9,280.9,179.7,115.0,64.13,7.17] ;+MIPS24
     nband eq 12: zerof = 1.d3*[2128.,2072.,1577.,1050.,674.9, $ ;+VVV YZJHKs
                                1577.,1050.,674.9,280.9,179.7,115.0,64.13]
     nband eq 13: zerof = 1.d3*[2128.,2072.,1577.,1050.,674.9, $ ;+VVV YZJHKs 
                                1577.,1050.,674.9,280.9,179.7,115.0,64.13,7.17] ;+MIPS24
  endcase
             
  f_nan    = !VALUES.F_NAN
  d_nan    = !VALUES.D_NAN
  if not keyword_set(xfov) then begin
     template_row = {$
                 MIR_NAME      :''   ,$
                 RA         :d_nan,$
                 DEC         :d_nan,$
                 IRMAG         :replicate(f_nan,nband),$
                 IRMAG_ERR     :replicate(f_nan,nband),$
;                 NIRPHOT_CAT   :replicate(-1,3),   $
                 SED_FLAG       :-1   ,$
                 SED_CHISQ_NORM:f_nan,$
                 SED_AV        :f_nan,$
                 SED_AV_ERR    :f_nan,$
                 SED_STAGE     :-1   $
                    }
  endif else begin
     template_row = {$
                 MIR_NAME      :''   ,$
                 RA         :d_nan,$
                 DEC         :d_nan,$
                 IRMAG         :replicate(f_nan,nband),$
                 IRMAG_ERR     :replicate(f_nan,nband),$
;                 NIRPHOT_CAT   :replicate(-1,3),   $
                 SED_FLAG       :-1   ,$
                 SED_CHISQ_NORM:f_nan,$
                 SED_AV        :f_nan,$
                 SED_AV_ERR    :f_nan,$
                 SED_STAGE     :-1,   $
                 XFOV          :0 $
                    }
  endelse 

;  template_row = null_structure(template_row)
  iSED = replicate(template_row, n_src)

  for i=0L,n_src-1 do begin   ;Populate info source-by-source
     
     restore,files[i]

     desig = strtrim(s.desig,2)
     MIR_NAME = desig   ;ASSUMING this is in a correct naming convention! But it's just a string, so shouldn't crash anything.

     iSED[i].MIR_NAME = MIR_NAME

     if not keyword_set(fk5) then begin
        euler,s.l,s.b,RA,DEC,2  ;ASSUMES S.l and S.b are Galactic coords
        iSED[i].RA = RA
        iSED[i].DEC = DEC
     endif else begin
        iSED[i].RA = s.l
        iSED[i].DEC = s.b
     endelse 

     ;Convert fluxes to magnitudes and populate photometry arrays
     bands_fit = where(s.valid ne 0,nfit)  ;Original, but see below for MAGS.
     bands_present = where(s.f ne 0,npresent)
     bands_lim = where(s.valid eq 3,nlim)
     ndata = nfit - nlim

     mags = replicate(f_nan,nband)
     emags = replicate(f_nan,nband)
;     mags[bands_fit] = -2.5*alog10(s.f[bands_fit]/zerof[bands_fit])
;     emags[bands_fit] = s.df[bands_fit]/s.f[bands_fit]
     mags[bands_present] = -2.5*alog10(s.f[bands_present]/zerof[bands_present])
     emags[bands_present] = s.df[bands_present]/s.f[bands_present]
     if nlim ne 0 then emags[bands_lim] = -99.99

     iSED[i].IRMAG = mags
     iSED[i].IRMAG_ERR = emags 

     iSED[i].SED_CHISQ_NORM = min(p.chi2)/ndata

     weights = exp(-1*p.chi2/2.)   ;Simplest possible weighting..
     iSED[i].SED_AV = total(weights*p.av)/total(weights)
     iSED[i].SED_AV_ERR = sqrt(total(weights*(p.av - iSED[i].SED_AV)^2)/total(weights))
  endfor 

                                ;Flag possible starburst galaxies and AGN
  mags = transpose(iSED.IRMAG[nband-5:nband-2,*]) ;IRAC mags only
  galflag_phot,mags,iSED.SED_AV,SED_FLAG,regionfile='ised_sedflag.reg', $
               ra=iSED.RA,dec=iSED.DEC
  iSED.SED_FLAG   = SED_FLAG 
  ;; if keyword_set(target_xpms) then begin
  ;;    ind_stellar = where(SED_FLAG eq 0,n_stellar)
  ;;    ind_pms = ind_stellar[where(ind_stellar ge n_yso,n_pms)]
  ;;    if n_pms ne 0 then begin
  ;;       iSED[n_yso:*].SED_FLAG = SED_FLAG[ind_pms]
  ;;    endif 
  ;; endif

                                ;Flag possible AGB stars
  ind_agb = where(ised.irmag[nband-2] - ised.irmag[nband-1] lt 2.2,nagb,complement=ind_notagb)
  iSED[ind_agb].SED_FLAG = 3
  print,nagb,' candidate AGB stars'
  makereg_xy,iSED[ind_agb].ra,iSED[ind_agb].dec,'ised_sedflag.reg',point='cross',color='green',size=9,/append

  ;Classify stages
  stageflag_r17,target_ysoc,SED_STAGE,ambiguous=0.67,regionfile='ised_stageflag.reg',sourcelist=sourcelist_ysoc,fk5=keyword_set(fk5)
  if not keyword_set(target_xpms) then $
     iSED.SED_STAGE = SED_STAGE $
  else begin
     iSED[0:n_yso-1].SED_STAGE = SED_STAGE
     iSED[n_yso:*].SED_STAGE = 3
     makereg_xy, iSED[n_yso:*].RA, iSED[n_yso:*].DEC, 'ised_stageflag.reg',color='dodgerblue',/append,group='Stage III'
  endelse 

                                ;Flag YSO candidates within ACIS FOV (OPTIONAL!)
  ;CURRENTLY picky about input file! Buggy?
  if keyword_set(xfov) then begin
     choose = 1
     n_poly = 1
     while choose le n_poly do begin
        ds9_polygonvertices,xfov,polygon_x,polygon_y,choose=choose,n_poly=n_poly
        o_poly=obj_new('IDLanROI', polygon_x, polygon_y)
        ind_xfov = where((o_poly->ContainsPoints(iSED.RA, iSED.DEC) EQ 1), num_in,NCOMPLEMENT=num_trimmed)
        print,num_in,' iSED sources fell inside X-ray FOV polygon '+strtrim(choose,2)+'/'+strtrim(n_poly,2)
        iSED[ind_xfov].XFOV = 1
        choose++
     endwhile 
  endif 

  ;Sort by R.A. and save
  rasrt = sort(iSED.RA)
  iSED = iSED[rasrt]
  
;  save,iSED,file='ised_r17.sav',/verbose
  if not keyword_set(target_xpms) then $
     mwrfits,iSED,'ised_r17.fits',/create $
  else $
     mwrfits,iSED,'ised_r17+xpms.fits',/create
     
  
end
