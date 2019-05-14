;Generate a fitter datafile from a custom Spitzer-UKIDSS matched
;catalog created by Pat Broos' 2019 workflow.


pro glimpseplus2data, glimpseplusfile, subset_type, data_out, subset=subset, satlim=satlim,vvv=vvv

;INPUT
;       GLIMPSEPLUSRILE 'string' -- Path to IDL save file containing
;                                   a special, UKIDSS-matched
;                                   GLIMPSE point-source list as a
;                                   structure of structures. E.G. '../counterparts/spitzer_ukidss_cat.sav'
;       SUBSET_TYPE     INTEGER  -- Select subset of sources in
;                                   SOURCE_LIST to write to
;                                   DATA_OUT. Choices are:
;                0: Output all sources.
;                1: Perform cut using a specified coordinate box.
;                2: Select sources from list of designations.
;
;OUTPUT
;       DATA_OUT        'string' -- Path to output fitter datafile to
;                                   write.
;
;KEYWORD
;       SUBSET=  FLOAT or 'string' -- Specifies subset according to
;                                     SUBSET_TYPE. If SUBSET_TYPE=1
;                                     then SUBSET = [l, b, dl, db]
;                                     decimal degrees specifying
;                                     coordinate box. If SUBSET_TYPE=2
;                                     then SUBSET is a string
;                                     specifying the path to a file
;                                     containing a list of GLIMPSE source
;                                     designations.
;       SATLIM=  FLOAT[3 or 5] -- Limiting 2MASS above which UKIDSS
;                                 counterpart photometry should NOT be
;                                 trusted due to suspected saturation
;                                 effects. DEFAULT = [11,12,10.5] for UKIDSS
;       /VVV   (switch) -- If set, assumed VVV
;                          data columns, not the UKIDSS default.
;
;NOTE: If SUBSET_TYPE is nonzero and SUBSET keyword is not set, then
;SUBSET will be specified interactively.
;
;VERSION HISTORY  
;  - Original: MSP wrote 15 Oct 2009
;  - PRODUCTION (v1.0)   8 May 2019 (seriously, 10 years??)
  
if n_params() lt 3 then begin
   print,'syntax -- glimpseplus2data, glimpseplusfile, subset_type, data_out, subset=[0],satlim=[11,12,10.5],/vvv'
   return
endif

;Set saturation limits by default
  if not keyword_set(satlim) then begin
     if not keyword_set(vvv) then $
        satlim = [11,12,10.5] $
     else $
        satlim = [12.5,13,11.5]
  endif    

;Load glimpse+ catalog structure
  guv = mrdfits(glimpseplusfile,1)

  nsrc = n_elements(guv)

;Determine subset for output
  case subset_type of
     0: begin                   ; Output ALL sources (SUBSET_TYPE = 0)
        whrsub = lindgen(nsrc)
        nsub = nsrc
     end 

     1: begin                   ; Output sources in specified coordinate box (SUBSET_TYPE = 1)
        if not keyword_set(subset) then begin
           subset = fltarr(4)
           read,prompt='Specify coordinate box: [l, b, dl, db] = ',subset
        endif else begin
           szsub = size(subset)
           if szsub[1] ne 4 or szsub[2] ne 4 then begin
              print,'Keyword SUBSET must be a 4-element floating-point vector or left undefined. RETURNING.'
              return
           endif
        endelse
        lminmax = [ subset[0]-0.5*subset[2], subset[0]+0.5*subset[2] ]
        bminmax = [ subset[1]-0.5*subset[3], subset[1]+0.5*subset[3] ]
        whrl = where(guv.Glon ge lminmax[0] and guv.Glon le lminmax[1],nwl)
        if nwl eq 0 then begin
           print,'WARNING: No sources found in longitude cut. RETURNING'
           return
        endif 
        whrb = where(guv[whrl].Glat ge bminmax[0] and guv[whrl].Glat le bminmax[1],nsub)
        if nsub eq 0 then begin
           print,'WARNING: No sources found in latitude cut! RETURNING'
           return
        endif 
        whrsub = whrl[whrb]
     end

     2: begin                   ; Output sources from designation list (SUBSET_TYPE = 2)
        if not keyword_set(subset) then begin
           subset = ''
           read,prompt='Specify file containing GLIMPSE source designations (GLLL.llll+BB.bbbb) = ',subset
        endif else begin
           szsub = size(subset)
           if szsub[0] ne 0 or szsub[2] ne 1 then begin
              print,'Keyword SUBSET must be a string or left undefined. RETURNING.'
              return
           endif
        endelse
        openr,1,subset
        nlines = numlines(subset)
        line = ''
        lines = strarr(nlines)
        for i=0L,nlines-1 do begin
           readf,1,line
           lines(i) = line
        endfor
        close,1
        linelen = strlen(lines[0])
        deslen = strlen(guv[0].GLIMPSE)
        if linelen ne deslen then begin
           print,'Designations in file '+subset+' MUST match designations in DESIG tag of input structure. RETURNING.'
           return
        endif 
        desigtr = strmid(guv.GLIMPSE,deslen-17,17)
        linestr = strmid(lines,linelen-17,17)
        keep = replicate(-1L,nlines)
        for i=0L,nlines-1 do begin
           whri = where(desigtr eq linestr[i])
           keep(i) = whri
        endfor
        whrkeep = where(keep ne -1,nsub)
        if nsub eq 0 then begin
           print,'WARNING: No sources found from '+subset+'! RETURNING.'
           return
        endif 
        whrsub = keep[whrkeep]
        ;sort in order of increasing l
        lsub = guv[whrsub].Glon
        lsort = sort(lsub)
        whrsub = whrsub[lsort]
     end

     else: begin
        print, 'SUBSET_TYPE must have value of 0, 1, or 2. RETURNING'
        return
     end 
  endcase

;Locate sources in subset
  guv = guv[whrsub]
  f = [ [guv.F_J_], [guv.F_H_], [guv.F_KS_], $
      [guv.F_3_6_], [guv.F_4_5_], [guv.F_5_8_], [guv.F_8_0_]]
  df = [ [guv.e_F_J_], [guv.e_F_H_], [guv.e_F_KS_], $
      [guv.e_F_3_6_], [guv.e_F_4_5_], [guv.e_F_5_8_], [guv.e_F_8_0_]]

;Set up output arrays
  ind_det = where(finite(f) and f gt -90.,complement=ind_nodet)
  keys = fix(f)
  keys[ind_det] = 1
  keys[ind_nodet] = 0
  f[ind_nodet] = 0.
  df[ind_nodet] = 0.

;New arrays to hold UKIDSS/VVV fluxes  
  f_nan = !VALUES.f_nan
  f_uvk = replicate(f_nan,nsub)
  df_uvk = replicate(f_nan,nsub)
  f_uvh = replicate(f_nan,nsub)
  df_uvh = replicate(f_nan,nsub)
  f_uvj = replicate(f_nan,nsub)
  df_uvj = replicate(f_nan,nsub)
  if keyword_set(vvv) then begin
     f_vz = replicate(f_nan,nsub)
     df_vz = replicate(f_nan,nsub)
     f_vy = replicate(f_nan,nsub)
     df_vy = replicate(f_nan,nsub)
  endif 
  
; WHERE A HIGH-QUALITY UKIDSS/VVV MAG IS AVAILABLE, Calculate its FLUX in mJy

  pperrmax = 4 + 16 + 64        ;Only Byte 0 bits are acceptable (no warnings)
  if not keyword_set(vvv) then $
     ind_uvk = where(finite(guv.k_1Apermag3) and guv.k_1Apermag3 ge satlim[2] $
                  and guv.k_1ppErrBits le pperrmax,n_uvk)
  ind_uvh = where(finite(guv.hApermag3) and guv.hApermag3 ge satlim[1] $
                  and guv.hppErrBits le pperrmax,n_uvh)
  ind_uvj = where(finite(guv.jApermag3) and guv.jApermag3 ge satlim[0] $
                  and guv.jppErrBits le pperrmax,n_uvj)

  if keyword_set(vvv) then begin
     ind_uvk = where(finite(guv.ksApermag3) and guv.ksApermag3 ge satlim[2] $
                     and guv.ksppErrBits le pperrmax,n_uvk)
     ind_vz = where(finite(guv.zApermag3) and guv.jApermag3 ge satlim[0] $
                     and guv.zppErrBits le pperrmax,n_uvz)
     ind_vy = where(finite(guv.yApermag3) and guv.jApermag3 ge satlim[0] $
                    and guv.yppErrBits le pperrmax,n_uvy)

     print,'Found ',[n_uvz,n_uvy,n_uvj,n_uvh,n_uvk], $
           ' high-quality VVV ZYJHKs counterparts, respectively.'
  endif else $
     print,'Found ',[n_uvj,n_uvh,n_uvk], $
           ' high-quality UKIDSS J, H, K, counterparts, respectively.'
     
;Convert from mag to mJy (zeropoints from Pickles & Depagne 2010, PASP, 122, 1437)
  if keyword_set(vvv) then $
     f0 = 1.e3*[1577.,1050.,674.9,2128.,2072.] $ ;J,H,Ks + Z,Y zeropoints in mJy -- 2MASS/VVV
  else $
     f0 = 1.e3*[1531,1006.,663.1] ;J,H,K zeropoints in mJy -- UKIDSS
  
  ;Fluxes (mJy)
  if not keyword_set(vvv) then $
     f_uvk = f0[2]*10^(-1*guv[ind_uvk].k_1Apermag3/2.5)
  f_uvh = f0[1]*10^(-1*guv[ind_uvh].hApermag3/2.5)
  f_uvj = f0[0]*10^(-1*guv[ind_uvj].jApermag3/2.5)
  
  ;Uncertainties (mJy)
  if not keyword_set(VVV) then $
     df_uvk = f_uvk*guv[ind_uvk].k_1Apermag3Err
  df_uvh = f_uvh*guv[ind_uvh].hApermag3Err
  df_uvj = f_uvj*guv[ind_uvj].jApermag3Err

  if keyword_set(vvv) then begin ;Add Ks, Z, Y from VVV
     f_uvk = f0[2]*10^(-1*guv[ind_uvk].ksApermag3/2.5)
     f_vz = f0[3]*10^(-1*guv[ind_vz].zApermag3/2.5)
     f_vy = f0[4]*10^(-1*guv[ind_vy].yApermag3/2.5)

     df_uvk = f_uvk*guv[ind_uvk].ksApermag3Err
     df_vz = f_vz*guv[ind_vz].zApermag3Err
     df_vy = f_vy*guv[ind_vy].yApermag3Err
  endif 
     
;Expand flux and flag arrays, populate with UKIDSS fluxes.
  if keyword_set(vvv) then n_band = 12 else n_band = 10
  fexpand = fltarr(nsub,n_band)
  fexpand[*,n_band-7:n_band-1] = f
  f = temporary(fexpand)
  dfexpand = fltarr(nsub,n_band)
  dfexpand[*,n_band-7:n_band-1] = df
  df = temporary(dfexpand)
  keysexpand = intarr(nsub,n_band)
  keysexpand[*,n_band-7:n_band-1] = keys
  keys = temporary(keysexpand)

  f[ind_uvk,n_band-8] = f_uvk
  df[ind_uvk,n_band-8] = df_uvk
  keys[ind_uvk,n_band-8] = 1
  keys[ind_uvk,n_band-5] = 0    ;Suppress 2MASS when UKIDSS/VVV present

  f[ind_uvh,n_band-9] = f_uvh
  df[ind_uvh,n_band-9] = df_uvh
  keys[ind_uvh,n_band-9] = 1
  keys[ind_uvh,n_band-6] = 0   

  f[ind_uvj,n_band-10] = f_uvj
  df[ind_uvj,n_band-10] = df_uvj
  keys[ind_uvj,n_band-10] = 1
  keys[ind_uvj,n_band-7] = 0
  
  if keyword_set(vvv) then begin
     f[ind_vz,0] = f_vz
     df[ind_vz,0] = df_vz
 ;    keys[ind_vz,0] = 1      ;We'll switch these on later, for YSO fitting!
     
     f[ind_vy,1] = f_vy
     df[ind_vy,1] = df_vy
 ;    keys[ind_vy,1] = 1
  endif
     
;Write output datafile
  openw,u,data_out,/get_lun
  for i=0L,nsub-1 do $
     if keyword_set(vvv) then $
        printf,u,$
               format='(A25,A5,2(1X,F9.5),12(1X,I1),24(1X,E11.4))' $
               ,guv[i].GLIMPSE,' ',guv[i].Glon,guv[i].Glat,keys[i,*] $
               ,f[i,0],df[i,0],f[i,1],df[i,1],f[i,2],df[i,2] $
               ,f[i,3],df[i,3],f[i,4],df[i,4],f[i,5],df[i,5] $
               ,f[i,6],df[i,6],f[i,7],df[i,7],f[i,8],df[i,8] $
               ,f[i,9],df[i,9],f[i,10],df[i,10],f[i,11],df[i,11] $
     else $
        printf,u,$
               format='(A25,A5,2(1X,F9.5),10(1X,I1),20(1X,E11.4))' $
               ,guv[i].GLIMPSE,' ',guv[i].Glon,guv[i].Glat,keys[i,*] $
               ,f[i,0],df[i,0],f[i,1],df[i,1],f[i,2],df[i,2] $
               ,f[i,3],df[i,3],f[i,4],df[i,4],f[i,5],df[i,5] $
               ,f[i,6],df[i,6],f[i,7],df[i,7],f[i,8],df[i,8] $
               ,f[i,9],df[i,9]
  close,u
  free_lun,u
  
  print,'Wrote ',strtrim(nsub,2),' sources to file ',data_out, $
        ' from ',strtrim(nsrc,2),' in catalog.'

end
