;Get single-band IR magnitudes of detected sources (not upper or lower
;limits) from a specified
;sedfitter data file.

pro magfromdata,data_in,band,mags,e_mags,nwav=nwav,mk=mk

;INPUT
;     DATA_IN    'string' -- Path to fitter datafile containing 2MASS, IRAC,
;                           and MIPS24 fluxes of sources.
;     BAND  integer  -- Number specifying the IR broad
;                           band for which magnitudes are desired. Must take
;                           on an integer value between 0 and 10. 
;
;OUTPUT
;     MAGS     float[N_SRC] -- Magnitude values of specified BAND.
;     E_MAGS     float[N_SRC] -- Uncertainties on MAGS
;
;KEYWORDs
;     NWAV=    integer  -- Specifies number of bands used in DATA_IN
;                          Allowed values are 7, 8, 10 (DEFAULT), 11,
;                          12, and 13.
;     /MK      (switch) -- Set if desired JHK mags are on the Mauna
;                          Kea system (e.g. UKIDSS photometry) instead
;                          of the (default) 2MASS
;                          system.  
;
;
;MSP wrote 7 March 2008
;  Added NWAV keyword April 2010
;  Force the MAGS output vector to contain one entry
;  per source, with NaN's for non-detections or limits.

;PRODUCTION v1.0  M. S. Povich 17 April 2019
; * Changed reading of input data file to conform with new column-based
;  standard formatting
; * Updated magnitude zero-points for more precise
; * Removed histogram plotting functionality and associated keywords. 
; * Added E_MAGS output  -- 3 May 2019
 
  
if n_params() lt 3 then begin
   print,'syntax: magfromdata, data_in, band, mags, e_mags, nwav=[10], /mk'
   return
endif


  if not keyword_set(nwav) then nwav = 10

  if band gt nwav then begin
     print,'BAND must have a value between 0 and 10 <= NWAV. RETURNING.'
     return
  endif 
  
;NAMES of selectable bands
      ;NOTE: Zero-point fluxes for ZYJHK(s) are taken from Table 3 of Pickles & Depagne (2010, PASP, 122, 1437), for Spitzer IRAC/MIPS we use IPAC values.
  if nwav le 11 then begin
;                  0   1   2     3   4   5       6       7     8        9       10
     bandnames = ['J','H','Ks','2J','2H','2Ks','[3.6]','[4.5]','[5.8]','[8.0]','[24]']
     zerof = 1.d3*[1577.,1050.,674.9,1577.,1050.,674.9,280.9,179.7,115.0,64.13,7.17]
     if keyword_set(mk) then begin
        bandnames[2] = 'K'
        zerof[0:2] = 1.d3*[1531,1006.,663.1]
     endif 
  endif else begin
;                 0   1    2   3   4    5    6     7     8        9      10      11     12
     bandnames = ['Z','Y','J','H','Ks','2J','2H','2Ks','[3.6]','[4.5]','[5.8]','[8.0]','[24]']
     zerof = 1.d3*[2128.,2072.,1577.,1050.,674.9,1577.,1050.,674.9,280.9,179.7,115.0,64.13,7.17]
     if keyword_set(mk) then begin
        bandnames[4] = 'K'
        zerof[2:4] = [1531,1006.,663.1]
     endif 
  endelse
     

;Read input datafile
  case 1 of             ;Dealing with multiple options for filters included in data file
     nwav eq 7: begin
        readcol, data_in, $
           format='A,F,F,I,I,I,I,I,I,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F', $
           name, l, b, $
           ff1,ff2,ff3,ff4,ff5,ff6,ff7, $
           f1,df1,f2,df2,f3,df3,f4,df4,f5,df5,f6,df6,f7,df7
        flag = [[ff1], [ff2], [ff3], [ff4], [ff5], [ff6],[ff7] ]
        f = [[f1], [f2], [f3], [f4], [f5], [f6],[f7] ]
        df = [[df1], [df2], [df3], [df4], [df5], [df6],[df7] ]
     end
     nwav eq 8: begin
        readcol, data_in, $
           format='A,F,F,I,I,I,I,I,I,I,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F', $
           name, l, b, $
           ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8, $
           f1,df1,f2,df2,f3,df3,f4,df4,f5,df5,f6,df6,f7,df7,f8,df8
        flag = [[ff1], [ff2], [ff3], [ff4], [ff5], [ff6], [ff7], [ff8] ]
        f = [[f1], [f2], [f3], [f4], [f5], [f6], [f7], [f8] ]
        df = [[df1], [df2], [df3], [df4], [df5], [df6], [df7], [df8] ]
        zerof = [1577.,1050.,674.9,280.9,179.7,115.0,64.13]
     end
     nwav eq 10: begin
        readcol, data_in, $
           format='A,F,F,I,I,I,I,I,I,I,I,I,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F', $
           name, l, b, $
           ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10, $
           f1,df1,f2,df2,f3,df3,f4,df4,f5,df5,f6,df6,f7,df7,f8,df8,f9,df9,f10,df10 
        flag = [[ff1], [ff2], [ff3], [ff4], [ff5], [ff6], [ff7], [ff8], [ff9], [ff10] ]
        f = [[f1], [f2], [f3], [f4], [f5], [f6], [f7], [f8], [f9], [f10] ]
        df = [[df1], [df2], [df3], [df4], [df5], [df6], [df7], [df8], [df9], [df10] ]
     end 
     nwav eq 12: begin
        readcol, data_in, $
           name, l, b, $
           format='A,F,F,I,I,I,I,I,I,I,I,I,I,I,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F', $
           ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12, $
           f1,df1,f2,df2,f3,df3,f4,df4,f5,df5,f6,df6,f7,df7,f8,df8,f9,df9,f10,df10,f11,df11,f12,df12 
        flag = [[ff1], [ff2], [ff3], [ff4], [ff5], [ff6], [ff7], [ff8], [ff9], [ff10], [ff11], [ff12] ]
        f = [[f1], [f2], [f3], [f4], [f5], [f6], [f7], [f8], [f9], [f10], [f11], [f12] ]
        df = [[df1], [df2], [df3], [df4], [df5], [df6], [df7], [df8], [df9], [df10], [df11], [df12] ]
     end
     nwav eq 11: begin
        readcol, data_in, $
           format='A,F,F,I,I,I,I,I,I,I,I,I,I,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F', $
           name, l, b, $
           ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11, $
           f1,df1,f2,df2,f3,df3,f4,df4,f5,df5,f6,df6,f7,df7,f8,df8,f9,df9,f10,df10,f11,df11 
        flag = [[ff1], [ff2], [ff3], [ff4], [ff5], [ff6], [ff7], [ff8], [ff9], [ff10], [ff11] ]
        f = [[f1], [f2], [f3], [f4], [f5], [f6], [f7], [f8], [f9], [f10], [f11] ]
        df = [[df1], [df2], [df3], [df4], [df5], [df6], [df7], [df8], [df9], [df10], [df11] ]
     end
     nwav eq 13: begin
        readcol, data_in, $
           format='A,F,F,I,I,I,I,I,I,I,I,I,I,I,I,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F', $
           name, l, b, $
           ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13, $
           f1,df1,f2,df2,f3,df3,f4,df4,f5,df5,f6,df6,f7,df7,f8,df8,f9,df9,f10,df10,f11,df11,f12,df12,f13,df13
        flag = [[ff1], [ff2], [ff3], [ff4], [ff5], [ff6], [ff7], [ff8], [ff9], [ff10], [ff11], f[12], f[13] ]
        f = [[f1], [f2], [f3], [f4], [f5], [f6], [f7], [f8], [f9], [f10], [f11], [f12] ]
        df = [[df1], [df2], [df3], [df4], [df5], [df6], [df7], [df8], [df9], [df10], [df11], [df12], [df13] ]
     end    
     else: begin
        print,'NWAV must have a value of 7, 8, 10, 11, 12, or 13. RETURNING.'
        return
     end 
  endcase 
  nlines = n_elements(name)

  flag = flag[*,band]
  f = f[*,band]
  df = df[*,band]
  ind_det = where(flag eq 1,nf)

  mags = replicate(!VALUES.F_NAN,nlines)
  e_mags = replicate(!VALUES.F_NAN,nlines)
  mags[ind_det] = -2.5*alog10(f[ind_det]/zerof[band])  
  e_mags[ind_det] = df[ind_det]/f[ind_det]
  
  print,strtrim(nf,1),'/',strtrim(nlines,1),' sources have detections in the '+bandnames[band]+' bandpass.'

end
