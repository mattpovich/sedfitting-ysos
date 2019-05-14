;Make a color cut from a fitter datafile containing a list of sources
;poorly-fit by stellar photospheres because of IR excess only in the single
;longest-wavelength IRAC band (either [8.0] or [5.8]).

pro malmcullav,data_in,data_good,data_bad,maxav=maxav,nwav=nwav

;INPUT
;     DATA_IN    'string' -- Path to fitter datafile containing 2MASS, IRAC,
;                           and [optionally] UKIDSS fluxes of sources.
;
;OUTPUT
;     DATA_GOOD    'string' -- Path to fitter datafile containing 2MASS, IRAC,
;                           and [optionally] UKIDSS fluxes of candidate YSOs to keep.
;     DATA_BAD    'string' -- Path to fitter datafile containing 2MASS, IRAC,
;                           and [optionally] UKIDSS fluxes of likely bad data
;                           (Malmquist bias) sources.
;
;KEYWORDS
;     MAXAV=       FLOAT -- Max Av used in fitting stellar photosphere
;                           models to data. Will be applied to estimate
;                           the [3.6]-[4.5] cutoff for reddened
;                           stellar photospheres (DEFAULT = 40.).
;
;     NWAV=    integer  -- Specifies number of bands used in DATA_IN
;                          Allowed values are 7, 8, 10 (DEFAULT), 11,
;                          12, and 13.
  

;CALLS
;   MAGFROMDATA
  
;MSP wrote 12 May 2009
;
;Added MAXAV keyword 12 Nov 2009
;Added UKIDSS keyword 24 August 2011
;Added /SILENT keyword and changed COUNTG,COUNTB to LONG 24 July 2018
;PRODUCTION VERSION v1.0
;   * Reconfigured code to call MAGFROM DATA instead of computing
;    magnitudes internally.
;   * Replaced UKIDSS with NWAV keyword.
;   * Removed SILENT keyword.
  
if n_params() lt 3 then begin
   print,'syntax -- malmcullav, data_in, data_good, data_bad, maxav=[40], nwav=[10]'
   return
endif

  if not keyword_set(maxav) then maxav = 40.


;Prepare output datafiles
  openw,u,data_good,/get_lun
  openw,v,data_bad,/get_lun
  countg = 0L
  countb = 0L

;Get IRAC magnitudes
  case 1 of
     nwav eq 7 or nwav eq 8: begin
        band1 = 3
        band2 = 4
        band3 = 5
        band4 = 6
     end
     nwav eq 10 or nwav eq 11: begin
        band1 = 6
        band2 = 7
        band3 = 8
        band4 = 9
     end
     nwav eq 12 or nwav eq 13: begin
        band1 = 8
        band2 = 9
        band3 = 10
        band4 = 11       
     end
     else: begin
        print,'NWAV must have a value of 7, 8, 10, 11, 12, or 13. RETURNING.'
        return
     end 
  endcase 

  magfromdata,data_in,band1,mag1,emag1,nwav=nwav
  magfromdata,data_in,band2,mag2,emag2,nwav=nwav
  magfromdata,data_in,band3,mag3,emag3,nwav=nwav
  magfromdata,data_in,band4,mag4,emag4,nwav=nwav

;NEW! Extinction correction
  filtkapremy = [0.56, 0.43] 
  filtkapremy = filtkapremy*23. ;scale to ISM at K
  kapvremy=220.
  maxa1 = maxav*filtkapremy[0]/kapvremy
  maxa2 = maxav*filtkapremy[1]/kapvremy

  nlines = file_lines(data_in)
  openr,w,data_in,/get_lun
  line = ''
  print,'Inspecting ',strtrim(nlines),' sources from file: ',data_in,'.'

  for i=0L,nlines-1 do begin

     readf,w,line
     
     if ~finite(mag1[i]) or ~finite(mag2[i]) then begin ;Check for bands 1+2
        printf,v,line       ;TOSS -- too many bad ones to cull out manually
        countb++
     endif else begin           ;Continue if source detected in both bands 1+2

        color21 = (mag1[i] - maxa1) - (mag2[i] - maxa2) ;Color corrected for max Av
        ecolor21 = sqrt(emag1[i]^2 + emag2[i]^2)
        
        if color21 gt ecolor21 then begin
           printf,u,line    ;KEEP For later visual SED inspection
           countg++
        endif else begin
           if ~finite(mag3[i]) or ~finite(mag4[i]) ne 1 then begin ;Check for bands 3+4
              printf,v,line ;TOSS to bad data sources
              countb++
           endif else begin     ;Continue if source detected in both bands 3+4
              
              color32 = mag2[i] - mag3[i]
              ecolor32 = sqrt(emag2[i]^2 + emag3[i]^2)

              color43 = mag3[i] - mag4[i]
              ecolor43 = sqrt(emag3[i]^2 + emag4[i]^2)

              if abs(color32) gt ecolor32 then begin
                 if color43 gt ecolor43 then begin ;MODIFICATION!
                    printf,u,line ;KEEP For later visual SED inspection
                    countg++
                 endif else begin
                    printf,v,line ;TOSS to bad data sources
                    countb++
                 endelse  
              endif else begin
                 if abs(color43) gt ecolor43 then begin
                    printf,v,line ;TOSS to bad data sources
                    countb++
                 endif else begin
                    printf,u,line ;KEEP For later visual SED inspection
                    countg++
                 endelse   
              endelse 
           endelse             
        endelse  
     endelse 
  endfor  
 
  print,'Wrote ',strtrim(countg,2),' sources to file: ',data_good,'.'
  print,'Wrote ',strtrim(countb,2),' sources to file: ',data_bad,'.'

  close,u
  close,v
  free_lun,u
  free_lun,v

end
