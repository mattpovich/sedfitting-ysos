;Read in a fitter datafile containing candidate YSOs and identify
;"strong" candidate unresolved 4.5 um excess objects (UGOs) based on
;the color-color criteria of Povich et al. 2011. Set these fluxes to
;upper limits and write out a new datafile.

pro ugoflag, data_in, data_out, nwav=nwav, regions=regions, coord=coord

;INPUT
;     DATA_IN    'string' -- Path to fitter datafile containing 2MASS, IRAC,
;                           and [optionally] UKIDSS fluxes of sources.
;
;OUTPUT
;     DATA_OUT    'string' -- Path to updated fitter datafile containing 2MASS, IRAC,
;                           and [optionally] UKIDSS fluxes.
;
;KEYWORD PARAMETERS
;     NWAV=    integer  -- Specifies number of bands used in DATA_IN
;                          Allowed values are 7, 8, 10 (DEFAULT), 11,
;                          12, and 13.
;     REGIONS=     SWITCH -- If set, outputs a ds9 regionfile called data_in+'_ugoc.reg'
;                            containing positions of the UGO sources.
;     COORD=      'string' -- Coordinate type passed to MAKEREG_XY if
;                             /REGIONS is set. DEFAULT='galactic'.
;
;CALLS
  ;  MAGFROMDATA, MAKEREG_XY
;  
;VERSION HISTORY
;     Original - M. S. Povich         August 2011
;     Refined color cuts - MSP November 2011
;
;PRODUCTION VERSION v1.0
;   * Reconfigured code to call MAGFROM DATA instead of computing
;    magnitudes internally.
;   * Replaced UKIDSS with NWAV keyword.

  
if n_params() lt 2 then begin
   print,'syntax - ugoflag, data_in, data_out, nwav=[10], /regions'
   return
endif

;Read input datafile
  nlines = file_lines(data_in)
  lines = strarr(nlines)
  openr,u,data_in,/get_lun
  line = '' 
  for j=0L,nlines-1 do begin
     readf,u,line
     lines(j) = line
  endfor 
  close,u
  free_lun,u

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
  
;  bands = indgen(nband)

  print,'Inspecting ',strtrim(nlines),' sources from file: ',data_in,'.'


                                ;Find the start, end position for each
                                ;COLUMN
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
  
  ;; flags = intarr(nlines,nwav)
  ;; for i=0,nwav-1 do flags[*,i] = fix(strmid(lines,colstart[i+3],colend[i+3]))
  ;; ind_cc = where(flags[*,band1] eq 1 and flags[*,band2] eq 1 and flags[*,band3] eq 1,ncc)
  ;; if ncc eq 0 then begin
  ;;    print,'No sources to be examined. RETURNING.'
  ;;    return
  ;; endif 

  ind_ugo = where(mag1 - mag2 gt 1.1 AND mag1 - mag2 gt (1.9/1.2)*(mag2 - mag3),nugo)
  if nugo eq 0 then begin
     print,'No sources are candidate UGOs! RETURNING.'
     return
  endif 

;  ind_update = ind_cc[ind_ugo]
  ind_update = ind_ugo
  print,'Updating ',nugo,' sources with strong 4.5um excess colors in file: '+data_in

  lineschange = lines[ind_update]
  strput,lineschange,'3',colstart[band2 + 3]

  strput,lineschange,'1.0000E+00',colstart[nwav + 3 + 2*band2 + 1]
  lines[ind_update] = lineschange
  if keyword_set(regions) then begin
     if not keyword_set(coord) then coord = 'galactic'
     l = double(strmid(lineschange,colstart[1],colend[1]-colstart[1]))
     b = double(strmid(lineschange,colstart[2],colend[2]-colstart[2]))

     makereg_xy,l,b,data_out+'_flagged.reg',coord=coord,point='diamond',color='green'
  endif 

  openw,u,data_out,/get_lun
  for i=0L,nlines-1 do printf,u,lines[i]
  close, u
  free_lun, u

   print,'Wrote ',strtrim(nlines,2),' sources to file: ',data_out,'.'

end
