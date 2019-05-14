;Turn ON unused photometry datapoints by changing their flags to 1

pro data_activate, data_in, data_out, nwav=nwav

  if n_params() lt 2 then begin
     print,'syntax -- data_activate, data_in, data_out, nwav=[10]'
     return
  endif

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

  ;Pull flags and fluxes from data lines
  flags = intarr(nlines,nwav)
  fluxes = fltarr(nlines,nwav)
  for j=0,nwav-1 do begin
     flags[*,j] = fix(strmid(lines,colstart[j+3],colend[j+3]))
     fluxes[*,j] = float(strmid(lines,colstart[2*j+3+nwav],colend[2*j+3+nwav]))
  endfor 

  ;Identify flags that need activation and set to 1
  ind_activate = where(flags eq 0 and fluxes ne 0,n_activate)
  flags[ind_activate] = 1

                                ;Put new flags back into data lines and write out new file
  openw,v,data_out,/get_lun
  for i=0L,nlines-1 do begin
     line = lines[i]
     strput,line,string(flags[i,*],format='('+strtrim(nwav)+'(I1,1X))'),colstart[3]
     printf,v,line
  endfor 
  close, v
  free_lun, v

   print,'Wrote ',strtrim(nlines,2),' sources to file: ',data_out,'.'


end
