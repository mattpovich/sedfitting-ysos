;Make a new file in sedfitter data format from an ASCII file make from
;the unformatted binary *.fitinfo output by the Python sedfitter.
;Same functionality as the old fits2data fortran utility which
;doesn't seem to have been implemented for the python sedfitter.

pro fitinfo2data, fitinfo_ascii, data_parent, data_out

;INPUT
;        FITINFO_ASCII      'string' -- Path to ASCII file produced
;                                       using the python
;                                       write_parameters() function.
;        DATA_PARENT        'string' -- Path to existing fitter
;                                       datafile that contains ALL of
;                                       the sources in FITINFO_ASCII
;                                       (must be same set OR superset)
;OUTPUT
;        DATA_OUT           'string' -- Path to new fitter datafile to
;                                       write.
;
;VERSION HISTORY
;   Original 24 July 2018 by MSP

  if n_params() lt 3 then begin
     print,'syntax -- fitinfo2data, fitinfo_ascii, data_parent, data_out'
     return
  endif

;Cleverly taking advantage of READCOL to pull only the rows we want
  readcol,fitinfo_ascii,format='A,I,I',source_name,ndata,nfits,/silent
  nsrc = n_elements(source_name)
                                ;Working around a weird bug for
                                ;READCOL format in the PMS model name case :-(
  namelen = strlen(source_name)
  if namelen[1] ne namelen[0] then begin
     ind_valid = where(namelen eq namelen[0],nsrc)
     source_name = source_name[ind_valid]
     ndata = ndata[ind_valid]
     nfits = nfits[ind_valid]
     endif
  print,nsrc,' sources in file: ',fitinfo_ascii ;added benefit of REPORTING how many sources we have, since python filter_output doesn't!

;Actually need the FULL lines and then to pull names from them!
  nlines = file_lines(data_parent)
  line = ''
  datalines = strarr(nlines)
  openr,u,data_parent,/get_lun
  for i=0L,nlines-1 do begin
     readf,u,line
     datalines[i] = line
  endfor

  dataname = strmid(strtrim(datalines,1),0,strpos(strtrim(datalines[0],1),' '))

  close,u
  free_lun,u
  
;Prepare output file
  openw,v,data_out,/get_lun

  count = 0L
  for i=0L,nsrc-1 do begin      ;loop over all sources fit
     ind_match = where(dataname eq source_name[i],nmatch)
     if nmatch eq 1 then begin
        printf,v,datalines[ind_match]
        count++
     endif
     ;BUG CHECK!
     if nmatch eq 0 then begin
        print,'WARNING! No match found for source ',source_name[i]
        print,'RETURNING. Check that FITINFO_ASCII is a subset of DATA_PARENT.'
        STOP
        return
     endif 
  endfor
  close,v
  free_lun,v
  print,'Wrote ',count,' lines to: ', data_out
  
end
