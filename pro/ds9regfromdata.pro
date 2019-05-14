;Make a DS9 regionfile from a fitter data file.

pro ds9regfromdata,datafile,regfile,coord=coord,color=color

;INPUT
;     DATAFILE    'string': fitter data file to read. Example: '/d/glimpse62/povich/G34/fitter_v3.0e/output_s/data_good'
;
;OUTPUT
;     REGFILE     'string': DS9 regionfile to write. Example: 'stellar_good.reg'
;KEYWORDS
;     COORD       'string': Coordinate type for regionfile. Default='galactic'
;     COLOR       'string': Color of the points in the regionfile.
;                      'blue' = well fit as stars (default)
;                      'red'  = well fit as YSOs
;                      'magenta' = poorly fit as YSOs
;
;VERSION HISTORY
  ;ORIGINAL B. Whitney wrote MAKEREG, lost in the mists of time
   ;PRODUCTION v1.0   M. S. Povich,
;  * Changed to program name to DS9REGFROMDATA
;  * Minor housekeeping and modernization to handle now-standard data file format where the first column is not longer fixed with but cannot contain spaces.
  
if n_params() lt 2 then begin
    print,'syntax: ds9regfromdata, datafile, regfile [, coord=, color=]'
    return
endif

  nlines = file_lines(datafile)
  line = ''
  xcoord = dblarr(nlines)
  ycoord = dblarr(nlines)
  openr,u,datafile,/get_lun
  for i=0L,nlines-1 do begin
     readf,u,line
     line = strtrim(line,2)
     xstart = strpos(line,' ')
     line = strtrim(strmid(line,xstart),1)
     xend = strpos(line,' ')
     xcoord[i] = double(strmid(line,0,xend))
     line = strtrim(strmid(line,xend),1)
     yend = strpos(line,' ')
     ycoord[i] = double(strmid(line,0,yend))
  endfor

  close,u
  free_lun,u

  if not keyword_set(coord) then coord='galactic'
  if not keyword_set(color) then color='blue'

  openw,v,regfile,/get_lun
  printf,v,'# Region file format: DS9 version 3.0'
  printf,v,'global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source'

  for i=0l,n_elements(xcoord)-1l do begin
     printf,v,coord+';circle('+strtrim(string(xcoord[i],format='(f10.5)'),2)+','+$
            strtrim(string(ycoord[i],format='(f10.5)'),2)+',6.5")'
  endfor

  close,v
  free_lun,v
  
end

