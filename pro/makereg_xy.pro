pro makereg_xy, x, y, regionfile, coord=coord, radius=radius, point=point, color=color, size=size, group=group, append=append

;+  
;PURPOSE
;Simple routine inspired by a script Leisa Townsley showed me. Takes a
;set of paired coordinate vectors (x,y) and outputs a ds9 regionfile
;with circles marking each point.
;
;INPUT
;      X,Y      double(npoints) -- Paired coordinates in the system
;                                  specified by the COORD keyword.
;
;OUTPUT
;      REGIONFILE      'string' -- Path to the output regionfile to
;                                  write.
;
;KEYWORDS
;      COORD=          'string' -- The coordinate system for X and
;                                  Y. This can be anything that ds9
;accepts. DEFAULT='fk5', other options include 'fk4', 'galactic', and
;'physical'. If you choose something invalid this routine will still
;run, but ds9 will choke on your regionfile. You can fix the
;regionfile in a text editor, if necessary.
;
;      RADIUS=         float -- The radius of the regionfile
;                               circles. DEFAULT=6.5 which is in WCS
;arcseconds. If COORDINATE is a WCS system, then RADIUS will be
;converted to degrees, otherwise, it will be kept as entered (pixels).
;
;      POINT=          'string' -- If set, ds9 region points with the
;                                  specified symbol
;(e.g. 'circle','box','diamond') will be used in place of circles, and
;the RADIUS keyword is ignored. CAVEAT EMPTOR: YOU are responsible for
;knowing the valid point types that ds9 uses.
;
;      COLOR=          'string' -- Color for regions. DEFAULT =
;                                  'green'
;      SIZE=           integer -- If POINT is set, specifies the
;                                 symbol size (DEFAULT = 11)
;      GROUP=          'string' -- Adds a group tag to the end of each
;                                  line
;      /APPEND         switch -- If set, checks for the existence of
;                                REGIONFILE, and appends data to the file.
;
;VERSION_HISTORY
;   Written by M. S. Povich, 16 May 2011
;     Added GROUP and APPEND keywords 3 Nov 2015
;     Added ability to assign different radii to each circle 1 May 2018
;-


if n_params() lt 3 then begin
   print,'syntax -- makereg_xy, x, y, regionfile, coord=[fk5], radius=[6.5], point=, color=[green], size=[11],'
   print,'                                  group=, /append'
   return
endif

  if not keyword_set(color) then color = 'green'
  if not keyword_set(coord) then coord = 'fk5'
  if not keyword_set(group) then group = ''

  if keyword_set(point) and ~keyword_set(size) then size = 11

  if not keyword_set(radius) then radius = 6.5
;  if coord eq 'fk5' or coord eq 'galactic' or coord eq 'fk4' or coord eq 'ecliptic' then radius = radius/3600.

  n = n_elements(x)
  if n_elements(y) ne n then begin
     print,'Input vector Y MUST have same length as X! Returning.'
     return
  endif
  if n_elements(radius) eq 1 then radius = replicate(radius,n)
  if n_elements(group) eq 1 then group = replicate(group,n)

  openw,u,regionfile,/get_lun,append=keyword_set(append)
  if not keyword_set(append) then begin
     printf,u,'# Region file format: DS9 version 4.1'
     printf,u,'global color='+color+' dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
;     printf,u,coord
  endif 
  if not keyword_set(point) then begin
     for i=0L,n-1 do printf,u,strtrim(coord,2)+';circle('+strtrim(x[i],2)+', '+strtrim(y[i],2)+', '+strtrim(radius[i],2)+'") # color='+color+' tag={'+group[i]+'}'
  endif else begin
     for i=0L,n-1 do printf,u,strtrim(coord,2)+';point('+strtrim(x[i],2)+', '+strtrim(y[i],2)+') # point='+point+' '+strtrim(size,2) + ' color='+color+' tag={'+group[i]+'}'      
  endelse 
  close,u
  free_lun,u

  print,'Wrote ',n,' regions to file ',regionfile,'.'
  






end
