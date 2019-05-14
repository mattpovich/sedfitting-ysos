;Based on a few lines of code in Pat Broos'
;ae_ds9_to_ciao_regionfile utility that output the vertex coordinates
;for ds9 polygon regions.
;NOTE -- This only works with coordinates saved in DECIMAL DEGREES!
;HANDLES MULTIPLE POLYGONS WITHIN SINGLE REGION FILE.

pro ds9_polygonvertices,ds9_filename,chosen_polygon_x,chosen_polygon_y,choose=choose, n_poly=n_poly

if not keyword_set(choose) then choose = 1
  
num_lines = file_lines(ds9_filename)
lines     = strarr(num_lines)
chosen_polygon_x = 0
chosen_polygon_y = 0
count = 0
  
;; Parse ds9 regionfile.
openr, region_unit, ds9_filename, /GET_LUN

    line = ''
    for i=0,num_lines-1 do begin
      readf, region_unit, line      
      
      ;; Try to determine if ds9 wrote the file.
      ;; Skip globals, and entries with "background" tag if desired.
      if                                        (strpos(line,'global')     NE -1) then continue
      
      if keyword_set(ignore_background_tag) AND (strpos(line,'background') NE -1) then continue
      
      ;; Strip off trailing comments
      comment_start = strpos(line,'#')
      if (comment_start GT 0) then line = strmid(line,0,comment_start)
      
      ;; Strip off leading "physical;" tags present in ds9 version 3.
      if (strmid(line,0,9) EQ 'physical;') then line = strmid(line,9)
      
      ;; Skip lines that are in celestial coordinates (e.g. "panda" regions written by AE's PSF hook screening).
;      if (strmid(line,0,5) EQ 'J2000') then continue
      
      ;; Save line.
      lines[i] = line

      first = strpos(line, 'polygon(')
      if (first NE -1) then begin
         count++
      ; Parse the polygon.
      ; Check for leading '-' or '+'.
         sign_str = strmid(line, 0, first)
         case 1 of
            (strpos(sign_str,'+') NE -1): sign='+'
            (strpos(sign_str,'-') NE -1): sign='-'
            else:                         sign=''
         endcase

         first   = first + 8
         last    = strpos(line, ')', first)
         coords  = strmid(line, first, (last-first))
         numbers = float(strsplit(coords, ',', /EXTRACT))
         ind       = 2*indgen(n_elements(numbers)/2)
         polygon_x = numbers[ind]
         polygon_y = numbers[ind+1]

         if count eq choose then begin
            chosen_polygon_x = polygon_x
            chosen_polygon_y = polygon_y 
         endif 

    ;     if not keyword_set(chosen_polygon_x) then chosen_polygon_x = polygon_x 
    ;     if not keyword_set(chosen_polygon_y) then chosen_polygon_y = polygon_y 
         
      endif
   endfor 
    close,region_unit
    free_lun, region_unit

    n_poly = count

 end 
