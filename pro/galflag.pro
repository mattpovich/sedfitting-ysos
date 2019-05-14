;Flag (1) candidate starburst galaxies and (2) AGN in a sample of YSOs
;that have been through the SED fitting procedure and processed
;through FITS2IDL_YSO.
;Strategy: deredden candidate YSOs using AV from the fits, set flags using simple color-mag criteria
;adapted from Gutermuth et al. (2008,2009).
;NOTE: This procedure can also be used on the output of stellar
;photosphere fits processed through FITS2IDL_STELLAR

pro galflag, target_dir, galflags, av, dav, sourcelist=sourcelist, noukidss=noukidss, regionfile=regionfile

;INPUT
;      TARGET_DIR        'string' --  Directory containing IDL save files of
;                                    fit parameters and sourcelist
;                                    produced by FITS2IDL_YSO.
;
;OUTPUT
;      GALFLAGS             INTEGER[nsrc] -- Source flags matching the order
;                                         of SOURCELIST.
;                                         0=candidate YSO,
;                                         1=candidate starburst galaxy
;                                         2= candidate AGN
;      AV                FLOAT[nsrc]   -- Visual extinction used in
;                                         dereddening
;                                         (probability-weighted MEAN
;                                         AV from ALL fits)
;      DAV               FLOAT[nsrc]   -- Weighted standard deviation
;                                         on AV.
;
;KEYWORDS
;     SOURCELIST=   'string' -- Name of file containing list of IDL save
;                                files to operate on (DEFAULT: 'sourcelist.txt')
;     NOUKIDSS=      SWITCH -- Set this switch if the fitter data DOES
;                              NOT CONTAIN BOTH
;                            UKIDSS photometry AND 2MASS
;                            photometry.
;     MAXAV=       FLOAT -- Max Av used in fitting stellar photosphere
;                           models to data. (DEFAULT = 0.).
;     REGIONFILE=   'string' -- Specifies path to write a DS9 regionfile showing all sources
;                               color-coded by classification: orange
;                               = YSOc, magenta = galc, blue = agnc
;
;CALLS
;
;VERSION HISTORY
;     Original - M. S. Povich         November 2011
;       Added DAV output parameter   May 2019  MSP
  
if n_params() lt 3 then begin
   print,'syntax - galflag, target_dir, flags, av, sourcelist=, /noukidss, regionfile='
   return
endif

  if not keyword_set(sourcelist) then sourcelist = 'sourcelist.txt'
  target = target_dir + '/'
  readcol,target + sourcelist,format='(A)',names
  files = target + names
  nsrc = n_elements(files)

;Convert fluxes to magnitudes, deredden, and make color cuts source-by-source 

  if keyword_set(noukidss) then begin
     zerof = [1594.,1024.,666.7,280.9,179.7,115.0,64.13]
     zeroe = [27.8,20.0,12.6,4.1,2.6,1.7,0.94,0.11]
     band1 = 3
     band2 = 4
     band3 = 5
     band4 = 6
  endif else begin
     zerof = [1594.,1024.,666.7,1594.,1024.,666.7,280.9,179.7,115.0,64.13]
     zeroe = [27.8,20.0,12.6,27.8,20.0,12.6,4.1,2.6,1.7,0.94,0.11]
     band1 = 6
     band2 = 7
     band3 = 8
     band4 = 9
  endelse 

  print,'Inspecting ',strtrim(nsrc),' sources in: ',target + sourcelist,'.'

  ; Extinction correction
  if keyword_set(noukidss) then $
     filtkapremy=[2.5, 1.55, 1., 0.56, 0.43, 0.43, 0.43] else $
        filtkapremy=[2.5, 1.55, 1., 2.5, 1.55, 1., 0.56, 0.43, 0.43, 0.43] ;Note numbers should not REALLY be the same for UKIDSS and 2MASS filters, but these values are not actualy used by this routine.
  filtkapremy=filtkapremy*23.   ;scale to ISM at K
  kapvremy=220.

  a1av = filtkapremy[band1]/kapvremy
  a2av = filtkapremy[band2]/kapvremy
  a3av = filtkapremy[band3]/kapvremy
  a4av = filtkapremy[band4]/kapvremy  

;Loop through all sources, get magnitudes, and make cuts.
  galflags = intarr(nsrc)
  av = fltarr(nsrc)
  dav = fltarr(nsrc)
  
  if keyword_set(regionfile) then begin
     colorkey = ['orange','magenta','dodgerblue']
     groupkey = ['YSOc','galc','AGNc']
     openw,u,regionfile,/get_lun
  endif

     ngalc = 0
     nagnc = 0

  for i=0L,nsrc-1 do begin

     restore,files[i]
     valid = s.valid

;DEREDDEN
     weights = exp(-1*p.chi2/2.)
     av[i] = total(weights*p.av)/total(weights)
     dav[i] = sqrt(total(weights*(p.av - av[i])^2)/total(weights))
     
;Make starburst galaxy and AGN cuts
     case 1 of

;4-band cuts
        valid[band1] eq 1 and valid[band2] eq 1 and valid[band3] eq 1 and valid[band4] eq 1: begin 
           f1 = 1.d-3*s.f[band1]   ;convert mJy to Jy
           f2 = 1.d-3*s.f[band2]
           f3 = 1.d-3*s.f[band3]
           f4 = 1.d-3*s.f[band4]
           mag1 = -2.5*alog10(f1/zerof[band1]) - a1av*av[i]
           mag2 = -2.5*alog10(f2/zerof[band2]) - a2av*av[i]
           mag3 = -2.5*alog10(f3/zerof[band3]) - a3av*av[i]
           mag4 = -2.5*alog10(f4/zerof[band4]) - a4av*av[i]
           
;Find candidate starburst galaxies NOTE this could reject PAH
;YSOs if MAXAV is not set!
           if (mag2 - mag3 lt (1.05/1.2)*(mag3 - mag4 - 1) and $
               mag2 - mag3 lt 1.05 and $
               mag3 - mag4 gt 1 and $
               mag2 gt 11.5) OR $
              ((mag1 - mag3) lt (1.5/2.)*(mag2 - mag4 - 1) and $
               mag2 - mag3 lt 1.5 and $
               mag2 - mag4 gt 1 and $
               mag2  gt 11.5) $
           then begin
              galflags[i] = 1
              ngalc++
           endif else begin
;Separate candidate AGN from candidate YSOs
              if (mag2 - mag4 gt 0.5 and $
                  mag2 gt 13.5 + (mag2 - mag4 - 2.3)/0.4 and $
                  mag2 gt 13.5) AND $
                 (mag2 gt 14. + (mag2 - mag4 - 0.5) or $
                  mag2 gt 14.5 - (mag2 - mag4 - 1.2)/0.3 or $
                  mag2 gt 14.5) $
              then begin
                 galflags[i] = 2
                 nagnc++
              endif 
           endelse  
        end 

;2-band (4.5 and 8.0) AGN cut
        valid[band2] eq 1 and valid[band4] eq 1 and (valid[band1] ne 1 or valid[band3] ne 1): begin 
           f2 = 1.d-3*s.f[band2] ;convert mJy to Jy
           f4 = 1.d-3*s.f[band4] ;convert mJy to Jy
           mag2 = -2.5*alog10(f2/zerof[band2]) - a2av*av[i]
           mag4 = -2.5*alog10(f4/zerof[band4]) - a4av*av[i]

           if (mag2 - mag4 gt 0.5 and $
                  mag2 gt 13.5 + (mag2 - mag4 - 2.3)/0.4 and $
                  mag2 gt 13.5) AND $
                 (mag2 gt 14. + (mag2 - mag4 - 0.5) or $
                  mag2 gt 14.5 - (mag2 - mag4 - 1.2)/0.3 or $
                  mag2 gt 14.5) $
           then begin
              galflags[i] = 2
              nagnc++
           endif 
        end

;Single-band (3.6) faint, exgal cut (send to GALc list)
;        valid[band1] eq 1 and (valid[band2] ne 1 or valid[band3] ne 1 or valid[band4] ne 1): begin
        valid[band1] eq 1 and (valid[band2] ne 1 or valid[band4] ne 1): begin
           f1 = 1.d-3*s.f[band1] ;convert mJy to Jy
           mag1 = -2.5*alog10(f1/zerof[band1]) - a1av*av[i]
           
           if mag1 gt 14.5 then begin
              galflags[i] = 1
              ngalc++
           endif 
        end 

;The ELSE case should not happen if the sources have already been run through MALMCULL. Print warning message and add source to GALC list.
        else: begin 
           print,'WARNING: Source ' + s.desig + ' was not detected at [3.6]! Flagging it as Galc.'
           galflags[i] = 1
           ngalc++
        end 
     endcase

     if keyword_set(regionfile) then $
        printf,u,'galactic;circle('+strtrim(s.l,2)+','+strtrim(s.b,2)+ $
               ',6.5") # color='+colorkey[galflags[i]]+' tag={'+groupkey[galflags[i]]+'}'
  endfor 
 
  
  print,'Flagged ',nsrc - ngalc - nagnc,' candidate YSOs.'
  print,'Flagged ',ngalc,' candidate starburst (PAH) galaxies.'
  print,'Flagged ',nagnc,' candidate AGN.'


;Create regionfile, if desired
  if keyword_set(regionfile) then begin
     close,u
     free_lun,u
     print,'Wrote '+strtrim(nsrc,2)+' sources to '+regionfile+'.'
  endif ;Done making regionfile

end

