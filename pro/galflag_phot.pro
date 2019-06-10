;Flag (1) candidate starburst galaxies and (2) AGN in a sample of YSOs
;that have been through the SED fitting procedure and processed
;through FITS2IDL_YSO.
;Strategy: deredden candidate YSOs using AV from the fits, set flags using simple color-mag criteria
;adapted from Gutermuth et al. (2008,2009).
;NOTE: This procedure can also be used on the output of stellar
;photosphere fits processed through FITS2IDL_STELLAR

pro galflag_phot, mags, av, galflags, regionfile=regionfile,ra=ra,dec=dec

;INPUT
;      MAGS       FLOAT[nsrc,4] --  IRAC 1234 mags of a source
;                                  sample. IMPORTANT: MISSING data
;                                  in MAGS must be assigned !VALUES.F_NAN.
;      AV                FLOAT[nsrc]   -- Visual extinction for sources, used in
;                                         dereddening.
;
;OUTPUT
;      GALFLAGS             INTEGER[nsrc] -- Source flags matching the order
;                                         of MAGS and AV.
;                                         0=candidate YSO,
;                                         1=candidate starburst galaxy
;                                         2= candidate AGN
;
;KEYWORD PARAMETERS
;     REGIONFILE=   'string' -- Specifies path to write a DS9 regionfile showing all sources
;                               color-coded by classification: orange
;                               = YSOc, magenta = galc, blue = agnc
;     RA=,DEC=    DOUBLE[nsrc] -- FK5 coordinates (decimal deg)
;                                    matching MAGS,
;                                    AV. ONLY required for making REGIONFILE.
  
;CALLS  
;
;VERSION HISTORY
;     Original GALFLAG - M. S. Povich         November 2011
;       Added DAV output parameter   May 2019  MSP
;
;PRODUCTION (v1.0)   MSP  9 May 2019
; * New! GALFLAG_PHOT, pass photometry and AV into routine
; rather than pulling from the save files and calculating internally.
  
if n_params() lt 3 then begin
   print,'syntax - galflag_phot, mags, av, galflags, regionfile=, ra=, dec='
   return
endif

  nsrc = n_elements(av)

  print,'GALFLAG_PHOT: Inspecting ',strtrim(nsrc),' sources.'

  ; Extinction correction
  filtkapremy = [0.56, 0.43, 0.43, 0.43]  ;IRAC bands, Indebetouw+2005
  filtkapremy=filtkapremy*23.   ;scale to ISM at K
  kapvremy=220.

  alamav = filtkapremy/kapvremy

;Loop through all sources, get magnitudes, and make cuts.
  galflags = intarr(nsrc)
  
  if keyword_set(regionfile) then begin
     if not keyword_set(ra) and not keyword_set(dec) then begin
        print,'WARNING! RA and DEC are required for REGIONFILE. RETURNING.'
        return
     endif 
     colorkey = ['orange','magenta','dodgerblue']
     groupkey = ['YSOc','galc','AGNc']
     openw,u,regionfile,/get_lun
  endif

     ngalc = 0
     nagnc = 0

  for i=0L,nsrc-1 do begin

;Make starburst galaxy and AGN cuts
     case 1 of
;4-band cuts
        finite(mags[i,0]) eq 1 and finite(mags[i,1]) eq 1 $
           and finite(mags[i,2]) eq 1 and finite(mags[i,3]) eq 1: $
           begin 
           mag1 = mags[i,0] - alamav[0]*av[i]
           mag2 = mags[i,1] - alamav[1]*av[i]
           mag3 = mags[i,2] - alamav[2]*av[i]
           mag4 = mags[i,3] - alamav[3]*av[i]
           
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
        finite(mags[i,1]) eq 1 and finite(mags[i,3]) eq 1 $
           and (finite(mags[i,0]) eq 0 or finite(mags[i,2]) eq 0): $
           begin 
           mag2 = mags[i,1] - alamav[1]*av[i]
           mag4 = mags[i,3] - alamav[3]*av[i]

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
        finite(mags[i,0]) eq 1 and $
           (finite(mags[i,1]) eq 0 or finite(mags[i,3]) eq 0): $
           begin
           mag1 = mags[i,0] - alamav[0]*av[i]
           
           if mag1 gt 14.5 then begin
              galflags[i] = 1
              ngalc++
           endif 
        end 

;The ELSE case should not happen if the sources have already been run through MALMCULL. Print warning message and add source to GALC list.
        else: begin 
;           print,'WARNING: Source ' + s.desig + ' was not detected at [3.6]! Flagging it as Galc.'
           print,'WARNING: Source not detected at [3.6]! Flagging it as Galc.'
           galflags[i] = 1
           ngalc++
        end 
     endcase

     if keyword_set(regionfile) then $
        printf,u,'fk5;circle('+strtrim(ra[i],2)+','+strtrim(dec[i],2)+ $
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

