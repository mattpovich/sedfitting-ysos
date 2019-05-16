;Flag candidate YSOs according to most probable SED stage: 1 = Stage
;0/I, 2 = Stage II/III, -1 = ambiguous

pro stageflag_r17, target_dir, stageflags, sourcelist=sourcelist, regionfile=regionfile, ambiguous=ambiguous, fk5=fk5

;INPUT
;      TARGET_DIR        'string' --  Directory containing IDL save files of
;                                    fit parameters and sourcelist
;                                    produced by FITS2IDL_YSO.
;
;OUTPUT
;      STAGEFLAGS             INTEGER[nsrc] -- Source flags matching the order
;                                         of SOURCELIST.
;                                         1=candidate protostar,
;                                         2=candidate disk-only
;                                         -1=ambiguous
;
;KEYWORDS
;     SOURCELIST=   'string' -- Name of file containing list of IDL save
;                                files to operate on (DEFAULT: 'sourcelist.txt')
;     REGIONFILE=   'string' -- Specifies path to write a DS9 regionfile showing all sources
;                               color-coded by classification: orange
;                               = YSOc, magenta = galc, blue = agnc
;     AMBIGUOUS=     FLOAT -- Probability (must be > 0.5) below which
;                             the Stage is considered ambiguous
;                             (DEFAULT = 0.67).
;     /FK5           SWITCH -- Set if the coordinates for the SED
;                              fitting results are celestial (default
;                              is Galactic). ONLY used in conjunction
;                              with REGIONFILE.

;CALLS
;
;VERSION HISTORY
;     Original - M. S. Povich         November 2011
;     ADAPTED to use the multiple model sets of Robitaille et
;     al. (2017), MSP August 2018

if n_params() lt 2 then begin
   print,'syntax - stageflag, target_dir, flags, sourcelist=, ambiguous=, regionfile='
   return
endif

  if not keyword_set(gas2dust) then gas2dust = 1.
  if keyword_set(ambiguous) then begin
     if ambiguous le 0.5 or ambiguous ge 1. then begin
        print,'Keyword AMBIGUOUS must have value between 0.5 and 1. SETTING TO ZERO, will have no effect!'
        ambiguous = 0
     endif 
  endif else ambiguous = 0.67
  
  if not keyword_set(sourcelist) then sourcelist = 'sourcelist.txt'
  target = target_dir + '/'
  readcol,target + sourcelist,format='(A)',names
  files = target + names
  nsrc = n_elements(files)

  print,'Finding stages for ',strtrim(nsrc),' sources in: ',target + sourcelist,'.'

 
;Loop through all sources and calculate Stages
  stageflags = intarr(nsrc)
  probstage = fltarr(2)
  
  if keyword_set(regionfile) then begin
     colorkey = ['goldenrod','red','yellow']
     groupkey = ['ambiguous','Stage 0/I','Stage II']
     if keyword_set(fk5) then coordkey = 'fk5' else coordkey = 'galactic'
     openw,u,regionfile,/get_lun
  endif

  for i=0L,nsrc-1 do begin

     restore,files[i]

     weights = exp(-1*p.chi2/2.) ;Tom doesn't like this, but I still don't see the problem using chi2 for RELATIVE probability!
     
     ;HARD-CODED for the 6 specific model sets I use! Tom DOES like this!
     ind_set1 = where(p.model_set eq 1,n_set1)
     probset1 = n_set1/1.e4
     if n_set1 ne 0 then weights[ind_set1] = weights[ind_set1]*probset1
     
     ind_set2 = where(p.model_set eq 2,n_set2)
     probset2 = n_set2/1.e4
     if n_set2 ne 0 then weights[ind_set2] = weights[ind_set2]*probset2

     ind_set14 = where(p.model_set eq 14,n_set14)
     probset14 = n_set14/1.e4
     if n_set14 ne 0 then weights[ind_set14] = weights[ind_set14]*probset14
     
     ind_set15 = where(p.model_set eq 15,n_set15)
     probset15 = n_set15/1.e4
     if n_set15 ne 0 then weights[ind_set15] = weights[ind_set15]*probset15

     ind_set16 = where(p.model_set eq 16,n_set16)
     probset16 = n_set16/4.e4
     if n_set2 ne 0 then weights[ind_set16] = weights[ind_set16]*probset16

     ind_set17 = where(p.model_set eq 17,n_set17)
     probset17 = n_set17/8.e4
     if n_set17 ne 0 then weights[ind_set17] = weights[ind_set17]*probset17

     ind1 = where(p.model_set gt 5,n1,complement=ind2,ncomplement=n2)
     if n1 ne 0 then probstage[0] = total(weights[ind1])
     if n2 ne 0 then probstage[1] = total(weights[ind2])

     probstage = probstage/total(probstage) ;normalize
     peak = max(probstage,stage_mostprob)
     if peak gt ambiguous then stageflags[i] = stage_mostprob + 1

     if keyword_set(regionfile) then $
        printf,u,coordkey+';circle('+strtrim(s.l,2)+','+strtrim(s.b,2)+ $
               ',6.5") # color='+colorkey[stageflags[i]]+' tag={'+groupkey[stageflags[i]]+'}'
  endfor 
 

;Reset flag values for ambiguous
  stageflags[where(stageflags eq 0)] = -1

  print,'Flagged ',strtrim(n_elements(where(stageflags eq -1)),2),' YSOs with ambiguous Stage.'
  print,'Flagged ',strtrim(n_elements(where(stageflags eq 1)),2),' YSOs as Stage 0/I'
  print,'Flagged ',strtrim(n_elements(where(stageflags eq 2)),2),' YSOs as Stage II/III'


;Create regionfile, if desired
  if keyword_set(regionfile) then begin
     close,u
     free_lun,u
     print,'Wrote '+strtrim(nsrc,2)+' sources to '+regionfile+'.'
  endif ;Done making regionfile


end

