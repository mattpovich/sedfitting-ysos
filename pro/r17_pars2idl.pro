pro r17_pars2idl,infile,target_dir,data_parent=data_parent,filter=filter

;Read the (expicitly-formatted) ASCII output of the python sedfitter
;write_parameters program into IDL structures. Creates one structure for each source 
;and saves in an IDL save
;file in the specified directory (TARGET_DIR).
;INNOVATION for use with the R17 model sets: collects parameters for
;ALL model sets that fit a given source in one place, including the
;factor N_fits/N_models to give relative probability of each model
;set, and facilitating comparisons of chi2 - chi2_best across model
;                                            sets as recommended in
;Robitaille et al. (2017).
;Also creates a list of the output filenames for guiding further automated
;analysis routines.

;INPUTS
;         INFILE     strarr(n_sets) -- Path(s) to parameter file(s) to be
;                                read, one for each model set fit
;                           EXAMPLE: ['pars_spsi_g4.txt', 'pars_sphi_g4.txt']
;                           NOTE this file MUST be created with the option of
;                           to output ALL parameters (NOT default).
;         TARGET_DIR 'string' --  Path to the output directory for the IDL save
;                       files. EXAMPLE: 'results_ysoc' or 'results_xysoc'
;
;KEYWORD PARAMETERS
;        DATA_PARENT    'string' -- Path to a common fitter data format
;                                   file used for all fits. If this
;                                   keyword is set, the output save
;                                   files will contain a source S
;                                   structure in addition to the
;                                   parameter P structure.
;         FILTER=         'string' -- Filter the model fits saved
;                                     using the chi2 criteria defined
;                                     in the SED fitter
;                                     manual. Examples:
;                                     'F1','N10','E2'.
;                                    NOTE: These filters start counting from
;                                    the BEST chi2 found in ANY model
;                                    set and cut across different
;                                    model sets!
;
;
;MSP adapted from READ_GENERIC 1 August 2018
;

if n_params() lt 2 then begin
    print,'syntax: r17_pars2idl, infile, target_dir [, data_parent=, filter=]'
    return
endif


;create output directory
  if file_test(target_dir) then exists = 1
  if keyword_set(exists) then begin
     check = ''
     read,check,prompt=target_dir+': Directory exists. Overwrite? y/[n]  '
     if check ne 'y' then begin
        print,'Aborting to avoid overwriting files.'
        return
     endif else file_delete,target_dir,/recursive
  endif 
  spawn,'mkdir '+target_dir
  target = target_dir + '/'

;initialize ASCII I/O
 
  if keyword_set(data_parent) then begin
     n_datalines = file_lines(data_parent)
     datalines = strarr(n_datalines)
     openr,w,data_parent,/get_lun
     line = '' 
     for j=0L,n_datalines-1 do begin
        readf,w,line
        datalines[j] = line
        cells = line.Split(' ')
        cells = cells[where(cells ne '',n_cols)]
        if j eq 0 then begin  ;What's wrong with initializing arrays this way?
           n_bands = n_cols/3 - 1
           datanames = strarr(n_datalines)
           l = fltarr(n_datalines)
           b = fltarr(n_datalines)
           valid = intarr(n_datalines,n_bands)
           flux = fltarr(n_datalines,n_bands)
           flux_error = fltarr(n_datalines,n_bands)
        endif
        datanames[j] = cells[0]
        l[j] = float(cells[1])
        b[j] = float(cells[2])
        valid[j,*] = fix(cells[3:3+n_bands-1])
        flux[j,*] = float(cells[3+n_bands:*:2])
        flux_error[j,*] = float(cells[4+n_bands:*:2])
     endfor 
     close,w
     free_lun,w
  endif 

  n_source = 0L                  ;RUNNING counter of UNIQUE sources processed
  openw,v,target_dir+'/sourcelist.txt',/get_lun  ;print sourcelist

  for i_set=0,n_elements(infile)-1 do begin
     n_srcinset = 0L
     n_head = 3       ;number of header lines in parameter file for each source
     h0 = file_lines(infile[i_set]) - n_head
     h = h0
     head = strarr(n_head)
     headline = ''
     infoline = ''
     modelline = ''
     openr,u,infile[i_set],/get_lun

     print,'Processing parameter file: '+infile[i_set]

     for i=0, n_head-1 do begin
        readf,u,headline
        head[i] = headline
     endfor 
     colnames = head[1].Split(' ')
     colnames = colnames[where(colnames ne '',n_colnames)]
     
     while keyword_set(h) do begin

        readf,u,infoline
        fit_info = infoline.Split(' ')
        fit_info = fit_info[where(fit_info ne '',n_info)]
        if n_info ne 3 then message,'BUG Alert! Should be 3 items in the last header line for each source in parfile.'
      
        source_name = fit_info[0]
        ndata = fix(fit_info[1])
        n_fits = long(fit_info[2])

                                ;Construct source structure (OPTIONAL)
        if keyword_set(data_parent) then begin
           ind_data = where(datanames eq source_name,n_datamatch)
           if n_datamatch eq 1 then begin
              s = {DESIG:source_name,L:l[ind_data],B:b[ind_data],VALID:valid[ind_data,*],$
                   F:flux[ind_data,*],DF:flux_error[ind_data,*]}
           endif else begin
              print,'WARNING! No unique match to source '+source_name+' found in '+data_parent +'!'
              print,'No S structure was created.'
              stop
           endelse 
        endif 

                                ;Create structure table for model
                                ;parameter SUPERSET
      ;This is wasteful of disk space, but saves my sanity!
        f_nan    = !VALUES.F_NAN
        template_row = {$
                       MODEL_SET        :-1   ,$
                       MODEL_NAME       :''   ,$
                       CHI2             :f_nan,$
                       AV               :f_nan,$
                       SCALE            :f_nan,$
                       STAR             :replicate(f_nan,2),$
                       DISK             :replicate(f_nan,5),$
                       ENVELOPE         :replicate(f_nan,2),$
                       DISK_RMIN        :f_nan,$
                       ENVELOPE_RMIN    :f_nan,$
                       CAVITY           :replicate(f_nan,3),$  ;IMPLEMENTED LATER
                       AMBIENT          :replicate(f_nan,2),$
                       SCATTERING       :f_nan,$
                       INCLINATION      :f_nan $
                       }

        pars = replicate(template_row, n_fits)
      
        for j=0L,n_fits-1 do begin
           readf,u,modelline
           modelcells = modelline.Split(' ')
           modelcells = modelcells[where(modelcells ne '',n_modcells)]
           pars[j].model_name = modelcells[1]
           pars[j].chi2 = modelcells[2]
           pars[j].av = modelcells[3]
           pars[j].scale = modelcells[4]
           pars[j].star = modelcells[5:6]
           pars[j].disk = modelcells[7:11]
         
           pars[j].scattering = modelcells[n_modcells-2]
           pars[j].inclination = modelcells[n_modcells-1]

         ;Check for various model parameter sets and handle each in turn
           case 1 of
              n_modcells eq 14: begin
                 if j eq 0 and h eq h0 then print,'sp--s-i model parameter set detected'
                 pars[j].model_set = 1
              end 
              n_modcells eq 15: begin
                 if j eq 0 and h eq h0 then print,'sp--h-i model parameter set detected'
                 pars[j].model_set = 2
                 pars[j].disk_rmin = modelcells[12]
              end
              n_modcells eq 18: begin
                 if j eq 0 and h eq h0 then print,'spu-smi model parameter set detected'
                 pars[j].model_set = 14
                 pars[j].envelope = modelcells[12:13]
                 pars[j].ambient = modelcells[14:15]
              end
              n_modcells eq 20: begin
                 if j eq 0 and h eq h0 then print,'spu-hmi model parameter set detected'
                 pars[j].model_set = 15
                 pars[j].envelope = modelcells[12:13]
                 pars[j].disk_rmin = modelcells[14]
                 pars[j].envelope_rmin = modelcells[15]
                 pars[j].ambient = modelcells[16:17]
              end
              n_modcells eq 21: begin ;spubsmi models
                 if j eq 0 and h eq h0 then print,'spubsmi model parameter set detected'
                 pars[j].model_set = 16
                 pars[j].envelope = modelcells[12:13]
                 pars[j].cavity = modelcells[14:16]
                 pars[j].ambient = modelcells[17:18]
              end
              n_modcells eq 23: begin ;spubhmi models
                 if j eq 0 and h eq h0 then print,'spubhmi model parameter set detected'
                 pars[j].model_set = 17
                 pars[j].envelope = modelcells[12:13]
                 pars[j].cavity = modelcells[14:16]
                 pars[j].disk_rmin = modelcells[17]
                 pars[j].envelope_rmin = modelcells[18]
                 pars[j].ambient = modelcells[19:20]
              end
              else: begin
                 print,'WARNING: Unexpected number of parameter columns in file: '+infile[i_select]+'! Returning'
                 return
              end
           endcase
        endfor 


;NEW! Check if this source was already present in fits of a previous
;model set.
        if not file_test(target_dir + '/' + source_name + '.dat') then begin
           n_source++ 
           p = temporary(pars)
           printf,v,source_name + '.dat'
        endif else begin        ;Handle pre-existing sources...
           restore,target_dir + '/' + source_name + '.dat'
           p = [p,pars]
        endelse 

;FILTER fits by chisq, if desired. NOTE: This works ACROSS ALL MODEL
;SETS, as recommended by R17.
        if keyword_set(filter) then begin
           ftype = strmid(filter,0,1)
           fnum = float(strmid(filter,1))
           chi2 = p.chi2
           chi2_best = min(chi2)
           case 1 of 
              ftype eq 'N': begin
                 fnum = fix(fnum)
                 srtchi2 = sort(chi2)
                 ind_good = srtchi2[0:fnum-1]
                 n_fits2 = fnum
              end 
              ftype eq 'C': begin
                 ind_good = where(chi2 lt fnum,n_fits2)
              end 
              ftype eq 'D': begin
                 ind_good = where(chi2 - chi2_best lt fnum,n_fits2)
              end 
              ftype eq 'E': begin
                 ind_good = where(chi2/ndata lt fnum,n_fits2)
              end 
              ftype eq 'F': begin
                 ind_good = where((chi2 - chi2_best)/ndata lt fnum,n_fits2)
              end 
              else: print,'Invalid FILTER type specified (N C D E or F required)! SKIPPING filtering on Chisq.' 
           endcase
           p = p[ind_good]
        endif 
        
        if isa(s) then save,s,p,file=target_dir + '/' + source_name + '.dat' $
        else save,p,file=target_dir + '/' + source_name + '.dat'
        n_srcinset++
        h = h - n_fits - 1
      
        if n_srcinset mod 100 eq 0 then begin
           print,n_srcinset,' sources from ' + infile[i_set] + ' completed.'
           print,n_source,' UNIQUE sources processed.'
        endif 
     endwhile 

     close,u
     free_lun,u
  endfor                        ;i_set

  close,v
  free_lun,v
  print,'Successfully generated parameter files for',n_source,' sources.'
  
end
