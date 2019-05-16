
;Plot a stamdard J-H versus H-K(s) color-color diagram, including
;loci of MS, red giant stars, and T Tauri stars
;plus reddening vectors

pro plot_nircc_rv,j,h,k,xr=xr,yr=yr,oplot=oplot,axes=axes,twomass=twomass

;INPUT
;    J,H,K          float[nsrc] -- three vectors containing the J, H,
;                                  and K magnitudes for a set of NSRC
;                                  sources. These must all be on the
;                                  same photometric system.
;KEYWORD PARAMETERS
;   XR=,YR=         float[2] -- Standard IDL XRANGE and YRANGE
;                               plotting range keywords.
;Set /OPLOT if you simply want to overplot annotations on an existing
;plot
;Set /AXES to make empty plot axes with appropriate labels.
;Set /twomass if plotting data taken using 2MASS K-short filter
;(OMIT FOR UKIDSS PHOTOMTERY). DO NOT MIX DIFFERENT PHOTOMETRIC
;SYSTEMS ON THE SAME PLOT!


;PRODUCTION VERSION (v1.0)   M. S. Povich 17 April 2019
;   * Implemented coyote graphics routines for better handling of plot
;   colors across platforms.
  
if n_params() lt 3 then begin
   print,'syntax -- plot_nircc_rv,j,h,k,xr=,yr=,/pub,/twomass'
   return
endif

  xthick = 2
  ythick = 2
  charthick = 2
  charsize = 1.6
  thick = 2
  if not keyword_set(xr) then xr = [-0.5,max(h-k)]
  
  if not keyword_set(oplot) then begin
     if not keyword_set(twomass) then xtit = 'H - K (mag)' else xtit = 'H - K!Ds!N (mag)'
     cgplot,h-k,j-h,ps=1,symsize=0.6,xr=xr,yr=yr,/xsty,/ysty $
            ,xtitle=xtit,ytitle='J - H (mag)' $
            ,xthick=xthick, ythick=ythick,charsize=charsize,th=thick,charth=charthick,/iso $
            ,nodata=keyword_set(axes),/window
  endif
   
  if not keyword_set(axes) then begin
;Intrinsic Colors -- Converted to 2MASS photometric system following
;                    Carpenter (2001, AJ, 121, 2851)
;Cool Giants:
      J_HG=[0.37,0.47,0.5,0.5,0.54,0.58,0.63,0.68,0.73,0.79,0.83,0.85,0.87,0.90,0.93,0.95,0.96,0.96]
      H_KG=[0.065,0.08,0.085,0.085,0.095,0.10,0.115,0.14,0.15,0.165,0.19,0.205,0.215,0.235,0.245,0.285,0.30,0.31]
      TwoMass_J_HG = 1.024*J_HG - 0.045
      TwoMass_H_KsG = 0.792*H_KG + 0.027
  
;Dwarves (Main Sequence) V:
      J_HMS=[-0.05,0.00,0.02,0.06,0.09,0.13,0.165,0.23,0.285,0.305,0.32,0.33,0.37,0.45,0.5,0.58,0.61,0.66,0.695,0.68,0.665,0.62,0.60,0.62,0.66]
      H_KMS=[-0.035,0.00,0.005,0.015,0.025,0.03,0.035,0.04,0.045,0.05,0.052,0.055,0.06,0.075,0.09,0.105,0.11,0.13,0.165,0.20,0.21,0.25,0.275,0.32,0.37]
      TwoMass_J_HMS = 1.024*J_HMS - 0.045
      TwoMass_H_KsMS = 0.792*H_KMS + 0.027
    
;Koorn Main Sequence Reddening Vector
      AV = float([0,5,10,15,20,25,30,35,40])
      pick = 0                  ;0 for hot end, 24 for cool end
      J_HRV=J_HMS[pick] + AV*(0.282-0.175)
      H_KRV=H_KMS[pick] + AV*(0.175-0.112)

;2Mass Main Sequence Reddening Vector
      TwoMass_J_HRV = TwoMass_J_HMS[pick] + 1.024*AV*(0.282-0.175)
      TwoMass_H_KsRV = TwoMass_H_KsMS[pick] + 0.792*AV*(0.175-0.112)

;This is the RV originating at the end of the giants line
      J_HRV2 = J_HG[16]+AV*(0.282-0.175)
      H_KRV2 = H_KG[16]+AV*(0.175-0.112)

;2MASS Giant Sequence Reddening Vectors
      TwoMass_J_HRV2 = TwoMass_J_HG[16] + 1.024*AV*(0.282-0.175)
      TwoMass_H_KsRV2 = TwoMass_H_KsG[16] + 0.792*AV*(0.175-0.112)

;CTTS locus
;From Meyer et. al. 1997 pg 290
      H_K_CTTS = [0.2,1.0]
      J_H_CTTS = 0.58*H_K_CTTS + 0.52

;2MASS CTTS Locus
      TwoMass_H_Ks_CTTS = 0.792*H_K_CTTS - 0.027
      TwoMass_J_H_CTTS = 0.58*TwoMass_H_Ks_CTTS + 0.52
   
;This is the RV originating at the end of the CTTS Locus
      J_H_RV3 = J_H_CTTS[1]+AV*(0.282-0.175)
      H_K_RV3 = H_K_CTTS[1]+AV*(0.175-0.112)

      TwoMass_J_H_RV3 = TwoMass_J_H_CTTS[1]+1.024*AV*(0.282-0.175)
      TwoMass_H_Ks_RV3 = TwoMass_H_Ks_CTTS[1]+0.792*AV*(0.175-0.112)

   ;Overplot loci and reddening vectors
      if not keyword_set(twomass) then begin
         cgoplot,H_KG,J_HG,linesty=0,color='orange',th=thick,/addcmd
         cgoplot,H_KMS,J_HMS,linesty=0,color='blue',th=thick,/addcmd
         cgoplot,H_K_CTTS,J_H_CTTS,color='brown',linestyle=2,th=thick,/addcmd
         cgoplot,H_KRV,J_HRV,color='blue',th=thick,ps=-1,linestyle=1,symsize=0.8,/addcmd
         cgoplot,H_KRV2,J_HRV2,color='orange',th=thick,ps=-1,linestyle=1,symsize=0.8,/addcmd
         cgoplot,H_K_RV3,J_H_RV3,color='brown',th=thick,ps=-1,linestyle=1,/addcmd
      endif else begin
         cgoplot,TwoMass_H_KsG,TwoMass_J_HG,linesty=0,color='orange',th=thick ,/addcmd
         cgoplot,TwoMass_H_KsMS,TwoMass_J_HMS,linesty=0,color='blue',th=thick,/addcmd
         cgoplot,TwoMass_H_Ks_CTTS,TwoMass_J_H_CTTS,color='brown',linestyle=2,th=thick,/addcmd
         cgoplot,TwoMass_H_KsRV,TwoMass_J_HRV,color='blue',th=thick,ps=-1,linestyle=1,symsize=0.8,/addcmd
         cgoplot,TwoMass_H_KsRV2,TwoMass_J_HRV2,color='orange',th=thick,ps=-1,linestyle=1,symsize=0.8,/addcmd
         cgoplot,TwoMass_H_Ks_RV3,TwoMass_J_H_RV3,color='brown',th=thick,ps=-1,linestyle=1,/addcmd
                                                                 
      endelse 

      cgtext,0.25,0.8,'A!DV!N interval (+) = 5 mag',charsize=charsize,co='blue',/norm,/addcmd

   endif  

END
