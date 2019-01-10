PRO stix_sim_clean_iterations, dmap, nitermax, gain, weights, u, v, xx, yy, $
    final_residual, component_values, components_index, $
    plot_intermed=plot_intermed, stop_message=stop_message

;+
; :project:
;   STIX imaging system simple simulation
;
; :description:
;   This procedure carries the iterative process of the CLEAN algorithm, to find the CLEAN components and construct the residual map
;   
; :inputs:
;   dmap: the dirty map with which we start
;   nitermax: the max number of iterations
;   gain: clean gain
;   weights: the weights associated with the visibilities
;   u and v: (u,v) coordinates of the visibilities
;   xx, yy: the (x,y) coordinates of the dirty map
;   
; :outputs:
;   final_residual: the final residual map
;   component_values: clean components
;   component_index: the indices in x and y directions of the position of the clean components
;   
; :keywords:
;   plot_intermed, default is 1: set to 1 to plot the intermediary residual maps during the process
;   stop_message, out: string containing the reason for which the iteration process stopped
;   
; :history:
;   2018/07/30, Sophie Musset (University of Minnesota), initial release 
;     
;-

  ;---------------------------------------------------------
  ; initialisation
  ;---------------------------------------------------------
  DEFAULT, plot_intermed, 1
  res_map = dmap          ; first residual map
  iter=0                  ; first iteration
  peakflux = max(res_map) ; first peak flux
  Ck = list()             ; list of clean components
  xk = list()             ; list of x indices for the clean components
  yk = list()             ; list of y indices for the clean components

  ;---------------------------------------------------------
  ; plot first residual map
  ;---------------------------------------------------------
  IF plot_intermed EQ 1 THEN ires = image(res_map, xx, yy, margin=0.1, axis_style=1, rgb_table=5, title='current residual map')

  ;---------------------------------------------------------
  ; iterations
  ;---------------------------------------------------------
  WHILE iter LT nitermax AND peakflux GT 0 DO BEGIN
    ; clean process continues until number of iteration reaches a max
    ; or until the maximum peak flux is negative
    
    IF n_elements(where(res_map EQ peakflux)) GT 1 THEN PRINT, 'several max'
    
    ; find the indices of the position of the peak
    peakpos = array_indices(res_map, where(res_map EQ peakflux))
    xxk = peakpos[0]
    yyk = peakpos[1]
    ; add the position to the list
    xk.add, xxk
    yk.add, yyk
    ; calculate the clean component
    comp = gain*peakflux
    ck.add, comp
    ; calculate the new residual map
    new_res = res_map
    FOR i=0, n_elements(xx)-1 DO BEGIN
      FOR j=0, n_elements(yy)-1 DO BEGIN
        new_res[i,j] = res_map[i,j] - comp*total(weights*cos( (u*(xx[i]-xx[xxk]) + v*(yy[j]-yy[yyk]) ) ) )
      ENDFOR
    ENDFOR
    ; update the current residual map with the new residual map
    res_map = new_res
    ; update the plot
    IF plot_intermed EQ 1 THEN ires = image(res_map, xx, yy, margin=0.1, axis_style=1, rgb_table=5, title='current residual map',/current)

    ; find the value of the new peak 
    peakflux = max(res_map)
    ; update the iteration counter
    iter = iter+1
  ENDWHILE

  ;---------------------------------------------------------
  ; find the reason for end of iteration process
  ;---------------------------------------------------------
  print, iter, ' iterations to CLEAN process'
  IF iter LT nitermax THEN BEGIN
    IF peakflux LE 0 THEN stop_message='negative' $
    ELSE stop_message = 'inknown reason for CLEAN algo to stop'
  ENDIF ELSE BEGIN
    stop_message = 'reached max number of iterations'
  ENDELSE
  
  ;---------------------------------------------------------
  ; define outputs
  ;---------------------------------------------------------
  final_residual = res_map
  component_values = ck.toarray()
  components_index = [[xk.toarray()],[yk.toarray()]]
  
  ;---------------------------------------------------------
  ; close plot window
  ;---------------------------------------------------------
  IF plot_intermed EQ 1 THEN ires.close
END


FUNCTION stix_sim_clean_map_from_vs, visibilities, weights=weights, nitermax=nitermax, cleangain=cleangain, $
  npix=npix, pixsize=pixsize, xc=xc, yc=yc, xvalues=xvalues, yvalues=yvalues, loud=loud, noplot=noplot

  ;+
  ; :project:
  ;   STIX imaging system simple simulation
  ;
  ; :description:
  ;   This function takes 30 complexe visibilities and create a cleaned map
  ;
  ; :inputs:
  ;   visibilities: array of 30 complex numbers (visibilities), typically calculated with the function stix_sim_calculate_visibilities
  ;
  ; :keywords:
  ;   weights:  in, default is 30 identical values: array of 30 values containing the weights to apply to each visibility (collimator)
  ;   nitermax: in, default is 30: number of max iteration of the clean process
  ;   cleanbeam: in, default is 0.05: clean beam (fraction) - this default value herited from RHESSI (check this)
  ;   npix:     in, default is 64: number of pixels in each dimension of the image
  ;   pixsize:  in, default is 4: pixel size in arcsec
  ;   xc:       in, default is 0: center of the map in arcsec
  ;   yc:       in, default is 0: idem
  ;   loud:     in, default is 0: set to 1 to print some messages in the process
  ;   noplot:   in, default is 0: set to 1 to not plot the final images
  ;   xvalues, out,            : array containing the x values for each pixel
  ;   yvalues, out,            : array containing the y values for each pixel
  ;
  ; :calls:
  ;   stix_sim_dirty_map_from_vs
  ;   stix_sim_clean_iterations
  ;   
  ; :history:
  ;   2018/07/30, Sophie Musset (University of Minnesota), initial release
  ;-

  ;-------------------------------------------------------------
  ; define default values and restore useful variables
  ;-------------------------------------------------------------

  DEFAULT, weights, replicate(1./n_elements(visibilities), n_elements(visibilities))
  DEFAULT, nitermax, 30
  DEFAULT, cleangain, 0.05
  DEFAULT, npix, 64
  DEFAULT, pixsize, 4
  DEFAULT, xc, 0
  DEFAULT, yc, 0
  DEFAULT, loud, 0
  DEFAULT, noplot, 0

  restore, 'C:\Users\SMusset\Documents\GitHub\STIX\stix_visibility_uv.sav'
  ; this restore the u and v variables

  ;-------------------------------------------------------------
  ; create the dirty image as a starting point
  ;-------------------------------------------------------------
  dmap = stix_sim_dirty_map_from_vs(visibilities, weights=weights, npix=npix, pixsize=pixsize, xc=xc, yc=yc, xvalues=xvalues, yvalues=yvalues)

  ;-------------------------------------------------------------
  ; do the CLEAN iteration process
  ;-------------------------------------------------------------
  stix_sim_clean_iterations, dmap, nitermax, cleangain, weights, u, v, xvalues, yvalues, final_residual, component_values, components_index

  ;-------------------------------------------------------------
  ; define the resolution corresponding to each visibility
  ;-------------------------------------------------------------
  resolut = 0.5/sqrt(u^2 + v^2)

  ;-------------------------------------------------------------
  ; calculate beam size
  ;-------------------------------------------------------------
  sigma = 0.45*sqrt(total(weights)/total(weights/(resolut^2)))
  IF loud EQ 1 THEN print, 'beam ', sigma, ' arcsec'
  
  ;-------------------------------------------------------------
  ; Create the convolved map
  ;-------------------------------------------------------------
  convolved_map = final_residual*0d
  FOR i=0, n_elements(xvalues)-1 DO BEGIN
    FOR j=0, n_elements(yvalues)-1 DO BEGIN
      rk = sqrt( (double(xvalues[i])-double(xvalues[reform(components_index[*,0])]))^2 + (double(yvalues[j])-double(yvalues[reform(components_index[*,1])]))^2 )
      convolved_map[i,j] = total( component_values * exp( -0.5*( rk/sigma )^2) )
    ENDFOR
  ENDFOR
  
  ;-------------------------------------------------------------
  ; Create the clean map
  ;-------------------------------------------------------------
  clean_map = convolved_map + final_residual

  ;-------------------------------------------------------------
  ; If specified, plot the dirty image, the residual image,
  ; the convolved image, and the final clean image, together
  ;-------------------------------------------------------------

  IF noplot NE 1 THEN BEGIN
    imag = image( dmap, xvalues, yvalues, margin=0.15, rgb_table=5, title='dirty map', xtitle='arcsec', ytitle='arcsec', axis_style=1, dimensions=[1600,1400], layout=[2,2,1])
    imag = image( final_residual, xvalues, yvalues, margin=0.15, rgb_table=5, title='residual map', xtitle='arcsec', ytitle='arcsec', axis_style=1, dimensions=[1600,1400], layout=[2,2,2], /current)
    imag = image( convolved_map, xvalues, yvalues, margin=0.15, rgb_table=5, title='convolved components', xtitle='arcsec', ytitle='arcsec', axis_style=1, dimensions=[1600,1400], layout=[2,2,3], /current)
    imag = image( clean_map, xvalues, yvalues, margin=0.15, rgb_table=5, title='clean map', xtitle='arcsec', ytitle='arcsec', axis_style=1, dimensions=[1600,1400], layout=[2,2,4], /current)
  ENDIF

  RETURN, clean_map
END
