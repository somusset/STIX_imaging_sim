FUNCTION stix_sim_dirty_map_from_vs, visibilities, weights=weights, npix=npix, pixsize=pixsize, xc=xc, yc=yc, xvalues=xvalues, yvalues=yvalues
  
  ;+
  ; :project:
  ;   STIX imaging system simple simulation
  ;
  ; :description:
  ;   This function takes 30 complexe visibilities and create a dirty map
  ;
  ; :inputs:
  ;   visibilities: array of 30 complex numbers (visibilities), typically calculated with the function stix_sim_calculate_visibilities
  ;   
  ; :keywords:
  ;   weights: in, default is 30 identical values: array of 30 values containing the weights to apply to each visibility (collimator)
  ;   npix:    in, default is 64: number of pixels in each dimension of the image
  ;   pixsize: in, default is 4: pixel size in arcsec
  ;   xc:      in, default is 0: center of the map in arcsec
  ;   yc:      in, default is 0: idem
  ;   xvalues, out,            : array containing the x values for each pixel
  ;   yvalues, out,            : array containing the y values for each pixel
  ;
  ; :history:
  ;   2018/07/30, Sophie Musset (University of Minnesota), initial release
  ;-
  
  ;-------------------------------------------------------------
  ; define default values for weight, several possibilities
  ;-------------------------------------------------------------
  
  ; DEFAULT, weights, replicate(1./n_elements(visibilities), n_elements(visibilities))
  ; DEFAULT, weights, exp( alog10([replicate(100.,3), replicate(10.*9./10.,3), replicate(10.*8./10.,3), $
  ;                    replicate(10.*7./10.,3), replicate(10.*6./10.,3), replicate(10.*5./10.,3), $
  ;                    replicate(10.*4./10.,3), replicate(10.*3./10.,3), replicate(10.*2./10.,3), $
  ;                    replicate(10.*1./10.,3)]) )
  
  pitches = [replicate(0.0380,3),replicate(0.0543,3),replicate(0.0777,3),$
             replicate(0.1112,3),replicate(0.1590,3),replicate(0.2275,3),$
             replicate(0.3254,3),replicate(0.4655,3),replicate(0.6659,3),$
             replicate(0.9526,3)]
  DEFAULT, weights, 1/pitches
  
  ;-------------------------------------------------------------
  ; define other default values and restore useful variables
  ;-------------------------------------------------------------
  
  DEFAULT, npix, 400
  DEFAULT, pixsize, 1
  DEFAULT, xc, 0
  DEFAULT, yc, 0

  restore, 'C:\Users\SMusset\Documents\GitHub\STIX\stix_visibility_uv.sav'
  ; this restore the u and v variables

  facto=2*!pi
  
  ;-------------------------------------------------------------
  ; calculate the x and y values of the image
  ;-------------------------------------------------------------

  xvalues = indgen(npix)*pixsize - (npix*pixsize/2 - xc)
  yvalues = indgen(npix)*pixsize - (npix*pixsize/2 - yc)

  ;-------------------------------------------------------------
  ; calculate the amplitude and phase of each visibitily
  ;-------------------------------------------------------------

  real = real_part(visibilities)
  imgn = imaginary(visibilities)
  amplitude = sqrt( real^2 + imgn^2 )
  phase = atan(visibilities, /phase)
  
  ;print, phase
  ;phase = dblarr(30)
  ;FOR k=0,29 DO BEGIN
  ;  phasee = stix_sim_complex_phase(visibilities[k])
  ;  phase[k] = phasee
  ;ENDFOR

  ;-------------------------------------------------------------
  ; construct the dirty image
  ;-------------------------------------------------------------

  ; initialisation
  dmap = dblarr(npix,npix)
  ; calculation of the value in each pixel
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap[i,j] = total(weights*amplitude*cos(facto*(u*xvalues[i]+v*yvalues[j])+phase))/total(weights)
    ENDFOR
  ENDFOR
  
  ;-------------------------------------------------------------
  ; do intermediate maps (for info plot)
  ;-------------------------------------------------------------  

  dmap1 = dblarr(npix,npix)
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap1[i,j] = total(weights[0]*amplitude[0]*cos(facto*(u[0]*xvalues[i]+v[0]*yvalues[j])+phase[0]))/total(weights[0])
    ENDFOR
  ENDFOR

  dmap2 = dblarr(npix,npix)
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap2[i,j] = total(weights[0:1]*amplitude[0:1]*cos(facto*(u[0:1]*xvalues[i]+v[0:1]*yvalues[j])+phase[0:1]))/total(weights[0:1])
    ENDFOR
  ENDFOR
  
  dmap4 = dblarr(npix,npix)
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap4[i,j] = total(weights[0:3]*amplitude[0:3]*cos(facto*(u[0:3]*xvalues[i]+v[0:3]*yvalues[j])+phase[0:3]))/total(weights[0:3])
    ENDFOR
  ENDFOR

  dmap8 = dblarr(npix,npix)
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap8[i,j] = total(weights[0:7]*amplitude[0:7]*cos(facto*(u[0:7]*xvalues[i]+v[0:7]*yvalues[j])+phase[0:7]))/total(weights[0:7])
    ENDFOR
  ENDFOR
  
  dmap16 = dblarr(npix,npix)
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap16[i,j] = total(weights[0:15]*amplitude[0:15]*cos(facto*(u[0:15]*xvalues[i]+v[0:15]*yvalues[j])+phase[0:15]))/total(weights[0:15])
    ENDFOR
  ENDFOR
  
  ;-------------------------------------------------------------
  ; for info, plot images with different number of visibilities
  ;-------------------------------------------------------------
  
  i=image(dmap1, xvalues, yvalues, margin=0.1, title='1 vis', layout=[5,1,1], axis_style=1, rgb_table=5, dimensions = [2100,500])
  i=image(dmap2, xvalues, yvalues, margin=0.1, title='2 vis', layout=[5,1,2], axis_style=1, rgb_table=5, /current)
  i=image(dmap4, xvalues, yvalues, margin=0.1, title='4 vis', layout=[5,1,3], axis_style=1, rgb_table=5, /current)
  i=image(dmap8, xvalues, yvalues, margin=0.1, title='8 vis', layout=[5,1,4], axis_style=1, rgb_table=5, /current)
  i=image(dmap16, xvalues, yvalues, margin=0.1, title='16 vis', layout=[5,1,5], axis_style=1, rgb_table=5, /current)

  ;-------------------------------------------------------------
  ; return the dirty image
  ;-------------------------------------------------------------
  
  RETURN, dmap
END