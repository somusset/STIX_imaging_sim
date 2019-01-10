FUNCTION stix_sim_calculate_visibilities, pixel_val, phase=phase

  ;+
  ; :project:
  ;   STIX imaging system simple simulation
  ;
  ; :description:
  ;   This function takes the values in four STIX pixels in 30 STIX collimators and return an array of 30 complex visibilities 
  ;   
  ; :inputs:
  ;   pixel_val: 30*4 array of pixel values
  ;   
  ; :keyword:
  ;   phase: default is 0, set to 1 to take into account the phase of the collimator  
  ;
  ; :Notes:
  ;   2018/07/30, Sophie Musset (University of MN), There is a remaining incertainty in the visibility calculation. It is mentionned in the 2016 paper
  ;                         that the phase direction of the collimator should be adapted in the phase factor...
  ;                         But when I do, the dirty image for a point source on axis is not the one that we expect.
  ;                         If I do not take the phase direction into consideration, it works fine, so that is the option retained here so far.
  ;                         Need more investigation.
  ;
  ; :history:
  ;   2018/07/30, Sophie Musset (UMN), initial release
  ;-

  ; initialisation of the option -- for testing purposes
  DEFAULT, phase, 0

  ;-----------------------------------------------------------------------------------------------
  ; initialisation
  ;-----------------------------------------------------------------------------------------------

  l = 0.22 ; largeur pixel in cm
  h = 0.46 ; longueur pixel in cm
  m1 = 4d/(!pi)^3*l*h*sin(!pi/4.)

  restore, 'stix_collimator_phases.sav'
  ; this restore the coll_phases variable

  vis_real_proxy = dblarr(30)
  vis_imgn_proxy = dblarr(30)
  
  ;-----------------------------------------------------------------------------------------------
  ; this calculates the values of C-A and D-B, where A, B, C, D are the values of the four pixels.
  ;-----------------------------------------------------------------------------------------------

  FOR k=0, 29 DO BEGIN
    vis_real_proxy[k] = pixel_val[k,2]-pixel_val[k,0] ; C-A
    vis_imgn_proxy[k] = pixel_val[k,3]-pixel_val[k,1] ; D-B
  ENDFOR

  ;-----------------------------------------------------------------------------------------------
  ; phase consideration -- for testing purposes
  ;-----------------------------------------------------------------------------------------------

  IF phase EQ 1 THEN vis_complex = complex(vis_real_proxy, vis_real_proxy)/(4*M1)*complex(cos(coll_phases*!pi/4.), sin(coll_phases*!pi/4.))$
  ELSE  vis_complex = complex(vis_real_proxy, vis_imgn_proxy)/(4*M1)*complex(cos(!pi/4.), sin(!pi/4.)) ; visibility in counts/s/cm2

  RETURN, vis_complex
END
