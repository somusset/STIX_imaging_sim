PRO stix_sim_pipeline, calculate_moire=calculate_moire, moire_file=moire_file, angles_sec=angles_sec, photon_flux=photon_flux, $
  row=row, npix=npix, pixsize=pixsize, nitermax=nitermax, stop=stop

  ;+
  ; :project:
  ;   STIX imaging system simple simulation 
  ;
  ; :description:
  ;   This procedure is a wrap-up that calculate the STIX pixel values, visibilities for a given incidence angle for a point source,
  ;   and create a dirty image and CLEAN it
  ;
  ; :keywords:
  ;   calculate_moire: in, default is 0, set to 1 to calculate the moire pattern and the values in the pixels
  ;   moire_file:      in, default is 'save_pixels_00.sav', name of the file that contain the moire/pixel information when we want to use a set of data already saved
  ;   angles_sec:      in, default is [0,0], angles of incidence on X and Y in arcsec
  ;   photon_flux      in, default is 2d6. This is the photon flux on detector before grids (flux in the detector will be about 1/4 of this input flux)
  ;   row:             in, default is 0. Set to 1 to select the bottom row of pixels, select the top row if not.
  ;   npix:            in, default is 400. This is the number of pixels in the resulting image.
  ;   pixsize:         in, default is 0.5. This is the pixel size in arcsec in the resulting image.
  ;   nitermax         in, default is 50. This is the maximum number of iterations in the clean process.
  ;   stop:            in, default is 0. Set to 1 to stop in the code
  ;
  ; :calls:
  ;   stix_sim_single_grid
  ;   stix_sim_calculate_visibilities
  ;   stix_sim_dirty_map_from_vs
  ;   stix_sim_clean_map_from_vs
  ; 
  ; :examples:
  ;   stix_sim_pipeline
  ;   stix_sim_pipeline, calculate_moire=1, angles_sec=[0., 5.]  
  ;   
  ; :history:
  ;   2019-01-10: Sophie Musset (University of Minnesota), initial release
  ; 
  ; :to be done:
  ;   2019-01-10: Sophie Musset (University of Minnesota), need to validate results, especially for sources with a non-zero incidence angle
  ;   
  ; :supporting documentation:
  ;   To know more details you will need to refer to a document called "STIX imaging simulation in IDL"
  ;-

  ;-----------------------------------------------------
  ; set default values
  ;-----------------------------------------------------

  DEFAULT, calculate_moire, 0
  DEFAULT, moire_file, 'save_pixels_00.sav'
  DEFAULT, angles_sec, [0.,0.] ; define incidence angles
  DEFAULT, photon_flux, 2d6
  DEFAULT, row, 0 ; can be 0 or 1 to select the top or the bottom row of pixels
  DEFAULT, npix, 400
  DEFAULT, pixsize, 0.5
  DEFAULT, nitermax, 50
  DEFAULT, stop, 0
  
  ;-----------------------------------------------------
  ; set hard-coded values
  ;-----------------------------------------------------

  l = 2.2 ; largeur pixel in mm
  h = 4.6 ; longueur pixel in mm
  eff_area = l*h*100 ; effective area in cm2

  ;-----------------------------------------------------
  ; read inputs
  ;-----------------------------------------------------

  ;define incidence angles in degree
  angles = angles_sec / 3600d ; in degree
  ; name angles to put in .sav name if created
  nameangles = strtrim(round(angles[0]),2)+'_'+strtrim(round(angles[1]),2)

  ;-----------------------------------------------------
  ; restore the STIX parameters
  ;-----------------------------------------------------
  
  restore, 'stix_grid_parameters.sav'
  ; this restore the grid_param and grid_names variables

  ;-----------------------------------------------------
  ; IF NEEDED calculate the pixel values for each STIX collimator
  ;-----------------------------------------------------

  IF calculate_moire EQ 1 THEN BEGIN
    pixels = dblarr(30,4)
    FOR k=0, 29 DO BEGIN
      r = stix_sim_single_grid(angles, reform(grid_param[k,*]), loud=1, /noplot, motif=motif, n_counts=photon_flux, plot_interm=0)
      pixels[k,*] = reform(r[*,row]) ; here we choose the first row of pixels
      ;print, 'flux ', total(r), ' counts/s'
    ENDFOR
    ;save it because it takes time
    save, pixels, filename='save_pixels_'+nameangles+'.sav'
    print, 'save_pixels_'+nameangles+'.sav', ' has been saved'
  ENDIF ELSE BEGIN
    restore, moire_file
  ENDELSE    
 
; THOSE LINES ARE TO PLOT GRID MOIRE PATTERN AND PIXEL INTENSITIES FOR TWO COLLIMATORS - THE PLOTS HAVE BEEN USED IN THE DOCUMENTATION    
;    r = stix_sim_single_grid(angles, reform(grid_param[7,*]), loud=1, motif=motif, n_counts=1d4, plot_interm=0)
;    r = stix_sim_single_grid(angles, reform(grid_param[22,*]), loud=1, motif=motif, n_counts=1d4, plot_interm=0)
  IF stop EQ 1 THEN stop

  ;-----------------------------------------------------
  ; calculate count flux
  ;-----------------------------------------------------

  flux = dblarr(30) ; will be in counts/s/cm2
  FOR k=0,29 DO flux[k] = total(reform(pixels[k,*])) /eff_area
  print, flux, 'counts/s/cm2'

  ;-----------------------------------------------------
  ; calculate the corresponding visibilities
  ;-----------------------------------------------------

  viz = stix_sim_calculate_visibilities(pixels)

  ;-----------------------------------------------------
  ; create dirty map and plot it
  ;-----------------------------------------------------
  
  dirty = stix_sim_dirty_map_from_vs(viz, npix=npix, pixsize=pixsize)
  i=image(dirty, rgb_table=5)
 
  IF stop EQ 1 THEN stop
  
  ;-----------------------------------------------------
  ; CLEAN the dirty map and plot results
  ;-----------------------------------------------------
  
  cleanm = stix_sim_clean_map_from_vs(viz, npix=npix, pixsize=pixsize, nitermax=nitermax)

  IF stop EQ 1 THEN stop

END
