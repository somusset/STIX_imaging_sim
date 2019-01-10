PRO stix_sim_calculate_grid_param, directory=directory, plotuvplane=plotuvplane

;+
; :project:
;   STIX imaging system simple simulation
;   
; :description:
;   This procedure calculate grid parameters for the 30 collimators of the STIX instrument. It produces 3 IDL sav files containing useful variables:
;   - 'stix_grid_parameters.sav' that contains variables 'grid_param', 'grid_names'. Grid param is an array of 30*4 values, containing for each collimator,
;   the pitch of the front grid in mm, the pitch of the rear grid in mm, the orientation of the front grid in degrees, the orientation of the rear grid in degrees.
;   grid_name is a string array of 30 values containing the names of the collimators (in the form '1a', etc)
;   - 'stix_visibility_uv.sav' that contains variables 'u' and 'v'. Those are the (u,v) values of each of the visibilities, in arcsec-1
;   - 'stix_collimator_phases.sav' that contains variable 'coll_phases', that contain the information about the phase direction of the collimator
; 
; :keyword:
;   directory: in, string containing the path to the directory in which to save the IDL save file produced
;   plotuvplane: in, set to 1 to plot the (u,v) plane coverage
;   
; :calls:
;   calcul_stix_grid_param
;   
; :example:
;   stix_sim_calculate_grid_param, dir='C:\Users\SMusset\Documents\GitHub\STIX\'
;   
; :history:
;   2018/07/30, SM (UMN), initial release
;-

;----------------------------------------------------------------------------------------
; Default values of keywords
;----------------------------------------------------------------------------------------

DEFAULT, directory, 'C:\Users\SMusset\Documents\GitHub\STIX\'
DEFAULT, plotuvplane, 0

;----------------------------------------------------------------------------------------
; declare STIX constants
;----------------------------------------------------------------------------------------

; STIX collimators have ten different pitches (corresponding to 10 different spatial resolutions)
pitches = [replicate(0.0380,3),replicate(0.0543,3),replicate(0.0777,3),$
           replicate(0.1112,3),replicate(0.1590,3),replicate(0.2275,3),$
           replicate(0.3254,3),replicate(0.4655,3),replicate(0.6659,3),$
           replicate(0.9526,3)]
           
; In order to reproduce the spiral response in the (u,v) plane, the orientation of the collimators are as follow:
angles = [150, 090, 030, 130, 070, 010, 110, 050, 170,$
          090, 030, 150, 070, 010, 130, 050, 170, 110,$
          030, 150, 090, 010, 130, 070, 170, 110, 050,$
          150, 090, 030]
          
; different phase orientation where introduced (even if not necessary).
;  They will be useful to estimate misalignment between grids and/or detectors
phases = [ 1,  1, -1, -1,  1, -1,  1,  1, -1, $
          -1, -1,  1, -1,  1,  1,  1, -1, -1, $
          -1,  1,  1, -1,  1, -1,  1, -1,  1, $
           1, -1, -1 ]

; collimator names: the number is linked to the pitch and the letter to the orientation
grid_names = ['1a', '1b', '1c', '2a', '2b', '2c', '3a', '3b', '3c', $
              '4a', '4b', '4c', '5a', '5b', '5c', '6a', '6b', '6c', $
              '7a', '7b', '7c', '8a', '8b', '8c', '9a', '9b', '9c', $
              '10a','10b','10c' ]

;----------------------------------------------------------------------------------------
; calculate grid parameters and save them in IDL save file
;----------------------------------------------------------------------------------------

; initialize a variable containing the grid parameters
grid_param = DBLARR(30,4) ; pitch of the front grid, pitch of the rear grid, orientation of the front grid, orientation of the rear grid         
; for each collimator, calculate the grid parameters by calling calcul_stix_grid_param with the pitch, orientation and phase for that collimator
FOR k=0,29 DO BEGIN
  res = calcul_stix_grid_param(pitches[k], angles[k], phases[k])
  grid_param[k,*] = res
ENDFOR

; save grid parameters and collimator name in an IDL save file
save, grid_param, grid_names, filename=directory+'stix_grid_parameters.sav'

;----------------------------------------------------------------------------------------
; calculate (u,v) coordinates of visibilities in arcsec-1 and save in IDL save file
;----------------------------------------------------------------------------------------

; declare constants
deg2rad = !pi/180d
rad2arcsec = 3600d*180d/!pi
D = 550d ; mm distance between grids

; calculate u and v
module = D/pitches ; amplitude of the visibilities in rad-1
u = module*cos(angles*deg2rad)/rad2arcsec ; u coordinate in arcsec-1
v = module*sin(angles*deg2rad)/rad2arcsec ; v coordinate in arcsec-1

; save in IDL save file
save, u, v, filename=directory+'stix_visibility_uv.sav'

;----------------------------------------------------------------------------------------
; save phase information in IDL save file
;----------------------------------------------------------------------------------------
coll_phases = phases
save, coll_phases, filename= directory+'stix_collimator_phases.sav'

;----------------------------------------------------------------------------------------
; plot the (u,v) plane coverage
;----------------------------------------------------------------------------------------
IF plotuvplane EQ 1 THEN s=scatterplot([u,-u],[v,-v], aspect_ratio=1, title='STIX (u,v) plane coverage', axis_style=2, dimensions=[900,900], sym_thick=2, sym_filled=1, sym_size=1.2, xtitle='arcsec-1', ytitle='arcsec-1')

END