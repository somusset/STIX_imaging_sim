FUNCTION calcul_stix_grid_param, pitch, orientation, phase

; pitch is the mean pitch used for the visibility, in mm
; orientation is the mean orientation for the visibility, in deg
; phase is the phase sens of the grid, is 1 or -1
; print, calcul_stix_grid_param(0.0380, 150., 1.)


DEFAULT, moire_width_mm, 8.8

du = 1./moire_width_mm/2.

deg2rad = !pi/180
rad2deg = 180/!pi

; calcul les composantes u,v de la visibilite
module = 1./pitch ; module de la visibilite
u = module*cos(orientation*deg2rad)
v = module*sin(orientation*deg2rad)

; calcul des composantes u pour la front and rear grids
uf = u - phase*du
ur = u + phase*du

modulef = sqrt( uf^2 + v^2)
moduler = sqrt( ur^2 + v^2)

pitchf = 1/modulef
pitchr = 1/moduler

orientationf = atan(v/uf)*rad2deg
orientationr = atan(v/ur)*rad2deg
IF orientationf LT 0 THEN orientationf=180.+orientationf
IF orientationr LT 0 THEN orientationr=180.+orientationr

RETURN, [pitchf, pitchr, orientationf, orientationr]

END