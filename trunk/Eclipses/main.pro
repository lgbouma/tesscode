PRO main
  fnums = mrdfits('../../trilegal/fnums.fits')	; File containing healpix numbers
  eclass = [    1, $ ; Planets
                1, $ ; EBs
                1, $ ; BEBs
                1  ] ; HEBs
  tile_wrapper, '../../trilegal/', fnums, $
	'all12.fits', $		; output file name. Note: fits files don't over-write, they add extensions
  	n_trial=1, $		; number of trials (don't exceed 10)
	eclass=eclass, $	; from above
	ps_only=0, $		; 1=only run postage stamps, 0=run ffis as well
	detmag=12		; If you want the sim to return a magnitude-limited catalog, set this to the limit.
				;  0=run the TESS model for detection
END
