PRO main
  fnums = mrdfits('fnums.fits')	; File containing healpix numbers
  eclass = [    1, $ ; Planets
                1, $ ; EBs
                1, $ ; BEBs
                1, $ ; HEBs
                0  ] ; BTPs
  tile_wrapper, '../../trilegal/', fnums, $
	'allpst5.fits', $		; output file name. Note: fits files don't over-write, they add extensions
  	n_trial=5, $		; number of trials (don't exceed 10)
	eclass=eclass, $	; from above
	ps_only=1, $		; 1=only run postage stamps, 0=run ffis as well
	detmag=0, $		; If you want the sim to return a magnitude-limited catalog, set this to the limit.
				;  0=run the TESS model for detection
        prf_file='psfs/dfrac_t75_f3p31.fits' ; PRF file produced by MATLAB
END
