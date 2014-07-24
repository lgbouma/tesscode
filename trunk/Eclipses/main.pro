PRO main
  fnums = mrdfits('../../trilegal/fnums.fits')
  ;fnums = 1000
  eclass = [    1, $ ; Planets
                1, $ ; EBs
                1, $ ; BEBs
                1  ] ; HEBs
  tile_wrapper, '../../trilegal/', fnums, 'run4.fits', n_trial=1, eclass=eclass
END
