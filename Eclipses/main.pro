PRO main
  fnums = mrdfits('../../trilegal/fnums.fits')
  ;fnums = 1000
  eclass = [    1, $ ; Planets
                0, $ ; EBs
                0, $ ; BEBs
                1  ] ; HEBs
  tile_wrapper, '../../trilegal/', fnums, 'planet1.fits', n_trial=1, eclass=eclass
END