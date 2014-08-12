PRO main
  fnums = mrdfits('../../trilegal/fnums.fits')
  ;fnums = 1000
  eclass = [    1, $ ; Planets
                0, $ ; EBs
                0, $ ; BEBs
                0  ] ; HEBs
  tile_wrapper, '../../trilegal/', fnums, 'level1f3p33.fits', n_trial=5, eclass=eclass
END
