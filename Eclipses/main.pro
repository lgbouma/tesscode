PRO main
  fnums = mrdfits('../../trilegal/fnums.fits')
  tile_wrapper, '../../trilegal/', fnums, 'ebnew.fits', n_trial=10
END