PRO fits2sav_wrapper, fstub
  fnames = file_search(fstub)
  restore, 'dartmouth_grid.sav'
  tt = mrdfits('tic_teff.fits');
  lfr = mrdfits('lfr.fits')
  numfil = n_elements(fnames)
  for ii=0, numfil-1 do begin
    if (strmatch(fstub, '*hp*')) then fits2sav, fnames[ii], ss, tt, nstar=nstar, dbl=1, dartcor=0, jlfr=0, suffix='raw'
    if (strmatch(fstub, '*dp*')) then fits2sav, fnames[ii], ss, tt, nstar=nstar, kmin=15., dartcor=0, jlfr=0, suffix='raw'
    if (strmatch(fstub, '*bk*')) then fits2sav, fnames[ii], ss, tt, nstar=nstar, tmin=21., dartcor=0, jlfr=0
  end
END
