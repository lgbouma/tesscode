PRO fits2sav_wrapper, fstub
  fnames = file_search(fstub)
  restore, 'dartmouth_grid.sav'
  tt = mrdfits('tic_teff.fits');
  lfr = mrdfits('lfr.fits')
  numfil = n_elements(fnames)
  for ii=0, numfil-1 do begin
    if (strmatch(fstub, '*hp*')) then fits2sav, fnames[ii], ss, tt, jlfr=lfr, nstar=nstar, dbl=1
    if (strmatch(fstub, '*dp*')) then fits2sav, fnames[ii], ss, tt, jlfr=lfr, nstar=nstar, kmin=15.
    if (strmatch(fstub, '*bk*')) then fits2sav, fnames[ii], ss, tt, jlfr=lfr, nstar=nstar, tmin=21.
  end
END
