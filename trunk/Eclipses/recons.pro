PRO recons, fstub, fname, dmax=dmax, icmax=icmax
  fnames = file_search(fstub)
  restore, 'dartmouth_grid.sav'
  tt = mrdfits('tic_teff.fits');
  lfr = mrdfits('lfr.fits')
  numfil = n_elements(fnames)
  numstar = lonarr(numfil)
  for ii=0, numfil-1 do begin
    fits2sav, fnames[ii], ss, tt, jlfr=lfr, nstar=nstar, dmax=dmax, icmax=icmax, dbl=1
    numstar[ii] = nstar
  end
  print, numfil, ' files contain ', total(numstar), ' stars within ', 10.^(dmax/5.+1.), ' pc.'
;  print, numfil, ' files contain ', total(numstar), ' stars brighter than Ic=', icmax
  nustar = replicate({starstruct}, total(numstar))
  idx0 = 0L
  for ii=0, numfil-1 do begin
    if (numstar[ii] gt 0) then begin
      thisfn = repstr(fnames[ii], '.fits', '.sav')
      restore, thisfn
      idx = idx0+lindgen(numstar[ii])
      nustar[idx] = star[where(star.mag.dm le dmax)]
;      nustar[idx] = star[where(star.mag.ic le icmax)]
      idx0 = idx0+numstar[ii]
    end
  end
  star = nustar
  save, star, filen=fname
END
