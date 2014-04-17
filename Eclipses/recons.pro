PRO recons, fstub
  fnames = file_search(fstub)
  numfil = n_elements(fnames)
  numstar = lonarr(numfil)
  for ii=0, numfil-1 do begin
    fits2sav, fnames[ii], nstar=nstar, dmax=0
    numstar[ii] = nstar
  end
  print, numfil, ' files contain ', total(numstar), ' stars within 10 pc.'
  nustar = replicate({starstruct}, total(numstar))
  idx0 = 0L
  for ii=0, numfil-1 do begin
    thisfn = repstr(fnames[ii], '.fits', '.sav')
    restore, thisfn
    idx = idx0+indgen(numstar[ii])
    nustar[idx] = star[where(star.coord.dm lt 0)]
    idx0 = idx0+numstar[ii]
  end
  save, nustar, filen='recons.sav'
END
