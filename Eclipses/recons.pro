PRO recons, fstub, dmax, fname
  fnames = file_search(fstub)
  numfil = n_elements(fnames)
  numstar = lonarr(numfil)
  for ii=0, numfil-1 do begin
    fits2sav, fnames[ii], nstar=nstar, dmax=dmax
    numstar[ii] = nstar
  end
  print, numfil, ' files contain ', total(numstar), ' stars within ', 10.^(dmax/5.), ' pc.'
  nustar = replicate({starstruct}, total(numstar))
  idx0 = 0L
  for ii=0, numfil-1 do begin
    if (numstar[ii] gt 0) then begin
      thisfn = repstr(fnames[ii], '.fits', '.sav')
      restore, thisfn
      idx = idx0+lindgen(numstar[ii])
      nustar[idx] = star[where(star.coord.dm le dmax)]
      idx0 = idx0+numstar[ii]
    end
  end
  star = nustar
  save, star, filen=fname
END