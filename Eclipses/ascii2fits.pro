PRO ascii2fits, fpre, fnums
  ;fnames = file_search(fstub)
  numfil = n_elements(fnums)
  numstar = intarr(numfil)
  for ii=0, numfil-1 do begin
    fname = fpre+string(fnums[ii], format='(I04)')+'.dat'
    readcol, fname, gc, logA, z, mini, logL, logT, logG, dm, av, $
        comp, bol, t, j, h, ks, kp, g, r, i, z, dd, mnow
    newtbl = [[gc], [logA], [z], [mini], [logL], [logT], [logG], [dm], [av], $
	[comp], [bol], [t], [j], [h], [ks], [kp], [g], [r], [i], [z], [dd], [mnow]]
    newname = repstr(fname, 'dat', 'fits')
    mwrfits, newtbl, newname
  end
END