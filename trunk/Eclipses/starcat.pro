PRO starcat, fstub, fname, dmax=dmax, icmax=icmax, rmax=rmax, homies=homies
  fnames = file_search(fstub)
;  restore, 'dartmouth_grid.sav'
;  tt = mrdfits('tic_teff.fits');
;  lfr = mrdfits('lfr.fits')
  numfil = n_elements(fnames)
;  numstar = lonarr(numfil)
;  for ii=0, numfil-1 do begin
;    fits2sav, fnames[ii], ss, tt, jlfr=lfr, nstar=nstar, dmax=dmax, icmax=icmax, homies=homies, dbl=1
;    numstar[ii] = nstar
;  end
;  print, numfil, ' files contain ', total(numstar), ' stars within ', 10.^(dmax/5.+1.), ' pc.'
;  print, numfil, ' files contain ', total(numstar), ' stars brighter than Ic=', icmax
;  nustar = replicate({starstruct}, 1E8)
  idx0 = 0L
  nparam=21
  nustar = fltarr(1e8, nparam)
  for ii=0, numfil-1 do begin
    ;if (numstar[ii] gt 0) then begin
      ;thisfn = repstr(fnames[ii], '.fits', '.sav')
      restore, fnames[ii]
      if (keyword_set(rmax)) then gd = where(star.r le rmax and star.mag.ic le icmax)
      if (keyword_set(dmax)) then gd = where(star.mag.dm le dmax)
      if (keyword_set(homies)) then gd = where(star.spl)
      if (gd[0] ne -1) then begin
        m = star[gd].m
        rad = star[gd].r
        teff = star[gd].teff
        v = star[gd].mag.v
        r = star[gd].mag.r
        ic = star[gd].mag.ic
        z = star[gd].mag.z
        j = star[gd].mag.j
        h = star[gd].mag.h
        k = star[gd].mag.k
        dm = star[gd].mag.dm
        av = star[gd].mag.dm
        ps = star[gd].pri + 2*star[gd].sec
        mv = star[gd].mag.mv
        mic = star[gd].mag.mic
        mj = star[gd].mag.mj
        icsys = star[gd].mag.icsys
        jsys = star[gd].mag.jsys
        mvsys = star[gd].mag.mvsys
        micsys = star[gd].mag.micsys
        mjsys = star[gd].mag.mjsys

        poop = [[m],[rad],[teff],[v],[r],[ic],[z],[j],[h],[k],$
	[dm],[av],[ps],[mv],[mic],[mj],[icsys],[jsys],[mvsys],[micsys],[mjsys]]

        idx = idx0+lindgen(n_elements(gd))
        nustar[idx,*] = poop
        idx0 = idx0+n_elements(gd)
        print, fnames[ii], ' on index ', idx0
      end
  end
  nustar = nustar[lindgen(idx0-1),*]
  mwrfits, nustar, fname
END
