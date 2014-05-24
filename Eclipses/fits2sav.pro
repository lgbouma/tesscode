PRO fits2sav, fname, nstar=nstar, imax=imax, dmax=dmax, dbl=dbl

  dat = mrdfits(fname, 0, h, /SILENT)
  if keyword_set(dbl) then dat = [dat, dat]

  ; Abs. mag ranges for binary properties
  mp_min = [0.0, 0.1, 0.6, 0.8, 1.0, 1.4]  
  mp_max = [     0.1, 0.6, 0.8, 1.0, 1.4, 99.]  
  nmp = n_elements(mp_min)
  mf = [0.22, 0.26, 0.34, 0.41, 0.50, 0.75]
  abar = [4.5, 5.3, 20.,  45.,  45., 350]
  psig = [0.5, 1.3, 2.0,  2.3,  2.3, 3.0]
  qgam = [4.0, 0.4, 0.35, 0.3, 0.3, -0.5]
  homf = [0.0, 3.9, 3.8, 3.7, 3.7, 3.7]
  qnorm = 0.9*(qgam+1.0)*mf/(1.0-0.1^(qgam+1.0))
;  print, qnorm
;  readcol, fname, gc, logA, z, mini, logL, logT, logG, dm, av, $
;	comp, bol, t, j, h, ks, kp, g, r, i, z, dd, mnow 
  gc   = dat[*,0]
  logA = dat[*,1] 
  z    = dat[*,2] 
  mini = dat[*,3] 
  logL = dat[*,4] 
  logT = dat[*,5] 
  logG = dat[*,6] 
  dm   = dat[*,7] 
  av   = dat[*,8]
  comp = dat[*,9]
  bol  = dat[*,10]
  t    = dat[*,11] 
  j    = dat[*,12] 
  h    = dat[*,13] 
  ks   = dat[*,14] 
  kp   = dat[*,15] 
  g    = dat[*,16] 
  r    = dat[*,17] 
  i    = dat[*,18] 
  z    = dat[*,19] 
  dd   = dat[*,20] 
  mnow = dat[*,21] 
 
  pris = where(comp eq 1)
  secs = where(comp eq 2)
  npri = n_elements(pris)
  star = replicate({starstruct}, npri)
  
  m1 = mini[pris]
  m2 = mini[secs]
  q = m2/m1
  
  star.logage = logA[pris]
  star.feh = z[pris]
  star.mini = mini[pris]
  star.teff = 10.^(logT[pris])
  star.logg = logG[pris]
  star.mag.dm = dm[pris]
  star.mag.av = av[pris]
  star.mag.g = g[pris]
  star.mag.v = g[pris] - 0.5784*(g[pris]-r[pris]) - 0.0038 ;Lupton 2005
  star.mag.r = r[pris]
  star.mag.i = i[pris]
  star.mag.ic = i[pris] - 0.3780*(i[pris]-z[pris]) - 0.3974 ; Lupton 2005
  star.mag.z = z[pris]
  star.mag.j = j[pris]
  star.mag.h = h[pris]
  star.mag.k = ks[pris]
  star.mag.t = t[pris]
  star.mag.icsys = star.mag.ic
  star.mag.jsys = star.mag.j
  star.mag.tsys = star.mag.t
  star.m = mnow[pris]
  star.r = sqrt(star.m)/sqrt(10.^star.logg/27542.3)

  idx0 = long(npri)

  for ii=0,nmp-1 do begin
    sind = where((randomu(seed, npri) lt qnorm[ii]*q^qgam[ii]) and (m1 gt mp_min[ii]) and (m1 le mp_max[ii]))   
    if (sind[0] ne 0) then begin 
      nsec = n_elements(sind)
      bin_star = replicate({starstruct}, nsec)
      bin_star.logage = logA[secs[sind]]
      bin_star.feh = z[secs[sind]]
      bin_star.mini = mini[secs[sind]]
      bin_star.teff = 10.^(logT[secs[sind]])
      bin_star.logg = logG[secs[sind]]
      bin_star.mag.dm = dm[secs[sind]]
      bin_star.mag.av = av[secs[sind]]
      bin_star.mag.g = g[secs[sind]]
      bin_star.mag.v = g[secs[sind]] - 0.5784*(g[secs[sind]]-r[secs[sind]]) - 0.0038 ;Lupton 2005
      bin_star.mag.r = r[secs[sind]]
      bin_star.mag.i = i[secs[sind]]
      bin_star.mag.ic = i[secs[sind]] - 0.3780*(i[secs[sind]]-z[secs[sind]]) - 0.3974 ; Lupton 2005
      bin_star.mag.z = z[secs[sind]]
      bin_star.mag.j = j[secs[sind]]
      bin_star.mag.h = h[secs[sind]]
      bin_star.mag.k = ks[secs[sind]]
      bin_star.mag.t = t[secs[sind]]
      bin_star.m = mnow[secs[sind]]
      bin_star.r = sqrt(bin_star.m)/sqrt(10.^bin_star.logg/27542.3)
    
      bin_inds = lindgen(nsec)+idx0
      ;print, min(bin_inds), max(bin_inds)
      idx0 = idx0 + nsec
      ; Set the companionship flags
      star[sind].pri = 1
      bin_star.sec = 1
      ; Cross-reference the indices
      star[sind].companion.ind = bin_inds 
      bin_star.companion.ind = sind
      ; Cross-reference the indices
      star[sind].companion.m = bin_star.m 
      bin_star.companion.m = star[sind].m
     
      ; Add up the fluxes
      star[sind].mag.icsys = -alog10(10.^(-1.*bin_star.mag.ic) + 10.^(-1.*star[sind].mag.ic))
      bin_star.mag.icsys = star[sind].mag.icsys
      star[sind].mag.jsys = -alog10(10.^(-1.*bin_star.mag.j) + 10.^(-1.*star[sind].mag.j))
      bin_star.mag.jsys = star[sind].mag.jsys
      star[sind].mag.tsys = -alog10(10.^(-1.*bin_star.mag.t) + 10.^(-1.*star[sind].mag.t))
      bin_star.mag.tsys = star[sind].mag.tsys
      
      ; Convert mean separation into mean period
      logpbar = alog10(365.25*abar[ii]^(1.5)*(star[sind].m*(1.0+q[sind]))^(-0.5))
      ;print, median(logpbar)
      ; Randomize the period
      logp = randomn(seed, nsec)*psig[ii] + logpbar
      bin_star.companion.p = 10.^logp
      ; Cross-reference p and a
      star[sind].companion.p = bin_star.companion.p
      bin_star.companion.a = (star[sind].m*(1.0+q[sind]))^(1./3.)*(star[sind].companion.p/365.25)^(2./3.)
      star[sind].companion.a = bin_star.companion.a
      angseps = star[sind].companion.a/(10.*10.^(star[sind].mag.dm/5.))
      ;  Inclination and phase of binary
      cosi = -1.0 + 2.0*randomu(seed, nsec)
      phi = !DPI*2.0*randomu(seed, nsec)
      angseps = angseps*sqrt(sin(phi)^2. + cosi^2.*cos(phi)^2.)
      star[sind].companion.sep = angseps
      bin_star.companion.sep = angseps

      ; Mark stars for triples and quadruples
      if (homf[ii] gt 0) then begin
        tq = homf[ii]^(-2.) + homf[ii]^(-1.)
        tqind = where((randomu(seed, nsec) lt tq))
        if (tqind[0] ne -1) then begin
          ntq = n_elements(tqind)
          ; For quadruples, BOTH pri and sec get split
          qind = where((randomu(seed, ntq) lt homf[ii]^(-1)), complement=tind)
          if (qind[0] ne -1) then begin
            star[sind[tqind[qind]]].spl = 1
            bin_star[tqind[qind]].spl = 1
          endif
          ; For triples, pri OR sec gets split
          if (tind[0] ne -1) then begin
            nt = n_elements(tind)
            tpind = where((randomu(seed, nt) lt 0.5), complement=tsind)
            if (tpind[0] ne -1) then star[sind[tqind[tind[tpind]]]].spl = 1
            if (tsind[0] ne -1) then bin_star[tqind[tind[tsind]]].spl = 1
          endif
        endif
      endif 


      star = struct_append(star, bin_star)
      delvar, bin_star
    endif
  end
  if (keyword_set(dmax)) then begin
	near = where(star.mag.dm le dmax)
        if (near[0] eq -1) then nstar = 0 else nstar = n_elements(near)
  end else if (keyword_set(imax)) then begin
	near = where(star.mag.ic le imax)
        if (near[0] eq -1) then nstar = 0 else nstar = n_elements(near)
  end else nstar = n_elements(star)
  newfname = repstr(fname, 'fits', 'sav')
  save, star, filen=newfname
END
