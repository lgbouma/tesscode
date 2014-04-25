PRO fits2sav, fname, nstar=nstar, dmax=dmax

  dat = mrdfits(fname, 0, h, /SILENT)
  dat = [dat, dat]

  ; Abs. mag ranges for binary properties
  mp_min = [0.0, 0.1, 0.6, 0.8, 1.0, 1.4]  
  mp_max = [     0.1, 0.6, 0.8, 1.0, 1.4, 99.]  
  nmp = n_elements(mp_min)
  mf = [0.22, 0.26, 0.34, 0.41, 0.50, 0.75]
  abar = [4.5, 5.3, 20.,  45.,  45., 350]
  psig = [0.5, 1.3, 2.0,  2.3,  2.3, 3.0]
  qgam = [4.0, 0.4, 0.35, 0.3, 0.3, -0.5]
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
  q = mini[secs]/mini[pris]
  
  star.logage = logA[pris]
  star.feh = z[pris]
  star.mini = mini[pris]
  star.teff = 10.^(logT[pris])
  star.logg = logG[pris]
  star.coord.dm = dm[pris]
  star.mag.av = av[pris]
  star.mag.g = g[pris]
  star.mag.r = r[pris]
  star.mag.i = i[pris]
  star.mag.z = z[pris]
  star.mag.j = j[pris]
  star.mag.h = h[pris]
  star.mag.k = ks[pris]
  star.mag.t = t[pris]
  star.m = mnow[pris]
  star.r = sqrt(star.m)/sqrt(10.^star.logg/27542.3)

  idx0 = npri

  for ii=0,nmp-1 do begin
    sind = where((randomu(seed, npri) lt qnorm[ii]*q^qgam[ii]) and (m1 gt mp_min[ii]) and (m1 le mp_max[ii]))    
    nsec = n_elements(sind)
    bin_star = replicate({starstruct}, nsec)
    bin_star.logage = logA[secs[sind]]
    bin_star.feh = z[secs[sind]]
    bin_star.mini = mini[secs[sind]]
    bin_star.teff = 10.^(logT[secs[sind]])
    bin_star.logg = logG[secs[sind]]
    bin_star.coord.dm = dm[secs[sind]]
    bin_star.mag.av = av[secs[sind]]
    bin_star.mag.g = g[secs[sind]]
    bin_star.mag.r = r[secs[sind]]
    bin_star.mag.i = i[secs[sind]]
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
    angseps = star[sind].companion.a/(10.*10.^(star[sind].coord.dm/5.))
    ; Inclination and phase of binary
    cosi = -1.0 + 2.0*randomu(seed, nsec)
    phi = !DPI*2.0*randomu(seed, nsec)
    angseps = angseps*sqrt(sin(phi)^2. + cosi^2.*cos(phi)^2.)
    star[sind].companion.sep = angseps
    bin_star.companion.sep = angseps

    star = struct_append(star, bin_star)
    delvar, bin_star

  end
  if (keyword_set(dmax)) then begin
	near = where(star.coord.dm le dmax)
        if (near[0] eq -1) then nstar = 0 else nstar = n_elements(near)
  endif else begin 
	nstar = n_elements(star)
  endelse
  newfname = repstr(fname, 'fits', 'sav')
  save, star, filen=newfname
END
