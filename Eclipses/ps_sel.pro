function ps_sel, tmag, teff, mass, rad, ph_p, minrad=minrad, per=per
  nstars = n_elements(tmag)
  if (keyword_set(minrad)) then minrad=minrad else minrad=2.0
  if (keyword_set(per)) then per=per else per=10.0
  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  ph_filt = dblarr(nfilt, nstars)
  ph_star = dblarr(nstars)
  geom_area = 69. ;74.6
  recipteff = 4000./teff
  for j=0, nfilt-1 do begin
    ph_filt[j,*] = ph_p[j,0] + ph_p[j,1]*recipteff + $
         ph_p[j,2]*recipteff^2. + ph_p[j,3]*recipteff^3.
  endfor
  ph_filt[where(ph_filt lt 0.0)] = 0.0
  ph_star = 1.0*10.^(-0.4*(tmag-10.))*total(ph_filt, 1)
  AU_IN_RSUN = 215.093990942
  REARTH_IN_RSUN = 0.0091705248
  a = mass^(1./3.) * (float(per)/365.25)^(2./3.)   ; in AU
  dur = rad * float(per) / (!DPI*a*AU_IN_RSUN)
  exptime = 2.*dur*24.*3600
  dep = (REARTH_IN_RSUN * float(minrad) / rad)^2.0
  sig = dep/7.0
  rn = 10.*sqrt(4.0*exptime/2.0)
  minphot = (1.+sqrt(1.+4.*sig^2.*rn^2.))/(2.*sig^2.)
  sel = where(ph_star gt (minphot/(exptime*geom_area)))
  return, sel
END
