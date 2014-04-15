pro add_planets, sstruct=sstruct, pstruct=pstruct, infile=infile, outfile=outfile, verbose=verbose, fov=fov, $
  plotname=plotname

  AU_IN_RSUN = 215.093990942D0          ; in solar radii
  REARTH_IN_RSUN = 0.0091705248         ; in solar radii
  MSUN_IN_MEARTH = 332946D0
  AU_DAY_IN_CM_S = 173145684D0
  RV_AMP = 0.6395 ; in m/s
  G_CM_S2 = 98.1 ; cm/s^2
 
  if (keyword_set(fov)) then fov=fov else fov=24.0 
  ccd_pix = 4096.0
  PIX_SCALE = fov*3600./ccd_pix


  if (keyword_set(sstruct)) then star = sstruct else restore, infile

  nstars = n_elements(star)

  if(keyword_set(dressing)) then begin
    hotstars = where(star.teff ge 4000.)
    nhotstars = total(star.teff ge 4000.)
    coolstars = where(star.teff lt 4000.)
    ncoolstars = total(star.teff lt 4000.)
  endif else begin
    hotstars = lindgen(nstars)
    nhotstars = nstars
    coolstars = -1
    ncoolstars = 0
  endelse

  if (keyword_set(verbose)) then verbose=1 else verbose=0     
 
     period_boundary = [0.8, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 245.0, 418.0]
     radius_boundary = [0.8, 1.25, 2.0, 4.0, 6.0, 22.0]
     planet_type = ['Earths', 'Super-Earths', 'Small Neptunes', 'Large Neptunes', 'Giants']

     rate_fressin = dblarr(11,5) ; period bin, radius bin
     rate_fressin[0,*] = [0.18, 0.17, 0.035, 0.004, 0.015]
     rate_fressin[1,*] = [0.61, 0.74, 0.18,  0.006, 0.067]
     rate_fressin[2,*] = [1.72, 1.49, 0.73,  0.11,  0.17]
     rate_fressin[3,*] = [2.70, 2.90, 1.93,  0.091, 0.18]
     rate_fressin[4,*] = [2.70, 4.30, 3.67,  0.29,  0.27]
     rate_fressin[5,*] = [2.93, 4.49, 5.29,  0.32,  0.23]
     rate_fressin[6,*] = [4.08, 5.29, 6.45,  0.49,  0.35]
     rate_fressin[7,*] = [3.46, 3.66, 5.25,  0.66,  0.71]
     rate_fressin[8,*] = [0.0,  6.54, 4.31,  0.43,  1.25]
     rate_fressin[9,*] = [0.0,  0.0,  3.09,  0.53,  0.94]
     rate_fressin[10,*] =[0.0,  0.0,  0.0,   0.24,  1.05]
     rate_fressin = rate_fressin/100.

     rate_dressing = rate_fressin ;just for now

     nplanets = total(round(rate_fressin * nhotstars)) + $
                total(round(rate_dressing * ncoolstars))
     ; Pre-allocate for speed
     ;planet = replicate(template_planet, total(nplanets))
     planet_rad = dblarr(nplanets)
     planet_per = dblarr(nplanets)
     planet_hid = lonarr(nplanets)
     idx0 = 0L
     for i=10,0,-1 do begin
        for j=4,0,-1 do begin
           binplanets = round(rate_fressin[i,j] * nhotstars)
           if (binplanets gt 0) then begin              
	      ;tmp_planet = replicate(template_planet, binplanets)
              hostids = hotstars[floor(double(nhotstars)*randomu(seed, binplanets))]
              if(j eq 0) then radpow = 0.0 else radpow = -1.7
              randomp, periods, -1.0, binplanets, range_x = [period_boundary[i], period_boundary[i+1]], seed=seed
              randomp, radii, radpow, binplanets, range_x = [radius_boundary[j], radius_boundary[j+1]], seed=seed
              ;tmp_planet.r = radii
              ;tmp_planet.p = periods
              idx = lindgen(binplanets) + idx0
 	      planet_per[idx] = periods
  	      planet_rad[idx] = radii
	      planet_hid[idx] = hostids
              ;planet[idx] = tmp_planet
              ;delvar, tmp_planet
              idx0 = max(idx) + 1
           endif
        endfor
     endfor
     if (keyword_set(dressing)) then begin
       for i=10,0,-1 do begin
        for j=4,0,-1 do begin
           binplanets = round(rate_dressing[i,j] * ncoolstars)
           if (binplanets gt 0) then begin
              ;tmp_planet = replicate(template_planet, binplanets)
              hostids = coolstars[floor(double(ncoolstars)*randomu(seed, binplanets))]
              if(j eq 0) then radpow = 0.0 else radpow = -1.7
              randomp, periods, -1.0, binplanets, range_x = [period_boundary[i], period_boundary[i+1]], seed=seed
              randomp, radii, radpow, binplanets, range_x = [radius_boundary[j], radius_boundary[j+1]], seed=seed
              ;tmp_planet.r = radii
              ;tmp_planet.p = periods
              idx = lindgen(binplanets) + idx0
              planet_per[idx] = periods
	      planet_rad[idx] = radii
              planet_hid[idx] = hostids
              ;print, i, j, binplanets, n_elements(tmp_planet), n_elements(idx), max(idx)/float(total(nplanets))
              ;planet[idx] = tmp_planet
              ;delvar, tmp_planet
              idx0 = max(idx) + 1
           endif
        endfor
     endfor
   endif

; Star properties to be randomized
;  pris = where(star.sec ne 1)
  star.dx = floor(10.*randomu(seed, nstars))
  star.dy = floor(10.*randomu(seed, nstars))
;  close = where(star.sec and (star.companion.sep/pixsc le 0.1))
; Random orbital orientation
  if (keyword_set(req)) then begin
    star.cosi = 0.0 ;star[pla].r/(AU_IN_RSUN * star[pla].planet.a) * (-1.0 + 2.0*randomu(seed, nplanets));
  endif else begin
    star.cosi = -1.0 + 2.0*randomu(seed, nstars)
  endelse

; Work out orbital distance and impact parameter
  allid = planet_hid
  planet_a = (star[allid].m)^(1./3.) * (planet_per/365.25)^(2./3.); in AU
  planet_s = (star[allid].r)^2.0 * (star[allid].teff/5777.0)^4. / (planet_a)^2. ; indicent flux wrt sun-earth value
; Equilibrium Temp.
  planet_teq = (star[allid].teff/5777.0)*sqrt(star[allid].r/(planet_a*AU_IN_RSUN))
; Impact parameter
  planet_b = (planet_a*AU_IN_RSUN / star[allid].r) * star[allid].cosi; assumes circular orbit
 

;  planet.m = 10.55*4.*!dpi/3.*(planet.r/3.9)^3. * $ ; assume MgSiO_3 composition (Seager 2007)
;     (1.0 + (1.0 - 3./5.*0.541)*(2./3.*!dpi*(planet.r/3.9)^2.)^0.541)
  ; Weiss & Marcy 2014:
  planet_m = dblarr(nplanets)
  lomass = where(planet_r lt 1.5)
  planet_m[lomass] = 0.440*(planet_r[lomass])^3. + 0.614*(planet_r[lomass])^4.
  himass = where(planet_r ge 1.5)
  planet_m[himass] = 2.69*(planet_r[himass])^0.93
; RV amplitude
  planet_k = RV_AMP*planet_p^(-1./3.) * planet_m * $ 
	sqrt(1.0-star[allid].cosi^2.) * (star[allid].m)^(-2./3.) 

	;2.0*!dpi*sqrt(1.0-star[pla].cosi^2.)*star[pla].planet.a * AU_DAY_IN_CM_S * planet_mass /  $
	;	(star[pla].planet.p * star[pla].m * MSUN_IN_MEARTH)
; Work out transit properties
  tra = where(abs(planet_b) lt 1.0)
  traid = planet_hid[tra]
  ntra = n_elements(tra)
  planet_eclip = replicate({structeclip}, ntra)
  planet_eclip.class=1
  planet_eclip.m1 = star[traid].m
  planet_eclip.m2 = planet_m[tra]
  planet_eclip.k = planet_k[tra]
  planet_eclip.r1 = star[traid].r
  planet_eclip.r2 = planet_rad[tra]*REARTH_IN_RSUN
  planet_eclip.teff1 = star[traid].teff
  planet_eclip.teff2 = planet_teq[tra]
  planet_eclip.a = planet_a[tra]
  planet_eclip.s = planet_s[tra]
  planet_eclip.p = planet_per[tra]
  planet_eclip.b = planet_b[tra]
  planet_eclip.hostid = planet_hid[tra]
  planet_eclip.dep1 = (planet_eclip.r2 / planet_eclip.r1 )^2.0
  planet_eclip.dep2 = (planet_eclip.teff2/planet_eclip.teff1)*(planet_eclip.r2/planet_eclip.r1 )^2.0
  planet_eclip.dur1 = planet_eclip.r1 * planet_eclip.p * sqrt(1.-(planet_eclip.b)^2.) / (!PI*planet_eclip.a*AU_IN_RSUN)
 ; planet_eclip.durpar = planet[tra].r * $
	REARTH_IN_RSUN * planet[tra].p / $
        sqrt(1.-(planet[tra].b)^2.) / $
        (!PI*planet[tra].a*AU_IN_RSUN)
  
  print, 'Created ', n_elements(tra), ' transiting planets out of ', nplanets, ' total.'

  if (keyword_set(outfile)) then save,filen=outfile, planet

  pstruct=planet
end
