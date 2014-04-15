pro add_ebs, sstruct=sstruct, pstruct=pstruct, infile=infile, outfile=outfile, verbose=verbose, fov=fov, $
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

  pris = where(star.pri eq 1)
  secs = where(star.sec eq 1)
  r1 = star[pris].r
  r2 = star[star[pris].companion.id].r
  m1 = star[pris].m
  m2 = star[star[pris].companion.id].m
  teff1 = star[pris].teff
  teff2 = star[star[pris].companion.id].teff
  a = star[pris].companion.a
  p = star[pris].companion.p
  cosi = star[pris].companion.cosi
  b1 = a*cosi/r1
  b2 = a*cosi/r2
   
  bin_ecl = (fabs(cosi) lt (r1+r2)/a)
  tot_ecl = (fabs(cosi) lt (r1-r2)/a)
  gra_ecl = (fabs(cosi) ge (r1-r2)/a) and (fabs(cosi) lt (r1+r2)/a)

  pdur14 = (p[bin_ecl]/!dpi)*asin((r1[bin_ecl]*sqrt(1.0-b1^2)+r2[bin_ecl])/a)
  pdur23 = (p[tot_ecl]/!dpi)*asin((r1[tot_ecl]*sqrt(1.0-b2^2)+r2[tot_ecl])/a)
  sdur14 = (p[bin_ecl]/!dpi)*asin((r2[bin_ecl]*sqrt(1.0-b2^2)+r1[bin_ecl])/a)
  sdur23 = (p[tot_ecl]/!dpi)*asin((r2[tot_ecl]*sqrt(1.0-b1^2)+r1[tot_ecl])/a)
  
  dur1 = pdur14/2.
  dur2 = sdur14/2.
  dur1[tot_ecl] = (pdur14[tot_ecl] + pdur23)/2.
  dur2[tot_ecl] = (sdur14[tot_ecl] + sdur23)/2.

  teff2phot, teff1, teff2, ph_p, phr1, phr2

  dep1 = dblarr(total(bin_ecl))
  dep2 = dblarr(total(bin_ecl))
  dep1[tot_ecl] = phr1*(r2[tot_ecl]/r1[tot_ecl])^2.
  dep2[tot_ecl] = phr2*(r1[tot_ecl]/r2[tot_ecl])^2.
  delt = (a[gra_ecl]*cosi[gra_ecl]) ; All defined on ph. 24-26 of Kopal (1979)
  ph1  = acos((delt^2. + r1[gra_ecl]^2. - r2[gra_ecl]^2.)/(2*r1[gra_ecl]*delt))
  ph2  = acos((delt^2. - r1[gra_ecl]^2. + r2[gra_ecl]^2.)/(2*r2[gra_ecl]*delt))
  da1  = r1[gra_ecl]^2.*(ph1-0.5*sin(2*ph1))
  da2  = r2[gra_ecl]^2.*(ph2-0.5*sin(2*ph2))
  dep1[gra_ecl] = phr1*(da1+da2)/(!dpi*r1^2.)
  dep2[gra_ecl] = phr2*(da1+da2)/(!dpi*r2^2.)

  planet_hid = lonarr(nplanets)
  idx0 = 0L
  hostids = hotstars[floor(double(nhotstars)*randomu(seed, binplanets))]
  idx = lindgen(binplanets) + idx0
  planet_per[idx] = periods
  planet_rad[idx] = radii
  planet_hid[idx] = hostids
  idx0 = max(idx) + 1
; Work out orbital distance and impact parameter
  allid = planet_hid
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
