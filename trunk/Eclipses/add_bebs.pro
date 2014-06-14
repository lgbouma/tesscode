pro add_bebs, star, bkgnd, estruct, frac, rad, ph_p, mult, $ ;input
	bebdil,  $ ;output
  	aspix=aspix, radmax=radmax

  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  sz_frac = size(frac)
  npts = sz_frac[1]*sz_frac[2]*sz_frac[4]
  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
  if (keyword_set(radmax)) then radmax=radmax else radmax=8.0
  if (keyword_set(sq_deg)) then sq_deg=sq_deg else sq_deg=13.54
  ; Background catalog contains 0.134 sq degrees of stars
  ; Radius of 0.134 sq degree circle in pixels
  radas = sqrt(sq_deg/!dpi)*3600.
  radpix = radas/aspix

  AU_IN_RSUN = 215.093990942D0          ; in solar radii
  REARTH_IN_RSUN = 0.0091705248         ; in solar radii
  MSUN_IN_MEARTH = 332946D0
  AU_DAY_IN_CM_S = 173145684D0
  RV_AMP = 1731446.8 ; in m/s
  G_CM_S2 = 98.1 ; cm/s^2
 
;  if (keyword_set(fov)) then fov=fov else fov=24.0
;  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
;  ccd_pix = 4096.0
;  PIX_SCALE = fov*3600./ccd_pix

  nstars = n_elements(star)
  nbks = n_elements(bkgnd)
  nebtot = 0L

  pris = where(bkgnd.pri eq 1)
  secs = where(bkgnd.sec eq 1)

  npri = n_elements(pris)
  r1 = bkgnd[pris].r
  r2 = bkgnd[bkgnd[pris].companion.ind].r
  m1 = bkgnd[pris].m
  m2 = bkgnd[bkgnd[pris].companion.ind].m
  teff1 = bkgnd[pris].teff
  teff2 = bkgnd[bkgnd[pris].companion.ind].teff
  tmag1 = bkgnd[pris].mag.t
  tsys = bkgnd[pris].mag.tsys
  tmag2 = bkgnd[bkgnd[pris].companion.ind].mag.t
  ars = bkgnd[pris].companion.a*AU_IN_RSUN
  a   = bkgnd[pris].companion.a
  p = bkgnd[pris].companion.p
  cosi = -1.0 + 2.0*randomu(seed, npri)
  b1 = ars*cosi/r1
  b2 = ars*cosi/r2

  for ii=0,mult-1 do begin
    ; Re-randomize the inclinations 
    cosi = -1.0 + 2.0*randomu(seed, npri)
    b1 = a*cosi/r1
    b2 = a*cosi/r2
    ; Where are the (non-contact) eclipsing systems? 
    bin_ecl = where((abs(cosi) lt (r1+r2)/ars) and (ars gt (r1+r2)))
    
    if (bin_ecl[0] ne -1) then begin
      neb = n_elements(bin_ecl)
      pdur14 = dblarr(neb)
      pdur23 = dblarr(neb)
      sdur14 = dblarr(neb) 
      sdur23 = dblarr(neb) 
      gress1 = dblarr(neb) 
      gress2 = dblarr(neb) 
      dep1 = dblarr(neb) 
      dep2 = dblarr(neb) 
      dur1 = dblarr(neb) 
      dur2 = dblarr(neb) 
      a1 = dblarr(neb)
      a2 = dblarr(neb)
      r1 = r1[bin_ecl]
      m1 = m1[bin_ecl]
      teff1 = teff1[bin_ecl]
      tmag1 = tmag1[bin_ecl]
      r2 = r2[bin_ecl]
      m2 = m2[bin_ecl]
      teff2 = teff2[bin_ecl]
      tmag2 = tmag2[bin_ecl]
      tsys = tsys[bin_ecl]
      a = a[bin_ecl]
      ars = ars[bin_ecl]
      p = p[bin_ecl]
      b1 = b1[bin_ecl]
      b2 = b2[bin_ecl]
      cosi = cosi[bin_ecl]
 
      tot_ecl = where((abs(cosi) lt (r1-r2)/ars) and (ars gt (r1+r2)))
      gra_ecl = where((abs(cosi) ge (r1-r2)/ars) $
        and (abs(cosi) lt (r1+r2)/ars) and (ars gt (r1+r2)))

      ; Same for all eclipse types
      pdur14 = (p[bin_ecl]/!dpi)*asin((r1[bin_ecl]*sqrt(1.0-b1^2)+r2[bin_ecl])/a)
      sdur14 = (p[bin_ecl]/!dpi)*asin((r2[bin_ecl]*sqrt(1.0-b2^2)+r1[bin_ecl])/a)
      phr1 = phot_ratio(teff1, teff2, tmag1, tmag2, ph_p) ; Flux ratios
      phr2 = 1.0-phr1
      ; Only for total eclipses
      if (tot_ecl[0] ne -1) then begin
        pdur23[tot_ecl] = (p[tot_ecl]/!dpi)*asin((r1[tot_ecl]*sqrt(1.0-b1^2)+r2[tot_ecl])/a)
        sdur23[tot_ecl] = (p[tot_ecl]/!dpi)*asin((r2[tot_ecl]*sqrt(1.0-b2^2)+r1[tot_ecl])/a)
        dur1[tot_ecl] = (pdur14[tot_ecl] + pdur23[tot_ecl])/2. ; Trapezoidal area
        dur2[tot_ecl] = (sdur14[tot_ecl] + sdur23[tot_ecl])/2.
        a1[tot_ecl] = (r2[tot_ecl]/r1[tot_ecl])^2.
        a2[tot_ecl] = (r1[tot_ecl]/r2[tot_ecl])^2.
      end
      if (gra_ecl[0] ne -1) then begin
        dur1[gra_ecl] = pdur14[gra_ecl]/2. ; FWHM of "V" shaped eclipse
        dur2[gra_ecl] = sdur14[gra_ecl]/2.
        delt = (a[gra_ecl]*cosi[gra_ecl]) ; All defined on ph. 24-26 of Kopal (1979)
        ph1  = acos((delt^2. + r1[gra_ecl]^2. - r2[gra_ecl]^2.)/(2*r1[gra_ecl]*delt))
        ph2  = acos((delt^2. - r1[gra_ecl]^2. + r2[gra_ecl]^2.)/(2*r2[gra_ecl]*delt))
        da1  = r1[gra_ecl]^2.*(ph1-0.5*sin(2*ph1))
        da2  = r2[gra_ecl]^2.*(ph2-0.5*sin(2*ph2))
        a1[gra_ecl] = (da1+da2)/(!dpi*r1[gra_ecl]^2.)
        a2[gra_ecl] = (da1+da2)/(!dpi*r2[gra_ecl]^2.)
      end
      dep1 = phr1*a1
      toodeep = where(a1 gt 1.0)
      if (toodeep[0] ne -1) then dep1[toodeep] = phr1[toodeep]
      dep2 = phr2*a2
      toodeep = where(a2 gt 1.0)
      if (toodeep[0] ne -1) then dep2[toodeep] = phr2[toodeep]
      eclipse_hid = pris[bin_ecl]
      
      ; RV amplitude
      ;eclipse_k = RV_AMP*eclipse_p^(-1./3.) * eclipse_m * $ 
      ;	sqrt(1.0-star[allid].cosi^2.) * (star[allid].m)^(-2./3.) 
      gress1 = (pdur14-pdur23)/2.0
      gress2 = (sdur14-sdur23)/2.0
      ; Work out transit properties
      eclip = replicate({eclipstruct}, neb)
      eclip.class=3
      eclip.m1 = m1[bin_ecl]
      eclip.m2 = m2[bin_ecl]
      eclip.k = RV_AMP*2.0*!dpi*a*m2 * $ 
	sqrt(1.0-cosi^2.)/(p*m1)
      eclip.r1 = r1
      eclip.r2 = r2
      eclip.teff1 = teff1
      eclip.teff2 = teff2
      eclip.a = a
      ;planet_eclip.s = s
      eclip.p = p
      eclip.b = b1
      eclip.cosi = cosi
      eclip.dep1 = dep1
      eclip.dep2 = dep2
      eclip.dur1 = dur1
      eclip.dur2 = dur2
      eclip.tsys = tsys
      print, 'Created ', neb, ' eclipsing binaries out of ', n_elements(pris), ' primaries.'
      if (ii gt 0) then estruct=struct_append(estruct, eclip) $
	else estruct = eclip
      nebtot = nebtot + neb
    end ; if ebs
  end ; mult loop

  keep = lonarr(nbebtot)
  ; Blend with target stars
  for ii=0, nebtot-1 do begin
    randomp, r, 1., nstars, range_x=[0., radpix]
    ; How many of these are primaries and fall within radmax?
    gd = where((r lt radmax) and (star.sec ne 1))
    if (gd[0] ne -1) then begin
      ; ID the brightest star
      brightt = min(star[gd].mag.t, ind)
      estruct[ii].hostid = gd[ind]
      estruct[ii].sep= r[gd[ind]] ; in pixels
      keep[ii] = 1 ; set the keep flag
    end
  end
  if (total(keep) gt 0) then begin
    keepers = where(keep)  
    estruct = estruct[keepers] 
    ;bkteff = estruct.teff1
    ;bkmagt = tsys[keepers]
    ;starteff = star[gd].teff
    ;starmagt = star[gd].mag.t
      
    ;phr1 = phot_ratio(bkteff, starteff, bkmagt, starmagt, ph_p) ; Flux ratios
    ;phr2 = 1.0-phr1
    ;phr = phr1/phr2
      
   ; bkrad  = r[gd]

  endif  else delvar, estruct
end
