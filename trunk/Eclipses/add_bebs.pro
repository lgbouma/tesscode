ipro add_bebs, star, bkgnd, estruct, frac, rad, ph_p, mult, aspix=aspix, radmax=radmax

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
  a = bkgnd[pris].companion.a
  p = bkgnd[pris].companion.p
  ;cosi = bkgnd[pris].companion.cosi
 
  for ii=0,mult-1 do begin
    ; Re-randomize the inclinations 
    cosi = -1.0 + 2.0*randomu(seed, npri)
    b1 = a*cosi/r1
    b2 = a*cosi/r2
    ; Where are the eclipsing systems? 
    bin_ecl = where((abs(cosi) lt (r1+r2)/a) and (a gt (r1+r2)/AU_IN_RSUN))
    
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
      tot_ecl = where(abs(cosi) lt (r1-r2)/a)
      gra_ecl = where((abs(cosi) ge (r1-r2)/a) and (abs(cosi) lt (r1+r2)/a))
      ; Same for all eclipse types
      pdur14 = (p[bin_ecl]/!dpi)*asin((r1[bin_ecl]*sqrt(1.0-b1^2)+r2[bin_ecl])/a)
      sdur14 = (p[bin_ecl]/!dpi)*asin((r2[bin_ecl]*sqrt(1.0-b2^2)+r1[bin_ecl])/a)
      teff2phot, teff1, teff2, ph_p, phr1, phr2 ; Flux ratios
      ; Only for total eclipses
      if (tot_ecl[0] ne -1) then begin
        pdur23[tot_ecl] = (p[tot_ecl]/!dpi)*asin((r1[tot_ecl]*sqrt(1.0-b2^2)+r2[tot_ecl])/a)
        sdur23[tot_ecl] = (p[tot_ecl]/!dpi)*asin((r2[tot_ecl]*sqrt(1.0-b1^2)+r1[tot_ecl])/a)
        dur1[tot_ecl] = (pdur14[tot_ecl] + pdur23[tot_ecl])/2. ; Trapezoidal area
        dur2[tot_ecl] = (sdur14[tot_ecl] + sdur23[tot_ecl])/2.
        dep1[tot_ecl] = phr1*(r2[tot_ecl]/r1[tot_ecl])^2.
        dep2[tot_ecl] = phr2*(r1[tot_ecl]/r2[tot_ecl])^2.
      end
      if (gra_ecl[0] ne -1) then begin
        dur1[gra_ecl] = pdur14[gra_ecl]/2. ; FWHM of "V" shaped eclipse
        dur2[gra_ecl] = sdur14[gra_ecl]/2.
        delt = (a[gra_ecl]*cosi[gra_ecl]) ; All defined on ph. 24-26 of Kopal (1979)
        ph1  = acos((delt^2. + r1[gra_ecl]^2. - r2[gra_ecl]^2.)/(2*r1[gra_ecl]*delt))
        ph2  = acos((delt^2. - r1[gra_ecl]^2. + r2[gra_ecl]^2.)/(2*r2[gra_ecl]*delt))
        da1  = r1[gra_ecl]^2.*(ph1-0.5*sin(2*ph1))
        da2  = r2[gra_ecl]^2.*(ph2-0.5*sin(2*ph2))
        dep1[gra_ecl] = phr1*(da1+da2)/(!dpi*r1^2.)
        dep2[gra_ecl] = phr2*(da1+da2)/(!dpi*r2^2.)
      end
      toodeep = where(dep1 gt 1.0)
      if (toodeep[0] ne -1) then dep1[toodeep] = 1.0
      toodeep = where(dep2 gt 1.0)
      if (toodeep[0] ne -1) then dep2[toodeep] = 1.0
      ; RV amplitude
      ;eclipse_k = RV_AMP*eclipse_p^(-1./3.) * eclipse_m * $ 
      ;	sqrt(1.0-star[allid].cosi^2.) * (star[allid].m)^(-2./3.) 
      gress1 = (pdur14-pdur23)/2.0
      gress2 = (sdur14-sdur23)/2.0
      ; Work out transit properties
      eclip = replicate({eclipstruct}, neb)
      eclip.class=2
      eclip.m1 = m1[bin_ecl]
      eclip.m2 = m2[bin_ecl]
      eclip.k = RV_AMP*2.0*!dpi*a[bin_ecl]*m2[bin_ecl]* $ 
	sqrt(1.0-cosi[bin_ecl]^2.)/(p[bin_ecl]*m1[bin_ecl])
      eclip.r1 = r1[bin_ecl]
      eclip.r2 = r2[bin_ecl]
      eclip.teff1 = teff1[bin_ecl]
      eclip.teff2 = teff2[bin_ecl]
      eclip.a = a[bin_ecl]
      ;planet_eclip.s = s
      eclip.p = p[bin_ecl]
      eclip.b = b1[bin_ecl]
      eclip.cosi = cosi[bin_ecl]
      eclip.dep1 = dep1
      eclip.dep2 = dep2
      eclip.dur1 = dur1
      eclip.dur2 = dur2
      print, 'Created ', neb, ' eclipsing binaries out of ', n_elements(pris), ' primaries.'
      if (ii gt 0) then estruct=struct_append(estruct, eclip) $
	else estruct = eclip
      nebtot = nebtot + neb
    end ; if ebs
  end ; mult loop

  keep = intarr(nbebtot)
  ; Blend with target stars
  for ii=0, nebtot-1 do begin
    randomp, r, 1., nstars, range_x=[0., radpix]
    ; How many of these are primaries and fall within radmax?
    gd = where((r lt radmax) and (star.sec ne 1))
    if (gd[0] ne -1) then begin
      ; ID the brightest star
      brightt = min(star[gd].mag.t, ind)
      estruct[ii].hostid = gd[ind]
      keep[0] = 1 ; set the keep flag
      bkteff = bkgnds[gd].teff
      bkmagt = bkgnds[gd].mag.tsys
      recipteff = 4000./bkteff
      ph_filt = dblarr(nfilt, n_elements(gd))
      bk_frac = dblarr(nfilt, n_elements(gd))
      ; Compute the fluxes in each sub-filter for each star
      for jj=0, nfilt-1 do begin
        ph_filt[jj,*] = ph_p[jj,0] + ph_p[jj,1]*recipteff + $
           ph_p[jj,2]*recipteff^2. + ph_p[jj,3]*recipteff^3.
      end
      ph_filt[where(ph_filt lt 0.0)] = 0.0
      ; Find brightest band; this becomes the "center"
      ph_max = max(ph_filt, ph_ind, dimension=1)
      ph_ind = ph_ind mod (nfilt)
      ; Identify prf pixel for each neighbor (slow?)
      rand_ind = randomu(seed, n_elements(gd))
      for kk=0, n_elements(gd)-1 do begin
        ;thisprf = reform(frac[*,*,fov[ii],*,*])
        ;thisrad =  reform(rad[*,*,fov[ii],*,ph_ind[kk]])
        ;minrad = min(abs(thisrad-r[gd[kk]]), rind)
        ;bk_frac[*,kk] = thisprf(rind+lindgen(nfilt)*npts)      
        thisprf = reform(frac[*,*,fov[ii],*,*])
        thisrad = reform( rad[*,*,fov[ii],*,ph_ind[kk]])
        inrad = where(thisrad lt radmax)
        rad_ind = inrad[floor(n_elements(inrad)*rand_ind[kk])]+lindgen(nfilt)*npts
        bk_frac[*,kk] = thisprf[rad_ind]
      end
      phr1 = phot_ratio(teff1, teff2, tmag1, tmag2, ph_p) ; Flux ratios
      estruct.dep1 = estruct.dep1*phr1
      eclip[ii].bk_ph = eclip[ii].bk_ph + $
        total(10.^(-0.4*(bkmagt-10.))*total(bk_frac*ph_filt, 1))
    end
  end
  if (total(keep) gt 0) then estruct = estruct(where[keep]) $
    else delvar, estruct
end
