PRO dilute_binary, eclip, star, frac, rad, ph_p, aspix=aspix
  ; How many filters in the prf?
  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  sz_frac = size(frac)
  npts = sz_frac[1]*sz_frac[2]*sz_frac[4]
  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
  ; Background catalog contains 0.134 sq degrees of stars
  ; Radius of 0.134 sq degree circle in pixels
  radmax = 8. ; only count flux up to 8 pixels away

  hostid = eclip.hostid
  binsys = where((star[hostid].pri or star[hostid].sec) and (star[hostid].companion.sep/aspix lt radmax))
  if (binsys[0] ne -1) then begin
    print, 'Diluting ' , n_elements(binsys), ' binaries'
    bintmag = star[star[hostid[binsys]].companion.ind].mag.t
    binteff = star[star[hostid[binsys]].companion.ind].teff
    binsep  = star[star[hostid[binsys]].companion.ind].companion.sep/aspix ; in pixels
    fov = eclip[binsys].coord.fov_ind

    recipteff = 4000./binteff
    ph_filt = dblarr(nfilt, n_elements(binsys))
    bin_frac = dblarr(nfilt, n_elements(binsys))
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
    for kk=0, n_elements(binsys)-1 do begin
      thisprf = reform(frac[*,*,fov[kk],*,*])
      thisrad =  reform(rad[*,*,fov[kk],*,ph_ind[kk]])
      minrad = min(abs(thisrad-binsep[kk]), rind)
      bin_frac[*,kk] = thisprf(rind+lindgen(nfilt)*npts)	
    end
    eclip[binsys].bk_ph = eclip[binsys].bk_ph + $
     10.^(-0.4*(bintmag-10.))*total(bin_frac*ph_filt, 1)
  end
END
