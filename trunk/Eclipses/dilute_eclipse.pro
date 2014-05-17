PRO dilute_eclipse, eclip, bkgnds, frac, rad, ph_p, aspix=aspix, sq_deg=sq_deg
  ; How many filters in the prf?
  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  sz_frac = size(frac)
  npts = sz_frac[1]*sz_frac[2]*sz_frac[4]
  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
  if (keyword_set(sq_deg)) then sq_deg=sq_deg else sq_deg=1.0
  ; Background catalog contains 0.134 sq degrees of stars
  ; Radius of 0.134 sq degree circle in pixels
  radas = sqrt(sq_deg/!dpi)*3600.
  radpix = radas/aspix
  print, 'Considering stars out to ', radas, ' arcsec or ', radpix, ' pixels'
  nbk = n_elements(bkgnds)
  radmax = 8. ; only count flux up to 8 pixels away
  
  neclip = n_elements(eclip)
  fov = eclip.coord.fov_ind

  for ii=0, neclip-1 do begin
    randomp, r, 1., nbk, range_x=[0., radpix]
    ; How many of these are primaries and fall within radmax?
    gd = where((r lt radmax) and (bkgnds.sec ne 1))
    if (gd[0] ne -1) then begin
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
      for kk=0, n_elements(gd)-1 do begin
        thisprf = reform(frac[*,*,fov[ii],*,*])
        thisrad =  reform(rad[*,*,fov[ii],*,ph_ind[kk]])
        minrad = min(abs(thisrad-r[gd[kk]]), rind)
        bk_frac[*,kk] = thisprf(rind+lindgen(nfilt)*npts)	
      end
      eclip[ii].dil = eclip[ii].dil + $
	total(10.^(-0.4*(bkmagt-10.))*total(bk_frac*ph_filt, 1))
    end
  end
END