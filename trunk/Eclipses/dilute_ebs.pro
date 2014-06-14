PRO dilute_ebs, eclip, star, frac, rad, ph_p, dx, dy, dilvec, $
	aspix=aspix, radmax=radmax
  ; as a function of npix in the target star aperture, calculate the 
  ; diluting photon flux and the depth of eclipse

  ; How many filters in the prf?
  sz_ph_p = size(ph_p)
  nfilt = sz_ph_p[1]
  sz_frac = size(frac)
  npts = sz_frac[1]*sz_frac[2]*sz_frac[4]
  imgsize = sqrt(sz_frac[4])
  if (keyword_set(aspix)) then aspix=aspix else aspix=21.1
  if (keyword_set(radmax)) then radmax=radmax else radmax=6.
  ; Background catalog contains 0.134 sq degrees of stars
  ; Radius of 0.134 sq degree circle in pixels

  ; allocate the appropriately-sized dilution vector
  dilvec = dblarr(n_elements(eclip), imgsize*imgsize/4)
  depvec = dblarr(n_elements(eclip), imgsize*imgsize/4)

  hostid = eclip.hostid
  nbin = n_elements(eclip)
  print, 'Diluting ' , n_elements(binsys), ' binaries'
  bintmag = eclip.tsys
  binteff = eclip.teff1
  binsep  = eclip.sep ; in pixels
  bintheta = 2.*!dpi*randomu(seed, nbin)
  binx = binsep*cos(bintheta)
  biny = binsep*sin(bintheta)

  binfov = eclip[binsys].coord.fov_ind

  recipteff = 4000./binteff
  ph_filt = dblarr(nfilt, nbin)
  bin_frac = dblarr(nfilt, n_elements(binsys))
  ; Compute the fluxes in each sub-filter for each star
  for jj=0, nfilt-1 do begin
    ph_filt[jj,*] = ph_p[jj,0] + ph_p[jj,1]*recipteff + $
       ph_p[jj,2]*recipteff^2. + ph_p[jj,3]*recipteff^3.
  end
  ph_filt[where(ph_filt lt 0.0)] = 0.0
  dilpix = dblarr(imgsize/2, imgsize/2)
  for kk=0, nbin-1 do begin
    pixdx = floor(binx[kk] + bindx[kk]/10. + 0.1)
    pixdy = floor(biny[kk] + bindy[kk]/10. + 0.1)
    indx  = round(10*(binx[kk] - pixdx + bindx[kk]/10.))
    indy  = round(10*(biny[kk] - pixdy + bindy[kk]/10.))
    xsel = indgen(imgsize/2) + imgsize/4 + pixdx
    ysel = indgen(imgsize/2) + imgsize/4 + pixdy
    thisprf = reform(frac[indx,indy,binfov[kk],*,0:(nfilt-1)])
    thisimg = reform(thisprf#ph_filt[*,kk], imgsize, imgsize)
    thisimgx  = thisimg[xsel,*]
    thisimgxy = thisimgx[*,ysel]
    dilpix = 10^(-0.4*(bintmag[kk]-10.))*thisimgxy
    dilvec[binsys[kk],*] = reform(dilpix, imgsize*imgsize/4) 
  end
 end
END
