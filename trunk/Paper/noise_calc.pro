PRO noise_calc 
  imag = findgen(151)/10.+5.
  ;print, imag
  npix_min=3
  npix_max=49
  geom_area = 63.9
  teff=3400.
  fov = 24.
  pixsc = fov*3600./4096.
  field_angle = fov/2.
  frac_fits = mrdfits('../ExpTimeCalc/frac24_105_3400k.fits')
  bk_fits = mrdfits('../ExpTimeCalc/fov24_1p0_bkgnd_param.fits')
  frac_3d = mean(frac_fits, dimension=1)
  frac_2d = mean(frac_3d, dimension=1)
  frac = frac_2d[2,*]
  elon=110.
  elat=30.
  noises = dblarr(n_elements(imag), npix_max-npix_min+1)
  dilution = dblarr(n_elements(imag), npix_max-npix_min+1)
  diln = dblarr(n_elements(imag))
  satn = dblarr(n_elements(imag))
  fov_ind = intarr(n_elements(imag)) + 2
  exptime = dblarr(n_elements(imag)) + 3600. ; what's the noise per hour?
  for ii=0, (npix_max-npix_min) do begin
     calc_noise, imag, exptime, noise, $
                 npix_aper=(ii+npix_min), $
                 frac_aper=frac[ii+npix_min-1], $
		 field_angle = field_angle, $
		 fov_ind = fov_ind, $
		 bk_p = bk_fits, $
                 teff=teff, $
                 e_pix_ro = 10.,$
                 subexptime=2., $
                 geom_area = geom_area, $
                 sys_lim = 60.0, $
                 pix_scale = pixsc, $
                 elon=elon, $
                 elat=elat, $
                 dilution=dil, $
                 e_star_sub=estar
  dilution[*,ii] = dil
    noises[*,ii] = (1.0+dil)*noise
    if (ii eq 0) then begin
	satn = (estar gt 150000.)
        print, noise
    end
  end
  minnoise = min(noises, ind, dimension=2)
  npix = ind / n_elements(imag) + npix_min
  print, median(npix)
  calc_noise, imag, exptime, noise, $
                 npix_aper=npix, $
                 frac_aper=frac[npix-1], $
		 field_angle = field_angle, $
		 fov_ind = fov_ind, $
		 bk_p = bk_fits, $
                 teff=teff, $
                 e_pix_ro = 10.,$
                 subexptime=2., $
                 geom_area = geom_area, $
                 sys_lim = 60.0, $
                 pix_scale = pixsc, $
                 elon=elon, $
                 elat=elat, $
                 dilution=diln, $
                 e_star_sub=estar, $
                 noise_star=shot_noise, $
	 	 noise_sky=bknd_noise, $
 		 noise_ro=read_noise, $
		 noise_sys=sys_noise
 
  fitsout = [[imag], [npix], [noise], [shot_noise], [bknd_noise], [read_noise], [sys_noise], [satn], [diln]]
  mwrfits, fitsout, 'noises.fits'
  
END


