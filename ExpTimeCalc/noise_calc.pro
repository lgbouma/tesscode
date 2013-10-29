  frac_fits = mrdfits('frac24_1p0.fits')
  noises = dblarr(n_elements(obs), npix_max)
  dilution = dblarr(n_elements(obs), npix_max)
  shot_noises = dblarr(n_elements(obs), npix_max)
  exptime = dblarr(n_elements(obs)) + 3600. ; what's the noise per hour?
  for ii=(npix_min-1),(npix_max-1) do begin
     thisfrac = frac_fits[*,*,*,ii]
     frac = thisfrac[star[obs].dx, star[obs].dy,fov_ind]
     calc_noise, star[obs].mag.i, exptime, noise, $
                 npix_aper=(ii+1), $
                 frac_aper=frac, $
                 teff=star[obs].teff, $
                 e_pix_ro = E_PIX_RO,$
                 subexptime=SUB_EXP_TIME, $
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
                 elon=star[obs].coord.elon, $
                 elat=star[obs].coord.elat, $
                 dilution=dil, $
                 e_star_sub=estar, $
                 noise_star=shot_noise
    dilution[*,ii] = dil
    noises[*,ii] = (1.0+dil)*noise
    shot_noises[*,ii] = shot_noise*1d6 ; in ppm

