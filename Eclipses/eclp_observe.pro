pro eclp_observe, sstruct=sstruct, pstruct=pstruct, sfile=sfile, pfile=pfile, outfile=outfile, $
	fov=fov, geomarea=geomarea, readnoise=readnoise, tranmin=tranmin, thresh=thresh, $
	nodil=nodil, ffi_len=ffi_len, $
	sys_limit=sys_limit, keep_ntra1=keep_ntra1, duty_cycle=duty_cycle, $
        bk_file=bk_file, sp_file=sp_file, ph_fits=ph_fits, prf_fits=prf_fits
	

 REARTH_IN_RSUN = 0.0091705248
;;;;;; basic parameters here

  DWELL_TIME = 13.66d0 ;27.4d0             ; days per field
  DOWNLINK_TIME = 16.0d0/24.0d0 ; correct for downlink time
  if (keyword_set(thresh)) then SNR_MIN=thresh else SNR_MIN = 7.0
  if (keyword_set(tra1nmin)) then NTRA_OBS_MIN = tra1nmin else NTRA_OBS_MIN = 2
  if (keyword_set(geomarea)) then GEOM_AREA=geomarea else GEOM_AREA=61.2; cm^2
  if (keyword_set(fov)) then fov=fov else fov=24.0
  if (keyword_set(sys_limit)) then sys_limit=sys_limit else sys_limit=60.0
  if (keyword_set(duty_cycle)) then duty_cycle=duty_cycle else duty_cycle=100.0
  if (keyword_set(ffi_len)) then ffi_len=ffi_len else ffi_len=30.
  apo_blank = (DWELL_TIME-DOWNLINK_TIME)*(1.0-duty_cycle/100.0)
 ;SYS_LIMIT = 60.0; ppm in 1 hour
  E_PIX_RO = 10.0 
  SUB_EXP_TIME = 2.0
  SATURATION = 150000.
  npix_max = 49
  npix_min = 3
  ; photon counts
  if (keyword_set(ph_file)) then ph_fits=mrdfits(ph_file) else ph_fits=0
  ; big frac file
  if (keyword_set(prf_file)) then prf_file=prf_file else prf_file='../ExpTimeCalc/bigfrac24_105_4700k.fits'
  prf_fits = mrdfits(prf_file)
  ;prf_fits = mean(mean(mrdfits(prf_file),dimension=1), dimension=1)
  ; background file
  if (keyword_set(bk_file)) then bk_fits = mrdfits(bk_file) else bk_fits=0
  ; spline fits to PSRR
  if (keyword_set(sp_file)) then begin
	sp_fits = mrdfits(sp_file) 
 	sp_fits[4,*] = sp_fits[4,*]-sp_fits[4,0]
 	sp_fits[5,*] = sp_fits[5,*]-sp_fits[5,0]
 	sp_fits[6,*] = sp_fits[6,*]-sp_fits[6,0]
 	sp_fits[7,*] = sp_fits[7,*]-sp_fits[7,0]
  endif
  ccd_pix = 4096.0
  PIX_SCALE = fov*3600./ccd_pix
  pix_scale_deg = fov/ccd_pix


  if keyword_set(sstruct) then star=sstruct  else restore, sfile
  if keyword_set(pstruct) then eclipse=pstruct else restore, pfile
  ;restore, infile
  nstars = n_elements(star)
  neclipses = n_elements(eclipse)
  print, 'Observing ', neclipses, ' eclipses around ', nstars,' stars.'
  
  ecid = eclipse.hostid

   print, 'Calculating number of eclipses'
  eclipse.neclp_obs1 = $
      n_eclip_blank(eclipse.p, DWELL_TIME, $
      2.0*double(star[ecid].npointings), periblank=DOWNLINK_TIME, apoblank=apo_blank, ein=e1)
  ; Flip the phase
  e2 = 1.0-e1
  eclipse.neclp_obs2 = $
      n_eclip_blank(eclipse.p, DWELL_TIME, $
      2.0*double(star[ecid].npointings), periblank=DOWNLINK_TIME, apoblank=apo_blank, ein=e1)

  print, 'Diluting FFIs'
  tra_ps = where(star[eclipse.hostid].ffi lt 1)
  if (tra1_ps[0] ne -1) then begin
    dur1_min = eclipse[tra_ps].dur1*24.0*60.0
    dur2_min = eclipse[tra_ps].dur2*24.0*60.0
    eclipse[tra_ps].dep1_eff = eclipse[tra_ps].dep1*dil_ffi(dur1_min, 1.0, ffis=nps)
    eclipse[tra_ps].dur1_eff = (nps*1.0)/(24.0*60.0)
    eclipse[tra_ps].dep2_eff = eclipse[tra_ps].dep2*dil_ffi(dur2_min, 1.0, ffis=nps)
    eclipse[tra_ps].dur2_eff = (nps*1.0)/(24.0*60.0)
  endif

  tra_ffi = where(star[eclipse.hostid].ffi gt 0)
  if (tra1_ffi[0] ne -1) then begin
    dur1_min = eclipse[tra_ffi].dur1*24.0*60.0
    dur2_min = eclipse[tra_ffi].dur2*24.0*60.0
    eclipse[tra_ffi].dep1_eff = eclipse[tra_ffi].dep1*dil_ffi(dur1_min, float(ffi_len), ffis=ffis)
    eclipse[tra_ffi].dur1_eff = (ffis*ffi_len)/(24.0*60.0)
    eclipse[tra_ffi].dep2_eff = eclipse[tra_ffi].dep2*dil_ffi(dur2_min, float(ffi_len), ffis=ffis)
    eclipse[tra_ffi].dur2_eff = (ffis*ffi_len)/(24.0*60.0)
  endif
; for each observed tra1nsiting eclipse, calculate snr
 
   obs = where((eclipse.neclp_obs1 gt 0) or (eclipse.neclp_obs2 gt 0))
   obsid = eclipse[obs].hostid
   nobs = n_elements(obs)
   print, 'Observing ', nobs, ' transits'
;  obs = indgen(n_elements(star))

  fov_ind = intarr(n_elements(obs))
  fov_ind[where((star[obsid].coord.fov_r ge 0.104*CCD_PIX) and $
		(star[obsid].coord.fov_r lt 0.365*CCD_PIX))] = 1
  fov_ind[where((star[obsid].coord.fov_r ge 0.365*CCD_PIX) and $
		(star[obsid].coord.fov_r lt 0.592*CCD_PIX))] = 2
  fov_ind[where(star[obsid].coord.fov_r  ge 0.592*CCD_PIX)] = 3
  field_angle = star[obsid].coord.fov_r / CCD_PIX * FOV 
 
  print, 'Stacking and sorting PRFs'
  ; ph_star is npix x nstar
  dx = floor(10*randomu(seed, nobs))
  dy = floor(10*randomu(seed, nobs))
  stack_prf, star[obsid].mag.t, star[obsid].teff, ph_fits, prf_fits, ph_star, dx=dx, dy=dy, fov_ind=fov_ind

  bin_sys = (star[obsid].pri or star[obsid].sec)
  bin_sep = star[obsid].companion.sep
  bin_tmag = star[star[obsid].companion.ind].mag.t
  bin_teff = star[star[obsid].companion.ind].teff
  bins = where(bin_sys)

  stack_prf, bin_tmag, bin_teff, ph_fits, prf_fits, ph_bin, dx=dx, dy=dy, fov_ind=fov_ind
 
  pix_sep = bin_sep[bins]/pix_scale   ; from arcsec to pixels
  r = sp_fits[fov_ind[bins],*]       ; distance (in pixels)
  di = sp_fits[fov_ind[bins]+4,*]    ; imag attenuation
  dibin = interpol(di, r, pix_sep)    ; interpolate over spline fit
  ph_bin = 10.0^(-0.4*dibin) * ph_bin ; attenuate the count rate
  
  print, 'Calculating noise'
  noises = dblarr(n_elements(obs), npix_max-npix_min+1)
  dilution = dblarr(n_elements(obs), npix_max-npix_min+1)
  shot_noises = dblarr(n_elements(obs), npix_max-npix_min+1)
  exptime = dblarr(n_elements(obs)) + 3600.
  for ii=0,(npix_max-npix_min) do begin
     calc_noise_filt, ph_star, exptime, noise, $
		 npix_aper=(ii+npix_min), $
                 field_angle=field_angle, $
  		 fov_ind=fov_ind, $
                 e_pix_ro = E_PIX_RO,$
		 subexptime=SUB_EXP_TIME, $
                 geom_area = GEOM_AREA, $
                 sys_lim = SYS_LIMIT, $
                 pix_scale = PIX_SCALE, $
		 bk_p = bk_fits, $
                 elon=star[obsid].coord.elon, $
                 elat=star[obsid].coord.elat, $
		 dilution=dil, $
		 e_star_sub=estar, $
		 noise_star=shot_noise, $
		 bin_sys= bin_sys, $
 		 bin_ph = ph_bin, $
  		 bin_sep = bin_sep
    dilution[*,ii] = dil
    if (keyword_set(nodil)) then noises[*,ii] = noise $
	else noises[*,ii] = dil*noise
    shot_noises[*,ii] = shot_noise*1d6
    ;print, median(shot_noise*1d6)
    ;if (ii eq 0) then star[obsid].sat = (estar gt SATURATION)
  end
  minnoise = min(noises, ind, dimension=2)
  star[obsid].npix = ind / n_elements(obs) + npix_min
  star[obsid].snr = 1.0/minnoise
  star[obsid].dil = dilution[ind]
  ; Calculate SNR in phase-folded lightcurve
  et1_folded = double(eclipse[obs].neclp_obs1) * $
    eclipse[obs].dur1_eff * 24.0 * 3600
  et2_folded = double(eclipse[obs].neclp_obs2) * $
    eclipse[obs].dur2_eff * 24.0 * 3600
  et1_eclp = eclipse[obs].dur1_eff * 24.0 * 3600
  et2_eclp = eclipse[obs].dur2_eff * 24.0 * 3600
  eclipse[obs].snr1 = eclipse[obs].dep1_eff / (minnoise) * sqrt(et1_folded/3600.)
  eclipse[obs].snr2 = eclipse[obs].dep2_eff / (minnoise) * sqrt(et2_folded/3600.)
  eclipse[obs].snr  = sqrt(eclipse[obs].snr1^2. + eclipse[obs].snr2^2.)
  eclipse[obs].snreclp1 = eclipse[obs].dep1_eff / (minnoise) * sqrt(et1_eclp/3600.)
  eclipse[obs].snreclp2 = eclipse[obs].dep2_eff / (minnoise) * sqrt(et2_eclp/3600.)
;  eclipse[obs].snrgress1 = eclipse[obs].snreclp1 * $
; 		sqrt(2. * eclipse[obs].neclp_obs1 * $
;		     REARTH_IN_RSUN * eclipse[obs].r / $
;		    (6.0*star[obsid].r*(1.0+eclipse[obs].b^2.)))
;   decide if it is 'detected'.
  det1 = where((eclipse.neclp_obs1 ge NTRA_OBS_MIN) and $
	      (eclipse.snr1 ge SNR_MIN))
  det2 = where((eclipse.neclp_obs2 ge NTRA_OBS_MIN) and $
	      (eclipse.snr2 ge SNR_MIN))
  det = where(((eclipse.neclp_obs1 + eclipse.neclp_obs2) ge NTRA_OBS_MIN) and $
	      (eclipse.snr ge SNR_MIN))
  print, 'Detected ', n_elements(det), ' eclipses.'
  eclipse[det1].det1 = 1
  eclipse[det2].det2 = 1
  eclipse[det].det = 1
;  detected = where(star.eclipse_hz.tra1 gt 0 and star.eclipse_hz.neclp1_obs ge NTRA_OBS_MIN and star.eclipse_hz.snr ge SNR_MIN)
;  star[detected].eclipse_hz.det = 1
  if keyword_set(sstruct) then sstruct=star 
  if keyword_set(pstruct) then pstruct=eclipse 
  if keyword_set(outfile) then save, filen=outfile, eclipse

end
