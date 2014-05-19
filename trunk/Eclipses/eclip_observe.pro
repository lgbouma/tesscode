pro eclip_observe, eclipse, star, frac, rad, ph_p, $
	aspix=aspix, geomarea=geomarea, readnoise=readnoise, $
        tranmin=tranmin, thresh=thresh, $
	nodil=nodil, ffi_len=ffi_len, saturation=saturation, $
	sys_limit=sys_limit, duty_cycle=duty_cycle
 REARTH_IN_RSUN = 0.0091705248
;;;;;; basic parameters here

  DWELL_TIME = 13.66d0 ;27.4d0             ; days per field
  DOWNLINK_TIME = 16.0d0/24.0d0 ; correct for downlink time
  if (keyword_set(thresh)) then SNR_MIN=thresh else SNR_MIN = 7.0
  if (keyword_set(tranmin)) then NTRA_OBS_MIN = tranmin else NTRA_OBS_MIN = 2
  if (keyword_set(geomarea)) then GEOM_AREA=geomarea else GEOM_AREA=61.2; cm^2
  if (keyword_set(sys_limit)) then sys_limit=sys_limit else sys_limit=60.0
  if (keyword_set(duty_cycle)) then duty_cycle=duty_cycle else duty_cycle=100.0
  if (keyword_set(ffi_len)) then ffi_len=ffi_len else ffi_len=30.
  if (keyword_set(saturation)) then saturation=saturation else saturation=150000.
  apo_blank = (DWELL_TIME-DOWNLINK_TIME)*(1.0-duty_cycle/100.0)
 ;SYS_LIMIT = 60.0; ppm in 1 hour
  SUB_EXP_TIME = 2.0
  npix_max = 49
  npix_min = 3

  nstars = n_elements(star)
  neclipses = n_elements(eclipse)
  print, 'Observing ', neclipses, ' eclipses around ', nstars,' stars.'
  
  ecid = eclipse.hostid
  print, 'Calculating number of eclipses'
  eclipse.neclip_obs1 = $
      n_eclip_blank(eclipse.p, DWELL_TIME, $
      2.0*double(eclipse.npointings), periblank=DOWNLINK_TIME, apoblank=apo_blank, ein=e1)
  ; Flip the phase
  e2 = 1.0-e1
  eclipse.neclip_obs2 = $
      n_eclip_blank(eclipse.p, DWELL_TIME, $
      2.0*double(eclipse.npointings), periblank=DOWNLINK_TIME, apoblank=apo_blank, ein=e2)

  print, 'Diluting FFIs'
  tra_ps = where(star[eclipse.hostid].ffi lt 1)
  if (tra_ps[0] ne -1) then begin
    dur1_min = eclipse[tra_ps].dur1*24.0*60.0
    dur2_min = eclipse[tra_ps].dur2*24.0*60.0
    eclipse[tra_ps].dep1_eff = eclipse[tra_ps].dep1*dil_ffi(dur1_min, 1.0, ffis=nps)
    eclipse[tra_ps].dur1_eff = (nps*1.0)/(24.0*60.0)
    eclipse[tra_ps].dep2_eff = eclipse[tra_ps].dep2*dil_ffi(dur2_min, 1.0, ffis=nps)
    eclipse[tra_ps].dur2_eff = (nps*1.0)/(24.0*60.0)
  endif

  tra_ffi = where(star[eclipse.hostid].ffi gt 0)
  if (tra_ffi[0] ne -1) then begin
    dur1_min = eclipse[tra_ffi].dur1*24.0*60.0
    dur2_min = eclipse[tra_ffi].dur2*24.0*60.0
    eclipse[tra_ffi].dep1_eff = eclipse[tra_ffi].dep1*dil_ffi(dur1_min, float(ffi_len), ffis=ffis)
    eclipse[tra_ffi].dur1_eff = (ffis*ffi_len)/(24.0*60.0)
    eclipse[tra_ffi].dep2_eff = eclipse[tra_ffi].dep2*dil_ffi(dur2_min, float(ffi_len), ffis=ffis)
    eclipse[tra_ffi].dur2_eff = (ffis*ffi_len)/(24.0*60.0)
  endif
; for each observed transiting eclipse, calculate snr
 
   obs = where((eclipse.neclip_obs1 gt 0) or (eclipse.neclip_obs2 gt 0))
   if (obs[0] ne -1) then begin
   obsid = eclipse[obs].hostid
   nobs = n_elements(obs)
   print, 'Observing ', nobs, ' transits'
;  obs = indgen(n_elements(star))

  print, 'Stacking and sorting PRFs'
  ; ph_star is npix x nstar
  dx = floor(10*randomu(seed, nobs))
  dy = floor(10*randomu(seed, nobs))
  stack_prf_eclip, star[obsid].mag.t, star[obsid].teff, ph_p, frac, star_ph, $
	dx=dx, dy=dy, fov_ind=eclipse[obs].coord.fov_ind
  bk_ph = eclipse[obs].bk_ph

  zodi_flux, eclipse[obs].coord.elat, aspix, zodi_ph
  eclipse[obs].zodi_ph = zodi_ph

  print, 'Calculating noise'
  noises = dblarr(n_elements(obs), npix_max)
  dilution = dblarr(n_elements(obs), npix_max)
  ;shot_noises = dblarr(n_elements(obs), npix_max-npix_min+1)
  exptime = dblarr(n_elements(obs)) + 3600.
  for ii=0,(npix_max-1) do begin
     calc_noise_eclip, star_ph, bk_ph, exptime, $
		 readnoise, sys_limit, noise, $
		 npix_aper=(ii+1), $
                 field_angle=eclipse[obs].coord.fov_r, $
		 subexptime=SUB_EXP_TIME, $
                 geom_area = GEOM_AREA, $
                 aspix=aspix, $
                 zodi_ph=zodi_ph, $
		 dilution=dil, $
		 e_tot_sub=estar, $
		 noise_star=shot_noise
    dilution[*,ii] = dil
    if (keyword_set(nodil)) then noises[*,ii] = noise $
	else noises[*,ii] = dil*noise
    ;shot_noises[*,ii] = shot_noise*1d6
    ;print, median(shot_noise*1d6)
    if (ii eq 0) then eclipse[obs].sat = (estar gt SATURATION)
  end
  noises = noises[*,(npix_min-1):(npix_max-1)]
  dilution = dilution[*,(npix_min-1):(npix_max-1)]
  minnoise = min(noises, ind, dimension=2)
  eclipse[obs].npix = ind / n_elements(obs) + npix_min
  eclipse[obs].dil = dilution[ind]
  ;eclipse[obs].star_ph = reform(star_ph[eclipse[obs].npix-1, *])
  ; Calculate SNR in phase-folded lightcurve
  et1_folded = double(eclipse[obs].neclip_obs1) * $
    eclipse[obs].dur1_eff * 24.0 * 3600
  et2_folded = double(eclipse[obs].neclip_obs2) * $
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
  det1 = where((eclipse.neclip_obs1 ge NTRA_OBS_MIN) and $
	      (eclipse.snr1 ge SNR_MIN))
  det2 = where((eclipse.neclip_obs2 ge NTRA_OBS_MIN) and $
	      (eclipse.snr2 ge SNR_MIN))
  det = where(((eclipse.neclip_obs1 + eclipse.neclip_obs2) ge NTRA_OBS_MIN) and $
	      (eclipse.snr ge SNR_MIN))
  print, 'Detected ', n_elements(det), ' eclipses.'
  if (det1[0] ne -1) then eclipse[det1].det1 = 1
  if (det2[0] ne -1) then eclipse[det2].det2 = 1
  if (det[0] ne -1)  then eclipse[det].det = 1
;  detected = where(star.eclipse_hz.tra gt 0 and star.eclipse_hz.neclp1_obs ge NTRA_OBS_MIN and star.eclipse_hz.snr ge SNR_MIN)
;  star[detected].eclipse_hz.det = 1
  end
end
