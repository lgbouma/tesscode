pro eclip_observe, eclipse, star, bk, deep, frac, rad, ph_p, $
	aspix=aspix, effarea=effarea, readnoise=readnoise, $
        tranmin=tranmin, thresh=thresh, $
	ps_len=ps_len, ffi_len=ffi_len, saturation=saturation, $
	sys_limit=sys_limit, duty_cycle=duty_cycle, $
        dwell_time=dwell_time, downlink=downlink, $
        subexptime=subexptime
 REARTH_IN_RSUN = 0.0091705248
;;;;;; basic parameters here
  if (keyword_set(ps_len)) then ps_len=ps_len else ps_len = 2.0
  if (keyword_set(thresh)) then SNR_MIN=thresh else SNR_MIN = 7.0
  if (keyword_set(tranmin)) then NTRA_OBS_MIN = tranmin else NTRA_OBS_MIN = 2
  if (keyword_set(effarea)) then effarea=effarea else effarea=61.2; cm^2
  if (keyword_set(sys_limit)) then sys_limit=sys_limit else sys_limit=60.0
  if (keyword_set(duty_cycle)) then duty_cycle=duty_cycle else duty_cycle=100.0
  if (keyword_set(ffi_len)) then ffi_len=ffi_len else ffi_len=30.
  if (keyword_set(subexptime)) then subexptime=subexptime else subexptime=2.
  if (keyword_set(saturation)) then saturation=saturation else saturation=150000.
  if (keyword_set(dwell_time)) then dwell_time=dwell_time else dwell_time=13.66d0 ; days per orbit
  if (keyword_set(downlink)) then downlink=downlink else downlink=16.0d/24.0d
  apo_blank = (dwell_time-downlink)*(1.0-duty_cycle/100.0)
  
  npix_max = 64
  npix_min = 3
  mask2d = intarr(16,16)
  mid = indgen(8)+4
  mask1 = mask2d
  mask2 = mask2d
  mask1[*,mid] = 1
  mask2[mid,*] = 1
  mask2d = mask1*mask2
  mask1d = reform(mask2d, 16*16)

  nstars = n_elements(star)
  neclipses = n_elements(eclipse)
  print, 'Observing ', neclipses, ' eclipses around ', nstars,' stars.'
  dx = floor(10*randomu(seed, neclipses))
  dy = floor(10*randomu(seed, neclipses))
  
  ecid = eclipse.hostid
  print, 'Calculating number of eclipses'
  eclipse.neclip_obs1 = $
      n_eclip(eclipse.p, dwell_time, $
      2.0*double(eclipse.npointings), periblank=downlink, apoblank=apo_blank, ein=e1)
 ; Flip the phase
  e2 = 1.0-e1
  eclipse.neclip_obs2 = $
      n_eclip(eclipse.p, DWELL_TIME, $
      2.0*double(eclipse.npointings), periblank=downlink, apoblank=apo_blank, ein=e2)

  print, 'Diluting FFIs'
  tra_ps = where(star[eclipse.hostid].ffi lt 1)
  if (tra_ps[0] ne -1) then begin
    dur1_min = eclipse[tra_ps].dur1*24.0*60.0
    dur2_min = eclipse[tra_ps].dur2*24.0*60.0
    eclipse[tra_ps].dep1_eff = eclipse[tra_ps].dep1*dil_ffi_eclip(dur1_min, float(ps_len), ffis=nps)
    eclipse[tra_ps].dur1_eff = (nps*ps_len)/(24.0*60.0)
    eclipse[tra_ps].dep2_eff = eclipse[tra_ps].dep2*dil_ffi_eclip(dur2_min, float(ps_len), ffis=nps)
    eclipse[tra_ps].dur2_eff = (nps*ps_len)/(24.0*60.0)
  endif

  tra_ffi = where(star[eclipse.hostid].ffi gt 0)
  if (tra_ffi[0] ne -1) then begin
    dur1_min = eclipse[tra_ffi].dur1*24.0*60.0
    dur2_min = eclipse[tra_ffi].dur2*24.0*60.0
    eclipse[tra_ffi].dep1_eff = eclipse[tra_ffi].dep1*dil_ffi_eclip(dur1_min, float(ffi_len), ffis=ffis)
    eclipse[tra_ffi].dur1_eff = (ffis*ffi_len)/(24.0*60.0)
    eclipse[tra_ffi].dep2_eff = eclipse[tra_ffi].dep2*dil_ffi_eclip(dur2_min, float(ffi_len), ffis=ffis)
    eclipse[tra_ffi].dur2_eff = (ffis*ffi_len)/(24.0*60.0)
  endif
; for each observed transiting eclipse, calculate snr
 
  obs = where((eclipse.neclip_obs1 + eclipse.neclip_obs2) gt NTRA_OBS_MIN)
  if (obs[0] ne -1) then begin
    obsid = eclipse[obs].hostid
    nobs = n_elements(obs)
    print, 'Observing ', nobs, ' transits'
    ;  obs = indgen(n_elements(star))

    print, 'Stacking and sorting PRFs'
    ; ph_star is npix x nstar
    stack_prf_eclip, star[obsid].mag.t, star[obsid].teff, ph_p, frac, star_ph, $
	dx=dx[obs], dy=dy[obs], fov_ind=eclipse[obs].coord.fov_ind, mask=mask1d
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
                 field_angle=eclipse[obs].coord.field_angle, $
		 subexptime=subexptime, $
                 geom_area = effarea, $
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
    for ss=0,nobs-1 do eclipse[obs[ss]].star_ph = star_ph[eclipse[obs[ss]].npix-1,ss]
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
  
; Dilute
    det = where(((eclipse.neclip_obs1 + eclipse.neclip_obs2) ge NTRA_OBS_MIN) and $
	      (eclipse.snr ge SNR_MIN))
    if (det[0] ne -1) then begin
      ;print, 'Calculating noise again'
      detid = eclipse[det].hostid
      ndet = n_elements(det)
      print, 'Re-observing ', ndet, ' transits'
      
      print, 'Stacking and sorting PRFs'
      ; ph_star is npix x nstar
      stack_prf_eclip, star[detid].mag.t, star[detid].teff, ph_p, frac, star_ph, $
	dx=dx[det], dy=dy[det], fov_ind=eclipse[det].coord.fov_ind, mask=mask1d, sind=sind
      ;bk_ph = eclipse[det].bk_ph

      dil_ph = dblarr(total(mask1d), ndet)
      
;     print, "Diluting with binary companions"
      bindil = where(eclipse[det].class eq 1)
      if (bindil[0] ne -1) then begin
      stack_prf_eclip, star[star[bindil].companion.ind].mag.t, star[detid].teff, ph_p, frac, dilvec, $
	dx=dx[det], dy=dy[det], fov_ind=eclipse[det].coord.fov_ind, mask=mask1d, sind=sind
;	dilute_binary, eclipse[det[bindil]], star, frac, rad, ph_p, $
;		dx[det[bindil]], dy[det[bindil]], dilvec, aspix=aspix, radmax=4.0
;        dil_ph[bindil] = dil_ph[bindil] + dilvec
      end
      print, "Diluting with other target stars"
      targdil = where(eclipse[det].class eq 1 or eclipse[det].class eq 2)
      if (targdil[0] ne -1) then begin
        dilute_eclipse_img, eclipse[det[targdil]], star, frac, ph_p, $
		dx[det[targdil]], dy[det[targdil]], dilvec, aspix=aspix, sq_deg=13.4, radmax=6.0
        dil_ph[targdil] = dil_ph[targdil] + dilvec
      end
      print, "Diluting with background stars"
      dilute_eclipse_img, eclipse[det], bk, frac, ph_p, $
		dx[det], dy[det], dilvec, aspix=aspix, sq_deg=0.134, radmax=4.0
      dil_ph = dil_ph + dilvec
      print, "Diluting with deep stars"
      dilute_eclipse_img, eclipse[det], deep, frac, ph_p, $
		dx[det], dy[det], dilvec, aspix=aspix, sq_deg=0.0134, radmax=2.0
      dil_ph = dil_ph + dilvec
   
      ; Sort into the same pixel order as target star flux
      for jj=0,ndet-1 do begin
        dil_ph[*,jj] = total(dil_ph[sind[*,jj],jj], /cumulative)
      end      

      zodi_flux, eclipse[det].coord.elat, aspix, zodi_ph
      eclipse[det].zodi_ph = zodi_ph

      noises = dblarr(n_elements(det), npix_max)
      dilution = dblarr(n_elements(det), npix_max)
      ;shot_noises = dblarr(n_elements(obs), npix_max-npix_min+1)
      exptime = dblarr(n_elements(det)) + 3600.
      for ii=0,(npix_max-1) do begin
        calc_noise_eclip, star_ph, dil_ph, exptime, $
		 readnoise, sys_limit, noise, $
		 npix_aper=(ii+1), $
                 field_angle=eclipse[det].coord.field_angle, $
		 subexptime=subexptime, $
                 geom_area = effarea, $
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
        if (ii eq 0) then eclipse[det].sat = (estar gt SATURATION)
      end
      noises = noises[*,(npix_min-1):(npix_max-1)]
      dilution = dilution[*,(npix_min-1):(npix_max-1)]
      minnoise = min(noises, ind, dimension=2)
      eclipse[det].npix = ind / n_elements(det) + npix_min
      eclipse[det].dil = dilution[ind]
      eclipse[det].snrhr = minnoise
      for ss=0,ndet-1 do eclipse[det[ss]].star_ph = star_ph[eclipse[det[ss]].npix-1,ss]
      ;eclipse[obs].star_ph = reform(star_ph[eclipse[obs].npix-1, *])
      ; Calculate SNR in phase-folded lightcurve
      et1_folded = double(eclipse[det].neclip_obs1) * $
      eclipse[det].dur1_eff * 24.0 * 3600
      et2_folded = double(eclipse[det].neclip_obs2) * $
      eclipse[det].dur2_eff * 24.0 * 3600
      et1_eclp = eclipse[det].dur1_eff * 24.0 * 3600
      et2_eclp = eclipse[det].dur2_eff * 24.0 * 3600
      eclipse[det].snr1 = eclipse[det].dep1_eff / (minnoise) * sqrt(et1_folded/3600.)
      eclipse[det].snr2 = eclipse[det].dep2_eff / (minnoise) * sqrt(et2_folded/3600.)
      eclipse[det].snr  = sqrt(eclipse[det].snr1^2. + eclipse[det].snr2^2.)
      eclipse[det].snreclp1 = eclipse[det].dep1_eff / (minnoise) * sqrt(et1_eclp/3600.)
      eclipse[det].snreclp2 = eclipse[det].dep2_eff / (minnoise) * sqrt(et2_eclp/3600.)

      det1 = where((eclipse.neclip_obs1 ge NTRA_OBS_MIN) and $
	      (eclipse.snr1 ge SNR_MIN))
      det2 = where((eclipse.neclip_obs2 ge NTRA_OBS_MIN) and $
	      (eclipse.snr2 ge SNR_MIN))
      det = where(((eclipse.neclip_obs1 + eclipse.neclip_obs2) ge NTRA_OBS_MIN) and $
	      (eclipse.snr ge SNR_MIN))
      if (det1[0] ne -1) then eclipse[det1].det1 = 1
      if (det2[0] ne -1) then eclipse[det2].det2 = 1
      if (det[0] ne -1)  then begin
        eclipse[det].det = 1
        print, 'Detected ', n_elements(det), ' eclipses.'
      endif
      ;  detected = where(star.eclipse_hz.tra gt 0 and star.eclipse_hz.neclp1_obs ge NTRA_OBS_MIN and star.eclipse_hz.snr ge SNR_MIN)
      ;  star[detected].eclipse_hz.det = 1
    end ;det if
  end ; obs if
end ; PRO
