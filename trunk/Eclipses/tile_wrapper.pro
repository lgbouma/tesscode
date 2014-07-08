PRO tile_wrapper, fpath, fnums, outname, eclip=eclip, n_trial=n_trial
  numfil = n_elements(fnums)

  ; User-adjustable settings (yes, that's you!)
  frac_file = 'bigfrac24_105_f3p33.fits' ; prf file 
  rad_file = 'bigrad24_105_f3p33.fits' ; radius file 
  ph_file = 'ph_T_filt.fits' ; photon fluxes for T=10 vs Teff
  cr_file = 'crnoise.fits' ; photon fluxes for T=10 vs Teff
  tic_file = 'tic_teff.fits'
  dart_file = 'dartmouth.sav'
  fov = 24.
  seg = 13
  skirt=6.
  effarea = 69.1 ; in cm^2. 
  aspix = 21.1 ; arcseconds per pixel
  readnoise= 10.0 ; in e- per subexposure
  subexptime = 2.0 ; sec in subexposure
  thresh = 7.0 ; detection threshold in phase-folded lightcurve
  tranmin = 2.0 ; minimum number of eclipses for detection
  min_depth=1D-6 ; minimum transit depth to retain from eclipses
  sys_limit=60. ; in ppm/hr
  if (keyword_set(n_trial)) then n_trial=n_trial else n_trial = 10 ; number of trials in this run
  ffi_len=30. ; in minutes
  ps_len=30. ; in minutes
  duty_cycle=100.+fltarr(numfil)
  nps = 100000
  saturation=150000.
  CCD_PIX = 4096.
  orbit_period = 13.66d0 ; days per orbit
  downlink = 16.0d0/24.0d0 ; downlink time in days
  eclass = [	0, $ ; Planets
	    	1, $ ; EBs
		1, $ ; BEBs
		1  ] ; HEBs

  ; Don't phuck with physics, though
  REARTH_IN_RSUN = 0.0091705248
  AU_IN_RSUN = 215.093990942
  nparam = 49 ; output table width

  ; Here we go!
  numtargets = lonarr(numfil)
  numbkgnd = lonarr(numfil)
  numdeeps = lonarr(numfil)
  numps = lonarr(numfil)
  ; Open the fits files
  frac_fits = mrdfits(frac_file)
  rad_fits = mrdfits(rad_file)/60. ; put into pixels
  ph_fits = mrdfits(ph_file)
  cr_fits = fltarr(100,64)
;  cr_fits = mrdfits(cr_file)
  restore, dart_file
  dartstruct = ss
  tic_fits = mrdfits(tic_file)
  ; Make random spherical coords
  u = randomu(seed, 1D7)
  v = randomu(seed, 1D7)
  phi = 2.*!dpi*u
  theta = acos(2.*v-1.)
  ang2pix_ring, 16, theta, phi, ipring
   
  totdet = 0L
  star_out = dblarr(1E7*n_trial,nparam)
  for ii=0, numfil-1 do begin
    ; Gather the .sav files
    print, 'Restoring files for tile ', fnums[ii]
    fname = fpath+'hp'+string(fnums[ii], format='(I04)')+'.sav'
    restore, fname
    targets = star
    numtargets[ii] = n_elements(targets)
    fname = fpath+'bk'+string(fnums[ii], format='(I04)')+'.sav'
    restore, fname
    bkgnds = star[where(star.mag.k gt 15)]
    numbkgnd[ii] = n_elements(bkgnds)
    fname = fpath+'dp'+string(fnums[ii], format='(I04)')+'.sav'
    restore, fname
    deeps = star[where(star.mag.t gt 21)]
    numdeeps[ii] = n_elements(deeps)
    delvar, star
   
    ; Choose which stars are postage stamps vs ffis
    targets.ffi = 1
    pri = where(targets.pri eq 1)
    selpri = ps_sel(targets[pri].mag.t, targets[pri].teff, targets[pri].m, targets[pri].r, ph_fits)
    if (selpri[0] ne -1) then begin 
      targets[pri[selpri]].ffi=0
      secffi = targets[pri[selpri]].companion.ind
      targets[secffi].ffi=0
      numps[ii] = numps[ii]+n_elements(selpri)
    end
    sing = where((targets.pri eq 0) and (targets.sec eq 0))
    selsing = ps_sel(targets[sing].mag.t, targets[sing].teff, targets[sing].m, targets[sing].r, ph_fits)
    if (selsing[0] ne -1) then begin 
      targets[sing[selsing]].ffi=0
      numps[ii] = numps[ii]+n_elements(selsing)
    end
    ecliplen_tot = 0L
    for jj=0,n_trial-1 do begin
      ; re-radomize the inclination
      targets.cosi = -1 + 2.0*randomu(seed, n_elements(targets))
      ; Add eclipses
      ecliplen =  make_eclipse(targets, bkgnds, eclip_trial, frac_fits, $
	rad_fits, ph_fits, dartstruct, tic_fits, eclass, min_depth=min_depth)
      if (ecliplen gt 0) then begin
        eclip_trial.trial = jj + 1
        ; Add coordinates to the eclipses
        thispix = where(ipring eq fnums[ii])
        ncoord = n_elements(thispix)
        coordind = lindgen(ecliplen) mod ncoord
        if (ecliplen gt ncoord) then $
          print, "Needed ", ecliplen, " coords but got ", ncoord, " on tile ", ii
        ;print, "Filling in tile ", ii
        glon = phi[thispix[coordind]]*180./!dpi
        glat = (theta[thispix[coordind]]-!dpi/2.)*180./!dpi
        ; Transform from galactic healpix to ecliptic observations
        euler, glon, glat, elon, elat, select=6
        euler, glon, glat, ra, dec, select=2
        eclip_trial.coord.elon = elon
        eclip_trial.coord.elat = elat
        eclip_trial.coord.ra = ra
        eclip_trial.coord.dec = dec
        eclip_trial.coord.glon = glon
        eclip_trial.coord.glat = glat
        eclip_trial.coord.healpix_n = fnums[ii]
      
        if (ecliplen_tot gt 0) then eclip = struct_append(eclip, eclip_trial) $
        else eclip = eclip_trial
        ecliplen_tot = ecliplen_tot + ecliplen
      end
    end
    ; Survey: figure out npointings and field angles
    eclip_survey, seg, fov, eclip, offset=skirt
    ; Determine FOV index
    fov_ind = intarr(n_elements(eclip))
    fov_ind[where((eclip.coord.fov_r ge 0.104*CCD_PIX) and $ ; was 104
                (eclip.coord.fov_r lt 0.365*CCD_PIX))] = 1   ; was 391
    fov_ind[where((eclip.coord.fov_r ge 0.365*CCD_PIX) and $
                (eclip.coord.fov_r lt 0.592*CCD_PIX))] = 2
    fov_ind[where(eclip.coord.fov_r  ge 0.592*CCD_PIX)] = 3
    eclip.coord.field_angle = eclip.coord.fov_r / CCD_PIX * fov
    eclip.coord.fov_ind=fov_ind

    ; Dilute
    ;print, "Diluting with binary companions"
    ;bindil = where(eclip.class eq 1)
    ;if (bindil[0] ne -1) then $
    ;  dilute_binary, eclip[bindil], targets, frac_fits, rad_fits, ph_fits, aspix=aspix, radmax=4.0
    ;print, "Diluting with other target stars"
    ;targdil = where(eclip.class eq 1 or eclip.class eq 2)
    ;if (targdil[0] ne -1) then $
    ;  dilute_eclipse, eclip, targets, frac_fits, rad_fits, ph_fits, aspix=aspix, sq_deg=13.4, radmax=4.0
    ;print, "Diluting with background stars"
    ;dilute_eclipse, eclip, bkgnds, frac_fits, rad_fits, ph_fits, aspix=aspix, sq_deg=0.134, radmax=4.0
    ;print, "Diluting with deep stars"
    ;dilute_eclipse, eclip, deeps, frac_fits, rad_fits, ph_fits, aspix=aspix, sq_deg=0.0134, radmax=2.0
    ; Observe
    eclip_observe, eclip, targets, bkgnds, deeps, $
       frac_fits, rad_fits, ph_fits, cr_fits, $
       aspix=aspix, effarea=effarea, sys_limit=sys_limit, $ ;infil=sp_name,outfile=spo_name
       readnoise=readnoise, thresh=thresh, tranmin=tranmin, ps_len=ps_len, $
       duty_cycle=duty_cycle[ii], ffi_len=ffi_len, saturation=saturation, $
       subexptime=subexptime, dwell_time=orbit_period, downlink=downlink
    det = where(eclip.det1 or eclip.det2 or eclip.det)
    if (det[0] ne -1) then begin
      detid = eclip[det].hostid
      ndet = n_elements(det)
      bins = targets[detid].pri + 2*targets[detid].sec
      tmp_star = [[eclip[det].trial], [targets[detid].mag.v], [targets[detid].mag.ic], $
                [targets[detid].mag.t], [targets[detid].mag.j], $
                [targets[detid].mag.h], [targets[detid].mag.k], [targets[detid].teff], $
                [eclip[det].coord.elon], [eclip[det].coord.elat], $
                [eclip[det].coord.glon], [eclip[det].coord.glat], $
                [eclip[det].p], [eclip[det].a], [eclip[det].s], [eclip[det].b], $
                [eclip[det].teff2], [eclip[det].m2], [eclip[det].r2], $
                [eclip[det].dep1], [eclip[det].dur1], [eclip[det].neclip_obs1], $
                [eclip[det].teff1], [eclip[det].m1], [eclip[det].r1], $ 
                [eclip[det].dep2], [eclip[det].dur2], [eclip[det].neclip_obs2], $
                [eclip[det].snreclp1], [eclip[det].snrgress1], $
                [eclip[det].snreclp2], [eclip[det].snrgress2], $
                [eclip[det].k], [eclip[det].snrhr], $
	        [eclip[det].star_ph], [eclip[det].bk_ph], [eclip[det].zodi_ph], $
                [eclip[det].npix], [eclip[det].dil], [targets[detid].ffi], [eclip[det].npointings] ,$
                [eclip[det].sat], [eclip[det].coord.fov_r], $
                [eclip[det].class], [eclip[det].sep], [eclip[det].tsys], $
                [bins], [targets[detid].companion.sep], [targets[targets[detid].companion.ind].mag.t]]
      idx = lindgen(ndet) + totdet
      star_out[idx,*] = tmp_star
      totdet = totdet+ndet
      ;stard = star[det]
      ;if (keyword_set(sav)) then save, filen=spo_name, stard
    end
  endfor
  if (totdet gt 0) then mwrfits, star_out[0:(totdet-1),*], outname
  print, numps, ' postage stamps assigned'
END
