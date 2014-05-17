PRO tile_wrapper, fpath, fnums, eclip=eclip
  frac_file = 'bigfrac24_105_f3p33.fits' ;'+psfstr+'.fits' 
  rad_file = 'bigrad24_105_f3p33.fits' ;'+psfstr+'.fits' 
  ph_file = 'ph_filt_t.fits' ;'+psfstr+'.fits' 
  fov = 24.
  seg = 13
  skirt=6.
  geomarea = 69.1 ;78.6 ;69.1 ;78.6 ;75.8 ;57.6 ;67.5
  aspix = 21.1
  readnoise= 10.0
  thresh = 7.0
  tranmin = 2.0
  sys_limit=60.
  n_trial = 5
  ffi_len=2. ; in minutes
  duty_cycle=100.
  nps = 100000
  REARTH_IN_RSUN = 0.0091705248
  CCD_PIX = 4096.

  ; Gather the .sav files
  numfil = n_elements(fnums)
  numtargets = lonarr(numfil)
  numbkgnd = lonarr(numfil)
  numdeeps = lonarr(numfil)
  ; Open the fits files
  frac_fits = mrdfits(frac_file)
  rad_fits = mrdfits(rad_file)/60. ; put into pixels
  ph_fits = mrdfits(ph_file)
  ; Make random spherical coords
  u = randomu(seed, 1D6)
  v = randomu(seed, 1D6)
  phi = 2.*!dpi*u
  theta = acos(2.*v-1.)
  ang2pix_ring, 16, theta, phi, ipring
  ;cgHistoPlot, ipring
  for ii=0, numfil-1 do begin
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
    for jj=0,n_trial-1 do begin
      ; Add eclipses
      make_eclipse, targets, eclip_trial, frac_fits, rad_fits, ph_fits
      eclip_trial.trial = jj + 1
      ; Add coordinates to the eclipses
      neclip = n_elements(eclip_trial)
      thispix = where(ipring eq fnums[ii])
      ncoord = n_elements(thispix)
      if (neclip gt ncoord) then begin
        print, "Need ", neclip, " coords but got ", ncoord, " on tile ", ii
      endif else begin
        ;print, "Filling in tile ", ii
        glon = phi[thispix[0:(neclip-1)]]*180./!dpi
        glat = (theta[thispix[0:(neclip-1)]]-!dpi)*180./!dpi
        ; Transform from galactic healpix to ecliptic observations
        euler, glon, glat, elon, elat, select=5
        euler, glon, glat, ra, dec, select=2
        eclip_trial.coord.elon = elon
        eclip_trial.coord.elat = elat
        eclip_trial.coord.ra = ra
        eclip_trial.coord.dec = dec
        eclip_trial.coord.glon = glon
        eclip_trial.coord.glat = glat
        eclip_trial.coord.healpix_n = fnums[ii]
      end
      if (jj gt 0) then eclip = struct_append(eclip, eclip_trial) $
      else eclip = eclip_trial
    end
    ; Survey
    eclip_survey, seg, fov, eclip, offset=skirt
    ; Determine FOV index
    fov_ind = intarr(n_elements(eclip))
    fov_ind[where((eclip.coord.fov_r ge 0.104*CCD_PIX) and $
                (eclip.coord.fov_r lt 0.365*CCD_PIX))] = 1
    fov_ind[where((eclip.coord.fov_r ge 0.365*CCD_PIX) and $
                (eclip.coord.fov_r lt 0.592*CCD_PIX))] = 2
    fov_ind[where(eclip.coord.fov_r  ge 0.592*CCD_PIX)] = 3
    eclip.coord.field_angle = eclip.coord.fov_r / CCD_PIX * FOV
    eclip.coord.fov_ind=fov_ind

    ; Dilute
    print, "Diluting with binary companions"
    dilute_binary, eclip[where(eclip.class eq 1)], targets, frac_fits, rad_fits, ph_fits, aspix=aspix
    print, "Diluting with other target stars"
    dilute_eclipse, eclip, targets, frac_fits, rad_fits, ph_fits, aspix=aspix, sq_deg=13.4
    print, "Diluting with background stars"
    dilute_eclipse, eclip, bkgnds, frac_fits, rad_fits, ph_fits, aspix=aspix, sq_deg=0.134
    print, "Diluting with deep stars"
    dilute_eclipse, eclip, deeps, frac_fits, rad_fits, ph_fits, aspix=aspix, sq_deg=0.0134
    stop
    ; Observe
    eclp_observe, sstruct=star, pstruct=eclip, $
       geomarea=geomarea, fov=fov, sys_limit=sys_limit, $ ;infil=sp_name,outfile=spo_name
        readnoise=readnoise, thresh=thresh, tranmin=tranmin, $
        prf_file=prf_file, bk_file=bk_file, sp_file=sp_file, ph_file=ph_file, $
        duty_cycle=duty_cycle, ffi_len=ffi_len
    det = where(eclip.det1 or eclip.det2 or eclip.det)
    detid = eclip[det].hostid
    ndet = n_elements(det)
    tmp_star = [[eclip.trial], [star[detid].mag.v], [star[detid].mag.ic], $
                [star[detid].mag.t], $
                [star[detid].mag.j], [star[detid].mag.k], $
                [star[detid].coord.elon], [star[detid].coord.elat], $
                [star[detid].coord.glon], [star[detid].coord.glat], $
                [star[detid].r], [star[detid].m], [star[detid].teff], $
                [planet[det].r], [planet[det].p], [planet[det].a], [planet[det].s], [planet[det].b], $
                [planet[det].dep], [planet[det].dur], [planet[det].ntra_obs], [star[detid].npointings], $
                [planet[det].snrtran], [planet[det].snrgress], $
                [planet[det].m], [planet[det].k], $
                [star[detid].npix], [star[detid].dil], $
                [star[detid].sat], [star[detid].coord.fov_r], [bins], $
                [star[detid].companion.sep], [star[star[detid].companion.ind].mag.ic], $
                [star[detid].ffi]]
    idx = lindgen(ndet) + totdet
    star_out[idx,*] = tmp_star
    totdet = totdet+ndet
    print, 'Detected ', ndet, ' planets on trial ', ii
    ;stard = star[det]
    ;if (keyword_set(sav)) then save, filen=spo_name, stard
  endfor
  mwrfits, star_out[0:(totdet-1),*], fname 
END
