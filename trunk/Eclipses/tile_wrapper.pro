PRO tile_wrapper, fstub, eclip=eclip
  filen = '../Sims/bs3_13x24.sav'
  prf_file = '../ExpTimeCalc/bigfrac24_105_f3p35.fits' ;'+psfstr+'.fits' 
  ph_file = '../ExpTimeCalc/ph_filt.fits' ;'+psfstr+'.fits' 
  fov = 24.
  seg = 13
  geomarea = 69.1 ;78.6 ;69.1 ;78.6 ;75.8 ;57.6 ;67.5
  readnoise= 10.0
  thresh = 7.0
  tranmin = 2.0
  sys_limit=60.
  n_trial = 25
  ffi_len=2. ; in minutes
  duty_cycle=100.
  nps = 100000
  REARTH_IN_RSUN = 0.0091705248

  ; Gather the .sav files
  fnames = file_search(fstub)
  numfil = n_elements(fnames)
  numstar = intarr(numfil)
  for ii=0, numfil-1 do begin
    restore, fnames[ii]
    numstar[ii] = n_elements(star)
    ; Add eclipses
    make_eclipse, sstruct=star, estruct=eclip
    ; Observe
    eclp_observe, sstruct=star, pstruct=eclip, $
       geomarea=geomarea, fov=fov, sys_limit=sys_limit, $ ;infil=sp_name,outfile=spo_name
        readnoise=readnoise, thresh=thresh, tranmin=tranmin, $
        prf_file=prf_file, bk_file=bk_file, sp_file=sp_file, ph_file=ph_file, $
        duty_cycle=duty_cycle, ffi_len=ffi_len
 
  end
END
