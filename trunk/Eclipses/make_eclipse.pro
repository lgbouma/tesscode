PRO make_eclipse, fstub, eclip=eclip
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
    eclp_observe, sstruct=star, pstruct=eclip
 
  end
    ; Add Planets
    add_planets, sstruct=star, pstruct=p_eclip
    ; Add EBs

    ; Add HEBs
  
    ; Add BEBs

    ; Combine the fluxes

    ; Dilute fluxes

END
