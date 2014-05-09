PRO make_eclipse, fstub
  ; Gather the .sav files
  fnames = file_search(fstub)
  numfil = n_elements(fnames)
  numstar = intarr(numfil)
  for ii=0, numfil-1 do begin
    restore, fnames[ii]
    numstar[ii] = n_elements(star)
    ; Add Planets

    ; Add EBs

    ; Add HEBs
  
    ; Add 
  end
  ;


END
