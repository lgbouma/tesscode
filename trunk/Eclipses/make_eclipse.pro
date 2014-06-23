function make_eclipse, sstruct, bkstruct, estruct, frac, rad, ph_p, dartstruct, tefftic, $
  eclass, min_depth=min_depth
  ecliplen = 0L
  if (keyword_set(min_depth)) then min_depth=min_depth else min_depth = 1D-5
  
  ; Add Planets to all target stars
  if (eclass[0]) then begin
    gd = add_planets(sstruct, p_eclip, frac, rad, ph_p, min_depth=min_depth)
    if (gd gt 0) then estruct = p_eclip
    ecliplen = ecliplen + gd
  endif

  ; Add EBs (identify EBs among target stars)
  if (eclass[1]) then begin
    gd = add_ebs(sstruct, eb_eclip, frac, rad, ph_p)
    if (gd gt 0) then begin
      if (ecliplen gt 0) then estruct = struct_append(estruct, eb_eclip) $
      else estruct = eb_eclip   
    endif
    ecliplen = ecliplen + gd
  endif
  
  ; Add BEBs (identify EBs among background stars, attach to target stars)
  if (eclass[2]) then begin
    gd = add_bebs(sstruct, bkstruct, beb_eclip, frac, rad, ph_p, 100)
    if (gd gt 0) then begin
      if (ecliplen gt 0) then estruct = struct_append(estruct, beb_eclip) $
      else estruct = beb_eclip   
    endif
    ecliplen = ecliplen + gd
  endif
  
  ; Add HEBs
  if (eclass[3]) then begin
    gd = add_hebs(sstruct, heb_eclip, frac, rad, ph_p, dartstruct, tefftic)
    if (gd gt 0) then begin
      if (ecliplen gt 0) then estruct = struct_append(estruct, heb_eclip) $
      else estruct = heb_eclip   
    endif
    ecliplen = ecliplen + gd
  endif

  return, ecliplen
END
