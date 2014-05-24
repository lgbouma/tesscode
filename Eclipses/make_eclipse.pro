PRO make_eclipse, sstruct, estruct, frac, rad, ph_p, min_depth=min_depth
  if (keyword_set(min_depth)) then min_depth=min_depth else min_depth = 1D-5
  ; Add Planets to all target stars
  add_planets, sstruct, p_eclip, frac, rad, ph_p
  estruct = p_eclip[where(p_eclip.dep1 gt min_depth)]
  ; Add EBs (identify EBs among target stars)

  ; Add HEBs - TODO

  ; Add BEBs (identify EBs among background stars, attach to target stars)
  
END
