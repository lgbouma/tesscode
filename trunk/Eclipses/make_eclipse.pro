PRO make_eclipse, sstruct, estruct, frac, rad, ph_p
  ; Add Planets to all target stars
  add_planets, sstruct, p_eclip, frac, rad, ph_p
  estruct = p_eclip
  ; Add EBs (identify EBs among target stars)

  ; Add HEBs - TODO

  ; Add BEBs (identify EBs among background stars, attach to target stars)
  
END
