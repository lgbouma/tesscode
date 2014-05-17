PRO starstruct__define

  template_star = {starstruct, $
                  r: 0.0, $     ; Radius (solar)
                  ;oldr: 0.0, $     ; Radius (solar)
                  m: 0.0, $     ; Mass (solar)
                  teff: 0.0, $  ; Effective T (Kelvin)
                  ;oldteff: 0.0, $  ; Effective T (Kelvin)
                  cosi: 0.0, $  ; cos inclination of planets
		  logage: 0.0, $
		  mini: 0.0, $
		  feh: 0.0, $
 	          logg: 0.0, $
                  npointings: 0, $ ; Number of TESS pointings this star gets
                  mag: {magstruct}, $
                  companion: {compstruct}, $
                  snr: 0.0, $      ; SNR per hour
                  pri: 0, $        ; Primary of binary?
                  sec: 0, $        ; Secondary of binary? 
                  spl: 0, $        ; Split (for triples and quadruples)
                  ffi: 0 $         ; Full-frame only (default is postage stamp)
                }

end
