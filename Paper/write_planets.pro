PRO write_planets, filen=filen
  nstars=1E6     
  template_planet = {$
                    n: 0, $   ; planets per star (set to one for now)
                    r: 0.0, $ ; can eventually replace by dblarr(5), $
                    m: 0.0, $ ; mass
                    p: 0.0, $ ; Period (days)
                    a: 0.0, $ ; Semimajor axis (AU)
                    s: 0.0, $ ; incident flux in units of Sun-->Earth flux
                    b: 0.0, $ ; Impact parameter (0-1)
                    k: 0.0, $ ; rv amplitude
                    tra: 0, $ ; Transit? (boolean)
                    dep: 0.0, $ ; Transit depth (0-1)
                    dur: 0.0, $ ; Transit duration (days)
                    ntra_obs: 0, $ ; Number of transits observed
                    det: 0, $  ; Detected?
                    snr: 0.0, $ ; SNR in phase-folded lightcurve
                    snrtran: 0.0, $ ; SNR per transit
                    durpar: 0.0, $    ; duration of in+engress (days)
                    snrgress: 0.0, $   ; SNR of ingress/egress
                    hostid: 0L $
                    }


     period_boundary = [0.8, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 245.0, 418.0]
     radius_boundary = [0.8, 1.25, 2.0, 4.0, 6.0, 22.0]
     planet_type = ['Earths', 'Super-Earths', 'Small Neptunes', 'Large Neptunes', 'Giants']

     rate_fressin = dblarr(11,5) ; period bin, radius bin
     rate_fressin[0,*] = [0.18, 0.17, 0.035, 0.004, 0.015]
     rate_fressin[1,*] = [0.61, 0.74, 0.18,  0.006, 0.067]
     rate_fressin[2,*] = [1.72, 1.49, 0.73,  0.11,  0.17]
     rate_fressin[3,*] = [2.70, 2.90, 1.93,  0.091, 0.18]
     rate_fressin[4,*] = [2.70, 4.30, 3.67,  0.29,  0.27]
     rate_fressin[5,*] = [2.93, 4.49, 5.29,  0.32,  0.23]
     rate_fressin[6,*] = [4.08, 5.29, 6.45,  0.49,  0.35]
     rate_fressin[7,*] = [3.46, 3.66, 5.25,  0.66,  0.71]
     rate_fressin[8,*] = [0.0,  6.54, 4.31,  0.43,  1.25]
     rate_fressin[9,*] = [0.0,  0.0,  3.09,  0.53,  0.94]
     rate_fressin[10,*] =[0.0,  0.0,  0.0,   0.24,  1.05]
     rate_fressin = rate_fressin/100.

     nplanets = total(round(rate_fressin * nstars))
     ; Pre-allocate for speed
     planet = replicate(template_planet, total(nplanets))
     idx0 = 0L
     for i=10,0,-1 do begin
        for j=4,0,-1 do begin
           binplanets = round(rate_fressin[i,j] * nstars)
           if (binplanets gt 0) then begin
              tmp_planet = replicate(template_planet, binplanets)
              tmp_planet.hostid = floor(double(nstars)*randomu(seed, binplanets))
              if (keyword_set(dressing)) then begin
                periods = fltarr(binplanets) + 10^((alog10(period_boundary[i])+alog10(period_boundary[i+1]))/2.0)
                radii   = fltarr(binplanets) + 10^((alog10(radius_boundary[j])+alog10(radius_boundary[j+1]))/2.0)
              endif else begin
                if(j eq 0) then radpow = 0.0 else radpow = -1.7
                randomp, periods, -1.0, binplanets, range_x = [period_boundary[i], period_boundary[i+1]], seed=seed
                randomp, radii, radpow, binplanets, range_x = [radius_boundary[j], radius_boundary[j+1]], seed=seed
              endelse
              tmp_planet.r = radii
              tmp_planet.p = periods
              idx = lindgen(binplanets) + idx0
              ;print, i, j, binplanets, n_elements(tmp_planet), n_elements(idx), max(idx)/float(total(nplanets))
              planet[idx] = tmp_planet
              delvar, tmp_planet
              idx0 = max(idx) + 1
           endif
        endfor
     endfor
     if (keyword_set(filen)) then mwrfits, [[planet.r], [planet.p]], filen
end
