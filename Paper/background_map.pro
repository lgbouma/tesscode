pro background_map,ps=ps

  !p.charsize=2
  !p.multi=[0,2,2,0,1]
  dlon=60 & dlat=30
  th=1

  if (keyword_set(ps)) then begin

     set_plot,'ps'
     device,filen='bgmap.eps',/encapsulated,xsize=8,ysize=4.25,/inches,$
           /color,bits_per_pixel=16,/isolatin,/helvetica
     !p.font=0
     !p.thick=4 & !x.thick=4 & !y.thick=4 & !p.charthick=1 & !p.charsize=0.9
     !y.margin=[4,2]; 4,2 is default
     !x.margin=[10,3]; 10,3 is default
     syms=0.5
     th=0.25
     
  endif

  lines_color = fsc_color('Charcoal')

  m=200
  n=100

  x = -2.0 + 4.0*dindgen(m)/double(m-1)
  y = -1.0 + 2.0*dindgen(n)/double(n-1)

  cts_zodi = dblarr(m,n)
  cts_stars = dblarr(m,n)
  
  for i=0,m-1 do begin
     for j=0,n-1 do begin
        z = sqrt(1.-(0.25*x[i])^2.-(0.5*y[j])^2.)
        r = (0.5*x[i])^2. + (y[j])^2.
        elon = 2.0*atan(z*x[i],2.*(2.*z^2. - 1.0))
        elat = asin(z*y[j])
        elon = elon*180./!PI & elat = elat*180./!PI
        euler, elon, elat, glon, glat, select=5
        cts_zodi[i,j] = calculate_counts_zodi(elon,elat)
        cts_stars[i,j] = calculate_counts_stars(glon,glat)
        if (r gt 1.0) then begin
           cts_zodi[i,j] = 0.0
           cts_stars[i,j] = 0.0
        endif
     endfor
  endfor

  print, min(cts_stars), max(cts_stars), min(cts_zodi), max(cts_zodi)

  z1=alog10(cts_stars)
  z2=alog10(cts_zodi)
  z=alog10(cts_stars+cts_zodi)
  ra = [alog10(1300.),alog10(45.)]

;  plotimage, z1, range=ra, /isotropic, $
;             imgxrange=[-180,180], imgyrange=[-90,90],$
;             xtit='ecliptic latitude', ytit='ecliptic longitude',$
;             title='unresolved stars'
;  aitoff_grid, dlon,dlat,color=lines_color,thick=th
  plotimage, z2, range=ra, /isotropic, $
             imgxrange=[-180,180], imgyrange=[-90,90],$
             xtit='Ecliptic Longitude', ytit='Ecliptic Latitude',$
             title='Zodiacal Light', font_size=4
  aitoff_grid, dlon,dlat,color=lines_color,thick=th
  plotimage, z, range=ra, /isotropic, $
             imgxrange=[-180,180], imgyrange=[-90,90],$
             xtit='Ecliptic Longitude', ytit='Ecliptic Latitude',$
             title='Unresolved Stars + Zodiacal Light', font_size=4
  aitoff_grid, dlon,dlat,color=lines_color,thick=th

;;;;;

  x = -2.0 + 4.0*dindgen(m)/double(m-1)
  y = -1.0 + 2.0*dindgen(n)/double(n-1)

  cts_zodi = dblarr(m,n)
  cts_stars = dblarr(m,n)
  
  for i=0,m-1 do begin
     for j=0,n-1 do begin
        z = sqrt(1.-(0.25*x[i])^2.-(0.5*y[j])^2.)
        r = (0.5*x[i])^2. + (y[j])^2.
        glon = 2.0*atan(z*x[i],2.*(2.*z^2. - 1.0))
        glat = asin(z*y[j])
        glon = glon*180./!PI & glat = glat*180./!PI
        euler, glon, glat, elon, elat, select=6
        cts_zodi[i,j] = calculate_counts_zodi(elon,elat)
        cts_stars[i,j] = calculate_counts_stars(glon,glat)
        if (r gt 1.0) then begin
           cts_zodi[i,j] = 0.0
           cts_stars[i,j] = 0.0
        endif
     endfor
  endfor

  z1=alog10(cts_stars)
  z2=alog10(cts_zodi)
  z=alog10(cts_stars+cts_zodi)

  plotimage, z1, range=ra, /isotropic, imgxrange=[-180,180], imgyrange=[-90,90],$
             xtit='Galactic Longitude', ytit='Galactic Latitude',$
             title='Unresolved Stars', font_size=4
  aitoff_grid, dlon,dlat,color=lines_color,thick=th
; plotimage, z2, range=ra, /isotropic, imgxrange=[-180,180], imgyrange=[-90,90],$
;             xtit='galactic longitude', ytit='galactic latitude',$
;             title='zodiacal light'
;  aitoff_grid, dlon,dlat,color=lines_color,thick=th
  plotimage, z, range=ra, /isotropic, imgxrange=[-180,180], imgyrange=[-90,90],$
             xtit='Galactic Longitude', ytit='Galactic Latitude',$
             title='Unresolved Stars + Zodiacal light', font_size=4
  aitoff_grid, dlon,dlat,color=lines_color,thick=th

;;;;

  for i=0,m-1 do begin
     for j=0,n-1 do begin
        z = sqrt(1.-(0.25*x[i])^2.-(0.5*y[j])^2.)
        r = (0.5*x[i])^2. + (y[j])^2.
        lon = 2.0*atan(z*x[i],2.*(2.*z^2. - 1.0))
        lat = asin(z*y[j])
        lon = lon*180./!PI & lat = lat*180./!PI
        euler, lon, lat, elon, elat, select=3
        euler, lon, lat, glon, glat, select=1
        cts_zodi[i,j] = calculate_counts_zodi(elon,elat)
        cts_stars[i,j] = calculate_counts_stars(glon,glat)
        if (r gt 1.0) then begin
           cts_zodi[i,j] = 0.0
           cts_stars[i,j] = 0.0
        endif
     endfor
  endfor

;  z1=alog10(cts_stars)
;  z2=alog10(cts_zodi)
;  z=alog10(cts_stars+cts_zodi)

;  plotimage, z1, range=ra, /isotropic, imgxrange=[-180,180], imgyrange=[-90,90],$
;             xtit='celestial longitude', ytit='celestial latitude',$
;             title='unresolved stars'
;  aitoff_grid, dlon,dlat,color=lines_color,thick=th
;  plotimage, z2, range=ra, /isotropic, imgxrange=[-180,180], imgyrange=[-90,90],$
;             xtit='celestial longitude', ytit='celestial latitude',$
;             title='zodiacal light'
;  aitoff_grid, dlon,dlat,color=lines_color,thick=th
;  plotimage, z, range=ra, /isotropic, imgxrange=[-180,180], imgyrange=[-90,90],$
;             xtit='celestial longitude', ytit='celestial latitude',$
;             title='unresolved stars + zodiacal light'
;  aitoff_grid, dlon,dlat,color=lines_color,thick=th

  if (keyword_set(ps)) then begin

     device,/close

     set_plot,'x'
     !p.font=-1
     !p.thick=1 & !x.thick=1 & !y.thick=1 & !p.charthick=1 & !p.charsize=1

  endif

end

function calculate_counts_zodi, elon, elat

  vmax = 23.3447
  dv = 1.148

  x = 90.0-abs(elat)
  vmag = vmax - dv*(x/90.0)^2.
  ;cts = 70.622 * 10.0^(-0.4*(vmag-22.8))
  cts = 76.957 * 10.0^(-0.4*(vmag-22.8))

  return, cts

end

function calculate_counts_stars, glon, glat

;  p = [18.9733D0,8.833D0,4.007D0,0.805D0]
  p = [18.66, 4.28, 0.523, 10.18, 0.459, -3.74]
  dlon = glon
  q = where(dlon gt 180.) & if(q[0] ne -1) then dlon[q] = 360.-dlon[q]
  dlon = abs(dlon)/180.0D0
  ;print, 'Using new model'
  ;p = bk_p[fov_ind,*]
  dlat = abs(glat)/90.0D0
  imag = p[0] + $
    p[1]*(1.0-exp(-dlon/p[2])) + $
    p[3]*(1.0-exp(-dlat/p[4])) + $
    p[5]*sqrt(dlon*dlat)

;  q = where(glon gt 180D0)
;  if (q[0] ne -1) then glon[q]=360D0-glon[q]
;  dglat = abs(glat)/40.0D0
;  dglon = (abs(glon)/180.D0)^(p[3])

;  imag = p[0] + p[1]*dglat + p[2]*dglon
 
  cts = 1.7D6 * 10.0^(-0.4*imag) * 73.0 * 20.^2. ; ph/s/pixel in TESS band
  q = where(cts lt 1.0)
  if (q[0] ne -1) then cts[q]=1.

  return, cts

end

;pro test
;  
;  n=1000
;  lon = 0.0
;  lat = 90.*dindgen(n)/double(n-1)
;
;  cts = calculate_counts_stars(lon,lat)
;  plot, lat, cts, /ylog, yra=[0.5,1E3]
;
;  for i=0,n-1 do begin
;     lon = double(i)*180.
;     cts = calculate_counts_stars(lon,lat)
;     oplot, lat, cts
;  endfor
;
;end
