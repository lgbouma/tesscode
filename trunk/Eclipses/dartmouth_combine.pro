PRO dartmouth_combine, fpath, starstruct=starstruct
 ages = [0.250, 0.299, 0.349, 0.400, 0.450, 0.500, 0.550, 0.599, 0.650, 0.699, $
         0.750, 0.800, 0.849, 0.900, 0.949, 1.000, 1.250, 1.500, 1.750, 2.000, $ 
         2.250, 2.500, 2.750, 3.000, 3.250, 3.500, 3.750, 4.000, 4.250, 4.500, $
	 4.750, 5.000, 5.500, 6.000, 6.500, 7.000, 7.500, 8.000, 8.500, 9.000, $
         9.500, 10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 13.00, 13.50, 14.00, 14.50, 15.00]
 nage = n_elements(ages)
 masses = (findgen(143) + 20.)/200.
 nmass = n_elements(masses)
 feh = [-1.0, -0.5, 0.0, 0.2]
 nfe = n_elements(feh)
; nmass = 265L
 starstruct = replicate({dartstruct}, nfe*nage*nmass)
 starstruct = reform(starstruct, nmass, nfe, nage) 
 for fi = 0,nfe-1 do begin
   if (feh[fi] lt 0) then fes='m' else fes='p' 
   for ai = 0,nage-1 do begin     
     fstub =  fpath +'a'+string(ages[ai]*1E3, format='(I05)')+'feh'+fes+string(abs(feh[fi]*10),'(I02)')+'*'
     fname = file_search(fstub)
     if ((n_elements(fname) gt 1) or (fname eq '')) then print, 'wrong file name from '+fstub else begin
       ;  print, 'Found '+fname+' from '+fstub
       ;  #EEP   M/Mo    LogTeff  LogG   LogL/Lo U       B       V       R       I       J       H       Ks      Kp      D51
       readcol, fname, eep, mass, logT, logG, logL, U, B, V, R, Ic, J, H, Ks, Kp, D51
       
 ;      starstruct[*,fi,ai].age = ages[ai]
 ;      starstruct[*,fi,ai].feh = feh[fi]
 ;      starstruct[*,fi,ai].m = mass[0:nmass-1]
 ;      starstruct[*,fi,ai].teff = 10.^logT[0:nmass-1]
 ;      starstruct[*,fi,ai].rad = sqrt(mass[0:nmass-1])/sqrt((10.^logG[0:nmass-1])/27542.3)
 ;      starstruct[*,fi,ai].logL = logL[0:nmass-1]
 ;      starstruct[*,fi,ai].u  = U[0:nmass-1]
 ;      starstruct[*,fi,ai].b  = B[0:nmass-1]
 ;      starstruct[*,fi,ai].v  = V[0:nmass-1]
 ;      starstruct[*,fi,ai].r  = R[0:nmass-1]
 ;      starstruct[*,fi,ai].ic = Ic[0:nmass-1]
 ;      starstruct[*,fi,ai].j  = J[0:nmass-1]
 ;      starstruct[*,fi,ai].h  = H[0:nmass-1]
 ;      starstruct[*,fi,ai].ks = Ks[0:nmass-1]
       
       starstruct[*,fi,ai].age = ages[ai]
       starstruct[*,fi,ai].feh = feh[fi]
       starstruct[*,fi,ai].m = masses
       starstruct[*,fi,ai].teff = 10.^(interpol(logT, mass, masses))
       starstruct[*,fi,ai].rad = interpol(sqrt(mass)/sqrt((10.^logG)/27542.3), mass, masses)
       starstruct[*,fi,ai].logL = interpol(logL, mass, masses)
       starstruct[*,fi,ai].u  = interpol(U, mass, masses)
       starstruct[*,fi,ai].b  = interpol(B, mass, masses)
       starstruct[*,fi,ai].v  = interpol(V, mass, masses)
       starstruct[*,fi,ai].r  = interpol(R, mass, masses)
       starstruct[*,fi,ai].ic = interpol(Ic, mass, masses)
       starstruct[*,fi,ai].j  = interpol(J, mass, masses)
       starstruct[*,fi,ai].h  = interpol(H, mass, masses)
       starstruct[*,fi,ai].ks = interpol(Ks, mass, masses)
      
       ;flen[fi,ai] = max(mass)

     end ; if file found
   end ; age loop
 end ; feh loop
END	

