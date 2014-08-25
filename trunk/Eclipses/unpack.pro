PRO unpack, fname

dat = mrdfits(fname)

trial = dat[*,0]	 ; 
dep1  = dat[*,19]	 ; 
dep2  = dat[*,25]	 ;  
dil   = dat[*,38]	 ; 
teff = dat[*,7]	 	 ; 

dep1_obs = dep1/dil	 ;
dep2_obs = dep2/dil	 ;

; Weed out whopping EBs and very hot stars/WDs
gd = where(dep1_obs lt 0.1 and dep2_obs lt 0.1 and teff lt 15000.)
dat = dat[gd,*]

trial  = dat[*,0]	 ; Trial number
vmag   = dat[*,1]	 ; Apparent V
icmag  = dat[*,2]	 ; Apparent Ic
tmag   = dat[*,3]	 ; Apparent TESS mag
jmag   = dat[*,4]	 ; Apparent 2MASS J
hmag   = dat[*,5]	 ; Apparent 2MASS H
kmag   = dat[*,6]	 ; Apparent 2MASS Ks
teff   = dat[*,7]	 ; Teff of target star
elon   = dat[*,8]	 ; ecliptic lon
elat   = dat[*,9]	 ; ecliptic lat
ra    = dat[*,10]	 ; RA
dec   = dat[*,11]	 ; DEC
p     = dat[*,12]	 ; Period [days]
a     = dat[*,13]	 ; Semimajor axis [AU]
s     = dat[*,14]	 ; Insolation [Sun->Earth]
b     = dat[*,15]	 ; Impact Parameter
teff2 = dat[*,16]	 ; Temperature of secondary body
m2    = dat[*,17]*332946. ; Mass of secondary body (Msun converted to Mearth)
r2    = dat[*,18]/0.00917 ; Radius of secondary (Rsun converted to Rearth)
dep1  = dat[*,19]	 ; Depth of primary eclipse (transit)
dur1  = dat[*,20]	 ; Duration of primary eclipse (transit)
necl1 = dat[*,21]	 ; Number of primary eclipses (transits) observed
teff1 = dat[*,22]	 ; Teff of primary
m1    = dat[*,23]	 ; Mass of primary
r1    = dat[*,24]	 ; Radius of primary (keep in Rsun)
dep2  = dat[*,25]	 ; Depth of secondary eclipse
dur2  = dat[*,26]	 ; Duration of secondary eclipse
necl2 = dat[*,27]	 ; Number of secondary eclipses
snr1  = dat[*,28]	 ; SNR of phase-folded primary eclipses
gress1= dat[*,29]	 ; SNR of folded primary ingress/egress
snr2  = dat[*,30]	 ; SNR of phase-folded secondary eclipses
gress2= dat[*,31]	 ; SNR of folded secondary ingress/egress
rvk   = dat[*,32]	 ; Radial velocity semi-amplitude
snrh  = dat[*,33]	 ; SNR per hour
starph= dat[*,34]	 ; photon flux from star
bkph  = dat[*,35]	 ; background photon flux
zodi  = dat[*,36]	 ; zodi photon flux
npix  = dat[*,37]	 ; Number of pixels in photometric aperture
dil   = dat[*,38]	 ; Dilution
ffi   = dat[*,39]	 ;
npts  = dat[*,40]	 ; 
sat   = dat[*,41]	 ; 
fovr  = dat[*,42]	 ; 
eclass= dat[*,43]	 ; 
hsep  = dat[*,44]	 ;
icsys = dat[*,45]	 ;
tsys  = dat[*,46]	 ;
jsys  = dat[*,47]	 ;
censhift1 = dat[*,48]	 ;
censhift2 = dat[*,49]	 ;
cenerr1   = dat[*,50]	 ;
cenerr2   = dat[*,51]	 ;
var  = dat[*,52]	 ;
bin  = dat[*,53]	 ; 
binsep = dat[*,54]	 ; 
bint   = dat[*,55]	 ; 
dep1_obs = dep1./dil	 ;
dep2_obs = dep2./dil	 ;
r2_obs = sqrt(dep1./dil).*r1/0.00917	 ;
r1_obs = sqrt(dep2./dil).*r2/0.00917	 ;
snr1f = snr1.*sqrt(necl1)	 	;
snr2f = snr2.*sqrt(necl2)	 	;
snrf = sqrt(snr1f^2. + snr2f^2.)	 ;

snrgress1 = snr1./sqrt(6*gress1.*dur1)	 ;
snrgress2 = snr2./sqrt(6*gress2.*dur2)	 ;

ddep = abs(dep1-dep2)	 ;
ddsnr = ddep/(1./snr1+1./snr2)	 ;

stop

END
