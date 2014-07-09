PRO sav2fits, fname, vmag=vmag, newfname=newfname

  restore, fname
  m = star.m
  rad = star.r
  teff = star.teff
  v = star.mag.v
  r = star.mag.r
  ic = star.mag.ic
  z = star.mag.z
  j = star.mag.j
  h = star.mag.h
  k = star.mag.k
  dm = star.mag.dm
  av = star.mag.dm
  ps = star.pri + 2*star.sec
  mv = star.mag.mv
  mic = star.mag.mic
  mj = star.mag.mj
  icsys = star.mag.icsys
  jsys = star.mag.jsys

  poop = [[m],[rad],[teff],[v],[r],[ic],[z],[j],[h],[k],[dm],[av],[ps],[mv],[mic],[mj],[icsys],[jsys]]
  if (keyword_set(vmag)) then begin
    gd = where(v le vmag)
    poop = poop[gd,*]
  end
  if (keyword_set(newfname)) then newfname=newfname else $ 
  newfname = repstr(fname, 'sav', 'fits')
  mwrfits, poop, newfname

END
