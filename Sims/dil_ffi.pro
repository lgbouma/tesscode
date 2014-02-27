function dil_ffi, dur_min, ffi_len, ffi_os=ffi_os, ffis=ffis
    ffi_os = randomu(seed, n_elements(dur_min)) ; how much of the first ffi is spent in transit?
    ;ffi_os = fltarr(n_elements(dur_min))
    nffi = float(floor((dur_min - ffi_len*ffi_os)/ffi_len)) ; these images are fully transited
    igrss = ffi_os ; this fraction of ingress image is transited
    egrss = (dur_min/ffi_len) - nffi - ffi_os ; this fraction of egress image is transitied
    ; Calculate snrs of 4 cases
    snr_ffi = nffi/sqrt(nffi)
    snr_ffi[where(nffi eq 0)] = 0 ; correct the divide by zero
    snr_igr = (nffi+igrss)/sqrt(nffi+1.)
    snr_egr = (nffi+egrss)/sqrt(nffi+1.)
    snr_both = (nffi+igrss+egrss)/sqrt(nffi+2.)
    ; Initial assumption
    dep_dil = 1. + fltarr(n_elements(dur_min))
    ffis = nffi
    add_igr = where((snr_igr gt snr_egr) and (snr_igr gt snr_ffi))
    if (add_igr[0] ne -1) then begin
      dep_dil[add_igr] = (nffi[add_igr]+igrss[add_igr])/(nffi[add_igr]+1.)
      ffis[add_igr] = nffi[add_igr] + 1.
    endif
    add_egr = where((snr_egr gt snr_igr) and (snr_egr gt snr_ffi))
    if (add_egr[0] ne -1) then begin
      dep_dil[add_egr] = (nffi[add_egr]+egrss[add_egr])/(nffi[add_egr]+1.)
      ffis[add_egr] = nffi[add_egr] + 1.
    endif
    add_both = where((snr_both gt snr_igr) and (snr_both gt snr_egr) and (snr_both gt snr_ffi))
    if (add_both[0] ne -1) then begin
      dep_dil[add_both] = (nffi[add_both]+igrss[add_both]+egrss[add_both])/(nffi[add_both]+2.)
      ffis[add_both] = nffi[add_both] + 2.
    endif
    noegr = where(dur_min lt ffi_len*ffi_os)
    if (noegr[0] ne -1) then begin
      print, n_elements(noegr)
      dep_dil[noegr] = dur_min[noegr]/ffi_len
      ffis[noegr] = 1.
    endif
    return, dep_dil
end

