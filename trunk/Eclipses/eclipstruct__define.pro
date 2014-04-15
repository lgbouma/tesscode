PRO eclipstruct__define

  eclip = {eclipstruct,  $
                    class: 0, $ ;1 planet around field star, 2=planet around binary, 3=heb, 4=beb/eb
                    r1: 0.0, $ ; radius of primary
                    r2: 0.0, $ ; radius of secondary
                    m1: 0.0, $ ; mass of primary
                    m2: 0.0, $ ; mass of secondary
                    k: 0.0, $
                    teff1: 0.0, $ ; temp of primary
                    teff2: 0.0, $ ; temp of secondary           
                    p: 0.0, $ ; Period (days)
                    a: 0.0, $ ; Semimajor axis (AU)
                    s: 0.0, $ ; Insolation of secondary
		    cosi: 0.0, $ ; -1 to 1
                    b: 0.0, $ ; Impact parameter (0-1)
                    dep1: 0.0, $ ; Eclipse depth (0-1)
                    dep1_eff: 0.0, $ ; Effective depth
                    dep2: 0.0, $ ; Eclipse depth (0-1)
                    dep2_eff: 0.0, $ ; Effective secondary eclipse depth
                    dur1: 0.0, $ ; Primary eclipse duration (days)
                    dur2: 0.0, $ ; Primary eclipse duration (days)
                    dur1_eff: 0.0, $ ; Effective duration
                    dur2_eff: 0.0, $ ; Effective duration
                    neclip_obs1: 0, $ ; Number of primary eclipses observed
                    neclip_obs2: 0, $ ; Number of secondary eclipses observed
                    snr: 0.0, $ ; SNR of primary eclipses in phase-folded lightcurve
                    snr1: 0.0, $ ; SNR of primary eclipses in phase-folded lightcurve
                    snr2: 0.0, $ ; SNR of primary eclipses in phase-folded lightcurve
                    snreclp1: 0.0, $ ; SNR per primary eclipse
                    snreclp2: 0.0, $ ; SNR per secondary eclipse
                    gress1: 0.0, $    ; duration of in+engress (days)
                    gress2: 0.0, $    ; duration of in+engress (days)
                    snrgress1: 0.0, $   ; SNR of ingress/egress
                    snrgress2: 0.0, $   ; SNR of ingress/egress
                    det: 0, $  ; Detected?
                    det1: 0, $  ; Detected primary?
                    det2: 0, $  ; Detected secondary?
                    hostid: 0L $
                    }

end
