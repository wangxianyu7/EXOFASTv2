pro derivepars, ss, logname=logname

G = ss.constants.GMsun/ss.constants.rsun^3*ss.constants.day^2
AU = ss.constants.au/ss.constants.rsun ;; R_sun
mjup = ss.constants.gmjupiter/ss.constants.gmsun ;; m_sun
rjup = ss.constants.rjupiter/ss.constants.rsun  ;; r_sun
mearth = ss.constants.gmearth/ss.constants.gmsun ;; m_sun
rearth = ss.constants.rearth/ss.constants.rsun  ;; r_sun
sigmaB = ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2 ;; Stefan-Boltzmann constant

for i=0L, n_elements(*ss.priors)-1 do begin
   prior = (*ss.priors)[i]
   if prior.gaussian_width ne 0d0 or prior.value[1] eq -1 then continue

   ;; the prior is linked to another variable -- get its value
   if prior.value[4] ne -1 then begin
      value = ss.(prior.value[0])[prior.value[1]].(prior.value[2])[prior.value[3]].(prior.value[4])[prior.value[5]].value
   endif else if prior.value[2] ne -1 then begin
      value = ss.(prior.value[0])[prior.value[1]].(prior.value[2])[prior.value[3]].value
   endif else begin
      value = ss.(prior.value[0])[prior.value[1]].value
   endelse
   
   ;; assign the value
   if prior.map[4] ne -1 then begin
      ss.(prior.map[0])[prior.map[1]].(prior.map[2])[prior.map[3]].(prior.map[4])[prior.map[5]].value = value
   endif else if prior.map[2] ne -1 then begin
      ss.(prior.map[0])[prior.map[1]].(prior.map[2])[prior.map[3]].value = value
   endif else if prior.map[0] ne -1 then begin
      ss.(prior.map[0])[prior.map[1]].value = value
   endif
endfor

for i=0L, ss.nstars-1 do begin
   ss.star[i].mstar.value = 10^ss.star[i].logmstar.value
   ss.star[i].rhostar.value = ss.star[i].mstar.value/(ss.star[i].rstar.value^3)*ss.constants.rhosun                                          ;; rho_sun
   ss.star[i].logg.value = alog10(ss.star[i].mstar.value/(ss.star[i].rstar.value^2)*ss.constants.gravitysun)                                 ;; cgs
   ss.star[i].lstar.value = 4d0*!dpi*ss.star[i].rstar.value^2*ss.star[i].teff.value^4*ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2 ;; lsun
   ss.star[i].fbol.value = (ss.star[i].lstar.value*ss.constants.lsun)/(4d0*!dpi*(ss.star[i].distance.value*ss.constants.pc)^2)                  ;; cgs

   ;; derive the age
   if ss.mist[i] and not ss.star[i].age.fit and ss.nstars eq 1 then begin
      for j=0L, ss.nsteps-1 do begin
         mistchi2 = massradius_mist(ss.star[i].eep.value[j],ss.star[i].mstar.value[j],$
                                    ss.star[i].initfeh.value[j],ss.star[i].age.value[j],$
                                    ss.star[i].teff.value[j],ss.star[i].rstar.value[j],$
                                    ss.star[i].feh.value[j],mistage=mistage,fitage=ss.star[i].age.fit,$
                                    tefffloor=ss.teffemfloor,fehfloor=ss.fehemfloor,$
                                    rstarfloor=ss.rstaremfloor, agefloor=ss.ageemfloor)
         
         ss.star[i].age.value[j] = mistage
      endfor
   endif else if ss.parsec[i] and not ss.star[i].age.fit and ss.nstars eq 1 then begin
      for j=0L, ss.nsteps-1 do begin
         parsecchi2 = massradius_parsec(ss.star[i].eep.value[j],ss.star[i].mstar.value[j],$
                                        ss.star[i].initfeh.value[j],ss.star[i].age.value[j],$
                                        ss.star[i].teff.value[j],ss.star[i].rstar.value[j],$
                                        ss.star[i].feh.value[j],parsec_age=parsec_age,fitage=ss.star[i].age.fit,$
                                        tefffloor=ss.teffemfloor,fehfloor=ss.fehemfloor,$
                                        rstarfloor=ss.rstaremfloor, agefloor=ss.ageemfloor)
         ss.star[i].age.value[j] = parsec_age
      endfor
   endif

   ;; derive the distance
   ;; if user fixes the distance to another star's distance, parallax won't be derived like this
;   if ss.star[i].distance.fit then ss.star[i].parallax.value = 1d3/ss.star[i].distance.value ;; mas
;   if ss.star[i].parallax.fit then ss.star[i].distance.value = 1d3/ss.star[i].parallax.value ;; pc
   ;; do this instead -- hack!
   if stdev(ss.star[i].distance.value) lt 1d-8 then ss.star[i].distance.value = 1d3/ss.star[i].parallax.value
   if stdev(ss.star[i].parallax.value) lt 1d-8 then ss.star[i].parallax.value = 1d3/ss.star[i].distance.value
   
   ss.star[i].absks.value = ss.star[i].appks.value - 2.5d0*alog10((ss.star[i].distance.value/10d0)^2)                       ;; mag

endfor   

for j=0, ss.ntel-1 do begin
   positive = where(ss.telescope[j].jittervar.value gt 0d0)
   ss.telescope[j].jitter.value[positive] = sqrt(ss.telescope[j].jittervar.value[positive])
endfor

;; find the extent of the data
minbjd = !values.d_infinity
maxbjd = -!values.d_infinity
for i=0, ss.ntel-1 do begin
   tmpminbjd = min((*ss.telescope[i].rvptrs).bjd,max=tmpmaxbjd)
   if tmpminbjd lt minbjd then minbjd = tmpminbjd
   if tmpmaxbjd gt maxbjd then maxbjd = tmpmaxbjd
endfor
for i=0, ss.ntran-1 do begin
   tmpminbjd = min((*ss.transit[i].transitptrs).bjd,max=tmpmaxbjd)
   if tmpminbjd lt minbjd then minbjd = tmpminbjd
   if tmpmaxbjd gt maxbjd then maxbjd = tmpmaxbjd
endfor
for i=0, ss.nastrom-1 do begin
   tmpminbjd = min((*ss.astrom[i].astromptrs).bjdtdb,max=tmpmaxbjd)
   if tmpminbjd lt minbjd then minbjd = tmpminbjd
   if tmpmaxbjd gt maxbjd then maxbjd = tmpmaxbjd
endfor

;; sanity check
if (maxbjd-minbjd) gt 1d5 then begin
   printandlog, 'WARNING: These data span more than 270 years!', ss.logname
   printandlog, 'Check for consistency between the timestamps in your input data files.', ss.logname
   printandlog, 'This will almost certainly not end well.', ss.logname
endif

;starbb36 = exofast_blackbody(ss.star[i].teff.value,replicate(3550d0/1d9,ss.nsteps),/wave)
;starbb45 = exofast_blackbody(ss.star[i].teff.value,replicate(4493d0/1d9,ss.nsteps),/wave)
starbb25 = exofast_blackbody(ss.star[i].teff.value,replicate(2500d0/1d9,ss.nsteps),/wave)
starbb50 = exofast_blackbody(ss.star[i].teff.value,replicate(5000d0/1d9,ss.nsteps),/wave)
starbb75 = exofast_blackbody(ss.star[i].teff.value,replicate(7500d0/1d9,ss.nsteps),/wave)

for i=0, ss.nplanets-1 do begin

   ;; derive the mass of the planet
   if ss.planet[i].logmp.fit then ss.planet[i].mpsun.value = 10^ss.planet[i].logmp.value $ ;; m_sun
   else ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value) ;; m_sun
   ss.planet[i].mp.value = ss.planet[i].mpsun.value/mjup ;; m_jupiter
   ss.planet[i].mpearth.value = ss.planet[i].mpsun.value/mearth ;; m_earth
   ss.planet[i].q.value = ss.planet[i].mpsun.value/ss.star[ss.planet[i].starndx].mstar.value                           ;; unitless

   ;; derive the radius of the planet   
   ss.planet[i].rpsun.value = ss.planet[i].p.value*ss.star[ss.planet[i].starndx].rstar.value ;; r_sun
   ss.planet[i].rp.value = ss.planet[i].rpsun.value/rjup ;; r_jupiter
   ss.planet[i].rpearth.value = ss.planet[i].rpsun.value/rearth ;; r_earth

   ;; derive period
   ss.planet[i].period.value = 10^ss.planet[i].logp.value  

   ;; derive eccentricity and argument of periastron based on the parameterization
   if ss.planet[i].secosw.fit and ss.planet[i].sesinw.fit then begin
      ss.planet[i].e.value = ss.planet[i].secosw.value^2 + ss.planet[i].sesinw.value^2
      ss.planet[i].omega.value = atan(ss.planet[i].sesinw.value,ss.planet[i].secosw.value)
   endif else if ss.planet[i].vcve.fit then begin

      ;; what sign of the quadratic solution?
      if ss.planet[i].sign.fit then begin
         ;; we fit for it
         sign = floor(ss.planet[i].sign.value)
      endif else begin
         ;; L implicitly defines it
         lsq = ss.planet[i].lsinw.value^2 + ss.planet[i].lcosw.value^2
         sign = (lsq lt 0.5d0)
      endelse

if 0 then begin
      if ss.planet[i].lsinw2.fit then begin
         ss.planet[i].sign.value = ss.planet[i].sign2.value
         ss.planet[i].lsinw.value = ss.planet[i].lsinw2.value
         flip = where(ss.planet[i].vcve.value gt 1d0 and ss.planet[i].lsinw2.value ge 0d0 and ~sign)
         if flip[0] ne -1 then begin
            ss.planet[i].lsinw.value[flip] = -ss.planet[i].lsinw2.value[flip]
            ss.planet[i].sign.value[flip] = ss.planet[i].sign2.value[flip] + 1d0
         endif
         flip = where(ss.planet[i].vcve.value le 1d0 and ss.planet[i].lsinw2.value lt 0d0 and sign)
         if flip[0] ne -1 then begin
            ss.planet[i].lsinw.value[flip] = -ss.planet[i].lsinw2.value[flip]
            ss.planet[i].sign.value[flip] = ss.planet[i].sign2.value[flip] - 1d0
         endif
      endif
endif else begin
;   ss.planet[i].sign.value = ss.planet[i].sign2.value
;   ss.planet[i].lsinw.value = ss.planet[i].lsinw2.value
endelse

      ;; based on L*cos(omega) and L*sin(omega)
      if ss.planet[i].lsinw.fit and ss.planet[i].lcosw.fit then begin
         ss.planet[i].omega.value = atan(ss.planet[i].lsinw.value, ss.planet[i].lcosw.value)
      endif

      ss.planet[i].e.value = vcve2e(ss.planet[i].vcve.value,omega=ss.planet[i].omega.value, sign=floor(ss.planet[i].sign.value))

   endif else if ss.planet[i].e.fit and ss.planet[i].omega.fit then begin
      ;; do nothing
   endif else if ss.planet[i].ecosw.fit and ss.planet[i].esinw.fit then begin
      ss.planet[i].e.value = sqrt(ss.planet[i].ecosw.value^2 + ss.planet[i].esinw.value^2)
      ss.planet[i].omega.value = atan(ss.planet[i].esinw.value,ss.planet[i].ecosw.value)
   endif else if ss.planet[i].qecosw.fit and ss.planet[i].qesinw.fit then begin
      ss.planet[i].e.value = (ss.planet[i].qecosw.value^2 + ss.planet[i].qesinw.value^2)^2
      ss.planet[i].omega.value = atan(ss.planet[i].qesinw.value,ss.planet[i].qecosw.value)
   endif

   zero = where(ss.planet[i].e.value eq 0d0,complement=nonzero)
   if zero[0] ne -1 then ss.planet[i].omega.value[zero] = !dpi/2d0 

   ; ss.planet[i].lambda.value = atan(ss.planet[i].lsinlambda.value,ss.planet[i].lcoslambda.value)
   ss.planet[i].lambda.value = atan(ss.planet[i].svsinicoslambda.value,ss.planet[i].svcosicoslambda.value)
   ss.planet[i].lambdadeg.value = ss.planet[i].lambda.value*180d0/!dpi
   ss.planet[i].vsini.value = (ss.planet[i].svsinicoslambda.value^2 + ss.planet[i].svcosicoslambda.value^2)^0.5
   ss.planet[i].omega.value = atan(ss.planet[i].lsinomega.value,ss.planet[i].lcosomega.value)
   ss.planet[i].omegadeg.value = ss.planet[i].omega.value*180d0/!dpi
   ss.planet[i].bigomega.value = atan(ss.planet[i].lsinbigomega.value,ss.planet[i].lcosbigomega.value)
   ss.planet[i].bigomegadeg.value = ss.planet[i].bigomega.value*180d0/!dpi
   ss.planet[i].esinw.value = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
   ss.planet[i].ecosw.value = ss.planet[i].e.value*cos(ss.planet[i].omega.value)
   ss.planet[i].sesinw.value = sqrt(ss.planet[i].e.value)*sin(ss.planet[i].omega.value)
   ss.planet[i].secosw.value = sqrt(ss.planet[i].e.value)*cos(ss.planet[i].omega.value)
   ss.planet[i].qesinw.value = (ss.planet[i].e.value)^0.25d0*sin(ss.planet[i].omega.value)
   ss.planet[i].qecosw.value = (ss.planet[i].e.value)^0.25d0*cos(ss.planet[i].omega.value)
   ss.planet[i].vcve.value = sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)

   ;; derive a/rstar, a
   ss.planet[i].arsun.value=(G*(ss.star[ss.planet[i].starndx].mstar.value+ss.planet[i].mpsun.value)*ss.planet[i].period.value^2/(4d0*!dpi^2))^(1d0/3d0) ;; (a1 + a2)/rsun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star[ss.planet[i].starndx].rstar.value                                                        ;; (a1 + a2)/rstar
   ss.planet[i].a.value = ss.planet[i].arsun.value*ss.constants.rsun/ss.constants.au                                           ;; AU

;   ;; estimate GR precession (eq 2, Jordan & Bakos, 2008)
;   ss.planet[i].omegagr.value = 7.78d0/(1d0-ss.planet[i].e.value^2)*ss.star[ss.planet[i].starndx].mstar.value*(0.05d0/ss.planet[i].a.value)*(1d0/ss.planet[i].period.value)

   ;; more precise estimate of GR precession (accounts for planet mass) (eq 1, Jordan & Bakos, 2008)
   n = (ss.constants.gmsun*(ss.star[ss.planet[i].starndx].mstar.value+ss.planet[i].mpsun.value)/(ss.planet[i].a.value*ss.constants.au)^3)^(0.5d0) ;; rad/s
   ss.planet[i].omegagr.value = 3d0*ss.constants.gmsun*ss.star[ss.planet[i].starndx].mstar.value*n/(ss.planet[i].a.value*ss.constants.au*ss.constants.c^2*(1d0-ss.planet[i].e.value^2)) * 180d0/!dpi*36525*86400d0 ;; deg/century

   ;; derive cosi, b, and chord (depending on which is fit)
   if ss.planet[i].chord.fit then begin
      ss.planet[i].b.value = sqrt((1d0+ss.planet[i].p.value)^2-ss.planet[i].chord.value^2)
      ss.planet[i].cosi.value = ss.planet[i].b.value/$
                                (ss.planet[i].ar.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value))
   endif else if ss.planet[i].b.fit then begin
      ss.planet[i].cosi.value = ss.planet[i].b.value/$
                                (ss.planet[i].ar.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value))
;;      ss.planet[i].chord.value = 
   endif else if ss.planet[i].cosi.fit then begin
      ss.planet[i].b.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*$
                             (1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value) ;; eq 7, Winn 2010
;;      ss.planet[i].chord.value = 
   endif

   ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   ss.planet[i].ideg.value = ss.planet[i].i.value*180d0/!dpi
   sini = sin(ss.planet[i].i.value)

   ss.planet[i].msini.value = ss.planet[i].mp.value*sini                                         ;; m_jupiter
   ss.planet[i].msiniearth.value = ss.planet[i].mpearth.value*sini                               ;; m_earth

   ;; derive RV semi-amplitude
   ss.planet[i].k.value = (2d0*!dpi*G/(ss.planet[i].period.value*(ss.star[ss.planet[i].starndx].mstar.value + ss.planet[i].mpsun.value)^2d0))^(1d0/3d0) * $
                          ss.planet[i].mpsun.value*sin(ss.planet[i].i.value)/sqrt(1d0-ss.planet[i].e.value^2)*$
                          ss.constants.rsun/ss.constants.meter/ss.constants.day ;; m/s

   ;; time of periastron
   ss.planet[i].phase.value=exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/pri)  
   phase2 = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/secondary)
   phasea = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/ascending)
   phased = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/descending)

   if ss.planet[i].tt.fit then begin
      ss.planet[i].tc.value = tc2tt(ss.planet[i].tt.value,$
                                    ss.planet[i].e.value,$
                                    ss.planet[i].i.value,$
                                    ss.planet[i].omega.value,$
                                    ss.planet[i].period.value,/tt2tc)                                    
   endif else if ss.planet[i].tt.derive then begin
      ss.planet[i].tt.value = tc2tt(ss.planet[i].tc.value,$
                                    ss.planet[i].e.value,$
                                    ss.planet[i].i.value,$
                                    ss.planet[i].omega.value,$
                                    ss.planet[i].period.value)
   endif

   ss.planet[i].tp.value = ss.planet[i].tc.value - ss.planet[i].period.value*(ss.planet[i].phase.value)
   ss.planet[i].ts.value = ss.planet[i].tc.value - ss.planet[i].period.value*(ss.planet[i].phase.value-phase2)
   ss.planet[i].ta.value = ss.planet[i].tc.value - ss.planet[i].period.value*(ss.planet[i].phase.value-phasea)
   ss.planet[i].td.value = ss.planet[i].tc.value - ss.planet[i].period.value*(ss.planet[i].phase.value-phased)

   if ss.planet[i].te.derive then begin
      ss.planet[i].te.value = tc2tt(ss.planet[i].ts.value,$
                                    ss.planet[i].e.value,$
                                    ss.planet[i].i.value,$
                                    ss.planet[i].omega.value,$
                                    ss.planet[i].period.value,/ts2te)
   endif

   ;; if given a prior, select the epoch closest to that prior
   if ss.planet[i].ts.prior ne 0d0 then begin
      nper = round((ss.planet[i].ts.prior - ss.planet[i].ts.value)/ss.planet[i].period.value)
      ss.planet[i].ts.value += nper*ss.planet[i].period.value
   endif else begin
      ;; otherwise select the epoch closest to Tc
      nper = round((ss.planet[i].tc.value - ss.planet[i].ts.value)/ss.planet[i].period.value)
      ss.planet[i].ts.value += nper*ss.planet[i].period.value
   endelse

   ;; it's possible tp,ts,ta,td could be split down the middle 
   ;; then the median would be meaningless -- correct that
   medper = median(ss.planet[i].period.value)
   ss.planet[i].tp.value = exofast_recenter(ss.planet[i].tp.value, medper)
   ss.planet[i].ts.value = exofast_recenter(ss.planet[i].ts.value, medper)
   ss.planet[i].ta.value = exofast_recenter(ss.planet[i].ta.value, medper)
   ss.planet[i].td.value = exofast_recenter(ss.planet[i].td.value, medper)

   ss.planet[i].teq.value = ss.star[ss.planet[i].starndx].teff.value*sqrt(1d0/(2d0*ss.planet[i].ar.value)) ;(f*(1d0-Ab))^(0.25d0)
   ss.planet[i].dr.value = ss.planet[i].ar.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value) ;; d/rstar = star-planet separation at transit

   Qp = 1d6 ;; tidal Q value   
   ss.planet[i].tcirc.value = 4d0*Qp/63d0/(ss.constants.day*365.25d9)*$
                              ((ss.planet[i].a.value*ss.constants.au)^3/(ss.constants.GMsun*ss.star[ss.planet[i].starndx].mstar.value))^(1d0/2d0)*$
                              (ss.planet[i].mpsun.value/ss.star[ss.planet[i].starndx].mstar.value)*$
                              (ss.planet[i].ar.value/ss.planet[i].p.value)^5d0*$
                              (1d0-ss.planet[i].e.value^2)^(13d0/2d0)/$
                              (1d0+6d0*ss.planet[i].e.value^2) ;; Adams & Laughlin, 2006, eq 2 (in Gyr)
;   ss.planet[i].tcirc.value = 1.6d0*ss.planet[i].mp.value*ss.star[ss.planet[i].starndx].mstar.value^(-3d0/2d0)*ss.planet[i].rp.value^(-5d0)*(ss.planet[i].a.value/0.05d0)^(13d0/2d0) ;; Adams & Laughlin, 2006, eq 3

   ss.planet[i].fave.value = ss.constants.sigmab*ss.star[ss.planet[i].starndx].teff.value^4/(ss.planet[i].ar.value*(1d0+ss.planet[i].e.value^2/2d0))^2/1d9    ;; 10^9 erg/s/cm^2

   ss.planet[i].b.value  = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)  ;; eq 7, Winn 2010
   ss.planet[i].bs.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0-ss.planet[i].esinw.value)  ;; eq 8, Winn 2010

   ;; approximate durations taken from Winn 2010 (close enough; these should only be used to schedule observations anyway)
   ss.planet[i].t14.value = ss.planet[i].period.value/!dpi*asin(sqrt((1d0+abs(ss.planet[i].p.value))^2 - ss.planet[i].b.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)
   ;; eqs 14, 16, Winn 2010
   t23 = ss.planet[i].period.value/!dpi*asin(sqrt((1d0-abs(ss.planet[i].p.value))^2 - ss.planet[i].b.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)

   ;; no transit, transit duration equation is undefined -- set to zero
   notransit = where(abs(ss.planet[i].b.value) gt 1d0+abs(ss.planet[i].p.value))
   if notransit[0] ne -1 then ss.planet[i].t14.value[notransit] = 0d0

   ;; grazing transit, the flat part of transit is zero
   grazing = where(abs(ss.planet[i].b.value) gt 1d0-abs(ss.planet[i].p.value))
   if grazing[0] ne -1 then t23[grazing] = 0d0

   ss.planet[i].tau.value = (ss.planet[i].t14.value-t23)/2d0
   ss.planet[i].tfwhm.value = ss.planet[i].t14.value-ss.planet[i].tau.value

   ss.planet[i].ptg.value = (ss.star[ss.planet[i].starndx].rstar.value+ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 + ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2) ;; eq 9, Winn 2010
   ss.planet[i].pt.value = (ss.star[ss.planet[i].starndx].rstar.value-ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 + ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2)  ;; eq 9, Winn 2010

   ;; approximate durations taken from Winn 2010 (close enough; these should only be used to schedule observations anyway)
   ss.planet[i].t14s.value = ss.planet[i].period.value/!dpi*asin(sqrt((1d0+abs(ss.planet[i].p.value))^2 - ss.planet[i].bs.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0-ss.planet[i].esinw.value) ;; eqs 14, 16, Winn 2010
   t23s = ss.planet[i].period.value/!dpi*asin(sqrt((1d0-abs(ss.planet[i].p.value))^2 - ss.planet[i].bs.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0-ss.planet[i].esinw.value)   

   ;; no transit => transit duration equation is undefined -- set to zero
   notransit = where(abs(ss.planet[i].bs.value) gt 1d0+abs(ss.planet[i].p.value))
   if notransit[0] ne -1 then ss.planet[i].t14s.value[notransit] = 0d0

   ;; grazing transit, the flat part of transit is zero
   grazing = where(abs(ss.planet[i].bs.value) gt 1d0-abs(ss.planet[i].p.value))
   if grazing[0] ne -1 then t23s[grazing] = 0d0

   ss.planet[i].taus.value = (ss.planet[i].t14s.value-t23s)/2d0
   ss.planet[i].tfwhms.value = ss.planet[i].t14s.value-ss.planet[i].taus.value

   ;; check for NaNs in the distributions -- this is a problem
   if (ss.planet[i].tfwhms.derive and (where(~finite(ss.planet[i].tfwhms.value)))[0] ne -1) or $
      (ss.planet[i].taus.derive and (where(~finite(ss.planet[i].taus.value)))[0] ne -1) or $
      (ss.planet[i].tfwhm.derive and (where(~finite(ss.planet[i].tfwhm.value)))[0] ne -1) or $
      (ss.planet[i].tau.derive and (where(~finite(ss.planet[i].tau.value)))[0] ne -1) then begin
      printandlog, 'There are NaNs in the derived transit/eclipse duration distributions. This is a bug. Please contact jason.eastman@cfa.harvard.edu ***and do not exit IDL***. It is likely we can recover your run and fix the bug. Exiting IDL will cause you to lose all results and make it far more difficult to diagnose the problem.', ss.logname
      stop
   endif

   ss.planet[i].psg.value = (ss.star[ss.planet[i].starndx].rstar.value+ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 - ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2) ;; eq 9, Winn 2010
   ss.planet[i].ps.value = (ss.star[ss.planet[i].starndx].rstar.value-ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 - ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2)  ;; eq 9, Winn 2010

   ss.planet[i].rhop.value = ss.planet[i].mpsun.value/(ss.planet[i].rpsun.value^3)*ss.constants.rhosun ;; cgs
   ss.planet[i].loggp.value = alog10(ss.planet[i].mpsun.value/ss.planet[i].rpsun.value^2*ss.constants.gravitysun) ;; cgs
   ss.planet[i].safronov.value = ss.planet[i].ar.value*ss.planet[i].q.value/ss.planet[i].p.value ;; unitless

   ;; depth != delta if grazing (ignore limb darkening)
   ss.planet[i].delta.value = ss.planet[i].p.value^2
   for k=0L, ss.nband-1L do begin
      tagnames = tag_names(ss.planet)
      match = (where(tagnames eq strupcase('depth_' + ss.band[k].name)))[0]
      if match eq -1 then begin
         printandlog, 'ERROR: no match for ' + ss.band[k].name + ' ...  this should not be possible'
         continue
      endif
      for j=0L, ss.nsteps-1L do begin
         
         u1 = ss.band[k].u1.value[j]
         u2 = ss.band[k].u1.value[j]

         exofast_occultquad_cel, abs(ss.planet[i].b.value[j]), u1, u2, ss.planet[i].p.value[j],mu1
         ss.planet[i].(match).value[j] = 1d0-mu1[0]

;         exofast_occultquad_cel, abs(ss.planet[i].b.value[j]), 0d0, 0d0, ss.planet[i].p.value[j],mu1
;         ss.planet[i].depth.value[j] = 1d0-mu1[0]

      endfor
   endfor

   ;; find the ideal Tc that minimizes the Tc uncertainty and the
   ;; covariance between Tc and Period
   minepoch = floor((minbjd - median(ss.planet[i].tc.value))/median(ss.planet[i].period.value))-1
   maxepoch =  ceil((maxbjd - median(ss.planet[i].tc.value))/median(ss.planet[i].period.value))+1
   bestepoch = (minepoch+maxepoch)/2d0
   bestepoch2 = (minepoch+maxepoch)/2d0
   mincovar = !values.d_infinity
   mincovar2 = !values.d_infinity
   burnndx = ss.burnndx
   goodchains = *(ss.goodchains)
   for epoch=minepoch, maxepoch do begin

      ;; trim the burn in and remove bad chains before computing covariances
      tt     = (reform(ss.planet[i].tt.value+epoch*ss.planet[i].period.value,$
                       ss.nsteps/ss.nchains,ss.nchains))[burnndx:*,goodchains]
      period = (reform(ss.planet[i].period.value,$
                       ss.nsteps/ss.nchains,ss.nchains))[burnndx:*,goodchains]

      tt = reform(tt,n_elements(tt))
      period = reform(period,n_elements(period))
      
      corr = correlate(transpose([[tt],[period]]))
      if abs(corr[0,1]) lt mincovar then begin
         mincovar = abs(corr[0,1])
         bestepoch = epoch
      endif

      ;; do the same for the secondary eclipse time
      te = (reform(ss.planet[i].te.value+epoch*ss.planet[i].period.value,$
                   ss.nsteps/ss.nchains,ss.nchains))[burnndx:*,goodchains]
      te = reform(te,n_elements(te))
      corr = correlate(transpose([[te],[period]]))
      if abs(corr[0,1]) lt mincovar2 then begin
         mincovars = abs(corr[0,1])
         bestepoch2 = epoch
      endif

   endfor

   ;; apply the best epoch to compute the optimal times
   ss.planet[i].tt.value = ss.planet[i].tt.value + bestepoch*ss.planet[i].period.value
   ss.planet[i].te.value = ss.planet[i].te.value + bestepoch2*ss.planet[i].period.value
   
   ;; convert to SSB frame
   for j=0L, ss.nsteps-1 do begin
      model_times = [ss.planet[i].tt.value[j],ss.planet[i].te.value[j],$
                     ss.planet[i].tc.value[j],ss.planet[i].ts.value[j]]
      
      observed_times = target2bjd(model_times,$
                                  inclination=ss.planet[i].i.value[j],$
                                  a=ss.planet[i].a.value[j],$
                                  tp=ss.planet[i].tp.value[j],$
                                  period=ss.planet[i].period.value[j],$
                                  e=ss.planet[i].e.value[j],$
                                  omega=ss.planet[i].omega.value[j],$
                                  q=1d0/ss.planet[i].q.value[j])
      
      ss.planet[i].t0.value[j] = observed_times[0] 
      ss.planet[i].te0.value[j] = observed_times[1] 
      ss.planet[i].tco.value[j] = observed_times[2] 
      ss.planet[i].tso.value[j] = observed_times[3] 

   endfor

   ;; blackbody eclipse depths

;   planetbb36 = exofast_blackbody(ss.planet[i].teq.value,replicate(3550d0/1d9,ss.nsteps),/wave)
;   x = ss.planet[i].p.value^2*planetbb36/starbb36
;   ss.planet[i].eclipsedepth36.value = x/(1d0+x)*1d6

;   planetbb45 = exofast_blackbody(ss.planet[i].teq.value,replicate(4493d0/1d9,ss.nsteps),/wave)
;   x = ss.planet[i].p.value^2*planetbb45/starbb45
;   ss.planet[i].eclipsedepth45.value = x/(1d0+x)*1d6

   planetbb25 = exofast_blackbody(ss.planet[i].teq.value,replicate(2500d0/1d9,ss.nsteps),/wave)
   x = ss.planet[i].p.value^2*planetbb25/starbb25
   ss.planet[i].eclipsedepth25.value = x/(1d0+x)*1d6

   planetbb50 = exofast_blackbody(ss.planet[i].teq.value,replicate(5000d0/1d9,ss.nsteps),/wave)
   x = ss.planet[i].p.value^2*planetbb50/starbb50
   ss.planet[i].eclipsedepth50.value = x/(1d0+x)*1d6

   planetbb75 = exofast_blackbody(ss.planet[i].teq.value,replicate(7500d0/1d9,ss.nsteps),/wave)
   x = ss.planet[i].p.value^2*planetbb75/starbb75
   ss.planet[i].eclipsedepth75.value = x/(1d0+x)*1d6

endfor

for i=0L, ss.nband-1 do begin  
   massfraction = ss.planet[0].mpsun.value/(ss.star[ss.planet[0].starndx].mstar.value + ss.planet[0].mpsun.value)
;   fluxfraction = ss.band[i].dilute.value
;   ss.band[i].phottobary.value = 1d0/(massfraction-fluxfraction)
   ss.band[i].eclipsedepth.value = ss.band[i].thermal.value + ss.band[i].reflect.value
endfor

end
