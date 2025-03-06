pro plotrv, ss, psname=psname, ndx=ndx, range=range, savfile=savfile

if n_elements(savfile) ne 0 then begin
   restore, savfile
   ss = mcmcss
endif

;; pick the best-fit index if not specified
if n_elements(ndx) eq 0 then begin
   if ss.nsteps eq 1 then ndx = 0 $
   else minchi2 = min(*ss.chi2,ndx)
endif

defsysv, '!GDL', exists=runninggdl
mydevice=!d.name
if keyword_set(psname) then begin
   set_plot, 'PS'
   aspect_ratio=1.5
   xsize=10.5
   ysize=xsize/aspect_ratio
   !p.font=0

   if runninggdl then psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.1.ps' $
   else psname0 = psname

   device, filename=psname0, /color, bits=24
   device, xsize=xsize,ysize=ysize
   loadct, 39, /silent
   colors = [0,159,95,254,223,31,207,111,191,47]
   red = 254
   symsize=0.5
   title = strarr(ss.nplanets)
   position1 = [0.23, 0.40, 0.95, 0.95]    ;; data plot
   position2 = [0.23, 0.20, 0.95, 0.40]    ;; residual plot
endif else begin
;   set_plot, 'X'
   colors= ['ffffff'x,'0000ff'x,'00ff00'x,'ff0000'x,'0080ff'x,$
            '800080'x,'00ffff'x,'ffff00'x,'80d000'x,'660000'x]
   red = '0000ff'x
   symsize = 1d0
   charsize = 1
   device,window_state=win_state
   symsize=1

   title = ss.planet.label
   position1 = [0.07, 0.22, 0.97, 0.95]    ;; data plot
   position2 = [0.07, 0.07, 0.97, 0.22]    ;; residual plot
   font=-1

;   if win_state[20] eq 1 then wset, 20 $
;   else window, 20, retain=2
endelse

symbols = [0,8,4,3,0,8,4,3]
fills = [1,1,1,1,0,0,0,0]
nsymbols = n_elements(symbols)
nfills = n_elements(fills)
ncolors = n_elements(colors)

allmindate = !values.d_infinity
allmaxdate = -!values.d_infinity

q = dblarr(ss.nplanets)
for i=0L, ss.nplanets-1 do begin
   q[i] = ss.star[ss.planet[i].starndx].mstar.value[ndx]/ss.planet[i].mpsun.value[ndx]
endfor

starrvs = [-1]
companionrvs = [-1]

legendndx = lonarr(ss.nplanets+1, ss.ntel)
rmtime = 0d0
for j=0L, ss.ntel-1 do begin
   rv = *(ss.telescope[j].rvptrs)
   mindate = min(rv.bjd,max=maxdate)
   if maxdate - mindate lt 0.5 then rmtime = rv.bjd
   if mindate lt allmindate then allmindate = mindate
   if maxdate gt allmaxdate then allmaxdate = maxdate
   
   if rv.planet eq -1 then legendndx[ss.nplanets,j] = 1B $
   else legendndx[rv.planet,j] = 1B

endfor
;t0 = (allmindate+allmaxdate)/2d0
t0 = ss.rvepoch

roundto = 10L^(strlen(strtrim(floor(allmaxdate-allmindate),2)));+1L) 
bjd0 = floor(allmindate/roundto)*roundto

allmindate = (allmindate-bjd0)*0.95 + bjd0
allmaxdate = (allmaxdate-bjd0)*1.05 + bjd0

;; if the timespan is > 100 years, there's got to be an error in the timestamps
if (allmaxdate-allmindate) gt 36500 then begin
   message, 'WARNING: RV data span > 100 years. Make sure your input timestamps are consistent. This may cause memory issues', /continue
endif

;; sample the pretty light curve 100 times per (minimum) period 
;; for the span of the data
cadence = min(ss.planet.period.value[ndx])/100d0
nsteps = (allmaxdate-allmindate)/cadence
prettytime = allmindate + (allmaxdate-allmindate)*dindgen(nsteps)/(nsteps-1.d0)

if n_elements(rmtime) gt 2 then begin
   rmnsteps = (max(rmtime)+0.5-min(rmtime)-0.5)/(1.0/60.0/24.0)
   rmprettytime = min(rmtime)-0.5 + (max(rmtime)-min(rmtime)+1)*dindgen(rmnsteps)/(rmnsteps-1.d0)
   prettytime = [prettytime, rmprettytime]
   prettytime = prettytime[sort(prettytime)]
endif




allprettymodel = prettytime*0d0
allprettymodel2 = prettytime*0d0

if not keyword_set(psname) then begin
   if win_state[20] eq 1 then wset, 20 $
   else window, 20, retain=2
   nrvfit = n_elements(where(ss.planet.fitrv))
   ny = ceil(sqrt(nrvfit))
   nx = ceil(nrvfit/double(ny))
   !p.multi = [0,nx,ny]
endif

;; calculate the residuals
for j=0, ss.ntel-1 do begin
   rv = *(ss.telescope[j].rvptrs)
   mintime = min(rv.bjd,max=maxtime)

   ;; subtract gamma, slope, and quadratic terms
;;   rv.residuals = rv.rv - (ss.telescope[j].gamma.value[ndx] + ss.star.slope.value[ndx]*(rv.bjd-t0) + ss.star.quad.value[ndx]*(rv.bjd-t0)^2)
   ; modelrv = (ss.telescope[j].gamma.value[ndx] + ss.star[0].slope.value[ndx]*(rv.bjd-t0) + ss.star[0].quad.value[ndx]*(rv.bjd-t0)^2)
   ; 
   rvtime = ((*(ss.telescope[j].rvptrs)).bjd)
   mintime = min(rvtime,max=maxtime)
   t0_each_rv = (mintime+maxtime)/2.d0
   modelrv = (ss.telescope[j].gamma.value + ss.telescope[j].srv.value*(rv.bjd-t0_each_rv) + ss.telescope[j].qrv.value*(rv.bjd-t0_each_rv)^2)
   modelrv = modelrv + ss.star[0].slope.value[ndx]*(rv.bjd-t0) + ss.star[0].quad.value[ndx]*(rv.bjd-t0)^2
;   exofast_forprint, rv.bjd, modelrv, textout=base+'.rv.trend.txt', format='(f0.10,x,f0.10)'
;   stop

   for i=0, ss.nplanets-1 do begin
      
      if ~ss.planet[i].fitrv then continue ;; skip unfit planets
      if rv.planet ne -1 then continue     ;; skip direct planet RVs

      ;; this clause is never executed
      if rv.planet eq i then begin
         print, 'rv.planet eq i', rv.planet, i
;; this needs to be debugged
         rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                            a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                            period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                            omega=ss.planet[i].omega.value[ndx]+!dpi,$
                            c=ss.constants.c/ss.constants.au*ss.constants.day,$
                            q=q[i])
         
         ;; calculate the RV model
         modelrv += exofast_rv(rvbjd,ss.planet[i].tp.value[ndx],ss.planet[i].period.value[ndx],$
                               0d0,ss.planet[i].K.value[ndx]*q[i],$
                               ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx]+!dpi,$
                               slope=0d0)
         
      endif else begin
         rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                            a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                            period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                            omega=ss.planet[i].omega.value[ndx]+!dpi,$
                            c=ss.constants.c/ss.constants.au*ss.constants.day,$
                            q=q[i], /primary)

         if ss.planet[i].rossiter then begin
            coeffs = quadld(ss.star[ss.planet[i].starndx].logg.value,$
                        ss.star[ss.planet[i].starndx].teff.value,$
                        ss.star[ss.planet[i].starndx].feh.value,'V')
            u1 = coeffs[0]
            u2 = coeffs[1]
         endif else begin
            u1 = 0d0 
            u2 = 0d0
         endelse
         band = ss.band[ss.telescope[j].bandndx]
         u1 = band.u1.value
         u2 = band.u2.value
         if ss.planet[i].svsinicoslambda.value eq 0d0 then begin
            this_lambda = 0d0
         endif else begin
            this_lambda = atan(ss.planet[i].svsinisinlambda.value, ss.planet[i].svsinicoslambda.value)
         endelse
         this_vsini = ss.planet[i].svsinicoslambda.value^2 + ss.planet[i].svsinisinlambda.value^2
         modelrv += exofast_rv(rvbjd,ss.planet[i].tp.value[ndx]+ss.telescope[j].rmttv.value,$
                               ss.planet[i].period.value[ndx],0d0,$
                               ss.planet[i].K.value[ndx],ss.planet[i].e.value[ndx],$
                               ss.planet[i].omega.value[ndx],slope=0,$
                               rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                              p=abs(ss.planet[i].p.value),vsini=this_vsini,$
                              lambda=this_lambda,$
                              vgamma=ss.star[ss.planet[i].starndx].vgamma.value, vzeta=ss.star[ss.planet[i].starndx].vzeta.value,$
                              vxi=ss.star[ss.planet[i].starndx].vxi.value, valpha=ss.star[ss.planet[i].starndx].valpha.value,$
                              u1=u1,u2=u2,deltarv=deltarv, exptime=ss.telescope[j].exptime, ninterp=ss.telescope[j].ninterp,$
                              srv=0, qrv=0)
      endelse
      
   endfor
   if tag_exist(ss.telescope[j],'detrend') then begin
      detrendadd = total((*ss.telescope[j].rvptrs).detrendadd*(replicate(1d0,n_elements((*ss.telescope[j].rvptrs).bjd))##(*ss.telescope[j].rvptrs).detrendaddpars.value),1)
      detrendmult = (1d0+total((*ss.telescope[j].rvptrs).detrendmult*(replicate(1d0,n_elements((*ss.telescope[j].rvptrs).bjd))##(*ss.telescope[j].rvptrs).detrendmultpars.value),1))
   endif else begin
      detrendadd = 0d0
      detrendmult = 1d0
   endelse
   
   modelrv += total((*ss.telescope[j].rvptrs).detrendadd*(replicate(1d0,n_elements((*ss.telescope[j].rvptrs).bjd))##(*ss.telescope[j].rvptrs).detrendaddpars.value),1)
   modelrv *= (1d0+total((*ss.telescope[j].rvptrs).detrendmult*(replicate(1d0,n_elements((*ss.telescope[j].rvptrs).bjd))##(*ss.telescope[j].rvptrs).detrendmultpars.value),1))
   ;; re-populate the residual array
   rv.residuals = rv.rv - modelrv
   *(ss.telescope[j].rvptrs) = rv
endfor

for i=0, ss.nplanets-1 do begin
   
   if ~ss.planet[i].fitrv then continue
   if rv.planet ne -1 then continue

   if rv.planet eq i then begin
;; this needs to be debugged

      ;; pretty model without quad, slope, or gamma
      prettymodel2 = exofast_rv(prettytime,ss.planet[i].tp.value[ndx],$
                                ss.planet[i].period.value[ndx],0d0,ss.planet[i].K.value[ndx]*q[i],$
                                ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx]+!dpi)
      allprettymodel2 += prettymodel2

   endif else begin


      if ss.planet[i].rossiter then begin
         coeffs = quadld(ss.star[ss.planet[i].starndx].logg.value,$
                     ss.star[ss.planet[i].starndx].teff.value,$
                     ss.star[ss.planet[i].starndx].feh.value,'V')
         u1 = coeffs[0]
         u2 = coeffs[1]
      endif else begin
         u1 = 0d0 
         u2 = 0d0
      endelse

      if ss.planet[i].svsinicoslambda.value eq 0d0 then begin
         this_lambda = 0d0
      endif else begin
         this_lambda = atan(ss.planet[i].svsinisinlambda.value, ss.planet[i].svsinicoslambda.value)
      endelse
      this_vsini = ss.planet[i].svsinicoslambda.value^2 + ss.planet[i].svsinisinlambda.value^2
      
      ;; pretty model without quad, slope, or gamma
      prettymodel = exofast_rv(prettytime,ss.planet[i].tp.value[ndx],$
                               ss.planet[i].period.value[ndx],0d0,ss.planet[i].K.value[ndx],$
                               ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx],$
                               slope=0,$
                               rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                               p=abs(ss.planet[i].p.value[ndx]),vsini=this_vsini,$
                               lambda=this_lambda,$
                               vgamma=ss.star[ss.planet[i].starndx].vgamma.value, vzeta=ss.star[ss.planet[i].starndx].vzeta.value,$
                               vxi =ss.star[ss.planet[i].starndx].vxi.value, valpha=ss.star[ss.planet[i].starndx].valpha.value,$
                               u1=u1,u2=u2,deltarv=prettydeltarv,srv=0,qrv=0)
                              ;  p=abs(ss.planet[i].p.value[ndx]),vsini=this_vsini[ndx],$
                              ;  lambda=this_lambda,vbeta=ss.star[ss.planet[i].starndx].vbeta.value,$
                              ;  vgamma=ss.star[ss.planet[i].starndx].vgamma.value[ndx], vzeta=ss.star[ss.planet[i].starndx].vzeta.value[ndx],$
                              ;  u1=u1,u2=u2,deltarv=prettydeltarv, exptime=ss.telescope[j].exptime[ndx], ninterp=ss.telescope[j].ninterp[ndx])
      allprettymodel += prettymodel
      ; print,this_lambda,this_vsini,ss.star[ss.planet[i].starndx].vgamma.value, ss.star[ss.planet[i].starndx].vzeta.value, ss.star[ss.planet[i].starndx].vxi.value, ss.star[ss.planet[i].starndx].valpha.value
      ; print,u1,u2
      ;print,this_lambda,
      if n_elements(psname) eq 1 then begin
         base = file_dirname(psname) + path_sep() + 'modelfiles' + path_sep() + file_basename(psname,'.rv.ps')
         exofast_forprint, prettytime, prettymodel, textout=base+'.prettymodelrv.planet.' + string(i,format='(i02)') + '.txt', format='(f0.10,x,f0.10)'
      endif
   endelse 

   ;; phase so primary is at 0.25
   prettyphase = (((prettytime - ss.planet[i].tc.value[ndx]) mod $
                   ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                  1.25d0) mod 1
   sorted = sort(prettyphase)

   allminrv = min(allprettymodel,max=allmaxrv)
   for j=0, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet ne -1 then continue

      err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])

      if ss.planet[i].svsinicoslambda.value eq 0d0 then begin
         this_lambda = 0d0
      endif else begin
         this_lambda = atan(ss.planet[i].svsinisinlambda.value, ss.planet[i].svsinicoslambda.value)
      endelse
      this_vsini = ss.planet[i].svsinicoslambda.value^2 + ss.planet[i].svsinisinlambda.value^2
      band = ss.band[ss.telescope[j].bandndx]
      u1 = band.u1.value
      u2 = band.u2.value
      modelrv = exofast_rv(rv.bjd,ss.planet[i].tp.value[ndx]+ss.telescope[j].rmttv.value,$
                           ss.planet[i].period.value[ndx],0d0,$
                           ss.planet[i].K.value[ndx],ss.planet[i].e.value[ndx],$
                           ss.planet[i].omega.value[ndx],slope=0,$
                           rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                           p=abs(ss.planet[i].p.value),vsini=this_vsini,$
                           lambda=this_lambda,$
                           vgamma=ss.star[ss.planet[i].starndx].vgamma.value, vzeta=ss.star[ss.planet[i].starndx].vzeta.value,$
                           vxi=ss.star[ss.planet[i].starndx].vxi.value, valpha=ss.star[ss.planet[i].starndx].valpha.value,$
                           u1=u1,u2=u2,deltarv=deltarv, exptime=ss.telescope[j].exptime, ninterp=ss.telescope[j].ninterp)
      mintime = min(rv.bjd,max=maxtime)
      minrv = min(rv.residuals-err+modelrv)
      maxrv = max(rv.residuals+err+modelrv)
      
      if minrv lt allminrv then allminrv = minrv
      if maxrv gt allmaxrv then allmaxrv = maxrv
   endfor
   xmin = 0d0
   xmax = 1d0

   if finite(range[0]) then xmin = range[0]
   if finite(range[1]) then xmax = range[1]
   if finite(range[2]) then allminrv = range[2]
   if finite(range[3]) then allmaxrv = range[3]
   xrange = [xmin,xmax]

   xtitle1='!3' + exofast_textoidl('Phase + (T_P - T_C)/P + 0.25',font=font)
   plot, [0], [0], xrange=xrange, yrange=[allminrv,allmaxrv], $
         ytitle='!3RV (m/s)', charsize=charsize, title=title[i],position=position1,xtickformat='(A1)'
   for j=0, ss.ntel-1 do begin 
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet ne -1 then continue
      err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])

      if ss.planet[i].svsinicoslambda.value eq 0d0 then begin
         this_lambda = 0d0
      endif else begin
         this_lambda = atan(ss.planet[i].svsinisinlambda.value, ss.planet[i].svsinicoslambda.value)
      endelse
      this_vsini = ss.planet[i].svsinicoslambda.value^2 + ss.planet[i].svsinisinlambda.value^2
      band = ss.band[ss.telescope[j].bandndx]
      u1 = band.u1.value
      u2 = band.u2.value
      modelrv = exofast_rv(rv.bjd,ss.planet[i].tp.value[ndx]+ss.telescope[j].rmttv.value,$
                           ss.planet[i].period.value[ndx],0d0,$
                           ss.planet[i].K.value[ndx],ss.planet[i].e.value[ndx],$
                           ss.planet[i].omega.value[ndx],slope=0,$
                           rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                           p=abs(ss.planet[i].p.value),vsini=this_vsini,$
                           lambda=this_lambda,$
                           vgamma=ss.star[ss.planet[i].starndx].vgamma.value, vzeta=ss.star[ss.planet[i].starndx].vzeta.value,$
                           vxi=ss.star[ss.planet[i].starndx].vxi.value, valpha=ss.star[ss.planet[i].starndx].valpha.value,$
                           u1=u1,u2=u2,deltarv=deltarv, exptime=ss.telescope[j].exptime, ninterp=ss.telescope[j].ninterp,t0=0d0,srv=0, qrv=0)
      time=(((rv.bjd-ss.planet[i].tc.value[ndx]) mod ss.planet[i].period.value[ndx])/$
            ss.planet[i].period.value[ndx]+1.25d0) mod 1
      plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
      oploterr, time, rv.residuals+modelrv, err, 8
      if keyword_set(psname) then begin
         base = file_dirname(psname) + path_sep() + 'modelfiles' + path_sep() + file_basename(psname,'.model')
         exofast_forprint, time, rv.residuals,modelrv, err, format='(f0.8,x,f0.6,x,f0.6,x,f0.6)', textout=base + '.rvphase.residuals.planet_'+ string(i,format='(i02)')+'.telescope_' + string(j,format='(i02)') + '.txt', /nocomment,/silent
      endif
   endfor

   if keyword_set(psname) then begin
      base = file_dirname(psname) + path_sep() + 'modelfiles' + path_sep() + file_basename(psname,'.model')
      exofast_forprint, prettyphase[sorted], prettymodel[sorted], format='(f0.8,x,f0.6)', textout=base + '.rvphase.model.planet_'+ string(i,format='(i02)') + '.txt', /nocomment,/silent
   endif
   
   oplot, prettyphase[sorted], prettymodel[sorted], color=red
   use = where(legendndx[ss.nplanets,*],nuse)
   if nuse gt 1 then exofast_legend, ss.telescope[use].label, color=colors[(indgen(ss.ntel) mod ncolors)[use]],/bottom,/right,psym=symbols[(indgen(ss.ntel) mod nsymbols)[use]], /useplotsym, charsize=0.5, fill=fills[(indgen(ss.ntel) mod nfills)[use]]

   ;; plot the residuals below
   ymin = !values.d_infinity
   ymax = -!values.d_infinity
   for j=0, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet ne -1 then continue
      err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
      ymin = min([(rv.residuals - err)*1.1,ymin])
      ymax = max([(rv.residuals + err)*1.1,ymax])
   endfor

   ;; make the plot symmetric about 0
   if ymin lt -ymax then ymax = -ymin
   if ymax gt -ymin then ymin = -ymax
   if finite(range[4]) then ymin = range[4]
   if finite(range[5]) then ymax = range[5]

   plot, [0],[0], position=position2, /noerase, $
         xrange=xrange, xtitle=xtitle1,$
         yrange=[ymin,ymax]/0.7d0, ytitle='O-C (m/s)', $
         /xstyle, /ystyle, yminor=2,yticks=2, ytickv=[ymin,0,ymax]
   oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  
   
   for j=0L, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet ne -1 then continue
      err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
      time=(((rv.bjd-ss.planet[i].tc.value[ndx]) mod ss.planet[i].period.value[ndx])/$
            ss.planet[i].period.value[ndx]+1.25d0) mod 1
      plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
      oploterr, time, rv.residuals, err, 8
   endfor

   ;; GDL can't do multi-page plots
   if runninggdl then begin
      device, /close
      psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.' + strtrim(i+2,2) + '.ps'
      device, filename=psname0, /color, bits=24, xsize=xsize,ysize=ysize
   endif

   ;; plot the RM (if modeled)
   if ~ss.planet[i].rossiter then continue

   period = ss.planet[i].period.value[ndx]

   ;; compute the transit duration
   sini = sin(ss.planet[i].i.value[ndx])
   t14 = period/!dpi*asin(sqrt((1d0+abs(ss.planet[i].p.value[ndx]))^2 - ss.planet[i].b.value[ndx]^2)/(sini*ss.planet[i].ar.value[ndx]))*sqrt(1d0-ss.planet[i].e.value[ndx]^2)/(1d0+ss.planet[i].esinw.value[ndx])
   xrange = [-t14,t14]*24d0 ;; hours
   ;; plot the residuals below

   ;; get the y limits of the plot main plot
   ymin = !values.d_infinity
   ymax = -!values.d_infinity
   for j=0, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet eq -1 then begin
         time = exofast_mod(rv.bjd - ss.planet[i].tc.value[ndx],period,/negative)*24d0
         inrange = where(time ge xrange[0] and time le xrange[1])
         if inrange[0] ne -1 then begin
            if ss.planet[i].svsinicoslambda.value eq 0d0 then begin
               this_lambda = 0d0
            endif else begin
               this_lambda = atan(ss.planet[i].svsinisinlambda.value, ss.planet[i].svsinicoslambda.value)
            endelse
            this_vsini = ss.planet[i].svsinicoslambda.value^2 + ss.planet[i].svsinisinlambda.value^2
            if ss.telescope[j].bandndx eq -1 then begin
               u1 = 0d0
            endif else begin
               band = ss.band[ss.telescope[j].bandndx]
               u1 = band.u1.value
               u2 = band.u2.value
            endelse
            modelrv = exofast_rv(rv.bjd[inrange],ss.planet[i].tp.value[ndx]+ss.telescope[j].rmttv.value,$
                                 ss.planet[i].period.value[ndx],0d0,$
                                 ss.planet[i].K.value[ndx],ss.planet[i].e.value[ndx],$
                                 ss.planet[i].omega.value[ndx],slope=0,$
                                 rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                                 p=abs(ss.planet[i].p.value),vsini=this_vsini,$
                                 lambda=this_lambda,$
                                 vgamma=ss.star[ss.planet[i].starndx].vgamma.value, vzeta=ss.star[ss.planet[i].starndx].vzeta.value,$
                                 vxi=ss.star[ss.planet[i].starndx].vxi.value, valpha=ss.star[ss.planet[i].starndx].valpha.value,$
                                 u1=u1,u2=u2,deltarv=deltarv, exptime=ss.telescope[j].exptime, ninterp=ss.telescope[j].ninterp,t0=0d0, srv=0, qrv=0)

            err = sqrt(rv.err[inrange]^2 + ss.telescope[j].jittervar.value[ndx])
            ymin = min([deltarv + rv.residuals[inrange] - err,ymin])
            ymax = max([deltarv + rv.residuals[inrange] + err,ymax])
         endif
      endif
   endfor


   xtitle1='!3' + exofast_textoidl('Time - T_C (hrs)',font=font)
   plot, [0], [0], xrange=xrange, yrange=[ymin,ymax], $
         ytitle='!3RV (m/s)', charsize=charsize, title=title[i],position=position1,xtickformat='(A1)',/xstyle
   for j=0, ss.ntel-1 do begin 
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet eq -1 then begin        
         
         time = exofast_mod(rv.bjd - ss.planet[i].tc.value[ndx], period,/negative)*24d0


         if ss.planet[i].svsinicoslambda.value eq 0d0 then begin
            this_lambda = 0d0
         endif else begin
            this_lambda = atan(ss.planet[i].svsinisinlambda.value, ss.planet[i].svsinicoslambda.value)
         endelse
         this_vsini = ss.planet[i].svsinicoslambda.value^2 + ss.planet[i].svsinisinlambda.value^2
         band = ss.band[ss.telescope[j].bandndx]
         u1 = band.u1.value
         u2 = band.u2.value
         modelrv = exofast_rv(rv.bjd,ss.planet[i].tp.value[ndx]+ss.telescope[j].rmttv.value,$
                              ss.planet[i].period.value[ndx],0d0,$
                              ss.planet[i].K.value[ndx],ss.planet[i].e.value[ndx],$
                              ss.planet[i].omega.value[ndx],slope=0,$
                              rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                              p=abs(ss.planet[i].p.value),vsini=this_vsini,$
                              lambda=this_lambda,$
                              vgamma=ss.star[ss.planet[i].starndx].vgamma.value, vzeta=ss.star[ss.planet[i].starndx].vzeta.value,$
                              vxi=ss.star[ss.planet[i].starndx].vxi.value, valpha=ss.star[ss.planet[i].starndx].valpha.value,$
                              u1=u1,u2=u2,deltarv=deltarv, exptime=ss.telescope[j].exptime, ninterp=ss.telescope[j].ninterp,t0=0d0, srv=0, qrv=0)

         plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
         oploterr, time, (deltarv + rv.residuals), rv.err, 8
         if keyword_set(psname) then begin
            base = file_dirname(psname) + path_sep() + 'modelfiles' + path_sep() + file_basename(psname,'.model')
            exofast_forprint, time, rv.residuals, deltarv, rv.err, format='(f0.8,x,f0.6,x,f0.6,x,f0.6)', textout=base + '.rm.residuals.planet_'+ string(i,format='(i02)')+'.telescope_' + string(j,format='(i02)') + '.txt', /nocomment,/silent
         endif
      endif
   endfor
   prettyrmtime = exofast_mod(prettytime - ss.planet[i].tc.value[ndx],period,/negative)*24d0
   sorted = sort(prettyrmtime)
   oplot, prettyrmtime[sorted], prettydeltarv[sorted], color=red

   if keyword_set(psname) then begin
      base = file_dirname(psname) + path_sep() + 'modelfiles' + path_sep() + file_basename(psname,'.model')
      exofast_forprint, prettyrmtime[sorted], prettydeltarv[sorted], format='(f0.8,x,f0.6)', textout=base + '.rm.model.planet_'+ string(i,format='(i02)') + '.txt', /nocomment,/silent
   endif
   
   use = where(legendndx[ss.nplanets,*],nuse)
   if nuse gt 1 then exofast_legend, ss.telescope[use].label, color=colors[(indgen(ss.ntel) mod ncolors)[use]],/bottom,/right,psym=symbols[(indgen(ss.ntel) mod nsymbols)[use]], /useplotsym, charsize=0.5, fill=fills[(indgen(ss.ntel) mod nfills)[use]]

   ;; plot the residuals below
   ymin = !values.d_infinity
   ymax = -!values.d_infinity
   for j=0, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)

      if rv.planet eq -1 then begin

         time = exofast_mod(rv.bjd - ss.planet[i].tc.value[ndx],period,/negative)*24d0
         inrange = where(time ge xrange[0] and time le xrange[1])
         
         if inrange[0] ne -1 then begin
            err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
            ymin = min([(rv.residuals[inrange] - err[inrange])*1.1,ymin])
            ymax = max([(rv.residuals[inrange] + err[inrange])*1.1,ymax])
         endif

      endif
   endfor

   ;; make the plot symmetric about 0
   if ymin lt -ymax then ymax = -ymin
   if ymax gt -ymin then ymin = -ymax
   plot, [0],[0], position=position2, /noerase, $
         xrange=xrange, xtitle=xtitle1,$
         yrange=[ymin,ymax]/0.7, ytitle='O-C (m/s)', $
         /xstyle, /ystyle, yminor=2,yticks=2, ytickv=[ymin,0,ymax]
   oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  
   
   for j=0L, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet ne -1 then continue
      err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])

      time = exofast_mod(rv.bjd - ss.planet[i].tc.value[ndx],period,/negative)*24d0
      plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
      oploterr, time, rv.residuals, err, 8
   endfor

   ;; GDL can't do multi-page plots
   if runninggdl then begin
      device, /close
      psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.RM.' + strtrim(i+2,2) + '.ps'
      device, filename=psname0, /color, bits=24, xsize=xsize,ysize=ysize
   endif

endfor

!p.multi=0

;; now plot all planets, unphased, including the slope and quadratic
;; terms (but removing detrending)
if not keyword_set(psname) then begin
   if win_state[21] eq 1 then wset, 21 $
   else window, 21, retain=2
endif else begin
   trend = (prettytime-t0)*ss.star[0].slope.value[ndx] + (prettytime-t0)^2*ss.star[0].quad.value[ndx]
   allprettymodel += trend
   exofast_forprint, prettytime, prettymodel, textout=base+'.prettymodelrv.trend.txt', format='(f0.10,x,f0.10)'
endelse

if n_elements(yrange) eq 0 then begin
;; scale for the unphased plot
   allminrv = min(allprettymodel,max=allmaxrv)
   for j=0, ss.ntel-1 do begin
      rvstr = *(ss.telescope[j].rvptrs)
      if rvstr.planet ne -1 then continue
      rv = rvstr.rv - ss.telescope[j].gamma.value[ndx]
      err = sqrt(rvstr.err^2 + ss.telescope[j].jittervar.value[ndx])
      if min(rv-err) lt allminrv then allminrv = min(rv-err)
      if max(rv+err) gt allmaxrv then allmaxrv = max(rv+err)
   endfor
endif else begin
   allminrv = yrange[0]
   allmaxrv = yrange[1]
endelse

xtitle2='!3' + exofast_textoidl('BJD_{TDB} - ' + string(bjd0,format='(i7)'),font=font)
plot, [0], [0], xrange=[allmindate,allmaxdate]-bjd0,/xstyle,$
      yrange=[allminrv,allmaxrv], $
      ytitle='!3RV (m/s)', position=position1,xtickformat='(A1)'
oplot, prettytime-bjd0, allprettymodel, color=red

for j=0, ss.ntel-1 do begin 
   rv = *(ss.telescope[j].rvptrs)
   if rv.planet ne -1 then continue
   err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
   plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]

   if tag_exist(ss.telescope[j],'detrend') then begin
      detrendadd = total((*ss.telescope[j].rvptrs).detrendadd*(replicate(1d0,n_elements((*ss.telescope[j].rvptrs).bjd))##(*ss.telescope[j].rvptrs).detrendaddpars.value),1)
      detrendmult = (1d0+total((*ss.telescope[j].rvptrs).detrendmult*(replicate(1d0,n_elements((*ss.telescope[j].rvptrs).bjd))##(*ss.telescope[j].rvptrs).detrendmultpars.value),1))
   endif else begin
      detrendadd = 0d0
      detrendmult = 1d0
   endelse
   oploterr, rv.bjd-bjd0, (rv.rv-ss.telescope[j].gamma.value[ndx]-detrendadd)/detrendmult, err, 8
endfor

use = where(legendndx[ss.nplanets,*],nuse)
if nuse gt 1 then exofast_legend, ss.telescope[use].label, color=colors[(indgen(ss.ntel) mod ncolors)[use]],/top,/right,psym=symbols[(indgen(ss.ntel) mod nsymbols)[use]], /useplotsym, charsize=0.5, fill=fills[(indgen(ss.ntel) mod nfills)[use]]

;; plot the residuals below
plot, [0],[0], position=position2, /noerase, $
      xrange=[allmindate,allmaxdate]-bjd0, xtitle=xtitle2,$
      yrange=[ymin,ymax]/0.7, ytitle='O-C (m/s)', $
      /xstyle, /ystyle, yminor=2,yticks=2, ytickv=[ymin,0,ymax]
oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  

for j=0L, ss.ntel-1 do begin
   rv = *(ss.telescope[j].rvptrs)
   if rv.planet ne -1 then continue
   err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
   plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
   oploterr, rv.bjd-bjd0, rv.residuals, err, 8
endfor

;; plot the companion RVs
nplanetrvs = 0L
for j=0L, ss.ntel-1 do begin
   rv = *(ss.telescope[j].rvptrs)
   if rv.planet ne -1 then nplanetrvs++
endfor

if nplanetrvs gt 0L then begin

   for i=0L, ss.nplanets-1 do begin
      
      prettymodel = exofast_rv(prettytime,ss.planet[i].tp.value[ndx],$
                               ss.planet[i].period.value[ndx],0d0,ss.planet[i].K.value[ndx]*q[i],$
                               ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx]+!dpi)
      
      ;; phase so primary is at 0.25
      prettyphase = (((prettytime - ss.planet[i].tc.value[ndx]) mod $
                      ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                     1.25d0) mod 1
      sorted = sort(prettyphase)
      
      
      allminphase = 1d0
      allmaxphase = 0d0
      allminrv = !values.d_infinity
      allmaxrv = -!values.d_infinity
      allminoc = !values.d_infinity
      allmaxoc = -!values.d_infinity
      
      for j=0L, ss.ntel-1 do begin
         rv = *(ss.telescope[j].rvptrs)
         if rv.planet ne i then continue
         
         ;; this needs to be debugged
         rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                            a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                            period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                            omega=ss.planet[i].omega.value[ndx]+!dpi,$
                            c=ss.constants.c/ss.constants.au*ss.constants.day,$
                            q=q[i])
         
         phasetime = ((( rvbjd - ss.planet[i].tc.value[ndx]) mod $
                       ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                      1.25d0) mod 1
         
         minphase = min(phasetime,max=maxphase)
         if minphase lt allminphase then allminphase = minphase
         if maxphase gt allmaxphase then allmaxphase = maxphase
         minrv = min(rv.rv,max=maxrv)
         if minrv lt allminrv then allminrv = minrv
         if maxrv gt allmaxrv then allmaxrv = maxrv
      endfor
      
;   xtitle1='!3' + exofast_textoidl('Phase + (T_P - T_C)/P + 0.25',font=font)
;   xtitle2='!3' + exofast_textoidl('BJD_{TDB} - ' + string(bjd0,format='(i7)'),font=font)
;   plot, [0], [0], xtitle=xtitle1, psym=8, ytitle='!3RV (m/s)', xrange=xrange, yrange=yrange
      
      plot, [0], [0], xrange=[allminphase,allmaxphase],/xstyle,$
            yrange=[allminrv,allmaxrv], $
            ytitle='!3RV (m/s)', position=position1,xtickformat='(A1)'
      
      oplot, prettyphase, prettymodel,  color=red
      
      use = where(legendndx[i,*],nuse)
      if nuse gt 1 then exofast_legend, ss.telescope[use].label, color=colors[(indgen(ss.ntel) mod ncolors)[use]],/bottom,/right,psym=symbols[(indgen(ss.ntel) mod nsymbols)[use]], /useplotsym, charsize=0.5, fill=fills[(indgen(ss.ntel) mod nfills)[use]]
      
      for j=0L, ss.ntel-1 do begin
         rv = *(ss.telescope[j].rvptrs)
         if rv.planet eq i then begin
            
            ;; this needs to be debugged
            rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                               a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                               period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                               omega=ss.planet[i].omega.value[ndx]+!dpi,$
                               c=ss.constants.c/ss.constants.au*ss.constants.day,$
                               q=q[i])
            
            ;; calculate the RV model
            modelrv = exofast_rv(rvbjd,ss.planet[i].tp.value[ndx],ss.planet[i].period.value[ndx],$
                                 0d0,ss.planet[i].K.value[ndx]*q[i],$
                                 ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx]+!dpi,$
                                 slope=0d0)
            
            ;; for residual plot
            if tag_exist(ss.telescope[j],'detrend') then begin
               detrendadd = total((*ss.telescope[j].rvptrs).detrendadd*(replicate(1d0,n_elements((*ss.telescope[j].rvptrs).bjd))##(*ss.telescope[j].rvptrs).detrendaddpars.value),1)
               detrendmult = (1d0+total((*ss.telescope[j].rvptrs).detrendmult*(replicate(1d0,n_elements((*ss.telescope[j].rvptrs).bjd))##(*ss.telescope[j].rvptrs).detrendmultpars.value),1))
            endif else begin
               detrendadd = 0d0
               detrendmult = 1d0
            endelse
            (*(ss.telescope[j].rvptrs)).residuals = (rv.rv-modelrv-ss.telescope[j].gamma.value[ndx]-detrendadd)/detrendmult
            minoc = min((*(ss.telescope[j].rvptrs)).residuals ,max=maxoc)
            
            if minoc lt allminoc then allminoc = minoc
            if maxoc gt allmaxoc then allmaxoc = maxoc
            
            phasetime = ((( rvbjd - ss.planet[i].tc.value[ndx]) mod $
                          ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                         1.25d0) mod 1
            
            plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
            oplot, phasetime, modelrv, psym=8
            
         endif
      endfor
      
      ;; **** this is a hack; do this better ****
      if finite(allminoc) and finite(allmaxoc) then begin
         ;; residual plot
         plot, [0],[0], position=position2, /noerase, $
               xrange=[allminphase,allmaxphase], xtitle=xtitle1,$
               yrange=[allminoc,allmaxoc]/0.7d0, ytitle='O-C (m/s)', $
               /xstyle, /ystyle, yminor=2,yticks=2, ytickv=[allminoc,0,allmaxoc]
         oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  
      endif
      
      for j=0L, ss.ntel-1 do begin
         rv = *(ss.telescope[j].rvptrs)
         
         if rv.planet eq i then begin
            
            rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                               a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                               period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                               omega=ss.planet[i].omega.value[ndx]+!dpi,$
                               c=ss.constants.c/ss.constants.au*ss.constants.day,$
                               q=q[i])
            phasetime = ((( rvbjd - ss.planet[i].tc.value[ndx]) mod $
                          ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                         1.25d0) mod 1
            plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
            oplot , phasetime, rv.residuals, psym=8
         endif
      endfor
      
   endfor
endif

if keyword_set(psname) then begin
   device, /close
   set_plot, mydevice
   exofast_fixps, psname
endif

end
