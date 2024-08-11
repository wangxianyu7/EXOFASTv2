pro resume_exofastv2, savfile, priorfile=priorfile, $
               prefix=prefix,$
               ;; data file inputs
               rvpath=rvpath, tranpath=tranpath, $
               astrompath=astrompath, dtpath=dtpath, $
               ;; SED model inputs
               fluxfile=fluxfile,mistsedfile=mistsedfile,$
               sedfile=sedfile,specphotpath=specphotpath,$
               noavprior=noavprior,$
               fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,$
               fehsedfloor=fehsedfloor, oned=oned,$
               ;; evolutionary model inputs
               yy=yy, nomist=nomist, parsec=parsec, $
               torres=torres, mannrad=mannrad,mannmass=mannmass, $
               teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
               rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
               ;; BEER model inputs
               fitthermal=fitthermal, fitellip=fitellip, $
               fitreflect=fitreflect, fitphase=fitphase, $
               fitbeam=fitbeam, derivebeam=derivebeam, $
               ;; star inputs
               nstars=nstars, starndx=starndx, $
               seddeblend=seddeblend, fitdilute=fitdilute, $
               ;; planet inputs
               nplanets=nplanets, $
               fittran=fittran, fitrv=fitrv, $
               rossiter=rossiter, fitdt=fitdt, $
               circular=circular, tides=tides, $ 
               alloworbitcrossing=alloworbitcrossing, $
               chen=chen, i180=i180, $
               ;; RV inputs
               fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch,$
               ;; transit inputs
               noclaret=noclaret, $
               ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs,$
               longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
               rejectflatmodel=rejectflatmodel,$
               noprimary=noprimary, requiresecondary=requiresecondary,$
               fitspline=fitspline, splinespace=splinespace, $
               fitramp=fitramp, fitwavelet=fitwavelet, $
               ;; reparameterization inputs
               fitlogmp=fitlogmp,$
               novcve=novcve, nochord=nochord, fitsign=fitsign, $
               fittt=fittt, earth=earth, $
               ;; plotting inputs
               transitrange=transitrange,rvrange=rvrange,$
               sedrange=sedrange,emrange=emrange, $
               ;; debugging inputs
               debug=debug, verbose=verbose, delay=delay, $
               ;; MCMC inputs
               maxsteps=maxsteps, nthin=nthin, maxtime=maxtime, $
               maxgr=maxgr, mintz=mintz, $
               dontstop=dontstop, $
               ntemps=ntemps, tf=tf, keephot=keephot, $
               randomfunc=randomfunc, seed=seed,$
               stretch=stretch, $
               nthreads=nthreads, $              
               ;; General inputs
               skiptt=skiptt, $
               usernote=usernote, $
               mksummarypg=mksummarypg,$
               nocovar=nocovar, $
               plotonly=plotonly, bestonly=bestonly, $
               badstart=badstart
               





;; if a virtual machine or runtime license, read the arguments from args.txt
if lmgr(/vm) or lmgr(/runtime) then begin
   ;; IDL_IDLBridge disabled in the virtual machine
   ;; no multi-threading without a license :(
   nthreads = 1L 

   par = command_line_args(count=numargs)
   if numargs eq 1 then begin
      argfile = par[0]
   endif else argfile = 'args.txt'

   if not file_test(argfile) then message, argfile + ', containing desired arguments to EXOFASTv2, does not exist'
   readargs, argfile, priorfile=priorfile, $
             prefix=prefix,$
             rvpath=rvpath, tranpath=tranpath, $
             astrompath=astrompath, dtpath=dtpath, $
             fluxfile=fluxfile,mistsedfile=mistsedfile,$
             sedfile=sedfile,specphotpath=specphotpath,$
             fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor, $
             fehsedfloor=fehsedfloor, oned=oned,$
             yy=yy, nomist=nomist, parsec=parsec, $
             torres=torres, mannrad=mannrad, mannmass=mannmass, $
             teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
             rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
             fitthermal=fitthermal, fitellip=fitellip, $
             fitreflect=fitreflect, fitphase=fitphase, $
             fitbeam=fitbeam, derivebeam=derivebeam, $
             nstars=nstars,starndx=starndx, $
             seddeblend=seddeblend,fitdilute=fitdilute, $
             nplanets=nplanets, $
             fittran=fittran, fitrv=fitrv, $
             rossiter=rossiter, fitdt=fitdt, $
             circular=circular, tides=tides, $ 
             alloworbitcrossing=alloworbitcrossing, $
             chen=chen, i180=i180, $
             fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch, $
             noclaret=noclaret, $
             ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
             longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
             rejectflatmodel=rejectflatmodel,$
             noprimary=noprimary, requiresecondary=requiresecondary,$
             fitspline=fitspline, splinespace=splinespace, $
             fitramp=fitramp, fitwavelet=fitwavelet, $              
             fitlogmp=fitlogmp,$
             novcve=novcve, nochord=nochord, fitsign=fitsign, $
             fittt=fittt, earth=earth, $             
             transitrange=transitrange,rvrange=rvrange,$
             sedrange=sedrange,emrange=emrange, $
             debug=debug, verbose=verbose, delay=delay,$
             maxsteps=maxsteps, nthin=nthin, maxtime=maxtime, $
             maxgr=maxgr, mintz=mintz, $
             dontstop=dontstop, $              
             ntemps=ntemps,tf=tf,keephot=keephot,$
             randomfunc=randomfunc, seed=seed,$
             stretch=stretch,$
             skiptt=skiptt, $
             usernote=usernote,$
             mksummarypg=mksummarypg,$
             nocovar=nocovar,$
             plotonly=plotonly, bestonly=bestonly,$
             logname=logname

endif


;; Load the save file, get the parameters, and set up the MCMC
restore, savfile
print, 'restored ', savfile
pars = str2parsarr(mcmcss)
npars = n_elements(pars[*,0])
pars = reform(pars,npars,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)
nfit = n_elements((*(mcmcss.tofit))[0,*])
previous_nsteps = mcmcss.nsteps/mcmcss.nchains
newsteps = maxsteps - previous_nsteps
combined_array = MAKE_ARRAY(npars, maxsteps, mcmcss.nchains)
newchi2 = dblarr(maxsteps, mcmcss.nchains)
nsteps_real_old = mcmcss.nsteps/mcmcss.nchains
oldchi2 = reform(*(mcmcss.chi2),nsteps_real_old,mcmcss.nchains)

newchi2[0:previous_nsteps-1, *] = oldchi2
chi2 = newchi2

combined_array[*, 0:previous_nsteps-1, *] = pars
combined_array[*, previous_nsteps:maxsteps-1, *] = dblarr(npars, newsteps, mcmcss.nchains)
pars = combined_array
resumendx = mcmcss.nsteps/mcmcss.nchains
chi2func = mcmcss.chi2func


;; default prefix for all output files (filename without extension)
if n_elements(prefix) eq 0 then prefix = 'fitresults/planet.'
basename = file_basename(prefix)

;; output to log file and the screen
logname = prefix + 'log'
file_delete, logname, /allow_nonexistent
if n_elements(usernote) ne 0 then printandlog, usernote, logname

;; name of the chi square function
; chi2func = 'exofast_chi2v2'

;; compile all routines now to keep output legible 
;; resolve_all doesn't interpret execute; it's also broken prior to IDL v6.4(?)
defsysv, '!GDL', exists=runninggdl  

;; default to NCORES threads, if we're running a full copy of IDL
if runninggdl or lmgr(/runtime) or lmgr(/vm) then begin
   nthreads=1
endif else if n_elements(nthreads) ne 1 then begin
   nthreads = !cpu.hw_ncpu
endif ;; else use the user's input   

if double(!version.release) ge 6.4d0 and ~lmgr(/vm) and ~lmgr(/runtime) and ~runninggdl then $
   resolve_all, resolve_either=[chi2func,'exofast_random','ramp_func'],skip_routines=['cggreek'],/cont,/quiet

;; output to log file too
logname = prefix + 'log'
file_delete, logname, /allow_nonexistent


exofast_path = getenv('EXOFAST_PATH')
;if (strpos(exofast_path,'//') ne -1) or 
if (strpos(exofast_path,'EXOFASTv2') eq -1) then begin
   printandlog, "ERROR: EXOFAST_PATH (" + exofast_path + ") not set properly",logname
   printandlog, "Typically, EXOFAST_PATH is something like ${HOME}/idl/EXOFASTv2/", logname
;   printandlog, "NOTE: '//' is not allowed as it is inconsistenly resolved", logname
   return
endif

;; if the directory doesn't exist, make it
dirname = file_dirname(prefix)
if dirname ne '.' then begin
   if ~file_test(dirname,/directory) then begin
      file_mkdir, dirname
   endif
endif

modeldirname = dirname + path_sep() + 'modelfiles'
if ~file_test(modeldirname,/directory) then begin
   file_mkdir, modeldirname
endif

;; some error checking on NTHREADS (multi-threading)
if nthreads gt !cpu.hw_ncpu then begin
   printandlog, "WARNING: Using more threads (" + strtrim(nthreads,2) + ") than physical cores (" + strtrim(!cpu.hw_ncpu,2) + "); this is likely to be inefficient", logname
;endif else if nthreads eq !cpu.hw_ncpu then begin
;   printandlog, "WARNING: Doing fit with " + strtrim(nthreads,2) + " threads. This may impact performance of other tasks (See NTHREADS argument).", logname
endif else begin
   printandlog, "Doing fit with " + strtrim(nthreads,2) + ' threads', logname
endelse

;; insert the commit id into the log
spawn, 'git -C $EXOFAST_PATH rev-parse HEAD', output, stderr
if stderr[0] eq '' then printandlog, "Using EXOFASTv2 commit " + output[0], logname

;; create the master structure 
;; keyword inheritance would be helpful here, but it might break multi-threading
ss = mkss(priorfile=priorfile, $
          prefix=prefix,$
          ;; data file inputs
          rvpath=rvpath, tranpath=tranpath, $
          astrompath=astrompath, dtpath=dtpath, $
          ;; SED model inputs
          fluxfile=fluxfile, mistsedfile=mistsedfile, $
          sedfile=sedfile, specphotpath=specphotpath,$
          noavprior=noavprior,$
          fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,$
          fehsedfloor=fehsedfloor, oned=oned,$
          ;; evolutionary model inputs
          yy=yy, nomist=nomist, parsec=parsec, $
          torres=torres, mannrad=mannrad, mannmass=mannmass,$
          teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
          rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
          ;; BEER model inputs
          fitthermal=fitthermal, fitellip=fitellip, $
          fitreflect=fitreflect, fitphase=fitphase,$
          fitbeam=fitbeam, derivebeam=derivebeam, $
          ;; star inputs
          nstars=nstars, starndx=starndx, $
          seddeblend=seddeblend, fitdilute=fitdilute, $
          ;; planet inputs
          nplanets=nplanets, $
          fittran=fittran,fitrv=fitrv,$
          rossiter=rossiter, fitdt=fitdt,$ 
          circular=circular, tides=tides, $
          alloworbitcrossing=alloworbitcrossing,$
          chen=chen, i180=i180,$
          ;; RV inputs
          fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch, $
          ;; transit inputs
          noclaret=noclaret,$
          ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
          longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
          rejectflatmodel=rejectflatmodel,$
          noprimary=noprimary, requiresecondary=requiresecondary,$
          fitspline=fitspline, splinespace=splinespace, $
          fitramp=fitramp, fitwavelet=fitwavelet, $            
          ;; reparameterization inputs
          fitlogmp=fitlogmp,$
          novcve=novcve, nochord=nochord, fitsign=fitsign, $
          fittt=fittt, earth=earth, $
          ;; plotting inputs
          transitrange=transitrange,rvrange=rvrange,$
          sedrange=sedrange,emrange=emrange, $
          ;; debugging inputs
          debug=debug, verbose=verbose, delay=delay, $
          ;; internal inputs
          chi2func=chi2func, $
          logname=logname)

if (size(ss))[2] ne 8 then begin
   badstart=1
   return
endif else badstart=0

if (size(ss))[2] ne 8 then begin
   badstart=1
   return
endif else badstart=0

npars = ss.npars
nfit = n_elements((*(ss.tofit))[0,*])

;; this is where threads were previously initialized

if n_elements(mintz) eq 0 then mintz = 1000d0
if n_elements(maxgr) eq 0 then maxgr = 1.01d0

; ;; by default, use 1 GB of RAM and set NTHIN so it converges
; if n_elements(maxsteps) eq 0 then begin
;    maxsteps = round(1024d0^3/(double(ss.nchains)*npars*8d0))
   
;    ;; set to 10x the limit expected from an idealized Gaussian case of N parameters
;    if n_elements(nthin) eq 0 then begin
;       a = [173881.7852504589d0,-3934.6300584134d0,512.5700554749d0]
;       nthin = ceil((a[0]+a[1]*nfit+a[2]*nfit^2)/(2d0*nfit)/maxsteps*10d0) > 1
;       printandlog, 'NTHIN set to ' + strtrim(nthin,2) + ' for the fit to converge.', logname
;       printandlog, 'This is an extremely rough estimate and varies wildly depending on the ',logname
;       printandlog, 'details of the fit and the accuracy of the starting conditions.',logname
;       printandlog, 'You should monitor the progress of the fit and increase NTHIN if necessary.', logname
;       printandlog, '', logname
;    endif
; endif 

; memrequired = double(ss.nchains)*double(maxsteps)*npars*8d0/(1024d0^3)
; printandlog, 'Fit will require ' + strtrim(memrequired,2) + ' GB of RAM for the final structure', logname
; if memrequired gt 2d0 then begin
;    printandlog, 'WARNING: this likely exceeds your available RAM and may crash after the end of a very long run. You likely want to reduce MAXSTEPS and increase NTHIN by the same factor. If you would like to proceed anyway, type ".con" to continue', logname
;    if ~lmgr(/vm) then stop
; endif
; printandlog, '', logname

if n_elements(nthin) eq 0 then nthin = 1L
if nthin lt 1L then nthin=1L

printandlog, 'MAXSTEPS set to ' + strtrim(maxsteps,2), logname
printandlog, 'NTHIN set to ' + strtrim(nthin,2), logname

; pars = str2pars(ss,scale=scale,name=name, angular=angular)

;; plot the data + starting guess
; modelfile = prefix + 'start'
ss.verbose=1B
; bestchi2 = call_function(chi2func, pars, psname=modelfile)
; ss.verbose = keyword_set(verbose)
; if ~finite(bestchi2) then begin
;    printandlog, 'Starting model is out of bounds; cannot recover. You must change the starting parameter(s) via the prior file.', logname
;    printandlog, 'Re-running starting model with /VERBOSE flag to identify the parameter', logname
;    ss.verbose=1B
;    bestchi2 = call_function(chi2func, pars, psname=modelfile)
;    printandlog, 'Starting model is out of bounds; cannot recover. You must change the starting parameter(s) via the prior file.', logname
;    badstart=1
;    return
; endif else begin
;    printandlog, 'The loglike of the starting model was ' + strtrim(-bestchi2/2d0,2), logname
; endelse
; if keyword_set(plotonly) then return

;; do it again for accurate timing 
;; after loading all the files into cache, not including plotting
; t0 = systime(/seconds)
; bestchi2 = call_function(chi2func, pars)
; modeltime = systime(/seconds)-t0

;; these are the starting values for all step parameters
; printandlog, 'These are the starting values for all fitted parameters', logname
; printandlog, 'and the amoeba stepping scale, which is roughly the range', logname
; printandlog, 'of parameter space it will explore around the starting value', logname
; printandlog, 'and is equal to 3x any Gaussian width. When priors are not ', logname
; printandlog, 'specified in ' + priorfile + ', a default guess is used.', logname
; printandlog, 'The parameter number is useful for tracking down unconstrained',logname
; printandlog, 'parameters', logname
; printandlog,'',logname 
; printandlog, '**************************************************************', logname
; printandlog, '*** IT IS WISE TO MAKE SURE THESE AGREE WITH YOUR          ***', logname
; printandlog, '*** EXPECTATION FROM THE PRIORFILE. IF NOT, YOUR           ***', logname
; printandlog, '*** STARTING PRIORS MAY NOT BE TRANSLATED CORRECTLY INTO   ***', logname
; printandlog, '*** THE FITTED PARAMETERIZATION BY                         ***', logname
; printandlog, '*** $EXOFAST_PATH/pars2step.pro. NOT ALL PARAMETER         ***', logname
; printandlog, '*** COMBINATIONS ARE ALLOWED/SUPPORTED.                    ***', logname
; printandlog, '*** WHEN IN DOUBT, SET PRIORS/CHANGE STARTING VALUES       ***', logname
; printandlog, '*** DIRECTLY IN THESE FITTED PARAMETERS.                   ***', logname
; printandlog, '**************************************************************', logname
; printandlog, '', logname
; printandlog, 'Par #      Par Name    Par Value       Amoeba Scale', logname
; for i=0, n_elements(name)-1 do printandlog, string(i, name[i], pars[i], scale[i], format='(i3,x,a15,x,f14.6,x,f14.6)'), logname
; printandlog, '', logname

; if ss.debug and ~lmgr(/vm) then begin
;    printandlog, 'program halted to give you time to inspect the priors. Type ".con" to continue', logname
;    stop
; end

; nmax = 1d5
; printandlog, 'It takes ' + strtrim(modeltime,2) + ' seconds to calculate a single model', logname
; printandlog, 'Beginning AMOEBA fit; this may take up to ' + string(modeltime*nmax/60d0,format='(f0.1)') + ' minutes if it takes the maximum allowed steps (' + strtrim(nmax,2) + ')', logname

; ;; do the AMOEBA fit
ss.amoeba = 1B
ss.delay =0
; best = exofast_amoeba(1d-5,function_name=chi2func,p0=pars,scale=scale,nmax=nmax)

; ss.delay = delay
; if best[0] eq -1 then begin
;    printandlog, 'ERROR: Could not find best combined fit; adjust your starting values and try again. You may want to set the /DEBUG keyword.', logname
;    return
; endif
; printandlog, 'Finished AMOEBA fit', logname
; save, best, filename=prefix + 'amoeba.idl'

;; update the parameter array with the chen-derived logks/rp (is this necessary?)
;best = str2pars(ss,scale=scale,name=name) 

;; try again?
;printandlog, 'restarting AMOEBA with chen enabled', logname
;printandlog, call_function(chi2func,best,modelrv=modelrv,modelflux=modelflux, psname=prefix + 'model'), logname
;best = exofast_amoeba(1d-8,function_name=chi2func,p0=best,scale=scale,nmax=nmax)
;printandlog, call_function(chi2func,best,modelrv=modelrv,modelflux=modelflux, psname=prefix + 'model'), logname

;; output the best-fit model fluxes/rvs
; bestchi2 = call_function(chi2func,best,modelrv=modelrv,modelflux=modelflux, psname=prefix + 'amoeba')
; printandlog, 'The best loglike found by AMOEBA was ' + strtrim(-bestchi2/2d0,2), logname
; printandlog, 'It should only be compared against the loglike of the same model with different starting points', logname

;; initialize the threads
if nthreads gt 1 then begin
   printandlog, 'Initializing ' + strtrim(nthreads,2) + ' threads', logname

   ;; load the stellar structure into the common block for each thread
   ;; (can't pass structures between threads, so we have to create it in each thread)
   thread_array = replicate(create_struct('obridge',obj_new("IDL_IDLBridge", output=''),$
                                          'j',-1L, 'm',-1L, 'k', -1L, $
                                          'newpars',dblarr(nfit),'status',0B, 'fac', 1d0,$
                                          'start', systime(/seconds)),nthreads)
   
   ;; get the current working directory (make sure all threads start there)
   cd, current=cwd

   for i=0L, nthreads-1 do begin
      ;; replicate copies the IDLBridge by reference, not by value!! reinitialize here
      thread_array[i].obridge = obj_new("IDL_IDLBridge", output='')
      
      ;; can't share a structure directly with a thread
      ;; must share components and create the structure within the thread
      ;; share every variable in memory to make it future proof
      help, output=helpoutput
      for j=0L, n_elements(helpoutput)-2 do begin
         ;; done with variables, we're done
         if helpoutput[j] eq 'Compiled Procedures:' then break
         
         ;; undefined variable; skip it
         if strpos(helpoutput[j],'UNDEFINED') eq 16 then continue

         entries = strsplit(helpoutput[j],/extract)

         ;; either a long variable name that spans two lines, or not a variable
         if strpos(helpoutput[j],'=') ne 26 then begin

            ;; not a variable; skip it
            if strpos(helpoutput[j+1],'=') ne 26 then continue

            ;; undefined variable, skip it
            if strpos(helpoutput[j+1],'UNDEFINED') eq 16 then continue

            ;; if it's not a lone (long) variable name, skip it
            if n_elements(entries) ne 1 then continue

            ;; this line is the name of a long variable name
            ;; skip the next line, which is its value
            j++

         endif else if n_elements(entries) gt 4 then begin
            ;; catch N-dimentional arrays and strings with spaces
            if strpos(helpoutput[j],'Array') ne 28 and strpos(helpoutput[j],'STRING') ne 16 then continue
         endif else if n_elements(entries) ne 4 then continue            
         ;; declare it in the thread
         ;; EXECUTE is ok here, since we can't use VM with IDLBridge anyway
         varname = entries[0]
         junk = execute("thread_array[i].obridge->setvar,'" + varname + "'," + varname)
;         print, 'setting ' + varname + ' to '
;         junk = execute('print, ' + strtrim(varname,2))
      endfor

      ;; disable NaN warnings inside each thread
      thread_array[i].obridge->setvar,'!except',0

      ;; make sure each thread is run from the current working directory
      thread_array[i].obridge->setvar,'cwd',cwd
      thread_array[i].obridge->execute,'cd, cwd'
      
      ;; compile all the codes in each thread so compilation messages don't pollute the screen
      if double(!version.release) ge 6.4d0 and ~lmgr(/vm) and ~lmgr(/runtime) and ~runninggdl then $
         thread_array[i].obridge->execute, "resolve_all, resolve_either=[chi2func,'exofast_random','ramp_func'], resolve_procedure=['exofastv2'],skip_routines=['cggreek'],/cont,/quiet"
      
      ; create the stellar stucture within each thread
      thread_array[i].obridge->execute,$
         'ss = mkss(priorfile=priorfile, prefix=prefix,'+$
         'rvpath=rvpath, tranpath=tranpath,'+$
         'astrompath=astrompath, dtpath=dtpath,'+$
         'fluxfile=fluxfile, mistsedfile=mistsedfile,'+$
         'sedfile=sedfile,specphotpath=specphotpath,'+$
         'noavprior=noavprior,'+$
         'fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,'+$
         'fehsedfloor=fehsedfloor, oned=oned,'+$
         'yy=yy, nomist=nomist, parsec=parsec,'+ $
         'torres=torres, mannrad=mannrad, mannmass=mannmass,'+$         
         'teffemfloor=teffemfloor, fehemfloor=fehemfloor,'+$
         'rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,'+$
         'fitthermal=fitthermal, fitellip=fitellip,'+ $
         'fitreflect=fitreflect,fitphase=fitphase,'+ $
         'fitbeam=fitbeam, derivebeam=derivebeam,'+ $
         'nstars=nstars,starndx=starndx,'+ $         
         'seddeblend=seddeblend, fitdilute=fitdilute,'+$
         'nplanets=nplanets,'+$
         'fittran=fittran, fitrv=fitrv,'+$
         'rossiter=rossiter,fitdt=fitdt,'+$
         'circular=circular,tides=tides,'+$
         'alloworbitcrossing=alloworbitcrossing,'+$
         'chen=chen, i180=i180,'+$
         'fitslope=fitslope, fitquad=fitquad,rvepoch=rvepoch,'+$
         'noclaret=noclaret,'+$
         'ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs,'+$
         'longcadence=longcadence,exptime=exptime,ninterp=ninterp,'+$ 
         'rejectflatmodel=rejectflatmodel,'+$
         'noprimary=noprimary, requiresecondary=requiresecondary,'+$
         'fitspline=fitspline, splinespace=splinespace,'+$
         'fitramp=fitramp, fitwavelet=fitwavelet,'+$
         'fitlogmp=fitlogmp,'+$
         'novcve=novcve, nochord=nochord, fitsign=fitsign,'+$
         'fittt=fittt, earth=earth,'+$
         'transitrange=transitrange,rvrange=rvrange,'+$
         'sedrange=sedrange,emrange=emrange,'+$
         'debug=debug, verbose=verbose,delay=delay,'+$
         '/silent,'+$
         'chi2func=chi2func,'+$
         'logname=logname)'
   endfor
endif

;; do the MCMC fit

resumendx = ss.nsteps/ss.nchains
best = pars[*,0,0]
bestchi2 = call_function(chi2func, best)
testpars = str2pars(ss,scale=scale,name=name, angular=angular)
if n_elements(testpars) ne n_elements(best) then begin
   print, 'Pealse check the starting values. The starting values are not consistent with the existing mcmc.idl.'
   stop
endif
print,'Running MCMC',
exofast_demcpt_multi, best, chi2func, pars, chi2=chi2,$
                        nthin=nthin,maxsteps=maxsteps, maxtime=maxtime, $
                        ntemps=ntemps, tf=tf, dontstop=dontstop, $
                        burnndx=burnndx, goodchains=goodchains, seed=seed, randomfunc=randomfunc, $
                        gelmanrubin=gelmanrubin, tz=tz, maxgr=maxgr, mintz=mintz, $
                        stretch=stretch, logname=logname, angular=angular, $
                        keephot=keephot, hotpars=hotpars, hotchi2=hotchi2, thread_array=thread_array,resumendx=resumendx

if pars[0] eq -1 then begin
   printandlog, 'MCMC Failed to find a stepping scale. This usually means one or more parameters are unconstrained by the data or priors.', logname
endif

bad = where(tz lt mintz or gelmanrubin gt maxgr,nbad)
if bad[0] ne -1 then begin
   printandlog, 'WARNING: The Gelman-Rubin statistic indicates ' + $
                  'the following parameters are not well-mixed', logname
   printandlog, '    Parameter   Rz     Tz', logname
   for i=0, nbad-1 do printandlog, string(name[bad[i]], gelmanrubin[bad[i]],tz[bad[i]], format='(a13,x,2(f0.4,x))'), logname
endif
printandlog, 'Synthesizing results; for long chains and/or many fitted parameters, this may take up to 15 minutes', logname

;; combine all chains
sz = size(pars)
npars = sz[1]
nsteps = sz[2]
nchains = sz[3]
pars = reform(pars,npars,nsteps*nchains)
chi2 = reform(chi2,nsteps*nchains)
minchi2 = min(chi2,bestndx)

printandlog, 'The best loglike found by MCMC was ' + strtrim(-minchi2/2d0,2), logname
printandlog, 'It should only be compared against the loglike of the same model with different starting points', logname
printandlog, '', logname
printandlog, 'Use BIC and AIC to compare different models', logname
bic = nfit*alog(ss.ndata) + minchi2
aic = 2d0*nfit + minchi2
printandlog, 'NDATA = ' + strtrim(ss.ndata,2),logname
printandlog, 'NFIT = ' + strtrim(nfit,2),logname
printandlog, 'BIC = ' + strtrim(bic,2),logname
printandlog, 'AIC = ' + strtrim(aic,2),logname
if minchi2 lt bestchi2 then begin
   printandlog, 'WARNING: MCMC found a better model that AMOEBA.', logname
   printandlog, 'Using mkprior to refine your starting values and refitting', logname
   printandlog, 'may result in a faster and more robust answer.', logname
   printandlog, 'Look at the chain plot before you trust this fit.', logname
endif


;; generate the model fit from the best MCMC values, not AMOEBA
bestamoeba = best
best = pars[*,bestndx]
modelfile = prefix + 'mcmc'
ss.verbose = 1
bestchi2 = call_function(chi2func,best,psname=modelfile, $
                         modelrv=modelrv, modelflux=modelflux)
ss.verbose = keyword_set(verbose)

;; make a new stellar system structure with only fitted and derived
;; parameters, populated by the pars array
;mcmcss = mcmc2str(pars, ss)
mcmcss = mkss(priorfile=priorfile, $
              prefix=prefix,$
              ;; data file inputs
              rvpath=rvpath, tranpath=tranpath, $
              astrompath=astrompath, dtpath=dtpath, $
              ;; SED model inputs
              fluxfile=fluxfile, mistsedfile=mistsedfile, $
              sedfile=sedfile, specphotpath=specphotpath,$
              noavprior=noavprior,$
              fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,$
              fehsedfloor=fehsedfloor, oned=oned,$
              ;; evolutionary model inputs
              yy=yy, nomist=nomist, parsec=parsec, $
              torres=torres, mannrad=mannrad, mannmass=mannmass,$
              teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
              rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
              ;; BEER model inputs
              fitthermal=fitthermal, fitellip=fitellip, $
              fitreflect=fitreflect, fitphase=fitphase,$
              fitbeam=fitbeam, derivebeam=derivebeam, $
              ;; star inputs
              nstars=nstars, starndx=starndx, $
              seddeblend=seddeblend, fitdilute=fitdilute, $
              ;; planet inputs
              nplanets=nplanets, $
              fittran=fittran,fitrv=fitrv,$
              rossiter=rossiter, fitdt=fitdt,$ 
              circular=circular, tides=tides, $
              alloworbitcrossing=alloworbitcrossing,$
              chen=chen, i180=i180,$
              ;; RV inputs
              fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch, $
              ;; transit inputs
              noclaret=noclaret,$
              ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
              longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
              rejectflatmodel=rejectflatmodel,$
              noprimary=noprimary, requiresecondary=requiresecondary,$
              fitspline=fitspline, splinespace=splinespace, $
              fitramp=fitramp, fitwavelet=fitwavelet, $            
              ;; reparameterization inputs
              fitlogmp=fitlogmp,$
              novcve=novcve, nochord=nochord, fitsign=fitsign, $
              fittt=fittt, earth=earth, $
              ;; plotting inputs
              transitrange=transitrange,rvrange=rvrange,$
              sedrange=sedrange,emrange=emrange, $
              ;; debugging inputs
              debug=debug, verbose=verbose, delay=delay, $
              ;; internal inputs
              nvalues=nsteps*nchains,$ 
              /silent, $
              chi2func=chi2func, $
              logname=logname, $
              best=best)

if (size(mcmcss))[2] ne 8 then return

mcmcss.nchains = nchains
mcmcss.burnndx = burnndx
*(mcmcss.goodchains) = goodchains
*(mcmcss.chi2) = chi2

pars2str, pars, mcmcss

;; populate residuals for the mcmcss file
for i=0L, mcmcss.ntran-1 do $
   (*mcmcss.transit[i].transitptrs).residuals = (*ss.transit[i].transitptrs).residuals
for i=0L, mcmcss.ntel-1 do $
   (*mcmcss.telescope[i].rvptrs).residuals = (*ss.telescope[i].rvptrs).residuals
   
;; derive all parameters
derivepars, mcmcss, logname=logname

spawn, 'git -C $EXOFAST_PATH rev-parse --short HEAD', output, stderr
if output[0] ne '' then versiontxt = ", created using EXOFASTv2 commit number " + output[0] $
else versiontxt = ''

;; output filenames
label = "tab:" + basename
caption = "Median values and 68\% confidence interval for " + basename + versiontxt
parfile = prefix + 'pdf.ps'
covarfile = prefix + 'covar.ps'
chainfile = prefix + 'chain.ps'
texfile = prefix + 'median.tex'
csvfile = prefix + 'median.csv'

exofast_plotdist_corner, mcmcss, pdfname=parfile, covarname=covarfile,nocovar=nocovar,logname=logname, csvfile=csvfile
;; save the chains for additional analysis 
;; wait until here because many things are populated in exofast_plotdist_corner
idlfile = prefix + 'mcmc.idl'

;; make a new prior file to start at the best fit found here
priorfileparts = strsplit(priorfile,'.',/extract,/preserve_null)
suffix = priorfileparts[n_elements(priorfileparts)-1]
if valid_num(suffix) then priorfileparts[n_elements(priorfileparts)-1] = strtrim(suffix+1L,2) $
else priorfileparts = [priorfileparts,'2']
priorfile2 = file_dirname(prefix) + path_sep() + file_basename(strjoin(priorfileparts,'.'))
mkprior2, mcmcss=mcmcss, priorfilename=priorfile2

;; GDL compatibility
if runninggdl then begin
   if keyword_set(keephot) and ntemps gt 1 then cmsave, mcmcss, hotpars, hotchi2, filename=idlfile $
   else cmsave, mcmcss, filename=idlfile
endif else begin
   if keyword_set(keephot) and ntemps gt 1 then save, mcmcss, hotpars, hotchi2, filename=idlfile $
   else save, mcmcss, filename=idlfile
endelse

exofast_latextab2, mcmcss, caption=caption, label=label,texfile=texfile
exofast_plotchains, mcmcss, chainfile=chainfile, logname=logname

if (total(mcmcss.ttvs) gt 0 or ~keyword_set(skiptt)) and mcmcss.ntran gt 0 then begin
   printandlog, 'The fit is done and can be interrupted without losing any results', logname
   printandlog, 'Now generating a table of the numerically solved times of ', logname
   printandlog, 'minimum projected separation, depth, and impact parameters for',logname
   printandlog, 'each transit file. This may take a while, but can be done at any',logname
   printandlog, 'point with the following command:', logname
   printandlog, "junk = exofast_gettt(filename='" + idlfile + "', filebase='" + prefix + "')", logname
   junk = exofast_gettt(mcmcss, filebase=prefix)

   if total(mcmcss.ttvs) gt 1 then begin
      ;; generate an O-C diagram for each planet
      readcol, prefix + 'transits.csv',label,planet,epoch,time,hierr,loerr, format='a,a,l,d,d,d', delimiter=',', comment='#',/silent
      err = (double(hierr) + double(loerr))/2d0
      telescope = label
      for i=0L, n_elements(label)-1 do telescope[i] = (strsplit(label[i],' UT ',/regex,/extract))[0]

      sorted = sort(planet)
      uniqplanets = planet[sorted[uniq(planet[sorted])]]
      for i=0L, n_elements(uniqplanets)-1 do begin
         match = where(planet eq uniqplanets[i]); and mcmcss.ttvs[*,i])
         omc, time[match], err[match], telescope=telescope[match], epsname=prefix + uniqplanets[i] + '.ttv.eps', $
              period=median(mcmcss.planet[i].period.value), t0=median(mcmcss.planet[i].tc.value),logname=logname
      endfor
   endif
   
   ;; generate a plot of the tdeltavs 
   ;; *** I think this breaks for multiple planets***
   if total(mcmcss.tdeltavs) gt 1 then begin
      plot_tdeltav, prefix + 'transits.csv'
   endif

endif

;; this makes a quick-look summary page, but requires gs, ps2pdf, grep,
;; convert and likely only works on linux
if keyword_set(mksummarypg) then begin
   mksummaryframe,idlfile=idlfile,base=prefix
endif

end