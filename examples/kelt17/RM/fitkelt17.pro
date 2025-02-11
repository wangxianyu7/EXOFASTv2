pro fitkelt17,debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt17','RM'])

;; use this initial run to refine priors with kelt17.priors as a
;; user-supplied guess and Parallel tempering to search widely
;exofastv2, nplanets=1, tranpath=path+'n20??????.*.dat',$
;           dtpath=path+'n201?????.KELT-17b.TRES.44000.fits',$
;           mistsedfile=path+'kelt17.sed',rvpath =path+'KELT-17b.*.rv',$
;           priorfile=path + 'kelt17.priors',prefix=path + 'fitresults/KELT-17b.pt.',$
;           debug=debug, verbose=verbose, maxsteps=5000,nthin=1, fitdt=[1], ntemps=8
;
;; regenerate the prior file based on the best fit found above with
;; parallel tempering
;mkprior,
;filename=path+'fitresults/KELT-17b.pt.mcmc.idl',priorfilename=path +
;'kelt17.priors.2'

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','kelt17','fitresults'])

;; rerun the fit starting at the previous best fit
exofastv2, nplanets=1, tranpath=path + 'n20??????.*.dat',$
           mistsedfile=path+'kelt17.sed',rvpath=path+'KELT-17b.*.rv',$
           priorfile=path + 'kelt17.priors',$
           prefix=outpath + 'KELT-17b.',rmbands=["notrm","V","V"],rossiter=[1],restorebest=1, $
           debug=debug, verbose=verbose, maxsteps=maxsteps,nthin=nthin, nthread=nthread

end
