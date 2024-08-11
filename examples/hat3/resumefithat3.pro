pro resumefithat3, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, $
             nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hat3'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','hat3','fitresults'])

;; Fit using the Torres relation
resume_exofastv2, outpath +'HAT-3b.Torres'+'.mcmc.idl', nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', $
           rvpath=path+'HAT-3b.*.rv',priorfile=path+'hat3.torres.priors', $
           prefix=outpath +'HAT-3b.Torres.', maxsteps=20000, $
           nthin=nthin, circular=[1], fitrv=[1],fittran=[1], $
           debug=debug,ntemps=5, verbose=verbose,  nthread=2, /KEEPHOT


end
