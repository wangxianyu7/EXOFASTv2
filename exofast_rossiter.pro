FUNCTION exofast_rossiter, time, inc_rad, ar, TPeriastron, period, e, omega_rad, p, u1, u2, lambda, vsini, $
         vgamma, vzeta, vxi, valpha, $ 
         exptime=exptime, ninterp=ninterp
    ; name: exofast_rossiter
    ; purpose: calculate the radial velocity anomaly due to the Rossiter-McLaughlin effect: Ohta2005, Hirano2010, Hirano2011
    ; inputs:
    ;       time: time array in BJD
    ;       inc_rad: inclination angle in radians
    ;       ar: a/Rs
    ;       tc: time of conjunction
    ;       period: orbital period
    ;       e: eccentricity
    ;       omega_rad: argument of periastron in radians
    ;       p: planet to star radius ratio
    ;       u1: limb darkening parameter 1
    ;       u2: limb darkening parameter 2
    ;       lambda: sky-projected angle between the stellar spin and the planetary orbit
    ;       vsini: projected stellar rotational velocity
    
    ;       vbeta: the Guassian dispersion of spectral lines in m/s, typically 2500 - 4500 m/s (see Hirano+11)
    ;       vgamma: the Lorentzian dispersion of spectral lines in m/s, typically 500-1500 m/s (see Hirano+11)
    ;       zeta: the macroturbulence dispersion in m/s, typically 2000 - 6500 km/s (see Hirano+11)

    ;       exptime: exposure time in minutes
    ;       ninterp: number of points to interpolate
    ;       ohta2005: set to 1 to use Ohta+2005 model
    ;       hirano2010: set to 1 to use Hirano+2010 model
    ;       hirano2011: set to 1 to use Hirano+2011 model
    ;       if none of the models are set, Hirano+2011 model is used by default
    ; outputs:
    ;       delta_rv: radial velocity anomaly in m/s
    
   
    
    COMPILE_OPT IDL2

    vsini = vsini / 1000D0
    vgamma = vgamma / 1000D0
    vzeta = vzeta / 1000D0
    vxi = vxi / 1000D0
    valpha = valpha / 1000D0
    ; Get the number of points
    npoints = N_ELEMENTS(time)
    IF N_ELEMENTS(ninterp) EQ 0 THEN ninterp = 1

    ; force supersampling to be an integer

    exptime = abs(median(ts_diff(time,1)))*24*60D0
    IF exptime LT 3 THEN BEGIN
        ninterp = 0
    ENDIF ELSE BEGIN
        ninterp = ROUND(exptime/2D0)
    ENDELSE


    ; Interpolate if ninterp is greater than 1
    IF ninterp GT 1 THEN BEGIN
        transitbjd = time # (DBLARR(ninterp) + 1D0) + $
                     ((DINDGEN(ninterp) / ninterp - (ninterp - 1D0) / (2D0 * ninterp)) / $
                     1440D0 * exptime) ## (DBLARR(npoints) + 1D0)
        modelflux = DBLARR(npoints, ninterp) + 1D0
    ENDIF ELSE BEGIN
        transitbjd = time
        modelflux = DBLARR(npoints) + 1D0
    ENDELSE

    ; Get the model flux using exofast_tran function
    tmpmodelflux = exofast_tran(time, inc_rad, ar, TPeriastron, period, e, omega_rad, p, u1, u2, 1)

    ; Calculate phase and true anomaly
    meananom = 2.D0 * !dpi * (1.D0 + (transitbjd - Tperiastron) / Period MOD 1)
    eccanom = exofast_keplereq(meananom, e)
    trueanom = 2.D0 * ATAN(SQRT((1.D0 + e) / (1.D0 - e)) * TAN(eccanom / 2.D0))
    

    cosi = COS(inc_rad)
    ; Calculate up, and vp values
    b = exofast_getb2(time, i=inc_rad, a=ar, tperiastron=TPeriastron, period=period, e=e,omega=omega_rad,z2=z, x2=x,y2=y)
    xp = x*cos(lambda) - y*sin(lambda)
    yp = x*sin(lambda) + y*cos(lambda)
    vp = vsini * xp
    ; Compute flux and ratios
  
    f = 1 - tmpmodelflux
    ratios = DBLARR(npoints)
    IF vsini lt 1 THEN BEGIN
        ; PRINT, 'Error: vsini cannot be zero'
        RETURN, f*0
    ENDIF

    max_group_size = 200     ; Maximum number of pairs per group
    delta_rv = DBLARR(npoints) ; Initialize result array
    
    group_count = CEIL(float(npoints) / float(max_group_size)) ; Number of groups
    
    ; Loop over each group
    FOR group = 0, group_count - 1 DO BEGIN
        start_idx = group * max_group_size
        end_idx = (group + 1) * max_group_size - 1
        end_idx = MIN([end_idx, npoints - 1])
        
        
        group_size = end_idx - start_idx + 1
        xy = DBLARR(group_size * 2) ; Prepare the xy array
        
        ; Fill xy array with interleaved xp and yp for this group
        xy[0:*:2] = STRING(xp[start_idx:end_idx])
        xy[1:*:2] = STRING(yp[start_idx:end_idx])
        zp = z[start_idx:end_idx]
        ; Parameters for this group
        pars = [u1, u2, vsini, p, vxi, vgamma, vzeta, valpha, 0, group_size]
        par_str = STRING(pars)
        
        ; Prepare the command
        run_input = [par_str, xy]
        cmd = '$EXOFAST_PATH/new_analytic7.exe ' + STRJOIN(run_input, ' ')
        ; Run the external program and capture the result
        SPAWN, cmd, result
        
        ; Convert result to delta_rv and store in the appropriate indices
        delta_rv_tmp = FLOAT(result) * 1000D0

        ; Check if any NaN exists in the result
        IF (WHERE(FINITE(delta_rv_tmp) EQ 0, count))[0] NE -1 THEN BEGIN
            PRINT, 'Error: NaN detected in result. Returning delta_rv = INF'
            RETURN, DBLARR(npoints) + 0d0
        ENDIF

        delta_rv[start_idx:end_idx] = delta_rv_tmp

        for i = 0, group_size - 1 DO BEGIN
        if zp[i] GT 0 THEN delta_rv[start_idx + i] = 0
        endfor

    ENDFOR
    RETURN, delta_rv
END


