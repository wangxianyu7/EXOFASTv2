FUNCTION Theta_sigma, sigma, zeta, cos_theta, sin_theta
    cos_component = EXP(-(!PI * zeta * cos_theta)^2 * sigma^2)
    sin_component = EXP(-(!PI * zeta * sin_theta)^2 * sigma^2)
    RETURN, 0.5 * (cos_component + sin_component)
END

FUNCTION integrated_J0, t, vsini, sigma, phi_space
    y = COS(2 * !PI * sigma * vsini * t * COS(phi_space))
    RETURN, TSUM(phi_space,y) / (2.*!PI)  
END

FUNCTION M_sigma, t, vsini, sigma, u1, u2, zeta
    J0_value = BESELJ(2 * !PI * vsini * sigma * t,0)
    numerator = (1 - u1 * (1 - SQRT(1 - t^2)) - u2 * (1 - SQRT(1 - t^2))^2)
    denominator = (1 - u1 / 3 - u2 / 6)
    exp_term = EXP(-!PI^2 * zeta^2 * sigma^2 * (1 - t^2)) + EXP(-!PI^2 * zeta^2 * sigma^2 * t^2)
    RETURN, numerator / denominator * exp_term * J0_value * t
END


FUNCTION integrated_M_sigma, sigma, vsini, u1, u2, zeta, t_space
    y = M_sigma(t_space, vsini, sigma, u1, u2, zeta)
    RETURN, TSUM(t_space,y) ;
END



FUNCTION compute_integral, costheta, sintheta, f, vp, u1, u2, vsini, vbeta, vgamma, zeta
    ;COMPILE_OPT idl2
    ; Prepare the necessary arrays and grid spaces
    t_space = FINDGEN(1001) / 1000
    sigma_space = FINDGEN(101) / 100 / (2 * vbeta)

    ; Pre-calculate constants outside the loop
    external_M = DBLARR(N_ELEMENTS(sigma_space))
    FOR j = 0, N_ELEMENTS(sigma_space) - 1 DO BEGIN
        external_M[j] = integrated_M_sigma(sigma_space[j], vsini, u1, u2, zeta, t_space)
    ENDFOR

    external_exp = EXP(-2 * !PI^2 * vbeta^2 * sigma_space^2 - 4 * !PI * vgamma * sigma_space)
	n = N_ELEMENTS(f)
    ratios = DBLARR(n)
	; Theta = DBLARR(n)
    FOR i = 0, n - 1 DO BEGIN
        ; Apply TSUM for integration where applicable
        numerators = external_exp * external_M * Theta_sigma(sigma_space, zeta, costheta[i], sintheta[i]) * SIN(2 * !PI * vp[i] * sigma_space) * sigma_space
        denominators = external_exp * external_M * (external_M - f[i] * Theta_sigma(sigma_space, zeta, costheta[i], sintheta[i]) * COS(2 * !PI * vp[i] * sigma_space)) * sigma_space^2
        ; Use TSUM for numerical integration
        numerator = TSUM(sigma_space, numerators)
        denominator = TSUM(sigma_space, denominators)
        ratios[i] = numerator / denominator
    ENDFOR
    RETURN, ratios
END





FUNCTION exofast_rossiter, time, inc_rad, ar, TPeriastron, period, e, omega_rad, p, u1, u2, lambda, vsini, $
         vbeta, vgamma, zeta, $ 
         exptime=exptime, ninterp=ninterp, rmmodel=rmmodel
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
    ; IF (N_ELEMENTS(ohta2005) EQ 0) AND (N_ELEMENTS(hirano2010) EQ 0) AND (N_ELEMENTS(hirano2011) EQ 0) THEN hirano2011 = 1
    IF rmmodel EQ 'ohta2005' THEN ohta2005 = 1
    IF rmmodel EQ 'hirano2010' THEN hirano2010 = 1
    IF rmmodel EQ 'hirano2011' THEN hirano2011 = 1

    vsini = vsini / 1000
    vbeta = vbeta / 1000
    vgamma = vgamma / 1000
    zeta = zeta / 1000

    ; Get the number of points
    npoints = N_ELEMENTS(time)
    IF N_ELEMENTS(ninterp) EQ 0 THEN ninterp = 1

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
    ; Calculate tp using exofast_getphase function
    ; tp = tc - period * exofast_getphase(e, omega_rad, /pri)
    
    ; Get the model flux using exofast_tran function
    tmpmodelflux = exofast_tran(time, inc_rad, ar, TPeriastron, period, e, omega_rad, p, u1, u2, 1)

    ; Calculate phase and true anomaly
    ;phase = exofast_getphase(e, omega_rad, /pri)
    ; Tperiastron = tc - period * phase
    meananom = 2.D0 * !dpi * (1.D0 + (transitbjd - Tperiastron) / Period MOD 1)
    eccanom = exofast_keplereq(meananom, e)
    trueanom = 2.D0 * ATAN(SQRT((1.D0 + e) / (1.D0 - e)) * TAN(eccanom / 2.D0))
    
    ; Calculate coordinates as seen from observer
    ;r = ar*(1d0-e^2)/(1d0+e*cos(trueanom))
    ;print,true
    ;x = - r* COS(trueanom + omega_rad)
    ;y = - r* SIN(trueanom + omega_rad) * COS(inc_rad)

    cosi = COS(inc_rad)
    ; Calculate up, and vp values
    z = exofast_getb2(time, i=inc_rad, a=ar, tperiastron=TPeriastron, period=period, e=e,omega=omega_rad,z2=depth, x2=xp,y2=yp)
    up = xp*cos(lambda) - yp*sin(lambda)
    vp = vsini * up

    ; Compute flux and ratios
    costheta = COS(trueanom)
    sintheta = SIN(trueanom)
    f = 1 - tmpmodelflux
    ratios = DBLARR(npoints)



    IF KEYWORD_SET(ohta2005) THEN BEGIN
        ; PRINT, 'Using ohta2005 model'
		; Ohta+2005 (OTS)
		delta_rv = - f / (1 - f) * vp * 1E3
    ENDIF ELSE IF KEYWORD_SET(hirano2010) THEN BEGIN
        ; PRINT, 'Using hirano2010 model'
        ;Hirano+2010
        sigma = vsini/1.31 ; Gaussian broadening kernel from Hirano+10. Hirano+10 found vsini/1.31 as an approximation that sometimes works
        vbeta_p_2 = vbeta^2
        vbeta_star_2 = vbeta^2 + sigma^2

        ; Calculate terms
        term1 = (2 * vbeta_star_2) / (vbeta_p_2 + vbeta_star_2)
        term2 = f * vp
        term3 = (1 - (vp^2) / (vbeta_p_2 + vbeta_star_2) + (vp^4) / (2 * (vbeta_p_2 + vbeta_star_2)^2))

        ; Calculate delta_v
        delta_rv = - (term1^(3/2)) * term2 * term3 * 1E3
    ENDIF ELSE IF KEYWORD_SET(hirano2011) THEN BEGIN
        ; PRINT, 'Using hirano2011 model'
        ; Call compute_integral to calculate ratios Hirano+ 2011
        ratios = compute_integral( costheta, sintheta, f, vp, u1, u2, vsini, vbeta, vgamma, zeta)
        ; Calculate delta_rv
        delta_rv = -f * ratios / (2 * !pi) * 1E3
    ENDIF ELSE BEGIN
        PRINT, 'Error: No model specified'
    ENDELSE
  
    RETURN, delta_rv
END


