FUNCTION computeGP, N, x, y, diag_, kernel_type, output_type, par1, par2, par3, par4, par5, $
    AUTO_GLUE=auto_glue, DEBUG=debug, VERBOSE=verbose
    ;print,'We are in the computeGP function'
    ; Enable debugging with ON_ERROR if DEBUG is set
    IF NOT(KEYWORD_SET(debug)) THEN ON_ERROR, 2

    ; Set default values if parameters are not defined
    N = LONG(N)
    par1 = DOUBLE(par1)
    par2 = DOUBLE(par2)
    par3 = DOUBLE(par3)
    par4 = DOUBLE(par4)
    par5 = DOUBLE(par5)
    kernel_type= LONG(kernel_type)
    ; 0: Real, 1: Complex, 2: SHO, 3: Matern32, 4: Matern52, 5: RotationTerm
    output_type= (SIZE(output_type, /TNAME) EQ 'UNDEFINED') ? 0L : LONG(output_type)
    ; 0: log likelihood, 1: prediction

    ; Generate default data for x, y, and yerr if not provided
    x = DOUBLE(x)
    y = DOUBLE(y)
    diag_ = DOUBLE(diag_)
    
    ; Check that x, y, and diag_ have the correct length
    IF (N NE SIZE(x, /N_ELEMENTS) OR N NE SIZE(y, /N_ELEMENTS) OR N NE SIZE(diag_, /N_ELEMENTS)) THEN $
        MESSAGE, 'x, y, and diag_ must all have length N'

    ; Prepare output array for the result
    result = DBLARR(N)

    ; Call the external function
    func = KEYWORD_SET(auto_glue) ? 'compute_GP_natural' : 'computeGP_wrapper'
    _ = CALL_EXTERNAL('/N/project/spinOrbit_Angles/codes/nv5/idl/lib/EXOFASTv2/celerite2IDL/computeGP.so', func, $
    N, x, y, diag_, result, kernel_type, output_type, par1, par2, par3, par4, par5, $
    VALUE=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], $
    /D_VALUE, /CDECL, AUTO_GLUE=auto_glue, VERBOSE=verbose, SHOW_ALL_OUTPUT=verbose)

    ; If output_type is 0 (log likelihood), return only the first element
    IF output_type EQ 0 THEN BEGIN
        RETURN, result[0]
    ENDIF ELSE BEGIN
        RETURN, result
    ENDELSE

END

