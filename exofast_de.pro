
; Mutation function with bounds
function mutation, XTemp, F, xMin, xMax
    m = n_elements(XTemp[*, 0]); checked, n par
    n = n_elements(XTemp[0, *]); checked, npop
    XMutationTmp = dblarr(m, n); checked
    seed = !NULL
    for i = 0, n - 1 do begin
        repeat begin
            r1 = fix(randomu(seed) * n)
            r2 = fix(randomu(seed) * n)
            r3 = fix(randomu(seed) * n)
        endrep until r1 ne i and r2 ne i and r3 ne i and r1 ne r2 and r1 ne r3 and r2 ne r3
        
        for j = 0, m - 1 do begin
            XMutationTmp[j,i] = XTemp[j, r1] + F * (XTemp[j, r2] - XTemp[j, r3])
			;print,XMutationTmp[j,i]
            ; Ensure the mutated value is within bounds
            if XMutationTmp[j,i] lt xMin[j] then XMutationTmp[j,i] = xMin[j]
            if XMutationTmp[j,i] gt xMax[j] then XMutationTmp[j,i] = xMax[j]
        endfor
    endfor
    return, XMutationTmp
end

; Crossover function with bounds
function crossover, XTemp, XMutationTmp, CR, xMin, xMax
    m = n_elements(XTemp[*, 0])
    n = n_elements(XTemp[0, *])
    XCorssOverTmp = dblarr(m, n)
    seed = !NULL
    for i = 0, n - 1 do begin
        for j = 0, m - 1 do begin
            r = randomu(seed)
            if r le CR then XCorssOverTmp[j, i] = XMutationTmp[j, i] $
            else XCorssOverTmp[j, i] = XTemp[j, i]
            ; Ensure the crossover value is within bounds
            if XCorssOverTmp[j,i] lt xMin[j] then XCorssOverTmp[j,i] = xMin[j]
            if XCorssOverTmp[j,i] gt xMax[j] then XCorssOverTmp[j,i] = xMax[j]
        endfor
    endfor
    return, XCorssOverTmp
end


; Selection procedure
pro selection, XTemp, XCorssOverTmp, fitnessVal, fitness_function
    m = n_elements(XTemp[*, 0])
    n = n_elements(XTemp[0, *])
    fitnessCrossOverVal = dblarr(n, 1)
    for i = 0, n - 1 do begin
        fitnessCrossOverVal[i, 0] = call_function(fitness_function, XCorssOverTmp[*, i])
        if fitnessCrossOverVal[i, 0] lt fitnessVal[i, 0] then begin
            XTemp[*, i] = XCorssOverTmp[*, i]
            fitnessVal[i, 0] = fitnessCrossOverVal[i, 0]
			
        endif
		;print,fitnessCrossOverVal[i, 0], fitnessVal[i, 0],'ssss'
    endfor
end

pro saveBest, fitnessVal
    m = n_elements(fitnessVal)
    tmp = 0
    for i = 1, m - 1 do begin
        if fitnessVal[tmp, 0] gt fitnessVal[i, 0] then tmp = i
    endfor
    ;print, 'Best fitness: ', fitnessVal[tmp, 0]
end


function exofast_de, npop, p0, scale, ngen, fitness_function, min_ptp
    ; reference: https://github.com/zhaozhiyong1989/Differential-Evolution/blob/master/DE.py
    if n_elements(fitness_function) eq 0 then begin
        print, 'A fitness function must be provided'
        return, -1
    endif
    ; Check bounds dimensions

    nrows = n_elements(p0)  ; Get the number of rows
	dim = n_elements(p0)
    ncols = 2
	bounds = dblarr(2, dim)
	for i = 0, dim - 1 do begin
        bounds[0, i] = p0[i] - scale[i]
		bounds[1, i] = p0[i] + scale[i]
	endfor


    xMin = bounds[0, *]
    xMax = bounds[1, *]

    ; Initialize population
    XTemp = dblarr(dim, npop)
    seed = 0
    for i = 0, npop - 1 do begin
        for j = 0, dim - 1 do begin
            XTemp[j,i] = xMin[j] + randomu(seed) * (xMax[j] - xMin[j])
        endfor
    endfor

    ; insert the initial guess
    for j = 0, dim - 1 do begin
        XTemp[j,0] = p0[j]
    endfor

   ; Calculate initial fitness
    fitnessVal = dblarr(npop,1)
    for i = 0, npop - 1 do begin
        fitnessVal[i, 0] = call_function(fitness_function, XTemp[*, i])
    endfor

    ; Evolution process
    for gen = 0, ngen do begin
        F = randomu(!NULL)*0.5 + 0.25
		CR = randomu(!NULL)*0.75 + 0.25
        XMutationTmp = mutation(XTemp, F, xMin, xMax)
        XCorssOverTmp = crossover(XTemp, XMutationTmp, CR, xMin, xMax)
        selection, XTemp, XCorssOverTmp, fitnessVal, fitness_function
        saveBest, fitnessVal
		print, 'DE: ngen', gen, ', delta chi2:',abs(max(fitnessVal) - min(fitnessVal)),'<',min_ptp, ', max chi2:', max(fitnessVal), ', min chi2:', min(fitnessVal)
        if abs(max(fitnessVal) - min(fitnessVal)) lt min_ptp then break
    endfor
	min_value = MIN(fitnessVal, index)
    return, XTemp[*,index]

END



