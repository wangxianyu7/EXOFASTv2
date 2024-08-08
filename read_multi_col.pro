FUNCTION read_multi_col, filename
    ; Open the file for reading
    openr, unit, filename, /get_lun

    ; Initialize an empty list to hold the data
    data_list = list()

    ; Read the file line by line
    while not EOF(unit) do begin
        ; Read a single line from the file
        line = ''
        readf, unit, line

        ; Split the line into components
        values = strsplit(line, ' ', /extract)
        
        ; Convert the components to numerical values
        if n_elements(values) GT 0 then begin
            ; Convert each value to double and add to a temporary array
            temp_values = dblarr(n_elements(values))
            for i = 0, n_elements(values)-1 do begin
                temp_values[i] = double(values[i])
            endfor
            ; Add the temporary array to the list
            data_list.add, temp_values
        endif
    endwhile

    ; Convert the list to a 2D array
    if data_list.count() EQ 0 then return, 0
    ncols = n_elements(data_list[0])
    nrows = data_list.count()
    result = dblarr(ncols, nrows)
    for i = 0, nrows-1 do begin
        result[*, i] = data_list[i]
    endfor

    ; Close the file
    free_lun, unit

    ; Return the result
    return, result
END
