pro getticsed
  ; Get command-line arguments
  args = command_line_args()
  ; idl -e getticsed -args "86396382"
  ; Ensure at least one argument is passed
  if n_elements(args) lt 1 then begin
    print, 'Error: No input argument provided. Usage: IDL> getticsed <input_file>'
    return
  endif

  ; Define the filenames based on the input argument
  input_file = args[0]
  prior_file = input_file + '.priors'
  sed_file = input_file + '.sed'

  ; Call the mkticsed procedure with appropriate arguments
  mkticsed, input_file, priorfile=prior_file, sedfile=sed_file

  ; Print messages for confirmation
  print, 'Called mkticsed with:'
  print, '  Input File: ', input_file
  print, '  Prior File: ', prior_file
  print, '  SED File: ', sed_file
end

