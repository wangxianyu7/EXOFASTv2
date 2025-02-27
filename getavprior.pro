;+
; NAME:
;   GETAVPRIOR
;
; PURPOSE:
;
;   Returns the maximum V-band extinction (Av) from Schlafly and
;   Finkbeiner (2011) to be used as an upper bound on EXOFASTv2 SED
;   fits.
;
;   If you use this code in your research, please cite:
;   Schlafly and Finkbeiner (2011)
;   http://adsabs.harvard.edu/abs/2011ApJ...737..103S
; 
; CALLING SEQUENCE:
;   maxav = getavprior(ra=ra,dec=dec,object=object,line=line)
;
; OPTIONAL INPUTS:
; 
;   RA     - The J2000 RA of the object in decimal degrees, from which to get the Av.
;   DEC    - The J2000 Dec of the object in decimal degrees, from which to get the Av.
;   OBJECT - The name of the object to be resolved by SIMBAD. Either
;            OBJECT or RA and DEC must be specified
; OPTIONAL OUTPUTS:
;   LINE   - The line in the prior file to impose an upper bound at
;            the maxav from the Schlegel maps.
; 
;  RESULT 
;   The maximum V-band extinction (Av)
;
; MODIFICATION HISTORY
; 
;  2018/04 -- Written, Jason Eastman (CfA)
;-
function getavprior, ra=ra, dec=dec, object=object, line=line

if keyword_set(object) then begin
   url = 'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr='  + object
endif else if n_elements(ra) eq 1 and n_elements(dec) eq 1 then begin
   url = 'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr='  + strtrim(ra,2) + '+'  + strtrim(dec,2) + '+equ+j2000'
endif else message, 'must specify either RA and DEC or OBJECT'

if float(!version.release) ge 8.5 then begin
   output =  strsplit(string(wget(url,/buffer)),string(10B),/extract)
endif else begin
   cmd = "curl '" + url + "'"
   spawn, cmd, output
endelse

if output[0] eq '' then begin
   message, 'Curl not installed? curl command failed: ' + cmd,/continue
   line = 'av 0 -1 0 99'
   maxav = 99d0
   return, 99d0
endif

if strpos(output[2],'Invalid object name') ne -1 then begin
   message, 'Invalid object name (try using RA and Dec instead)',/continue
   line = 'av 0 -1 0 99'
   maxav = 99d0
   return, maxav
endif

match = where(strtrim(output,2) eq '</maxValueSandF>')
maxav = strtrim((strsplit(output[match-1]*3.1d0*1.5d0*0.87d0,'(',/extract))[0],2)
line = 'av 0 -1 0 ' + maxav

return, double(maxav)

end
