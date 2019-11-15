;+
; NAME:
;      HLINE
;     
; PURPOSE:
;      Draw a horizontal line on a pre-existing plot window.
;
; CALLING SEQUENCE:
;      HLINE, VAL
;
; INPUTS:
;
;      VAL: The y-value where the vertical line should be drawn
;
; KEYWORD PARAMETERS:
;
;      All keyword parameters are passed to OPLOT.
;
; SIDE EFFECTS:
;
;      Causes a horizontal line to appear on your screen.
;
; RESTRICTIONS:
;
;      This program won't do anything else. Sorry, them's the 
;      restrictions.
;
; EXAMPLE:
;
;      Draw a horizontal line at x = 4
;      IDL> plot, findgen(10)
;      IDL> hline, 4
;
; MODIFICATION HISTORY:
; Written sometime in 2003 by JohnJohn
;-

pro hline, val,_extra=extra, xlog=xlog, min=min, max=max
nv = n_elements(val)
if 1-keyword_set(min) then if keyword_set(xlog) then min = 1d-20 else $
  min = !x.crange[0]
if 1-keyword_set(max) then if keyword_set(xlog) then max = 1d20 else $
  max = !x.crange[1]
for i = 0, nv-1 do oplot,[min,max],fltarr(2)+val[i],_extra=extra 
end
