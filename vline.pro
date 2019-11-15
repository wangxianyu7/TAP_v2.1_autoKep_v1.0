;+
; NAME:
;      VLINE
;     
; PURPOSE:
;      Draw a vertical line on a pre-existing plot window.
;
; CALLING SEQUENCE:
;      VLINE, VAL
;
; INPUTS:
;
;      VAL: The x-value or array of x-values where the vertical
;      line(s) should be drawn
;
; KEYWORD PARAMETERS:
;
;      All keyword parameters are passed to OPLOT.
;
; SIDE EFFECTS:
;
;      Causes a vertical line to appear on your screen.
;
; RESTRICTIONS:
;
;      This program won't do anything else. Sorry, them's the 
;      restrictions.
;
; EXAMPLE:
;
;      Draw a vertical line at x = 0
;      IDL> plot, findgen(10)
;      IDL> vline, 5
;
; MODIFICATION HISTORY:
; Written sometime in 2003 by JohnJohn
;-

pro vline, val,_extra=extra,ylog=ylog, min=min, max=max
nv = n_elements(val)
if 1-keyword_set(min) then if keyword_set(ylog) then min = 1d-20 else $
  min = !y.crange[0]
if 1-keyword_set(max) then if keyword_set(ylog) then max = 1d20 else $
  max = !y.crange[1]

for i = 0,nv-1 do oplot,fltarr(2)+val[i],[min,max],_extra=extra 
end
