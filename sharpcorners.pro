pro sharpcorners, _REF_EXTRA=_extra
;+
; NAME:
;       SHARPCORNERS
;
; PURPOSE:
;       If you look closely, especially if you make your axes thick,
;       the corners of your plots are not sharp.  (Try example below.)
;       This procedure makes sharp corners on your plot.
;
; CALLING SEQUENCE:
;       SHARPCORNERS
;
; Graphics Keywords: 
;       [, COLOR=value] [, LINESTYLE={0|1|2|3|4|5}] [, /NOCLIP] [, 
;       THICK=value]
;
; INPUTS:
;       None.
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       Plot corners are now sharp!
;
; EXAMPLE:
;       IDL> plot, findgen(5), XTHICK=3, YTHICK=3
;       IDL> sharpcorners, THICK=3
;
; MODIFICATION HISTORY:
;   18 Oct 2002  Written by Tim Robishaw, Berkeley
;-

corners = [0,0,0,1,1,1,1,0,0,0,0,1]

for i = 0, 4 do $
  plots, [!x.window[corners[2*i  ]], !x.window[corners[2*(i+1)  ]]], $
         [!y.window[corners[2*i+1]], !y.window[corners[2*(i+1)+1]]], $
         /NORMAL, _EXTRA=_extra

end; sharpcorners
