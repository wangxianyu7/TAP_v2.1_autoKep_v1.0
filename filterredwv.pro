function filterredwv,x,sigma_r,sigma_w,zeropad=zeropad

;+
; NAME:
;       FILTERREDWV
; PURPOSE:
;       Returns the 1/f component from a data vector x that is assumed
;       to be described for parameters sigma_r and sigma_w as 1/f
;       noise plus white noise as detailed by Carter & Winn (2009)
;
; CALLING SEQUENCE:
;       result = filterredwv(x, sigma_r, sigma_w)
;
; INPUTS:
;       X - Data vector, must be a power of two
;       SIGMA_R - Red noise amplitude parameter (not equal to
;                 1/f-component RMS)
;       SIGMA_W - White noise amplitude parameter (approximately equal
;                 to white-component RMS)
;
; OUTPUTS:
;       RESULT - 1/f component with same number of element as for input X.
;
; NOTES:
;       Refer to the paper by Carter & Winn (2009) for theory and
;       details.  In the notation of that paper, gamma=1 for this
;       algorithm.
;
; REVISION HISTORY:
;       Written,    J. Carter               September, 2009
;-   

gamma=1
d=x
 els = n_elements(d)
  pow = ceil(alog(els)/alog(2.))
  if (2^pow ne els and keyword_set(zeropad)) then begin
;;     print, 'inside if'
     diff = 2.^pow-els
     left = floor(diff/2.)
     right = diff-left
     if diff gt 1 then x = [replicate(0d,left),d,replicate(0d,right)] else $
        x = [d,replicate(0d,right)] 
  endif else x = d
  
J = alog(double(n_elements(x)))/alog(2d0)
;;help, J
if (abs(J-fix(J)) ne 0) then message,'Data length must be a power of two'
J = fix(J)

info = WV_FN_DAUBECHIES(2,wavelet,scaling,ioff,joff)
wv = wv_dwt(x,wavelet,scaling,ioff,joff)

sm2 = sigma_r^2*(gamma eq 1 ? 1.0/(2.0*alog(2.0)) : 2-2^gamma)+sigma_w^2

wv[0] *= 1-(sigma_w^2)/sm2

k = 1.
for i=1,J,1 do begin
   sm2 = sigma_r^2*2^(-gamma*i*(1.0d))+sigma_w^2
   for m=0.,2^(i-1)-1,1 do begin
      wv[k] *= 1-sigma_w^2/sm2
      k+=1
   endfor
endfor

redcomp = wv_dwt(wv,wavelet,scaling,ioff,joff,/inverse)

if keyword_set(zeropad) then begin
   if n_elements(left) ne 0 then begin
      if left ne 0 and right ne 0 then begin
         x = d
                                ;whitecomp = whitecomp(left:(els+left))
         redcomp = redcomp(left:(els+left))
      endif
   endif
endif
d=0L
return,redcomp

end
