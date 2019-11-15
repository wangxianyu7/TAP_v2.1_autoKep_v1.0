function ekepler,m,e
  ;; real*8 e,e0,eps,m,ms,pi2,f0,f1,f2,f3,d1,d2,d3
  ;; c This routine solves Kepler's equation for E as a function of (e,M)
  ;; c using the procedure outlined in Murray & Dermott:
  eps=1.d-10
  pi2=2.*acos(-1.d0)
  ms=(m mod pi2)
  d3=1.d10
  e0=ms+e*0.85d0*sin(ms)/abs(sin(ms))
  while(max(abs(d3)) gt eps) do  begin
     f3=e*cos(e0)
     f2=e*sin(e0)
     f1=1.-f3
     f0=e0-ms-f2
     d1=-f0/f1
     d2=-f0/(f1+0.5*d1*f2)
     d3=-f0/(f1+d2*0.5*(f2+d2*f3/3.))
     e0=e0+d3
  endwhile
  ekep=e0+m-ms
  return,ekep
end

pro kepler_ma,m,e,f
  nm=n_elements(m)
  ekep=dblarr(nm)
  i=long64(0)
  if(e ne 0.d0) then begin
     ekep=ekepler(m,e)
     f=2.d0*atan(sqrt((1.d0+e)/(1.d0-e))*tan(0.5d0*ekep))
  endif else begin
     f=m
  endelse
  nm0=where(m eq 0.d0)
  if(total(nm0) ge 0) then f(nm0)=0.d0
  
  nm0=0L
  ekep = 0L
  i = 0L
end

