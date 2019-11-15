function TAP_MPFit_function,param,time=time,dfit=dfit,finefit=finefit,flux=flux $
                            ,redfit=redfit,tdat=tdat,fdat=fdat,efit=efit,rebin=rebin,smooth=smooth

  if n_elements(flux) eq 0 then flux = time*0+1d0
  inc = param[1]
  if inc gt 90d0 then inc = 90d0 - (inc mod 90d0)
   
  tap_transite,time,[param[0],inc,1d0/param[2],param[3:n_elements(param)-1]],finefit;,/trend
  
  if rebin eq 0 then dfit = finefit
  if rebin eq 1 then dfit = ((dblarr(smooth[1])+smooth[1]^(-1d0))#reform(finefit,smooth[1],n_elements(tdat)))[0,*]
  dfit *= poly(tdat-min(tdat),param[9:10])
  
  if keyword_set(redfit) then begin
     redfit = (filterredwv((fdat-dfit),param[11],param[12],/zeropad))[0:n_elements(fdat)-1]
  endif
  inc = 0L
  if keyword_set(redfit) then return,(fdat-dfit-redfit) else return,(fdat-dfit)
end

function tap_diff,ptr
  return,max(*ptr)-min(*ptr)
end

function tap_mm,ptr,adjust=adjust,shift=shift
  if 1-keyword_set(shift) then shift = 0d0
  if 1-keyword_set(adjust) then adjust = 0d0
  return,[min(*ptr)-shift-tap_diff(ptr)*adjust,max(*ptr)-shift+tap_diff(ptr)*adjust]
end

function tap_colors
  colors= {white:  -1 ,$
           red:     0 ,$
           blue:    1 ,$  
           black:   2 ,$
           green:   4 ,$
           orange:  5 ,$
           gray:    6 ,$
           violet:  7 ,$
           dkblue:  8 ,$
           dkgreen: 9 ,$
           dkorange:10,$
           red2:    11,$
           red3:    12,$
           green1:  13,$
           green2:  14,$
           green3:  15 $
          }
  tvlct,0,0,0,colors.black       ; black
  tvlct,255,255,255,colors.white ; white
  tvlct,255,0,0,colors.red       ; red
  tvlct,155,0,0,colors.red2      ; red
  tvlct,055,0,0,colors.red3      ; red
  tvlct,0,180,100,colors.green   ; green
  tvlct,0,220,100,colors.dkgreen ; green
  tvlct,0,110,220,colors.blue    ; blue
  tvlct,0,110,420,colors.dkblue  ; blue
  tvlct,250,112,0,colors.orange  ; orange
  tvlct,255,220,0,colors.dkorange; orange
  tvlct,120,120,120,colors.gray
  tvlct,131,61,172,colors.violet ; indigo_violet_1
  tvlct,0,220,80,colors.green1   ; green
      tvlct,0,250,60,colors.green2   ; green
  tvlct,0,160,40,colors.green3   ; green

return,colors

end


pro tap::execute_mcmc,restart=restart
  self.message = 'Running Markov Chain Monte Carlo'
  if keyword_set(restart) then self.message = 'Resuming Markov Chain Monte Carlo'
  self->message
  
  transit = (*self.transits)[0]->get()
  transit.diff = self.diff
  (*self.transits)[0]->set,transit
  transit = 0L
  
  if keyword_set(restart) then $
     x = tapmcmc(base_parent=self,input_transits=self.transits,restart=restart) else $
        x = tapmcmc(base_parent=self,input_transits=self.transits)
  
                                ;x = obj_new('TAPmcmc',base_parent=self, input_transits=self.transits)
  x->start
  
  self.message = 'MCMC Complete.'
  self->message
  
  ptr_free,self.mcmc_stuff
  self.mcmc_stuff = ptr_new(x->info())
  cd,current=thisdir
  (*self.mcmc_stuff).savefile =  thisdir+'/'+ (*self.mcmc_stuff).savefile
  widget_control,self.mcmc_fld[1],set_value = (*self.mcmc_stuff).savefile
  widget_control,self.buttons[0],sensitive=1
  x->destroy
  self->loadmcmc
end

pro TAP_Event,event
  Widget_Control, event.id, Get_UValue=info
  Call_Method, info.method, info.object, event
end


function tap::userquery,message
 D = dialog_message(message,/question,/center)
 if strcmp(d,'Yes') then return,1
 return,0
end


function TAP::stripspecial,string
  return,strjoin(strsplit(string,"_",/extract),"\_")
end

pro tap::remake_t,noreset=noreset
  transit = (*self.transits)[self.active_transit-1]->get() 
  case transit.rebin of
     0: t = ptr_new(transit.transit_tday)
     1: begin 
        if transit.t_int[0] eq 0 then begin
           transit.t_int[0] = median(transit.transit_tday[1:n_elements(transit.transit_tday)-1]-transit.transit_tday[0:n_elements(transit.transit_tday)-2]) 
           transit.t_int[1] = 1
        endif
        
        t1 = dblarr(transit.t_int[1])+1d0
        t2 =(findgen(transit.t_int[1])+1-5d-1*(transit.t_int[1]+1d0))*(transit.t_int[0]/transit.t_int[1])/1440d0
        t3 = (dblarr(n_elements(transit.transit_tday))+1)
   
        t = ptr_new((reform((t1)#transit.transit_tday +t2#t3,transit.t_int[1]*n_elements(transit.transit_tday),1))[*,0]) 
        t1 = 0
        t2 = 0
        t3 = 0
     end
  endcase
  ptr_free,transit.model_tfine,transit.model_ffine
  transit.model_tfine = ptr_new(*t)
  transit.model_ffine = ptr_new((*t)*0d0)
  ptr_free,t
  (*self.transits)[self.active_transit-1]->set,transit
  transit = 0L
  if 1-keyword_set(noreset) then begin
     widget_control,self.settings[2],set_value=((*self.transits)[self.active_transit-1]->get()).t_int[0]
     widget_control,self.settings[3],set_value=round(((*self.transits)[self.active_transit-1]->get()).t_int[1])
  endif
end

pro TAP::ButtonEvent,event
  widget_control, event.id, GET_UVALUE= uvalue
  widget_control, /Hourglass
  
  case uvalue.value of
     'create_plots': if self.create_plots then self.create_plots = 0 else self.create_plots = 1 
     'create_mcmcascii':  if self.create_mcmcascii then self.create_mcmcascii = 0 else self.create_mcmcascii = 1
     'diffset': begin
        self.diff = event.value
        self->lcplot
     end
     'Save Current Setup Button': begin
        path = dialog_pickfile(dialog_parent=(*self.bases)[0],title='Select / Type save file')
        if path ne '' then begin
           widget_control,self.transitfile2_fld[1],set_value = path
        endif 
     end
     'Delete Active':begin        
        transit = (*self.transits)[self.active_transit-1]->get()
        for k=0,n_elements(transit.params)-1,1 do begin
           ptr_free,transit.params[k].mcmc_chain
           ptr_free,transit.params[k].refined_mcmc_chain
        endfor
        ptr_free,transit.model_t,transit.model_f,transit.model_tfine,transit.model_ffine,transit.mcmc_files
        transit = 0L
        (*self.transits)[self.active_transit-1]->destroy

        if self.num_transits eq 1 then (*self.transits)=0L else $
           if self.active_transit eq self.num_transits then (*self.transits)=(*self.transits)[0:self.active_transit-2] else $
              if self.active_transit eq 1 then (*self.transits)=(*self.transits)[self.active_transit:self.num_transits-1] else $
                 (*self.transits)=[(*self.transits)[0:self.active_transit-2],(*self.transits)[self.active_transit:self.num_transits-1]]
        
        self.num_transits-=1
        self.active_transit = max([1,self.active_transit-1])  
        
        ;widget_control,(*self.bases)[2],MAP=0
        self->setup,'multi'
        widget_control,(*self.bases)[2],MAP=1
        self->setup,'cross lc links'
        self->setup,'fit'
        
        self->lcplot
        self.message = 'light curve deleted.'
        self->message
     end
     'Saved Setup File Button': begin
        path = dialog_pickfile(dialog_parent=(*self.bases)[0],title='Select Save File',filter='TAP_setup.idlsav')
        if path ne '' then begin
           widget_control,self.transitfile1_fld[1],set_value = path
           widget_control,self.buttons[1],sensitive=1
        endif 
     end
     'Load Setup Button': begin
        self.message = 'Loading existing setup...'
        self->message
        widget_control,self.transitfile1_fld[1],get_value=path
        restore,path
        ptr_free,self.transits
        self.transits = ptr_new(TAP_state)
        TAP_state = 0L

        self.diff = ((*self.transits)[0]->get()).diff
  

        self.active_transit =  n_elements(*self.transits) 

        self.num_transits = n_elements(*self.transits)
        self->setup,'fit'      

        self->setup,'multi'      
        self->setup,'adjustparams'
        self->setup,'cross lc links'
        self->lcplot
        self.numcol=2
        self->prepcol
        resume = 0
        if strcmp('-1',(*((*self.transits)[0]->get()).mcmc_files)[0]) eq 0 then $
           if n_elements(*((*self.transits)[0]->get()).mcmc_files) lt (((*self.transits)[0]->get()).mcmc_params)[0] then $
              resume = self->userquery('TAP has detected an imcomplete MCMC execution.  Continue it?')
        self.message = 'Load complete.'
        self->message
        if resume then self->execute_mcmc,restart='/'+strjoin((strsplit(path,'/',/extract))[0:n_elements(strsplit(path,'/',/extract))-2],'/')+'/'  else begin
           transit = (*self.transits)[0]->get()
           *transit.mcmc_files = '-1'
           (*self.transits)[0]->set,transit
           transit = 0L
        endelse
     end
     'Transit File Button': begin
        path = dialog_pickfile(dialog_parent=(*self.bases)[0],title='Select Transit File',/must_exist)
        if path ne '' then begin
           widget_control,self.transitfile_fld[1],set_value = path
           widget_control,self.buttons[2],sensitive=1
        endif
     end
     'FileType': self.filetype = event.value
     'RebinType': begin
        widget_control,self.setup_smooth,sensitive=0
        case event.value of
           'None': self.setup_rebin = 0
          ; 'Rebin to data cadence': self.setup_rebin = 1
           'Rebin to "Input Integration"': begin
              self.setup_rebin = 1
              widget_control,self.setup_smooth,sensitive=1
           end
        endcase
     end
     'RebinTypeA': begin
        transit = (*self.transits)[self.active_transit-1]->get() 
        case event.value of
           'None': begin
              transit.rebin = 0
              widget_control,self.settings[2],sensitive=0
              widget_control,self.settings[3],sensitive=0
           end
           'Rebin to "Input Integration"': begin
              transit.rebin = 1
              widget_control,self.settings[2],sensitive=1
              widget_control,self.settings[3],sensitive=1
;              widget_control,self.setup_smooth,sensitive=1
           end
        endcase
        (*self.transits)[self.active_transit-1]->set,transit
        transit = 0L
        self->remake_t
        self->updatemodel
        self->lcplot
        
     end
     'settint': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        self.smooth_val[uvalue.set] = *event.value
;event.value
     end
     'settint2': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 1
        transit = (*self.transits)[self.active_transit-1]->get() 
        transit.t_int[uvalue.set] = *event.value
        (*self.transits)[self.active_transit-1]->set,transit
        transit = 0L
        self->remake_t,/noreset
        self->updatemodel
        self->lcplot
        
;event.value
     end
     'Load Transit Button': begin
        self->loadtransit,event
        ;for i=1,n_elements(*self.bases)-1,1 do widget_control,(*self.bases)[i],MAP=0
        ;self->setup,'load'
        ;widget_control,(*self.bases)[1],MAP=1
     end
     'Clear Data Path': begin
        widget_control,self.transitfile_fld[1],set_value = ''
        widget_control,self.buttons[2],sensitive=0
     end
     'Clear Setup Path': begin
        widget_control,self.transitfile1_fld[1],set_value = ''
        widget_control,self.buttons[1],sensitive=0
     end
     'ActiveTransit': begin
        self.active_transit = event.value
        widget_control,self.settings[0],set_value=string(self.active_transit,format='(i2.2)')+": "+((*self.transits)[self.active_transit-1]->get()).fname
        widget_control,self.settings[1],set_value=((*self.transits)[self.active_transit-1]->get()).rebin   
        widget_control,self.settings[2],set_value=((*self.transits)[self.active_transit-1]->get()).t_int[0]
        widget_control,self.settings[3],set_value=round(((*self.transits)[self.active_transit-1]->get()).t_int[1])
        if ((*self.transits)[self.active_transit-1]->get()).rebin eq 1 then begin
           widget_control,self.settings[2],sensitive=1
           widget_control,self.settings[3],sensitive=1
        endif else begin
           widget_control,self.settings[2],sensitive=0
           widget_control,self.settings[3],sensitive=0
        endelse
        self->lcplot
        self->setup,'adjustparams'
        self->setup,'cross lc links'
        self->prepcol
     end
     'Setup Cross LC Locks': begin
        centertlb,(*self.extra_windows)[0]
        widget_control,(*self.extra_windows)[0],/realize
     end
     'Manual Parameter Adjustment': begin
        self->setup,'adjustparams'
        centertlb,(*self.extra_windows)[1]
        widget_control,(*self.extra_windows)[1],/realize
     end
     'Adjust Limits and Locks': begin
        self->setup,'adjustlimlocks'
        centertlb,(*self.extra_windows)[2]
        widget_control,(*self.extra_windows)[2],/realize
     end 
     'MCMC Parameters': begin
        self->setup,'mcmcparams'
        centertlb,(*self.extra_windows)[3]
        widget_control,(*self.extra_windows)[3],/realize
     end
     'Gaussian Priors': begin
        self->setup,'gaussianpriors'
        centertlb,(*self.extra_windows)[4]
        widget_control,(*self.extra_windows)[4],/realize
     end
     'Execute MCMC': begin
        self->execute_mcmc
     end
     'MCMC File Button': begin
        path = dialog_pickfile(dialog_parent=self.tap_base,title='Select MCMC SETUP file',/must_exist,filter='TAP_setup.idlsav')
        if path ne '' then begin
           widget_control,self.mcmc_fld[1],set_value = path
           widget_control,self.buttons[0],sensitive=1
        endif
     end
     'Load MCMC Button': begin
        self->loadmcmc
    
     end
     'Clear MCMC Path': begin
        widget_control,self.mcmc_fld[1],set_value=''
        widget_control,self.buttons[0],sensitive=0
                                ; self->setup,'inference'
     end 
     'xlock_freeall': begin
        widget_control, event.id, GET_UVALUE= uvalue
        for i=0,self.num_transits-1,1 do begin
           transit = (*self.transits)[i]->get()
           transit.params[uvalue.param].set = i+1
           widget_control,self.fld[i,uvalue.param],set_value = round(transit.params[uvalue.param].set)
           (*self.transits)[i]->set,transit
           transit = 0L
        endfor
    
        self->updatemodel
     end
     'xlock_lockall': begin
        widget_control, event.id, GET_UVALUE= uvalue
        for i=0,self.num_transits-1,1 do begin
           transit = (*self.transits)[i]->get()
           transit.params[uvalue.param].set = 1
           widget_control,self.fld[i,uvalue.param],set_value = round(transit.params[uvalue.param].set)
           (*self.transits)[i]->set,transit
           transit = 0L
        endfor
    
        self->updatemodel
     end
     else: print,'Unknown Button Event "'+uvalue.value+'"'
  endcase
end

pro TAP::QuitWidget,event
  widget_control, event.id, GET_UVALUE= uvalue

  ;help,/heap
;  stop
  widget_control,(*self.extra_windows)[uvalue.wid],map=0
  if uvalue.wid eq 0 then self.fld[*,*] = 0
;     for i=0,n_elements(self.fld[*,0])-1,1 do $
;        for j=0,n_elements(self.fld[0,*])-1 do if self.fld[i,j] ne 0 then begin
;     widget_control,self.fld[i,j],/destroy
;     self.fld[i,j] = 0
;  endif
  widget_control,(*self.extra_windows)[uvalue.wid],/destroy
  
  case uvalue.wid of
     0: self->setup,'cross lc links'
     1: self->setup,'adjustparams'
     2: self->setup,'adjustlimlocks'
     3: self->setup,'mcmcparams'
     4: self->setup,'gaussianpriors'
     else: print,uvalue.wid+' is an unknown setup'
  endcase 
  uvalue = 0L
end

;; pro TAP::QuitLinks,event
;;   widget_control,(*self.extra_windows)[0],/destroy
;;   self->setup,'cross lc links'
;; end

;; pro TAP::QuitManualAdjust,event
;;   widget_control,(*self.extra_windows)[1],/destroy
;;   self->setup,'adjustparams'
;; end

;; pro TAP::QuitAdjustLL,event
;;   widget_control,(*self.extra_windows)[2],/destroy
;;  ; self->setup,'adjustlimlocks'
;; end

;; pro tap::quitmcmcparams,event
;;   widget_control,(*self.extra_windows)[3],/destroy
;; end

pro TAP::MCMC_inference
  self.message = 'Conducting Bayesian inference.'
  self->message

  for i=0,n_elements(*self.transits)-1,1 do begin
     self.message = 'Transit '+string(i+1,format='(i2.2)')+' of '+string(self.num_transits,format='(i2.2)')
     self->message
     self.active_transit = i+1
     transit = (*self.transits)[i]->get()
     for k=0,n_elements(transit.params)-1,1 do begin
        bigarr = *transit.params[k].refined_mcmc_chain
        if strcmp(transit.params[k].param,'Inclination') then if (where(bigarr ge 90d0))[0] ne -1 then bigarr[where(bigarr ge 90)] = 90d0 - (bigarr[where(bigarr ge 90d0)] mod 90d0)
        
        sorted = bigarr[sort(bigarr)]
        range = n_elements(bigarr)/100d0
        transit.params[k].mcmc_val = [sorted[50d0*range],$
                                      sorted[84.135d0*range]-sorted[50d0*range],$
                                      sorted[50d0*range]-sorted[15.865d0*range]]
        if transit.params[k].fixed then transit.params[k].mcmc_val[1:2] = -1d0
        transit.params[k].value =  transit.params[k].mcmc_val[0]
        
        sorted = 0L
        range = 0L
        bigarr = 0L
     endfor
     (*self.transits)[i]->set,transit
     self.message = 'Combined '+string(self.mcmc_complete,format='(i2.2)')+' MCMC chains.'
     if transit.mcmc_params[0] gt 1 then self->message
     transit=0L
     self->updatemodel
     self.numcol=4
     self->prepcol
     self->lcplot
  endfor
  
  if self.create_ascii then self->create_ascii
  if self.create_plots then self->create_plots
  

  widget_control,self.mcmc_fld[1],get_value=path
  tap_state = (*self.transits)
  for i=0,n_elements(tap_state)-1,1 do begin
     transit = tap_state[i]->get()
     for j=0,n_elements(transit.params)-1,1 do begin
        *transit.params[j].mcmc_chain = -1
        *transit.params[j].refined_mcmc_chain = -1        
     endfor
     tap_state[i]->set,transit     
  endfor
  save,TAP_state,filename=path
;;   for i=0,n_elements(tap_state)-1,1 do begin
;;      st = TAP_state[i]->get()
;;      for j=0,n_elements(st.params)-1,1 do begin
;;         ptr_free,st.params[j].mcmc_chain
;;         ptr_free,st.params[j].refined_mcmc_chain
;;      endfor
;;      ptr_free,st.model_t,st.model_f,st.model_tfine,st.model_ffine,st.mcmc_files
;;      TAP_state[i]->destroy
;;      st = 0L
;;   endfor 
  TAP_state = 0L
  

  self.message = 'MCMC inference complete.'
  self->message
end







pro TAP::plot_open, filename, $
                      LANDSCAPE=landscape, $
                      XSIZE=xsize, $
                      YSIZE=ysize, $
                      INCHES=inches, $
                      COLOR=color, $
                      ENCAPSULATED=encapsulated, $
                      BITS_PER_PIXEL=bits_per_pixel, $
                      _REF_EXTRA=_extra

  self.base_plot = !d.name
  set_plot, 'PS', COPY=keyword_set(COLOR), INTERPOLATE=keyword_set(COLOR)
  
  device, FILENAME  = filename, $
          LANDSCAPE = keyword_set(LANDSCAPE), $
          XSIZE     = xsize, $
          YSIZE     = ysize, $
          XOFFSET   = xoffset, $
          YOFFSET   = yoffset, $
          INCHES    = keyword_set(INCHES), $
          COLOR     = keyword_set(COLOR), $
          BITS_PER_PIXEL = 8, $
          ENCAPSULATED   = keyword_set(ENCAPSULATED), $
          _EXTRA    = _extra
end




pro TAP::create_plots

  self.message = 'Creating MCMC plots...'
  self->message
  
  tempP = !p
  tempX = !x
  tempY = !y

  !p.charsize=1.1
  !p.charthick=4
  !x.thick = 5
  !y.thick = 5
  !p.thick= 5
  !p.font = 0
  
  ;; first plot--ALL LCs together
  tap_readcol,self.plot_dir+"ascii_phased_data.ascii",p,f,rf,t,m,r,format='(d,d,d,d,d,d)'
  
  ysize=4
  self->plot_open,self.plot_dir+'plot_alltransit_phased_lc.eps',xsize=6,ysize=ysize,/inches,/encapsulated,/color
  yr = [min(f)-(15d0*stdev(r)),max(f)+(3d0*stdev(r))]
  plot,p[sort(p)]*24d0,f[sort(p)],psym=8,color=(*self.colors).black,background=(*self.colors).white,yrange=yr, $
       /ys,/xs,xrange=mm(p*24d0),symsize=.5,xthick=4,ythick=4,xtitle='Hours from Mid Transit',ytitle='Relative Flux'
                                ; oplot,p[sort(p)],m[sort(p)],color=(*self.colors).blue,thick=4
  oplot,p[sort(p)]*24d0,r[sort(p)]+(min(f)-(8d0*stdev(r))),color=(*self.colors).black,psym=8,symsize=.5
                                ; hline,(min(f)-(8d0*stdev(r))),color=(*self.colors).blue,thick=4
  
  sharpcorners,thick=4,color=(*self.colors).black
  self->plot_close

  
  tap_readcol,self.plot_dir+"ascii_phased_model.ascii",t2,p2,mf2,format='(d,d,d)'
  model = interpol(mf2[sort(p2)],p2[sort(p2)],p[sort(p)])
  ysize=4
  self->plot_open,self.plot_dir+'plot_alltransit_phase_modandresid.eps',xsize=6,ysize=ysize,/inches,/encapsulated,/color
  yr = [min(f)-(15d0*stdev(r)),max(f)+(3d0*stdev(r))]
  plot,p[sort(p)]*24d0,model+r[sort(p)],psym=8,color=(*self.colors).black,background=(*self.colors).white,yrange=yr, $
       /ys,/xs,xrange=mm(p*24d0),symsize=.5,xthick=4,ythick=4,xtitle='Hours from Mid Transit',ytitle='Relative Flux'
  oplot,p[sort(p)]*24d0,model,color=(*self.colors).blue,thick=5
  oplot,p[sort(p)]*24d0,r[sort(p)]+(min(f)-(8d0*stdev(r))),color=(*self.colors).black,psym=8,symsize=.5
  hline,(min(f)-(8d0*stdev(r))),color=(*self.colors).blue,thick=5
  
  sharpcorners,thick=4,color=(*self.colors).black
  self->plot_close
 
  
  !p.multi=[0,1,1]
  ysize = min([8,3d0 + 1d0*n_elements(*self.transits)])
  self->plot_open,self.plot_dir+'plot_alltransit_lightcurve.eps',xsize=6,ysize=ysize,/inches,/encapsulated,/color
  
  
  for i=0,n_elements(*self.transits)-1,1 do $
     (*self.plot_windows)[0].xrange = [$
     min([(*self.plot_windows)[0].xrange[0],min(((*self.transits)[i]->get()).transit_tday- $
                                                (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit')  $
                                                                                           eq 1)].value-(((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))]),$
                                                                 max([(*self.plot_windows)[0].xrange[1],$
                                                                      max(((*self.transits)[i]->get()).transit_tday-$
                                                                          (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value $
                                                                           - (((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))])]
  



  tot = self.num_transits-1
  lci = (*self.transits)[0]->get()
  lcf = (*self.transits)[self.num_transits-1]->get()
 ; yrange = [.995*min(lci.transit_flux)-tot*.01,1.00*max(lcf.transit_flux)+((tot)*.019)]
  diff = self.diff
  if diff eq 0 then diff = 10d0*max([lci.residuals,lcf.residuals])
  yrange = [min(lci.transit_flux)-tot*(diff/2d0)-diff,max(lcf.transit_flux)+((tot+.5)*diff)]

  for i=0,self.num_transits-1,1 do begin
   lc = (*self.transits)[i]->get() 
  ;   diff = 5d0*max(lc.residuals)
     midt = lc.params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
     if i eq 0 then $
        plot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,yrange=yrange,$
color=(*self.colors).black,background=(*self.colors).white,xrange=(*self.plot_windows)[0].xrange*24d0,$
             /xs,/ys,xtitle='Hours from Mid Transit',ytitle='Relative Flux + Constant',title=title,/nodata
     
     oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,color=(*self.colors).black
     
     oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(diff*i),color=lc.modcol,thick=2
     oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)*(i+1)),psym=8,symsize=.6,color=(*self.colors).black
     hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=5
     oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=5
     oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=5
     xyouts,!x.crange[0]+.05d0*(!x.crange[1]-!x.crange[0]),1d0+diff*i+diff/6,$
            string(i+1,format='(i2.2)')+': '+lc.fname,color=(*self.colors).black,charsize=1,charthick=2,align=0
     lc = 0L   
     midt = 0L
  endfor
  lci = 0L
  lcf = 0L
  tot = 0L
  yrange = 0L
  self->plot_close
  ;spawn,'open plot_alltransit_lightcurve.eps'

  self->plot_open,self.plot_dir+'plot_alltransit_mcmc_system.eps',xsize=7,ysize=8,/inches,/encapsulated,/color
  !y.tickname=' '   
  !p.multi = [0,2,3]
  bin =30
  
  xtitle= ['Days','Degrees','Parameter Value','Parameter Value','Days since Mid Transit']
  order = [0,1,2,3,4]

  for j=0,n_elements(order)-1,1 do begin
     yr = [0,0]
     all_lock = 1
  
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        if transit.params[order[j]].mcmc_val[1] ge 0 then begin
                                ;  xmin = floor(min(*transit.params[where(strcmp(((*self.transits)[0]->get()).params.param,'Mid Transit') eq 1)].refined_mcmc_chain))    
           modify = 0d0
           if strcmp(transit.params[order[j]].param,'Mid Transit') then modify = transit.params[order[j]].mcmc_val[0]
           plothist,(*transit.params[order[j]].refined_mcmc_chain)-modify,tx,ty,$
                    bin=tap_diff(( (*self.transits)[0]->get()).params[order[j]].refined_mcmc_chain)/bin,/noplot
           txr = tap_mm(transit.params[order[j]].refined_mcmc_chain,shift=modify,adjust=.05)
           if all_lock eq 1 then xr = txr else xr = [min([xr[0],txr[0]]),max([xr[1],txr[1]])]
           if strcmp(transit.params[order[j]].param,'Inclination') then xr[1] = min([90d0,xr[1]]) 
           
                                ;  stop
           yr = [0,max([yr[1],max(ty)])]
                                ;  print,yr
           tx = 0L
           ty = 0L
           transit = 0L
           r = 0L
           all_lock = 0
        endif
     endfor
     ;; plot the hist
     first = 1
     yr[1]*= 1.1d0
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        
        if all_lock then begin
           yr = [0,2]
           if first then begin
              !x.tickname=' '
              plot,[0,2],[0,2],/nodata,title=transit.params[order[j]].param,background=(*self.colors).white,$
                   color=(*self.colors).black,thick=4,xthick=4,ythick=4,charthick=4
              xyouts,1,.9,'Fixed',color=(*self.colors).black,charsize=3d,charthick=3d,align=0.5
              !x.tickname=''
              first = 0L
           endif
           if j eq 0 then xyouts,!x.crange[1]-.3d0*(!x.crange[1]-!x.crange[0]),yr[1]-(.05*(i+1)*yr[1]),transit.fname,color=transit.color,charsize=0.8d0
           
        endif else begin
           modify = 0d0
           if strcmp(transit.params[order[j]].param,'Mid Transit') then modify = transit.params[order[j]].mcmc_val[0]
           if transit.params[order[j]].mcmc_val[1] ge 0 then begin
              if first eq 1 then begin
                 plothist,(*transit.params[order[j]].refined_mcmc_chain)-modify,$
                          bin=tap_diff(( (*self.transits)[0]->get()).params[order[j]].refined_mcmc_chain)/bin,$
                          color=transit.color, $
                          background=(*self.colors).white,axiscolor=(*self.colors).black, $
                          xtitle=xtitle[j],title=transit.params[order[j]].param, $
                          xrange=xr,/xs, $
                          thick=2*n_elements(*self.transits),xticks=3,xthick=2,ythick=2,yrange=yr,/ys 
                 first = 0
              endif else $
                 plothist,(*transit.params[order[j]].refined_mcmc_chain)-modify, $
                          bin=tap_diff(( (*self.transits)[0]->get()).params[order[j]].refined_mcmc_chain)/bin,$
                                ;  bin=tap_diff(transit.params[order[j]].refined_mcmc_chain)/bin,$
                          color=transit.color,/overplot, $
                          background=(*self.colors).white,axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                          title=transit.params[order[j]].param,xrange=xr,/xs, $
                          thick=2*n_elements(*self.transits)-1.2*i,xticks=3,xthick=2,ythick=2,yrange=yr,/ys
              if j eq 0 then xyouts,!x.crange[1]-.3d0*(!x.crange[1]-!x.crange[0]),yr[1]-(.05*(i+1)*yr[1]),transit.fname,color=transit.color,charsize=0.8d0
           endif 
        endelse
              transit = 0L
              
     endfor
     
  endfor 
  !y.tickname=''   
  self->plot_close
  





  ;;order = [0,1]
  
  ;;for j=0,n_elements(order)-1,1 do $
  ;;   for k=0,n_elements(*self.transits)-1,1 do begin
  ;;   transit  = (*self.transits)[k]->get() 
  ;;   if 
  ;;  endfor
    

  for i=0,n_elements(*self.transits)-1,1 do begin
     !y.tickname=' '
  !p.charsize = 1.4d0
  
     
     ;; plot system variables
     transit  = (*self.transits)[i]->get() 
     self->plot_open,self.plot_dir+'plot_'+transit.fname+'_mcmc_params1.eps',xsize=7,ysize=8,/inches,/encapsulated,/color
     !p.multi = [0,2,3]
     bin =30


  xmin = floor(min(*transit.params[where(strcmp(((*self.transits)[0]->get()).params.param,'Mid Transit') eq 1)].refined_mcmc_chain))    
  
     xtitle= ['Days','Degrees','Parameter Value','Parameter Value','Days since '+sigfig(xmin,8)]
     order = [0,1,2,3,4]
     for j=0,n_elements(order)-1,1 do begin
        if transit.params[order[j]].fixed eq 1 then begin
           !x.tickname=' '
           plot,[0,2],[0,2],/nodata,title=transit.params[order[j]].param,background=(*self.colors).white,$
                color=(*self.colors).black,thick=4,xthick=4,ythick=4,charthick=4
           xyouts,1,1,'Fixed',color=(*self.colors).black,charsize=4d,charthick=3d,align=0.5
           !x.tickname=''
        endif else begin
           modify = 0d0
           if strcmp(transit.params[order[j]].param,'Mid Transit') then modify = xmin
           xrange = tap_mm(transit.params[order[j]].refined_mcmc_chain,adjust=.1,shift=modify)
           if order[j] eq 1 then xrange[1] = min([xrange[1],90d0])
           plothist,(*transit.params[order[j]].refined_mcmc_chain)-modify, $
                    bin=tap_diff(transit.params[order[j]].refined_mcmc_chain)/bin, $
                    color=(*self.colors).black,background=(*self.colors).white, $
                    axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                    title=transit.params[order[j]].param, $
                    xrange=xrange,/xs, $
                    thick=4,xticks=3,xthick=4,ythick=4,charthick=4
           vline,transit.params[order[j]].mcmc_val[0]-modify,color=(*self.colors).blue,thick=4       
           vline,transit.params[order[j]].mcmc_val[0]-modify+transit.params[order[j]].mcmc_val[1], $
                 color=(*self.colors).blue,thick=4,linestyle=1
           vline,transit.params[order[j]].mcmc_val[0]-modify-transit.params[order[j]].mcmc_val[2], $
                 color=(*self.colors).blue,thick=4,linestyle=1
        endelse
     endfor
     self->plot_close

     self->plot_open,self.plot_dir+'plot_'+transit.fname+'_mcmc_params2.eps',xsize=7,ysize=2.7,/inches,/encapsulated,/color
     !p.multi = [0,2,1]
     !p.charsize *= (1d0/2)
     
     xtitle= ['Parameter Value','Parameter Value']
     order = [7,8]
     for j=0,n_elements(order)-1,1 do begin
        if transit.params[order[j]].fixed eq 1 then begin
           !x.tickname=' '
           plot,[0,2],[0,2],/nodata,title=transit.params[order[j]].param,background=(*self.colors).white,$
                color=(*self.colors).black,thick=4,xthick=4,ythick=4,charthick=4
           xyouts,1,1,'Fixed',color=(*self.colors).black,charsize=4d,charthick=3d,align=0.5
           !x.tickname=''
        endif else begin
           modify = 0d0
           if strcmp(transit.params[order[j]].param,'Mid Transit') then modify = xmin
           plothist,(*transit.params[order[j]].refined_mcmc_chain)-modify, $
                    bin=tap_diff(transit.params[order[j]].refined_mcmc_chain)/bin, $
                    color=(*self.colors).black,background=(*self.colors).white, $
                    axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                    title=transit.params[order[j]].param, $
                    xrange=tap_mm(transit.params[order[j]].refined_mcmc_chain,adjust=.1,shift=modify),/xs, $
                    thick=4,xticks=3,xthick=4,ythick=4,charthick=4
           vline,transit.params[order[j]].mcmc_val[0]-modify,color=(*self.colors).blue,thick=4       
           vline,transit.params[order[j]].mcmc_val[0]-modify+transit.params[order[j]].mcmc_val[1], $
                 color=(*self.colors).blue,thick=4,linestyle=1
           vline,transit.params[order[j]].mcmc_val[0]-modify-transit.params[order[j]].mcmc_val[2], $
                 color=(*self.colors).blue,thick=4,linestyle=1
        endelse
     endfor
     self->plot_close
     !p.charsize *= 2d0



     self->plot_open,self.plot_dir+'plot_'+transit.fname+'_mcmc_params3.eps',xsize=7,ysize=8,/inches,/encapsulated,/color
     !p.multi = [0,2,3]
     xtitle= ['Parameter Value','Parameter Value','Parameter Value','Parameter Value','Parameter Value','Parameter Value']
     order = [5,6,9,10,11,12]
     for j=0,n_elements(order)-1,1 do begin
        if transit.params[order[j]].fixed eq 1 then begin
           !x.tickname=' '
           plot,[0,2],[0,2],/nodata,title=transit.params[order[j]].param,background=(*self.colors).white,$
                color=(*self.colors).black,thick=4,xthick=4,ythick=4,charthick=4
           xyouts,1,1,'Fixed',color=(*self.colors).black,charsize=4d,charthick=3d,align=0.5
           !x.tickname=''
        endif else begin
           modify = 0d0
           if strcmp(transit.params[order[j]].param,'Mid Transit') then modify = xmin
           plothist,(*transit.params[order[j]].refined_mcmc_chain)-modify, $
                    bin=tap_diff(transit.params[order[j]].refined_mcmc_chain)/bin, $
                    color=(*self.colors).black,background=(*self.colors).white, $
                    axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                    title=transit.params[order[j]].param, $
                    xrange=tap_mm(transit.params[order[j]].refined_mcmc_chain,adjust=.1,shift=modify),/xs, $
                    thick=4,xticks=3,xthick=4,ythick=4,charthick=4
           vline,transit.params[order[j]].mcmc_val[0]-modify,color=(*self.colors).blue,thick=4       
           vline,transit.params[order[j]].mcmc_val[0]-modify+transit.params[order[j]].mcmc_val[1], $
                 color=(*self.colors).blue,thick=4,linestyle=1
           vline,transit.params[order[j]].mcmc_val[0]-modify-transit.params[order[j]].mcmc_val[2], $
                 color=(*self.colors).blue,thick=4,linestyle=1
        endelse
     endfor
     self->plot_close
     !y.tickname=''      


     ;spawn,'open '+self.plot_dir+'plot_'+transit.fname+'_mcmc_LD.eps'
     
     
     
  endfor
  !y.tickname=''    
                                ; stop
  
  
  
  
  
  !p = tempp
  !x = tempx
  !y = tempy
  tempp = 0L
  tempx = 0L
  tempy = 0L

 ; stop
end


pro TAP::plot_close
  device, /CLOSE_FILE
  set_plot, self.base_plot
end


function TAP::makeline,array
  line = ' & '
  if 0 then begin
     print,array
     stop
  endif
  if array[1] eq -1 then begin
     i=1
     while (sigfig(array[0],i) ne array[0]) and i lt 15 do i++
     line += sigfig(array[0],i)
     line += '\tablenotemark{a}' 
  endif else begin
     if abs(array[0]) ge 1 then begin
        line += sigfig(array[0],strlen(sigfig(array[1],2))-2+strlen(strtrim(floor(array[0]),2)))
        line +=  ' $^{+'+sigfig(array[1],2)+'}_{-'+sigfig(array[2],2)+'}$'
     endif else begin
     ;   print,'needs fix'
        i = 0d0
        while (abs(array[0])*(10^i) lt 1) do i++
        if array[0] lt 0 then i++
        line += sigfig(array[0],strlen(sigfig(array[1],2))-2-i+strlen(strtrim(floor(array[0]),2)))
        line +=  ' $^{+'+sigfig(array[1],2)+'}_{-'+sigfig(array[2],2)+'}$'
     endelse   
  endelse
  if 0 then print,line
  array = 0L
  return,line
end  

pro TAP::create_ascii
  self.message = 'Writing MCMC results to ascii files...'
  self->message

  openw,strlun,self.plot_dir+'MCMC_tables.tex',/get_lun,width=2500,bufsize=0
  printf,strlun,'\documentclass[preprint,10pt]{aastex}'
  printf,strlun,'\usepackage{graphicx}'
  printf,strlun,'\usepackage{epstopdf}'
  printf,strlun,'\begin{document}'
  printf,strlun,""
  printf,strlun,"\begin{deluxetable}{lc}"  
  printf,strlun,"\tablewidth{0pt}"
  printf,strlun,"\tablecaption{Overview of TAP MCMC Parameters}"
  printf,strlun,"\tablehead{"
  printf,strlun,"\colhead{Parameter} & \colhead{Value}}"
  printf,strlun,"\startdata"
  printf,strlun,"TAP version & "+self.version+'\\'  
  printf,strlun,"TAPmcmc version & "+((*self.transits)[0]->get()).mcmc_version+'\\'
  printf,strlun,'\hline'
  printf,strlun,"MCMC Chains & "+string(self.mcmc_complete,format='(i)')+'\\'
  printf,strlun,"Chain Length & "+string((((*self.transits)[self.active_transit-1])->get()).mcmc_params[1],format='(i)')+'\\'
  printf,strlun,"Total Inference Links & "+string(10d0*(n_elements(*((*self.transits)[0]->get()).params[0].refined_mcmc_chain)), $
                                                  format='(i)')+'\\ \hline'
  rebin=[-1]
  texp = [-1]
  texp1 = [-1]
  for i=0,n_elements(*self.transits)-1,1 do begin
     rebin = [rebin,((*self.transits)[i]->get()).rebin]
     texp  = [texp, ((*self.transits)[i]->get()).t_int[0]]
     texp1  = [texp1, ((*self.transits)[i]->get()).t_int[1]]
  endfor
  
  rebin_txt1 = ''
  rebin_txt2 = ''
  rebin_tbl = 0
  if robust_sigma(rebin[1:self.num_transits]) eq 0 then begin
     case rebin[1] of 
        0: rebin_txt = 'No Re-sampling'
        1: begin
           rebin_txt = 'Resample and Rebin Mandel \& Agol'
           case robust_sigma(texp[1:self.num_transits]) of
              0: case robust_sigma(texp1[1:self.num_transits]) of
                 0: rebin_txt1 = string(texp[1],format='(d9.6)')+' minutes, '+string(texp1[1],format='(i2.2)')+' Samples'
              endcase
              else: begin
                 rebin_txt1 = 'Table \ref{tbl:rebin} values.'
                 rebin_tbl = 1
              end
           endcase
        end
     endcase
  endif else begin
     rebin_txt = 'Varied, see Table \ref{tbl:rebin}'
     rebin_tbl=1
  endelse
  
  printf,strlun,"Long Int Mode: & "+rebin_txt+'\\'
  if strcmp(rebin_txt1,'') eq 0 then printf,strlun,"Resample technique & "+rebin_txt1+'\\' 

  printf,strlun,'\hline'

  for i=0,n_elements(*self.transits)-1,1 do printf,strlun,"Transit "+ $
     string(i+1,format='(i2.2)')+' & '+self->stripspecial(((*self.transits)[i]->get()).fname)+'\\'
  printf,strlun,"\enddata"  
  printf,strlun,"\tablecomments{TAP MCMC "+self.version+"}"
  printf,strlun,"\label{tbl:mcmcpar}"
  printf,strlun,"\end{deluxetable}"
  printf,strlun,""


  ;stop
  
  if rebin_tbl then begin
     printf,strlun,""
     printf,strlun,"\begin{deluxetable}{lcc}"  
     printf,strlun,"\tablewidth{0pt}"
     printf,strlun,"\tablecaption{Overview of TAP MCMC Rebin Parameters}"
     printf,strlun,"\tablehead{"
     printf,strlun,"\colhead{Transit} & \colhead{Mode} & \colhead{Technique}}"
     printf,strlun,"\startdata"
     for i=0,n_elements(*self.transits)-1,1 do begin
        transit = string(i+1,format='(i2.2)')+': '+self->stripspecial(((*self.transits)[i]->get()).fname)
        case rebin[i+1] of
           0: begin
              mode = 'No Re-sampling'
              technique = '\nodata'
           end
           1: begin
              mode = 'Resample and Rebin Mandel \& Agol'
              technique = string(texp[i+1],format='(d9.6)')+' minutes, '+string(texp1[i+1],format='(i2.2)')+' Samples'
           end
        endcase
        printf,strlun,transit+' & '+mode+' & '+technique+' \\'
     endfor
   ;  printf,strlun,"\hline"
     

   ;  stop
     

     printf,strlun,"\enddata"  
     printf,strlun,"\tablecomments{TAP MCMC "+self.version+"}"
     printf,strlun,"\label{tbl:rebin}"
     printf,strlun,"\end{deluxetable}"
     printf,strlun,""
  endif
  mode = 0L
  technique = 0L
  transit = 0L
  rebin_txt = 0L
  rebin_txt1 = 0L
  rebin = 0L
  texp =  0L
  
  
  if self.create_plots then begin
     printf,strlun,'\begin{figure*}[h]'
     printf,strlun,'\centering'
     printf,strlun,'\includegraphics[]{'+self.plot_dir+'plot_alltransit_phase_modandresid.eps}'
     printf,strlun,'\caption{TAP MCMC '+self.version+'}'
     printf,strlun,'\end{figure*}'


     printf,strlun,'\begin{figure*}[h]'
     printf,strlun,'\centering'
     printf,strlun,'\includegraphics[]{'+self.plot_dir+'plot_alltransit_phased_lc.eps}'
     printf,strlun,'\caption{TAP MCMC '+self.version+'}'
     printf,strlun,'\end{figure*}'


     printf,strlun,'\begin{figure*}[h]'
     printf,strlun,'\centering'
     printf,strlun,'\includegraphics[]{'+self.plot_dir+'plot_alltransit_lightcurve.eps}'
     printf,strlun,'\caption{TAP MCMC '+self.version+'}'
     printf,strlun,'\end{figure*}'
  endif

  write_table = [-1,-1]
  while write_table[1] lt n_elements(*self.transits)-1 do begin

     format = ''
     fixed = 1
     
     write_table[0] = write_table[1]+1
     write_table[1] = min([write_table[0]+2,n_elements(*self.transits)-1])
     for i=write_table[0],write_table[1],1 do format = format+'c'
     printf,strlun,""
     printf,strlun,"\begin{deluxetable}{l"+format+"}"  
                                ; if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
     printf,strlun,"\tablewidth{0pt}"
     printf,strlun,"\tablecaption{Wavelet basis red noise MCMC analysis}"
     printf,strlun,"\tablehead{"
     printf,strlun,"\colhead{Parameter}"
     printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{Value} \\"
     for i=write_table[0],write_table[1],1 do printf,strlun,$
        "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,15)+" }"
     printf,strlun,"}"
     printf,strlun,"\startdata"
     for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
        line = string(((*self.transits)[0]->get()).params[i].param,format='(A15)')
        for j=write_table[0],write_table[1],1 do begin
           line += self->makeline(((*self.transits)[j]->get()).params[i].mcmc_val) 
        endfor
        line += '\\'
        printf,strlun,line  
     endfor
     
     printf,strlun,"\enddata"
     if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
     
     printf,strlun,"\tablecomments{TAP MCMC "+self.version+"}"
     printf,strlun,"\label{tbl:tapmcmc1}"
     printf,strlun,"\end{deluxetable}"
     
  endwhile
  
  write_table = [-1,-1]
  while write_table[1] lt n_elements(*self.transits)-1 do begin
     
     format = ''
     fixed = 1
     
     write_table[0] = write_table[1]+1
     write_table[1] = min([write_table[0]+5,n_elements(*self.transits)-1])
     for i=write_table[0],write_table[1],1 do format = format+'c'
     printf,strlun,""
     ;; PARAMETER LOCKS
     printf,strlun,"\begin{deluxetable}{l"+format+"}"    ;;if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
     printf,strlun,"\tablewidth{0pt}"
     printf,strlun,"\tablecaption{Multi Curve MCMC Parameter Lock Matrix}"
     printf,strlun,"\tablehead{"
     printf,strlun,"\colhead{Parameter}"
     printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{Jump Rate (Requested)} \\"
     for i=write_table[0],write_table[1],1 do printf,strlun,$
  "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,5)+"... }"
     printf,strlun,"}"
     printf,strlun,"\startdata"
     for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
        line = string(((*self.transits)[0]->get()).params[i].param,format='(A15)')
        for j=write_table[0],write_table[1],1 do begin
           if ((*self.transits)[j]->get()).params[i].mcmc_val[1] eq -1 then $
              line += '& \nodata' else $
                 line += ' & ' + string((((*self.transits)[j]->get()).params[i].jumpct/((*self.transits)[j]->get()).params[i].jumptot),format='(d4.2)') +$
                         ' ('+string(((*self.transits)[j]->get()).params[i].accept,format='(d4.2)')+')'
        endfor
        line += '\\'
        printf,strlun,line  
     endfor
     
     printf,strlun,"\enddata"
  ;   if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
     
     printf,strlun,"\tablecomments{TAP MCMC "+self.version+".  Jump rates (with request).}"
     printf,strlun,"\label{tbl:tapmcmc2}"
     printf,strlun,"\end{deluxetable}"
     
  endwhile


  
  write_table = [-1,-1]
  while write_table[1] lt n_elements(*self.transits)-1 do begin
     
     format = ''
     fixed = 1
     
     write_table[0] = write_table[1]+1
     write_table[1] = min([write_table[0]+5,n_elements(*self.transits)-1])
     for i=write_table[0],write_table[1],1 do format = format+'c'
     printf,strlun,""
     ;; PARAMETER LOCKS
     printf,strlun,"\begin{deluxetable}{l"+format+"}"    ;;if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
     printf,strlun,"\tablewidth{0pt}"
     printf,strlun,"\tablecaption{Multi Curve MCMC Parameter Lock Matrix}"
     printf,strlun,"\tablehead{"
     printf,strlun,"\colhead{Parameter}"
     printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{MCMC Parameter Set} \\"
     for i=write_table[0],write_table[1],1 do printf,strlun,$
  "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,5)+"... }"
     printf,strlun,"}"
     printf,strlun,"\startdata"
     for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
        line = string(((*self.transits)[0]->get()).params[i].param,format='(A15)')
        for j=write_table[0],write_table[1],1 do begin
           line += ' & ' + string(((*self.transits)[j]->get()).params[i].set,format='(i3)')
           if ((*self.transits)[j]->get()).params[i].mcmc_val[1] eq -1 then line += '\tablenotemark{a}' 
        endfor
        line += '\\'
        printf,strlun,line  
     endfor
     
     printf,strlun,"\enddata"
     if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
     
     printf,strlun,"\tablecomments{TAP MCMC "+self.version+".  Transits with the same values for any parameter row are locked together in the MCMC analysis.}"
     printf,strlun,"\label{tbl:tapmcmc2}"
     printf,strlun,"\end{deluxetable}"
     
     endwhile




  printf,strlun,"\end{document}"
  close,strlun
  
  openw,strlun,self.plot_dir+'ascii_phased_data.ascii',/get_lun,width=2500,bufsize=0
  printf,strlun,'#  Phase [Days from Tmid]    Flux_corr      Raw_Flux  Orig_T_days     model_flux  Residual'
  printf,strlun,"#   tap_readcol,'"+self.plot_dir+"ascii_phased_data.ascii',p,f,rf,t,m,r,format='(d,d,d,d,d,d)'"
  
  compare = 'Airmass Y-int'
  compare2 = 'Airmass Slope'
  for i=0,self.num_transits-1,1 do begin
     lc = (*self.transits)[i]->get() 
     midt = lc.params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
     fit = ptr_new(poly((lc.transit_tday-min(lc.transit_tday)),[lc.params[where(strcmp(lc.params.param,compare) eq 1)].value,lc.params[where(strcmp(lc.params.param,compare2) eq 1)].value]))
     for j=0,n_elements(lc.transit_tday)-1,1 do begin
        
        printf,strlun,(lc.transit_tday[j]-midt),lc.transit_flux[j]/(*fit)[j],lc.transit_flux[j],lc.transit_tday[j],(*lc.model_f)[j],lc.transit_flux[j]-(*lc.model_f)[j],format='(d20.8,d20.8,d20.8,d25.8,d20.8,d20.8)'
     endfor
     ptr_free,fit
  endfor
  close,strlun
  compare = 0L
  compare2 = 0L
  


  
  transit = (*self.transits)[0]->get()
  tpar = transit.params.value
  tpar[9] = 1d0
  tpar[10] = 0d0
;;  dum = TAP_MPFit_function(tpar,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,redfit=redfit,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int)
  dum = TAP_MPFit_function(tpar,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int)

  openw,strlun,self.plot_dir+'ascii_phased_model.ascii',/get_lun,width=2500,bufsize=0
  printf,strlun,'# Phased from '+transit.fname+' Parameters'
  printf,strlun,'# T_days   Phase [hours from Tmid]     Model Flux'
  printf,strlun,"#   tap_readcol,'"+self.plot_dir+"ascii_phased_model.ascii',t,p,mf,format='(d,d,d)'"
  
  midt = transit.params[where(strcmp(transit.params.param,'Mid Transit') eq 1)].value - (transit.epoch*transit.params[0].value)
  for j=0,n_elements(ffit)-1,1 do printf,strlun,(*transit.model_tfine)[j],((*transit.model_tfine)[j]-midt),ffit[j],format='(d20.8,d12.8,d12.8)'


  close,strlun
  midt=0L
  dum = 0L
  ffit=0L
  dfit=0L
  transit = 0L
  tpar = 0

  
  
  
  if self.create_mcmcascii then begin
     self.message = '   ...Writing MCMC chains to ascii files'
     self->message
     
     for i=0,n_elements(*self.transits)-1,1 do begin
        transit = (*self.transits)[i]->get()
        openw,strlun,self.plot_dir+'ascii_'+transit.fname+'_MCMC.ascii',/get_lun,width=2500,bufsize=0
        
        head = '## TAP '+self.version+' Analysis of '+string(transit.mcmc_params[0],format='(i2.2)')+' MCMC Chains of '+transit.fname
        printf,strlun,head
        ;; head = '##   '+string(self.num_combined,format='(i)')+' MCMC chains combined.'
        ;; printf,strlun,head
        head = '##   Burn-in stripped, every 10th link saved, leaving '+ $
               string(n_elements(*transit.params[0].refined_mcmc_chain),format='(i)')+ $
               ' links in this file, representative of '$
               +string(n_elements(*transit.params[0].refined_mcmc_chain)*10d0,format='(i)')+ ' total MCMC links.'
        printf,strlun,head
        head = "## tap_readcol,'"+self.plot_dir+"ascii_"+transit.fname+"_MCMC.ascii',P,i,adr,rdr,tmid,mu1,mu2,ecc,omega,yint,slope,sigr,sigw,format='d,d,d,d,d,d,d,d,d,d,d,d,d'"
        printf,strlun,head
        
        head = '##     Period             Inclination            a/Rs              Rp/Rs                    '+$
               'Mid-Transit          LD_u1_linear       LD_u2_quadratic       Eccentricity          Omega   '+$
               'Linear_Y-int        Linear_Slope           Sigma_R             Sigma_W'
        printf,strlun,head
        
        format='(d19.15,d20.15,d20.15,d20.15,d30.15,d20.15,d20.15,d20.15,d20.15,d20.15,d20.15,d20.15,d20.15)'
        for j=0d0,n_elements(*transit.params[0].refined_mcmc_chain)-1d0,1d0 do begin
           printf,strlun,(*transit.params[0].refined_mcmc_chain)[j], $
                  (*transit.params[1].refined_mcmc_chain)[j], $          
                  (*transit.params[2].refined_mcmc_chain)[j], $
                  (*transit.params[3].refined_mcmc_chain)[j], $          
                  (*transit.params[4].refined_mcmc_chain)[j], $
                  (*transit.params[5].refined_mcmc_chain)[j], $          
                  (*transit.params[6].refined_mcmc_chain)[j], $
                  (*transit.params[7].refined_mcmc_chain)[j], $          
                  (*transit.params[8].refined_mcmc_chain)[j], $
                  (*transit.params[9].refined_mcmc_chain)[j], $          
                  (*transit.params[10].refined_mcmc_chain)[j], $
                  (*transit.params[11].refined_mcmc_chain)[j], $          
                  (*transit.params[12].refined_mcmc_chain)[j], $
                  format=format 
        endfor
        
        close,strlun
        transit = 0L
        
     endfor
  endif
  
  
end

pro TAP::LoadTransit,event
  widget_control,self.transitfile_fld[1],get_value=path

  case self.filetype of
     'IDL Save File': BEGIN
        self.num_transits+=1
        self.active_transit+=1
        
        restore,path
        self->preptransit,input=lc,type='struct'
     end
     'ASCII File': Begin
        tap_readcol,path,temp1,temp2,format='d,d'
        
        bounds = where(temp1 eq -1 and temp2 eq -1)
        if bounds[0] eq -1 then begin
           ntransit = 1
           bounds=[-1,n_elements(temp1)]
        endif else begin
           ntransit = n_elements(where(temp1 eq -1 and temp2 eq -1))+1
           bounds=[-1,bounds,n_elements(temp1)]
        endelse
        self.message = 'Detected '+string(ntransit,format='(i2.2)')+' light curves.'
        self->message
        for i=0,ntransit-1,1 do begin       
           self.num_transits+=1
           self.active_transit+=1
           
           t1 = temp1[bounds[i]+1:bounds[i+1]-1]
           t2 = temp2[bounds[i]+1:bounds[i+1]-1]
     
           struct= replicate({hjd:0d0, f: 0d0, e:0d0},n_elements(t1))
           struct.hjd=t1[sort(t1)]
           struct.f = t2[sort(t1)]
           
           self->preptransit,input=struct,type='text'
           line=0L
           struct = 0L
           wc = 0L
        endfor
        temp1=0L
        temp2=0L
     end
  endcase
  self->setup,'multi'
  self->setup,'cross lc links'
  self->setup,'fit'

  self.message = 'light curve(s) loaded.'
  self->message

  ;; self->setup,'inference'
end

pro TAP::OrganizeTransits,input  
  if self.num_transits eq 1 then *self.transits = obj_new('transit') else *self.transits = [*self.transits,obj_new('transit')]

  transit = ptr_new(self->prepTransitStruct(input))
  (*transit).rebin = self.setup_rebin
  (*transit).t_int = self.smooth_val
  (*transit).epoch = 0
  (*self.transits)[self.active_transit-1]->set,*transit

  input = 0L
  
  ptr_free,transit
  self->setupmodel
end

pro TAP::setupmodel
  transit = (*self.transits)[self.active_transit-1]->get()
 
  self->setparinfo
  (*self.parinfo)[0].fixed = 1 
  x = mpfit('TAP_MPFit_function',$
            transit.params.value,$
            functargs={time:*transit.model_tfine ,$
                       flux:*transit.model_ffine ,$
                       tdat:transit.transit_tday, $   
                       fdat:transit.transit_flux, rebin:transit.rebin, smooth:transit.t_int $
                      },parinfo=(*self.parinfo),quiet=1)
  
  transit.params.value = x
  x = 0L
  (*self.transits)[self.active_transit-1]->set,transit
  
  transit = 0L
  self->updateModel
  end


pro TAP::adjustslider,event
  widget_control, event.id, GET_UVALUE= uvalue
  transit = (*self.transits)[self.active_transit-1]->get()
  transit.params[where(strcmp(transit.params.param,uvalue.value) eq 1)].value = event.value
  (*self.transits)[self.active_transit-1]->set,transit
  transit = 0L
  self->updatemodel
end

pro TAP::updateModel
  transit = (*self.transits)[self.active_transit-1]->get()
  links = transit.params.set
  params= transit.params.value
  transit=0L

  for i=0,self.num_transits-1,1 do begin
     transit = (*self.transits)[i]->get()
     if (where(transit.params.set eq links))[0] ne -1 then begin
        transit.params[where(transit.params.set eq links)].value = params[where(transit.params.set eq links)]
        redfit = dblarr(n_elements(transit.rednoise))
;;        dum = TAP_MPFit_function(transit.params.value,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,redfit=redfit,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int)
        dum = TAP_MPFit_function(transit.params.value,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int)
        transit.rednoise = redfit
        *transit.model_f = dfit*1d0
        transit.residuals = transit.transit_flux - dfit
        *transit.model_ffine = ffit
        redfit=0L
        ffit = 0L
        dfit = 0L
        dum = 0L
        
        if abs((transit.params[4].value - median(transit.transit_tday))-transit.epoch*transit.params[0].value) gt .5*transit.params[0].value then begin
           epoch = fillarr(1,-5000,5000)
           transit.epoch = epoch[where(abs((transit.params[4].value - median(transit.transit_tday))-epoch*transit.params[0].value) eq min(abs((transit.params[4].value - median(transit.transit_tday))-epoch*transit.params[0].value)))]
           epoch = 0L
        endif
        
        (*self.transits)[i]->set,transit
     endif
     transit = 0L
  endfor
  links = 0L
  params = 0L
  self->lcplot
  self.numcol=2
  self->PrepCol
end

pro TAP::adjustxlock_event,event
  widget_control, event.id, GET_UVALUE= uvalue
  if n_elements(*event.value) eq 0 then (*event.value) = 0
  
  case uvalue.type of
     'val':begin
        transit = (*self.transits)[uvalue.value]->get()
        transit.params[uvalue.param].set = (*event.value)
        (*self.transits)[uvalue.value]->set,transit
        transit=0L
        self.active_transit = uvalue.value+1
        self->updatemodel
     end
     else: print,"unknown adjustxlock_event"
  endcase
end

pro TAP::SetParInfo
  ptr_free,self.parinfo

  self.parinfo = ptr_new(replicate({fixed: 0, $
                                    limited: [0,0], $
                                    limits: [0.,0.], $
                                    relstep: 0.d0},$
                                   13))
  
  (*self.parinfo).fixed = ((*self.transits)[self.active_transit-1]->get()).params.fixed
  (*self.parinfo).limited = ((*self.transits)[self.active_transit-1]->get()).params.limited
  (*self.parinfo).limits = ((*self.transits)[self.active_transit-1]->get()).params.limits
  
  transit = (*self.transits)[self.active_transit-1]->get()
  if transit.rebin gt 0 then begin
     (*self.parinfo)[5].fixed = 1 
     (*self.parinfo)[6].fixed = 1 
     (*self.parinfo)[9].fixed = 1 
     (*self.parinfo)[10].fixed = 1 
  endif
  transit = 0L
  
end

function tap::window
  return,(*self.plot_windows)[0]
end

function TAP::PrepTransitStruct,input
  
  scale = (input.hjd)[1]-(input.hjd)[0]
                                ; t_final = ptr_new(makearr(((mm(double(input.hjd)))[1]-(mm(double(input.hjd)))[0])*(24d0*60d0), mm(double(input.hjd))+[-10d0*scale,10d0*scale]))
  t_final = ptr_new(fillarr((1d0/(24d0*60d0)),mm(double(input.hjd))+[-5d0*scale,5d0*scale]))
  
  case self.setup_rebin of
     0: t = ptr_new(input.hjd)
     1: t = ptr_new((reform((dblarr(self.smooth_val[1])+1)#input.hjd $
                            +(findgen(self.smooth_val[1])+1-5d-1*(self.smooth_val[1]+1d0))*$
                            (self.smooth_val[0]/self.smooth_val[1])/1440d0 $
                            #(dblarr(n_elements(input.hjd))+1),self.smooth_val[1]*n_elements(input.hjd),1))[*,0])
  endcase  


  transit = {transit_tday: input.hjd  , $
             transit_flux: input.f    , $ 
             params: replicate({param: '', limits: [0d0,0d0], fixed: 0, limited: [0,0], $
                                value: 0d0, accept: 0.44d0, beta: 44d-2, set:self.active_transit, $
                                b_last: 0d0, sval: 2d0, prior: [0d0,0d0,0d0], $
                                mcmc_chain:ptr_new([-1]), jumpct: 0d0, jumptot: 0d0, $
                                curr_link: 0d0, new_link: 0d0, mcmc_val: dblarr(3), $
                                refined_mcmc_chain: ptr_new(), runval: dblarr(3)},13) ,$
             epoch: 0,$
             new_redl: 0d0,$
             curr_redl: 0d0,$
             mcmc_links: 0d0,$
             model_t: ptr_new(*t_final) ,$
             model_f: ptr_new(*t_final*0d0),$
             model_tfine: ptr_new(*t),$
             model_ffine: ptr_new((*t)*0d0),$
             rebin: 0L, $
             t_int: dblarr(2),$
             residuals: dblarr(n_elements(input.hjd)) ,$
             rednoise:  dblarr(n_elements(input.hjd)) ,$
             color: 0L ,$
             modcol: 0L ,$
             redcol: 0L ,$
             fname: '',$
             mcmc_params: [5,1d6,5d6],$
             mcmc_complete: 0L ,$
             mcmc_files: ptr_new('-1') ,$
             mcmc_version: '' ,$
             diff: 0d0 $
            } 
  
  widget_control,self.transitfile_fld[1],get_value=path
  transit.color = (*self.colors).(self.active_transit+2)
  transit.modcol = (*self.colors).blue
  transit.redcol = (*self.colors).red
  if n_elements(strsplit(path,'/')) gt 1 then $
     transit.fname = (strsplit(((strsplit(path,'/',/extract))[n_elements(strsplit(path,'/',/extract))-1]),'.',/extract))[0] else $
        transit.fname = (strsplit(((strsplit(path,'\',/extract))[n_elements(strsplit(path,'\',/extract))-1]),'.',/extract))[0]

;;print,transit.fname
  
  
  for i=1,n_elements(*self.label_indices)-1,1 do transit.params[i-1].param = strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)
  
  input = 0L

;;  compare = ''
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].min =
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].max = 
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].locked = 
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited = 

  compare = 'Period'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = self.init_period
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [.5d0,15d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited =  [0,0]

  compare = 'Inclination'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 88d0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [80d0, 100d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited =  [0,0]


  compare = 'a/R*'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 9d0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [1d0,60d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited = [0,0]


  compare = 'Rp/R*'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = .1d0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [.001d0,.2d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited = [0,0]

  
  compare = 'Mid Transit'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = transit.transit_tday[where(min(transit.transit_flux) eq transit.transit_flux)]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [min(transit.transit_tday),max(transit.transit_tday)]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited = [0,0]


  compare = 'Linear LD'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = .2d0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [0d0, 1d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited =  [1,1]

  compare = 'Quad LD'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = .2d0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [-1d0, 1d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited =  [1,1]

  compare = 'Eccentricity'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 0d0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [0d0,.1d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 1
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited = [0,0]

  compare = 'Omega'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 0d0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [0d0,3.14d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 1
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited =  [0,0]

  compare = 'Airmass Y-int'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 1d0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [.9999d0,1.0001d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited =  [0,0]


  compare = 'Airmass Slope'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 0d0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [-0.001d0,0.001d0]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited =  [0,0]


  compare = 'Sigma Red'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 1d-4
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [0d0,.5d-1]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited = [1,1]

  compare = 'Sigma White'
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 1d-3
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limits = [0d0,1d-1]
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].fixed = 0
  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited = [1,1]

  compare=0L
 ; print_struct,transit.params
  scale = 0L
  ptr_free,t_final,t
  return,transit
end

pro TAP::Message,event
  widget_control,self.message_window,/append,$
                 set_value='('+curr_date(format='hh:mm:ss')+') '+self.message
  ;print,self.message
  self.message=''
end

pro TAP::PrepTransit,event,input=input,type=type
  case type of
     'struct' :BEGIN
        input.hjd = input[sort(input.hjd)].hjd
        input.f = input[sort(input.hjd)].f
        input.e = input[sort(input.hjd)].e
        self->organizetransits,input
        self.message='Imported existing Transit Structure'
        self->message
     END
     'text' : BEGIN
        self->organizetransits,input
 ;       self.message='Imported ASCII Transit File'
     END
  endcase
  


  self->lcplot
end


;function TAP::PrepModel,event
;  setp=ptr_new([-1, 0.d0])
;  
;  return,model
;end

pro TAP::MenuMap,event
  Widget_Control, event.id, Get_UValue=info
  
  case info.type of
     'main': begin
        for i=1,n_elements(*self.bases)-1,1 do widget_control,(*self.bases)[i],MAP=0
        case info.value of
          ; 'General': begin
          ;    widget_control,self.general_base,/map
          ; end
           'Load Transit': begin
             ; self->setup,'load'
              widget_control,   (*self.bases)[1],/map
           end
           'Manage Transits': begin
             ; self->setup,'multi'
              widget_control,   (*self.bases)[2],/map
           end
           'Fit': begin
             ; self->setup,'fit'
              widget_control,   (*self.bases)[3],/map
           end
           'MCMC Inference': begin
             ; self->setup,'inference'
              widget_control,   (*self.bases)[4],/map
           end
           else:  print,'Unknown MenuMap Event "'+info.value+'"'
        endcase
     end
  endcase
end

pro TAP::start
  centertlb,self.tap_base
  widget_control,self.tap_base,/realize
  widget_control,self.tap_base, XOffset=0, YOffset=0

  widget_control, (*self.plot_windows)[0].window, Get_Value=wid
  (*self.plot_windows)[0].w_id=wid
  
  window, xsize=(*self.plot_windows)[0].x $
          , ysize=(*self.plot_windows)[0].y $
          , /pixmap,/free
  (*self.plot_windows)[0].pix_window = !d.window
  self->lcplot
end

pro TAP::destroy
  self.message = 'Destroying Widget'
  self->message

  for i=0,n_elements(self.buttons)-1,1 do if (self.buttons)[i] then widget_control,(self.buttons)[i],/destroy
  for i=0,n_elements(*self.slider)-1,1 do  if (*self.slider)[i].id then widget_control,(*self.slider)[i].id,/destroy
  for i=0,n_elements(self.settings)-1,1 do if self.settings[i] then widget_control,self.settings[i],/destroy
  for i=0,n_elements(*self.bases)-1,1 do if (*self.bases)[i] then widget_control,(*self.bases)[i],/destroy
  if self.tap_base then widget_control,self.tap_base,/destroy
  ptr_free,self.colors
  ptr_free,self.slider
  ptr_free,self.plot_windows
  ptr_free,self.parinfo
  
  ptr_free,self.mcmc_stuff

  for i=0,n_elements(*self.extra_windows)-1,1 do if (*self.extra_windows)[i] then widget_control,(*self.extra_windows)[i],map=0
  for i=0,n_elements(self.fld[*,0])-1,1 do $
     for j=0,n_elements(self.fld[0,*])-1 do if self.fld[i,j] ne 0 then begin
     widget_control,self.fld[i,j],/destroy
     self.fld[i,j] = 0
  endif
  for i=0,n_elements(*self.extra_windows)-1,1 do if (*self.extra_windows)[i] then widget_control,(*self.extra_windows)[i],/destroy


  for i=0,self.num_transits-1,1 do begin
     transit = (*self.transits)[i]->get()   
     for k=0,n_elements(transit.params)-1,1 do begin
        ptr_free,transit.params[k].mcmc_chain
        ptr_free,transit.params[k].refined_mcmc_chain
     endfor
     ptr_free,transit.model_t,transit.model_f,transit.model_tfine,transit.model_ffine,transit.mcmc_files
     transit = 0L
     (*self.transits)[i]->destroy
  endfor
  
  ptr_free,self.transits
 
  ptr_free,self.extra_windows
  ptr_free,self.label_indices
  ptr_free,self.label
  ptr_free,self.menus
  ptr_free,self.bases

  obj_destroy,self
end


pro TAP_Cleanup,event

end


pro tap::PrepCol,runval=runval
  cols = self.numcol
  
  if keyword_set(runval) then begin
     cols = 4
     col2 = ['Value',string(((*self.transits)[runval-1]->get()).params.runval[0],format='(d17.8)')]
     col3 = ['85.1%',string(((*self.transits)[runval-1]->get()).params.runval[1],format='(d9.5)')]
     col4 = ['15.9%',string(((*self.transits)[runval-1]->get()).params.runval[2],format='(d9.5)')] 
     
     badformat = where(strcmp(strmid(col3,0,1),'*'))
     if badformat[0] ne -1 then col3[badformat] = string(((*self.transits)[runval-1]->get()).params[badformat].runval[1],format='(e9.2)')
     badformat = 0L
     badformat = where(strcmp(strmid(col4,0,1),'*'))
     if badformat[0] ne -1 then col4[badformat] = string(((*self.transits)[runval-1]->get()).params[badformat].runval[2],format='(e9.2)')
     badformat = 0L
     
  endif else begin
     if cols gt 1 then begin
        col2 = ['Value',string(((*self.transits)[self.active_transit-1]->get()).params.value,format='(d17.8)')]
        if cols gt 2 then begin
           col3 = ['85.1%',string(((*self.transits)[self.active_transit-1]->get()).params.mcmc_val[1],format='(d9.5)')]
           col4 = ['15.9%',string(((*self.transits)[self.active_transit-1]->get()).params.mcmc_val[2],format='(d9.5)')]

           badformat = where(strcmp(strmid(col3,0,1),'*'))
           if badformat[0] ne -1 then col3[badformat] = string(((*self.transits)[self.active_transit-1]->get()).params[badformat].mcmc_val[1],format='(e9.2)')
           badformat = 0L
           badformat = where(strcmp(strmid(col4,0,1),'*'))
           if badformat[0] ne -1 then col4[badformat] = string(((*self.transits)[self.active_transit-1]->get()).params[badformat].mcmc_val[2],format='(e9.2)')
           badformat = 0L
        endif
     endif
  endelse
  
  
  line = ''
  case cols of 
     1:  for i=0,n_elements(*self.label)-1,1 do $
        widget_control,(*self.label)[i],set_value=string((*self.label_indices)[i],format='(A14)') 
     2:  for i=0,n_elements(*self.label)-1,1 do $
        widget_control,(*self.label)[i],set_value=string((*self.label_indices)[i],col2[i],format='(A14,A17)')
     4: begin
        if (where(1d0*col3[1:n_elements(col3)-1] eq -1))[0] ne -1 then col3[where(1d0*col3[1:n_elements(col3)-1] eq -1)+1] = string('Fixed',format='(A8)')
        if (where(1d0*col4[1:n_elements(col4)-1] eq -1))[0] ne -1 then col4[where(1d0*col4[1:n_elements(col4)-1] eq -1)+1] = string('Fixed',format='(A8)')
        for i=0,n_elements(*self.label)-1,1 do $
           widget_control,(*self.label)[i],set_value=string((*self.label_indices)[i],col2[i],col3[i],col4[i],format='(A14,A17,A10,A10)')
     end
     else: print,'DEBUG: PrepCol unknown # of cols:',cols
     
  endcase
  
  cols=0L
  col1=0L
  col2=0L
  col3=0L
  col4=0L
  
end

pro tap::setup,set_this
  case set_this of 
     'load': begin
        ys = 80

        (*self.bases)[1] = widget_base((*self.bases)[0] ,$
                                        /column,$
                                        MAP=0)

        bigrow =  widget_base( (*self.bases)[1] ,$
                           /column,$
                           /Base_align_center,frame=1)
        row = widget_base( bigrow,$
                           /ROW,$
                           /Base_align_center,frame=0)
        
        button = widget_button(row,$
                               font=self.buttonfont,$
                               xsize=100,$
                               value='Transit File',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Transit File Button'})
    
        path = coyote_field2(row,$
                             labelfont=self.buttonfont,$
                             FIELDFONT=self.textfont,$
                             TITLE=':',$
                             VALUE = self.transitpath,$
                             UVALUE='Transit File Field',$
                             XSIZE=65,$
                             TEXTID=textid)
        self.transitfile_fld = [path,textid]
        
        clear = widget_button(row,$
                              FONT=self.buttonfont,$
                              VALUE='Clear',$
                              uvalue={object:self, $
                                      method:'ButtonEvent', $
                                      value:'Clear Data Path'})
        
        mrow = widget_base(bigrow,$
                           /ROW,$
                           /Base_align_center)
        

        if 0 then begin
           row = widget_base(mrow,frame=1,$
                             /column,ysize=ys,$
                             /Base_align_left)
           bg = cw_bgroup(row,$
                          ['ASCII File','IDL Save File'],$
                          /column,$
                          LABEL_top='File Type:',$
                          /RETURN_NAME,$
                          /NO_RELEASE,$
                          UVALUE={object:self, method:'ButtonEvent', value:'FileType'},$      
                          FONT=self.buttonfont,$
                          /EXCLUSIVE,$
                          SET_VALUE=0)   
        endif
        self.filetype = 'ASCII File'
        

        row = widget_base(mrow,frame=1,$
                          /column,ysize=ys,$
                          /Base_align_left)
        
                                ; row = widget_base(mrow,$
                                ;                   /column,frame=1,$
                                ;                   /Base_align_center)
        
                                ; radio = widget_base(row,column=1,/nonexclusive)
                                ;  label = widget_label(row,value='Long integrations?')
      ;  options = ['None','Rebin to data cadence','Rebin to "Input Integration"']
        options = ['None','Rebin to "Input Integration"']
        bg = cw_bgroup(row,$
                       options,$
                       /column,$
                       LABEL_top='Long Integrations?',$
                       /RETURN_NAME,$
                       /NO_RELEASE,$
                       set_value =self.setup_rebin,$
                       UVALUE={object:self, method:'ButtonEvent', value:'RebinType'},$      
                       FONT=self.buttonfont,$
                       /EXCLUSIVE)
                                ;self.setup_rebin = 0
        options = 0L
        
        if self.setup_rebin eq 1 then x = 1 else x = 0
        self.setup_smooth = widget_base(mrow,frame=1,$
                                        /column,$
                                        /Base_align_left, sensitive=x,ysize=ys)
        
        lbl = widget_label(self.setup_smooth,value='Input Integration')
        
        row = widget_base(self.setup_smooth,/row)
        fld = coyote_field2(row,$
                            TITLE='Minutes:',$
                            UVALUE={object:self, method:'ButtonEvent', $
                                    value:'settint',set:0},$
                            decimal=4, digits=8, $
                            VALUE=self.smooth_val[0],$
                            XSIZE=10,/doubleValue,event_pro='TAP_event',$
                            textid=textid)    
        
     ;   lbl = widget_label(self.setup_smooth,value='Samples per Integration')
        
        fld = coyote_field2(row,$
                            TITLE='N_samp:',$
                            UVALUE={object:self, method:'ButtonEvent', $
                                    value:'settint',set:1},$
                            decimal=4, digits=8, $
                            VALUE=round(self.smooth_val[1]),$
                            XSIZE=5,/integerValue,event_pro='TAP_event',$
                            textid=textid)    
        
        
        row = widget_base(mrow,frame=0,$
                          /column,$
                          /Base_align_left)
        
                                ;button = widget_button(radio,value='Long integrations',$
                                ;uname = ,$ 
         ;                      uvalue={object:self,$
          ;                             method:'setup_rebin',$
           ;                            value:0 ,$
            ;                           type: ''$
             ;                         })
        ;widget_control,button,set_button=self.setup_rebin
        

       ; row = widget_base(mrow,$
       ;                   /ROW,$
       ;                   /Base_align_center,ysize=40)
        
        self.buttons[2] =  widget_button(row,$
                                font=self.buttonfont,$
                                xsize=100,$
                                value='Load Transit',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Load Transit Button'},sensitive=0)
       ; widget_control,   (*self.bases)[1],/map
              

        row = widget_base( (*self.bases)[1] ,$
                           /ROW,$
                           /Base_align_center,frame=1)
        
        button = widget_button(row,$
                               font=self.buttonfont,$
                               xsize=100,$
                               value='Setup File',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Saved Setup File Button'})
    
        path = coyote_field2(row,$
                             TITLE=':',$
                             VALUE = self.transitpath1,$
                             UVALUE='Transit Setup Field',$
                             XSIZE=65,$
                             TEXTID=textid)
        self.transitfile1_fld = [path,textid]
        
        clear = widget_button(row,$
                              FONT=self.buttonfont,$
                              VALUE='Clear',$
                              uvalue={object:self, $
                                      method:'ButtonEvent', $
                                      value:'Clear Setup Path'})
        
        self.buttons[1] =  widget_button(row,$
                                font=self.buttonfont,$
                                xsize=110,$
                                value='Load Setup',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Load Setup Button'},sensitive=0)
      ;  widget_control,   (*self.bases)[1],/map

        if 0 then begin
        row = widget_base( (*self.bases)[1] ,$
                           /ROW,$
                           /Base_align_center,frame=1)
        button = widget_button(row,$
                               font=self.buttonfont,$
                               xsize=100,$
                               value='Save Setup',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Save Current Setup Button'})
        
        path = coyote_field2(row,$
                             TITLE=':',$
                             VALUE = self.transitpath2,$
                             UVALUE='Transit Save Setup Field',$
                             XSIZE=65,$
                             TEXTID=textid)
        self.transitfile2_fld = [path,textid]
        
        clear = widget_button(row,$
                              FONT=self.buttonfont,$
                              VALUE='Clear',$
                              uvalue={object:self, $
                                      method:'ButtonEvent', $
                                      value:'Clear Data Path'})
        
        button =  widget_button(row,$
                                font=self.buttonfont,$
                                xsize=110,$
                                value='Save Setup',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Save Setup'})
     endif              
     end
     'multi': begin
        tframe = 0
        if self.num_transits ne 0 then sens=1 else sens = 0
        (*self.bases)[2] = widget_base( (*self.bases)[0],$
                                        /column,$
                                        MAP=0,sensitive=sens,frame=tframe)
        
        base = widget_base((*self.bases)[2],/row,frame=tframe)
        col = widget_base(base,/column,frame=1)
        lbl = widget_label(col,value='General')

        if self.num_transits lt 2 then sens1 = 0 else sens1 = 1
        slider = widget_slider(col,minimum=1,maximum=max([2,self.num_transits]), $
                               value=max([1,self.active_transit]), $
                               uvalue={object:self, method:'ButtonEvent', $
                                       value:'ActiveTransit'},/drag,title='Active Transit',sensitive=sens1)

        slider = cw_fslider(col,title='Plot scaling',min=0,max=.02,value=self.diff,$
                            uvalue={object:self, method:'buttonevent',value:'diffset'},drag=1,$
                            /double, /edit,xsize=150)
        
        col=widget_base(base,/column,frame=1,/base_align_center)
        
        self.settings[0] = widget_label(col,value='settings:',/dynamic_resize)
        
        
        if self.num_transits gt 0 then widget_control,self.settings[0],set_value=string(self.active_transit,format='(i2.2)')+": "+((*self.transits)[self.active_transit-1]->get()).fname+' settings:'
        
        row = widget_base(col,/row)
        col2 = widget_base(row,frame=1,/col)
      options = ['None','Rebin to "Input Integration"']
        self.settings[1] = cw_bgroup(col2,$
                                     options,$
                                     /column,$
                                     LABEL_top='Long Integrations?',$
                                     /RETURN_NAME,$
                                     /NO_RELEASE,$
                                     set_value =0,$
                                     UVALUE={object:self, method:'ButtonEvent', value:'RebinTypeA'},$      
                                     FONT=self.buttonfont,$
                                     /EXCLUSIVE)
                                ;self.setup_rebin = 0
        options = 0L
        if self.num_transits gt 0 then widget_control,self.settings[1],set_value=((*self.transits)[self.active_transit-1]->get()).rebin        
        col2 = widget_base(row,frame=1,/col)
        
        lbl = widget_label(col2,value='Input Integration')
        
      ;  row = widget_base(,/row)
        self.settings[2] = coyote_field2(col2,$
                            LABELFONT=self.buttonfont,$
                            FIELDFONT=self.textfont,$
                            TITLE='Minutes:',$
                            UVALUE={object:self, method:'ButtonEvent', $
                                    value:'settint2',set:0},$
                                         decimal=4, digits=8, $
                                         VALUE=self.smooth_val[0],$
                                         XSIZE=10,/doubleValue,event_pro='TAP_event',$
                            textid=textid)    
        
        
     ;   lbl = widget_label(self.setup_smooth,value='Samples per Integration')
        
        self.settings[3] = coyote_field2(col2,$
                            LABELFONT=self.buttonfont,$
                            FIELDFONT=self.textfont,$
                            TITLE='N_samp:',$
                            UVALUE={object:self, method:'ButtonEvent', $
                                    value:'settint2',set:1},$
                            decimal=4, digits=8, $
                            VALUE=round(self.smooth_val[1]),$
                            XSIZE=5,/integerValue,event_pro='TAP_event',$
                            textid=textid)    
        
        
        button = widget_button(col,xsize=150,value='Delete Active Transit',$
                               uvalue={object:self,$
                                       method:'ButtonEvent',$
                                       value:'Delete Active'})

        
        col=widget_base(base,/column,frame=1)    
        

        

lbl = widget_label(col,value='Inter-Transit Settings:')

        button =  widget_button(col,$
                                font=self.buttonfont,$
                                xsize=60,$
                                value='Set Links',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Setup Cross LC Locks'})    
             
        
     end
     'oldmulti': begin
        if self.num_transits ne 0 then sens=1 else sens = 0
        
   ;     if (*self.bases)[2] ne 0 then widget_control,(*self.bases)[2],/destroy
                                ;    if (*self.bases)[2] eq 0 then begin
        (*self.bases)[2] = widget_base( (*self.bases)[0],$
                                        /column,$
                                        MAP=0,sensitive=sens)
        
        row = widget_base( (*self.bases)[2] ,$
                           /column)
        
        button =  widget_button(row,$
                                font=self.buttonfont,$
                                xsize=200,$
                                value='Setup Cross LC Locks',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Setup Cross LC Locks'})
        if self.num_transits gt 0 then begin
           transits='01: '+((*self.transits)[0]->get()).fname
           for i=1,n_elements(*self.transits)-1,1 do transits = [transits,string(i+1,format='(i2.2)')+': '+((*self.transits)[i]->get()).fname]
           bg = cw_bgroup(row,$
                          transits,$
                          ROW=n_elements(transits)/4d0,$
                          LABEL_LEF='Set Active Transit:',$
                          /RETURN_NAME,$
                          /NO_RELEASE,$
                          UVALUE={object:self, method:'ButtonEvent', value:'ActiveTransit'},$      
                          FONT=self.buttonfont,$
                          /EXCLUSIVE,/scroll,$
                          x_scroll_size=500, y_scroll_size=50,$ ;,ysize=20,$
                          
                          SET_VALUE=self.active_transit-1)
           transits = 0L
        endif
                                ; endif
     end
     'fit': begin
        if self.num_transits ne 0 then sens=1 else sens = 0
     ;   if (*self.bases)[3] ne 0 then widget_control,(*self.bases)[3],/destroy
        (*self.bases)[3] = widget_base( (*self.bases)[0],$
                                        /column,$
                                        MAP=0,sensitive=sens)
        
        
        row1_base = widget_base(  (*self.bases)[3] ,$
                                  /ROW)
        
        col1_base = widget_base(row1_base,$
                                /COLUMN,$
                                /BASE_ALIGN_LEFT,$
                                FRAME=1)
        
        label = widget_label(col1_base,$
                             value='Manual Parameter Adjustments and Setup',$
                             font=self.buttonfont,$
                             /align_left)
        
;;; COL 1: manual adjust!        
        
        
        button = widget_button(col1_base,$
                               font  = self.buttonfont,$
                               value = 'Manual Parameter Adjustment',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value: 'Manual Parameter Adjustment'},$
                               /align_center)
        
        
        
        button = widget_button(col1_base,$
                               font=self.buttonfont,$
                               value='Adjust Limits and Locks',$
                               uvalue= {object:self,$
                                        method: 'ButtonEvent',$
                                        value: 'Adjust Limits and Locks'},$
                               /align_center)
        
        
        
        
;;; COL1 END
        
                                ; col2_base = widget_base(row1_base,$
                                ;                         /COLUMN,$
                                ;                         /BASE_ALIGN_center,$
                                ;                         FRAME=2)
                                ; 
                                ; label = widget_label(col2_base,$
                                ;                      value='Levenberg-Markwart',$
                                ;                      font=self.buttonfont,$
                                ;                      /align_left)
  ;;; COL 2: automatic LM fit...   
        
        
  col3_base = widget_base(row1_base,$
                          /column,$
                          /BASE_ALIGN_left,$
                          FRAME=2)
  
  label = widget_label(col3_base,$
                       value='Markov Chain Monte Carlo',$
                       font=self.buttonfont,$
                       /align_center)
  ;;; col 3 MCMC
        

  button = widget_button(col3_base,$
                         value = 'MCMC Parameters',$
                         uvalue={object:self, $
                                 method:'ButtonEvent', $
                                 value: 'MCMC Parameters'},$
                         /align_center,$
                         xsize=110.,$
                                ; ysize=30.,$
                         sensitive=1)
  
  button = widget_button(col3_base,$
                         value = 'Gaussian Priors',$
                         uvalue={object:self, $
                                 method:'ButtonEvent', $
                                 value: 'Gaussian Priors'},$
                         /align_center,$
                         xsize=110.,$
                                ; ysize=30.,$
                         sensitive=1) 

  button = widget_button(col3_base,$
                         value = 'Execute Chain',$
                         uvalue={object:self, $
                                 method:'ButtonEvent', $
                                 value: 'Execute MCMC'},$
                         /align_center,$
                         xsize=110.,$
                                ; ysize=30.,$
                         sensitive=1)
  
        
        
     end
     'inference': begin
        (*self.bases)[4] = widget_base( (*self.bases)[0],$
                                        /column,$
                                        MAP=0)
        row1_base = widget_base(  (*self.bases)[4] ,$
                                  /column)
        
        
        col1_base = widget_base(row1_base,$
                                /COLUMN,$
                                /BASE_ALIGN_LEFT,$
                                FRAME=1)
        label = widget_label(col1_base,value='1) Select options for output:'$
                             ,font=self.buttonfont, /align_left)
        row = widget_base(col1_base,/row,frame=0,/base_align_center)
        radio = widget_base(row,column=1,/nonexclusive)
        button = widget_button(radio,value='Create .eps Plots',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'create_plots'})
        widget_control,button,set_button=self.create_plots
        
       row = widget_base(col1_base,/row,frame=0,/base_align_center)
        radio = widget_base(row,column=1,/nonexclusive)
        button = widget_button(radio,value='Output MCMC chains to ascii',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'create_mcmcascii'})
        widget_control,button,set_button=self.create_mcmcascii
        




        col1_base = widget_base(row1_base,$
                                /COLUMN,$
                                /BASE_ALIGN_LEFT,$
                                FRAME=1)
        label = widget_label(col1_base,value='2) Load a compatible TAP_setup.idlsav file'$
                             ,font=self.buttonfont, /align_left)
        row = widget_base(col1_base,/row,frame=0,/base_align_center)
        
        button = widget_button(row,$
                               font=self.buttonfont,$
                               xsize=150,$
                               value='MCMC Save File',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'MCMC File Button'})
        
        if self.mcmc_fld[1] ne 0 then widget_control,self.mcmc_fld[1],get_value=pval else pval= ''
        path = coyote_field2(row,$
                             labelfont=self.buttonfont,$
                             FIELDFONT=self.textfont,$
                             TITLE=':',$
                             VALUE = pval,$
                             UVALUE='MCMC File Field',$
                             XSIZE=70,$
                             TEXTID=textid)
        self.mcmc_fld = [path,textid]
       
        pval = 0L
        
        self.buttons[0] =  widget_button(row, font=self.buttonfont, xsize=40, value='Load',$
                                         sensitive=0, uvalue={object:self, method:'ButtonEvent', $
                                                              index: 0, value:'Load MCMC Button'})
        
        clear = widget_button(row,FONT=self.buttonfont,VALUE='Clear', xsize=40,$
                              uvalue={object:self, $
                                      method:'ButtonEvent', $
                                      value:'Clear MCMC Path'})
        
        
        
        
        
        
     end
     'cross lc links': begin
        (*self.extra_windows)[0] = widget_base(title='MCMC Multi Chain Parameter Sets',/column,frame=3,uname='links')
        XManager, 'TAP' $
                  , (*self.extra_windows)[0] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[0] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:0, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[0],frame=1,/row)
        for i=0,self.num_transits-1,1 do begin
           column = widget_base(work_base,frame=1,/column)
           values = round(((*self.transits)[i]->get()).params.set)
         ;  help,values
           label=widget_label(column,font=self.buttonfont $
                              , Value=((*self.transits)[i]->get()).fname,/align_center)   
           for j=0,n_elements(((*self.transits)[i]->get()).params)-1,1 do begin
              row = widget_base(column,frame=0,/row,/align_center)
              if i eq 0 then begin
                 button = widget_button(row,value='Lock All',$
                                       uvalue={object:self,  method:'ButtonEvent', value:'xlock_lockall',param:j})
                 button = widget_button(row,value='Free',$
                                       uvalue={object:self,  method:'ButtonEvent', value:'xlock_freeall',param:j})
                 label=widget_label(row,FONT=self.buttonfont,Value=(*self.label_indices)[j+1],xsize=80,/align_center)
                 
              endif
              self.fld[i,j] = coyote_field2(row,$
                                            LABELFONT=self.buttonfont,$
                                            FIELDFONT=self.textfont,$
                                            TITLE='',$
                                            /integervalue, $
                                            UVALUE={object:self, method:'adjustxlock_event',value:i,$
                                                    type: 'val', param:j $
                                                   },$
                                            VALUE=values[j],$
                                            /positive,$
                                            XSIZE=2,$
                                            scr_ysize=30,$
                                            event_pro='TAP_event',$
                                            textid=textid)
           endfor
           values = 0L
        endfor
                     
     end
     'gaussianpriors': begin
        (*self.extra_windows)[4] = widget_base(title='MCMC Gaussian Priors',/column,frame=3,uname='mcmcgaussianpriors')
        XManager, 'TAP' $    
                  , (*self.extra_windows)[4] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[4] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:4, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[4],frame=0,/column,/base_align_center)

        transit = (*self.transits)[self.active_transit-1]->get()
        
        for i=0,n_elements(transit.params)-1,1 do begin
           row = widget_base(work_base,frame=0,/row)
           label = widget_label(row,font=self.buttonfont,value=transit.params[i].param,xsize=80,/align_right)
           fld = coyote_field2(row,$
                               LABELFONT=self.buttonfont,$
                               FIELDFONT=self.textfont,$
                               TITLE='Value:',$
                               UVALUE={object:self, method:'adjustLL_event',value:i,$
                                       type: 'prior_val'$
                                      },$
                               VALUE=((transit.params[i]).prior)[1],$
                               XSIZE=15,$
                               /doubleValue,event_pro='TAP_event',$
                               textid=textid)
           fld = coyote_field2(row,$
                               LABELFONT=self.buttonfont,$
                               FIELDFONT=self.textfont,$
                               TITLE='Value:',$
                               UVALUE={object:self, method:'adjustLL_event',value:i,$
                                       type: 'prior_sig'$
                                      },$
                               VALUE=((transit.params[i]).prior)[2],$
                               XSIZE=15,$
                               /doubleValue,event_pro='TAP_event',$
                               textid=textid)
           radio = widget_base(row,column=1,/nonexclusive)
           button = widget_button(radio,value='Enable',$
                                ;uname = ,$ 
                                  uvalue={object:self,$
                                          method:'adjustLL_event',$
                                          value: i,$
                                          type: 'prior_penalize'$
                                         })
           widget_control,button,set_button= (transit.params[i].prior)[0]
          
        endfor
       ; transit = 0L
     end
     'mcmcparams': begin
        (*self.extra_windows)[3] = widget_base(title='MCMC Chain Parameters',/column,frame=3,uname='mcmcparams')
        XManager, 'TAP' $    
                  , (*self.extra_windows)[3] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[3] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:3, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[3],frame=0,/column,/base_align_center)

        base      = widget_base(work_base,frame=1,/column)
        label     = widget_label(base,FONT=self.buttonfont,Value='Number of Chains:',/align_left)
        
        fld = coyote_field2(base,$
                            LABELFONT=self.buttonfont,$
                            FIELDFONT=self.textfont,$
                            TITLE=' ',$
                            UVALUE={object:self, method:'adjustLL_event',value:'numchain',$
                                    type: 'mcmcsetup'$
                                   },$
                            VALUE=(((*self.transits)[self.active_transit-1])->get()).mcmc_params[0],$
                            XSIZE=10,$
                            
                                ;  format='(G10.2)',$
                            /doubleValue, $
                            ;/integervalue,$
                            decimal = 0,$
                            event_pro='TAP_event',$
                            textid=textid)
        
        
        ;; chain length
         ;;;  'chainmin':  self.mcmc_chainlength[0] = *event.value
         ;;;  'chainmax':  self.mcmc_chainlength[1] = *event.value
        
       ; base      = widget_base(work_base,frame=1,/column)
  label     = widget_label(base,FONT=self.buttonfont,Value='Chain Length:',/align_left)
  
  fld = coyote_field2(base,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Minimum Links:',$
                      UVALUE={object:self, method:'adjustLL_event',value:'chainmin',$
                              type: 'mcmcsetup'$
                             },$
                      VALUE=(((*self.transits)[self.active_transit-1])->get()).mcmc_params[1],$
                      XSIZE=10,$
                      decimal=0,$
                                ;  format='(G10.2)',$
                      /doubleValue,event_pro='TAP_event',$
                      textid=textid)
  
  if 0 then begin
  fld = coyote_field2(base,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Maximum Links:',$
                      UVALUE={object:self, method:'adjustLL_event',value:'chainmax',$
                              type: 'mcmcsetup'$
                             },$
                      VALUE=(((*self.transits)[self.active_transit-1])->get()).mcmc_params[2],$
                      XSIZE=10,$
                      decimal=0,$
                                ;  format='(G10.2)',$
                      /doubleValue,event_pro='TAP_event',$
                      textid=textid)

endif
     end
     'adjustparams': begin
        (*self.extra_windows)[1] = widget_base(title='Manual Parameter Adjustment',/column,frame=3,uname='links')
        XManager, 'TAP' $
                  , (*self.extra_windows)[1] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[1] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:1, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[1],frame=0,/column,/base_align_center)
        
        
        label = widget_label(work_base,value='LC: '+((*self.transits)[self.active_transit-1]->get()).fname)



        row = widget_base(work_base,/column,/align_center,frame=1)
        
        label = widget_label(row,$
                             value='System Parameters',$
                             font=self.buttonfont,$
                             /align_left)
        
        col_params = widget_base(row,$
                                 frame=0,$
                                 column=3,$
                                 /align_left)
        for i=0,4,1 do begin
        ;   help,/heap
           (*(self.slider))[i].id =$
              CW_Fslider(col_params, $
                         title=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         min=((*self.transits)[self.active_transit-1]->get()).params[i].limits[0], $
                         max=((*self.transits)[self.active_transit-1]->get()).params[i].limits[1],$
                         value = ((*self.transits)[self.active_transit-1]->get()).params[i].value,$
                         format='(G20.15)',/double, /edit,$
                         uname=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         uvalue={object:self, method:'AdjustSlider', value:((*self.transits)[self.active_transit-1]->get()).params[i].param},$
                         drag=1)
         ;  help,/heap
         ;  stop
        endfor
        for i=7,8,1 do begin
           (*(self.slider))[i].id =$
              CW_Fslider(col_params, $
                         title=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         min=((*self.transits)[self.active_transit-1]->get()).params[i].limits[0], $
                         max=((*self.transits)[self.active_transit-1]->get()).params[i].limits[1],$
                         value = ((*self.transits)[self.active_transit-1]->get()).params[i].value,$
                         format='(G20.15)',/double, /edit,$
                         uname=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         uvalue={object:self, method:'AdjustSlider', value:((*self.transits)[self.active_transit-1]->get()).params[i].param},$
                         drag=1)
           endfor
        row = widget_base(work_base,/column,/align_left,frame=1)
        
        label = widget_label(row,$
                             value='Quadratic Limb Darkening',$
                             font=self.buttonfont,$
                             /align_left)
        
        col_params = widget_base(row,$
                                 frame=0,$
                                 column=3,$
                                 /align_left)
        
           for i=5,6,1 do begin
           (*(self.slider))[i].id =$
              CW_Fslider(col_params, $
                         title=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         min=((*self.transits)[self.active_transit-1]->get()).params[i].limits[0], $
                         max=((*self.transits)[self.active_transit-1]->get()).params[i].limits[1],$
                         value = ((*self.transits)[self.active_transit-1]->get()).params[i].value,$
                         format='(G20.15)',/double, /edit,$
                         uname=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         uvalue={object:self, method:'AdjustSlider', value:((*self.transits)[self.active_transit-1]->get()).params[i].param},$
                         drag=1)
           endfor
           row = widget_base(work_base,/column,/align_left,frame=1)
        
        label = widget_label(row,$
                             value='Data Specific Parameters',$
                             font=self.buttonfont,$
                             /align_left)
        
        col_params = widget_base(row,$
                                 frame=0,$
                                 row=2,$
                                 /align_left)
        
           for i=9,12,1 do begin
           (*(self.slider))[i].id =$
              CW_Fslider(col_params, $
                         title=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         min=((*self.transits)[self.active_transit-1]->get()).params[i].limits[0], $
                         max=((*self.transits)[self.active_transit-1]->get()).params[i].limits[1],$
                         value = ((*self.transits)[self.active_transit-1]->get()).params[i].value,$
                         format='(G20.15)',/double, /edit,$
                         uname=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         uvalue={object:self, method:'AdjustSlider', value:((*self.transits)[self.active_transit-1]->get()).params[i].param},$
                         drag=1)
        endfor
        end
     'adjustlimlocks' : begin
        ;;   'Adjust Limits and Locks': begin
        ;;      self->setup,'adjustlimlocks'
        ;;      centertlb,(*self.extra_windows)[2]
        ;;      widget_control,(*self.extra_windows)[2],/realize
        ;;   end
        (*self.extra_windows)[2] = widget_base(title='MCMC Limits and Locks' , $
                                               /column,frame=3,uname='links')
        XManager, 'TAP' $
                  , (*self.extra_windows)[2] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[2] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:2, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[2],frame=0,/column)
        label = widget_label(work_base,value='LC: ' + $
                             ((*self.transits)[self.active_transit-1]->get()).fname,/align_center)
        col = widget_base(work_base,/column,/align_center,frame=0)

        transit = (*self.transits)[self.active_transit-1]->get()
        
        label = widget_label(col,value='System Parameters',font=self.buttonfont,/align_left)
        col_params = widget_base(col,frame=1,/column,/align_left)
        for i=0,4,1 do self->llrow,transit,i,col_params
        for i=7,8,1 do self->llrow,transit,i,col_params
        label = widget_label(col,value='Quadratic Limb Darkening',font=self.buttonfont,/align_left)
        col_params = widget_base(col,frame=1,/column,/align_left)
        for i=5,6,1 do   self->llrow,transit,i,col_params
        label = widget_label(col,value='Data Specific Parameters',font=self.buttonfont,/align_left)
        col_params = widget_base(col,frame=1,/column,/align_left)
        for i=9,12,1 do self->llrow,transit,i,col_params
                           
        transit = 0L
        
        
     end
     else: print,'unknown SETUP: '+set_this
  endcase
end


pro TAP::AdjustLL_event,event
  widget_control, event.id, GET_UVALUE= uvalue
 ; print,*event.value
 ; print,uvalue.type
  transit = (*self.transits)[self.active_transit-1]->get()
  case uvalue.type of
     'min': begin
        if n_elements(*event.value) gt 0 then begin
           transit.params[uvalue.value].limits[0] = *event.value
           for i=0,self.num_transits-1,1 do begin
              if i ne self.active_transit-1 then begin
                 trans = (*self.transits)[i]->get()
                 if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                    trans.params[uvalue.value].limits[0] = *event.value
                    (*self.transits)[i]->set,trans
                 endif
                 trans = 0L
              endif
           endfor
        endif
     end
     'max': begin
        if n_elements(*event.value) gt 0 then begin
           transit.params[uvalue.value].limits[1] = *event.value
           for i=0,self.num_transits-1,1 do begin
              if i ne self.active_transit-1 then begin
                 trans = (*self.transits)[i]->get()
                 if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                    trans.params[uvalue.value].limits[1] = *event.value
                    (*self.transits)[i]->set,trans
                 endif
                 trans = 0L
              endif
           endfor
        endif
     end
     'lock': begin
        if transit.params[uvalue.value].fixed eq 0 then $
           transit.params[uvalue.value].fixed = 1 else $
              transit.params[uvalue.value].fixed = 0
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].fixed = transit.params[uvalue.value].fixed
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'limit': begin
        if transit.params[uvalue.value].limited[0] eq 0 then $
           transit.params[uvalue.value].limited = [1,1] else $
              transit.params[uvalue.value].limited = [0,0]
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].limited = transit.params[uvalue.value].limited
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'mcmcacceptrate': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        transit.params[uvalue.value].accept = (*event.value)
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].accept = transit.params[uvalue.value].accept
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'mcmcsetup': begin
        transit = 0L
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        case uvalue.value of
           'numchain':  begin
              for i=0,n_elements(*self.transits)-1,1 do begin
                 tran = (*self.transits)[i]->get()
                 tran.mcmc_params[0] = *event.value
                 (*self.transits)[i]->set,tran
                 if self.active_transit-1 eq i then transit = tran
                 tran = 0L
              endfor
           end
           'chainmin':   begin
              for i=0,n_elements(*self.transits)-1,1 do begin
                 tran = (*self.transits)[i]->get()
                 tran.mcmc_params[1] = *event.value
                 (*self.transits)[i]->set,tran
                 if self.active_transit-1 eq i then transit = tran
                 tran = 0L
              endfor
           end
           'chainmax':    begin
              for i=0,n_elements(*self.transits)-1,1 do begin
                 tran = (*self.transits)[i]->get()
                 tran.mcmc_params[2] = *event.value
                 (*self.transits)[i]->set,tran
                 if self.active_transit-1 eq i then transit = tran
                 tran = 0L
              endfor
           end
        endcase
     end
     'prior_val': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        transit.params[uvalue.value].prior[1] = (*event.value)
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].prior = transit.params[uvalue.value].prior
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'prior_sig': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        transit.params[uvalue.value].prior[2] = (*event.value)
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].prior = transit.params[uvalue.value].prior
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'prior_penalize': begin
        if transit.params[uvalue.value].prior[0] eq 0 then $
           transit.params[uvalue.value].prior[0] = 1 else $
              transit.params[uvalue.value].prior[0] = 0
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].prior = transit.params[uvalue.value].prior
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     else: begin
        print,'UNKNOWN adjust LL event'+uvalue.type
        stop
     end
  endcase
  (*self.transits)[self.active_transit-1]->set,transit
  transit=0L
                                ; endif
end


pro TAP::LLrow,transit,i,col_params
  row = widget_base(col_params,frame=0,/row,/base_align_left)
  label=widget_label(row,FONT=self.buttonfont, $
                     Value=transit.params[i].param+':', $
                     xsize=90,/align_right)
  fld = coyote_field2(row,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Min',$
                      UVALUE={object:self, method:'AdjustLL_event', $
                              value:i, type: 'min'},$
                      decimal=10, digits=20, $
                      VALUE=transit.params[i].limits[0],$
                      XSIZE=15,/doubleValue,event_pro='TAP_event',$
                      textid=textid)    
  fld = coyote_field2(row,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Max',$
                      UVALUE={object:self, method:'AdjustLL_event', $
                              value:i, type: 'max'},$
                      decimal=10, digits=20, $
                      VALUE=transit.params[i].limits[1],$
                      XSIZE=15,/doubleValue,event_pro='TAP_event',$
                      textid=textid)    
  radio = widget_base(row,column=1,/nonexclusive)
  button = widget_button(radio,value='Lock',$ 
                         uvalue={object:self, method:'AdjustLL_event',$
                                 value: i, type: 'lock'})
  widget_control,button,set_button=transit.params[i].fixed
  fld = coyote_field2(row,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='MCMC Accept Rate:',$
                      UVALUE={object:self, method:'AdjustLL_event', $
                              value:i, type: 'mcmcacceptrate' },$
                      VALUE=transit.params[i].accept,$
                      XSIZE=5,$
                      decimal=2,$
                      digits=3,$
                                ;  format='(G10.2)',$
                      /doubleValue,event_pro='TAP_event',$
                      textid=textid)
  radio = widget_base(row,column=1,/nonexclusive)
  button = widget_button(radio,value='Apply Limits to Fitting.',$
                         uvalue={object:self, method:'AdjustLL_event',$
                                 value: i, type: 'limit'})
  widget_control,button,set_button=transit.params[i].limited[0]
end

pro tap::widget_setup
  self.tap_base = widget_base(title='Transit Analysis Package '+self.version,/column)
  
  XManager, 'TAP' $
            , self.tap_base $
            , /no_block $
            , cleanup = 'TAP_cleanup'
  

  quit_button = widget_button(self.tap_base,$
                              font = self.buttonfont,$
                              value = 'Quit',$
                              uvalue={object:self, method:'Quit'})
 
  message = '('+ curr_date(format='hh:mm:ss yyyymmdd') +') TAP Tools '+self.version
  self.message_window = widget_text(self.tap_base, $
                                    font = self.textfont, $
                                    value = message, $
                                    /scroll, $
                                    ysize=4)
  
  main_base = widget_base(self.tap_base,$
                          /row,$
                          frame=5, event_pro='TAP_AllEvents')
  



  col1_base = widget_base(main_base,$
                          /COLUMN,$
                          /BASE_ALIGN_RIGHT,$
                          /FRAME)

  ysize=22.
  xsize=320.
  for i=0,n_elements(*self.label_indices)-1,1 do begin
     (*self.label)[i] = widget_label(col1_base,$
                                     font=self.textfont, $
                                     value = ' ', $
                                     /align_left,$
                                     ysize = ysize,$
                                     xsize = xsize)
  endfor
  self.numcol = 1
  self->PrepCol

  col2_base = widget_base(main_base,$
                          /column,$
                          /base_align_left,$
                          /frame)

  (*self.plot_windows)[0].x = 400d0
  (*self.plot_windows)[0].y = 350d0
  
  (*self.plot_windows)[0].window = widget_draw(col2_base $
                                               , xsize=(*self.plot_windows)[0].x $
                                               , ysize=(*self.plot_windows)[0].y $
                                               , uvalue='Plot Window 1')
  
 

  ;; menu
  menubar = widget_base(self.tap_base,$
                        /row)
  
  row = widget_base(menubar,$
                    /row,$
                    /toolbar,$
                    /exclusive,$
                    /base_align_center)
  
  for i=0,n_elements(*self.menus)-1 do begin
     button = widget_button(row,$
                            value=' '+(*self.menus)[i]+' ',$
                            uvalue={object:self, method:'menumap', value:(*self.menus)[i], type:'main'},$
                            /no_release,$
                            font=self.buttonfont)
     if i eq 0 then widget_control, button, /SET_BUTTON
  endfor

  ;; workspace:
  (*self.bases)[0]  = widget_base(self.tap_base,$
                                  frame=5)
  

  self->setup,'load'
  self->setup,'multi'
  self->setup,'fit'
  self->setup,'inference'
end

pro TAP::lcplot,event
  !p.multi = [0,1,1]
  device,decomposed=0
  
  wset,(*self.plot_windows)[0].pix_window
  if self.num_transits ne 0 then begin
     (*self.plot_windows)[0].xrange = [1d4,-1d4]
     for i=0,n_elements(*self.transits)-1,1 do $
        (*self.plot_windows)[0].xrange = [$
        min([(*self.plot_windows)[0].xrange[0],min(((*self.transits)[i]->get()).transit_tday- $
                                                   (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit')  $
                                                                                              eq 1)].value-(((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))]),$
        max([(*self.plot_windows)[0].xrange[1],$
             max(((*self.transits)[i]->get()).transit_tday-$
                 (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value $
                  - (((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))])]
     
     
     tot = self.num_transits-1
     lci = (*self.transits)[0]->get()
     lcf = (*self.transits)[self.num_transits-1]->get()
     diff = self.diff
     if diff eq 0 then diff = 10d0*max([lci.residuals,lcf.residuals])
     
     
     yrange = [min(lci.transit_flux)-tot*(diff/2d0)-diff,max(lcf.transit_flux)+((tot+.5)*diff)]
     
     title = string(self.active_transit,format='(i2.2)')+": "+((*self.transits)[self.active_transit-1]->get()).fname
     
     for i=0,self.num_transits-1,1 do begin
        lc = (*self.transits)[i]->get() 
                                ;   diff = 5d0*max(lc.residuals)
        midt = lc.params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
        if i eq 0 then $
           plot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,$
                yrange=yrange,color=(*self.colors).black,background=(*self.colors).white,$
                xrange=(*self.plot_windows)[0].xrange*24d0,/xs,/ys,xtitle='Hours from Mid Transit',$
                ytitle='Relative Flux + Constant',title=title,/nodata
        
        oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,color=lc.color
        
        oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(diff*i),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)*(i+1)),psym=8,symsize=.6,color=lc.color
        hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=2
        lc = 0L   
        midt = 0L
        
        
        if i+1 eq self.active_transit then begin
           plotsym,7,3.0,/fill,thick=5
           oplot,[(*self.plot_windows)[0].xrange[0]*24d0],[1d0+(diff*i)],psym=8,color=(*self.colors).blue,symsize=1
           plotsym,6,3.0,/fill,thick=5
           oplot,[(*self.plot_windows)[0].xrange[1]*24d0],[1d0+(diff*i)],psym=8,color=(*self.colors).blue,symsize=1
           plotsym,0,1.2,/fill
        endif
     endfor
     lci = 0L
     lcf = 0L
     tot = 0L
     yrange = 0L
  endif else plot,[0],[0],/nodata,color=(*self.colors).gray,background=(*self.colors).gray
  wset,(*self.plot_windows)[0].w_id
  device,copy=[0,0,(*self.plot_windows)[0].x,(*self.plot_windows)[0].y,0,0,(*self.plot_windows)[0].pix_window]
  
end

pro TAP::quit,event
  self->destroy
end

pro tap::loadmcmc
  self.message = 'Restoring Savefile...'
  self->message  
  tap_state = 0L
  widget_control,self.mcmc_fld[1],get_value=path
 ; obj_destroy,self.restoreSaved
  self.restoreSaved = obj_new('IDL_Savefile',path)
;  restore,path
 ; stop
  
  if strcmp(strmid(path,strlen(path)-16,16),'TAP_setup.idlsav') then self.plot_dir = strmid(path,0,strlen(path)-16) else $
     if n_elements(strsplit(path,'/')) gt 1 then self.plot_dir = '/'+strjoin((strsplit(path,'/',/extract))[0:n_elements(strsplit(path,'/'))-2],'/')+'/' $
     else self.plot_dir = '\'+strjoin((strsplit(path,'\',/extract))[0:n_elements(strsplit(path,'\'))-2],'\')+'\'
  
  if n_elements(strsplit(path,'/')) le 1 then self.backslash = 1
  
  if where(self.restoresaved->names() eq 'TAP_STATE') eq -1 then begin
     self.message = 'Incompatible savefile.'
     self->message
  endif  else begin
     restore,path
     if (tap_state[0]->get()).mcmc_complete eq 0 then begin
        self.message = 'Savefile contains no MCMC chains.'
        self->message
        for i=0,n_elements(tap_state)-1,1 do begin
           st = TAP_state[i]->get()
           for j=0,n_elements(st.params)-1,1 do begin
              ptr_free,st.params[j].mcmc_chain
              ptr_free,st.params[j].refined_mcmc_chain
           endfor
           ptr_free,st.model_t,st.model_f,st.model_tfine,st.model_ffine,st.mcmc_files
           TAP_state[i]->destroy
           st = 0L
        endfor
        
     endif else begin
                                ;     self.message = 'Restoring Savefile...'
                                ;     self->message  
        
        ptr_free,self.transits
        self.transits = ptr_new(TAP_state)
        TAP_state = 0L
        self.active_transit = 1
        self.num_transits = n_elements(*self.transits)
        for i=0,self.num_transits-1,1 do begin
           transit = (*self.transits)[i]->get()   
           transit.params.jumpct = 0
           transit.params.jumptot= 0   
           (*self.transits)[i]->set,transit
           transit=0L
        endfor

        self.mcmc_complete =0
        transit = (*self.transits)[0]->get()   
        new = 0
        for i=0,n_elements((*transit.mcmc_files))-1,1 do begin
           if strcmp((*transit.mcmc_files)[i],'-1') ne 1 then begin
              file = (*transit.mcmc_files)[i]
              if self.backslash then strreplace,file,'/','\'
              self.mcmc_complete++
              self->addtofull,new,self->loadintocurr(self.plot_dir+file)
              file = 0
              new++
           endif
        endfor
        transit = 0L
        self->setup,'adjustparams'
        self->setup,'cross lc links'
 
        self->mcmc_inference

        self->setup,'multi'
        self->setup,'fit'
  
        
        self.diff = ((*self.transits)[0]->get()).diff
     endelse
  endelse
  obj_destroy,self.restoreSaved
  
end

pro tap::addtofull,cclear,go
  if cclear eq 0 then clearit=1 else clearit = 0
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()  
     for j=0,n_elements(transit.params)-1,1 do begin
        range = [.1*n_elements(*transit.params[j].mcmc_chain),n_elements(*transit.params[j].mcmc_chain)-1]
        if clearit then begin
           ptr_free,transit.params[j].refined_mcmc_chain
           if go then transit.params[j].refined_mcmc_chain = ptr_new((*transit.params[j].mcmc_chain)[range[0]:range[1]])
        endif else if go then *transit.params[j].refined_mcmc_chain = [*transit.params[j].refined_mcmc_chain,(*transit.params[j].mcmc_chain)[range[0]:range[1]]]  
     endfor
     
     if n_elements(strsplit(transit.fname,'/')) gt 1 then transit.fname = (strsplit(((strsplit(transit.fname,'/',/extract))[n_elements(strsplit(transit.fname,'/',/extract))-1]),'.',/extract))[0]  
     if n_elements(strsplit(transit.fname,'\')) gt 1 then transit.fname = (strsplit(((strsplit(transit.fname,'\',/extract))[n_elements(strsplit(transit.fname,'\',/extract))-1]),'.',/extract))[0]  

     
     (*self.transits)[i]->set,transit
     transit = 0L
  endfor
end

function tap::loadintocurr,filename
  if file_test(filename) then begin
     restore,filename
     for i=0,n_elements(tap_state)-1,1 do begin
        tt = (*self.transits)[i]->get()  
        st = TAP_state[i]->get()
        for j=0,n_elements(tt.params)-1,1 do begin
           ptr_free,tt.params[j].mcmc_chain
           tt.params[j].mcmc_chain = ptr_new(*st.params[j].mcmc_chain)
           tt.params[j].jumpct += st.params[j].jumpct
           tt.params[j].jumptot += st.params[j].jumptot
           ptr_free,st.params[j].mcmc_chain
           ptr_free,st.params[j].refined_mcmc_chain
           
        endfor
        (*self.transits)[i]->set,tt
        ptr_free,st.model_t,st.model_f,st.model_tfine,st.model_ffine,st.mcmc_files
        TAP_state[i]->destroy
        tt = 0L
        st = 0L
     endfor
     return,1
  endif else return,0
end

function TAP::INIT,input_ascii=input_ascii,smooth=smooth,period=period
  !p.font = 0;
  plotsym,0,1.2,/fill

  self.diff = 0
  self.buttonfont = ''
  self.textfont=''
  self.create_plots = 1
  self.create_ascii = 1
  self.create_mcmcascii=0
  
  self.smooth_val = [29.4244d0,4d0]

  device,decomposed=0
 
  ptr_free,self.colors,self.slider

  obj_destroy,self.restoreSaved

  self.colors = ptr_new(tap_colors())
  self.version = 'v2.104'
  self.slider = ptr_new(replicate({id: 0},22))

  ptr_free,self.plot_windows
  self.plot_windows = ptr_new(replicate({x: 0d, y: 0d, w_id: 0L, pw_id: 0L, window:0L, pix_window: 0L, xrange: [0d0,0d0], yrange: [0d0,0d0]},1))

  ptr_free,self.label_indices,self.label
  self.label_indices = ptr_new(['Parameters  ','Period:','Inclination:','a/R*:','Rp/R*:','Mid Transit:','Linear LD:','Quad LD:','Eccentricity:','Omega:','Airmass Y-int:','Airmass Slope:','Sigma Red:','Sigma White:'])
  self.label = ptr_new(lonarr(14));

  ptr_free,self.menus
  ptr_free,self.bases
  self.menus = ptr_new(['Load Transit','Manage Transits','Fit','MCMC Inference'])
  self.bases = ptr_new(lonarr(5))
  

  ptr_free,self.extra_windows
  self.extra_windows = ptr_new(lonarr(5))

  ptr_free,self.transits
  self.transits = ptr_new(-1)
  temp = obj_new('transit')
  temp->destroy
  
  if keyword_set(period) then self.init_period = period else self.init_period = 3d0
  
  if keyword_set(smooth) then begin
     self.setup_rebin=1
     self.smooth_val = smooth
  endif
  
  self->widget_setup
  ptr_free,self.mcmc_stuff
  
  self.mcmc_stuff = ptr_new()
 
  if keyword_set(input_ascii) then begin
     self.filetype = 'ASCII File'
     widget_control,self.transitfile_fld[1],set_value =input_ascii
  endif
  
  for i=1,n_elements(*self.bases)-1,1 do widget_control,(*self.bases)[i],MAP=0

  widget_control,(*self.bases)[1],MAP=1        
  
  Widget_Control,self.tap_base, set_UVALUE=self 
  self->start
  if keyword_set(input_ascii) then self->loadtransit,0
  return,1
end


pro tap__define
  struct = {tap,$
            $ ;; GENERAL STUFF
            version: '',$ 
            colors: ptr_new(), $
            buttonfont: '',$
            textfont: '',$
            restoreSaved: obj_new(),$
            $ ;; MAIN WIDGET stuff
            tap_base: 0L, $
            message_window: 0L, $
            message: '',$
            plot_windows: ptr_new(),$
            diff: 0d0,$
            create_plots: 0L,$
            create_ascii: 0L,$
            create_mcmcascii: 0L,$
            $ ;; EXTRA WIDGET WINDOWS
            extra_windows: ptr_new(),$
            $ ;; WIDGET pieces
            slider: ptr_new(),$
            ;buttons: ptr_new(),$
            label_indices: ptr_new(),$
            label: ptr_new(),$
            menus: ptr_new(),$
            bases: ptr_new(),$
            settings: lonarr(5),$
            fld: lonarr(50,13),$
            numcol: 0L,$
            $ ;; TRANSIT + MODEL stuff
            setup_rebin: 0L,$
            setup_smooth: 0L,$
            smooth_val: dblarr(2),$
            num_transits: 0d,$
            active_transit: 0d,$
          ;;  transit_colors: ptr_new(),$
            transits: ptr_new(),$
            parinfo: ptr_new(),$
            backslash: 0,$
            base_plot: '',$
            $ ;; PATHS:    
            init_period: 0d0,$
            base_path: '',$
            filetype: '',$
            transitpath:'',$
            transitpath1:'',$
            transitpath2:'',$
            transitfile_fld:[0,0], $
            transitfile1_fld:[0,0], $
            transitfile2_fld:[0,0], $
            mcmc_fld:[0,0],$
            plot_dir: '' ,$
            buttons: lonarr(5),$
            mcmc_stuff: ptr_new(),$
            mcmc_complete: 0,$
            $ ;; MCMC parameters
            $ ;; END
            exists: 0L $
}
end

pro tap,_REF_EXTRA=_extra
  tap = obj_new('tap',_EXTRA=_extra)
end





