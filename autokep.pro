
pro autokep_event,event
  Widget_Control, event.id, Get_UValue=info
  Call_Method, info.method, info.object, event
end





function aklast,array,num
  return,array[n_elements(array)-num]
end

function akgapless,array
  test = where(array eq 1)
 ; print,test
  maxgap=0
  if test[0] ne -1 then begin
     gapless=0
     for i=0,n_elements(test)-2,1 do begin
        if test[i+1] eq test[i]+1 then gapless++ else gapless=0
        maxgap = max([maxgap,gapless]) 
     endfor
  endif
  test = 0L
  gapless=0L
  return,maxgap
end

pro autokep::manualdetect,event,second=second
  found_transit = 0
  user_reject = 0
  sig = robust_sigma((*self.all_lc).flux)
  
  (*self.plot_windows)[1].bg = (*self.colors).gray
  (*self.plot_windows)[2].bg = (*self.colors).gray
  self->plot,''
  (*self.plot_windows)[1].bg = (*self.colors).white
  (*self.plot_windows)[2].bg = (*self.colors).white
  
  wset,(*self.plot_windows)[0].w_id
  
  plot,[0],[0],$
       xr=(*self.plot_windows)[0].xrange,/xs,$
       yrange=(*self.plot_windows)[0].yrange,/ys,psym=8,symsize=.6,color=(*self.colors).black,$
       background=(*self.plot_windows)[0].bg,xmargin=[8,3],ymargin=[3,1],xthick=1.2,ythick=1.2,$
       xtitle='Days since '+sigfig(floor(min((*self.all_lc).time)),8),xticks=5,yticks=2,/noerase,/nodata
  
  if 1-keyword_set(second) then self.message = 'Click before and after a transit in the white plot' else begin
     self.message = 'Click before and after a transit ONE epoch away'
  endelse
  self->message
  
  ;print,(*self.plot_windows)[0].xrange,(*self.plot_windows)[0].yrange
  
  cursor, xpos, ypos, /up, /data
  vline,xpos,color=(*self.colors).violet,thick=2
  ;print,xpos
  cursor, xpos1, ypos, /up, /data
  vline,xpos1,color=(*self.colors).violet,thick=2
  ;print,xpos1
  self.active_bound=[min([xpos,xpos1]),max([xpos,xpos1])]+floor(min((*self.all_lc).time))
  self->setactivelc,/noplot
  
  transit = where((*self.active_lc).flux le 1d0-3d0*sig)
  if transit[0] ne -1 then begin
     (*self.active_lc)[transit].transit = 1  
     if akgapless((*self.active_lc).transit) gt 3 then $
        if total(aklast((*self.active_lc).transit,5)) eq 0 then begin
        self.message = '...possible transit signature detected'
        self->message
        midpt = median((*self.active_lc)[transit].time)
        store = self.active_bound
        x =  akgapless((*self.active_lc).transit)
        self.t_duration = (x+4)*self.kep_dt
        self.duration = self.t_duration+(20*self.kep_dt)
        self.w_duration = self.duration
        
        self.active_bound[0] = midpt - self.duration/2d0
        self.active_bound[1] = midpt + self.duration/2d0
        self->setactivelc,/noplot
                                ;   transit = where((*self.active_lc).time gt midpt - (x/2d0)*self.kep_dt and $
                                ;               (*self.active_lc).time lt midpt + (x/2d0)*self.kep_dt)
        (*self.active_lc)[where((*self.active_lc).time ge midpt - self.t_duration/2 and $
                                (*self.active_lc).time le midpt + self.t_duration/2)].transit = 1
        (*self.active_lc)[where((*self.active_lc).time lt midpt - self.t_duration/2 or $
                                (*self.active_lc).time gt midpt + self.t_duration/2)].oot = 1
        (*self.plot_windows)[0].bg = (*self.colors).gray
        (*self.plot_windows)[2].bg = (*self.colors).gray
        self->plot,''
        (*self.plot_windows)[0].bg = (*self.colors).white
        (*self.plot_windows)[2].bg = (*self.colors).white
        self->updateslider
        
        if self->userquery('Is the desired transit centered in the white background plot?') then  found_transit = 1 else user_reject = 1
     endif
  endif
  transit = 0L
  
  if found_transit eq 0 and user_reject eq 0 then begin
     midpt = median((*self.active_lc).time)
     store = self.active_bound
     self.duration = max((*self.active_lc).time)-min((*self.active_lc).time)
     self.t_duration = .5*self.duration
     self.w_duration = self.duration
     
     self.active_bound[0] = midpt - self.duration/2d0
     self.active_bound[1] = midpt + self.duration/2d0
     self->setactivelc,/noplot
                                ;   transit = where((*self.active_lc).time gt midpt - (x/2d0)*self.kep_dt and $
                                ;               (*self.active_lc).time lt midpt + (x/2d0)*self.kep_dt)
     (*self.active_lc)[where((*self.active_lc).time ge midpt - self.t_duration/2 and $
                             (*self.active_lc).time le midpt + self.t_duration/2)].transit = 1
     (*self.active_lc)[where((*self.active_lc).time lt midpt - self.t_duration/2 or $
                             (*self.active_lc).time gt midpt + self.t_duration/2)].oot = 1
     (*self.plot_windows)[0].bg = (*self.colors).gray
     (*self.plot_windows)[2].bg = (*self.colors).gray
     self->plot,''
     (*self.plot_windows)[0].bg = (*self.colors).white
     (*self.plot_windows)[2].bg = (*self.colors).white
     self->updateslider
     
     
     if self->userquery('No automatic detection, add anyway?') then  found_transit = 1
  endif
  
  if found_transit then begin
     self->addtransit,midpt
     
     if self.num_transits eq 1 then begin
        self.message = 'Transit located.  Search for 2nd transit.'
        self->message
        widget_control,self.extra_bases[0],sensitive=1 
     endif else widget_control,self.extra_bases[1],sensitive=1
     
  endif
  
  self->plot,''
end

function autokep::userquery,message
 D = dialog_message(message,/question,/center)
 if strcmp(d,'Yes') then return,1
 return,0
end


pro autokep::autodetect,event,clear=clear

 ; self.active_trf = 1d0
  sig = robust_sigma((*self.all_lc).flux)
  
  if keyword_set(clear) then begin
     self.message = 'Auto Detect...'
     self->message
     for i=0,self.num_transits-1,1 do  (*self.transits)[i]->destroy ;obj_destroy,(*self.transits)[i]
     
     ptr_free,self.transits
     self.transits = ptr_new()
     self.num_transits = 0 
     ptr_free,self.midts
     self.midts = ptr_new(-1)
     
     self.active_bound=[min((*self.all_lc).time),0]
     self.active_bound[1] = self.active_bound[0]+1d0
  endif else begin
     self.message = 'Auto Detect 2nd...'
     self->message    
     self.active_bound[0] = self.active_bound[1] 
     self.active_bound[1] = self.active_bound[0]+1d0
  endelse
  
  self->setactivelc,/noplot
  
  no_t = 1
  nothing = 0
  found_transit = 0
  while no_t do begin
     if ((max((*self.all_lc).time) - self.active_bound[0]) lt 1) or (max((*self.all_lc).time) - self.active_bound[1] lt 0) then begin
        no_t = 0
        nothing = 1
     endif else begin
        if (*self.active_lc)[0].time ne -1 then begin
           transit = where((*self.active_lc).flux le 1d0-3d0*sig)
           if transit[0] ne -1 then begin
              (*self.active_lc)[transit].transit = 1  
              if akgapless((*self.active_lc).transit) gt 3 then $
                 if total(aklast((*self.active_lc).transit,5)) eq 0 then begin
                 self.message = '...possible transit signature detected'
                 self->message
                 midpt = median((*self.active_lc)[transit].time)
                 store = self.active_bound
                 x =  akgapless((*self.active_lc).transit)
                 self.t_duration = (x+4)*self.kep_dt
                 self.duration = self.t_duration+(20*self.kep_dt)
        self.w_duration = self.duration

                 self.active_bound[0] = midpt - self.duration/2d0
                 self.active_bound[1] = midpt + self.duration/2d0
                 self->setactivelc,/noplot
                                ;   transit = where((*self.active_lc).time gt midpt - (x/2d0)*self.kep_dt and $
                                ;               (*self.active_lc).time lt midpt + (x/2d0)*self.kep_dt)
                 (*self.active_lc)[where((*self.active_lc).time ge midpt - self.t_duration/2 and $
                                         (*self.active_lc).time le midpt + self.t_duration/2)].transit = 1
                 (*self.active_lc)[where((*self.active_lc).time lt midpt - self.t_duration/2 or $
                                         (*self.active_lc).time gt midpt + self.t_duration/2)].oot = 1
                 (*self.plot_windows)[0].bg = (*self.colors).gray
                 (*self.plot_windows)[2].bg = (*self.colors).gray
                 self->plot,''
                 (*self.plot_windows)[0].bg = (*self.colors).white
                 (*self.plot_windows)[2].bg = (*self.colors).white
                   self->updateslider

                 if self->userquery('Is the desired transit centered in the white background plot?') then begin
                    found_transit = 1
                    no_t = 0
                 endif else begin
                    self.message = '...continuing search'
                    self->message
                    self.active_bound[0] = store[1]
                    self.active_bound[1] = self.active_bound[0]+1d0
                    self->setactivelc,/noplot
                 endelse
                 self->plot,'active_lc'
                 
              endif
              transit=0L
           endif
           self->plot,'active_lc'    
        endif
       
        if found_transit then begin ;; if transit found
           self->addtransit,midpt
           
           if self.num_transits eq 1 then begin
              self.message = 'Transit located.  Search for 2nd transit.'
              self->message
              widget_control,self.extra_bases[0],sensitive=1 
           endif else widget_control,self.extra_bases[1],sensitive=1
           
        endif else begin ;; no transit found
           self.active_bound[1]+=(self.kep_dt*5d0)
           self->setactivelc,/noplot
        endelse
     endelse
  endwhile
  
  self.active_trf=-1d0
  
end

pro autokep::setepoch
  ptr_free,self.epoch
  zmidt = max((*self.midts))
  tmm = mm((*self.all_lc).time)
 
  epocht = zmidt
  if self.period ne 0 then begin
     ct = epocht[0] - self.period
     if ct ge tmm[0] then $
        while ct ge tmm[0] do begin
        ;print,ct,tmm[0]
        epocht = [epocht, ct]
        ct-=self.period 
     endwhile
     
     ct = epocht[0] + self.period
     if ct le tmm[1] then $
        while ct le tmm[1] do begin
        epocht = [epocht, ct]
        ct+=self.period 
     endwhile
  endif
  self.epoch = ptr_new(epocht[sort(epocht)])
  
  epocht=0L
  ct = 0L
  zmidt= 0L
  tmm = 0L
end



pro autokep::plotepoch
  y = [!y.crange[0],!y.crange[1],!y.crange[1],!y.crange[0]]
 ; self.duration = self.t_duration*6d0
  zmidt = (*self.midts)[0]
  dur = self.duration
  
  self->setepoch
  epoch = *self.epoch-floor(min((*self.all_lc).time))
  for i=0,n_elements(*self.epoch)-1,1 do begin
     polyfill,[(epoch[i]-(dur/2)),(epoch[i]-(dur/2)),$
               (epoch[i]+(dur/2)),(epoch[i]+(dur/2))],y,color=(*self.colors).yellow
     xyouts,epoch[i],(*self.plot_windows)[0].yrange[0]+1d-1*((*self.plot_windows)[0].yrange[1]-(*self.plot_windows)[0].yrange[0]),string(i+1,format='(i0.0)'),charsize=1,charthick=1,color=(*self.colors).black,align=-1
  endfor
  y = 0L
  zmidt=0L
  dur=0L
  epoch = 0L
end

pro autokep::plotall,event
  for i=0,self.num_transits-1,1 do begin
     lc = (*self.transits)[i]->get() 
     oplot,[lc.lc.time-floor(min((*self.all_lc).time))],[lc.lc.flux],psym=8,$
           symsize=.8,color=(*self.colors).dkgreen
     lc = 0L
  endfor
end

pro autokep::plotwrite,hours=hours,midt=midt  
  if self.duration ne self.w_duration then begin
     if 1-keyword_set(hours) then modif = 1d0 else modif = 24d0
     if 1-keyword_set(midt) then subt = floor(min((*self.all_lc).time)) else subt = ((*self.transits)[self.active_transit]->get()).midt
     vline,[(((*self.transits)[self.active_transit]->get()).midt+(self.w_duration/2d0))-subt]*modif,color=(*self.colors).violet,thick=2
     vline,[(((*self.transits)[self.active_transit]->get()).midt-(self.w_duration/2d0))-subt]*modif,color=(*self.colors).violet,thick=2
  endif
end


pro autokep::plotactive,event,offset=offset,hours=hours,midt=midt
  if 1-keyword_set(offset) then offset=0
  if 1-keyword_set(hours) then modif = 1d0 else modif = 24d0
  if 1-keyword_set(midt) then subt = floor(min((*self.all_lc).time)) else subt = ((*self.transits)[self.active_transit]->get()).midt
  if (*self.active_lc)[0].time ne -1 then begin
     oplot,[(*self.active_lc).time-subt]*modif,[(*self.active_lc).flux]+offset,psym=8,$
           symsize=.6,color=(*self.colors).green
     transit = (where((*self.active_lc).transit))
     if transit[0] ne -1 then $
        oplot,[((*self.active_lc).time)[transit]-subt]*modif,$
              [((*self.active_lc).flux)[transit]]+offset,psym=8,symsize=.6,color=(*self.colors).red  
     transit = 0L
     oot = (where((*self.active_lc).oot))
     if oot[0] ne -1 then $
        oplot,[((*self.active_lc).time)[oot]-subt]*modif,$
              [((*self.active_lc).flux)[oot]]+offset,psym=8,symsize=.6,color=(*self.colors).blue
     oot=0L
  endif
end

pro autokep::setactive,event
  lc = (*self.transits)[self.active_transit]->get()
  self.active_bound=[lc.midt-self.duration/2d0,lc.midt+self.duration/2d0]
  self->setactivelc,/noplot      
  (*self.active_lc)[where((*self.active_lc).time ge lc.midt - self.t_duration/2 and $
                          (*self.active_lc).time le lc.midt + self.t_duration/2)].transit = 1
  (*self.active_lc)[where((*self.active_lc).time lt lc.midt - self.t_duration/2 or $
                          (*self.active_lc).time gt lc.midt + self.t_duration/2)].oot = 1
  
  self->updateslider,/nodur
  lc = 0L
end


pro autokep::addtransit,midt

  ;help,/heap
  ;stop
  
  transit = {lc: *self.active_lc,$
             midt: midt,$
             epoch: 0d0,$
             oot_order: -1,$
             corr_lcf: dblarr(n_elements(*self.active_lc))}

  if self.num_transits eq 0 then begin
     ptr_free,self.transits
     self.transits = ptr_new(obj_new('transit'))
  endif else *self.transits = [*self.transits,obj_new('transit')]
  (*self.transits)[self.num_transits]->set,transit
  transit = 0L
  if self.num_transits eq 0 then *self.midts = midt else *self.midts = [*self.midts,midt]
  sort = sort(*self.midts)
  *self.midts = (*self.midts)[sort]  
  *self.transits = (*self.transits)[sort]
  sort = 0L
  self.num_transits+=1
  self.active_transit = self.num_transits-1
  if self.num_transits gt 1 then widget_control, self.sliders[3],sensitive=1
  if self.num_transits gt 0 then widget_control, self.bases[4],sensitive=1

  widget_control, self.sliders[4],sensitive=1 
  widget_control, self.sliders[5],sensitive=1


;  help,/heap
 ; stop  


  self->setepoch
  self->plot,'transits'
  self->updatefits,event,/init
  self->updateslider

end

pro autokep::map_epochs,event
  for i=0,self.num_transits-1,1 do begin
   ;  transit = (*self.transits)[i]->get()
   ;  ptr_free,transit.oot_trend
   ;  transit = 0L
;     obj_destroy,(*self.transits)[i]
 (*self.transits)[i]->destroy
  endfor
  ptr_free,self.transits
  self.transits = ptr_new()
  
  sig = robust_sigma((*self.all_lc).flux)
  self.num_transits=0
  
  self->setepoch
  epoch = *self.epoch
  for i=0,n_elements(epoch)-1,1 do begin
     self.active_bound=[epoch[i]-self.duration/2d0,epoch[i]+self.duration/2d0]
     self->setactivelc,/noplot      
     if ((*self.active_lc).flux)[0] ne -1 and n_elements((*self.active_lc)) gt 1 then begin
        transit = where((*self.active_lc).flux le 1d0-3d0*sig)
        if transit[0] ne -1 then begin
           (*self.active_lc)[transit].transit = 1  
           midpt = median((*self.active_lc)[transit].time)
           
           epoch[i] = midpt
           x =  akgapless((*self.active_lc).transit)
           if x ge 3 then begin
              self.t_duration = (x+4)*self.kep_dt
              self.duration = max([self.duration,self.t_duration+(20*self.kep_dt)])
              self.w_duration = self.duration
              
              self.active_bound[0] = midpt - self.duration/2d0
              self.active_bound[1] = midpt + self.duration/2d0
              self->setactivelc,/noplot
              
              (*self.active_lc)[where((*self.active_lc).time ge epoch[i] - self.t_duration/2 and $
                                      (*self.active_lc).time le epoch[i] + self.t_duration/2)].transit = 1
              (*self.active_lc)[where((*self.active_lc).time lt epoch[i] - self.t_duration/2 or $
                                      (*self.active_lc).time gt epoch[i] + self.t_duration/2)].oot = 1
              self->updateslider  

              self->addtransit,midpt     
              
              self->plot,''
           endif

           midpt=0L
           x = 0L
        endif
     endif
;stop
;     print,*self.epoch
  endfor
;  print,self.period
  
     self.period = median(abs(epoch[0:n_elements(epoch)-1] - epoch[1:n_elements(epoch)-1]))
     *self.epoch =  epoch
                                ; stop
     epoch = 0L
  
end

pro autokep::setactivelc,noplot=noplot
  lc_area = where((*self.all_lc).time ge self.active_bound[0] and $
                  (*self.all_lc).time le self.active_bound[1])
  
  if lc_area[0] eq -1 then *self.active_lc = {time: -1, flux: -1, err: -1} else begin
     lc = replicate({time: 0d0, flux: 0d0, err: 0d0, transit: 0d0, oot: 0d0, finalrange: 0d0},n_elements(lc_area))
     lc.time = ((*self.all_lc).time)[lc_area]
     lc.flux = ((*self.all_lc).flux)[lc_area]
     lc.err  = ((*self.all_lc).err)[lc_area]
     *self.active_lc = lc
     lc = 0L
  endelse
  lc_area = 0L
  
  if 1-keyword_set(noplot) then self->plot,'active_lc'
end

pro autokep::buttonevent,event
  widget_control, event.id, GET_UVALUE= uvalue
  widget_control, /Hourglass
  
  case uvalue.value of 
     'Fits File Button': begin
        path = dialog_pickfile(dialog_parent=self.bases[0],title='Select Fits File',/must_exist)
        if path ne '' then begin
           widget_control,self.loadfile_fld[1],set_value = path
        endif
     end
     'Clear Data Path':  widget_control,self.loadfile_fld[1],set_value = ''
     'Load all Fits': begin
        spawn,'ls *.fits',temp
        for i=0,n_elements(temp)-1,1 do begin
           path = temp[i]
           if strcmp((*self.fits_list)[0],'none') then *self.fits_list[0] = path[0] else $
              *self.fits_list = [*self.fits_list, path[0]]
            
        file = mrdfits(path[0],1)
        lc =  replicate({time: 0d0, flux: 0d0, err: 0d0},n_elements(where(finite(file.ap_corr_flux))))
        lc.time  = file[where(finite(file.ap_corr_flux))].barytime
        lc.flux = file[where(finite(file.ap_corr_flux))].ap_corr_flux
        lc.err  = file[where(finite(file.ap_corr_flux))].ap_corr_err 
        
        lc.err /= median(lc.flux)
        lc.flux /= median(lc.flux)
        
        if (*self.all_lc)[0].time eq -1 then new_lc = lc else begin
           new_lc = replicate({time: 0d0, flux: 0d0, err: 0d0},n_elements(lc)+n_elements(*self.all_lc))
           new_lc.time = [lc.time,(*self.all_lc).time]
           new_lc.flux = [lc.flux,(*self.all_lc).flux]
           new_lc.err = [lc.err,(*self.all_lc).err]
        endelse
        sort = sort(new_lc.time)
        new_lc.time = (new_lc.time)[sort]
        new_lc.flux = (new_lc.flux)[sort]
        new_lc.err = (new_lc.err)[sort]
        
        *self.all_lc = new_lc
        new_lc = 0L
        lc = 0L
        sort = 0L

        self->plot,''  
        endfor
        temp=0L
     end
     'Load Fits File': begin
        widget_control,self.loadfile_fld[1],get_value=path
        if strcmp((*self.fits_list)[0],'none') then *self.fits_list[0] = path[0] else $
           *self.fits_list = [*self.fits_list, path[0]]
        
        file = mrdfits(path[0],1)
        lc =  replicate({time: 0d0, flux: 0d0, err: 0d0},n_elements(where(finite(file.ap_corr_flux))))
        lc.time  = file[where(finite(file.ap_corr_flux))].barytime
        lc.flux = file[where(finite(file.ap_corr_flux))].ap_corr_flux
        lc.err  = file[where(finite(file.ap_corr_flux))].ap_corr_err 
        
        lc.err /= median(lc.flux)
        lc.flux /= median(lc.flux)
        
        if (*self.all_lc)[0].time eq -1 then new_lc = lc else begin
           new_lc = replicate({time: 0d0, flux: 0d0, err: 0d0},n_elements(lc)+n_elements(*self.all_lc))
           new_lc.time = [lc.time,(*self.all_lc).time]
           new_lc.flux = [lc.flux,(*self.all_lc).flux]
           new_lc.err = [lc.err,(*self.all_lc).err]
        endelse
        sort = sort(new_lc.time)
        new_lc.time = (new_lc.time)[sort]
        new_lc.flux = (new_lc.flux)[sort]
        new_lc.err = (new_lc.err)[sort]
        
        *self.all_lc = new_lc
        new_lc = 0L
        lc = 0L
        sort = 0L

        self->plot,''        
     end
     'apply_corr': if self.apply_corr then self.apply_corr = 0 else self.apply_corr = 1
     'secondary': if self.exp_second then self.exp_second = 0 else self.exp_second = 1
     'totap': if self.exp_to_tap then self.exp_to_tap = 0 else self.exp_to_tap = 1
     'Manual Detect': self->manualdetect,event
     'Manual Detect 2': begin
        self->manualdetect,event,/second
        if self.num_transits ge 2 then begin
           self.period = (*self.midts)[1]-(*self.midts)[0]
           self->plot,''
        endif
     end
     'Auto Detect': self->autodetect,event,/clear
     'Auto Detect 2': begin
        self->autodetect,event
        if self.num_transits ge 2 then begin
           self.period = (*self.midts)[1]-(*self.midts)[0]
           self->plot,''
        endif
     end
     'map': begin
        self->map_epochs,event
     end
     'adjmap': begin
        self->adj_map,event
     end
     'uq': self.yon = uvalue.set
     'init_fits': self->updatefits,event,/init,/all
     'exp_fname':begin
        if n_elements(*event.value) eq 0 then val = '' else val = *event.value
        self.exp_fname=val
        val = 0
     end
     'export':begin
        self->export
        if self.exp_to_tap then self->quit
     end
     else: print,'Unknown Button Event "'+uvalue.value+'"'
  endcase
end

pro autokep::export
  filename = self.exp_fname + '.ascii'
  self.message = 'Exporting LCs to '+filename
  
  format = '(d20.10,d20.10,d20.10)'
  openw,lun,filename,/get_lun,bufsize=0,width=250
  for i=0,self.num_transits-1,1 do begin
     if i ne 0 then printf,lun,-1,-1,-1,format=format
     lc = (*self.transits)[i]->get() 
     good = where(lc.lc.time ge (lc.midt-(self.w_duration/2d0)) and lc.lc.time le (lc.midt+(self.w_duration/2d0)))
     for j=0,n_elements(good)-1,1 do begin
        if self.apply_corr then $ 
           printf,lun,lc.lc[good[j]].time,lc.lc[good[j]].flux/lc.corr_lcf[good[j]],lc.lc[good[j]].err/lc.corr_lcf[good[j]],format=format else $
              printf,lun,lc.lc[good[j]].time,lc.lc[good[j]].flux,lc.lc[good[j]].err,format=format 
     endfor
  endfor
  close,/all

  

end

function autokep_colors
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
           green3:  15,$
           yellow: 16$
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
  tvlct,255,255,0,colors.yellow

  return,colors
  
end

pro autokep::message,event
  widget_control,self.message_window,/append,$
                 set_value='('+curr_date(format='hh:mm:ss')+') '+self.message
  self.message=''
end

function akmm,arr,modif=modif
  if 1-keyword_set(modif) then modif = 0d0
  val = mm(arr)
  valr = val[1]-val[0]
  val[1] += valr*modif
  val[0] -= valr*modif
  valr = 0L
  return,val
end

pro autokep::menumap,event
  Widget_Control, event.id, Get_UValue=info
  case info.type of
     'main': begin
        for i=1,n_elements(where(self.bases ne 0))-1,1 do widget_control,self.bases[i],MAP=0
        case info.value of
           '1) Load Fits': widget_control,self.bases[1],/map
           '2) Detect Transits': widget_control,self.bases[2],/map      
           '3) Adjust Transits': widget_control,self.bases[3],/map      
           '4) Export Transits': widget_control,self.bases[4],/map

        endcase
     end
  endcase
end



pro autokep::quitextra,event
  Widget_Control, event.id, Get_UValue=info
  widget_control,self.extra_windows[info.wid],/destroy
  info = 0L
end

pro autokep::updateslider,event,nodur=nodur
  if 1-keyword_set(nodur) then begin
     widget_control,self.sliders[1],set_value=[self.t_duration,self.t_duration/2d0,self.t_duration*2d0]
     widget_control,self.sliders[2],set_value=[self.duration,self.duration/2d0,self.duration*2d0]
  endif
  widget_control,self.sliders[3],set_value=[self.active_transit+1]
  widget_control,self.sliders[3],set_slider_max=max([self.num_transits,2])

  widget_control,self.sliders[6],set_value=[self.w_duration,self.t_duration,self.duration]

  if self.num_transits ge 1 then begin
     if ((*self.transits)[self.active_transit]->get()).oot_order ge 1 then begin
        widget_control,self.sliders[4],set_value = ((*self.transits)[self.active_transit]->get()).oot_order
  endif 
  endif
end

pro autokep::updatetransits,event
  for i=0,self.num_transits-1,1 do begin
     lc = (*self.transits)[i]->get() 
     self.active_bound=[lc.midt-self.duration/2d0,lc.midt+self.duration/2d0]
     self->setactivelc,/noplot      
     (*self.active_lc)[where((*self.active_lc).time ge lc.midt - self.t_duration/2 and $
                             (*self.active_lc).time le lc.midt + self.t_duration/2)].transit = 1
     (*self.active_lc)[where((*self.active_lc).time lt lc.midt - self.t_duration/2 or $
                             (*self.active_lc).time gt lc.midt + self.t_duration/2)].oot = 1
     
     lcn  = {lc: *self.active_lc,$
             midt: lc.midt,$
             epoch: lc.epoch,$
             oot_order: lc.oot_order,$
             corr_lcf: dblarr(n_elements(*self.active_lc))}

     if lcn.oot_order ne -1 then begin
        contf,lcn.lc.flux,cont,nord=lcn.oot_order,plot=0,frac=1,$
              mask=where(lcn.lc.oot),sbin=n_elements(where(lcn.lc.oot))
        lcn.corr_lcf = cont
     endif
     
     lc = 0L
     (*self.transits)[i]->set,lcn
     lcn=0L 
  endfor
  self->setactive,event
  self->plot,''
end

pro autokep::adjustslider,event
  widget_control, event.id, GET_UVALUE= uvalue
  case uvalue.value of
     'period': begin
        self.period = event.value
     end
     'duration2': begin 
        self.duration = event.value
        self.w_duration = self.duration
       ; self->setepoch
       ; self->plotepoch
     end
     'duration': begin 
        self.duration = event.value
        self.w_duration = self.duration
        self->updatetransits
        self->setepoch
     end
     't_duration': begin 
        self.t_duration = event.value
        self->updatetransits
        self->setepoch
     end
     'w_duration': begin 
        self.w_duration = event.value
        ;self->updatetransits
        ;self->setepoch
        self->plot,''
     end
     'active_t': begin
        self.active_transit = event.value-1
        
        self->setactive,event
        self->updateslider,event
     end
     'fitorder': begin
        if uvalue.doall then begin
           for i=0,self.num_transits-1,1 do begin 
              lc = (*self.transits)[i]->get()
              lc.oot_order = event.value
              (*self.transits)[i]->set,lc
              lc=0L
           endfor
           self->updatefits,/all,/holdplot
           self->setactive
        endif else begin
           lc = (*self.transits)[self.active_transit]->get()
           lc.oot_order = event.value
           (*self.transits)[self.active_transit]->set,lc
           lc=0L
           self->updatefits
        endelse
        self->plot,''
     end
  endcase

;  (*self.active_lc)[where((*self.active_lc).time ge (*self.epoch)[i] - self.t_duration/2 and $
;                          (*self.active_lc).time le epoch[i] + self.t_duration/2)].transit = 1
;  (*self.active_lc)[where((*self.active_lc).time lt epoch[i] - self.t_duration/2 or $
;                          (*self.active_lc).time gt epoch[i] + self.t_duration/2)].oot = 1
  

;  self->updateslider
 
  self->plot,''
end

pro autokep::adj_map,event
  self.extra_windows[0] = widget_base(title="Adjust Map Parameters",/column)
  XManager, 'autokep' $
            , self.extra_windows[0] $
            , /no_block   
  quit_button = widget_button(self.extra_windows[0],$
                              value = 'Quit',$
                              uvalue={object:self, method:'QuitExtra', wid: 0})
  
  self.sliders[0] = cw_fslider(self.extra_windows[0],title='Period',min=self.period/4d0,max=self.period*4d0,$
                               value = self.period,format='(g15.5)',/double,/edit,$
                               uname='period',uvalue={object:self, method:'AdjustSlider', $
                                                      value:'period'},$
                               drag=1)
  
  

  self.sliders[7] = cw_fslider(self.extra_windows[0],title='Duration',min=self.duration/2d0,max=self.duration*2d0,$
                               value = self.duration,format='(g15.5)',/double,/edit,$
                               uname='duration',uvalue={object:self, method:'AdjustSlider', $
                                                      value:'duration2'},$
                               drag=1)
  

  widget_control,self.extra_windows[0],/realize
end

pro autokep::updatefits,event,init=init,all=all,holdplot=holdplot
  if keyword_set(all) then begin
     act = self.active_transit
     for i=0,self.num_transits-1,1 do begin
        lc = (*self.transits)[i]->get() 
        
        if keyword_set(init) or lc.oot_order eq -1 then lc.oot_order = 1
                                ;*lc.oot_trend = polyfit(lc.lc[where(lc.lc.oot)].time,lc.lc[where(lc.lc.oot)].flux,lc.oot_order)
                                ;lc.corr_lcf = poly(lc.lc.time,*lc.oot_trend)
        contf,lc.lc.flux,cont,nord=lc.oot_order,plot=0,frac=1,$
              mask=where(lc.lc.oot),sbin=n_elements(where(lc.lc.oot))
        lc.corr_lcf = cont
        
        (*self.transits)[i]->set,lc
        
        ln=0L 
        if 1-keyword_set(holdplot) then begin
           self.active_transit = i
           self->setactive,event
           self->plot,''
        endif
     endfor
     self.active_transit = act
     act = 0L
     self->setactive,event
     self->plot,''       
  endif else begin
     lc = (*self.transits)[self.active_transit]->get() 
     
     if keyword_set(init) or lc.oot_order eq -1 then lc.oot_order = 1
     
     ;*lc.oot_trend = polyfit(lc.lc[where(lc.lc.oot)].time,lc.lc[where(lc.lc.oot)].flux,lc.oot_order)
     ;lc.corr_lcf = poly(lc.lc.time,*lc.oot_trend)
     contf,lc.lc.flux,cont,nord=lc.oot_order,plot=0,frac=1,$
           mask=where(lc.lc.oot),sbin=n_elements(where(lc.lc.oot))
     lc.corr_lcf = cont
     
     (*self.transits)[self.active_transit]->set,lc
     
     lc=0L 
     self->plot,'activeonly'
  endelse
  self->updateslider
end


pro autokep::plottrend
  if self.num_transits gt 0 then begin
     lc = (*self.transits)[self.active_transit]->get() 
     if lc.oot_order gt 0 then begin
        oplot,lc.lc.time-floor(min((*self.all_lc).time)),lc.corr_lcf,thick=2,color=(*self.colors).black
     endif
     lc = 0L
  endif
end

pro autokep::plot,type
  plot = lonarr(3)
  case type of
     'nolc': begin
        for i=0,2,1 do begin
           wset,(*self.plot_windows)[i].pix_window
           plot,[0],[0],/nodata,color=(*self.colors).gray,background=(*self.colors).gray
           wset,(*self.plot_windows)[i].w_id
           device,copy=[0,0,(*self.plot_windows)[i].x,(*self.plot_windows)[i].y,0,0,(*self.plot_windows)[i].pix_window]
           plot=[0,0,0]
        endfor
     end
     'full_lc': plot = [1,0,0]
     'active_lc': plot = [1,1,0]
     'activeonly': plot = [0,1,0]
     'transits': plot = [1,0,1]
     else: plot = [1,1,1]
  endcase
  
  if plot[0] then begin
     !p.multi = [0,1,1]
     wset,(*self.plot_windows)[0].pix_window
     (*self.plot_windows)[0].xrange = mm((*self.all_lc).time-floor(min((*self.all_lc).time)))
     (*self.plot_windows)[0].yrange = akmm((*self.all_lc).flux,modif=0.05d0)
     plot,[0],[0],$
          xr=(*self.plot_windows)[0].xrange,/xs,$
          yrange=(*self.plot_windows)[0].yrange,/ys,psym=8,symsize=.6,color=(*self.colors).black,$
          background=(*self.plot_windows)[0].bg,xmargin=[8,3],ymargin=[3,1],xthick=1.2,ythick=1.2,$
          xtitle='Days since '+sigfig(floor(min((*self.all_lc).time)),8),xticks=5,yticks=2,/nodata
     if self.period ne 0 then self->plotepoch
     oplot,(*self.all_lc).time-floor(min((*self.all_lc).time)),(*self.all_lc).flux,psym=8,symsize=.6,color=(*self.colors).black
     if self.num_transits ne 0 then self->plotall
     if total(self.active_bound) ne 0 then begin
        self->plotactive
     endif

     wset,(*self.plot_windows)[0].w_id
     device,copy=[0,0,(*self.plot_windows)[0].x,(*self.plot_windows)[0].y,0,0,(*self.plot_windows)[0].pix_window]
  endif 

  if plot[1] then begin
     if total(self.active_bound) ne 0 then begin
        !p.multi = [0,1,1]
        wset,(*self.plot_windows)[1].pix_window
        (*self.plot_windows)[1].xrange = mm((*self.active_lc).time-floor(min((*self.all_lc).time)))
        (*self.plot_windows)[1].yrange =(*self.plot_windows)[0].yrange
        if self.active_trf ge 0 then  (*self.plot_windows)[1].xrange[1] = (*self.plot_windows)[1].xrange[0]+self.active_trf
        if   (*self.plot_windows)[1].xrange[0] ne   (*self.plot_windows)[1].xrange[1] then begin
           plot,[0],[0],$
                xr=(*self.plot_windows)[1].xrange,/xs,$
                yrange=(*self.plot_windows)[1].yrange,/ys,psym=8,symsize=.6,color=(*self.colors).black,$
                background=(*self.plot_windows)[1].bg,xmargin=[8,3],ymargin=[3,1],xthick=1.2,ythick=1.2,$
                xtitle='Days since '+sigfig(floor(min((*self.all_lc).time)),8),xticks=3,yticks=4,/nodata
           self->plotactive
           self->plottrend
           self->plotwrite
        endif
        wset,(*self.plot_windows)[1].w_id
        device,copy=[0,0,(*self.plot_windows)[1].x,(*self.plot_windows)[1].y,0,0,(*self.plot_windows)[1].pix_window]
     endif
  endif
  
  if plot[2] then begin
     if self.num_transits gt 0 then begin
        !p.multi = [0,1,1]
        device,decomposed=0
        
        wset,(*self.plot_windows)[2].pix_window
        (*self.plot_windows)[2].xrange = [1d4,-1d4]
        for i=0,self.num_transits-1,1 do begin
           transit = (*self.transits)[i]->get()
           (*self.plot_windows)[2].xrange[0] = [min([(*self.plot_windows)[2].xrange[0],min(transit.lc.time)-(transit.midt)])]
           (*self.plot_windows)[2].xrange[1] = [max([(*self.plot_windows)[2].xrange[1],max(transit.lc.time)-(transit.midt)])]
           transit = 0L
        endfor

        tot = self.num_transits-1
        lci = (*self.transits)[0]->get()
        lcf = (*self.transits)[self.num_transits-1]->get()
        
        sepr = 8d0*robust_sigma(lci.lc[where(lci.lc.oot)].flux)
        
        yrange = [min(lci.lc.flux)-sepr,max(lcf.lc.flux)+((tot)*(3*sepr))+sepr]
        self->setepoch
        
        for i=0,self.num_transits-1,1 do begin
           lc = (*self.transits)[i]->get() 
           
           midt = (*self.midts)[i]
           epoch = *self.epoch
           val = where(abs(epoch-midt) eq min(abs(epoch-midt)))


           ;print,abs(epoch-midt)
           ;print,val+1
           midt1 = (*self.midts)[i]
           midt = lc.midt
           if midt1 ne midt then begin
              ;print,midt1,midt
              stop
           endif
           if i eq 0 then $
              plot,(lc.lc.time-midt)*24d0,lc.lc.flux+(3*sepr*i),psym=8,symsize=.6,yrange=yrange,color=(*self.colors).black,background=(*self.plot_windows)[2].bg,xrange=(*self.plot_windows)[2].xrange*24d0,/xs,/ys,xtitle='Hours from Mid Transit',ytitle='Rel. Flux + Constant',xmargin=[9,3],ymargin=[3,1] else $
                 oplot,(lc.lc.time-midt)*24d0,lc.lc.flux+(3*sepr*i),psym=8,symsize=.6,color=(*self.colors).black
           xyouts,(((*self.plot_windows)[2].xrange)[0]+0.1d0*(abs((*self.plot_windows)[2].xrange[1]-(*self.plot_windows)[2].xrange[0])))*24d0,1d0+(3*sepr*i)-sepr,string(val+1,format='(i0.0)'),color=(*self.colors).black,charsize=0.80
           lc = 0L
           midt = 0L
           epoch = 0L
           val = 0L
           if i eq self.active_transit then self->plotactive,offset=(3*sepr*i),/hours,/midt
        endfor
           self->plotwrite,/hours,/midt
           
        wset,(*self.plot_windows)[2].w_id
        device,copy=[0,0,(*self.plot_windows)[2].x,(*self.plot_windows)[2].y,0,0,(*self.plot_windows)[2].pix_window]
      ;  stop
     endif
     
  endif
  plot=0L
end


pro autokep::quit,even
  self.message = 'Destroying Widget'
  self->message

  
  self.buttons = 0L


  ptr_free,self.colors
  ptr_free,self.fits_list
  ptr_free,self.plot_windows
  ptr_free,self.menus
  ptr_free,self.all_lc
  ptr_free,self.active_lc
  for i=0,self.num_transits-1,1 do (*self.transits)[i]->destroy
  ptr_free,self.transits
  
  export=0L
  if self.exp_to_tap then begin
     export = 1L
     period = self.period
     ascii = self.exp_fname+'.ascii'
     smooth = [29.4244d0,5]
  endif
  
  widget_control,self.autokep_base,/destroy
  obj_destroy,self
  
  if export then tap,input_ascii=ascii,smooth=smooth,period=period
  ascii=0L
  smooth=0L
  period=0L

end

pro autokep::start
  centertlb,self.autokep_base
  widget_control,self.autokep_base,/realize

  for i=0,n_elements(*self.plot_windows)-1,1 do begin
     widget_control, (*self.plot_windows)[i].window, Get_Value=wid
     (*self.plot_windows)[i].w_id=wid
     
     window, xsize=(*self.plot_windows)[i].x $
             , ysize=(*self.plot_windows)[i].y $
             , /pixmap,/free
     (*self.plot_windows)[i].pix_window = !d.window
  endfor
  self->plot,'nolc'
  widget_control,self.bases[1],/map
  
end


pro autokep::widget_setup
  self.autokep_base = widget_base(title='autoKep Kepler LC Prep '+self.version,/column)
  
  XManager, 'autokep', self.autokep_base, /no_block 
  
  quit_button = widget_button(self.autokep_base,$
                              value = 'Quit',$
                              uvalue={object:self, method:'Quit'})
  
  message = '('+ curr_date(format='hh:mm:ss yyyymmdd') +') autoKep '+self.version
  self.message_window = widget_text(self.autokep_base, $
                                    value = message, $
                                    /scroll, $
                                    ysize=3)
  message = 0L
  
  main_base = widget_base(self.autokep_base,/column,$
                          frame=1, event_pro='autoKep_AllEvents')
  
  row_base = widget_base(main_base,/row)
  (*self.plot_windows)[0].x = 583d0
  (*self.plot_windows)[0].y = 100d0
  (*self.plot_windows)[0].window = widget_draw(row_base $
                                               , xsize=(*self.plot_windows)[0].x $
                                               , ysize=(*self.plot_windows)[0].y $
                                               , uvalue='Plot Window 1')
  row_base = widget_base(main_base,/row)
  (*self.plot_windows)[1].x = 290d0
  (*self.plot_windows)[1].y = 250d0
  (*self.plot_windows)[1].window = widget_draw(row_base $
                                               , xsize=(*self.plot_windows)[1].x $
                                               , ysize=(*self.plot_windows)[1].y $
                                               , uvalue='Plot Window 2')
  (*self.plot_windows)[2].x = 290d0
  (*self.plot_windows)[2].y = 250d0
  (*self.plot_windows)[2].window = widget_draw(row_base $
                                               , xsize=(*self.plot_windows)[2].x $
                                               , ysize=(*self.plot_windows)[2].y $
                                               , uvalue='Plot Window 3')


  menubar = widget_base(self.autokep_base,/row)
  
  row = widget_base(menubar,/row,/toolbar,/exclusive,/base_align_center)
  
  for i=0,n_elements(*self.menus)-1 do begin
     button = widget_button(row,$
                            value=' '+(*self.menus)[i]+' ',$
                            uvalue={object:self, method:'menumap', value:(*self.menus)[i], type:'main'},$
                            /no_release)
     if i eq 0 then widget_control, button, /SET_BUTTON
  endfor

  self.bases[0]  = widget_base(self.autokep_base,$
                               frame=1)
  
  self->setup,'load_fits'
  self->setup,'detect_transits' 
  self->setup,'adjust_transits'
  self->setup,'export_transits'

end

pro autokep::setup,set_this
  case set_this of 
     'load_fits': begin
        self.bases[1] = widget_base(self.bases[0] ,$
                                     /column, MAP=0,frame=0)
       ; workcol = widget_base(self.bases[1],/column,frame=1)
        row = widget_base(self.bases[1],/row,frame=0)

        button = widget_button(row,$
                               xsize=100,$
                               value='Fits File',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Fits File Button'})

        path = coyote_field2(row, TITLE=':', VALUE = self.loadpath,$
                             UVALUE='Fits File Field',XSIZE=65, TEXTID=textid)
        self.loadfile_fld = [path,textid]
        
        clear = widget_button(row, VALUE='Clear',$
                              uvalue={object:self, $
                                      method:'ButtonEvent', $
                                      value:'Clear Data Path'})
        
        row = widget_base(self.bases[1],/row,frame=0,/base_align_center)
        
        button = widget_button(row,xsize=100, value='Load Fits File',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Load Fits File'})
        label = widget_label(row,value="-or-")
        
        button = widget_button(row,xsize=100, value='Load all .fits',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Load all Fits'})

     end
     'detect_transits': begin
        self.bases[2] = widget_base(self.bases[0] ,$
                                    /column, MAP=0,frame=0)
        row = widget_base(self.bases[2],/row,frame=0)

        ;; FIRST TRANSIT
        base = widget_base(row,/column,frame=1)
        label = widget_label(base,value='1) First Transit')
        
        button = widget_button(base,xsize=100, value='Auto Detect',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Auto Detect'})
        
        button = widget_button(base,xsize=100, value='Manual Detect',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Manual Detect'})
        
        self.extra_bases[0] = widget_base(row,/column,frame=1,sensitive=0)
        label = widget_label(self.extra_bases[0],value='2) Second Transit')
        
        button = widget_button(self.extra_bases[0],xsize=100, value='Auto Detect',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Auto Detect 2'})
        
        button = widget_button(self.extra_bases[0],xsize=100, value='Manual Detect',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Manual Detect 2'})


        self.extra_bases[1] = widget_base(row,/column,frame=1,sensitive=0)
        label = widget_label(self.extra_bases[1],value='3) Map Remaining Transits')
        

        button = widget_button(self.extra_bases[1],xsize=100, value='Adj. Yellow Maps',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'adjmap'})
        
        button = widget_button(self.extra_bases[1],xsize=100, value='Map Transits',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'map'})
        
     end
     'adjust_transits': begin
        self.bases[3] = widget_base(self.bases[0] ,$
                                    /column, MAP=0,frame=0)
        row = widget_base(self.bases[3],/row,frame=0)
      
        ;; Overall Params
        base = widget_base(row,/column,frame=1)
        self.sliders[1] =  cw_fslider(base,title='Duration [Days]',min=self.t_duration/2,max=self.t_duration*2,$
                               value = self.duration,format='(g15.5)',/double,/edit,$
                               uname='t_duration',uvalue={object:self, method:'AdjustSlider', $
                                                        value:'t_duration'},$
                                      drag=1)
        
        self.sliders[2] =  cw_fslider(base,title='OOT [Days]',min=self.duration/2,max=self.duration*2,$
                               value = self.duration,format='(g15.5)',/double,/edit,$
                               uname='duration',uvalue={object:self, method:'AdjustSlider', $
                                                        value:'duration'},$
                                      drag=1)
        base = widget_base(row,/column,frame=1)
        label = widget_label(base, value='Transit')
        self.sliders[3] = widget_slider(base,/vertical,minimum=1,maximum=max([self.num_transits,2]),$
                                        value=self.active_transit+1,uvalue={object:self, method:'AdjustSlider', $
                                                                            value:'active_t'},/drag,$
                                        ysize=120,sensitive=0)
        base = widget_base(row,/column,frame=1)
        label = widget_label(base, value='Fit OOT Trends')
        
;        button = widget_button(base,xsize=100, value='Initial Fits',$
 ;                              uvalue={object:self, $
  ;                                     method:'ButtonEvent', $
   ;                                    value:'init_fits'})
        label = widget_label(base, value='active order:')
                
        self.sliders[4] = widget_slider(base,minimum=1,maximum=5,$
                                        value=1,uvalue={object:self, method:'AdjustSlider', $
                                                        value:'fitorder',doall:0},/drag,$
                                        xsize=50,sensitive=0)
        label = widget_label(base, value='all order:')
        
        self.sliders[5] = widget_slider(base,minimum=1,maximum=5,$
                                        value=1,uvalue={object:self, method:'AdjustSlider', $
                                                        value:'fitorder',doall:1},/drag,$
                                        xsize=50,sensitive=0)
     end
     'export_transits': begin
        self.bases[4] = widget_base(self.bases[0] ,$
                                    /column, MAP=0,frame=0,sensitive=0)
        col =  widget_base(self.bases[4],/col,frame=0)
        row = widget_base(col,/row,frame=0)
        
        base = widget_base(row,/column,frame=1)
        label = widget_label(base, value='1) Write Region:')
        

        self.sliders[6] =  cw_fslider(base,title='Write Region [Days]',min=self.t_duration,max=self.duration,$
                                      value = self.w_duration,format='(g15.5)',/double,/edit,$
                                      uname='w_duration',uvalue={object:self, method:'AdjustSlider', $
                                                                 value:'w_duration'},$
                                      drag=1)
        
        base = widget_base(row,/column,frame=1)
        label = widget_label(base, value='2) Export Options:')
        radio = widget_base(base,column=1,/nonexclusive)
        button = widget_button(radio,value='Apply OOT Corrections',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'apply_corr'})
        widget_control,button,set_button=self.apply_corr
        
     ;   radio = widget_base(base,column=1,/nonexclusive)
     ;   button = widget_button(radio,value='Write Secondary Eclipse',$ 
     ;                          uvalue={object:self, method:'buttonevent',$
     ;                                  value: 'secondary'})
     ;   widget_control,button,set_button=self.exp_second
          
        radio = widget_base(base,column=1,/nonexclusive)
        button = widget_button(radio,value='Export Transits to TAP',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'totap'})
        widget_control,button,set_button=self.exp_to_tap
        base = widget_base(row,/column,frame=1)
        label = widget_label(base, value='3) Filename Root (.ascii appended):')
        path = coyote_field2(base, TITLE='', VALUE = self.exp_fname,$
                             UVALUE={object:self, method:'buttonevent',value:'exp_fname'},XSIZE=33, TEXTID=textid,event_pro='autokep_event')
        
        row = widget_base(col,/row,frame=0,/align_center,/base_align_center)  
        button = widget_button(row,xsize=200,ysize=50,value='Export Transits',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'export'})
        
     end
     else: begin
        print,set_this+' is an unknown setup!'
        stop
     end
  endcase
end


function autokep::init
  !p.font = 0
  plotsym,0,1.2,/fill
  device,decomposed=0
 
  ptr_free,self.colors
  self.colors = ptr_new(autokep_colors())
  
  self.version = 'v1.01'
  self.apply_corr = 1

  ptr_free,self.plot_windows
  self.plot_windows = ptr_new(replicate({x: 0d, y: 0d, w_id: 0L, pw_id: 0L, window:0L, pix_window: 0L, xrange: [0d0,0d0], yrange: [0d0,0d0], bg:(*self.colors).white},3))
  
  ptr_free,self.menus
  self.menus = ptr_new(['1) Load Fits','2) Detect Transits','3) Adjust Transits','4) Export Transits'])
  ptr_free,self.midts
  self.midts = ptr_new(-1)
  self.exp_fname = 'autoKepExport'
  
  ptr_free,self.fits_list
  self.fits_list=ptr_new('none')
  ptr_free,self.all_lc,self.active_lc
  self.all_lc = ptr_new({time: -1, flux: -1, err: -1})
  self.active_lc = ptr_new({time: -1, flux: -1, err: -1})
  self.active_trf = -1
  
  self.kep_dt = 29.4d0/(60d0*24d0)
  self.duration = 24d0*self.kep_dt
  self.t_duration = 12d0*self.kep_dt
  self.w_duration = self.duration
  
  self->widget_setup
 ; widget_control,self.bases[1],MAP=1        
    
  Widget_Control,self.autokep_base,set_UVALUE=self 
  self->start
  
  return,1
end

pro autokep__define
  struct = {autokep,$
            $ ;; GENERAL STUFF
            version: '',$ 
            colors: ptr_new(), $
            kep_dt: 0d0,$
            $ ;; fits files
            fits_list: ptr_new(),$
            loadpath: '',$
            loadfile_fld: dblarr(2),$
            $ ;; MAIN WIDGET stuff
            autokep_base: 0L, $
            message_window: 0L, $
            message: '',$
            plot_windows: ptr_new(),$
            buttons: lonarr(5),$
            sliders: lonarr(8),$
            extra_bases: lonarr(5),$
            extra_windows: lonarr(5),$
            menus: ptr_new(),$
            bases: lonarr(5),$
            $ ;; DATA
            all_lc: ptr_new(),$
            active_lc: ptr_new(),$
            active_bound: dblarr(2),$
            active_trf: 0d0,$
            transits: ptr_new(),$
            active_transit: 0,$
            num_transits: 0,$
            yon: 0L,$
            $ ;; LC search
            midts: ptr_new(-1),$
            period: 0d0,$
            duration: 0d0,$
            t_duration: 0d0,$
            w_duration: 0d0,$
            exp_to_tap: 0,$
            exp_second: 0,$
            exp_fname: '',$
            apply_corr: 0,$
            epoch: ptr_new(),$
            exists: 0 $
           }
end
  
pro autokep
  autokep = obj_new('autokep')
end
