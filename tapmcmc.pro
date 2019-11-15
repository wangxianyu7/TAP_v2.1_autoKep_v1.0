
pro tapmcmc::destroy,event
  ptr_free,self.need_adjust
  ptr_free,self.adjust
  ptr_free,self.dec
  ptr_free,self.inc
  ptr_free,self.phi
  ptr_free,self.currlikes
  ptr_free,self.newlikes
  ptr_free,self.jumps
  ptr_free,self.sets
  
  obj_destroy,self.progbarobj1
  obj_destroy,self.progbarobj2
  
  obj_destroy,self.rand_obj
  obj_destroy,self.rand_obj2

  ptr_free,self.plot_windows

  for i=0,n_elements(*self.widget_bases)-1,1 do if (*self.widget_bases)[i] then widget_control,(*self.widget_bases)[i],/destroy
  ptr_free,self.widget_bases
  
  ptr_free,self.colors
  widget_control,self.mcmc_base,/destroy
  obj_destroy,self
 end

pro tapmcmc::storelink,event,pick=pick 
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.params[pick].jumptot++
     if (*self.jumps)[i] then begin
        transit.params.curr_link = transit.params.new_link
        transit.params[pick].jumpct++
        (*self.currlikes)[i] = (*self.newlikes)[i]
                                ;transit.curr_redl = transit.new_redl
     endif
     (*self.transits)[i]->set,transit
     transit = 0L
  endfor
  (*self.jumps)*=0d0
  
  if (self.iter mod 10) eq 0 then begin
     for i=0,n_elements(*self.transits)-1,1 do begin
        transit = (*self.transits)[i]->get()    
        for j=0,n_elements(transit.params)-1,1 do $
           *transit.params[j].mcmc_chain = [*transit.params[j].mcmc_chain,transit.params[j].curr_link]
                             
        (*self.transits)[i]->set,transit
        transit = 0L
     endfor
     self->updatemod
     
     if (self.iter mod 1000) eq 0 then begin
        self->curr_analyze
        self.base_widget->prepcol,runval=self.runval
     endif
     
     if (self.iter mod 100) eq 0 then self->plot
  endif
end

pro TAPmcmc::lcplot,event
  
  !p.multi = [0,1,1]
                                ;device,decomposed=0
  wset,(*self.plot_windows)[0].pix_window
  (*self.plot_windows)[0].xrange = [1d4,-1d4]
  
  if self.which_plot ge self.num_transits then begin
     for i=0,n_elements(*self.transits)-1,1 do  (*self.plot_windows)[0].xrange = [min([(*self.plot_windows)[0].xrange[0],min(((*self.transits)[i]->get()).transit_tday-(((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value-(((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))]),$
                                                                               max([(*self.plot_windows)[0].xrange[1],$
                                                                                    max(((*self.transits)[i]->get()).transit_tday-$
                                                                                        (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))])]
  

  tot = self.num_transits-1
  lci = (*self.transits)[0]->get()
  lcf = (*self.transits)[self.num_transits-1]->get()
  diff = ((*self.transits)[0]->get()).diff
  if diff eq 0 then diff = 5d0*max([lci.residuals,lcf.residuals])

  if self.phased then begin
     yrange = [min(lci.transit_flux)-(diff/2d0)-diff,max(lcf.transit_flux)+diff]
     
   ;  diff = 0L
     
     for i=0,self.num_transits-1,1 do begin
        lc = (*self.transits)[i]->get() 
                                ;   diff = 5d0*max(lc.residuals)
        midt = lc.params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
        if i eq 0 then $
           plot,(lc.transit_tday-midt)*24d0,lc.transit_flux,psym=8,symsize=.6,yrange=yrange,color=lc.color,background=(*self.colors).white,xrange=(*self.plot_windows)[0].xrange*24d0,/xs,/ys,xtitle='Hours from Mid Transit',ytitle='Relative Flux + Constant' else $
              oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux,psym=8,symsize=.6,color=lc.color
        
   ;     oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)),psym=8,symsize=.6,color=lc.color
    ;    hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=2
    ;    oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=2
    ;    oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=2
        lc = 0L   
        midt = 0L
     endfor
     lci = 0L
     lcf = 0L
     tot = 0L
     yrange = 0L
     
  endif else begin 
     yrange = [min(lci.transit_flux)-tot*(diff/2d0)-diff,max(lcf.transit_flux)+((tot+.5)*diff)]
     
     for i=0,self.num_transits-1,1 do begin
        
        lc = (*self.transits)[i]->get() 
                                ;   diff = 5d0*max(lc.residuals)
        midt = lc.params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
        if i eq 0 then $
           plot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,yrange=yrange,color=lc.color,background=(*self.colors).white,xrange=(*self.plot_windows)[0].xrange*24d0,/xs,/ys,xtitle='Hours from Mid Transit',ytitle='Relative Flux + Constant' else $
              oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,color=lc.color
        
     oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(diff*i),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)*(i+1)),psym=8,symsize=.6,color=lc.color
        hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=2
        lc = 0L   
        midt = 0L
     endfor
     lci = 0L
     lcf = 0L
     tot = 0L
     yrange = 0L
     ;; plot,((*self.transits)[0]->get()).transit_tday-((*self.transits)[0]->get()).params[where(strcmp(((*self.transits)[0]->get()).params.param,'Mid Transit') eq 1)].value,((*self.transits)[0]->get()).transit_flux,psym=8,symsize=.4,background=(*self.colors).white,color= ((*self.transits)[0]->get()).color,xrange=(*self.plot_windows)[0].xrange
     
  endelse
endif  else begin
   (*self.plot_windows)[0].xrange = [min([(*self.plot_windows)[0].xrange[0],min(((*self.transits)[self.which_plot]->get()).transit_tday-(((*self.transits)[self.which_plot]->get()).params[where(strcmp(((*self.transits)[self.which_plot]->get()).params.param,'Mid Transit') eq 1)].value-(((*self.transits)[self.which_plot]->get()).params[0].value*((*self.transits)[self.which_plot]->get()).epoch)))]),$
                                     max([(*self.plot_windows)[0].xrange[1],$
                                          max(((*self.transits)[self.which_plot]->get()).transit_tday-$
                                              (((*self.transits)[self.which_plot]->get()).params[where(strcmp(((*self.transits)[self.which_plot]->get()).params.param,'Mid Transit') eq 1)].value - (((*self.transits)[self.which_plot]->get()).params[0].value*((*self.transits)[self.which_plot]->get()).epoch)))])]
   i=0
   lc = (*self.transits)[self.which_plot]->get()
   diff = ((*self.transits)[0]->get()).diff
   if diff eq 0 then diff = 5d0*max(lc.residuals)
   midt = lc.params[where(strcmp(lc.params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
   yrange = [min(lc.transit_flux)-(diff/2d0)-diff,max(lc.transit_flux)+(.5*diff)]

   plot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,yrange=yrange,color=(*self.colors).black,background=(*self.colors).white,xrange=(*self.plot_windows)[0].xrange*24d0,/xs,/ys,xtitle='Hours from Mid Transit',ytitle='Relative Flux + Constant',title=string(self.which_plot+1,format='(i2.2)')+': '+lc.fname
   oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),color=lc.color,psym=8,symsize=.6
     oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(diff*i),color=lc.modcol,thick=2
   oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)*(i+1)),psym=8,symsize=.6,color=lc.color
   hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=2
   oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=2
   oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=2
   lc = 0L
endelse


  wset,(*self.plot_windows)[0].w_id
  device,copy=[0,0,(*self.plot_windows)[0].x,(*self.plot_windows)[0].y,0,0,(*self.plot_windows)[0].pix_window]
  
end

pro TAPmcmc::oldlcplot,event
  !p.multi = [0,1,1]
  device,decomposed=0
  wset,(*self.plot_windows)[0].pix_window
  

  (*self.plot_windows)[0].xrange = [1d4,-1d4]
  for i=0,n_elements(*self.transits)-1,1 do  (*self.plot_windows)[0].xrange = [min([(*self.plot_windows)[0].xrange[0],min(((*self.transits)[i]->get()).transit_tday-(((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value-(((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))]),$
                                                                               max([(*self.plot_windows)[0].xrange[1],$
                                                                                    max(((*self.transits)[i]->get()).transit_tday-$
                                                                                        (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))])]
  
  tot = self.num_transits-1
  lci = (*self.transits)[0]->get()
  lcf = (*self.transits)[self.num_transits-1]->get()
  yrange = [.99*min(lci.transit_flux)-tot*.005,1.00*max(lcf.transit_flux)+((tot)*.015)]
  
  for i=0,self.num_transits-1,1 do begin
     lc = (*self.transits)[i]->get() 
 
     midt = lc.params[where(strcmp(lc.params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
    
     
     if i eq 0 then begin
        plot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(.015*i),psym=8,symsize=.6,yrange=yrange,color=(*self.colors).black,background=(*self.colors).white,xrange=(*self.plot_windows)[0].xrange*24d0,/xs,/ys,/nodata,xtitle='Hours since Mid Transit',ytitle='Relative Flux + Constant'
        oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(.015*i),psym=8,symsize=.6,color=lc.color
     endif else $
           oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(.015*i),psym=8,symsize=.6,color=lc.color
     
     oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(.015*i),color=lc.modcol,thick=2
     oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+(.005*(i+1)),psym=8,symsize=.6,color=lc.color
     hline,yrange[0]+(.005*(i+1)),color=lc.modcol,thick=2
     oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(.015*i)+lc.rednoise,color=lc.redcol,thick=2
     oplot,(lc.transit_tday-midt)*24d0,yrange[0]+(.005*(i+1))+lc.rednoise,color=lc.redcol,thick=2
     oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(.015*i),psym=8,symsize=.6,color=lc.color

     lc = 0L   
     midt = 0L
  endfor

  lci = 0L
  lcf = 0L
  tot = 0L
  yrange = 0L
  wset,(*self.plot_windows)[0].w_id
  device,copy=[0,0,(*self.plot_windows)[0].x,(*self.plot_windows)[0].y,0,0,(*self.plot_windows)[0].pix_window]
 ; stop
  
end

pro tapmcmc::Message,event
  self.message = '('+curr_date(format='hh:mm:ss')+') '+self.message
  widget_control, (*self.widget_bases)[0],/append,$
                  set_value=self.message

  self.message=''
end

pro tapmcmc::setuprun,throw=throw
  if 1-keyword_set(throw) then throw=0d
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     tags = n_elements(transit.params.param)
     for j=0,tags-1,1 do begin
        transit.params[j].jumpct = 0d0
        transit.params[j].jumptot = 0d0

        if transit.params[5].value + transit.params[6].value ge 1 then transit.params[6].value = 9d-1 - transit.params[5].value
        if transit.params[5].value + transit.params[6].value le 0 then transit.params[6].value = -0.99d0*transit.params[5].value
      ;  stop
        ptr_free,transit.params[j].mcmc_chain
        transit.params[j].mcmc_chain = ptr_new(transit.params[j].value)
        if transit.params[j].limited[0] then $
           if (*transit.params[j].mcmc_chain)[0] lt transit.params[j].limits[0] then $
              (*transit.params[j].mcmc_chain)[0] = transit.params[j].limits[0]      
        if transit.params[j].limited[1] then $
           if (*transit.params[j].mcmc_chain)[0] gt transit.params[j].limits[1] then $
              (*transit.params[j].mcmc_chain)[0] = transit.params[j].limits[1]      
         
        transit.params[j].curr_link = (*transit.params[j].mcmc_chain)[0]
     endfor
     transit.new_redl = 0L
     transit.curr_redl = 0L
     (*self.transits)[i]->set,transit
     transit = 0L
     tags = 0L
  endfor
  
  self.jumpcount*=0
  self.jumptot*=0
  self.iter = 0d0
  self->likelihood,/current
  self->updatemod
  self->plot
end

pro tapmcmc::updatemod
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.params.value = transit.params.curr_link
     rednoise = dblarr(n_elements(transit.transit_tday))
     dum = TAP_MPFit_function(transit.params.value,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,redfit=rednoise,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int)
     transit.residuals=transit.transit_flux - dfit
     dum = 0L
     transit.rednoise=rednoise
     *transit.model_f = dfit
     transit.residuals = transit.transit_flux - dfit
     *transit.model_ffine = ffit
     redfit=0L
     ffit = 0L
     dfit = 0L

     (*self.transits)[i]->set,transit
     rednoise = 0L
     transit = 0L
  endfor
end

pro tapmcmc::guessbetas
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.params.beta = transit.params.value/1d3
     transit.params.beta = [6d-5,0.07,0.05,0.0004,0.0004,0.04,0.05,0.0001,0.0001,0.0001,1.8d-05,0.0001,7d-05]
     transit.params[4].beta = 1d-2
     if (where(transit.params.beta eq 0))[0] ne -1 then transit.params[where(transit.params.beta eq 0)].beta = 1d-4
     (*self.transits)[i]->set,transit
     transit = 0L
  endfor
end

pro tapmcmc::plot
  xtitle= 'Parameter Value'
  histplot = -1
  bin = 25
  maxthick = 6

  if (self.iter mod 1000) eq 0 and self.iter gt 250 then begin
     if (self.iter mod 10000) eq 0 then self.plot+=1
     if (self.iter mod 2500) eq 0 then begin
        if self.which_plot gt self.num_transits then self.which_plot = -1
        self.which_plot +=1
     endif
     if self.plot ne 13 then begin
        k = 0
        while k lt 14 do begin
           if self.plot ge 12 then self.plot = 0
           if (self.jumpcount)[self.plot] ge 20 and self.plot ne 4 then k = 14 else begin
              k++
              self.plot +=1
           endelse
        endwhile
        k = 0L
     endif
     histplot = self.plot
  endif
  

  if histplot ge 0 then begin
     xr = [1d50,-1d50]
     yr = [0,0]
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        r = [.1*n_elements(*transit.params[self.plot].mcmc_chain),n_elements(*transit.params[self.plot].mcmc_chain)-1]
        xr = [min([xr[0],min((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]])]),max([xr[1],max((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]])])]
        transit = 0L
        r = 0L
     endfor
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        r = [.1*n_elements(*transit.params[self.plot].mcmc_chain),n_elements(*transit.params[self.plot].mcmc_chain)-1]
        if stdev((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]]) ne 0 then begin
           plothist,(*transit.params[self.plot].mcmc_chain)[r[0]:r[1]],tx,ty,bin=(max(xr)-min(xr))/bin,/noplot
                                ;  stop
           yr = [0,max([yr[1],max(ty)])]
        endif
      ;  print,yr
        tx = 0L
        ty = 0L
        transit = 0L
        r = 0L
     endfor
    ; print,yr
     yr[1]*=1.1
     !p.multi=[0,1,1]
     window = self.base_widget->window()
     wset,window.pix_window
     first = 0
     for i=0,n_elements(*self.transits)-1,1 do begin
        transit = (*self.transits)[i]->get()
        ; added by Xian-Yu
        ;print, ((*self.transits)[i]->get()).color
        if transit.params[self.plot].prior[0] eq 1 then extra = ' (Gaussian Penalty)' else extra=''
        
        r = [.1*n_elements(*transit.params[self.plot].mcmc_chain),$
             n_elements(*transit.params[self.plot].mcmc_chain)-1]
        if first eq 0 then begin
           if stdev((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]]) ne 0 then begin
              plothist,(*transit.params[self.plot].mcmc_chain)[r[0]:r[1]], $
                       bin=(max(xr)-min(xr))/bin,color=((*self.transits)[i]->get()).color, $
                       xtitle=xtitle, $background=(*self.colors).white,axiscolor=(*self.colors).black,$
                       title=((*self.transits)[i]->get()).params[self.plot].param+extra,xrange=xr,/xs, $
                       thick=maxthick,xticks=3,xthick=2,ythick=2,yrange=yr,/ys
              first = 1
           endif
        endif else begin
           if stdev((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]]) ne 0 then begin
              plothist,(*transit.params[self.plot].mcmc_chain)[r[0]:r[1]], $
                       bin=(max(xr)-min(xr))/bin,color=((*self.transits)[i]->get()).color, $
                       xtitle=xtitle, $background=(*self.colors).white,axiscolor=(*self.colors).black,$
                       title=((*self.transits)[i]->get()).params[self.plot].param,/overplot,xrange=xr, $
                       /xs, xticks=3,thick=maxthick*(n_elements(*self.transits)-i)/n_elements(*self.transits),yrange=yr,/ys
           endif
        endelse
        if transit.params[self.plot].prior[0] eq 1 then begin
           vline,transit.params[self.plot].prior[1],thick=3,color=((*self.transits)[i]->get()).color
           vline,transit.params[self.plot].prior[1]+transit.params[self.plot].prior[2],thick=3,color=((*self.transits)[i]->get()).color,linestyle=2
           vline,transit.params[self.plot].prior[1]-transit.params[self.plot].prior[2],thick=3,color=((*self.transits)[i]->get()).color,linestyle=2
        endif
        
        
        transit = 0L
        r = 0L
     endfor
     xr = 0L
     yr = 0L
     wset,window.w_id
     device,copy=[0,0,window.x,window.y,0,0,window.pix_window]
     window = 0L
  endif
  bin = 0L
  histplot = 0L
  
  
  self->lcplot
                            
  
  
;stop
end

pro tapmcmc::likelihood,current=current,new=new
  likelihood = 0L
  if keyword_set(new) then *self.newlikes *= 0d0 else  *self.currlikes *= 0d0
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     if keyword_set(current) then $
        params = transit.params.curr_link else $
           params=transit.params.new_link
     dum = TAP_MPFit_function(params,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,redfit=rednoise,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int)
     
     
     likelihood = waveletlike((transit.transit_flux-dfit),params[11],params[12],/zeropad)
                                ;transit = 0L
                                ;transit = (*self.transits)[i]->get()
     
     if (where(transit.params.prior[0] eq 1))[0] ne -1 then begin
        these = where(transit.params.prior[0] eq 1)
        for j=0,n_elements(these)-1,1 do likelihood -= ((params[these[j]]-transit.params[these[j]].prior[1])^2d0/(transit.params[these[j]].prior[2])^2d0)
        these = 0L
     endif
     
     if keyword_set(new) then begin
        ;transit.new_redl = likelihood 
        (*self.newlikes)[i] = likelihood
        if i eq 0 then self.new_redl = likelihood else self.new_redl += likelihood
     endif else begin
        ;transit.curr_redl = likelihood
        (*self.currlikes)[i] = likelihood
        ;if finite(transit.curr_redl) eq 0 then transit.curr_redl = -1d5
        if finite((*self.currlikes)[i]) eq 0 then  (*self.currlikes)[i] = -1d5
        if i eq 0 then self.curr_redl = likelihood else self.curr_redl += likelihood
     endelse
    ; (*self.transits)[i]->set,transit
     transit = 0L 
     likelihood = 0L
  endfor  
end

pro tapmcmc::testheap,title
  print,title
  title=0L
  help,/heap
 stop
  
end

pro tapmcmc::burnbetas,event
  if self.testheap then self->testheap,'burnbetas 0'
     
  lim1 = 1.15d0
  lim2 = .85d0

  coarse_adjust = 1
  if keyword_set(fine) then coarse_adjust = 0
  fine_adjust=2
  if coarse_adjust ne 0 or fine_adjust ne 0 then self->guessbetas
  self->setuprun
  initial_lock = 0d
  stabilized   = 0d 
  self.message = ' Stabilizing Betas'
  if 1-keyword_set(fine) then   self->message
 ; if self.testheap then self->testheap,'burnbetas 1'

  if coarse_adjust then begin
     tot_iter = 2d3
     while not initial_lock do begin
        tot_iter--
      ;  if self.testheap then self->testheap,'pre newlink'
        self->newlink
      ;  if self.testheap then self->testheap,'post newlink'
        if (self.iter mod 500d) eq 0 then begin
           for i=0,n_elements(*self.transits)-1,1 do begin
              ptr_free,self.adjust,self.dec,self.inc
              transit = (*self.transits)[i]->get()
              self.adjust = ptr_new(transit.params.beta*0d0)
              self.dec = ptr_new(where((transit.params.jumpct/transit.params.jumptot) gt transit.params.accept*1.25d0))
              self.inc = ptr_new(where((transit.params.jumpct/transit.params.jumptot) le transit.params.accept*0.75d0))
              fixit = ptr_new(where((transit.params.fixed) eq 1))     
              transit = 0L
              if (*self.dec)[0] ge 0 then (*self.adjust)[(*self.dec)] = 1.d0
              if (*self.inc)[0] ge 0 then (*self.adjust)[(*self.inc)] = -1.d0
              if (*fixit)[0] ge 0 then (*self.adjust)[(*fixit)] = 0.d0
              if (where(*self.adjust ne 0))[0] ge 0 then begin
                 ;;self.beta_adjusttype= 'crude' ; crude adjustment
                 
                 ;;if self.testheap then self->testheap,'burnbetas: Pre crude adhyst'
                 self->AdjustBeta,'crude',i
                 ;;if self.testheap then self->testheap,'burnbetas: Post crude adhyst'
                                  
                 
                                ;  print,'  adjusting:', (*self.adjust)
                                ;  print,'  NEW BETAS:', (*self.betas_store)
              endif else begin
                 initial_lock=1
              endelse
              ptr_free,fixit
           endfor
        endif
        if tot_iter le 0 then initial_lock = 1
      ;;; MORE
     endwhile

     tot_iter = 0L
     self.message = '  Coarse Stabilization Done'
     self->message
   ;  if self.testheap then self->testheap,'burnbetas 2: post coarse'


 ;;; MORE?
  endif else initial_lock = 1
 
  
  while fine_adjust gt 0 do begin
     if self.testheap then self->testheap,'burnbetas 3: fineloop'
     
     while stabilized lt n_elements(*self.transits) do begin 
       ; if self.testheap then self->testheap,'burnbetas 3.1: pre newlink'
     
        self->newlink       
        stabilized = 0
        for i=0,n_elements(*self.transits)-1,1 do begin
           cs = 0
           transit = (*self.transits)[i]->get()
           check = where(transit.params.fixed eq 0)
           if self.debug then begin
              if i eq 0 then begin
                 if (self.iter mod 250d0) eq 0 or self.iter lt 3 then begin
                    print,'-1? : ',(where((transit.params[check].accept*(lim1+(.1*(max([self.iter-6d3,0])/6d3))) lt (transit.params[check].jumpct/transit.params[check].jumptot)) or $
                                          (transit.params[check].accept*(lim2-(.1*(max([self.iter-6d3,0])/6d3))) gt (transit.params[check].jumpct/transit.params[check].jumptot))))
                    print,'1d4< ',self.iter       
                    print,transit.params[check].param
                    print,transit.params[check].jumptot
                    print,'-1?  ',(where(transit.params[check].jumptot le (200d0-min([150d0,(max([self.iter-6d3,0])/250d0)]))))
                 
                 endif
              endif
           endif
           cs =  (where((transit.params[check].accept*(lim1+(.1*(max([self.iter-6d4,0])/6d3))) lt (transit.params[check].jumpct/transit.params[check].jumptot)) or (transit.params[check].accept*(lim2-(.1*(max([self.iter-6d3,0])/6d3)))  gt (transit.params[check].jumpct/transit.params[check].jumptot))))
           cs2 = where(transit.params[check].jumptot le (200d0-min([150d0,(max([self.iter-6d3,0])/250)])))
           if self.debug and (self.iter mod 500) eq 0 then print,cs
           if cs[0] eq -1 and  self.iter gt 1d4 and cs2[0] eq -1 then begin
              if self.debug then if i eq 0 then begin
                 print,'GOOD!: fine_adj(1)',fine_adjust
                 print,i,stabilized
                                ;   stop
              endif
              stabilized += 1
              if fine_adjust eq 1 and i eq 0 then begin
                 self.message = '  Fine Stabilization Done'
                 self->message
              endif
           endif 
           if cs[0] ne -1 then begin
              ptr_free,self.adjust,self.need_adjust
              self.adjust = ptr_new(transit.params.beta*0d0)
              self.need_adjust = ptr_new(where(((transit.params.jumpct/transit.params.jumptot) - transit.params.accept)^2d0 gt $
                                               (transit.params.sval*((transit.params.accept*(1-transit.params.accept))/transit.params.jumptot)) and $
                                               transit.params.jumptot gt 2d2))
              if (*self.need_adjust)[0] ne -1 then (*self.adjust)[(*self.need_adjust)] = 1.d0
              
          
              if (where((*self.adjust) eq 1))[0] ne -1 then begin
                 if self.debug then print,stabilized
                 self->AdjustBeta,'fine',i
                                ;self->message
              endif             ;else if  stabilized += 1
           endif   
           transit = 0L
           check= 0L
           
        endfor
     endwhile
     stabilized = 0L
     fine_adjust -= 1
  endwhile
  if self.testheap then self->testheap,'burnbetas 5: post fine'
  
  self.message=' Betas have stabilized'
  self->message
                                ;stop
end

pro tapmcmc::adjustbeta,type,active
  lim1 = 1.2d0
  lim2 = 0.8d0
  transit = (*self.transits)[active]->get()
  if self.debug then begin
     print,active
     print,'  rates:    ', (transit.params.jumpct/transit.params.jumptot)
     print,'   tots:    ', (transit.params.jumptot)
     print,'  betas:    ', transit.params.beta
     print,'  change:   ', *self.adjust
     print,'  iters:    ', self.iter, .44*(lim1+(.1*(max([self.iter-6d3,0])/6d3))),.44*(lim2-(.1*(max([self.iter-6d3,0])/1d4))),200d0-min([150d0,(max([self.iter-6d3,0])/250)])
     print,''
  endif
  case type of 
     'crude': begin
        if (where((*self.adjust) eq  1))[0] ge 0 then transit.params[where((*self.adjust) eq  1)].beta *= 1.5d0
        if (where((*self.adjust) eq -1))[0] ge 0 then transit.params[where((*self.adjust) eq -1)].beta /= 1.5d0
     end
     'fine':begin
        ptr_free,self.phi
        self.phi= ptr_new(dblarr(n_elements(transit.params.beta))*0d)
        area = where(transit.params.jumpct/transit.params.jumptot lt .1)
        if area[0] ne -1 then (*self.phi)[area] = .5d0
        area = where((0.5d*(transit.params.accept)) lt (transit.params.jumpct/transit.params.jumptot))
        if area[0] ne -1 then (*self.phi)[area] = 1d0
        area = where(((0.2d*transit.params.accept lt (transit.params.jumpct/transit.params.jumptot)) and $
                      ((0.5d*transit.params.accept) ge (transit.params.jumpct/transit.params.jumptot))))
        if area[0] ne -1 then  (*self.phi)[area] = 1.5d0
        area = where(((0.1d*transit.params.accept lt (transit.params.jumpct/transit.params.jumptot)) and $
                      ((0.2d*transit.params.accept) ge (transit.params.jumpct/transit.params.jumptot))))
        if area[0] ne -1 then   (*self.phi)[area] = 2d0
        
        transit.params[where((*self.adjust) eq 1)].beta = (transit.params.beta*((((transit.params.jumpct+1)/transit.params.jumptot)/transit.params.accept)^(*self.phi)))[where((*self.adjust) eq 1)]
        

        inc = 0
        dec = 0
        
        for i=0,n_elements(transit.params.beta)-1,1 do begin
           if (*self.adjust)[i] then begin
              inc = -1                                                                                                ; assume decrease
              if ((transit.params[i].jumpct/transit.params[i].jumptot)/(transit.params[i].accept)) le 1 then inc = 1 ; increase
              
              case transit.params[i].b_last of
                 0: transit.params[i].b_last = inc
                 -1: if inc eq 1 then transit.params[i].sval++
                 1: if inc eq -1 then transit.params[i].sval++
              endcase
              transit.params[i].b_last = inc     
           endif
        endfor
        transit.params[where(*self.adjust eq 1)].jumpct = 0d
        transit.params[where(*self.adjust eq 1)].jumptot = 0d
     end
  endcase
  (*self.transits)[active]->set,transit
  transit=0L
end


pro tapmcmc::newlink
  ;if self.testheap then self->testheap,'newlink0'
          
  fail = 1
  pick = ptr_new()
  while fail do begin 
     fail = 0
     all_lock = 0
     ptr_free,pick
     pick = ptr_new(sort(self.rand_obj2->getrandomnumbers(13,/uniform)))
     linklock = ptr_new([-1])
     linkval = ptr_new([-1])
     for i=0,n_elements(*self.transits)-1,1 do begin
        if fail eq 0 then begin
           transit = (*self.transits)[i]->get()
           transit.params.new_link = transit.params.curr_link
           if transit.params[(*pick)[0]].fixed ne 1 then begin 
              if (where(*linklock eq transit.params[(*pick)[0]].set))[0] ne -1 then begin
                 transit.params[(*pick)[0]].new_link = (*linkval)[where(*linklock eq transit.params[(*pick)[0]].set)]
              endif else begin
                 transit.params[(*pick)[0]].new_link += (self.rand_obj->getrandomnumbers(1,/normal,/double)*(transit.params[(*pick)[0]].beta))
                 *linklock = [*linklock,transit.params[(*pick)[0]].set]
                 *linkval = [*linkval,transit.params[(*pick)[0]].new_link]
              endelse
           endif else all_lock+=1

           if transit.params[5].new_link + transit.params[6].new_link ge 1d0 or transit.params[5].new_link + transit.params[6].new_link le 0d0 then fail = 1 else $
              if transit.params[(*pick)[0]].limited[0] eq 1 then if transit.params[(*pick)[0]].new_link le transit.params[(*pick)[0]].limits[0] then fail=1 else $
                 if transit.params[(*pick)[0]].limited[1] eq 1 then if transit.params[(*pick)[0]].new_link gt transit.params[(*pick)[0]].limits[1] then fail=1
           ;; transit.params.curr_link = transit.params.new_link
           if fail eq 0 then (*self.transits)[i]->set,transit
           transit = 0L
        endif
     endfor
     ptr_free,linklock,linkval

     if all_lock eq n_elements(*self.transits) then fail = 1

  endwhile
  
  self.jumptot[(*pick)[0]]++
  self->likelihood,/new
  tjump = self.rand_obj2->getrandomnumbers(1,/uniform)
  
                                ; stop
  sets = (*self.sets)[*,(*pick)[0]]
  uniqsets = uniq(sets[sort(sets)])
  for i=0,n_elements(uniqsets)-1,1 do begin
     set = (sets[sort(sets)])[uniqsets[i]]
     prob = exp((total((*self.newlikes)[where(sets eq set)])-total((*self.currlikes)[where(sets eq set)]))/2d0)
     if tjump le min([prob,1d0]) then begin
        (*self.jumps)[where(sets eq set)] = 1
     endif
  endfor
  sets = 0L
  uniqsets = 0L
  if max(*self.jumps) then self.jumpcount[(*pick)[0]]++

  self.iter++
  prob = 0L
  tjump = 0L

  self->StoreLink,pick=(*pick)[0]
  ptr_free,pick
end

pro tapmcmc::savechain

  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.mcmc_complete = 1
     (*self.transits)[i]->set,transit
     transit = 0L
  endfor
  
  tap_state = (*self.transits)
  
;  stop
  save,TAP_state,filename='MCMC_chains/'+self.save_header+'_'+string(self.active_chain,format='(i2.2)')+'of'+string(self.num_chains,format='(i2.2)')+'.idlsav'
  TAP_state = 0L

  restore,'TAP_setup.idlsav'
  adjust= TAP_state[0]->get()
  if adjust.mcmc_complete eq 0 then *adjust.mcmc_files='MCMC_chains/'+self.save_header+'_'+string(self.active_chain,format='(i2.2)')+'of'+string(self.num_chains,format='(i2.2)')+'.idlsav' else $
     *adjust.mcmc_files=[*adjust.mcmc_files,'MCMC_chains/'+self.save_header+'_'+string(self.active_chain,format='(i2.2)')+'of'+string(self.num_chains,format='(i2.2)')+'.idlsav']
  adjust.mcmc_complete++
  TAP_State[0]->set,adjust
  adjust=0L
  save,TAP_state,filename='TAP_setup.idlsav'
  
end

function tapmcmc::info
  return,{savefile: self.save_header+'/TAP_setup.idlsav',$
          version: self.version}
end

pro tapmcmc::start,event
  centertlb,self.mcmc_base
  widget_control,self.mcmc_base,/realize
  
  widget_control, (*self.plot_windows)[0].window, Get_Value=wid
  (*self.plot_windows)[0].w_id=wid
   window,xsize=(*self.plot_windows)[0].x, $
         ysize=(*self.plot_windows)[0].y, $
         /pixmap,/free
  (*self.plot_windows)[0].pix_window = !d.window
  
  self.message = self.save_header
  self->message
  
  self->ExecuteMCMC
end

pro tapmcmc::runchain,event
  if self.initial_run eq 1 then begin
     self->BurnBetas
     if self.testheap then self->testheap,'postbb'
     
     self->SetupRun,throw=0d
     if self.testheap then self->testheap,'postbb post setup'
     
  endif else self->setuprun
  
  if self.testheap then self->testheap,'begin init run'
  for i=1d0,self.chain_length[0],1d0 do begin
     if (self.iter mod 1000d0) eq 0 then begin
        widget_control,self.progbar1,set_value=((self.iter/self.chain_length[0]))
        widget_control,self.progbar2,set_value=((self.iter+(self.chain_length[0]*(self.active_chain-1)))/(self.num_chains*self.chain_length[0]))
     endif
     self->newlink
  endfor
  if self.testheap then self->testheap,'end init run'

  self->updatemod
  widget_control,self.progbar1,set_value=0d0
  widget_control,self.progbar2,set_value=((self.iter+(self.chain_length[0]*(self.active_chain-1)))/(self.num_chains*self.chain_length[0]))
end

pro tapmcmc::ExecuteMCMC,event
  self->lcplot
  ;; run the initial loops:
  
  self.initial_run = 1
  set = n_elements(*((*self.transits)[0]->get()).mcmc_files)
  if strcmp((*((*self.transits)[0]->get()).mcmc_files)[0],'-1') eq 0 then set++
  
  for i=set,self.num_chains,1 do begin
     self.active_chain = i
     widget_control,self.progbar1,set_value=((self.iter/self.chain_length[0]))
     widget_control,self.progbar2,set_value=((self.iter+(self.chain_length[0]*(self.active_chain-1)))/(self.num_chains*self.chain_length[0]))
     self.message = 'Running Chain ' + string(self.active_chain,format='(i2.2)') + ' of ' + string(self.num_chains,format='(i2.2)')
     self->message
     self->RunChain
     self->savechain
  endfor
  
  self.initial_run = 0
  
  self.message = 'MCMC Run Complete'
  self->message
  
  cd,'../'
end


pro tapmcmc::curr_analyze
  self.runval++
  if self.runval gt self.num_transits then self.runval = 1
  transit = (*self.transits)[self.runval-1]->get()
  for j=0,n_elements(transit.params)-1,1 do begin
     num = n_elements(*transit.params[j].mcmc_chain)
     sorted = ((*transit.params[j].mcmc_chain)[(.1*num):(num-1)])[sort((*transit.params[j].mcmc_chain)[(.1*num):(num-1)])]
     range = n_elements(sorted)/100d0
     transit.params[j].runval = [sorted[50d0*range],$
                                 sorted[84.135d0*range]-sorted[50d0*range],$
                                 sorted[50d0*range]-sorted[15.865d0*range]]
     if transit.params[j].fixed then transit.params[j].runval[1:2] = -1d0
     
     num = 0L
     range = 0L
     sorted = 0L
  endfor
  (*self.transits)[self.runval-1]->set,transit
  transit = 0L
end

function tapmcmc::INIT,$
   base_parent=base_parent,$
   input_transits=input_transits,$
   restart=restart, version=version
  
  
  if keyword_set(restart) then begin
     self.save_header = restart
     cd,self.save_header
     self.save_header = (strsplit(self.save_header,'/',/extract))[n_elements(strsplit(self.save_header,'/'))-1]
;     print,self.save_header
  endif else begin
     
     ;Modified by Xian-Yu, 20191115
     ;self.save_header = 'd:\'+'TAPmcmc_'+curr_date(format='yyyymmdd_hhmm')
     self.save_header = 'd:\'+'TAPmcmc_'+curr_date(format='yyyymmdd_hhmm')

     spawn,'mkdir '+self.save_header
     cd,self.save_header
     print,'The results folder:'+self.save_header
     spawn,'mkdir MCMC_chains'
  endelse
  
  ptr_free,self.colors
  self.colors = ptr_new(tap_colors())
  obj_destroy,self.base_widget
  self.base_widget = base_parent
  !EXCEPT = 0
  self.debug = 0
  self.testheap = 0
  self.version = version
  
  self.mcmc_base =  widget_base(frame=1,/column,title='TAP MCMC '+self.version) 
  
  work_base = widget_base(self.mcmc_base,/column,frame=1)
  Xmanager, 'TAPmcmc',$
            self.mcmc_base,$
            /no_block  
  
  ptr_free,self.plot_windows
  self.plot_windows = ptr_new(replicate({x: 0d, y: 0d, w_id: 0L, pw_id: 0L, window:0L, pix_window: 0L, xrange: [0,0], yrange: [0,0]},1))
  (*self.plot_windows)[0].x = 300d0
  (*self.plot_windows)[0].y = 300d0
  (*self.plot_windows)[0].window = $
     widget_draw(work_base,$
                 xsize=(*self.plot_windows)[0].x,$
                 ysize=(*self.plot_windows)[0].y,$
                 uvalue = 'mcmc plot 1',/align_center,frame=1)
  
  obj_destroy,self.progbarobj1,self.progbarobj2
  self.progbar1 = cw_progress(work_base,obj_ref=self.progbarobj1,/red,ysize=10d,xsize=(*self.plot_windows)[0].x)                 
  self.progbar2 = cw_progress(work_base,obj_ref=self.progbarobj2,/red,ysize=10d,xsize=(*self.plot_windows)[0].x)
  
  ptr_free,self.widget_bases
  self.widget_bases = ptr_new(lonarr(2))
  
  message = '('+ curr_date(format='hh:mm:ss yyyymmdd') +') TAP MCMC '+self.version
  (*self.widget_bases)[0] = widget_text(work_base, $
                                        font = textfont, $
                                        value = message, $
                                        /scroll, $
                                        ysize=5)
  message = 0L


  ptr_free,self.transits
  self.transits = input_transits
  self.num_transits = n_elements(*self.transits)

  self.which_plot = self.num_transits  
  
  transit = (*self.transits)[0]->get()
  transit.mcmc_version = self.version
  (*self.transits)[0]->set,transit
  transit = 0L

  TAP_state = (*self.transits)
  save,TAP_state,filename='TAP_setup.idlsav'
  tap_state = 0L
  
  self.num_chains = ((*self.transits)[0]->get()).mcmc_params[0]
  self.chain_length =((*self.transits)[0]->get()).mcmc_params[1:2]
  
  if 1-keyword_set(restart) then begin
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        ptr_free,transit.params.mcmc_chain
        for j=0,n_elements(transit.params)-1,1 do begin
           ptr_free,transit.params[j].mcmc_chain
           transit.params[j].mcmc_chain=ptr_new()           
        endfor        
        transit.mcmc_complete = 0
        (*self.transits)[i]->set,transit
        transit = 0L
     endfor
  endif else self.restart = restart
 
  ptr_free,self.sets
  self.sets = ptr_new(lonarr(self.num_transits,n_elements(((*self.transits)[0]->get()).params)))
  for i=0,self.num_transits-1,1 do (*self.sets)[i,*] =  ((*self.transits)[i]->get()).params.set
  ptr_free,self.jumps,self.newlikes,self.currlikes
  self.jumps = ptr_new(lonarr(self.num_transits))
  self.newlikes = ptr_new(dblarr(self.num_transits))
  self.currlikes = ptr_new(dblarr(self.num_transits))
  

  self.rand_obj = obj_new('tap_RandomNumberGenerator')
  self.rand_obj2 = obj_new('tap_RandomNumberGenerator')
  
  return,1
end



pro tapmcmc__define
  struct = { tapmcmc, $
             colors: ptr_new(),$
             version: '',$
             debug: 0L,$
             testheap: 0L,$
             restart: '',$
             $ ;; widget stuff
             base_widget: obj_new(),$
             mcmc_base: 0L, $
             widget_bases: ptr_new(),$
             message_window: 0L, $
             message: '',$
             plot_windows: ptr_new(),$
             phased: 0L,$
             phasect: 0L,$
             progbar1: 0L,$
             progbar2: 0L,$
             progbarobj1: obj_new(),$
             progbarobj2: obj_new(),$
             $ ;; RaNDOM Numers!
             rand_obj: obj_new(),$
             rand_obj2: obj_new(),$
             $ 
             sets: ptr_new(),$
             jumps: ptr_new(),$
             newlikes: ptr_new(),$
             currlikes: ptr_new(),$
             $
             initial_run: 0,$
             active_chain: 0,$
             num_chains: 0,$
             chain_length: dblarr(2),$
             $
             jump:0L ,$
             jumpcount: dblarr(13),$
             jumptot:   dblarr(13),$
             $
             adjust: ptr_new(),$
             need_adjust:  ptr_new(),$
             dec: ptr_new(),$
             inc: ptr_new(),$
             which_plot: 0,$
             $
             save_header: '',$
             $
             runval: 0L,$
             plot: 0d0,$
             iter: 0d0,$
             phi: ptr_new(),$
             num_transits: 0,$
             curr_redl: 0d0,$
             new_redl: 0d0,$
             transits: ptr_new()  $
           }
  
  




end

function tapmcmc,base_parent=base_parent,input_transits=input_transits,only_version=only_version,_REF_EXTRA=_extra
  version = 'v2.1'
  if keyword_set(only_version) then return,{version: version} else $
     return,obj_new('tapmcmc',base_parent=base_parent,input_transits=input_transits,version=version,_EXTRA=_extra)
end
