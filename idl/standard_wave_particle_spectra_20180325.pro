pro standard_wave_particle_spectra_20180325
  ;; example crib sheet to load fff data for Y. Shen
  thm_init
  ;;setups to make the plot better
  device,decomposed=0
  loadct,39
  !p.background=255
  !p.color=0
  !p.charsize = 0.8
  !p.thick = 3.0


  tplot_options,'xmargin',[25.,12.]
  tplot_options,'ymargin',[3,3]
  tplot_options,'xticklen',0.05
  tplot_options,'yticklen',0.03
  tplot_options,'xthick',2
  tplot_options,'ythick',2
  tplot_options,'ygap',0.25   ;eliminate the vertical spacing between panels

  axisthick = 4.0
  charthick = 4.0
  linethick = 2.0
  charsize = 3.0
  symsize = 1.0
  landscape=0 ; set landscape = 1 if you want to make any postscripts generated in landscape

  encapsulated=0 ; set encapsulated = 1 if you want any postscripts generated to be eps

  savepath='/Applications/harris/spedas/idl/' ;;change to your local directory
  probes=['a']
  ;probes=['a']
  ;;various ranges of conjunction intervals
  tranges=[['2018-03-25/20:40','2018-03-25/21:00']]
  ; , ['2014-08-21/09:00','2014-08-21/11:10']]
  ntot=n_elements(tranges[0,*])
  for i=0.,ntot-1. do begin
    tr=reform(tranges[*,i])
    timespan,time_double(tr(0)), time_double(tr(1))-time_double(tr(0)),/sec
    sname=time_string(tr[0],/date_only,format=2)
    for j=0., n_elements(probes)-1. do begin
      sc=probes[j]
      thm_load_state, probe = sc, trange = tr,/get_support_data
      thm_cotrans, probe=sc, 'th'+sc+'_state_pos','th'+sc+'_state_pos_sm',in_coord='gei',out_coord='gsm'
      MLT_L_lat,'th'+sc+'_state_pos_sm',MLT0=MLT0,L0=L0,lat0=lat0
      get_data, 'th'+sc+'_state_pos_sm', data = sc_pos
      
;      store_data, 'L', data={x: sc_pos.x, y: L0}
      store_data, 'R', data={x: sc_pos.x, y: sqrt(sc_pos.y[*,0]^2+sc_pos.y[*,1]^2+sc_pos.y[*,2]^2)/6371.0}
      
      store_data, 'MLT', data={x: sc_pos.x, y: MLT0}
      store_data, 'Lat', data={x: sc_pos.x, y: lat0*180.0/!pi}

      ;      options,'L',ytitle='L',ysubtitle='[Re]'
      ;      options,'MLT',ytitle='MLT',ysubtitle='[hr]'
      thm_load_fgm, probe = sc, trange = tr, level=2, coord='dsl', data='fgs'
      thm_load_fgm, probe = sc, trange = tr, level=2, coord='dsl', data='fgl'
      thm_load_fgm, probe = sc, trange = tr, level=2, coord='gsm', data='fgl'
      thm_load_fgm, probe = sc, trange = tr, level=2, coord='gsm', data='fgs'
      thm_load_fft, probe=sc, trange=tr, level=2,datatype='fff*'
      tvectot, 'th'+sc+'_fgs_dsl',tot='Bmag'
      get_data, 'Bmag', data = bmag
    
      fce_eq = 28.*abs(bmag.y)
      store_data,'Pressure_B',data={x:bmag.x,y:(bmag.y)^2/(2*1.26)*1e-3}
      

            store_data,strjoin('th'+sc+'_fce_eq'),data={x:bmag.x,y:fce_eq},dlim={colors:'w',thick:'2',linestyle:'0'}
      ;      store_data,'th'+sc+'_2fce_eq',data={x:bmag.x,y:fce_eq*2.0},dlim={colors:0,thick:'2',linestyle:'2'}
      ;      store_data,'th'+sc+'_3fce_eq',data={x:bmag.x,y:fce_eq*3.0},dlim={colors:0,thick:'1',linestyle:'0'}
      ;      store_data,'th'+sc+'_15fce_eq',data={x:bmag.x,y:fce_eq*1.5},dlim={colors:0,thick:'1',linestyle:'2'}
      ;      store_data,'th'+sc+'_25ce_eq',data={x:bmag.x,y:fce_eq*2.5},dlim={colors:0,thick:'1',linestyle:'2'}
            store_data,'th'+sc+'_fce_eq_01',data={x:bmag.x,y:fce_eq*0.1},dlim={colors:'k',thick:'1.5',linestyle:'0'}
      store_data,'th'+sc+'_fce_eq_half',data={x:bmag.x,y:fce_eq*0.5},dlim={colors:'r',thick:'1.5',linestyle:'0'}
      ;            store_data,'th'+sc+'_fcp_eq',data={x:bmag.x,y:fce_eq/1837},dlim={colors:1,thick:'2',linestyle:'0'}
      ;      store_data,'th'+sc+'_fche_eq',data={x:bmag.x,y:fce_eq/1837/4.0},dlim={colors:1,thick:'1',linestyle:'4'}
      ;      store_data,'th'+sc+'_fco_eq',data={x:bmag.x,y:fce_eq/1837/16.0},dlim={colors:1,thick:'1',linestyle:'1'}
      store_data,strjoin('th'+sc+'_fLHR'),data={x:bmag.x,y:fce_eq/1837.*43.},dlim={colors:'r',thick:'1.5',linestyle:'0'}
      


      tplot_names,'th'+sc+'_ff?_*_scm?',names=bvar
      tplot_names,'th'+sc+'_ff?_*_e?c??',names=evar
      get_data,bvar[0],data=thx_b1
      get_data,bvar[1],data=thx_b2
      store_data,'th'+sc+'_fff_32_bsum',data={x:thx_b1.x,y:thx_b1.y+thx_b2.y,v:thx_b1.v}
      get_data,evar[0],data=thx_e1,dlim=dlim
      get_data,evar[1],data=thx_e2
      store_data,'th'+sc+'_fff_32_esum',data={x:thx_e1.x,y:thx_e1.y+thx_e2.y,v:thx_e1.v}
      
      options,['th'+sc+'_fff_32_bsum','th'+sc+'_fff_32_esum'],spec=1
      options,strjoin('th'+sc+'_fff_32_bsum'),zrange=[1e-9,1e-2]
      options,strjoin('th'+sc+'_fff_32_esum'),/zlog,zrange=[1e-12,1e-5]
      
      get_data,strjoin('th'+sc+'_fff_32_bsum'),data=b_sum
      get_data,strjoin('th'+sc+'_fff_32_esum'),data=e_sum
      
      time_clip,strjoin('th'+sc+'_fff_32_bsum'),time_double(tr(0)),time_double(tr(1)),newname='th'+sc+'_fff_32_bsum_clip'
      tinterpol_mxn, 'th'+sc+'_fce_eq','th'+sc+'_fff_32_esum', newname ='th'+sc+'_fce_eq_int'
      get_data,'th'+sc+'_fce_eq_int',data=fce
      
      
      whistler_pow1=dblarr(n_elements(b_sum.x))
      for i=0,n_elements(b_sum.x)-1 do begin
        if i eq 0 then begin
          ind=where(b_sum.y[i,*] gt 1e-4 and b_sum.v[*] lt 3000.0 and b_sum.v[*] gt 100.0)
        endif else begin
          ind=[]
          ind=where(b_sum.y[i,*] gt 1e-4 and b_sum.v[*] lt 3000.0 and b_sum.v[*] gt 100.0)
        endelse
        if ind[0] le 0 then begin
          whistler_pow1[i]=0.0
        endif else begin
          whistler_pow1[i]=sqrt(1.5*total(e_sum.y[i,ind]*(e_sum.v[ind+1]-e_sum.v[ind]),2,/nan))*1000.
        endelse
      endfor
      
      whistler_pow2=dblarr(n_elements(b_sum.x))
      for i=0,n_elements(b_sum.x)-1 do begin
        if i eq 0 then begin
          ind=where(b_sum.y[i,*] gt 1e-4 and b_sum.v[*] lt 3000.0 and b_sum.v[*] gt 50.0)
        endif else begin
          ind=[]
          ind=where(b_sum.y[i,*] gt 1e-4 and b_sum.v[*] lt 3000.0 and b_sum.v[*] gt 50.0)
        endelse
        if ind[0] le 0 then begin
          whistler_pow2[i]=0.0
        endif else begin
          whistler_pow2[i]=sqrt(1.5*total(e_sum.y[i,ind]*(e_sum.v[ind+1]-e_sum.v[ind]),2,/nan))*1000.
        endelse
      endfor
      
      ech_pow=dblarr(n_elements(e_sum.x))
      for i=0,n_elements(e_sum.x)-1 do begin
        if i eq 0 then begin
          ind=where(e_sum.v[*] ge fce.y[i] and e_sum.v[*] lt fce.y[i]*1.6)
        endif else begin
          ind=[]
          ind=where(e_sum.v[*] ge fce.y[i] and e_sum.v[*] lt fce.y[i]*1.6)
        endelse
        if ind[0] le 0 or total(e_sum.y[i,ind]) lt total(e_sum.y[i,(ind[0]-n_elements(ind)+1):ind[0]]) then begin
          ech_pow[i]=0.0
        endif else begin
          ech_pow[i]=sqrt(total(e_sum.y[i,ind]*(e_sum.v[ind+1]-e_sum.v[ind]),2,/nan))*1000.
        endelse
      endfor

      
      
      ffactor=thx_e1.v[7:31]-thx_e1.v[6:30] ; including frequeny bins between 100 and 4000 Hz
      wfactor=thx_e1.v[4:31]-thx_e1.v[3:30]
      num=n_elements(thx_e1.x)
      
      n_factors=dblarr(num,n_elements(ffactor))
      w_factors=dblarr(num,n_elements(wfactor))
      for i=0,num-1 do begin
        n_factors[i,*]=ffactor
        w_factors[i,*]=wfactor
      endfor
      store_data,'thx_integrate_e',data={x:thx_e1.x[0:num-1],y:sqrt(1.5*total((thx_e1.y[0:num-1,7:31]+thx_e2.y[0:num-1,7:31])*n_factors,2,/nan))*1000.}
      store_data,'thx_integrate_e_2',data={x:thx_e1.x[0:num-1],y:sqrt(1.5*total((thx_e1.y[0:num-1,4:31]+thx_e2.y[0:num-1,4:31])*w_factors,2,/nan))*1000.}
;      time_clip,strjoin('thx_integrate_e'),thx_e1.x[0],thx_e1.x[num-1],newname=strjoin('thx_integrate_e_clip')
      get_data,'thx_integrate_e',data=e_tot
      get_data,'thx_integrate_e_2',data=e_tot_2
      root_e1=dblarr(n_elements(b_sum.x))
      root_e2=dblarr(n_elements(b_sum.x))
      for i=0,n_elements(b_sum.x)-1 do begin
        if (e_tot.y[i]-whistler_pow1[i]) gt 0.0 then begin
          root_e1[i]=e_tot.y[i]-whistler_pow1[i]
        endif else begin
          root_e1[i]=0.0
        endelse
      endfor
      for i=0,n_elements(b_sum.x)-1 do begin
        if (e_tot_2.y[i]-whistler_pow2[i]) gt 0.0 then begin
          root_e2[i]=e_tot_2.y[i]-whistler_pow2[i]
        endif else begin
          root_e2[i]=0.0
        endelse
      endfor
;      store_data,'integrate_tot',data={x:b_sum.x,y:[[root_e2],[root_e1]]}
      store_data,'integrate_tot',data={x:b_sum.x,y:root_e1}
      options, 'integrate_tot',yrange=[0.0,6],ytitle='TDS Ew',ylog=0,ysubtitle='[mV/m]',colors='k',labels='';,labels=['100-4000 Hz!cexcluding whistler!cand ECH power'],LABFLAG=1

      store_data,'integrate_tot_2',data={x:b_sum.x,y:root_e2}
      options, 'integrate_tot_2',yrange=[0.0,6],ytitle='TDS Ew',ylog=0,ysubtitle='[mV/m]',colors='r',labels='';,labels=['100-4000 Hz!cexcluding whistler!cand ECH power'],LABFLAG=1
      
     
      ;      tplot,'thx_integrate_e_clip'
      ;      get_data,'thx_integrate_e_clip',data=int_esum
      ;      openw,lun3,'THA_20140827_0600_0625_Int_E2_Vm_g4Hz.txt',/get_lun
      ;      for w=0,n_elements(int_esum.y)-1 do begin
      ;        printf,lun3, int_esum.y(w), Format='(E12.4)'
      ;      endfor
      ;      Free_lun, lun3
      ;
      ;      stop
      
;      ;here we get electric field and drift data
;
;      mysc='th'+sc
;      sclet=sc
;      thm_load_state,probe=sclet,/get_supp
;      thm_load_fit,probe=sclet,trange=tr,suff='_dsl' ; this gives both FGM and EFI (spin-fits)
;     
;
;
;      ;; first get Ex/y offsets
;
;      trange4offset=['2018-03-25/20:40:00','2018-03-25/21:00:00']
;      Exyzoffsets = spd_tplot_average(mysc+'_efs_dsl',trange4offset)
;      Exoffset = Exyzoffsets[0]
;      Eyoffset = Exyzoffsets[1]
;
;      get_data,'th'+sclet+'_efs_0_dsl',data=thx_efs_dsl
;
;      thx_efs_dsl.y(*,0)=thx_efs_dsl.y(*,0)-Exoffset
;      thx_efs_dsl.y(*,1)=thx_efs_dsl.y(*,1)-Eyoffset
;      
;
;      store_data,'th'+sclet+'_efs_new_dsl',data={x:thx_efs_dsl.x,y:thx_efs_dsl.y}
;      options,'th'+sclet+'_efs_new_dsl','colors',['g','b','r']
;      options,'th'+sclet+'_efs_new_dsl','databar',0
;
;      ; convert to GSM
;      ; 
;      thm_cotrans,mysc+'_efs_new',in_suffix='_dsl',out_suffix='_gsm'
;      ;
;      ; interpolate across gaps (when Bz is small and E*B=0 cannot be implemented resulting in gaps)
;      ;
;      tdeflag,mysc+'_efs_new_gsm','linear',/overwrite
;      options,'th'+sclet+'*gsm','databar',0
;      options,strjoin('th'+sc+'_'+'efs_new_gsm'),colors=['g','b','r'],labels=['Ex','Ey','Ez'],LABFLAG=1, ytitle='efs_gsm', ysubtitle='[mV/m]'
;      
;      
;      
      datatype5='peir'

      ;load particle data
      thm_part_load,probe=sc,trange=tr,datatype=datatype5

      ;Note ESA background removal is now enabled by default.
      ;Use esa_bgnd=0 keyword to thm_part_products to disable background removal

      ;Note spacecraft potential must be loaded by the user.  If none is found, none will be used in products generation.
      thm_part_products,probe=sc,datatype=datatype5,trange=tr,outputs='moments', $
        mag_name='th'+sc+'_fgs',pos_name='th'+sc+'_state_pos'
      thm_cotrans,strjoin('th'+sc+'_'+datatype5+'_velocity'),in_coord='dsl',out_coord='gsm',out_suffix='_gsm'
      options,strjoin('th'+sc+'_'+datatype5+'_velocity_gsm'),yrange=[-500,500],colors=['k','b','r'],labels=['Vx','Vy','Vz'],LABFLAG=1,ytitle='Vi gsm'
      options,strjoin('th'+sc+'_'+'fgs_gsm'),colors=['k','b','r'],labels=['Bx','By','Bz'],LABFLAG=1, ytitle='B gsm', ysubtitle='[nT]'

      store_data, strjoin('th'+sc+'_fff_B'),data=['th'+sc+'_fff_32_bsum','th'+sc+'_fce_eq_half','th'+sc+'_fce_eq_01']
      store_data, strjoin('th'+sc+'_fff_E'),data=['th'+sc+'_fff_32_esum','th'+sc+'_fce_eq_half','th'+sc+'_fce_eq_01']
      options, strjoin('th'+sc+'_fff_B th'+sc+'_fff_E'),yrange=[10,4.6e3],/ylog,ystyle=1
      options, strjoin('th'+sc+'_fff_B'),/zlog,zrange=[1e-9,1e-2],ytitle='Frequency', ysubtitle='[Hz]',ztitle='B FFF !c[(nT)!u2!n/Hz]'
      options, strjoin('th'+sc+'_fff_E'),/zlog,zrange=[1e-12,1e-5],ytitle='Frequency', ysubtitle='[Hz]',ztitle='E FFF !c[(V/m)!u2!n/Hz]'
      time_stamp,/off
      
      


      ;tplot,['th'+sc+'_fff_B','th'+sc+'_fff_E'],title='TH-'+strupcase(sc)
      ;stop ;;remove/comment stop to continue

      datatype='peer'

      thm_load_fit,probe=sc,coord='dsl',trange=tr ;this loads the fluxgate mag data(spacecraft generated spin-fit)

      ;load particle data
;      thm_part_load,probe=sc,trange=tr,datatype=datatype
      thm_part_load,probe=sc,trange=tr,datatype='peeb'

      ;;;Note ESA background removal is now enabled by default.
      ;;;Use esa_bgnd=0 keyword to thm_part_products to disable background removal

;      ;;;Note spacecraft potential must be loaded by the user.  If none is found, none will be used in products generation.
;      thm_part_products,probe=sc,datatype=datatype,trange=tr,outputs='energy',$
;        mag_name='th'+sc+'_fgl_dsl',pos_name='th'+sc+'_state_pos',pitch=[0,22.5],suffix='_para',units='eflux';,energy=[10,1000]
;      get_data,'th'+sc+'peer_eflux_energy_para',data=peer_eflux
;      stop
;        
;      thm_part_products,probe=sc,datatype=datatype,trange=tr,outputs='pa',$
;        mag_name='th'+sc+'_fgs',pos_name='th'+sc+'_state_pos',suffix='_fast_10_1keV',units='eflux',energy=[10,1000]
      
      
      
      thm_part_products,probe=sc,datatype='peeb',trange=tr,outputs='energy',mag_name='th'+sc+'_fgl_dsl',pos_name='th'+sc+'_state_pos',pitch=[0,22.5],suffix='_para_1',units='eflux'
      options, strjoin('tha_peeb_eflux_energy_para_1'), ytitle='peeb eflux!c0-22.5deg', ysubtitle='[eV]',/ylog, ztitle='[eV/cm!u2!n-sr-s-eV]',zrange=[1e3,8e7],/zlog,labflag=1
      thm_part_products,probe=sc,datatype='peeb',trange=tr,outputs='energy',mag_name='th'+sc+'_fgl_dsl',pos_name='th'+sc+'_state_pos',pitch=[22.5,45],suffix='_para_2',units='eflux'
      options, strjoin('tha_peeb_eflux_energy_para_2'), ytitle='peeb!c22.5-45deg!ceflux', ysubtitle='[eV]',/ylog, ztitle='[eV/cm!u2!n-sr-s-eV]',zrange=[1e3,8e7],/zlog,labflag=1
      thm_part_products,probe=sc,datatype='peeb',trange=tr,outputs='energy',mag_name='th'+sc+'_fgl_dsl',pos_name='th'+sc+'_state_pos',pitch=[45,67.5],suffix='_para_3',units='eflux'
      options, strjoin('tha_peeb_eflux_energy_para_3'), ytitle='peeb!c45-67.5deg!ceflux', ysubtitle='[eV]',/ylog, ztitle='[eV/cm!u2!n-sr-s-eV]',zrange=[1e3,8e7],/zlog,labflag=1
      thm_part_products,probe=sc,datatype='peeb',trange=tr,outputs='energy',mag_name='th'+sc+'_fgl_dsl',pos_name='th'+sc+'_state_pos',pitch=[67.5,90],suffix='_para_4',units='eflux'
      options, strjoin('tha_peeb_eflux_energy_para_4'), ytitle='peeb!c67.5-90deg!ceflux', ysubtitle='[eV]',/ylog, ztitle='[eV/cm!u2!n-sr-s-eV]',zrange=[1e3,8e7],/zlog,labflag=1
      
      get_data,'tha_peeb_eflux_energy_para_1',data=brst_eflux_1,dlim=dlim
      get_data,'tha_peeb_eflux_energy_para_2',data=brst_eflux_2,dlim=dlim
      get_data,'tha_peeb_eflux_energy_para_3',data=brst_eflux_3,dlim=dlim
      get_data,'tha_peeb_eflux_energy_para_4',data=brst_eflux_4,dlim=dlim
      eflux1=brst_eflux_1.y
      eflux2=brst_eflux_2.y
      eflux3=brst_eflux_3.y
      eflux4=brst_eflux_4.y
      efluxe=brst_eflux_1.v
      efluxt=brst_eflux_1.x
;      save,filename=savepath+'20180325_204828_205408_peeb_padist.sav',efluxt,eflux1,eflux2,eflux3,eflux4,eflux,efluxe
;      stop
      Openr,lun1,strjoin(savepath+'Ew_precflux_Dsdcorrected.txt'),/Get_LUN
      file1=strjoin(savepath+'Ew_precflux_Dsdcorrected.txt')
      nrows1 = File_Lines(file1)
      Ew=Dblarr(nrows1)
      prec=Dblarr(nrows1)
      data=Dblarr(2,nrows1)
      readf,lun1,data
      Ew=data(0,*)
      prec=data(1,*)
      free_lun,lun1
      ind=where(Ew lt 0.5)
      prec[ind]=0.0
      
      store_data,'Ew',data={x:brst_eflux_1.x,y:reform(Ew,n_elements(Ew))}
      store_data,'prec_flux',data={x:brst_eflux_1.x,y:reform(prec,n_elements(prec))}
      options,'Ew',yrange=[0.0,6],ytitle='EHs Ew',ylog=0,ysubtitle='[mV/m]',colors='k',labels=''
      options,'prec_flux',colors='r',ylog=0,yrange=[0,1.0],ytitle='precipitating eflux!c ', ysubtitle='[ergs/cm!u2!n/s]'
      
      
      
      tvectot, 'th'+sc+'_fgl_dsl',tot='Bmag'
      split_vec,'th'+sc+'_fgl_dsl'
      time_clip,'th'+sc+'_fgl_dsl_x',brst_eflux_1.x[0], brst_eflux_1.x[-1],newname='th'+sc+'_fgl_dsl_x_clip'
      time_clip,'th'+sc+'_fgl_dsl_z',brst_eflux_1.x[0], brst_eflux_1.x[-1],newname='th'+sc+'_fgl_dsl_z_clip'
      time_clip,'Bmag',brst_eflux_1.x[0], brst_eflux_1.x[-1],newname='Bmag_clip'
      store_data,'Bm',data=['Bmag','th'+sc+'_fgl_dsl_x_clip','th'+sc+'_fgl_dsl_z_clip']
      options,'Bm',colors=['k','b','r'],yrange=[1,1],ytitle='B',ysubtitle='[nT]',ylog=0,labels=['|B|','Bx','Bz'],labflag=1
      tplot_multiaxis,['tha_peeb_eflux_energy_para_1','Ew'],'prec_flux',2
      tlimit,[brst_eflux_1.x[0], brst_eflux_1.x[-1]]
      tplot_multiaxis,['tha_peeb_eflux_energy_para_1','Ew'],'prec_flux',2
      popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_HGI_peeb_prec_updated',encapsulated=encapsulated,land=landscape,xsize=4.5,ysize=3.5,units=cm
      tplot_multiaxis,['tha_peeb_eflux_energy_para_1','Ew'],'prec_flux',2
      pclose
      stop
      
      
      
      get_data,'tha_peeb_eflux_energy_para_1',data=brst_eflux,dlim=dlim
      brstt=brst_eflux.x
      brsty=brst_eflux.y
      brste=brst_eflux.v
      time_clip,'integrate_tot_2',brst_eflux.x[0],brst_eflux.x[-1],newname='integrate_tot_int'
      get_data,'integrate_tot_int',data=Ew_int
      thEwt=Ew_int.x
      thEw=Ew_int.y
      save,filename=savepath+'THA_Ew_20180325_204828_205408_Ew_50_4k.sav',thEwt,thEw
      stop
;      save,filename=savepath+'20180325_204828_205408_peeb_para.sav',brstt,brsty,brste
      integrate_eflux=make_array(n_elements(brst_eflux.x))
      for i=0,n_elements(brst_eflux.x)-1 do begin
        integrate_eflux[i]=total(brst_eflux.y[i,7:27]*(brst_eflux.v[8:28]-brst_eflux.v[7:27]),/nan)/6.2415e11*0.4783
      endfor
      store_data,'integrate_eflux',data={x:brst_eflux.x,y:integrate_eflux}
      options,'integrate_eflux',colors='r',ylog=0,yrange=[1,1],ytitle='total eflux', ysubtitle='[ergs/cm!u2!n/s]'
      thm_part_products,probe=sc,datatype='peeb',trange=tr,outputs='pa',mag_name='th'+sc+'_fgl_dsl',pos_name='th'+sc+'_state_pos',energy=[50,500],suffix='_50_500eV',units='eflux';,energy=[10,1000]
      options, strjoin('tha_peeb_eflux_pa_50_500eV'), ytitle='peeb!c0.05-0.5keV!ceflux', ysubtitle='[eV]',ylog=0, zrange=[1e5,1e7], ztitle=''  
      thm_part_products,probe=sc,datatype='peeb',trange=tr,outputs='pa',mag_name='th'+sc+'_fgl_dsl',pos_name='th'+sc+'_state_pos',energy=[500,1000],suffix='_500_1keV',units='eflux';,energy=[10,1000]
      options, strjoin('tha_peeb_eflux_pa_500_1keV'), ytitle='peeb!c0.5-1keV!ceflux', ysubtitle='[eV]',ylog=0,zrange=[1e6,5e7], ztitle='' 
      thm_part_products,probe=sc,datatype='peeb',trange=tr,outputs='pa',mag_name='th'+sc+'_fgl_dsl',pos_name='th'+sc+'_state_pos',energy=[1000,5000],suffix='_1k_5keV',units='eflux';,energy=[10,1000]
      options, strjoin('tha_peeb_eflux_pa_1k_5keV'), ytitle='peeb!c1-5keV!ceflux', ysubtitle='[eV]',ylog=0,zrange=[1e6,5e7], ztitle=''   
      thm_part_products,probe=sc,datatype='peeb',trange=tr,outputs='pa',mag_name='th'+sc+'_fgl_dsl',pos_name='th'+sc+'_state_pos',energy=[5000,15000],suffix='_5k_15keV',units='eflux';,energy=[10,1000]
      options, strjoin('tha_peeb_eflux_pa_5k_15keV'), ytitle='peeb!c5-15keV!ceflux', ysubtitle='[eV]',ylog=0,zrange=[1e6,5e7], ztitle=''
      tvectot, 'th'+sc+'_fgl_dsl',tot='Bmag'
      time_clip,'Bmag',time_double('2018-03-25/20:50'), time_double('2018-03-25/20:54'),newname='Bmag_clip'
      options,'Bmag_clip',colors='k',yrange=[1,1],ytitle='|B|',ysubtitle='[nT]',ylog=0
      tplot_multiaxis,['tha_peeb_eflux_energy_para','Bmag_clip','integrate_tot','tha_peeb_eflux_pa_50_500eV','tha_peeb_eflux_pa_500_1keV','tha_peeb_eflux_pa_1k_5keV','tha_peeb_eflux_pa_5k_15keV'],'integrate_eflux',3
;      tlimit,[time_double('2018-03-25/20:50'), time_double('2018-03-25/20:54')]
      popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_HGI_PAD_peeb',encapsulated=encapsulated,land=landscape,xsize=6.0,ysize=5.5,units=cm
      tplot_multiaxis,['tha_peeb_eflux_energy_para','Bmag_clip','integrate_tot','tha_peeb_eflux_pa_50_500eV','tha_peeb_eflux_pa_500_1keV','tha_peeb_eflux_pa_1k_5keV','tha_peeb_eflux_pa_5k_15keV'],'integrate_eflux',3
;      tlimit,[time_double('2018-03-25/20:50'), time_double('2018-03-25/20:54')]
      pclose   
      stop
      
      thm_part_products,probe=sc,datatype=datatype,trange=tr,outputs='energy',$
        mag_name='th'+sc+'_fgs',pos_name='th'+sc+'_state_pos',units='eflux',suffix='_omni';,energy=[10,1000]

      
      get_data,strjoin('th'+sc+'_'+datatype+'_eflux_'+'energy_para'),data=eflux_para
      
      stop
      get_data,strjoin('th'+sc+'_'+datatype+'_eflux_'+'energy_omni'),data=eflux_omni
      options, strjoin('th'+sc+'_'+datatype+'_eflux_'+'energy_omni'), ytitle='peer !comni eflux', ysubtitle='[eV]',/ylog, ztitle='[eV/cm!u2!n-sr-s-eV]',zrange=[1e5,1e9],/zlog,labflag=1
      
      ;time period of 06:12:18 seems to be most intense injection, corresponding to index of 463
      ;      mat_para=eflux_para.y[463,*]  ;31 energy bins
      ;      mat_energy=eflux_para.v  ;31 energy bins
      ;      mat_anti=eflux_anti.y[463,*] ; 31 energy bins
      ;      save,filename=savepath+'eflux_para_anti_20140827.sav',mat_energy,mat_para,mat_anti



      ;8-64eV 9-84eV 10-111eV 11-145eV 12-191eV 13-251eV 14-330eV 15-435eV 16-572eV 17-753eV 18-992eV 19-1304eV 20-1717eV 21-2261eV
      ;22-2977eV 23-3918eV 24-5157eV 25-6789eV 26-8936eV 27-11763eV 28-15485eV 29-20383eV 30-26831eV
      ;for i=8,30,2
;      energyname=['111eV','191eV','330eV','572eV','992eV','1717eV','2977eV','5157eV','8936eV','15485eV']
;
;      tvlct,r,g,b,/get
;      index=indgen(10)+1
;      r[index[0]]=0 & g[index[0]]=0 & b[index[0]]=139
;      r[index[1]]=0 & g[index[1]]=23.2 & b[index[1]]=255
;      r[index[2]]=0 & g[index[2]]=162.3 & b[index[2]]=255
;      r[index[3]]=0 & g[index[3]]=255 & b[index[3]]=209
;      r[index[4]]=0 & g[index[4]]=255 & b[index[4]]=70
;      r[index[5]]=70 & g[index[5]]=255 & b[index[5]]=0
;      r[index[6]]=207 & g[index[6]]=255 & b[index[6]]=0
;      r[index[7]]=255 & g[index[7]]=223 & b[index[7]]=0
;      r[index[8]]=255 & g[index[8]]=174 & b[index[8]]=0
;      r[index[9]]=255 & g[index[9]]=90 & b[index[9]]=0
;      ;     r[index[11]]=255 & g[index[11]]=0 & b[index[11]]=0
;      tvlct, r, g, b
;
;      for i=0,n_elements(energyname)-1 do begin
;        k=i*2+10
;        store_data,strjoin('energy_'+num2str(i)),data={x:eflux_omni.x,y:eflux_omni.y[*,k]}, $
;          limits = {colors:i+1,labels:energyname[i], labflag:1, ylog:1,thick:'1.5'}
;        if i eq 0 then begin
;          lines=strjoin('energy_'+num2str(i))
;        endif else begin
;          lines=[lines,strjoin('energy_'+num2str(i))]
;        endelse
;      endfor
;      store_data,'energy_lines',data=lines
;      options, 'energy_lines', ytitle='peer !comni eflux', ysubtitle='[eV/cm!u2!n-sr-s-eV]',yrange=[1e4,1e8],/ylog,labflag=1



      ;here we examine spacecraft potential data and the corresponding density
      
      ; Get moms from peir and peer correcting sc pot accordingly - this gives omni spectra too!
      ;  have to add esa background removal
      sclet=sc
      thm_load_mom,probe=sclet; (reads L1 file, calibrates on the fly)
      thm_load_esa_pkt,probe=sclet ; reads L0 data, all modes, but is fast
      

;      store_data,'th'+sc+'_peer_en_counts_pot',data='th'+sc+'_peer_en_counts th'+sc+'_pxxm_pot'
;      ;
;      store_data,'th'+sc+'_pexm_density',data='th'+sc+'_peim_density th'+sc+'_peem_density'
;      options,'th'+sc+'_pexm_density','colors',['b','r'],labels=['Ni','Ne'],LABFLAG=1
;      ;
;      ylim,'th'+sc+'_pe?m_density',0.05,5,1
;      ylim,'th'+sc+'_p???_en_counts*',6,30000,1
     

      calc," 'tha_pxxm_pot2' = ('tha_pxxm_pot'+1.)*0.5 "
;      store_data,'th'+sc+'_peer_en_counts_pot2',data='th'+sc+'_peer_en_counts th'+sc+'_pxxm_pot2'
;      store_data,'th'+sc+'_peeb_en_counts_pot2',data='th'+sc+'_peeb_en_counts th'+sc+'_pxxm_pot2'
;      ylim,'th'+sc+'_p???_en_counts*',6,30000,1
;
;;      options,['th'+sc+'_peer_en_counts_pot','th'+sc+'_peer_en_counts_pot2','th'+sc+'_peeb_en_counts_pot2'],yrange=[7,200],/ylog
;;      tplot,['th'+sc+'_peer_en_counts_pot','th'+sc+'_peer_en_counts_pot2','th'+sc+'_peeb_en_counts_pot2']
;;      stop
;      
      thm_part_products,probe=sclet,datatype='peer',outputs='moments energy',sc_pot_name='th'+sc+'_pxxm_pot2',mag_name='th'+sclet+'_fgs',pos_name='th'+sclet+'_state_pos',esa_bgnd_remove=1
;      store_data,'th'+sc+'_peer_eflux_energy_pot',data='th'+sc+'_peer_eflux_energy th'+sc+'_pxxm_pot2'
;      ylim,'th'+sc+'_pee?_eflux_e*',6,30000,1
;      thm_part_products,probe=sclet,datatype='peir',outputs='moments energy',sc_pot_name='th'+sc+'_pxxm_pot2'
;      ylim,'th'+sc+'_pei?_eflux_e*',6,25000,1
;      store_data,'th'+sc+'_pexr_density',data='th'+sc+'_peir_density th'+sc+'_peer_density'
;      ylim,'th'+sc+'_pe?r_density',0.05,2,1
;      options,'th'+sc+'_pexr_density','colors',['b','r']
;      options,'th'+sc+'_pexr_density','labels',['Ni','Ne']
;      options,'th'+sc+'_pexr_density','LAGFLAG',1
;      options,'*','databar',0.
      
      esa_datatype = 'peer'
      sst_datatype = 'psef'

;This will automatically load the required particle data and interpolate
;between the two instruments to produce the combined product.
;The original ESA and SST data will be passed out via the last two keywords.
combined = thm_part_combine(probe=sc, trange=tr, $
  esa_datatype=esa_datatype, sst_datatype=sst_datatype, $
  orig_esa=esa, orig_sst=sst)

;Default midpoint energies for interpolation above ESA energies are(in eV):
;Ions: [25000,26000.,28000.,30000.0,34000.0,41000.0,53000.0,67400.0,95400.0,142600.,207400.,297000.,421800.,654600.,1.13460e+06,2.32980e+06,4.00500e+06]
;Electrons: [27000,28000.,29000.,30000.0, 31000.0,41000.0,52000.0,65500.0,93000.0,139000.,203500.,293000.,408000.,561500.,719500.]

thm_part_products, dist_array=combined, outputs='energy'
options,'th'+sc+'_pterf_eflux_energy',ytitle='pterf !comni eflux!cenergy',ysubtitle='[eV]',yrange=[6,2e5],ztitle='[eV/cm!u2!n-sr-s-eV]',zrange=[1e3,8e7],/zlog,labflag=1
      
      
      
      get_data,'th'+sc+'_peer_avgtemp',data=Te
;      get_data,'th'+sc+'_peir_avgtemp',data=Ti
      
      
      ;Vsc
      get_data, 'th'+sclet+'_pxxm_pot2', data = Vdata, index = index
      Vdata.y[*] = -Vdata.y[*]
      store_data, 'th'+sclet+'_pxxm_-sc_pot', data = Vdata
      
   
      
      
      ; Now add the scpot derived density based on Lecture 13 crib:
      
      ;thm_scpot2dens_opt_n,probe=sclet,data='peer' ; you can modify this to input Te, scpot yourself... but here use as is
      
      ;SM position
      
      thm_load_state, probe = sclet,trange = tr, coord = 'sm', suffix = '_sm_coordinates'
      
      get_data, 'th'+sclet+'_state_pos_sm_coordinates', data = tmp
      
      store_data, 'th'+sclet+'_R', data = {x:tmp.x, y:sqrt(tmp.y[*, 0]^2+tmp.y[*, 1]^2+tmp.y[*, 2]^2)/6371.2}, $
        dlim = {colors:['k'], labels:['R'], ysubtitle:'[km]', labflag:1, constant:0, ytitle:'th'+sclet+'_R'}
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;  ;load data
      ;  thm_load_esa, probe = sclet, level = 2, datatype = 'peer_vthermal'
        ;Vth
        get_data, 'th'+sclet+'_peer_vthermal', data = Vthdata, index = index
      
      
      
        ;R
        tinterpol_mxn, 'th'+sclet+'_R', 'th'+sclet+'_pxxm_pot2', newname = 'th'+sclet+'_R_int'
        get_data, 'th'+sclet+'_R_int', data = R, index = index1
      
      
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;Vth correction
        factor=1.0
        Vthdata.y = (Vthdata.y*(1-(-1/(1+exp(-2*(R.y-4)))+1))+10^4*(-1/(1+exp(-2*(R.y-4)))+1))*factor
        store_data, 'th'+sclet+'_vthermal2', data = Vthdata
      
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;scpot-density conversion
      ;  Case sclet Of
      ;    'a': Ne_scpot = (10^(Vdata.y/22.0)*20000.+10^(Vdata.y/ 5.0)*160000.+10^(Vdata.y/2.0)*15000000.0+10^(Vdata.y/0.4)*1500000000000.0)/Vthdata.y
      ;    'b': Ne_scpot = (10^(Vdata.y/25.5)*16500.                          +10^(Vdata.y/2.0)*3000000.0+10^(Vdata.y/0.4)*1500000000000.0)/Vthdata.y
      ;    'c': Ne_scpot = (10^(Vdata.y/26.0)*22000.                          +10^(Vdata.y/2.0)*3000000.0+10^(Vdata.y/0.25)*1e14)/Vthdata.y
      ;    'd': Ne_scpot = (10^(Vdata.y/25.5)*20000.+10^(Vdata.y/ 5.0)*30000.+10^(Vdata.y/2.0)*10000000.0+10^(Vdata.y/0.2)*5000000000000000.0)/Vthdata.y
      ;    'e': Ne_scpot = (10^(Vdata.y/25.5)*20000.+10^(Vdata.y/ 5.0)*30000.+10^(Vdata.y/2.0)*10000000.0+10^(Vdata.y/0.2)*3000000000000000.0)/Vthdata.y
      ;  Endcase
      
      Ne_scpot = (10^(Vdata.y/22.0)*20000.+10^(Vdata.y/ 5.0)*160000.+10^(Vdata.y/2.0)*15000000.0+10^(Vdata.y/0.4)*1500000000000.0)/Vthdata.y
      store_data, 'th'+sclet+'_peer_density_npot', data = {x:Vdata.x, y:Ne_scpot}, dlim = {colors:['k'], labels:['Ne_scpot'], ysubtitle:'[cm-3]', labflag:1, constant:0, ylog:1}
      
;      store_data,'th'+sclet+'_pexr_density_pot',data='th'+sc+'_peir_density th'+sc+'_peer_density th'+sc+'_peer_density_npot'
;      ylim,'th'+sc+'_pe?r_dens*',0.01,1,1
;      options,'th'+sc+'_pexr_density*','colors',['b','r','k']
;      options,'th'+sc+'_pexr_density*','labels',['Ni','Ne','Ne_scpot']
;      options,'th'+sc+'_pexr_density*','labflag',1
;;      tplot,['th'+sc+'_pexr_density_pot','th'+sc+'_peer_eflux_energy_pot']
;;      
;;      tplot_apply_databar
;;      stop
      ;
      fp_over_fce=8989.0*sqrt(Ne_scpot)/fce_eq
      store_data,'fpefce',data={x:Vdata.x, y:fp_over_fce}
      options,'fpefce',colors=['k'],ytitle='fpe/fce',yrange=[1,100],/ylog
      options,'th'+sclet+'_peer_density_npot',yrange=[0.1,0.4],ylog=0
;      
;      
;      store_data,'Pressure_plasma',data={x:Vdata.x,y:1.38*1.16*Ne_scpot*(Te.y+Ti.y)*1e-4}
;      store_data,'Pressures',data=['Pressure_B','Pressure_plasma']
;      options,'Pressures',colors=['b','r'],labels=['magnetic P','plasma P'],LABFLAG=1, ytitle='Pressures', ysubtitle='[nPa]'
;;      
;;      
;;
;;      
;;
;;      tplot, ['th'+sc+'_'+'fgl_gsm','th'+sc+'_'+'efs_new_gsm',strjoin('th'+sc+'_'+datatype5+'_velocity'),'th'+sclet+'_peer_density_npot','Pressures','energy_lines','thx_integrate_e_clip_new'],var_label=['Lat','MLT','L'],title='TH-'+strupcase(sc)
;;      tlimit,[time_double('2018-03-25/20:50'),time_double('2018-03-25/20:55')]
;;      popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_closer',encapsulated=encapsulated,land=landscape,xsize=8,ysize=10.5,units=cm
;;      tplot, ['th'+sc+'_'+'fgl_gsm','th'+sc+'_'+'efs_new_gsm',strjoin('th'+sc+'_'+datatype5+'_velocity'),'th'+sclet+'_peer_density_npot','Pressures','energy_lines','thx_integrate_e_clip_new'],var_label=['Lat','MLT','L'],title='TH-'+strupcase(sc)
;;      tlimit,[time_double('2018-03-25/20:50'),time_double('2018-03-25/20:55')]
;;      pclose
;;      stop
;      
      tplot, ['th'+sc+'_'+'fgs_gsm','th'+sc+'_'+datatype5+'_velocity_gsm','fpefce','th'+sc+'_pterf_eflux_energy','th'+sc+'_fff_B','th'+sc+'_fff_E','integrate_tot'],var_label=['Lat','MLT','R'],title='TH-'+strupcase(sc)
      popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_HGI',encapsulated=encapsulated,land=landscape,xsize=4.0,ysize=4.5,units=cm

      tplot, ['th'+sc+'_'+'fgs_gsm','th'+sc+'_'+datatype5+'_velocity_gsm','fpefce','th'+sc+'_pterf_eflux_energy','th'+sc+'_fff_B','th'+sc+'_fff_E','integrate_tot'],var_label=['Lat','MLT','R'],title='TH-'+strupcase(sc)

      pclose
      stop

      tr2 = [time_double('2018-03-25/20:40'), time_double('2018-03-25/21:00')]
      timespan, tr2
      thm_load_state, probe = sc, trange = tr2,/get_support_data
      thm_load_fgm, probe = sc, trange = tr2, level=2, coord='dsl', data='fgs'
      thm_load_fit,probe=sc,coord='dsl',trange=tr2
      ;;load the wave burst mode E-field data
      thm_efi_clean_efw,probe=sc, trange=tr2
      ;;load the wave burst mode B-field data
      thm_load_scm,level=1,/get_suppport_data,probe=sc, trange=tr2,datatype='scw'
      bwtrange=thm_jbt_get_btrange('th'+sc+'_scw', nb=nbw, tind=tindw)
      ;burst data from 06:11:34 to 06:11:42 in 8 seconds
      nbt=1
      

      if (nbw ge 1) and (bwtrange[0] ne -1) then begin
        ;clip to the burst
        time_clip, strjoin('th'+sc+'_efw_clean_dsl'),bwtrange[0,nbt],bwtrange[1,nbt],newname=strjoin('th'+sc+'_efw_clip')
        time_clip, strjoin('th'+sc+'_scw'),bwtrange[0, nbt],bwtrange[1, nbt],newname=strjoin('th'+sc+'_scw_clip')
        tsmooth2, strjoin('th'+sc+'_fgs_dsl'),20,newname=strjoin('th'+sc+'_fgs_dsl_sm')
        time_clip, strjoin('th'+sc+'_fgs_dsl_sm'),bwtrange[0,nbt],bwtrange[1,nbt],newname=strjoin('th'+sc+'_fgs_dsl_clip')
        time_clip, strjoin('th'+sc+'_state_pos'),tr2[0],tr2[1],newname=strjoin('th'+sc+'_state_pos_clip')
        tinterpol_mxn, strjoin('th'+sc+'_state_pos_clip'),strjoin('th'+sc+'_fgs_dsl_clip'), newname = 'th'+sc+'_state_pos_int'



        thm_fac_matrix_make, strjoin('th'+sc+'_fgs_dsl_clip'), pos_var_name=strjoin('th'+sc+'_state_pos_int'),other_dim='Rgeo', $
          newname=strjoin('th'+sc+'_fgs_fac_mat')
        
        tvector_rotate, strjoin('th'+sc+'_fgs_fac_mat'), strjoin('th'+sc+'_efw_clip'), $
          newname=strjoin('th'+sc+'_efw_full_fac'), error=error1
        options,strjoin('th'+sc+'_efw_full_fac'),colors=['k','b','r'],ytitle='efw fac !c0-8192 Hz',ysubtitle='[mV/m]',labels=['Ex','Ey','Ez'],LABFLAG=1
        time_clip, strjoin('th'+sc+'_efw_full_fac'),time_double('2018-03-25/20:52:32'),time_double('2018-03-25/20:52:35'),newname=strjoin('th'+sc+'_efw_full_fac_clip')
        get_data,strjoin('th'+sc+'_efw_full_fac_clip'),data=efw_full
        Ewvt=efw_full.y
        store_data,'E_tot',data={x:efw_full.x,y:sqrt(efw_full.y[*,0]^2+efw_full.y[*,1]^2+efw_full.y[*,2]^2)}
        store_data,'E_perp',data={x:efw_full.x,y:sqrt(efw_full.y[*,0]^2+efw_full.y[*,1]^2)}
        store_data,'E_para',data={x:efw_full.x,y:efw_full.y[*,2]}
        wav_data,'E_perp',trange=['2018-03-25/20:52:32','2018-03-25/20:52:35'],maxpoints=2.0^17 
        wav_data,'E_para',trange=['2018-03-25/20:52:32','2018-03-25/20:52:35'],maxpoints=2.0^17
        wav_data,'E_tot',trange=['2018-03-25/20:52:32','2018-03-25/20:52:35'],maxpoints=2.0^17
        options,strjoin('E_perp_wv_pow'),zrange=[5e-4,5e-1],/zlog,ytitle='th'+sc+' !cefw Eperp!u2!n !cwavelet',ysubtitle='[Hz]',ztitle=['[(mV)^2/Hz]'],LABFLAG=1 
        options,strjoin('E_para_wv_pow'),zrange=[1e-3,5e-1],/zlog,ytitle='efw Ez',ysubtitle='[Hz]',ztitle=['[(mV)^2/Hz]'],LABFLAG=1 
        options,strjoin('E_tot_wv_pow'),zrange=[1e-3,5e-1],/zlog,ytitle='th'+sc+' !cefw |E|!u2!n !cwavelet',ysubtitle='[Hz]',ztitle=['[(mV)^2/Hz]'],LABFLAG=1
 
        
        
        tvector_rotate, strjoin('th'+sc+'_fgs_fac_mat'), strjoin('th'+sc+'_scw_clip'), $
          newname=strjoin('th'+sc+'_scw_full_fac'), error=error2
        options,strjoin('th'+sc+'_scw_full_fac'),colors=['k','b','r'],ytitle='tha !cscw fac !c0-4096 Hz',ysubtitle='[nT]',labels=['Bx','By','Bz'],LABFLAG=1
         time_clip, strjoin('th'+sc+'_scw_full_fac'),time_double('2018-03-25/20:52:32'),time_double('2018-03-25/20:52:35'),newname=strjoin('th'+sc+'_scw_full_fac_clip')
        get_data,strjoin('th'+sc+'_scw_full_fac_clip'),data=scw_full
        Bwvt=scw_full.y
        store_data,'B_y',data={x:scw_full.x,y:scw_full.y[*,1]}
        store_data,'B_z',data={x:scw_full.x,y:scw_full.y[*,2]}
        store_data,'B_tot',data={x:scw_full.x,y:sqrt(scw_full.y[*,0]^2+scw_full.y[*,1]^2+scw_full.y[*,2]^2)}
        wav_data,'B_y',trange=['2018-03-25/20:52:32','2018-03-25/20:52:35'],maxpoints=2.0^17
        wav_data,'B_z',trange=['2018-03-25/20:52:32','2018-03-25/20:52:35'],maxpoints=2.0^17
;        wav_data,'B_tot',trange=['2018-03-25/20:52:32','2018-03-25/20:52:35'],maxpoints=2.0^17
        options,strjoin('B_y_wv_pow'),zrange=[5.e-8,1e-5],/zlog,ytitle='scw By',ysubtitle='[Hz]',ztitle=['[nT^2/Hz]'],LABFLAG=1
        options,strjoin('B_z_wv_pow'),zrange=[5.e-8,1e-5],/zlog,ytitle='scw Bz',ysubtitle='[Hz]',ztitle=['[nT^2/Hz]'],LABFLAG=1
;        options,strjoin('B_tot_wv_pow'),zrange=[5.e-8,1e-5],/zlog,ytitle='th'+sc+' !cscw |B|!u2!n !cwavelet',ysubtitle='[Hz]',ztitle=['[nT^2/Hz]'],LABFLAG=1
        
          
        get_data,strjoin('th'+sc+'_efw_clip'),data=efw_clip,limits=lim,dlimits=dlim
        Ew=efw_clip.y
        nfft=1024
        filtere=digital_filter(5./8192.,4000.0/8192.,50,nfft*2,/double)     ;5-4096 Hz
        for i=0,2 do Ew[*,i]=convol(Ew[*,i],filtere,/center)
        store_data,strjoin('th'+sc+'_efw_dsl_hp'),data={x:efw_clip.x,y:Ew},limits=lim,dlimits=dlim
        

;        tsmooth2,strjoin('th'+sc+'_efw_clip'),1024.0,newname=strjoin('th'+sc+'_efw_clip_sm')
;        dif_data,'th'+sc+'_efw_clip','th'+sc+'_efw_clip_sm',newname=strjoin('th'+sc+'_efw_dsl_hp')
        ;       tinterpol_mxn, strjoin('th'+sc+'_fgs_fac_mat'),strjoin('th'+sc+'_efw_dsl_hp'), newname = 'th'+sc+'_fgs_fac_mat_int'
        tvector_rotate, strjoin('th'+sc+'_fgs_fac_mat'), strjoin('th'+sc+'_efw_dsl_hp'), $
          newname=strjoin('th'+sc+'_efw_fac'), error=error1
        options,strjoin('th'+sc+'_efw_fac'),colors=['k','b','r'],ytitle='efw fac !c5-4000 Hz',ysubtitle='[mV/m]',labels=['Ex','Ey','Ez'],LABFLAG=1
        ;    tvector_rotate, strjoin('th'+probes+'_fgs_fac_mat'), strjoin('th'+probes+'_efw_clip_sm'), $
        ;      newname=strjoin('th'+probes+'_efw_clip_sm_fac'), error=error3

        
        
        get_data,strjoin('th'+sc+'_scw_clip'),data=scw_clip,limits=lim,dlimits=dlim
        Bw=scw_clip.y
        filterb=digital_filter(5./4096,4000.0/4096,50,nfft,/double)
        for i=0,2 do Bw[*,i]=convol(Bw[*,i],filterb,/center)
        store_data,strjoin('th'+sc+'_scw_dsl_hp'),data={x:scw_clip.x,y:Bw},limits=lim,dlimits=dlim
        
;        tsmooth2,strjoin('th'+sc+'_scw_clip'),1024.0,newname=strjoin('th'+sc+'_scw_clip_sm')
;        dif_data,'th'+sc+'_scw_clip','th'+sc+'_scw_clip_sm',newname=strjoin('th'+sc+'_scw_dsl_hp')
        tvector_rotate, strjoin('th'+sc+'_fgs_fac_mat'), strjoin('th'+sc+'_scw_dsl_hp'), $
          newname=strjoin('th'+sc+'_scw_fac'), error=error2
        options,strjoin('th'+sc+'_scw_fac'),colors=['k','b','r'],ytitle='scw fac!c5-4000 Hz',ysubtitle='[nT]',labels=['Bx','By','Bz'],LABFLAG=1
        ;    tvector_rotate, strjoin('th'+probes+'_fgs_fac_mat'), strjoin('th'+probes+'_scw_clip_sm'), $
        ;      newname=strjoin('th'+probes+'_scw_clip_sm_fac'), error=error4

        ;      tplot,[strjoin('th'+sc+'_efw_fac'),strjoin('th'+sc+'_scw_fac')]
        ;      tlimit,[bwtrange[0,nbt],bwtrange[1,nbt]]
        ;      stop

;        ;; do the wave polarization analysis
;        twavpol,strjoin('th'+sc+'_efw_fac'),nopfft=4096,steplength=floor(2048.0)
;        options,strjoin('th'+sc+'_efw_fac_powspec'),spec=1,yrange=[10.0,2000.],/ylog,/zlog,ystyle=1
;
;
;        twavpol,strjoin('th'+sc+'_scw_fac'),nopfft=4096.0,steplength=floor(2048.0)
;        options,strjoin('th'+sc+'_scw_fac_powspec'),spec=1,yrange=[10.0,2000.],/ylog,/zlog,ystyle=1


        

;        ; Get clipped time data for dispersion analysis in Matlab
;        time_clip, strjoin('th'+sc+'_efw_fac'),time_double('2018-03-25/20:52:32'),time_double('2018-03-25/20:52:35'),newname=strjoin('th'+sc+'_efw_fac_clip')
        get_data,'th'+sc+'_efw_full_fac_clip',data=efw_clip_fac
        store_data,'efw_perp',data={x:efw_clip_fac.x,y:sqrt(efw_clip_fac.y[*,0]^2+efw_clip_fac.y[*,1]^2)}
        get_data,'efw_perp',data=efw_perp
;        time_clip, strjoin('th'+sc+'_scw_fac'),time_double('2018-03-25/20:52:32'),time_double('2018-03-25/20:52:35'),newname=strjoin('th'+sc+'_scw_fac_clip')
        get_data,'th'+sc+'_scw_full_fac_clip',data=scw_clip_fac
        store_data,'scw_perp',data={x:scw_clip_fac.x,y:sqrt(scw_clip_fac.y[*,0]^2+scw_clip_fac.y[*,1]^2)}
        get_data,'scw_perp',data=scw_perp
;        perp_E=efw_perp.y
;        perp_B=scw_perp.y
;        save,filename=savepath+'Wave_E_B_20180325_205232_205235.sav',perp_E,perp_B
;        stop
        
        
        
        ; Get data
        get_data,'th'+sc+'_efw_fac',lim=lime,dlim=dlime,data=efw
        get_data,'th'+sc+'_scw_fac',lim=limb,dlim=dlimb,data=scw


        ; Calculate bandpass filtered, time domain efw
        ;ndata = n_elements(efw.x)       ;jmm, 22-may-2012 efw was interpolated earlier
        E=efw.y
        filter=digital_filter(100./8192,4000./8192,50,nfft,/double)
        for i=0,2 do E[*,i]=convol(E[*,i],filter,/center)
        store_data,strjoin('th'+sc+'_efw_fac_bandpass'),data={x:efw.x,y:E},limits=lime,dlimits=dlim3  
        options,strjoin('th'+sc+'_efw_fac_bandpass'),colors=['k','b','r'],ytitle='efw fac !c0.1-4.0 kHz',ysubtitle='[mV/m]',labels=['Ex','Ey','Ez'],LABFLAG=1

        E_low=efw.y
        filter=digital_filter(5./8192,100./8192,50,nfft,/double)
        for i=0,2 do E_low[*,i]=convol(E_low[*,i],filter,/center)
        store_data,strjoin('th'+sc+'_efw_fac_bandpass_low'),data={x:efw.x,y:E_low},limits=lime,dlimits=dlime
        options,strjoin('th'+sc+'_efw_fac_bandpass_low'),colors=['k','b','r'],ytitle='efw fac !c5-100 Hz',ysubtitle='[mV/m]',labels=['Ex','Ey','Ez'],LABFLAG=1
        
        
        store_data,strjoin('th'+sc+'_efw_fac_x'),data={x:efw.x,y:E[*,0]}
        options,strjoin('th'+sc+'_efw_fac_x'),colors=['k'],ytitle='th'+sc+' !cefw fac x!c0.1-4 kHz',ysubtitle='[mV/m]',labels=['Ex'],LABFLAG=1
        store_data,strjoin('th'+sc+'_efw_fac_y'),data={x:efw.x,y:E[*,1]}
        options,strjoin('th'+sc+'_efw_fac_y'),colors=['b'],ytitle='th'+sc+' !cefw fac y!c0.1-4 kHz',ysubtitle='[mV/m]',labels=['Ey'],LABFLAG=1
        store_data,strjoin('th'+sc+'_efw_fac_z'),data={x:efw.x,y:E[*,2]}
        options,strjoin('th'+sc+'_efw_fac_z'),colors=['r'],ytitle='th'+sc+' !cefw fac z!c0.1-4 kHz',ysubtitle='[mV/m]',labels=['Ez'],LABFLAG=1
        E_r=dblarr(n_elements(E[*,2]))
        for i=0,n_elements(E[*,2])-1 do begin
          if abs(E[i,2]) gt 4.0 and abs(E[i,0]^2+E[i,1]^2) gt 1.0 then begin
            E_r[i]=abs(E[i,2])/sqrt(E[i,0]^2+E[i,1]^2)
          endif else begin
            E_r[i]=!values.f_nan
          endelse
        endfor
        store_data,strjoin('E_r'),data={x:efw.x,y:E_r}
        options,strjoin('E_r'),colors=['k'],ytitle='Epara/Eperp'

;        tplot,[strjoin('th'+sc+'_scw_fac_powspec'),strjoin('th'+sc+'_efw_fac_powspec'),strjoin('th'+sc+'_scw_fac'),strjoin('th'+sc+'_efw_fac'),strjoin('th'+sc+'_efw_fac_bandpass'),strjoin('th'+sc+'_efw_fac_bandpass_low')]
;        tlimit,[bwtrange[0,nbt],bwtrange[1,nbt]]
;        popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_burst_1',encapsulated=encapsulated,land=landscape,xsize=8,ysize=7.5,units=cm
;        tplot,[strjoin('th'+sc+'_scw_fac_powspec'),strjoin('th'+sc+'_efw_fac_powspec'),strjoin('th'+sc+'_scw_fac'),strjoin('th'+sc+'_efw_fac'),strjoin('th'+sc+'_efw_fac_bandpass'),strjoin('th'+sc+'_efw_fac_bandpass_low')]
;        tlimit,[bwtrange[0,nbt],bwtrange[1,nbt]]
;        pclose
;        stop
        
;
;                tplot,[strjoin('th'+sc+'_scw_fac'),strjoin('th'+sc+'_efw_fac'),strjoin('th'+sc+'_efw_fac_bandpass'),strjoin('th'+sc+'_efw_fac_bandpass_low')]
;                tlimit,[time_double('2018-03-25/20:52:27.0'), time_double('2018-03-25/20:52:27.5')]
;                popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_burst_1_zoom',encapsulated=encapsulated,land=landscape,xsize=7,ysize=6.5,units=cm
;                tplot,[strjoin('th'+sc+'_scw_fac'),strjoin('th'+sc+'_efw_fac'),strjoin('th'+sc+'_efw_fac_bandpass'),strjoin('th'+sc+'_efw_fac_bandpass_low')]
;                tlimit,[time_double('2018-03-25/20:52:27.0'), time_double('2018-03-25/20:52:27.5')]
;                pclose

;                tplot,[strjoin('th'+sc+'_efw_full_fac'),strjoin('E_perp_wv_pow'),strjoin('E_para_wv_pow'),strjoin('B_perp_wv_pow'),strjoin('B_para_wv_pow')]
;                tlimit,[time_double('2018-03-25/20:52:33.5'), time_double('2018-03-25/20:52:34.0')]
;                popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_burst_2_zoom_wavelet',encapsulated=encapsulated,land=landscape,xsize=7,ysize=6.5,units=cm
;                tplot,[strjoin('th'+sc+'_efw_full_fac'),strjoin('E_perp_wv_pow'),strjoin('E_para_wv_pow'),strjoin('B_perp_wv_pow'),strjoin('B_para_wv_pow')]
;                tlimit,[time_double('2018-03-25/20:52:33.5'), time_double('2018-03-25/20:52:34.0')]
;                pclose
;                stop
                
                tplot,[strjoin('E_para_wv_pow'),strjoin('B_z_wv_pow'),strjoin('B_y_wv_pow')]
                tlimit,[time_double('2018-03-25/20:52:33.5'), time_double('2018-03-25/20:52:34.0')]
                popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_burst_2_zoom_wavelet',encapsulated=encapsulated,land=landscape,xsize=4.5,ysize=3,units=cm
                tplot,[strjoin('E_para_wv_pow'),strjoin('B_z_wv_pow'),strjoin('B_y_wv_pow')]
                tlimit,[time_double('2018-03-25/20:52:33.5'), time_double('2018-03-25/20:52:34.0')]
                pclose
                stop

                tplot,[strjoin('th'+sc+'_scw_fac'),strjoin('th'+sc+'_efw_fac'),strjoin('th'+sc+'_efw_fac_bandpass_low'),strjoin('th'+sc+'_efw_fac_bandpass')]
                tlimit,[time_double('2018-03-25/20:52:33.5'), time_double('2018-03-25/20:52:34.0')]
                popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_burst_2_zoom',encapsulated=encapsulated,land=landscape,xsize=4.5,ysize=4.0,units=cm
                tplot,[strjoin('th'+sc+'_scw_fac'),strjoin('th'+sc+'_efw_fac'),strjoin('th'+sc+'_efw_fac_bandpass_low'),strjoin('th'+sc+'_efw_fac_bandpass')]
                tlimit,[time_double('2018-03-25/20:52:33.5'), time_double('2018-03-25/20:52:34.0')]
                pclose
                stop


        
;        tplot,[strjoin('th'+sc+'_efw_fac_x'),strjoin('th'+sc+'_efw_fac_y'),strjoin('th'+sc+'_efw_fac_z'),strjoin('E_r')]
;        tlimit,[time_double('2018-03-25/20:52:27.0'), time_double('2018-03-25/20:52:27.5')]
;        popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_burst_1_zoom_EH',encapsulated=encapsulated,land=landscape,xsize=7,ysize=6.5,units=cm
;        tplot,[strjoin('th'+sc+'_efw_fac_x'),strjoin('th'+sc+'_efw_fac_y'),strjoin('th'+sc+'_efw_fac_z'),strjoin('E_r')]
;        tlimit,[time_double('2018-03-25/20:52:27.0'), time_double('2018-03-25/20:52:27.5')]
;        pclose 
;        stop

;        tplot,[strjoin('th'+sc+'_efw_fac_x'),strjoin('th'+sc+'_efw_fac_y'),strjoin('th'+sc+'_efw_fac_z'),strjoin('E_r')]
;        tlimit,[time_double('2018-03-25/20:52:33.5'), time_double('2018-03-25/20:52:34.0')]
;        popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_burst_2_zoom_EH',encapsulated=encapsulated,land=landscape,xsize=7,ysize=6.5,units=cm
;        tplot,[strjoin('th'+sc+'_efw_fac_x'),strjoin('th'+sc+'_efw_fac_y'),strjoin('th'+sc+'_efw_fac_z'),strjoin('E_r')]
;        tlimit,[time_double('2018-03-25/20:52:33.5'), time_double('2018-03-25/20:52:34.0')]
;        pclose
        tplot,[strjoin('th'+sc+'_efw_fac_bandpass')]
        tlimit,[time_double('2018-03-25/20:52:33.768'), time_double('2018-03-25/20:52:33.770')]
        popen,savepath+'THEMIS'+strupcase(sc)+'_20180325_burst_2_zoom_EH_blow',encapsulated=encapsulated,xsize=5,ysize=2.8,units=cm
        tplot,[strjoin('th'+sc+'_efw_fac_bandpass')]
        tlimit,[time_double('2018-03-25/20:52:33.768'), time_double('2018-03-25/20:52:33.770')]
        pclose
        stop
      endif


    endfor ;;loop for each probe

  endfor ;;loop for different events
  ;stop
end