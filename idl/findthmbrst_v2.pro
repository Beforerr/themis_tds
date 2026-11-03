;+
  ; :Author: qlpdl
  ;-
pro findthmbrst_v2
  
  thm_init       ; initiate THEMIS data
  
  ;;;;;;;;;;;;;;;;;;;;;new code setups to make the plot better
  device,decomposed=0
  loadct,39
  !p.background=255
  !p.color=0
  !p.charsize = 1.3
  !p.thick = 1.0
  tplot_options,'xmargin',[25.,12.]
  tplot_options,'ymargin',[3,3]
  tplot_options,'xticklen',0.05
  tplot_options,'yticklen',0.03
  tplot_options,'xthick',2
  tplot_options,'ythick',2
  tplot_options,'ygap',0.25   ;eliminate the vertical spacing between panels
  ;;;;;;;;;;;;;;;;;;;;;
  
  ;start time to examine
  time1='2017-08-31/00:00'
  ;end time to examine , my own preference is to check per day, you can adjust this time span
  time2='2017-09-01/00:00'
  CD, '\Users\qlpdl\OneDrive\Desktop\Spedas\SpedasWork' ; this folder needs to be created beforehand by you
    openw,1,strjoin(strsplit(strmid(time1,0,10),'-',/EXTRACT))+'.txt'     ;open a file with time as name to write data into
  tr2 = [time_double(time1), time_double(time2)]    ; set the numeric time range that functions use
  timespan, tr2     ; set up time range for analysis
  ;;load the wave burst mode E-field data0
  sc='a'; three spacecraft 'a', 'd', 'e' in the near-Earth magnetosphere, 'b', 'c' are in lunar orbit
  nfft=1024     ;size of fft window
  thm_efi_clean_efw,probe=sc, trange=tr2     ;preparing wave burst electric fields data
  ewtrange=thm_jbt_get_btrange('th'+sc+'_efw', nb=nbw, tind=tindw)   ;get burst collection data during the entire intervals
  for i=0,nbw-1 do begin       ;loop over burst data collections
    if ewtrange[0,i] lt time_double(time1) or ewtrange[1,i] gt time_double(time2) then continue    ; bad data intervals
    if (ewtrange[1,i]-ewtrange[0,i]) lt 1.0 then continue    ; empty data collection
    time_clip,'th'+sc+'_efw_clean_fac',ewtrange[0,i],ewtrange[1,i],newname=strjoin('th'+sc+'_efw_clip') ; clip data to specific burst times
    get_data,strjoin('th'+sc+'_efw_clip'),data=efw_clip
    efw_clip_len_check=efw_clip.y[*,2]
    ind_len=where(finite(efw_clip_len_check))
    if n_elements(ind_len) lt 1024 then continue
    options,'th'+sc+'_efw_clip',colors=['k','b','r'],labels=['Ex','Ey','Ez'],labflag=1    ; set up plotting keywords
    wav_data,'th'+sc+'_efw_clip',trange=[time_string(ewtrange[0,i]),time_string(ewtrange[1,i])],maxpoints=2.0^18     ;wavelet analyzing the burst data
    options,strjoin('th'+sc+'_efw_clip_wv_pow'),zrange=[1e-4,5e-1],/zlog,ytitle='Ez wavelet',ysubtitle='[Hz]',ztitle=['[(mV/m)^2/Hz]'],LABFLAG=1 ; set options of wavelet spectra for plotting
    tplot,['th'+sc+'_efw_clip','th'+sc+'_efw_clip_wv_pow']  ;plot both time series and wavelet time-frequency spectra
    tlimit,[ewtrange[0,i],ewtrange[1,i]]     ;set/limit a time window to display
        
    tlimit   ;here is what you can do with cursor (one click for start and one click for end on the plot) to zoom in details of the spectra
    ;tlimit,/full   can be used to zoom out the entire data period
    
    check:  ;here is a loop

 stop        ; stop here to decide if data is worth looking into and save
    read, B, prompt='enter:'         ; this is interactive, type 1 to select the event, during selection process, if error occurred, type -1 to do again, type anything else (e.g., 0) to jump out to proceed loop over i
    case 1 of
      strcmp(B,'1'): begin          ;if 1, then look into data, analyze event, and save TDS amplitudes, time, MLT, R; this needs to be updated for additional parameters
        tcur=timerange(/current)    ;get time span of the data of interest
        makepng,strjoin(STRSPLIT(time_string(ewtrange[0,i]), '-/:', /EXTRACT))+'_'+strjoin(STRSPLIT(time_string(ewtrange[1,i]), '-/:', /EXTRACT)) ; save png figures for record
        get_data,'th'+sc+'_efw_clip',data=efw_z     ; time series data, efw_z.x is time (use time_string to check time data, time_double can be used to convert back to numerical, efw_z.y is 3-dimension electric field data
        f=0.5/(efw_z.x[1]-efw_z.x[0])      ; Nyquist frequency of the data (4096Hz), data sampling rate is often 8192, sometimes 16384 (8192Hz) 
        filtere=digital_filter(100./f,3000.0/f,50,nfft,/double)     ;100-3000 Hz, adjustable for individial cases, 50 is typically good for suppression of Gibbs spectra leakage/noise 
        Ew=efw_z.y[*,2]    ; this is one dimension Ez data, more code needed to be added here to ensure we have B-field aligned E data (parallel electric field)
        Ew=convol(Ew,filtere,/center)   ;convolution for filtering 
        diff1=min(abs(efw_z.x-tcur[0]),index1)   ;automatically select start index
        diff2=min(abs(efw_z.x-tcur[1]),index2)    ;automatically select end index
        maxEz=max(abs(Ew[index1:index2]))  ; output this TDS maximum amplitude parameter
        
        ;; new code added here to get survey-mode power spectra density derived root-mean-square intensity for frequencies above 50 Hz
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        del_data,'th'+sc+'_'+'f*'
        tr3=[ewtrange[0,i]-30.0, ewtrange[0,i]+30.0]  ;extending the burst time interval to be within 2 min
        thm_load_fgm, probe = sc, trange = tr3, level=2, coord='dsl', data='fgs'
        thm_load_fgm, probe = sc, trange = tr3, level=2, coord='gsm', data='fgs'  ;get the background magnetic field data in gsm coordinates
        thm_load_fbk, probe=sc, trange=tr3, level=2,datatype='f*'
        tplot_names,'th'+sc+'_'+'fb_e?c12',names=fb_null
        if STRCMP(fb_null, '', 1, /FOLD_CASE) then begin
          copy_data,strjoin('th'+sc+'_'+'fb_edc34'),'th'+sc+'_fb_edc'
        endif else begin
          copy_data,strjoin('th'+sc+'_'+'fb_edc12'),'th'+sc+'_fb_edc'
        endelse
        get_data,'th'+sc+'_fb_edc',data=fb_edc   ;get FBK electric field spectra data
        if size(fb_edc, /TYPE) ne 8 then goto, check
        store_data,'Ew_fbk',data= {x:fb_edc.x,y:sqrt(3.0*total(fb_edc.y[*,1:3]^2*(!pi)^2/8.0,2,/nan))} ; ~ 30Hz--1 kHz integration
        get_data,'Ew_fbk',data=Ew_fbk
        rms_Ew=max(Ew_fbk.y,/nan)  ; maximum observed root-mean-square intensity within 1 min of the burst interval
        get_data,strjoin('th'+sc+'_'+'fgs_gsm'),data=fgs_gsm
        mB_gsm=mean(fgs_gsm.y,dimension=1,/nan)   ;get mean value of local magnetic field
        mBx_gsm=mB_gsm[0]
        mBy_gsm=mB_gsm[1]
        mBz_gsm=mB_gsm[2]      
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        tinterpol_mxn,'th'+sc+'_state_pos_gsm','th'+sc+'_efw_clean_fac',newname='th'+sc+'_state_pos_gsm_int' ; interpolate spacecraft position data in gsm coordinates
        time_clip,strjoin('th'+sc+'_state_pos_gsm_int'),ewtrange[0,i],ewtrange[1,i],newname=strjoin('th'+sc+'_state_pos_gsm_clip') ; clip data to specific wave burst intervals
        thm_cotrans,strjoin('th'+sc+'_state_pos_gsm_clip'),strjoin('th'+sc+'_state_pos_sm'),probe=sctlimi,in_coord='gsm',out_coord='sm' ; coordinate transform data into SM coordinates
        tvectot,strjoin('th'+sc+'_state_pos_sm'),tot='pos_tot'    ;calculate the radial distance of the spacecraft in km
        get_data,'pos_tot',data=pos_tot    ; obtain the radial distance
        R=mean(pos_tot.y/6371.0)   ; convert to Earth's radii
        MLT_L_Lat,strjoin('th'+sc+'_state_pos_sm'),MLT0=mlt,L0=L,lat0=lat     ;calculate magnetic local time (mlt) of the spacecraft
        mmlt=mean([mlt[0],mlt[-1]])    ;get mean value of mlt
        
        ;;new code added here to get the mean value of magnetic latitude, position
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        mlat=mean([lat[0],lat[-1]])*180.0/!pi         ; in deg
        get_data,strjoin('th'+sc+'_state_pos_gsm_clip'),data=pos_gsm
        mpos_gsm=mean(pos_gsm.y/6371.0,dimension=1,/nan)   ;get mean value of position in gsm in Re unit
        mposx_gsm=mpos_gsm[0]
        mposy_gsm=mpos_gsm[1]
        mposz_gsm=mpos_gsm[2] 
        ; now also get the ion flow velocity corresponding to TDS observations
        datatype1='peir'
        ;;load particle data
        thm_part_load,probe=sc,trange=tr,datatype=datatype1
        thm_part_products,probe=sc,datatype=datatype1,trange=tr3,outputs='moments',mag_name=strjoin('th'+sc+'_fgs'),pos_name='th'+sc+'_state_pos'
        tplot_names,strjoin('th'+sc+'_peir_velocity'),names=vel_null
        if STRCMP(vel_null, '', 1, /FOLD_CASE) then begin
          mvx_gsm=0.0
          mvy_gsm=0.0
          mvz_gsm=0.0
        endif else begin
          thm_cotrans,strjoin('th'+sc+'_'+datatype1+'_velocity'),in_coord='dsl',out_coord='gsm',out_suffix='_gsm'
          get_data,strjoin('th'+sc+'_'+datatype1+'_velocity_gsm'),data=v_gsm
          mv_gsm_x=v_gsm.y[*,0]
          mv_gsm_y=v_gsm.y[*,1]
          mv_gsm_z=v_gsm.y[*,2]
          ind=where(abs(mv_gsm_x) eq max(abs(mv_gsm_x),/nan)) ; get index of array where Vx is maximum
          mvx_gsm=mv_gsm_x[ind[0]]
          mvy_gsm=mv_gsm_y[ind[0]]
          mvz_gsm=mv_gsm_z[ind[0]]
        endelse    
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
        ;;;;;;;;;;;;;;;;;; the following writing code has been changed
        printf, 1, FORMAT='(A20,A25,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',time_string(ewtrange[0,i]),time_string(ewtrange[1,i]),ewtrange[1,i]-ewtrange[0,i],R,mposx_gsm,mposy_gsm,mposz_gsm,mmlt,mlat,maxEz,rms_Ew,mBx_gsm,mBy_gsm,mBz_gsm,mvx_gsm,mvy_gsm,mvz_gsm   ;save these parameters to file
        ;;;;;;;;;;;;;;;;;;
        end 
      strcmp(B,'-1'): goto,check     ;if you want to adjust time window of plots or some mistakes are made, you can start over (you could erase the data line your just saved and do it again)
      else:
    endcase

  endfor
  stop
  FREE_LUN, 1  ;close file

;start: ewtrange[0,nbw[i]], end: ewtrange[1,nbw[i]] when loop over nbw (the number of burst intervals)

end