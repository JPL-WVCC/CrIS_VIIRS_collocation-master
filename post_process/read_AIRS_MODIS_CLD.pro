pro read_AIRS_MODIS,year,st_day,ed_day

 cpu, tpool_nthreads=1

mod_taubin=[0.,0.3,1.3,3.6,9.4,23,60,200]
mod_ctpbin=([50,180,310,440,560,680,800,1100])  
nbin_tau=n_elements(mod_taubin)-1  
nbin_ctp=n_elements(mod_ctpbin)-1
outdir='/raid15/qyue/VIIRS/MODIS/201601/MODIS_Cloud_onAIRS/'
indexdir='/raid15/qyue/VIIRS/MODIS/201601/Index/'

AIRSL1dir0='/archive/AIRSOps/airs/gdaac/v5/'
;CrISL2dir='/peate_archive/.data7/Ops/sndr/production/v01_23_01/2015/06/01/l2climcaps/'
;variable_list=['obs_time_tai93','lon','lat','air_temp','air_temp_qc','spec_hum','spec_hum_qc','air_pres_h2o','air_pres',$
;               'land_frac','sol_zen','air_pres_nsurf','air_pres_h2o_nsurf','air_pres_lay','air_temp_dof','h2o_vap_dof',$
;               'rel_hum','rel_hum_qc','cld_frac','cld_frac_qc','cld_top_pres','cld_top_pres_qc','mw_surf_class']
;content_list=['for_cld_frac_tot','for_cld_frac_tot_qc','for_cld_top_pres_tot','for_cld_top_pres_tot_qc']
;clist=['fg_air_temp','fg_h2o_vap_mol_lay','fg_surf_air_temp']
;mollist=['h2o_vap_mol_lay','h2o_vap_mol_lay_qc']

cld_list1=['Cloud_Top_Temperature','Cloud_Top_Pressure','Cloud_Optical_Thickness','Cloud_Effective_Radius','Cloud_Water_Path']
cld_list2=['Cloud_Mask','Cloud_Phase_Cloud_Top_Properties','Cloud_Phase_Optical_Properties']

MODISL2dir='/raid15/qyue/VIIRS/MODIS/201601/CLDPROP/CLDPROP_L2_MODIS_Aqua/2016/'

for iday=st_day,ed_day do begin

AIRSL1dir=AIRSL1dir0+'2016/01/'+string(iday,format='(i2.2)')+'/airibrad/'

mode=string(year,format='(i4.4)')+'01'+strtrim(string(iday),2)
mode1=string(year,format='(i4.4)')+string(iday,format='(i3.3)')

;VIIRSL2dir='/raid15/qyue/VIIRS/VIIRS/201501/CLDPROP/CLDPROP_L2_VIIRS_SNPP/2015/'+string(iday,format='(i3.3)')+'/'
spawn,'ls '+MODISL2dir+string(iday,format='(i3.3)')+'/CLDPROP_L2_MODIS_Aqua.A'+mode1+'*.nc',MODIScldfiles

for ifile=0,23 do begin
    indexfile=indexdir+'IND_AIRS_MODISMOD_'+mode+'_'+strcompress(string(ifile),/remove_all)+'.nc'
    ierr=read_ncdf(indexfile,index)
    ngran=10
   print,'index_file',indexfile
    result=size(index.dy_size,/dimension)
    dysize=index.dy_size
    ysize2=total(reform(index.dy_size,result[0]*result[1]),/cumulative,/integer)   ;cumulatively the end index of dy in each CrIS pixel

   okfile=-1
   for igran=0,ngran-1 do begin
      spawn,'ls '+AIRSL1dir+'AIRS*.'+string(ifile*ngran+igran+1,format='(i3.3)')+'.L1B*.hdf',AIRSL1file 
      okfile=okfile+1
      ierr=read_airs_swath(AIRSL1file[0],3,data1,content_list=['Latitude','Longitude','Time']) 
      if okfile eq 0 then begin
         AIRS_lat1=data1.latitude
         AIRS_lon1=data1.longitude
         AIRS_tai931=data1.time
      endif
      if okfile GE 1 then begin
         AIRS_lat1=concat(AIRS_lat1,data1.latitude,dimension=1)
         AIRS_lon1=concat(AIRS_lon1,data1.longitude,dimension=1)
         AIRS_tai931=concat(AIRS_tai931,data1.time,dimension=1)
      endif
   endfor
       
   okfile=-1
   if ifile eq 0 then begin
   MODISfiles=MODIScldfiles[ifile*12:(ifile+1)*12]
   ngran=13
   endif
   if ifile eq 23 then begin
      MODISfiles=MODIScldfiles[ifile*12-1:(ifile+1)*12-1]
      ngran=13
   endif
   if ifile gt 0 and ifile lt 23 then begin
      MODISfiles=MODIScldfiles[ifile*12-1:(ifile+1)*12]
      ngran=14
   endif

   for iVgran=0,ngran-1 do begin 

          okfile=okfile+1 
          modisfile=modisfiles[ivgran] ;#of ok VIIRS files matched to this CrIS granule        
          ierr=read_nc4(modisfile,data,groupname='geophysical_data',$
              content_list=[cld_list1,cld_list2])
          if okfile eq 0 then begin
             ctp=data.cloud_top_pressure*0.1 
             cod=data.cloud_optical_thickness*0.01
             Re=data.cloud_effective_radius*0.01
             ctt=data.cloud_top_temperature*0.008+100.
             CWP=data.cloud_water_path
             cphasetop=data.Cloud_Phase_Cloud_Top_Properties
             cphaseopt=data.Cloud_Phase_Optical_Properties
             cmaskt=data.cloud_mask
             data=0        
          endif
          if okfile GE 1 then begin  
             ctp=concat(ctp,data.cloud_top_pressure*0.1,dimension=1)
             cod=concat(cod,data.cloud_optical_thickness*0.01,dimension=1)
             Re=concat(Re,data.cloud_effective_radius*0.01,dimension=1)
             ctt=concat(ctt,data.cloud_top_temperature*0.008+100.,dimension=1)
             CWP=concat(cwp,data.cloud_water_path,dimension=1)
             cphasetop=concat(cphasetop,data.Cloud_Phase_Cloud_Top_Properties,dimension=1)
             cphaseopt=concat(cphaseopt,data.Cloud_Phase_Optical_Properties,dimension=1)
             cmaskt=concat(cmaskt,data.cloud_mask,dimension=2)
             data=0 
          endif
   endfor
    
   if okfile eq -1 then continue
   cmask=reform(cmaskt[0,*,*])*(-1)
   
   find=where(reform(cmaskt[0,*,*].bitget(0)) eq 0,nfind)
   if nfind GE 1 then cmask[find]=-1
   find=where(reform(cmaskt[0,*,*].bitget(2)) eq 0 and reform(cmaskt[0,*,*].bitget(1)) eq 0,nfind)
   if nfind GE 1 then cmask[find]=0 ;cloudy
   find=where(reform(cmaskt[0,*,*].bitget(2)) eq 0 and reform(cmaskt[0,*,*].bitget(1)) eq 1,nfind)
   if nfind GE 1 then cmask[find]=1 ;uncertain
   find=where(reform(cmaskt[0,*,*].bitget(2)) eq 1 and reform(cmaskt[0,*,*].bitget(1)) eq 0,nfind)
   if nfind GE 1 then cmask[find]=2 ;possibly clear
   find=where(reform(cmaskt[0,*,*].bitget(2)) eq 1 and reform(cmaskt[0,*,*].bitget(1)) eq 1,nfind)
   if nfind GE 1 then cmask[find]=3 ;Confidently clear
   
   daymask=reform(cmaskt[0,*,*])*(-1) & glintmask=daymask & icemask=daymask & sfcmask=daymask
   find=where(reform(cmaskt[0,*,*].bitget(3)) eq 0,nfind)
   if nfind GE 1 then daymask[find]=0 ;night
    find=where(reform(cmaskt[0,*,*].bitget(3)) eq 1,nfind)
   if nfind GE 1 then daymask[find]=1 ;day 
   find=where(reform(cmaskt[0,*,*].bitget(4)) eq 0,nfind)
   if nfind GE 1 then glintmask[find]=0 ;sunglint
    find=where(reform(cmaskt[0,*,*].bitget(4)) eq 1,nfind)
   if nfind GE 1 then glintmask[find]=1 ;not sunglint 
    find=where(reform(cmaskt[0,*,*].bitget(5)) eq 0,nfind)
   if nfind GE 1 then icemask[find]=0 ;ice/snow
    find=where(reform(cmaskt[0,*,*].bitget(5)) eq 1,nfind)
   if nfind GE 1 then icemask[find]=1 ;unfrozen
   find=where(reform(cmaskt[0,*,*].bitget(7)) eq 0 and reform(cmaskt[0,*,*].bitget(6)) eq 0,nfind)
   if nfind GE 1 then sfcmask[find]=0 ;water
   find=where(reform(cmaskt[0,*,*].bitget(7)) eq 0 and reform(cmaskt[0,*,*].bitget(6)) eq 1,nfind)
   if nfind GE 1 then sfcmask[find]=1 ;coast
   find=where(reform(cmaskt[0,*,*].bitget(7)) eq 1 and reform(cmaskt[0,*,*].bitget(6)) eq 0,nfind)
   if nfind GE 1 then sfcmask[find]=2 ;desert
   find=where(reform(cmaskt[0,*,*].bitget(7)) eq 1 and reform(cmaskt[0,*,*].bitget(6)) eq 1,nfind)
   if nfind GE 1 then sfcmask[find]=3 ;land 

   cmaskt=0
   
   find=where(ctp lt 1,nfind)
   if nfind GE 1 then ctp[find]=!values.f_nan
   find=where(cod lt 0,nfind)
   if nfind GE 1 then COD[find]=!values.f_nan  
   find=where(Re lt 0,nfind)
   if nfind GE 1 then Re[find]=!values.f_nan 
   find=where(ctt lt 0,nfind)
   if nfind GE 1 then ctt[find]=!values.f_nan  
   find=where(cwp lt 0,nfind)
   if nfind GE 1 then cwp[find]=!values.f_nan

   MODIS_ISCCP_Hist1=fltarr(nbin_tau,nbin_ctp,result[0],result[1])*!values.f_nan 
   MODIS_CLD_Mask1=fltarr(5,result[0],result[1])*!values.f_nan
   MODIS_CLDTOP_Phase1=fltarr(5,result[0],result[1])*!values.f_nan
   MODIS_OPT_Phase1=fltarr(5,result[0],result[1])*!values.f_nan
   MODIS_CTP_moment1=fltarr(4,result[0],result[1])*!values.f_nan
   MODIS_COD_moment1=fltarr(4,result[0],result[1])*!values.f_nan
   MODIS_CTT_moment1=fltarr(4,result[0],result[1])*!values.f_nan
   MODIS_CWP_moment1=fltarr(4,result[0],result[1])*!values.f_nan
   MODIS_Re_moment1=fltarr(4,result[0],result[1])*!values.f_nan   
   MODIS_Re_ISCCPHist1=fltarr(nbin_tau,nbin_ctp,result[0],result[1])*!values.f_nan
   MODIS_SFC_Mask1=fltarr(10,result[0],result[1])*!values.f_nan
   
   dysize=index.dy_size
   findcol=where(dysize GE 1,nfindcol)
   col = findcol mod result[0]
   row = (findcol / result[0]) mod result[1]
 ;  frame = findcol / (result[1]*result[0])
   for i=0,nfindcol-1 do begin
       time=systime(/seconds)
       if i eq 0 and findcol[0] eq 0 then begin
            dy=index.dy[0:ysize2[0]-1] & dx=index.dx[0:ysize2[0]-1]
       endif else begin
            dy=index.dy[ysize2[findcol[i-1]]:ysize2[findcol[i]]-1]
            dx=index.dx[ysize2[findcol[i-1]]:ysize2[findcol[i]]-1]
       endelse      
       re1=Re[dx,dy]
       for im=0,4 do begin
           find=where(cmask[dx,dy] eq im-1,nfind)
           MODIS_CLD_Mask1[im,col[i],row[i]]=nfind
           find=where(cphaseopt[dx,dy] eq im,nfind)
           MODIS_OPT_Phase1[im,col[i],row[i]]=nfind
           find=where(cphasetop[dx,dy] eq im,nfind)
           if im eq 4 then find=where(cphasetop[dx,dy] eq 6,nfind)
           MODIS_CLDTOP_Phase1[im,col[i],row[i]]=nfind
       endfor
;set up sfc mask:
       for im=0,1 do begin
           find=where(daymask[dx,dy] eq im,nfind)
           MODIS_SFC_Mask1[im,col[i],row[i]]=nfind
       endfor
       for im=0,1 do begin
           find=where(glintmask[dx,dy] eq im,nfind)
           MODIS_SFC_Mask1[im+2,col[i],row[i]]=nfind
       endfor
       for im=0,1 do begin
           find=where(icemask[dx,dy] eq im,nfind)
           MODIS_SFC_Mask1[im+4,col[i],row[i]]=nfind
       endfor
       for im=0,3 do begin
           find=where(sfcmask[dx,dy] eq im,nfind)
           MODIS_SFC_Mask1[im+6,col[i],row[i]]=nfind
       endfor
       for itau=0,nbin_tau-1 do begin
           for ictp=0,nbin_ctp-1 do begin
               find=where(Ctp[dx,dy] Gt mod_ctpbin[ictp] and ctp[dx,dy] LE mod_ctpbin[ictp+1] and COD[dx,dy] Gt mod_taubin[itau] and COD[dx,dy] LE mod_taubin[itau+1],nfind)  
               MODIS_ISCCP_Hist1[itau,ictp,col[i],row[i]]=nfind
               find=where(Ctp[dx,dy] Gt mod_ctpbin[ictp] and ctp[dx,dy] LE mod_ctpbin[ictp+1] and COD[dx,dy] Gt mod_taubin[itau] and COD[dx,dy] LE mod_taubin[itau+1] and finite(re1) ,nfind)
               if nfind GE 1 then MODIS_Re_ISCCPHist1[itau,ictp,col[i],row[i]]=mean(re1[find],/nan)
           endfor
       endfor
       MODIS_CTP_moment1[*,col[i],row[i]]=moment(ctp[dx,dy],/nan)
       MODIS_COD_moment1[*,col[i],row[i]]=moment(COD[dx,dy],/nan)
       MODIS_CTT_moment1[*,col[i],row[i]]=moment(CTT[dx,dy],/nan)
       MODIS_Re_moment1[*,col[i],row[i]]=moment(Re[dx,dy],/nan)
       MODIS_CWP_moment1[*,col[i],row[i]]=moment(CWP[dx,dy],/nan)


   endfor
   for igran=0,9 do begin
      AIRS_lat=AIRS_lat1[*,igran*135:(igran+1)*135-1]
      AIRS_lon=AIRS_lon1[*,igran*135:(igran+1)*135-1]
      AIRS_tai93=AIRS_tai931[*,igran*135:(igran+1)*135-1]
      MODIS_ISCCP_Hist=MODIS_ISCCP_Hist1[*,*,*,igran*135:(igran+1)*135-1]
      MODIS_CLD_Mask=MODIS_CLD_Mask1[*,*,igran*135:(igran+1)*135-1]
      MODIS_CLDTOP_Phase=MODIS_CLDTOP_Phase1[*,*,igran*135:(igran+1)*135-1]
      MODIS_OPT_Phase=MODIS_OPT_Phase1[*,*,igran*135:(igran+1)*135-1]
      MODIS_CTP_moment=MODIS_CTP_moment1[*,*,igran*135:(igran+1)*135-1]
      MODIS_COD_moment=MODIS_COD_moment1[*,*,igran*135:(igran+1)*135-1]
      MODIS_CTT_moment=MODIS_CTT_moment1[*,*,igran*135:(igran+1)*135-1]
      MODIS_CWP_moment=MODIS_CWP_moment1[*,*,igran*135:(igran+1)*135-1]
      MODIS_Re_moment=MODIS_Re_moment1[*,*,igran*135:(igran+1)*135-1]
      MODIS_Re_ISCCPHist=MODIS_Re_ISCCPHist1[*,*,*,igran*135:(igran+1)*135-1]
      MODIS_SFC_Mask=MODIS_SFC_MAsk1[*,*,igran*135:(igran+1)*135-1]
      save,AIRS_lon,AIRS_lat,AIRS_tai93,MODIS_SFC_Mask,MODIS_ISCCP_Hist,MODIS_Re_ISCCPHist,$
        MODIS_CLD_Mask,MODIS_CTP_moment,MODIS_COD_moment,MODIS_CLDTOP_Phase,MODIS_OPT_Phase,MODIS_CTT_moment,MODIS_Re_moment,MODIS_CWP_moment,$
        /compress,filename=outdir+'MODIS_Cloud_onAIRS_'+mode+'_G'+string(ifile*10+igran+1,format='(i3.3)')+'.sav'
  endfor
   
endfor

endfor

stop

end


       
           
