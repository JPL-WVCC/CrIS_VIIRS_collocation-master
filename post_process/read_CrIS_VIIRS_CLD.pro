pro read_cris_viirs,year,st_day,ed_day

 cpu, tpool_nthreads=1

mod_taubin=[0.,0.3,1.3,3.6,9.4,23,60,200]
mod_ctpbin=([50,180,310,440,560,680,800,1100])  
nbin_tau=n_elements(mod_taubin)-1  
nbin_ctp=n_elements(mod_ctpbin)-1
outdir='/raid15/qyue/VIIRS/VIIRS/201601/VIIRS_Cloud_onCrIS/'
indexdir='/raid15/qyue/VIIRS/VIIRS/201601/Index/'

CrISL1dir0='/peate_archive/.data1/Ops/snpp/gdisc/2/'
;CrISL2dir='/peate_archive/.data7/Ops/sndr/production/v01_23_01/2015/06/01/l2climcaps/'
;variable_list=['obs_time_tai93','lon','lat','air_temp','air_temp_qc','spec_hum','spec_hum_qc','air_pres_h2o','air_pres',$
;               'land_frac','sol_zen','air_pres_nsurf','air_pres_h2o_nsurf','air_pres_lay','air_temp_dof','h2o_vap_dof',$
;               'rel_hum','rel_hum_qc','cld_frac','cld_frac_qc','cld_top_pres','cld_top_pres_qc','mw_surf_class']
;content_list=['for_cld_frac_tot','for_cld_frac_tot_qc','for_cld_top_pres_tot','for_cld_top_pres_tot_qc']
;clist=['fg_air_temp','fg_h2o_vap_mol_lay','fg_surf_air_temp']
;mollist=['h2o_vap_mol_lay','h2o_vap_mol_lay_qc']

cld_list1=['Cloud_Top_Temperature','Cloud_Top_Pressure','Cloud_Optical_Thickness','Cloud_Effective_Radius','Cloud_Water_Path']
cld_list2=['Cloud_Mask','Cloud_Phase_Cloud_Top_Properties','Cloud_Phase_Optical_Properties']

VIIRSL2dir='/raid15/qyue/VIIRS/VIIRS/201601/CLDPROP/CLDPROP_L2_VIIRS_SNPP/2016/'  ;001/'

for iday=st_day,ed_day do begin

CrISL1dir=CrISL1dir0+'2016/01/'+string(iday,format='(i2.2)')+'/crisl1b/'

mode=string(year,format='(i4.4)')+'01'+strtrim(string(iday),2)
mode1=string(year,format='(i4.4)')+string(iday,format='(i3.3)')

;VIIRSL2dir='/raid15/qyue/VIIRS/VIIRS/201501/CLDPROP/CLDPROP_L2_VIIRS_SNPP/2015/'+string(iday,format='(i3.3)')+'/'
spawn,'ls '+VIIRSL2dir+string(iday,format='(i3.3)')+'/CLDPROP_L2_VIIRS_SNPP.A'+mode1+'*.nc',VIIRScldfiles


for ifile=0,239 do begin
    indexfile=indexdir+'IND_CrIS_VIIRSMOD_'+mode+'_'+strcompress(string(ifile),/remove_all)+'.nc'
    ierr=read_ncdf(indexfile,index)
    ngran=3
   print,'index_file',indexfile
    result=size(index.dy_size,/dimension)
    dysize=index.dy_size
    ysize2=total(reform(index.dy_size,result[0]*result[1]*result[2]),/cumulative,/integer)   ;cumulatively the end index of dy in each CrIS pixel

   okfile=-1
   igran=ifile+1
   spawn,'ls '+CrISL1dir+'SNDR*.g'+string(igran,format='(i3.3)')+'.L1B*.nc',crisL1file
   if crisL1file[0] eq '' then continue
   ierr=read_ncdf(crisL1file[0],data1,variable_list=['lat','lon','obs_time_tai93']) 
   CrIS_lat=data1.lat
   CrIS_lon=data1.lon
   CrIS_tai93=data1.obs_time_tai93
   
   if ifile eq 0 then begin
       ngran=2 & viirsfiles=viirscldfiles[ifile:ifile+1]
   endif 
   if ifile eq 239 then begin
       ngran=2 & viirsfiles=viirscldfiles[ifile-1:ifile]
   endif
   if ifile gt 0 and ifile lt 239 then begin
       ngran=3 & viirsfiles=viirscldfiles[ifile-1:ifile+1]
   endif
   
   okfile=-1
   for iVgran=0,ngran-1 do begin 

          okfile=okfile+1 
          viirsfile=viirsfiles[ivgran] ;#of ok VIIRS files matched to this CrIS granule        
          ierr=read_nc4(viirsfile,data,groupname='geophysical_data',$
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

   VIIRS_ISCCP_Hist=fltarr(nbin_tau,nbin_ctp,result[0],result[1],result[2])*!values.f_nan 
   VIIRS_CLD_Mask=fltarr(5,result[0],result[1],result[2])*!values.f_nan
   VIIRS_CLDTOP_Phase=fltarr(5,result[0],result[1],result[2])*!values.f_nan
   VIIRS_OPT_Phase=fltarr(5,result[0],result[1],result[2])*!values.f_nan
   VIIRS_CTP_moment=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_COD_moment=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_CTT_moment=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_CWP_moment=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_Re_moment=fltarr(4,result[0],result[1],result[2])*!values.f_nan   
   VIIRS_Re_ISCCPHist=fltarr(nbin_tau,nbin_ctp,result[0],result[1],result[2])*!values.f_nan
   VIIRS_SFC_Mask=fltarr(10,result[0],result[1],result[2])*!values.f_nan
   
   dysize=index.dy_size
   findcol=where(dysize GE 1,nfindcol)
   col = findcol mod result[0]
   row = (findcol / result[0]) mod result[1]
   frame = findcol / (result[1]*result[0])
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
           VIIRS_CLD_Mask[im,col[i],row[i],frame[i]]=nfind
           find=where(cphaseopt[dx,dy] eq im,nfind)
           VIIRS_OPT_Phase[im,col[i],row[i],frame[i]]=nfind
           find=where(cphasetop[dx,dy] eq im,nfind)
           if im eq 4 then find=where(cphasetop[dx,dy] eq 6,nfind)
           VIIRS_CLDTOP_Phase[im,col[i],row[i],frame[i]]=nfind
       endfor
;set up sfc mask:
       for im=0,1 do begin
           find=where(daymask[dx,dy] eq im,nfind)
           VIIRS_SFC_Mask[im,col[i],row[i],frame[i]]=nfind
       endfor
       for im=0,1 do begin
           find=where(glintmask[dx,dy] eq im,nfind)
           VIIRS_SFC_Mask[im+2,col[i],row[i],frame[i]]=nfind
       endfor
       for im=0,1 do begin
           find=where(icemask[dx,dy] eq im,nfind)
           VIIRS_SFC_Mask[im+4,col[i],row[i],frame[i]]=nfind
       endfor
       for im=0,3 do begin
           find=where(sfcmask[dx,dy] eq im,nfind)
           VIIRS_SFC_Mask[im+6,col[i],row[i],frame[i]]=nfind
       endfor       
       for itau=0,nbin_tau-1 do begin
           for ictp=0,nbin_ctp-1 do begin
               find=where(Ctp[dx,dy] Gt mod_ctpbin[ictp] and ctp[dx,dy] LE mod_ctpbin[ictp+1] and COD[dx,dy] Gt mod_taubin[itau] and COD[dx,dy] LE mod_taubin[itau+1],nfind)  
               VIIRS_ISCCP_Hist[itau,ictp,col[i],row[i],frame[i]]=nfind
               find=where(Ctp[dx,dy] Gt mod_ctpbin[ictp] and ctp[dx,dy] LE mod_ctpbin[ictp+1] and COD[dx,dy] Gt mod_taubin[itau] and COD[dx,dy] LE mod_taubin[itau+1] and finite(re1) ,nfind)
               if nfind GE 1 then VIIRS_Re_ISCCPHist[itau,ictp,col[i],row[i],frame[i]]=mean(re1[find],/nan)
           endfor
       endfor
       VIIRS_CTP_moment[*,col[i],row[i],frame[i]]=moment(ctp[dx,dy],/nan)
       VIIRS_COD_moment[*,col[i],row[i],frame[i]]=moment(COD[dx,dy],/nan)
       VIIRS_CTT_moment[*,col[i],row[i],frame[i]]=moment(CTT[dx,dy],/nan)
       VIIRS_Re_moment[*,col[i],row[i],frame[i]]=moment(Re[dx,dy],/nan)
       VIIRS_CWP_moment[*,col[i],row[i],frame[i]]=moment(CWP[dx,dy],/nan)


   endfor
   save,CrIS_lon,CrIS_lat,CrIS_tai93,VIIRS_SFC_Mask,VIIRS_ISCCP_Hist,VIIRS_Re_ISCCPHist,$
        VIIRS_CLD_Mask,VIIRS_CTP_moment,VIIRS_COD_moment,VIIRS_CLDTOP_Phase,VIIRS_OPT_Phase,VIIRS_CTT_moment,VIIRS_Re_moment,VIIRS_CWP_moment,$
        /compress,filename=outdir+'VIIRS_Cloud_onCrIS_'+mode+'_G'+string(igran,format='(i3.3)')+'.sav'
   
endfor

endfor

stop

end


       
           
