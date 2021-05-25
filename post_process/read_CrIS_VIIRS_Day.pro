pro read_cirs_viirs,year,day

mod_taubin=[0.,0.3,1.3,3.6,9.4,23,60,200]
mod_ctpbin=([50,180,310,440,560,680,800,1100])  
nbin_tau=n_elements(mod_taubin)-1  
nbin_ctp=n_elements(mod_ctpbin)-1

indexdir='/home/qyue/VIIRS/CrIS_VIIRS_collocation-master/'

CrISL2dir='/peate_archive/.data7/Ops/sndr/production/v01_23_01/2015/06/01/l2climcaps/'
;variable_list=['obs_time_tai93','lon','lat','air_temp','air_temp_qc','spec_hum','spec_hum_qc','air_pres_h2o','air_pres',$
;               'land_frac','sol_zen','air_pres_nsurf','air_pres_h2o_nsurf','air_pres_lay','air_temp_dof','h2o_vap_dof',$
;               'rel_hum','rel_hum_qc','cld_frac','cld_frac_qc','cld_top_pres','cld_top_pres_qc','mw_surf_class']
;content_list=['for_cld_frac_tot','for_cld_frac_tot_qc','for_cld_top_pres_tot','for_cld_top_pres_tot_qc']
;clist=['fg_air_temp','fg_h2o_vap_mol_lay','fg_surf_air_temp']
;mollist=['h2o_vap_mol_lay','h2o_vap_mol_lay_qc']

VIIRSL2dir='/raid15/qyue/VIIRS/VIIRS/20150601/'

spawn,'ls '+VIIRSL2dir+'CLDPROP_L2_VIIRS_SNPP.*.nc',VIIRScldfiles


for ifile=0,23 do begin
    indexfile=indexdir+'sample'+strcompress(string(10*ifile),/remove_all)+'.nc'
    ierr=read_ncdf(indexfile,index)
    ngran=10
   print,'index_file',indexfile
    result=size(index.dy_size,/dimension)
stop

    ysize2=reform(total(reform(index.dy_size,result[0]*result[1]*result[2]),/cumulative),result[0],result[1],result[2])   ;cumulatively the end index of dy in each CrIS pixel
 
   okfile=-1   
   for iVgran=0,ngran-1 do begin 

          okfile=okfile+1 
          viirsfile=viirscldfiles[ifile*ngran+ivgran] ;#of ok VIIRS files matched to this CrIS granule        
          ierr=read_nc4(viirsfile,data,groupname='geophysical_data',$
              content_list=['Cloud_Top_Pressure','Cloud_Optical_Thickness','Cloud_Mask','Cloud_Effective_Radius'])
          ierr=read_nc4(viirsfile,data1,groupname='geolocation_data',$
              content_list=['latitude','longitude'])
          if okfile eq 0 then begin
             lat=data1.latitude
             lon=data1.longitude
             ctp=data.cloud_top_pressure*0.1 
             cod=data.cloud_optical_thickness*0.01
             Re=data.cloud_effective_radius*0.01
             cmaskt=data.cloud_mask
             data=0 & data1=0       
          endif
          if okfile GE 1 then begin  
             lat=concat(lat,data1.latitude,dimension=1)
             lon=concat(lon,data1.longitude,dimension=1)    
             ctp=concat(ctp,data.cloud_top_pressure*0.1,dimension=1)
             cod=concat(cod,data.cloud_optical_thickness*0.01,dimension=1)
             Re=concat(Re,data.cloud_effective_radius*0.01,dimension=1)
             cmaskt=concat(cmaskt,data.cloud_mask,dimension=2)
             data=0 & data1=0
          endif
   endfor
   if okfile eq -1 then continue
   cmask=reform(cmaskt[0,*,*])*(-1)
   find=where(reform(cmaskt[0,*,*].bitget(2)) eq 0 and reform(cmaskt[0,*,*].bitget(1)) eq 0,nfind)
   if nfind GE 1 then cmask[find]=0 ;cloudy
   find=where(reform(cmaskt[0,*,*].bitget(2)) eq 0 and reform(cmaskt[0,*,*].bitget(1)) eq 1,nfind)
   if nfind GE 1 then cmask[find]=1 ;uncertain
   find=where(reform(cmaskt[0,*,*].bitget(2)) eq 1 and reform(cmaskt[0,*,*].bitget(1)) eq 0,nfind)
   if nfind GE 1 then cmask[find]=2 ;possibly clear
   find=where(reform(cmaskt[0,*,*].bitget(2)) eq 1 and reform(cmaskt[0,*,*].bitget(1)) eq 1,nfind)
   if nfind GE 1 then cmask[find]=3 ;Confidently clear
   cmaskt=0
   find=where(ctp lt 1,nfind)
   if nfind GE 1 then ctp[find]=!values.f_nan
   find=where(cod lt 0,nfind)
   if nfind GE 1 then COD[find]=!values.f_nan  
   find=where(Re lt 0,nfind)
   if nfind GE 1 then Re[find]=!values.f_nan 

   VIIRS_ISCCP_Hist=fltarr(nbin_tau,nbin_ctp,result[0],result[1],result[2])*!values.f_nan 
   VIIRS_COD_Hist=fltarr(nbin_tau,result[0],result[1],result[2])*!values.f_nan
   VIIRS_CTP_Hist=fltarr(nbin_CTP,result[0],result[1],result[2])*!values.f_nan
   VIIRS_CLD_Mask=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_CTP_moment=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_COD_moment=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_Re_ISCCPHist=fltarr(nbin_tau,nbin_ctp,result[0],result[1],result[2])*!values.f_nan
   VIIRS_lat=fltarr(result[0],result[1],result[2])*!values.f_nan
   VIIRS_lon=fltarr(result[0],result[1],result[2])*!values.f_nan
   
   dysize=index.dy_size
   for k=0,result[0]-1 do begin
       for j=0,result[1]-1 do begin
           for i=0,result[2]-1 do begin
               if dysize[k,j,i] lT 1 then continue
               offset1=k+9*(j+30*i) 
               if offset1 GE 1 then begin
                     dy=index.dy[ysize2[offset1-1]:ysize2[offset1]-1]
                     dx=index.dx[ysize2[offset1-1]:ysize2[offset1]-1]
               endif else begin
                     dy=index.dy[0:ysize2[offset1]-1]
                     dx=index.dx[0:ysize2[offset1]-1]
               endelse
               re1=Re[dx,dy]
               for im=0,3 do begin
                   find=where(cmask[dx,dy] eq im,nfind)
                   VIIRS_CLD_Mask[im,k,j,i]=nfind
               endfor
               for itau=0,nbin_tau-1 do begin
                   find=where(COD[dx,dy] Gt mod_taubin[itau] and COD[dx,dy] LE mod_taubin[itau+1],nfind)  
                   VIIRS_COD_Hist[itau,k,j,i]=nfind
                   for ictp=0,nbin_ctp-1 do begin
                       find=where(Ctp[dx,dy] Gt mod_ctpbin[ictp] and ctp[dx,dy] LE mod_ctpbin[ictp+1] and COD[dx,dy] Gt mod_taubin[itau] and COD[dx,dy] LE mod_taubin[itau+1],nfind)  
                       VIIRS_ISCCP_Hist[itau,ictp,k,j,i]=nfind
                       find=where(Ctp[dx,dy] Gt mod_ctpbin[ictp] and ctp[dx,dy] LE mod_ctpbin[ictp+1] and COD[dx,dy] Gt mod_taubin[itau] and COD[dx,dy] LE mod_taubin[itau+1] and finite(re1) ,nfind)
                       if nfind GE 1 then VIIRS_Re_ISCCPHist[itau,ictp,k,j,i]=mean(re1[find],/nan)
                   endfor
               endfor
               for ictp=0,nbin_ctp-1 do begin
                   find=where(Ctp[dx,dy] GE mod_ctpbin[ictp] and ctp[dx,dy] LE mod_ctpbin[ictp+1],nfind)  
                   VIIRS_Ctp_Hist[ictp,k,j,i]=nfind
               endfor
               VIIRS_CTP_moment[*,k,j,i]=moment(ctp[dx,dy],/nan)
               VIIRS_COD_moment[*,k,j,i]=moment(COD[dx,dy],/nan)
               VIIRS_lat[k,j,i]=mean(lat[dx,dy],/nan)
               VIIRS_lon[k,j,i]=mean(lon[dx,dy],/nan)
           endfor
       endfor
   endfor
   save,VIIRS_lon,VIIRS_lat,VIIRS_ISCCP_Hist,VIIRS_Re_ISCCPHist,VIIRS_COD_Hist,VIIRS_CTP_Hist,VIIRS_CLD_Mask,VIIRS_CTP_moment,VIIRS_COD_moment,$
        filename='VIIRS_Cloud_onCrIS_20150601_F'+string(ifile,format='(i3.3)')+'.sav'
   
endfor

stop

end


       
           
