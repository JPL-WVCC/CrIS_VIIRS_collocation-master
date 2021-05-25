pro read_cirs_viirs,year,day

mod_taubin=[0.,0.3,1.3,3.6,9.4,23,60,200]
mod_ctpbin=([50,180,310,440,560,680,800,1100])  
nbin_tau=n_elements(mod_taubin)-1  
nbin_ctp=n_elements(mod_ctpbin)-1

indexfile='/home/qyue/VIIRS/CrIS_VIIRS_collocation-master/sample.nc'

CrISL2dir='/peate_archive/.data7/Ops/sndr/production/v01_23_01/2015/06/01/l2climcaps/'
;variable_list=['obs_time_tai93','lon','lat','air_temp','air_temp_qc','spec_hum','spec_hum_qc','air_pres_h2o','air_pres',$
;               'land_frac','sol_zen','air_pres_nsurf','air_pres_h2o_nsurf','air_pres_lay','air_temp_dof','h2o_vap_dof',$
;               'rel_hum','rel_hum_qc','cld_frac','cld_frac_qc','cld_top_pres','cld_top_pres_qc','mw_surf_class']
;content_list=['for_cld_frac_tot','for_cld_frac_tot_qc','for_cld_top_pres_tot','for_cld_top_pres_tot_qc']
;clist=['fg_air_temp','fg_h2o_vap_mol_lay','fg_surf_air_temp']
;mollist=['h2o_vap_mol_lay','h2o_vap_mol_lay_qc']

VIIRSL2dir='/raid15/qyue/VIIRS/VIIRS/20150601/'

ierr=read_ncdf(indexfile,index)

ngran=240
result=size(index.dy_size,/dimension)

ysize2=reform(total(reform(index.dy_size,result[0]*result[1]*result[2]),/cumulative),result[0],result[1],result[2])   ;cumulatively the end index of dy in each CrIS pixel

spawn,'ls '+VIIRSL2dir+'VNP03MOD*.nc',VIIRSgeos

;VIIRS L2 files not complete, this is to find out which file is ready and the dimension of each file
gran=-1 & viirsL2files='' & row=0 & col=0
for ifile=0,n_elements(VIIRSgeos)-1 do begin
      ;VIIRSgeo=ncdf_open(VIIRSgeos[ifile],/nowrite)
      ncdf_list,VIIRSgeos[ifile],/dimensions,out=dim,/quiet
      row=[row,float(strmid(dim[9],strlen(dim[9])-4,4))]
     ; col=[col,float(strmid(dim[10],strlen(dim[10]-4,4))]
      VIIRScldfile=VIIRSL2dir+'VNPCLDPROP.A2015152.'+strmid(VIIRSgeos[ifile],strpos(VIIRSgeos[ifile],'.A2015152.')+10,4)+'*.nc'
      spawn,'ls '+VIIRScldfile,cldfile
      if cldfile[0] eq '' then continue
      gran=[gran,ifile] ;saved VIIRS granules that I have
      viirsL2files=[viirsL2files,cldfile[0]]    
endfor

;record # of rows cumulatively, first element is zero, row dimension is one element more than granule numbers 
row=total(row,/cumulative)
gran=gran[1:*] & viirsL2files=viirsL2files[1:*]
rowbnd1=row[gran-1]-1 & rowbnd2=row[gran]-1  ;this keeps the row bounds of available VIIRS granules


;loop over all CrIS Granules 
;first check whether CrIS L2 file is ready
;then calculate the dy and dysize for each granule
for igran=0,ngran-1 do begin
   spawn,'ls '+CrISL2dir+'SNDR*.g'+string(igran+1,format='(i3.3)')+'.L2_RET*.nc',crisL2file
   if crisL2file[0] eq '' then continue
 
   dysize=index.dy_size[*,*,45l*igran:45l*(igran+1)-1] ;this is the number of VIIRS pixels
   if igran eq 0 then dy=index.dy[0:max(ysize2[*,*,0:45l*(igran+1)-1])-1] else $
   dy=index.dy[max(ysize2[*,*,45l*(igran-1)]):max(ysize2[*,*,45l*igran])-1]  
 
   okfile=-1
   for iVgran=0,n_elements(gran)-1 do begin
       find=where(dy GE rowbnd1[ivgran] and dy LE rowbnd2[ivgran],nfind)
      
       if nfind GE 1 then begin
          okfile=okfile+1 & print,viirsL2files[ivgran] ;#of ok VIIRS files matched to this CrIS granule        
          ierr=read_nc4(viirsL2files[ivgran],data,groupname='geophysical_data',$
              content_list=['Cloud_Top_Pressure','Cloud_Optical_Thickness','Cloud_Mask'])
          ierr=read_nc4(viirsL2files[ivgran],data1,groupname='geolocation_data',$
              content_list=['latitude','longitude'])
          if okfile eq 0 then begin
             bnd1=rowbnd1[ivgran]
             lat=data1.latitude
             lon=data1.longitude
             ctp=data.cloud_top_pressure*0.1 
             cod=data.cloud_optical_thickness*0.01
             cmaskt=data.cloud_mask
             data=0 & data1=0       
          endif
          if okfile GE 1 then begin  
             bnd2=rowbnd2[ivgran] 
             lat=concat(lat,data1.latitude,dimension=1)
             lon=concat(lon,data1.longitude,dimension=1)    
             ctp=concat(ctp,data.cloud_top_pressure*0.1,dimension=1)
             cod=concat(cod,data.cloud_optical_thickness*0.01,dimension=1)
             cmaskt=concat(cmaskt,data.cloud_mask,dimension=2)
             data=0 & data1=0
          endif
       endif
   endfor
   if okfile eq -1 then continue
   print,'CrIS gran',igran
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

   VIIRS_ISCCP_Hist=fltarr(nbin_tau,nbin_ctp,9,30,45)*!values.f_nan 
   VIIRS_COD_Hist=fltarr(nbin_tau,9,30,45)*!values.f_nan
   VIIRS_CTP_Hist=fltarr(nbin_CTP,9,30,45)*!values.f_nan
   VIIRS_CLD_Mask=fltarr(4,9,30,45)*!values.f_nan
   VIIRS_CTP_moment=fltarr(4,9,30,45)*!values.f_nan
   VIIRS_COD_moment=fltarr(4,9,30,45)*!values.f_nan
   VIIRS_lat=fltarr(9,30,45)*!values.f_nan
   VIIRS_lon=fltarr(9,30,45)*!values.f_nan
   
   find=where(dy GE bnd1 and dy LE bnd2,nfind) 
   dysizegood=dysize*0 & dysizegood[find]=dysize[find]
   for k=0,8 do begin
       for j=0,29 do begin
           for i=0,44 do begin
               if dysizegood[k,j,i] lT 1 then continue
               offset1=k+9*(j+30*(i+45l*igran)) 
                  if offset1 GE 1 then begin
                     dy=index.dy[ysize2[offset1-1]:ysize2[offset1]-1]-bnd1
                     dx=index.dx[ysize2[offset1-1]:ysize2[offset1]-1]
                  endif else begin
                     dy=index.dy[0:ysize2[offset1]-1]-bnd1
                     dx=index.dx[0:ysize2[offset1]-1]
                  endelse
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
   save,VIIRS_lon,VIIRS_lat,VIIRS_ISCCP_Hist,VIIRS_COD_Hist,VIIRS_CTP_Hist,VIIRS_CLD_Mask,VIIRS_CTP_moment,VIIRS_COD_moment,$
        filename='VIIRS_Cloud_onCrIS_20150601_G'+string(igran+1,format='(i3.3)')+'.sav'
endfor

stop

end


       
           
