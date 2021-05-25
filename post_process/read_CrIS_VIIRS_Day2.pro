pro read_cirs_viirs,year,st_day,ed_day

 cpu, tpool_nthreads=1

mod_taubin=[0.,0.3,1.3,3.6,9.4,23,60,200]
mod_ctpbin=([50,180,310,440,560,680,800,1100])  
nbin_tau=n_elements(mod_taubin)-1  
nbin_ctp=n_elements(mod_ctpbin)-1

indexdir='/raid15/qyue/VIIRS/VIIRS/201501/Index/'

CrISL1dir0='/peate_archive/.data1/Ops/snpp/gdisc/2/'
CrISL2dir='/peate_archive/.data7/Ops/sndr/production/v01_23_01/2015/06/01/l2climcaps/'
;variable_list=['obs_time_tai93','lon','lat','air_temp','air_temp_qc','spec_hum','spec_hum_qc','air_pres_h2o','air_pres',$
;               'land_frac','sol_zen','air_pres_nsurf','air_pres_h2o_nsurf','air_pres_lay','air_temp_dof','h2o_vap_dof',$
;               'rel_hum','rel_hum_qc','cld_frac','cld_frac_qc','cld_top_pres','cld_top_pres_qc','mw_surf_class']
;content_list=['for_cld_frac_tot','for_cld_frac_tot_qc','for_cld_top_pres_tot','for_cld_top_pres_tot_qc']
;clist=['fg_air_temp','fg_h2o_vap_mol_lay','fg_surf_air_temp']
;mollist=['h2o_vap_mol_lay','h2o_vap_mol_lay_qc']

VIIRSL2dir='/raid15/qyue/VIIRS/VIIRS/201501/CLDPROP/'

for iday=st_day,ed_day do begin

CrISL1dir=CrISL1dir0+'2015/01/'+string(iday,format='(i2.2)')+'/crisl1b/'

mode=string(year,format='(i4.4)')+'01'+string(iday,format='(i2.2)')
mode1=string(year,format='(i4.4)')+string(iday,format='(i3.3)')

VIIRSL2dir='/raid15/qyue/VIIRS/VIIRS/201501/CLDPROP/CLDPROP_L2_VIIRS_SNPP/2015/'+string(iday,format='(i3.3)')+'/'
spawn,'ls '+VIIRSL2dir+'CLDPROP_L2_VIIRS_SNPP.A'+mode1+'*.nc',VIIRScldfiles


for ifile=0,23 do begin
    indexfile=indexdir+'IND_CrIS_VIIRSMOD_'+mode+'_'+strcompress(string(10*ifile),/remove_all)+'.nc'
    ierr=read_ncdf(indexfile,index)
    ngran=10
   print,'index_file',indexfile
    result=size(index.dy_size,/dimension)
    dysize=index.dy_size
    ysize2=total(reform(index.dy_size,result[0]*result[1]*result[2]),/cumulative,/integer)   ;cumulatively the end index of dy in each CrIS pixel

   okfile=-1
   for igran=0,ngran-1 do begin
       ;spawn,'ls '+CrISL2dir+'SNDR*.g'+string(ifile*ngran+igran+1,format='(i3.3)')+'.L2_RET*.nc',crisL2file
       spawn,'ls '+CrISL1dir+'SNDR*.g'+string(ifile*ngran+igran+1,format='(i3.3)')+'.L1B*.nc',crisL2file
       okfile=okfile+1
       ierr=read_ncdf(crisL2file,data1,variable_list=['lat','lon']) 
          if okfile eq 0 then begin
             CrIS_lat=data1.lat
             CrIS_lon=data1.lon
          endif
          if okfile GE 1 then begin
             CrIS_lat=concat(CrIS_lat,data1.lat,dimension=2)
             CrIS_lon=concat(CrIS_lon,data1.lon,dimension=2)
          endif
   endfor   
   okfile=-1   
   for iVgran=0,ngran-1 do begin 

          okfile=okfile+1 
          viirsfile=viirscldfiles[ifile*ngran+ivgran] ;#of ok VIIRS files matched to this CrIS granule        
          ierr=read_nc4(viirsfile,data,groupname='geophysical_data',$
              content_list=['Cloud_Top_Pressure','Cloud_Optical_Thickness','Cloud_Mask','Cloud_Effective_Radius'])
          if okfile eq 0 then begin
             ctp=data.cloud_top_pressure*0.1 
             cod=data.cloud_optical_thickness*0.01
             Re=data.cloud_effective_radius*0.01
             cmaskt=data.cloud_mask
             data=0 & data1=0       
          endif
          if okfile GE 1 then begin  
    
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
   VIIRS_CLD_Mask=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_CTP_moment=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_COD_moment=fltarr(4,result[0],result[1],result[2])*!values.f_nan
   VIIRS_Re_ISCCPHist=fltarr(nbin_tau,nbin_ctp,result[0],result[1],result[2])*!values.f_nan

   
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
       for im=0,3 do begin
           find=where(cmask[dx,dy] eq im,nfind)
           VIIRS_CLD_Mask[im,col[i],row[i],frame[i]]=nfind
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

   endfor
   save,CrIS_lon,CrIS_lat,VIIRS_ISCCP_Hist,VIIRS_Re_ISCCPHist,VIIRS_CLD_Mask,VIIRS_CTP_moment,VIIRS_COD_moment,$
        /compress,filename='VIIRS_Cloud_onCrIS_'+mode+'_F'+string(ifile,format='(i3.3)')+'.sav'
   
endfor

endfor

stop

end


       
           
