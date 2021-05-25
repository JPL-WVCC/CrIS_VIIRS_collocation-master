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




for ifile=0,0 do begin
    indexfile=indexdir+'IND_CrIS_VIIRSMOD_'+mode+'_'+strcompress(string(10*ifile),/remove_all)+'.nc'
    indexfile2=indexdir+'test_IND_CrIS_VIIRSMOD_20150115_1.nc'
    ierr=read_ncdf(indexfile,index)
    ngran=10
   print,'index_file',indexfile
    result=size(index.dy_size,/dimension)
    dysize=index.dy_size
    ysize2=total(reform(index.dy_size,result[0]*result[1]*result[2]),/cumulative,/integer)   ;cumulatively the end index of dy in each CrIS pixel

   ierr=read_ncdf(indexfile2,index1)
   result1=size(index1.dy_size,/dimension)
   spawn,'ls '+CrISL1dir+'SNDR*.g002.L1B*.nc',crisL2file
   ierr=read_ncdf(crisL2file,data1,variable_list=['lat','lon'])
   cris_lat1=data1.lat
   cris_lon1=data1.lon
   data1=0
   ierr=read_nc4('/raid15/qyue/VIIRS/VIIRS/201501/VNP03MOD/VNP03MOD.A2015015.0006.001.2017259084650.nc',data,groupname='geolocation_data',content_list=['latitude','longitude'])
   ctp1=data.latitude
   cod1=data.longitude
   data=0
   ysize21=total(reform(index1.dy_size,result1[0]*result1[1]*result1[2]),/cumulative,/integer)   ;cumulatively the end index of dy in each CrIS pixel
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
          ierr=read_nc4(viirsfile,data,groupname='geolocation_data',$
              content_list=['latitude','longitude'])
          if okfile eq 0 then begin
             ctp=data.latitude 
             cod=data.longitude
             data=0 & data1=0       
          endif
          if okfile GE 1 then begin  
    
             ctp=concat(ctp,data.latitude,dimension=1)
             cod=concat(cod,data.longitude,dimension=1)
              data=0 & data1=0
          endif
   endfor
   


 ;  find=where(cris_lon lt 0,nfind)
 ;  if nfind GE 1 then cris_lon[find]=360+cris_lon[find]
 ;  find=where(cod lt 0,nfind)
 ;  if nfind GE 1 then cod[find]=360+cod[find] 
   dysize=index.dy_size
   findcol=where(dysize GE 1,nfindcol)
   col = findcol mod result[0]
   row = (findcol / result[0]) mod result[1]
   frame = findcol / (result[1]*result[0])

;   cgmap_set,/stereographic, (60+75)/2,(150+210)/2,/continents,limit=[60,150,75,210],/grid
   !p.multi=[0,3,2]
   for i=0,nfindcol-1 do begin
       time=systime(/seconds)
       if frame[i] EQ 88 then begin

       if row[i] eq 0 or row[i] eq 14 or row[i] eq 29 then begin
 
       if i eq 0 and findcol[0] eq 0 then begin
            dy=index.dy[0:ysize2[0]-1] & dx=index.dx[0:ysize2[0]-1]
       endif else begin
            dy=index.dy[ysize2[findcol[i-1]]:ysize2[findcol[i]]-1]
            dx=index.dx[ysize2[findcol[i-1]]:ysize2[findcol[i]]-1]
       endelse 
       print,col[i],row[i],frame[i]     
       lon=cod[dx,dy] & lat=ctp[dx,dy]
       limit=[min(lat),min(lon),max(lat),max(lon)]
       if col[i] eq 0 then cgplot,[min(lon),max(lon)],[min(lat),max(lat)],/nodata,xtitle='Lon',ytitle='Lat',charsize=2,$
xrange=[min(cris_lon[*,row[i],frame[i]])-0.6,max(cris_lon[*,row[i],frame[i]])+0.6],yrange=[min(cris_lat[*,row[i],frame[i]])-0.3,max(cris_lat[*,row[i],frame[i]])+0.3]
       cgplot,lon,lat,psym=3,/overplot,symsize=0.1
       cgplot,cris_lon[col[i],row[i],frame[i]],cris_lat[col[i],row[i],frame[i]],psym=9,color='red',/overplot
 
       endif
        
      endif
   endfor

  ; find=where(cris_lon1 lt 0,nfind)
  ; if nfind GE 1 then cris_lon1[find]=360+cris_lon1[find]
  ; find=where(cod1 lt 0,nfind)
  ; if nfind GE 1 then cod1[find]=360+cod1[find] 
   dysize1=index1.dy_size
   findcol=where(dysize1 GE 1,nfindcol)
   col = findcol mod result1[0]
   row = (findcol / result1[0]) mod result1[1]
   frame = findcol / (result1[1]*result1[0])

;   cgmap_set,/stereographic, (60+75)/2,(150+210)/2,/continents,limit=[60,150,75,210],/grid

   for i=0,nfindcol-1 do begin
       time=systime(/seconds)
       if frame[i] EQ 43 then begin

       if row[i] eq 0 or row[i] eq 14 or row[i] eq 29 then begin
 
       if i eq 0 and findcol[0] eq 0 then begin
            dy=index1.dy[0:ysize21[0]-1] & dx=index1.dx[0:ysize21[0]-1]
       endif else begin
            dy=index1.dy[ysize21[findcol[i-1]]:ysize21[findcol[i]]-1]
            dx=index1.dx[ysize21[findcol[i-1]]:ysize21[findcol[i]]-1]
       endelse 
       print,col[i],row[i],frame[i]     
       lon=cod1[dx,dy] & lat=ctp1[dx,dy]
       limit=[min(lat),min(lon),max(lat),max(lon)]
       if col[i] eq 0 then cgplot,[min(lon),max(lon)],[min(lat),max(lat)],/nodata,xtitle='Lon',ytitle='Lat',charsize=2,$
xrange=[min(cris_lon1[*,row[i],frame[i]])-0.6,max(cris_lon1[*,row[i],frame[i]])+0.6],yrange=[min(cris_lat1[*,row[i],frame[i]])-0.3,max(cris_lat1[*,row[i],frame[i]])+0.3]
       cgplot,lon,lat,psym=3,/overplot,symsize=0.1
       cgplot,cris_lon1[col[i],row[i],frame[i]],cris_lat1[col[i],row[i],frame[i]],psym=9,color='red',/overplot
 
       endif
        
      endif
   endfor
   
endfor

endfor

stop

end


       
           
