pro read_CrISL2_VIIRS,year,month,st_day,ed_day

syear=string(year,format='(i4.4)')
smonth=string(month,format='(i2.2)')

rootj='/peate_archive/NPPOps/sndr/production/v01_33_00/'
rootc='/peate_archive/NPPOps/sndr/production/v01_33_00/'

viirsdir='/home/qyue/VIIRS/CrIS_VIIRS_collocation-master/post_process/VIIRS_Cloud_onCrIS/'
ECMdir='/peate_sci/eff/ecmwfAnalysis/cris/'

ngran=10
nbin_tau=7 & nbin_ctp=7

cal=0

if cal eq 1 then begin

variable_list=['obs_time_tai93','lon','lat','air_temp','air_temp_qc','spec_hum','spec_hum_qc','mw_surf_class','air_pres_h2o','air_pres',$
               'land_frac','sol_zen','air_pres_nsurf','air_pres_h2o_nsurf','air_pres_lay','rel_hum','rel_hum_qc'];,'o3_mmr','o3_mmr_qc']
content_list=['for_cld_frac_tot','for_cld_frac_tot_qc','for_cld_top_pres_tot','for_cld_top_pres_tot_qc']
clist=['fg_air_temp','fg_h2o_vap_mol_lay']
jlist=['nn_air_temp','nn_h2o_vap_mol_lay']
mollist=['h2o_vap_mol_lay','h2o_vap_mol_lay_qc']

for iday=st_day,ed_day do begin
    print,iday
    sday=string(iday,format='(i2.2)')
        
    SIPSJ=rootj+syear+'/'+smonth+'/'+sday+'/snpp_l2_chart_ret/'
    SIPSC=rootc+syear+'/'+smonth+'/'+sday+'/snpp_l2_climcaps_ret/'

    restore,ECMdir+syear+'/'+smonth+'/ecmwf_analysis_cris_'+syear+'-'+smonth+'-'+sday+'.dat'
    
    q=xe.ph2o*1000. & t=xe.ptemp 
    
    for ifile=0,23 do begin
        restore,viirsdir+'VIIRS_Cloud_onCrIS_'+syear+smonth+sday+'_F'+string(ifile,format='(i3.3)')+'.sav'
        result=size(VIIRS_ISCCP_Hist,/dimensions)
        if result[4] ne 450 then begin
             print,result,'VIIRS_Cloud_onCrIS_'+syear+smonth+sday+'_F'+string(ifile,format='(i3.3)')+'.sav', ' not right'
             continue
        endif
        print,'VIIRS_Cloud_onCrIS_'+syear+smonth+sday+'_F'+string(ifile,format='(i3.3)')+'.sav'
;        CrIS_lon,CrIS_lat,VIIRS_ISCCP_Hist,VIIRS_Re_ISCCPHist,VIIRS_CLD_Mask,VIIRS_CTP_moment,VIIRS_COD_moment
        VIIRS_ISCCP_Hist=reform(total(VIIRS_ISCCP_Hist,3,/nan),nbin_tau,nbin_ctp,30,450)
        VIIRS_ISCCP_Histpro=reform(transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[0,1,4,2,3]),nbin_tau,nbin_ctp,30*450l*66)

    
        for igran=0,ngran-1 do begin
           granid=ifile*ngran+igran+1
           spawn,'ls '+SIPSJ+'SNDR*.g'+string(granid,format='(i3.3)')+'.L2_CHART_RET*.nc',filej
           spawn,'ls '+SIPSC+'SNDR*.g'+string(granid,format='(i3.3)')+'.L2_CLIMCAPS_RET*.nc',filec 
           if filej[0] eq '' or filec[0] eq '' then message,'stop missing'
           ierr=read_ncdf(filej[0],dataj,variable_list=variable_list)
           ierr=read_nc4(filej[0],dataj1,groupname='aux',content_list=[content_list,jlist])
           ierr=read_nc4(filej[0],dataj2,groupname='mol_lay',content_list=mollist)
           find=where(dataj1.for_cld_frac_tot_qc gt 1,nfind) 
           if nfind GE 1 then dataj1.for_cld_frac_tot[find]=!values.f_nan
           find=where(dataj1.for_cld_top_pres_tot_qc gt 1,nfind) 
           if nfind GE 1 then dataj1.for_cld_top_pres_tot[find]=!values.f_nan      
           jcf=dataj1.for_cld_frac_tot & jctp=dataj1.for_cld_top_pres_tot
       
           ierr=read_ncdf(filec[0],datac,variable_list=variable_list)
           ierr=read_nc4(filec[0],datac1,groupname='aux',content_list=[content_list,clist])
           ierr=read_nc4(filec[0],datac2,groupname='mol_lay',content_list=mollist)
           find=where(datac1.for_cld_frac_tot_qc gt 1,nfind) 
           if nfind GE 1 then datac1.for_cld_frac_tot[find]=!values.f_nan
           find=where(datac1.for_cld_top_pres_tot_qc gt 1,nfind) 
           if nfind GE 1 then datac1.for_cld_top_pres_tot[find]=!values.f_nan
           ccf=datac1.for_cld_frac_tot & cctp=datac1.for_cld_top_pres_tot

       ;     jtai=reform(dataj.obs_time_tai93,30*45)
       ;    jlon=reform(dataj.lon,30*45)
           jlat=reform(dataj.lat,30*45)


;           ecbuffer=!null
;           ECMWF_INTERP_llt,tai,lon,lat,ecbuffer

      ;  if ecbuffer eq !null then continue
           EC_qx=reform(q[34:99,*,*,granid-1]) & EC_Tx=reform(t[34:99,*,*,granid-1]) 
           find=where(EC_qx lt 0,nfind) & if nfind GE 1 then EC_qx[find]=!values.f_nan
           find=where(EC_tx lt 0,nfind) & if nfind GE 1 then EC_tx[find]=!values.f_nan
           

           pres_h2o=dataj.air_pres_h2o/100. & pres_t=dataj.air_pres/100. & pres_lay=dataj.air_pres_lay/100.
 
           tempq=reform(dataj2.h2o_vap_mol_lay,100,30*45)*1.0e-4 & tempt=reform(dataj.air_temp,100,30*45)  
           tempqa=reform(dataj1.nn_h2o_vap_mol_lay,100,30*45)*1.0e-4  & tempta=reform(dataj1.nn_air_temp,100,30*45) 
           find=where(tempq GE 1.0e30,nbad) & if nbad GE 1 then tempq[find]=!values.f_nan 
           find=where(tempqa GE 1.0e30,nbad) & if nbad GE 1 then tempqa[find]=!values.f_nan
           find=where(tempt GE 1.0e30,nbad) & if nbad GE 1 then tempt[find]=!values.f_nan 
           find=where(tempta GE 1.0e30,nbad) & if nbad GE 1 then tempta[find]=!values.f_nan

            tempq1=reform(datac2.h2o_vap_mol_lay,100,30*45)*1.0e-4   & tempt1=reform(datac.air_temp,100,30*45)
           tempqa1=reform(datac1.fg_h2o_vap_mol_lay,100,30*45)*1.0e-4  & tempta1=reform(datac1.fg_air_temp,100,30*45)
           find=where(tempq1 GE 1.0e30,nbad) & if nbad GE 1 then tempq1[find]=!values.f_nan 
           find=where(tempqa1 GE 1.0e30,nbad) & if nbad GE 1 then tempqa1[find]=!values.f_nan 
           find=where(tempt1 GE 1.0e30,nbad) & if nbad GE 1 then tempt1[find]=!values.f_nan 
           find=where(tempta1 GE 1.0e30,nbad) & if nbad GE 1 then tempta1[find]=!values.f_nan

           for ipix=0,30*45l-1 do begin
            result=prof_interpol(tempq[*,ipix],pres_lay,pres_t,method=2,fillvalue=!values.f_nan)
            tempq[*,ipix]=cd_to_mmr(pres_t[*],tempt[*,ipix],result,jlat[ipix])
            result=prof_interpol(tempqa[*,ipix],pres_lay,pres_t,method=2,fillvalue=!values.f_nan)
            tempqa[*,ipix]=cd_to_mmr(pres_t[*],tempta[*,ipix],result,jlat[ipix])
            result=prof_interpol(tempq1[*,ipix],pres_lay,pres_t,method=2,fillvalue=!values.f_nan)
            tempq1[*,ipix]=cd_to_mmr(pres_t[*],tempt1[*,ipix],result,jlat[ipix])
            result=prof_interpol(tempqa1[*,ipix],pres_lay,pres_t,method=2,fillvalue=!values.f_nan)
            tempqa1[*,ipix]=cd_to_mmr(pres_t[*],tempta1[*,ipix],result,jlat[ipix])
           endfor

           tempq=reform(tempq,100,30,45) & tempqa=reform(tempqa,100,30,45)
           tempq1=reform(tempq1,100,30,45) & tempqa1=reform(tempqa,100,30,45)

           if igran eq 0 then begin
           lf=dataj.land_frac & solzen=dataj.sol_zen 
           j_cf=dataj1.for_cld_frac_tot & j_ctp=dataj1.for_cld_top_pres_tot
           c_cf=datac1.for_cld_frac_tot & c_ctp=datac1.for_cld_top_pres_tot
           j_sh=tempq[34:99,*,*] & j_sh_qc=dataj.spec_hum_qc
           j_rh=dataj.rel_hum & j_rh_qc=dataj.rel_hum_qc
           j_t=dataj.air_temp[34:99,*,*] & j_t_qc=dataj.air_temp_qc[34:99,*,*]
           j_sh_fg=tempqa[34:99,*,*] & j_t_fg=tempta[34:99,*,*]
           c_sh=tempq1[34:99,*,*] & c_sh_qc=datac.spec_hum_qc
           c_rh=datac.rel_hum & c_rh_qc=datac.rel_hum_qc
           c_t=datac.air_temp[34:99,*,*] & c_t_qc=datac.air_temp_qc[34:99,*,*]
           c_sh_fg=tempqa1[34:99,*,*] & c_t_fg=tempta1[34:99,*,*]
           EC_T=EC_Tx & EC_q=EC_qx
           endif
           if igran gt 0 then begin          
           lf=[[lf],[dataj.land_frac]] & solzen=[[solzen],[dataj.sol_zen]]
           j_cf=[[j_cf],[dataj1.for_cld_frac_tot]] & j_ctp=[[j_ctp],[dataj1.for_cld_top_pres_tot]] 
           c_cf=[[c_cf],[datac1.for_cld_frac_tot]] & c_ctp=[[c_ctp],[datac1.for_cld_top_pres_tot]]
 
           j_sh=[[[j_sh]],[[tempq[34:99,*,*]]]] & j_sh_qc=[[[j_sh_qc]],[[dataj.spec_hum_qc]]]
           j_rh=[[[j_rh]],[[dataj.rel_hum]]] & j_rh_qc=[[[j_rh_qc]],[[dataj.rel_hum_qc]]]
           j_t=[[[j_t]],[[dataj.air_temp[34:99,*,*]]]] & j_t_qc=[[[j_t_qc]],[[dataj.air_temp_qc[34:99,*,*]]]]
           j_sh_fg=[[[j_sh_fg]],[[tempqa[34:99,*,*]]]]
           j_t_fg=[[[j_t_fg]],[[tempta[34:99,*,*]]]]
           c_sh=[[[c_sh]],[[tempq1[34:99,*,*]]]] & c_sh_qc=[[[c_sh_qc]],[[datac.spec_hum_qc]]]
           c_rh=[[[c_rh]],[[datac.rel_hum]]] & c_rh_qc=[[[c_rh_qc]],[[datac.rel_hum_qc]]]
           c_t=[[[c_t]],[[datac.air_temp[34:99,*,*]]]] & c_t_qc=[[[c_t_qc]],[[datac.air_temp_qc[34:99,*,*]]]]
           c_sh_fg=[[[c_sh_fg]],[[tempqa1[34:99,*,*]]]]
           c_t_fg=[[[c_t_fg]],[[tempta1[34:99,*,*]]]]
           EC_T=[[[EC_T]],[[EC_Tx]]] & EC_q=[[[EC_q]],[[EC_qx]]]
           endif
      
           EC_Tj=EC_T & VIIRS_Histj=VIIRS_ISCCP_Histpro
           find=where(j_t_qc gt 1,nfind) 
           if nfind GE 1 then begin
               j_t[find]=!values.f_nan & EC_Tj[find]=!values.f_nan & VIIRS_Histj[*,*,find]=0
               
           endif
           EC_qj=EC_q
           find=where(j_sh_qc gt 1,nfind) 
           if nfind GE 1 then begin
               j_sh[find]=!values.f_nan & EC_qj[find]=!values.f_nan
           endif          
       
           EC_Tc=EC_T & VIIRS_Histc=VIIRS_ISCCP_Histpro
           find=where(c_t_qc gt 1,nfind) 
           if nfind GE 1 then begin
               c_t[find]=!values.f_nan & EC_Tc[find]=!values.f_nan & VIIRS_Histc[*,*,find]=0
           endif
           EC_qc=EC_q
           find=where(c_sh_qc gt 1,nfind) 
           if nfind GE 1 then begin
               c_sh[find]=!values.f_nan & EC_qc[find]=!values.f_nan
           endif 
      
           find=where(EC_T ne EC_T,nfind)
           if nfind GE 1 then begin
               VIIRS_ISCCP_Histpro[*,*,find]=0
           endif
 
           
        endfor ;all 10 granules 
 
           EC_T_ISCCPx=reform(total(total(transpose(rebin(EC_T,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp)
           EC_q_ISCCPx=reform(total(total(transpose(rebin(EC_q,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp)  
           EC_Tj_ISCCPx=reform(total(total(transpose(rebin(EC_Tj,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp) 
           EC_Tc_ISCCPx=reform(total(total(transpose(rebin(EC_Tc,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp)
           EC_qj_ISCCPx=reform(total(total(transpose(rebin(EC_qj,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp)
           EC_qc_ISCCPx=reform(total(total(transpose(rebin(EC_qc,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp) 
           j_T_ISCCPx=reform(total(total(transpose(rebin(j_T,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp) 
           j_q_ISCCPx=reform(total(total(transpose(rebin(j_sh,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp)
           c_T_ISCCPx=reform(total(total(transpose(rebin(c_T,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp) 
           c_q_ISCCPx=reform(total(total(transpose(rebin(c_sh,66,30,450,nbin_tau,nbin_ctp,/sample),[0,3,4,1,2])* $
                      transpose(rebin(VIIRS_ISCCP_Hist,nbin_tau,nbin_ctp,30,450,66,/sample),[4,0,1,2,3]),4,/nan),4,/nan),66,nbin_tau*nbin_ctp)
           VIIRS_ISCCPxj=reform(transpose(total(total(reform(VIIRS_Histj,nbin_tau,nbin_ctp,66,30,450),4,/nan),4,/nan),[2,0,1]),66,nbin_tau*nbin_ctp)
           VIIRS_ISCCPxc=reform(transpose(total(total(reform(VIIRS_Histc,nbin_tau,nbin_ctp,66,30,450),4,/nan),4,/nan),[2,0,1]),66,nbin_tau*nbin_ctp)
           VIIRS_ISCCPx=reform(transpose(total(total(reform(VIIRS_ISCCP_Histpro,nbin_tau,nbin_ctp,66,30,450),4,/nan),4,/nan),[2,0,1]),66,nbin_tau*nbin_ctp)
    
           if ifile eq 0 and iday eq st_day then begin
              EC_T_ISCCP=EC_T_ISCCPx & EC_q_ISCCP=EC_q_ISCCPx
              EC_Tj_ISCCP=EC_Tj_ISCCPx & EC_qj_ISCCP=EC_qj_ISCCPx
              EC_Tc_ISCCP=EC_Tc_ISCCPx & EC_qc_ISCCP=EC_qc_ISCCPx
              j_T_ISCCP=j_T_ISCCPx & j_q_ISCCP=j_q_ISCCPx
              c_T_ISCCP=c_T_ISCCPx & c_q_ISCCP=c_q_ISCCPx
              VIIRS_ISCCP=VIIRS_ISCCPx
              VIIRS_ISCCPj=VIIRS_ISCCPxj
              VIIRS_ISCCPc=VIIRS_ISCCPxc
           endif else begin
              EC_T_ISCCP=[[[EC_T_ISCCP]],[[EC_T_ISCCPx]]] & EC_q_ISCCP=[[[EC_q_ISCCP]],[[EC_q_ISCCPx]]]
              EC_Tj_ISCCP=[[[EC_Tj_ISCCP]],[[EC_Tj_ISCCPx]]] & EC_qj_ISCCP=[[[EC_qj_ISCCP]],[[EC_qj_ISCCPx]]]
              EC_Tc_ISCCP=[[[EC_Tc_ISCCP]],[[EC_Tc_ISCCPx]]] & EC_qc_ISCCP=[[[EC_qc_ISCCP]],[[EC_qc_ISCCPx]]]
              j_T_ISCCP=[[[j_T_ISCCP]],[[j_T_ISCCPx]]] & j_q_ISCCP=[[[j_q_ISCCP]],[[j_q_ISCCPx]]]
              c_T_ISCCP=[[[c_T_ISCCP]],[[c_T_ISCCPx]]] & c_q_ISCCP=[[[c_q_ISCCP]],[[c_q_ISCCPx]]]
              VIIRS_ISCCP=[[[VIIRS_ISCCP]],[[VIIRS_ISCCPx]]]
              VIIRS_ISCCPj=[[[VIIRS_ISCCPj]],[[VIIRS_ISCCPxj]]]
              VIIRS_ISCCPc=[[[VIIRS_ISCCPc]],[[VIIRS_ISCCPxc]]]
           endelse
      endfor
 endfor

 VIIRS_ISCCP=total(VIIRS_ISCCP,3,/nan)  ;66,9
 VIIRS_ISCCPj=total(VIIRS_ISCCPj,3,/nan)  ;66,9
 VIIRS_ISCCPc=total(VIIRS_ISCCPc,3,/nan)  ;66,9
 EC_T_ISCCP=total(EC_T_ISCCP,3,/nan)/VIIRS_ISCCP
 EC_q_ISCCP=total(EC_q_ISCCP,3,/nan)/VIIRS_ISCCP 
 EC_Tj_ISCCP=total(EC_Tj_ISCCP,3,/nan)/VIIRS_ISCCPj
 EC_qj_ISCCP=total(EC_qj_ISCCP,3,/nan)/VIIRS_ISCCPj 
 EC_Tc_ISCCP=total(EC_Tc_ISCCP,3,/nan)/VIIRS_ISCCPc
 EC_qc_ISCCP=total(EC_qc_ISCCP,3,/nan)/VIIRS_ISCCPc
 j_T_ISCCP=total(j_T_ISCCP,3,/nan)/VIIRS_ISCCPj
 j_q_ISCCP=total(j_q_ISCCP,3,/nan)/VIIRS_ISCCPj  
 c_T_ISCCP=total(c_T_ISCCP,3,/nan)/VIIRS_ISCCPc
 c_q_ISCCP=total(c_q_ISCCP,3,/nan)/VIIRS_ISCCPc


save, pres_h2o,VIIRS_ISCCP,VIIRS_ISCCPj,VIIRS_ISCCPc,EC_T_ISCCP,EC_q_ISCCP,EC_Tj_ISCCP,EC_qj_ISCCP,EC_Tc_ISCCP,EC_qc_ISCCP,j_T_ISCCP,j_q_ISCCP,c_T_ISCCP,c_q_ISCCP,$
      /compress,filename=viirsdir+'SNPP_MeanProfile_ISCCP_'+syear+smonth+'.sav'

endif

restore,viirsdir+'SNPP_MeanProfile_ISCCP_'+syear+smonth+'.sav'

j_T_ISCCP[65,*]=!values.f_nan
j_q_ISCCP[65,*]=!values.f_nan
c_T_ISCCP[65,*]=!values.f_nan
c_q_ISCCP[65,*]=!values.f_nan

 VIIRS_ISCCP=reform(VIIRS_ISCCP,66,nbin_tau,nbin_ctp)  ;66,9
 VIIRS_ISCCPj=reform(VIIRS_ISCCPj,66,nbin_tau,nbin_ctp)  ;66,9
 VIIRS_ISCCPc=reform(VIIRS_ISCCPc,66,nbin_tau,nbin_ctp)  ;66,9
 EC_T_ISCCP=reform(EC_T_ISCCP,66,nbin_tau,nbin_ctp)
 EC_q_ISCCP=reform(EC_q_ISCCP,66,nbin_tau,nbin_ctp) 
 EC_Tj_ISCCP=reform(EC_Tj_ISCCP,66,nbin_tau,nbin_ctp)
 EC_qj_ISCCP=reform(EC_qj_ISCCP,66,nbin_tau,nbin_ctp) 
 EC_Tc_ISCCP=reform(EC_Tc_ISCCP,66,nbin_tau,nbin_ctp)
 EC_qc_ISCCP=reform(EC_qc_ISCCP,66,nbin_tau,nbin_ctp)
 j_T_ISCCP=reform(j_T_ISCCP,66,nbin_tau,nbin_ctp)
 j_q_ISCCP=reform(j_q_ISCCP,66,nbin_tau,nbin_ctp)  
 c_T_ISCCP=reform(c_T_ISCCP,66,nbin_tau,nbin_ctp)
 c_q_ISCCP=reform(c_q_ISCCP,66,nbin_tau,nbin_ctp)

yieldj=float(VIIRS_ISCCPj)/VIIRS_ISCCP*100.
yieldc=float(VIIRS_ISCCPc)/VIIRS_ISCCP*100.



psname=viirsdir+'SNPP_MeanProfile_ISCCP_'+syear+smonth

cgps_open,psname+'.ps'
cgdisplay,1200,1200
pos=cglayout([7,7],OXMargin=[1, 1], OYMargin=[6, 14], XGap=2, YGap=2)
charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.25 : cgDefCharsize()*0.25

taulabel = ['< 0.3','0.3-1.3', '1.3-3.6', '3.6-9.4', '9.4-23.', '23.-60.','>60']
ctpbin=strcompress(string([180,310,440,560,680,800,1100]),/remove_all)

thick=5
xrange=[-5,5]
cgarrow,-0.03,0.04,-0.03,0.95,linestyle=0,color='black',thick=3,/normal
cgarrow,-0.03,0.04,1.,0.04,linestyle=0,color='black',thick=3,/normal
cgtext,-0.08,0.9,textoidl('CTP'),/normal,charsize=charsize*2
cgtext,0.95,-0.01,textoidl('\tau'),/normal,charsize=charsize*2

cgtext,0.5,0.9,'CHART',color='red',/normal,charsize=charsize*2
cgtext,0.6,0.9,'CLIMCAPS',color='blue',/normal,charsize=charsize*2

for ictp=0,nbin_ctp-1 do begin
    for itau=0,nbin_tau-1 do begin

    position=pos[*,itau+nbin_tau*ictp]
   
    if itau eq 0 then ytitle='CTP '+ctpbin[ictp] else ytitle=''
    if ictp eq nbin_ctp-1 then xtitle='T (K)' else xtitle=''
    if ictp eq 0 then title=textoidl('\tau')+taulabel[itau] else title=''

    cgplot,EC_T_ISCCP[*,itau,ictp],pres_h2o,color='black',thick=thick,yrange=[1100,100],title=title,$
     charsize=charsize,xrange=xrange,position=position,/noerase,background=background,xtitle=xtitle,ytitle=ytitle,/nodata
    
    cgplot,j_T_ISCCP[*,itau,ictp]-EC_T_ISCCP[*,itau,ictp],pres_h2o,color='red',thick=thick,/overplot
    cgplot,EC_Tj_ISCCP[*,itau,ictp]-EC_T_ISCCP[*,itau,ictp],pres_h2o,color='red',thick=thick,/overplot,linestyle=2
    cgplot,c_T_ISCCP[*,itau,ictp]-EC_T_ISCCP[*,itau,ictp],pres_h2o,color='blue',thick=thick,/overplot
    cgplot,EC_Tc_ISCCP[*,itau,ictp]-EC_T_ISCCP[*,itau,ictp],pres_h2o,color='blue',thick=thick,/overplot,linestyle=2
    
    endfor

endfor
cgText, 0.5, 0.95, /Normal,'Temperature (K): CrIMSS-ECMWF',Charsize=cgDefCharsize()*1.25, Alignment=0.5

erase

thick=5
xrange=[-5,5]
cgarrow,-0.03,0.04,-0.03,0.95,linestyle=0,color='black',thick=3,/normal
cgarrow,-0.03,0.04,1.,0.04,linestyle=0,color='black',thick=3,/normal
cgtext,-0.08,0.9,textoidl('CTP'),/normal,charsize=charsize*2
cgtext,0.95,-0.01,textoidl('\tau'),/normal,charsize=charsize*2

cgtext,0.5,0.9,'CHART',color='red',/normal,charsize=charsize*2
cgtext,0.6,0.9,'CLIMCAPS',color='blue',/normal,charsize=charsize*2

for ictp=0,nbin_ctp-1 do begin
    for itau=0,nbin_tau-1 do begin

    position=pos[*,itau+nbin_tau*ictp]
   
    if itau eq 0 then ytitle='CTP '+ctpbin[ictp] else ytitle=''
    if ictp eq nbin_ctp-1 then xtitle='T (K)' else xtitle=''
    if ictp eq 0 then title=textoidl('\tau')+taulabel[itau] else title=''

    cgplot,EC_T_ISCCP[*,itau,ictp],pres_h2o,color='black',thick=thick,yrange=[1100,100],title=title,$
     charsize=charsize,xrange=xrange,position=position,/noerase,background=background,xtitle=xtitle,ytitle=ytitle,/nodata
    
    cgplot,j_T_ISCCP[*,itau,ictp]-EC_Tj_ISCCP[*,itau,ictp],pres_h2o,color='red',thick=thick,/overplot
    cgplot,c_T_ISCCP[*,itau,ictp]-EC_Tc_ISCCP[*,itau,ictp],pres_h2o,color='blue',thick=thick,/overplot
    
    endfor

endfor
cgText, 0.5, 0.95, /Normal,'Temperature (K): CrIMSS-ECMWF(QC)',Charsize=cgDefCharsize()*1.25, Alignment=0.5

erase

xrange=[-60,60]
cgarrow,-0.03,0.04,-0.03,0.95,linestyle=0,color='black',thick=3,/normal
cgarrow,-0.03,0.04,1.,0.04,linestyle=0,color='black',thick=3,/normal
cgtext,-0.08,0.9,textoidl('CTP'),/normal,charsize=charsize*2
cgtext,0.95,-0.01,textoidl('\tau'),/normal,charsize=charsize*2

cgtext,0.5,0.9,'CHART',color='red',/normal,charsize=charsize*2
cgtext,0.6,0.9,'CLIMCAPS',color='blue',/normal,charsize=charsize*2

for ictp=0,nbin_ctp-1 do begin
    for itau=0,nbin_tau-1 do begin

    position=pos[*,itau+nbin_tau*ictp]
   
    if itau eq 0 then ytitle='CTP '+ctpbin[ictp] else ytitle=''
    if ictp eq nbin_ctp-1 then xtitle='Q (g/kg)' else xtitle=''
    if ictp eq 0 then title=textoidl('\tau')+taulabel[itau] else title=''

    cgplot,EC_q_ISCCP[*,itau,ictp],pres_h2o,color='black',thick=thick,yrange=[1100,100],title=title,$
     charsize=charsize,xrange=xrange,position=position,/noerase,background=background,xtitle=xtitle,ytitle=ytitle,/nodata
    
    cgplot,(j_q_ISCCP[*,itau,ictp]-EC_q_ISCCP[*,itau,ictp])/EC_q_ISCCP[*,itau,ictp]*100,pres_h2o,color='red',thick=thick,/overplot
    cgplot,(EC_qj_ISCCP[*,itau,ictp]-EC_q_ISCCP[*,itau,ictp])/EC_q_ISCCP[*,itau,ictp]*100,pres_h2o,color='red',thick=thick,/overplot,linestyle=2
    cgplot,(c_q_ISCCP[*,itau,ictp]-EC_q_ISCCP[*,itau,ictp])/EC_q_ISCCP[*,itau,ictp]*100,pres_h2o,color='blue',thick=thick,/overplot
    cgplot,(EC_qc_ISCCP[*,itau,ictp]-EC_q_ISCCP[*,itau,ictp])/EC_q_ISCCP[*,itau,ictp]*100,pres_h2o,color='blue',thick=thick,/overplot,linestyle=2
    
    endfor

endfor
cgText, 0.5, 0.95, /Normal,'MMR (%): (CrIMSS-ECMWF)/ECMWF',Charsize=cgDefCharsize()*1.25, Alignment=0.5

erase

xrange=[-60,60]
cgarrow,-0.03,0.04,-0.03,0.95,linestyle=0,color='black',thick=3,/normal
cgarrow,-0.03,0.04,1.,0.04,linestyle=0,color='black',thick=3,/normal
cgtext,-0.08,0.9,textoidl('CTP'),/normal,charsize=charsize*2
cgtext,0.95,-0.01,textoidl('\tau'),/normal,charsize=charsize*2

cgtext,0.5,0.9,'CHART',color='red',/normal,charsize=charsize*2
cgtext,0.6,0.9,'CLIMCAPS',color='blue',/normal,charsize=charsize*2

for ictp=0,nbin_ctp-1 do begin
    for itau=0,nbin_tau-1 do begin

    position=pos[*,itau+nbin_tau*ictp]
   
    if itau eq 0 then ytitle='CTP '+ctpbin[ictp] else ytitle=''
    if ictp eq nbin_ctp-1 then xtitle='Q (g/kg)' else xtitle=''
    if ictp eq 0 then title=textoidl('\tau')+taulabel[itau] else title=''

    cgplot,EC_q_ISCCP[*,itau,ictp],pres_h2o,color='black',thick=thick,yrange=[1100,100],title=title,$
     charsize=charsize,xrange=xrange,position=position,/noerase,background=background,xtitle=xtitle,ytitle=ytitle,/nodata
    
    cgplot,(j_q_ISCCP[*,itau,ictp]-EC_qj_ISCCP[*,itau,ictp])/EC_qj_ISCCP[*,itau,ictp]*100,pres_h2o,color='red',thick=thick,/overplot
    cgplot,(c_q_ISCCP[*,itau,ictp]-EC_qc_ISCCP[*,itau,ictp])/EC_qc_ISCCP[*,itau,ictp]*100,pres_h2o,color='blue',thick=thick,/overplot
    
    endfor

endfor
cgText, 0.5, 0.95, /Normal,'MMR (%): (CrIMSS-ECMWF(QC))/ECMWF(QC)',Charsize=cgDefCharsize()*1.25, Alignment=0.5

cgps_close
cgps2pdf,psname+'.ps',psname+'.pdf'
spawn,'rm '+psname+'.ps'


stop

end

               
           
                
