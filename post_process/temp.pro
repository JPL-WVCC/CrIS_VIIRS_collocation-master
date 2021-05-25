pro test

rootc='/peate_archive/NPPOps/sndr/production/v01_23_01/'
variable_list=['obs_time_tai93','lon','lat','air_temp','air_temp_qc','cld_frac','cld_frac_qc','cld_top_pres','cld_top_pres_qc']
ngran=10

positions = cgLayout([2,2], OYMargin=[10,12], YGap=4, OXMargin=[3,3],xgap=2)
for ifile=0,22 do begin
    filename='VIIRS_Cloud_onCrIS_20150601_F'+string(ifile,format='(i3.3)')+'.sav'
    restore,filename
    CF=reform((VIIRS_CLD_Mask[0,*,*,*]*1.+VIIRS_CLD_Mask[1,*,*,*]*0.5+VIIRS_CLD_Mask[2,*,*,*]*0.25)/total(VIIRS_CLD_Mask,1,/nan))
    result=size(CF,/dimensions)
    hist=total(VIIRS_ISCCP_Hist,3,/nan)/transpose(rebin(total(total(total(VIIRS_ISCCP_Hist,3,/nan),1,/nan),1,/nan),30,450,7,7,/sample),[2,3,0,1])*100.
    for igran=0,ngran-1 do begin
        gran=string(ngran*ifile+igran+1,format='(i3.3)')
        spawn,'ls '+rootc+'2015/06/01/l2climcaps/SNDR*g'+gran+'.L2_RET*.nc',filec
        ierr=read_ncdf(filec[0],datac,variable_list=variable_list)
        find=where(datac.cld_frac_qc gt 1,nfind) 
        if nfind GE 1 then datac.cld_frac[find]=!values.f_nan
        find=where(datac.cld_top_pres_qc gt 1,nfind) 
        if nfind GE 1 then datac.cld_top_pres[find]=!values.f_nan
        if igran eq 0 then begin
           temp1= datac.cld_frac
           temp2= datac.cld_top_pres
        endif else begin
           temp1= concat(temp1,datac.cld_frac,dimension=3)
           temp2= concat(temp2,datac.cld_top_pres,dimension=3)
        endelse
    endfor
    ccf=(total(temp1,/nan,1))
    find=where(ccf GE 1,nfind) & if nfind GE 1 then ccf[find]=1

    cctp=total(temp1*temp2,1,/nan)/total(temp1,/nan,1)/100
    find=where(ccf ne ccf,nfind) 
    if nfind GE 1 then begin
       temp11=temp1[0,*,*,*] & temp22=temp2[0,*,*,*]
       cctp[find]=temp22[find] & ccf[find]=temp11[find]
    endif
   
    lon=reform(cris_lon[4,*,*]) & lat=reform(cris_lat[4,*,*])
    CF=mean(CF[*,*,*],/nan,dimension=1)
    ctp=mean(reform(VIIRS_CTP_moment[0,*,*,*]),/nan,dimension=1)
    ctpstd=mean(reform(VIIRS_CTP_moment[1,*,*,*]),/nan,dimension=1)

    psname=filename
    cgps_open,filename+'.ps'    
    p = positions[*,0]
    cgmap_set,0,180,/continents,position=p
    cgloadct,13,/brewer,/reverse
    cgcontour,CF,lon,lat,levels=indgen(10)*0.1,/fill,/overplot
    yspace = 4.0*(!D.Y_CH_SIZE)/!D.Y_SIZE    
    p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.02]
    cgColorbar, divisions=10, Range=[0,1], Position=p, $
       Charsize=cgDefCharsize()*0.65, title='VIIRS Mean CF on CrIS'

    p = positions[*,1]
    cgmap_set,0,180,/continents,position=p,/noerase
    cgloadct,33
    cgcontour,CTP,lon,lat,levels=indgen(10)*100+100,/fill,/overplot
    p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.02]
    cgColorbar, divisions=10, Range=[0,1000], Position=p, $
       Charsize=cgDefCharsize()*0.65, title='VIIRS Mean CTP on CrIS'

    yspace=5.0*(!D.Y_CH_SIZE)/!D.Y_SIZE
    p = positions[*,2]
    cgmap_set,0,180,/continents,position=p,/noerase,title='CrIS CLIMCAPS CF',Charsize=cgDefCharsize()*0.65
    cgloadct,13,/brewer,/reverse
    cgcontour,reform(cCF[4,*,*]),lon,lat,levels=indgen(10)*0.1,/fill,/overplot
    p = [0.1, p[1]-yspace-0.02, 0.425, p[1]-yspace]
    cgColorbar, divisions=10, Range=[0,1], Position=p, $
       Charsize=cgDefCharsize()*0.65

    p = positions[*,3]
    cgmap_set,0,180,/continents,position=p,/noerase,title='CrIS CLIMCAPS CTP',Charsize=cgDefCharsize()*0.65
    cgloadct,33
    cgcontour,reform(cCTP[4,*,*]),lon,lat,levels=indgen(10)*100.+100,/fill,/overplot
    p = [0.585, p[1]-yspace-0.02, 0.915, p[1]-yspace]
    cgColorbar, divisions=10, Range=[0,1000], Position=p, $
       Charsize=cgDefCharsize()*0.65
   
    cgps_close
    cgps2pdf,psname+'.ps',psname+'.pdf'
    spawn,'rm '+psname+'.ps'    
            
endfor
stop

end
