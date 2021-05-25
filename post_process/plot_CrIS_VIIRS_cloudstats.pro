pro plot_viirs_cloud

mod_taubin=[0.,0.3,1.3,3.6,9.4,23,60,200.]
mod_ctpbin=[50,180,310,440,560,680,800,1100]
ecfbin=[0,0.1,0.2,0.4,0.6,0.8,1.0]
nbin_tau=n_elements(mod_taubin)
nbin_ctp=n_elements(mod_ctpbin)
nbin_ecf=n_elements(ecfbin)
dir1='/raid15/qyue/VIIRS/VIIRS/201601/VIIRS_Cloud_onCrIS/'
dir2='/raid15/qyue/VIIRS/MODIS/201601/MODIS_Cloud_onAIRS/'
dir3=dir2

spawn,'ls '+dir1+'VIIRS_Cloud_onCrIS_2016011_G*.sav',myfiles1
spawn,'ls '+dir2+'MODIS_Cloud_onAIRS_2016011_G*.sav',myfiles2
spawn,'ls '+dir3+'MYD06_Cloud_onAIRS_2016011_G*.sav',myfiles3

ngran=n_elements(myfiles1)-1  ;-1 because of issues in modis airs collocation

sc1=13 & sc2=16 ;zero-based sc line number (0~29, sc2>=sc1)
sc11=sc1*3 & sc22=(sc2+1)*3-1  ;convert to 0~89
for i=0,ngran-1 do begin
    print,i
    myfile1=myfiles1[i] & myfile2=myfiles2[i] & myfile3=myfiles3[i]
    restore,myfile1 & restore,myfile2 & restore,myfile3
    if i eq 0 then begin
       airs_long=airs_lon[sc11:sc22,*] & airs_lati=airs_lat[sc11:sc22,*] 
       cris_long=cris_lon[*,sc1:sc2,*] & cris_lati=cris_lat[*,sc1:sc2,*]
       cmodis_isccp_hist=modis_isccp_hist[*,*,sc11:sc22,*]
       cmyd06_isccp_hist=myd06_isccp_hist[*,*,sc11:sc22,*] 
       cviirs_isccp_hist=viirs_isccp_hist[*,*,*,sc1:sc2,*]
       cmodis_cldtop_phase=modis_cldtop_phase[*,sc11:sc22,*]
       cviirs_cldtop_phase=viirs_cldtop_phase[*,*,sc1:sc2,*]
       cmodis_opt_phase=modis_opt_phase[*,sc11:sc22,*]
       cmyd06_opt_phase=myd06_opt_phase[*,sc11:sc22,*]
       cviirs_opt_phase=viirs_opt_phase[*,*,sc1:sc2,*]
       cmodis_cld_mask=modis_cld_mask[*,sc11:sc22,*]
       cmyd06_cld_mask=myd06_cld_mask[*,sc11:sc22,*]
       cviirs_cld_mask=viirs_cld_mask[*,*,sc1:sc2,*]
       
       cmodis_sfc_mask=modis_sfc_mask[*,sc11:sc22,*]
       cmyd06_sfc_mask=myd06_sfc_mask[*,sc11:sc22,*]
       cviirs_sfc_mask=viirs_sfc_mask[*,*,sc1:sc2,*]

       cmodis_ctp_moment=modis_ctp_moment[*,sc11:sc22,*]
       cmyd06_ctp_moment=myd06_ctp_moment[*,sc11:sc22,*]
       cviirs_ctp_moment=viirs_ctp_moment[*,*,sc1:sc2,*]
       cmodis_cod_moment=modis_cod_moment[*,sc11:sc22,*]
       cmyd06_cod_moment=myd06_cod_moment[*,sc11:sc22,*]
       cviirs_cod_moment=viirs_cod_moment[*,*,sc1:sc2,*]
       cmodis_re_moment=modis_re_moment[*,sc11:sc22,*]
       cmyd06_re_moment=myd06_re_moment[*,sc11:sc22,*]
       cviirs_re_moment=viirs_re_moment[*,*,sc1:sc2,*] 
    endif
    if i GE 1 then begin
       airs_long=concat(airs_long,airs_lon[sc11:sc22,*],dimension=1) & airs_lati=concat(airs_lati,airs_lat[sc11:sc22,*],dimension=1)
       cris_long=concat(cris_long,cris_lon[*,sc1:sc2,*],dimension=2) & cris_lati=concat(cris_lati,cris_lat[*,sc1:sc2,*],dimension=2) 

       cmodis_isccp_hist=concat(cmodis_isccp_hist,modis_isccp_hist[*,*,sc11:sc22,*],dimension=3) 
       cmyd06_isccp_hist=concat(cmyd06_isccp_hist,myd06_isccp_hist[*,*,sc11:sc22,*],dimension=3)
       cviirs_isccp_hist=concat(cviirs_isccp_hist,viirs_isccp_hist[*,*,*,sc1:sc2,*],dimension=4)
       cmodis_cldtop_phase=concat(cmodis_cldtop_phase,modis_cldtop_phase[*,sc11:sc22,*],dimension=2) 
       cviirs_cldtop_phase=concat(cviirs_cldtop_phase,viirs_cldtop_phase[*,*,sc1:sc2,*],dimension=3) 
       cmodis_opt_phase=concat(cmodis_opt_phase,modis_opt_phase[*,sc11:sc22,*],dimension=2) 
       cmyd06_opt_phase=concat(cmyd06_opt_phase,myd06_opt_phase[*,sc11:sc22,*],dimension=2) 
       cviirs_opt_phase=concat(cviirs_opt_phase,viirs_opt_phase[*,*,sc1:sc2,*],dimension=3) 
       cmodis_cld_mask=concat(cmodis_cld_mask,modis_cld_mask[*,sc11:sc22,*],dimension=2) 
       cmyd06_cld_mask=concat(cmyd06_cld_mask,myd06_cld_mask[*,sc11:sc22,*],dimension=2) 
       cviirs_cld_mask=concat(cviirs_cld_mask,viirs_cld_mask[*,*,sc1:sc2,*],dimension=3) 
       
       cmodis_sfc_mask=concat(cmodis_sfc_mask,modis_sfc_mask[*,sc11:sc22,*],dimension=2) 
       cmyd06_sfc_mask=concat(cmyd06_sfc_mask,myd06_sfc_mask[*,sc11:sc22,*],dimension=2) 
       cviirs_sfc_mask=concat(cviirs_sfc_mask,viirs_sfc_mask[*,*,sc1:sc2,*],dimension=3)

       cmodis_ctp_moment=concat(cmodis_ctp_moment,modis_ctp_moment[*,sc11:sc22,*],dimension=2) 
       cmyd06_ctp_moment=concat(cmyd06_ctp_moment,myd06_ctp_moment[*,sc11:sc22,*],dimension=2) 
       cviirs_ctp_moment=concat(cviirs_ctp_moment,viirs_ctp_moment[*,*,sc1:sc2,*],dimension=3) 
       cmodis_cod_moment=concat(cmodis_cod_moment,modis_cod_moment[*,sc11:sc22,*],dimension=2) 
       cmyd06_cod_moment=concat(cmyd06_cod_moment,myd06_cod_moment[*,sc11:sc22,*],dimension=2) 
       cviirs_cod_moment=concat(cviirs_cod_moment,viirs_cod_moment[*,*,sc1:sc2,*],dimension=3) 
       cmodis_re_moment=concat(cmodis_re_moment,modis_re_moment[*,sc11:sc22,*],dimension=2) 
       cmyd06_re_moment=concat(cmyd06_re_moment,myd06_re_moment[*,sc11:sc22,*],dimension=2) 
       cviirs_re_moment=concat(cviirs_re_moment,viirs_re_moment[*,*,sc1:sc2,*],dimension=3) 
       
    endif
   ; print,'myd06',n_elements(where(myd06_cld_mask[0,*,*] ne myd06_cld_mask[0,*,*])),$
   ;       'viirs',n_elements(where(viirs_cld_mask[0,*,*,*] ne viirs_cld_mask[0,*,*,*])),$
   ;       'MODIS',n_elements(where(modis_cld_mask[0,*,*] ne modis_cld_mask[0,*,*]))
endfor


npixel1=n_elements(airs_long) & npixel2=n_elements(cris_long)

airs_long=reform(airs_long,npixel1) & airs_lati=reform(airs_lati,npixel1)
cris_long=reform(cris_long,npixel2) & cris_lati=reform(cris_lati,npixel2)


modiscnt=reform(cmodis_sfc_mask[2,*,*]+cmodis_sfc_mask[3,*,*])
myd06cnt=reform(cmyd06_sfc_mask[2,*,*]+cmyd06_sfc_mask[3,*,*])
viirscnt=reform(cviirs_sfc_mask[2,*,*,*]+cviirs_sfc_mask[3,*,*,*])

cmodis_sfc_mask=reform(cmodis_sfc_mask/transpose(rebin(modiscnt,sc22-sc11+1,135l*ngran,10,/sample))*100,10,npixel1)
cmyd06_sfc_mask=reform(cmyd06_sfc_mask/transpose(rebin(myd06cnt,sc22-sc11+1,135l*ngran,10,/sample))*100,10,npixel1)
cviirs_sfc_mask=reform(cviirs_sfc_mask/transpose(rebin(viirscnt,9,(sc2-sc1+1),45l*ngran,10,/sample))*100,10,npixel2)

cmodis_isccp_hist=reform(cmodis_isccp_hist,nbin_tau-1,nbin_ctp-1,npixel1)
cmyd06_isccp_hist=reform(cmyd06_isccp_hist,nbin_tau-1,nbin_ctp-1,npixel1)
cviirs_isccp_hist=reform(cviirs_isccp_hist,nbin_tau-1,nbin_ctp-1,npixel2)

cmodis_cldtop_phase=reform(cmodis_cldtop_phase,5,npixel1)
cviirs_cldtop_phase=reform(cviirs_cldtop_phase,5,npixel2)
cmodis_opt_phase=reform(cmodis_opt_phase,5,npixel1)
cmyd06_opt_phase=reform(cmyd06_opt_phase,5,npixel1)
cviirs_opt_phase=reform(cviirs_opt_phase,5,npixel2)

cmodis_cld_mask=reform(cmodis_cld_mask,5,npixel1)
cmyd06_cld_mask=reform(cmyd06_cld_mask,5,npixel1)
cviirs_cld_mask=reform(cviirs_cld_mask,5,npixel2)

cmodis_ctp_mean=reform(cmodis_ctp_moment[0,*,*],npixel1)
cmyd06_ctp_mean=reform(cmyd06_ctp_moment[0,*,*],npixel1)
cviirs_ctp_mean=reform(cviirs_ctp_moment[0,*,*,*],npixel2)
cmodis_cod_mean=reform(cmodis_cod_moment[0,*,*],npixel1)
cmyd06_cod_mean=reform(cmyd06_cod_moment[0,*,*],npixel1)
cviirs_cod_mean=reform(cviirs_cod_moment[0,*,*,*],npixel2)

ncon=5
;airs_isccp_hist1=fltarr(nbin_ecf-1,nbin_ctp-1,ncon)
;airs_isccp_hist2=fltarr(nbin_ecf-1,nbin_ctp-1,ncon)
;cris_isccp_hist1=fltarr(nbin_ecf-1,nbin_ctp-1,ncon)
;cris_isccp_hist2=fltarr(nbin_ecf-1,nbin_ctp-1,ncon)
;airs_clr1=fltarr(ncon) & airs_clr2=airs_clr1 & cris_clr1=airs_clr1 & cris_clr2=airs_clr1
myd06_isccp=fltarr(nbin_tau-1,nbin_ctp-1,ncon)
modis_isccp=fltarr(nbin_tau-1,nbin_ctp-1,ncon)
viirs_isccp=fltarr(nbin_tau-1,nbin_ctp-1,ncon)
myd06_ophase=fltarr(5,ncon) & modis_ophase=myd06_ophase & viirs_ophase=myd06_ophase
myd06_cphase=fltarr(5,ncon) & modis_cphase=myd06_cphase & viirs_cphase=myd06_cphase
myd06_cmask=fltarr(5,ncon) & modis_cmask=myd06_cmask & viirs_cmask=myd06_cmask

for icon=0,ncon-1 do begin
    case icon of
    0: begin
        flag1=finite(airs_long) and cmodis_sfc_mask[1,*] GE 90
        flag2=finite(cris_long) and cviirs_sfc_mask[1,*] GE 90
       end
    1: begin
        flag1=cmodis_sfc_mask[9,*] GE 70 and cmodis_sfc_mask[1,*] GE 90
        flag2=cviirs_sfc_mask[9,*] GE 70  and cviirs_sfc_mask[1,*] GE 90;land
       end
    2: begin
        flag1=cmodis_sfc_mask[6,*] GE 90 and cmodis_sfc_mask[1,*] GE 90
        flag2=cviirs_sfc_mask[6,*] GE 90  and cviirs_sfc_mask[1,*] GE 90;ocean
       end
    3: begin
        flag1=cmodis_sfc_mask[4,*] GE 50 and cmodis_sfc_mask[1,*] GE 90
        flag2=cviirs_sfc_mask[4,*] GE 50 and cviirs_sfc_mask[1,*] GE 90;ice SFC
       end
    4: begin
        flag1=abs(airs_lati) LE 30 and cmodis_sfc_mask[1,*] GE 90 ;tropics
        flag2=abs(cris_lati) LE 30 and cviirs_sfc_mask[1,*] GE 90
       end
    endcase
;    for iecf=0,nbin_ecf-2 do begin
;        for ictp=0,nbin_ctp-2 do begin
;           find1=where(airs_ecf1 Gt ecfbin[iecf] and airs_ecf1 le ecfbin[iecf+1] and $
 ;                airs_ctp1 Gt mod_ctpbin[ictp] and airs_ctp1 le mod_ctpbin[ictp+1] and flag,nfind1)
;           airs_isccp_hist1[iecf,ictp,icon]=nfind1
;           find2=where(airs_ecf2 Gt ecfbin[iecf] and airs_ecf2 le ecfbin[iecf+1] and $
;                 airs_ctp2 Gt mod_ctpbin[ictp] and airs_ctp2 le mod_ctpbin[ictp+1] and flag,nfind2)
;           airs_isccp_hist2[iecf,ictp,icon]=nfind2
;           find1=where(cris_ecf1 Gt ecfbin[iecf] and cris_ecf1 le ecfbin[iecf+1] and $
;                 cris_ctp1 Gt mod_ctpbin[ictp] and cris_ctp1 le mod_ctpbin[ictp+1] and flag,nfind1)
;           cris_isccp_hist1[iecf,ictp,icon]=nfind1
;           find2=where(cris_ecf2 Gt ecfbin[iecf] and cris_ecf2 le ecfbin[iecf+1] and $
;                 cris_ctp2 Gt mod_ctpbin[ictp] and cris_ctp2 le mod_ctpbin[ictp+1] and flag,nfind2)
;           cris_isccp_hist2[iecf,ictp,icon]=nfind2
 ;       endfor
 ;   endfor
;    airs_isccp_hist1[*,*,icon]=airs_isccp_hist1[*,*,icon]/total(airs_isccp_hist1[*,*,icon],/nan)*100
;    airs_isccp_hist2[*,*,icon]=airs_isccp_hist2[*,*,icon]/total(airs_isccp_hist2[*,*,icon],/nan)*100
;    cris_isccp_hist1[*,*,icon]=cris_isccp_hist1[*,*,icon]/total(cris_isccp_hist1[*,*,icon],/nan)*100
;    cris_isccp_hist2[*,*,icon]=cris_isccp_hist2[*,*,icon]/total(cris_isccp_hist2[*,*,icon],/nan)*100

    find1=where(flag1) & find2=where(flag2)
    myd06_isccp[*,*,icon]=total(cmyd06_isccp_hist[*,*,find1],3,/nan)/total(cmyd06_isccp_hist[*,*,find1],/nan)*100.
    modis_isccp[*,*,icon]=total(cmodis_isccp_hist[*,*,find1],3,/nan)/total(cmodis_isccp_hist[*,*,find1],/nan)*100.
    viirs_isccp[*,*,icon]=total(cviirs_isccp_hist[*,*,find2],3,/nan)/total(cviirs_isccp_hist[*,*,find2],/nan)*100.
    myd06_ophase[*,icon]=total(cmyd06_opt_phase[*,find1],2,/nan)/total(cmyd06_opt_phase[*,find1],/nan)*100
    modis_ophase[*,icon]=total(cmodis_opt_phase[*,find1],2,/nan)/total(cmodis_opt_phase[*,find1],/nan)*100
    viirs_ophase[*,icon]=total(cviirs_opt_phase[*,find2],2,/nan)/total(cviirs_opt_phase[*,find2],/nan)*100
    
    modis_cphase[*,icon]=total(cmodis_cldtop_phase[*,find1],2,/nan)/total(cmodis_cldtop_phase[*,find1],/nan)*100
    viirs_cphase[*,icon]=total(cviirs_cldtop_phase[*,find2],2,/nan)/total(cviirs_cldtop_phase[*,find2],/nan)*100

    myd06_cmask[*,icon]=total(cmyd06_cld_mask[*,find1],2,/nan)/total(cmyd06_cld_mask[*,find1],/nan)*100
    modis_cmask[*,icon]=total(cmodis_cld_mask[*,find1],2,/nan)/total(cmodis_cld_mask[*,find1],/nan)*100
    viirs_cmask[*,icon]=total(cviirs_cld_mask[*,find2],2,/nan)/total(cviirs_cld_mask[*,find2],/nan)*100
  ;  airs_clr1[icon]=n_elements(where(airs_ecf1[find] le 0.01))/float(n_elements((airs_ecf1[find])))*100
  ;  airs_clr2[icon]=n_elements(where(airs_ecf2[find] le 0.01))/float(n_elements((airs_ecf2[find])))*100
  ;  cris_clr1[icon]=n_elements(where(cris_ecf1[find] le 0.01))/float(n_elements((cris_ecf1[find])))*100
  ;  cris_clr2[icon]=n_elements(where(cris_ecf2[find] le 0.01))/float(n_elements((cris_ecf2[find])))*100
endfor


psname=dir1+'Stats_Cloud_CrIS_AIRS_200601_Nadir'

cgps_open,psname+'.ps'
charsize= (!D.Name EQ 'PS') ? cgDefCharsize()*0.5 : cgDefCharsize()*0.75

mapoff=0
if mapoff eq 0 then begin
cgmap_set,/continents,/grid,/noerase
cgplot,airs_long,airs_lati,psym=3,color='red',/overplot
cgplot,airs_long[where(cmodis_sfc_mask[8,*] GE 50)],airs_lati[where(cmodis_sfc_mask[8,*] GE 50)],psym=3,color='black',/overplot
cgplot,airs_long[where(cmodis_sfc_mask[9,*] GE 70)],airs_lati[where(cmodis_sfc_mask[9,*] GE 70)],psym=3,color='green',/overplot
cgplot,airs_long[where(cmodis_sfc_mask[4,*] GE 50)],airs_lati[where(cmodis_sfc_mask[4,*] GE 50)],psym=3,color='blue',/overplot
cgmap_continents,/coasts,color='pink'
al_legend,['MODIS all data','desert','frozen sfc','land'],color=['red','black','blue','green'],linestyle=replicate(0,4),/left,box=0,charsize=charsize

erase
cgmap_set,/continents,/grid,/noerase
cgplot,airs_long,airs_lati,psym=3,color='red',/overplot
cgplot,airs_long[where(cmyd06_sfc_mask[8,*] GE 50)],airs_lati[where(cmyd06_sfc_mask[8,*] GE 50)],psym=3,color='black',/overplot
cgplot,airs_long[where(cmyd06_sfc_mask[9,*] GE 70)],airs_lati[where(cmyd06_sfc_mask[9,*] GE 70)],psym=3,color='green',/overplot
cgplot,airs_long[where(cmyd06_sfc_mask[4,*] GE 50)],airs_lati[where(cmyd06_sfc_mask[4,*] GE 50)],psym=3,color='blue',/overplot
cgmap_continents,/coasts,color='pink'

al_legend,['MYD06 all data','desert','frozen sfc','land'],color=['red','black','blue','green'],linestyle=replicate(0,4),/left,box=0,charsize=charsize

erase

cgmap_set,/continents,/grid,/noerase
cgplot,cris_long,cris_lati,psym=3,color='red',/overplot
cgplot,cris_long[where(cviirs_sfc_mask[8,*] GE 50)],cris_lati[where(cviirs_sfc_mask[8,*] GE 50)],psym=3,color='black',/overplot
cgplot,cris_long[where(cviirs_sfc_mask[9,*] GE 70)],cris_lati[where(cviirs_sfc_mask[9,*] GE 70)],psym=3,color='green',/overplot
cgplot,cris_long[where(cviirs_sfc_mask[4,*] GE 50)],cris_lati[where(cviirs_sfc_mask[4,*] GE 50)],psym=3,color='blue',/overplot
cgmap_continents,/coasts,color='pink'

al_legend,['VIIRS all data','desert','frozen sfc','land'],color=['red','black','blue','green'],linestyle=replicate(0,4),/left,box=0,charsize=charsize

erase
endif

pbar=[0.125, 0.98, 0.95, 0.99]
alg_titles=['MYD06','MODIS','VIIRS','AIRS V7','C-AIRS','C-SNPP-FSR','C-SNPP-NSR']
titles=['01/01/2016','Land','Ocean','Ice/Snow SFC','Tropics']
pos=cglayout([4,2],oymargin=[10,12],ygap=8,oxmargin=[2,2],xgap=6)
for icon=0,ncon-1 do begin
     for ipos=0,2 do begin
        case ipos of
        0: data=myd06_isccp[*,*,icon]
        1: data=modis_isccp[*,*,icon]
        2: data=viirs_isccp[*,*,icon]
       ; 3: data=airs_isccp_hist1[*,*,icon]
       ; 4: data=airs_isccp_hist2[*,*,icon]
       ; 5: data=cris_isccp_hist1[*,*,icon]
       ; 6: data=cris_isccp_hist2[*,*,icon]
        endcase
        if ipos Le 2 then begin
           xtitle='COD' & xarrays=mod_taubin
        endif else begin
           xtitle='ECF' & xarrays=ecfbin
        endelse
        yarrays=mod_ctpbin
        axis_format={xticks:n_elements(xarrays)-1,xTickname:strcompress(string(xarrays,format=('(f5.1)'))),yticks:nbin_ctp-1,ytickname:reverse(strcompress(string(mod_ctpbin)))}
        scaledata=cgscalevector(data,0,255,minvalue=0,maxvalue=15)   
        cgloadct,33,rgb_table=palette
        cgimage,reverse(scaledata,2),/axes,ytitle='CTP (hPa)',xtitle=xtitle,title=alg_titles[ipos],position=pos[*,ipos],charsize=charsize*0.6,$
               axkeywords=axis_format,/noerase,palette=palette
 ;       histcolors=bytscl(data,min=0,max=15,top=10)
 ;       cgloadct,33,ncolors=10
 ;       for itau=0,n_elements(xarrays)-2 do begin
 ;           for ictp=0,nbin_ctp-2 do begin
 ;               cgPolygon,[xarrays[itau],xarrays[itau],xarrays[itau+1],xarrays[itau+1],xarrays[itau]],$
  ;                        [yarrays[ictp],yarrays[ictp+1],yarrays[ictp+1],yarrays[ictp],yarrays[ictp]],$
  ;                        /normal,color=histcolors[itau,ictp]
  ;          endfor
  ;      endfor
    endfor
    cgloadct,33
    cgcolorbar,divisions=10,range=[0,15],title=titles[icon]+' Cld_Stats %',TLocation='top',position=pbar
    erase
endfor

;bar plot
!p.multi=0
colors=['black','red','blue','salmon','skyblue']
xtickv=[3,10,15,19,23,27,32,38,44] & xtickname=['No Mask','Phase Clr','Cld','U.Cld','U.Clr','Clr','Liquid','Ice','Mix/Uncert']
for icon=0,ncon-1 do begin
    cgplot,[0,48],[0,100],/nodata,xtitle='cloud phase types',ytitle='percentage',title=titles[icon]+' cloud phase stats',charsize=charsize,/noerase,xtickv=xtickv,xtickname=xtickname,xticks=n_elements(xtickv)-1
    ;plot no cld mask
    xarrays=indgen(7)+0.5
    for i=0,5 do begin
      case i of
      0:data=myd06_cmask[0,icon]
      1:data=modis_cmask[0,icon]
      2:data=viirs_cmask[0,icon]
      3:data=myd06_ophase[0,icon]
      4:data=modis_ophase[0,icon]
      5:data=viirs_ophase[0,icon]
      endcase
      if i le 2 then cgpolygon,[xarrays[i],xarrays[i],xarrays[i+1],xarrays[i+1],xarrays[i]],[0,data,data,0,0],/data,color=colors[i] $
      else cgpolygon,[xarrays[i],xarrays[i],xarrays[i+1],xarrays[i+1],xarrays[i]],[0,data,data,0,0],/data,color=colors[i-2],/fill
    endfor
    ;plot clr
    xarrays=indgen(6)+7.5
    for i=0,4 do begin
      case i of
      0:data=myd06_ophase[1,icon]
      1:data=modis_ophase[1,icon]
      2:data=viirs_ophase[1,icon]
      3:data=modis_cphase[0,icon]
      4:data=viirs_cphase[0,icon]
      endcase
      cgpolygon,[xarrays[i],xarrays[i],xarrays[i+1],xarrays[i+1],xarrays[i]],[0,data,data,0,0],/data,color=colors[i],/fill 
    endfor
 ;   cgplot,[xarrays[0]],[airs_clr1[icon]],psym=9,color='black',/overplot
 ;   cgplot,[xarrays[1]],[airs_clr2[icon]],psym=9,color='red',/overplot
 ;   cgplot,[xarrays[3]],[cris_clr1[icon]],psym=9,color='blue',/overplot
 ;   cgplot,[xarrays[4]],[cris_clr2[icon]],psym=9,color='magenta',/overplot
     ;plot cloud mask
    for j=0,3 do begin
      xarrays=indgen(4)+13.5+4*j
      for i=0,2 do begin
        case i of
         0:data=myd06_cmask[j+1,icon]
         1:data=modis_cmask[j+1,icon]
         2:data=viirs_cmask[j+1,icon]
        endcase
        cgpolygon,[xarrays[i],xarrays[i],xarrays[i+1],xarrays[i+1],xarrays[i]],[0,data,data,0,0],/data,color=colors[i] 
      endfor   
    endfor
    ;plot water
    xarrays=indgen(6)+29.5
    for i=0,4 do begin 
      case i of
      0:data=myd06_ophase[2,icon]
      1:data=modis_ophase[2,icon]
      2:data=viirs_ophase[2,icon]
      3:data=modis_cphase[1,icon]
      4:data=viirs_cphase[1,icon]
      endcase
      cgpolygon,[xarrays[i],xarrays[i],xarrays[i+1],xarrays[i+1],xarrays[i]],[0,data,data,0,0],/data,color=colors[i],/fill 
    endfor
    ;plot ice
    xarrays=indgen(6)+35.5
    for i=0,4 do begin 
      case i of
      0:data=myd06_ophase[3,icon]
      1:data=modis_ophase[3,icon]
      2:data=viirs_ophase[3,icon]
      3:data=modis_cphase[2,icon]
      4:data=viirs_cphase[2,icon]
      endcase
      cgpolygon,[xarrays[i],xarrays[i],xarrays[i+1],xarrays[i+1],xarrays[i]],[0,data,data,0,0],/data,color=colors[i],/fill 
    endfor
    ;plot uncertain and mix
    xarrays=indgen(6)+41.5
    cgaxis,yaxis=1,yrange=[0,4],/save,ytitle='Percentage for Mix and Uncert.',charsize=charsize
    for i=0,4 do begin 
      case i of
      0:data=myd06_ophase[4,icon]
      1:data=modis_ophase[4,icon]
      2:data=viirs_ophase[4,icon]
      3: begin
         data=modis_cphase[3,icon]+modis_cphase[4,icon] & data2=modis_cphase[3,icon]
         end
      4: begin
         data=viirs_cphase[3,icon]+viirs_cphase[4,icon] & data2=viirs_cphase[3,icon]
         end
      endcase
      cgpolygon,[xarrays[i],xarrays[i],xarrays[i+1],xarrays[i+1],xarrays[i]],[0,data,data,0,0],/data,color=colors[i],/fill 
      if i GE 3 then cgplot,[xarrays[i],xarrays[i+1]],[data2,data2],color='black',thick=8,/overplot
    endfor
    al_legend,['MYD06 Optical Phase','MODIS Optical Phase','VIIRS Optical Phase','MODIS CLDTOP Phase','VIIRS CLDTOP Phase'],$
      color=colors,/right,box=0,linestyle=replicate(0,5),thick=6,charsize=charsize
;    al_legend,['AIRS V7 Clr','CLIMCAPS-Aqua Clr','CLIMCAPS-SNPP FSR Clr','CLIMCAPS-SNPP NSR Clr'],color=['black','red','blue','magenta'],$
;     /left,box=0,psym=replicate(9,4),charsize=charsize

erase
endfor

cgps_close
cgps2pdf,psname+'.ps',psname+'.pdf'
spawn,'rm '+psname+'.ps'

stop
end
