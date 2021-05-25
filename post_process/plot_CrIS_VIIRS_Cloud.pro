pro plot_CrIS_VIIRS_Cloud,year,month,day,hour,gran1,gran2

;;; datadir='/home/qyue/VIIRS/CrIS_VIIRS_collocation-master/post_process/VIIRS_Cloud_onCrIS/'
datadir='/home/leipan/tmp/VIIRS_Cloud_onCrIS/'

mode=string(year,format='(i4.4)')+string(month,format='(i2.2)')+string(day,format='(i2.2)')+'_F'+string(hour,format='(i3.3)')

restore,datadir+'VIIRS_Cloud_onCrIS_'+mode+'.sav'

lon=reform(cris_lon[*,*,(gran1-1)*45:(gran2)*45-1],9*long(30)*(gran2-gran1+1)*45)
lat=reform(cris_lat[*,*,(gran1-1)*45:(gran2)*45-1],9*long(30)*(gran2-gran1+1)*45)

CF0=reform((VIIRS_CLD_MASK[0,*,*,(gran1-1)*45:(gran2)*45-1])/total(VIIRS_CLD_MASK[*,*,*,(gran1-1)*45:(gran2)*45-1],1,/nan)*100.,9*long(30)*(gran2-gran1+1)*45)
CF1=reform((VIIRS_CLD_MASK[1,*,*,(gran1-1)*45:(gran2)*45-1])/total(VIIRS_CLD_MASK[*,*,*,(gran1-1)*45:(gran2)*45-1],1,/nan)*100.,9*long(30)*(gran2-gran1+1)*45)
CF2=reform((VIIRS_CLD_MASK[2,*,*,(gran1-1)*45:(gran2)*45-1])/total(VIIRS_CLD_MASK[*,*,*,(gran1-1)*45:(gran2)*45-1],1,/nan)*100.,9*long(30)*(gran2-gran1+1)*45)
CF3=reform((VIIRS_CLD_MASK[3,*,*,(gran1-1)*45:(gran2)*45-1])/total(VIIRS_CLD_MASK[*,*,*,(gran1-1)*45:(gran2)*45-1],1,/nan)*100.,9*long(30)*(gran2-gran1+1)*45)

COD0=reform((VIIRS_COD_Moment[0,*,*,(gran1-1)*45:(gran2)*45-1]),9*long(30)*(gran2-gran1+1)*45)
COD1=reform((VIIRS_COD_Moment[1,*,*,(gran1-1)*45:(gran2)*45-1]),9*long(30)*(gran2-gran1+1)*45)
COD2=reform((VIIRS_COD_Moment[2,*,*,(gran1-1)*45:(gran2)*45-1]),9*long(30)*(gran2-gran1+1)*45)
COD3=reform((VIIRS_COD_Moment[3,*,*,(gran1-1)*45:(gran2)*45-1]),9*long(30)*(gran2-gran1+1)*45)
CTP0=reform((VIIRS_CTP_Moment[0,*,*,(gran1-1)*45:(gran2)*45-1]),9*long(30)*(gran2-gran1+1)*45)
CTP1=reform((VIIRS_CTP_Moment[1,*,*,(gran1-1)*45:(gran2)*45-1]),9*long(30)*(gran2-gran1+1)*45)
CTP2=reform((VIIRS_CTP_Moment[2,*,*,(gran1-1)*45:(gran2)*45-1]),9*long(30)*(gran2-gran1+1)*45)
CTP3=reform((VIIRS_CTP_Moment[3,*,*,(gran1-1)*45:(gran2)*45-1]),9*long(30)*(gran2-gran1+1)*45)

psname=datadir+'VIIRS_Cloud_onCrIS_'+mode+'_Cloud'
cgps_open,psname+'.ps'

positions = cgLayout([2,2], OYMargin=[4,6], YGap=14, OXMargin=[4,2],xgap=6)
Charsize=1

cgloadct,33
 
p=positions[*,0]
zcolors = BytScl(CF0, Min=0, Max=150)
centerlon=(min(lon)+max(lon))*0.5 & centerlat=(min(lat)+max(lat))*0.5
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[0,100],divisions=10,tlocation='TOP',tickinterval=10,title='Cloudy CF',/Discrete,position=p,charsize=charsize*0.9

 
p=positions[*,1]
zcolors = BytScl(CF1, Min=0, Max=100)
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[0,100],divisions=10,tlocation='TOP',tickinterval=10,title='Uncertain CF',/Discrete,position=p,charsize=charsize*0.9

p=positions[*,2]
zcolors = BytScl(CF2, Min=0, Max=100)
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[0,100],divisions=10,tlocation='TOP',tickinterval=10,title='Possibly Clear Fraction',/Discrete,position=p,charsize=charsize*0.9

p=positions[*,3]
zcolors = BytScl(CF3, Min=0, Max=100)
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[0,100],divisions=10,tlocation='TOP',tickinterval=10,title='Confidently Clear Fraction',/Discrete,position=p,charsize=charsize*0.9
erase


p=positions[*,0]
zcolors = BytScl(COD0, Min=0, Max=100)
centerlon=(min(lon)+max(lon))*0.5 & centerlat=(min(lat)+max(lat))*0.5
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[0,100],divisions=10,tlocation='TOP',tickinterval=10,title='Mean VIIRS COD on CrIS',/Discrete,position=p,charsize=charsize*0.9

 
p=positions[*,1]
zcolors = BytScl(COD1, Min=0, Max=80)
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[0,80],divisions=10,tlocation='TOP',tickinterval=8,title='VIIRS COD Variance on CrIS',/Discrete,position=p,charsize=charsize*0.9

p=positions[*,2]
zcolors = BytScl(COD2, Min=-8, Max=8)
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[-8,8],divisions=10,tlocation='TOP',tickinterval=1.6,title='VIIRS COD Skewness on CrIS',/Discrete,position=p,charsize=charsize*0.9

p=positions[*,3]
zcolors = BytScl(COD3, Min=-3, Max=57)
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[-3,57],divisions=10,tlocation='TOP',tickinterval=6,title='VIIRS COD Kurtosis on CrIS',/Discrete,position=p,charsize=charsize*0.9
erase
p=positions[*,0]
zcolors = BytScl(CTP0, Min=200, Max=1000)
centerlon=(min(lon)+max(lon))*0.5 & centerlat=(min(lat)+max(lat))*0.5
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[200,1000],divisions=10,tlocation='TOP',tickinterval=80,title='Mean VIIRS CTP on CrIS',/Discrete,position=p,charsize=charsize*0.9

 
p=positions[*,1]
zcolors = BytScl(CTP1, Min=0, Max=500)
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[0,500],divisions=10,tlocation='TOP',tickinterval=50,title='VIIRS CTP Variance on CrIS',/Discrete,position=p,charsize=charsize*0.9

p=positions[*,2]
zcolors = BytScl(CTP2, Min=-8, Max=10)
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[-8,10],divisions=10,tlocation='TOP',tickinterval=1.8,title='VIIRS CTP Skewness on CrIS',/Discrete,position=p,charsize=charsize*0.9

p=positions[*,3]
zcolors = BytScl(CTP3, Min=-3, Max=97)
cgplot,lon,lat,position=p,/nodata,charsize=charsize,/noerase
cgmap_set,/cylindrical,centerlat,centerlon,limit=[min(lat),min(lon),max(lat),max(lon)],position=p,/noerase,/continents
cgplots,lon,lat,psym=16,color=zcolors,symsize=0.2
cgmap_grid,color='charcoal'
yspace = 1.28*(!D.Y_CH_SIZE)/!D.Y_SIZE
p = [p[0], p[3]+yspace, p[2], p[3]+yspace+0.01] & print,p
cgcolorbar,range=[-3,97],divisions=10,tlocation='TOP',tickinterval=10,title='VIIRS CTP Kurtosis on CrIS',/Discrete,position=p,charsize=charsize*0.9
erase

cgps_close

cgps2pdf,psname+'.ps',psname+'.pdf'
spawn,'rm '+psname+'.ps'

stop

end


