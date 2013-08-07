
;---------------------------------------------------------------------------
;           GNSPI algorithm for FILLING THE SLC-OFF GAP OF ETM+ IMAGES
;                           Using TM input images
;        Can process whole ETM+ scene using block strategy
;            Developed by Xiaolin Zhu,email: zhu.381@osu.edu
;             Department of Geography,The Ohio State University                  
;                         update date:  3/17/2013
;                     Copyright belong to Xiaolin Zhu
; Please cite the reference:
; Zhu, X., Liu, D. and Chen, J. 2012. A new geostatistical approach for filling
;      gaps in Landsat ETM+ SLC-off images, Remote Sensing of Environment,
;      124,49-60. 
;---------------------------------------------------------------------------


;function for open the file

Pro GetData,ImgData = ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,$
    FileName = FileName,Map_info = map_Info, Fid = Fid
    Filter = ['all file;*.*']
    Envi_Open_File,FileName,R_Fid = Fid
    Envi_File_Query,Fid,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type
    map_info = envi_get_map_info(fid=Fid)
    dims = [-1,0,ns - 1 ,0,nl - 1]
    case Data_Type Of
        1:ImgData = BytArr(ns,nl,nb)    ;  BYTE  Byte
        2:ImgData = IntArr(ns,nl,nb)    ;  INT  Integer
        3:ImgData = LonArr(ns,nl,nb)    ;  LONG  Longword integer
        4:ImgData = FltArr(ns,nl,nb)    ;  FLOAT  Floating point
        5:ImgData = DblArr(ns,nl,nb)    ;  DOUBLE  Double-precision floating
        6:ImgData = COMPLEXARR(ns,nl,nb); complex, single-precision, floating-point
        9:ImgData = DCOMPLEXARR(ns,nl,nb);complex, double-precision, floating-point
        12:ImgData = UINTARR(ns,nl,nb)   ; unsigned integer vector or array
        13:ImgData = ULONARR(ns,nl,nb)   ;  unsigned longword integer vector or array
        14:ImgData = LON64ARR(ns,nl,nb)   ;a 64-bit integer vector or array
        15:ImgData = ULON64ARR(ns,nl,nb)   ;an unsigned 64-bit integer vector or array
    EndCase
    For i = 0,nb-1 Do Begin
       Dt = Envi_Get_Data(Fid = Fid,dims = dims,pos=i)
       ImgData[*,*,i] = Dt[*,*]
    EndFor
End

; estimate semi-variogram function
PRO vgram2, data, locs, bintype, binparms, vargram, vardist, count
sz = size(locs)
if sz(0) eq 1 then d = 1 else d = sz(2)
N = n_elements(data)

;construct bin boundaries
if bintype eq 1 then begin
nbins = binparms(0) & binsize = binparms(1) & minbin = binparms(2)
bounds = binsize*findgen(nbins+1) + minbin
endif else if bintype eq 2 then begin
nbins = n_elements(binparms)-1 & bounds = binparms
endif

vargram = fltarr(nbins)
vardist = fltarr(nbins)
count = fltarr(nbins)
maxdist = 0.
for i=1,N-1 do begin
for j=0,i-1 do begin
if d eq 1 then pdist = abs(locs(i) - locs(j)) else pdist =norm(locs(*,i) - locs(*,j))
if pdist gt maxdist then maxdist = pdist
subs = where(bounds le pdist, cnt)
if cnt gt 0 and cnt le nbins then begin
sub = max(subs)
vardist(sub) = vardist(sub) + pdist
vargram(sub) = vargram(sub) + .5*(data(i)-data(j))^2
count(sub) = count(sub) + 1.
endif
endfor
endfor

for i=0,nbins-1 do begin
if count(i) gt 0 then begin
vardist(i) = vardist(i)/count(i)
vargram(i) = vargram(i)/count(i)
endif else vardist(i) = .5*(bounds(i) + bounds(i+1))
endfor

END

; exponential function of variogram

PRO gfunct, X, A, F, pder  ;fuction with nugget
  bx = EXP(A[1] * X)
  F = A[0] * bx + A[2]

  IF N_PARAMS() GE 4 THEN $
    pder = [[bx], [A[0] * X * bx], [replicate(1.0, N_ELEMENTS(X))]]
END

PRO gfunct2, X, A, F, pder   ;fuction without nugget
  bx = EXP(A[1] * X)
  F = A[0]*(1.0-bx)
    IF N_PARAMS() GE 3 THEN $
    pder = [[bx], [A[0] * X * bx]]
END

;-------------------------------------
; ordinary kriging function
;-------------------------------------
FUNCTION O_kriging,X,Y,Z,X0,Y0,sill,range,nugget

  N_E=N_ELEMENTS(X)

  D=FLTARR(N_E,N_E)
  for iok=0,N_E-1,1 do begin
    for jok=0,N_E-1,1 do begin
     D[iok,jok]=((X[iok]-X[jok])^2+(Y[iok]-Y[jok])^2)^0.5
    endfor
  endfor

  Ch=sill*exp(-3.0*D/range)
  ind_d0=where(D eq 0)
  Ch[ind_d0]=sill+nugget

  cs=fltarr(1,N_E)
  for iok=0,N_E-1,1 do begin
    Ds=((X[iok]-x0)^2+(Y[iok]-y0)^2)^0.5
    cs[0,iok]=sill*exp(-3.0*Ds/range)
  end

  Ch_plus=fltarr(N_E+1,N_E+1)+1
  Ch_plus[N_E,N_E]=0
  Ch_plus[0:N_E-1,0:N_E-1]=Ch
;  print,'ch_plus',Ch_plus
  cs_plus=fltarr(1,N_E+1)+1
  cs_plus[0,0:N_E-1]=cs
;   print,'cs_plus',cs_plus
  lamta=INVERT(Ch_plus)##cs_plus
;  print,'weight',lamta[0:N_E-1]
  result=fltarr(2)
  result[0]=total(lamta[0:N_E-1]*Z)
  result[1]=(sill-transpose(cs_plus)##INVERT(Ch_plus)##cs_plus)^0.5*1.96
  return, result
end





;##################################################################
;                  main program
;##################################################################

pro GNSPI

 t0=systime(1)                  ;the initial time of program running

 ;please set the following parameters
;----------------------------------------------------------------------
 sample_size=20                 ;set the sample size of sample pixels
 size_wind=12                   ;set the maximum window size
 class_num=4                    ;set the estimated number of classes
 num_series=1                   ;set the number of images in the time-series except the input image
 DN_min=0.0                       ;set the range of DN value of the image,If byte, 0 and 255
 DN_max=1.0
 patch_long=500                ;set the size of block,if process whole ETM scene, set 1000
 temp_file='G:\temp'            ;set the temporary file location
;------------------------------------------------------------------------


 ;open the SLC-off ETM+ image
  FileName1 = Dialog_PickFile(title = 'Open the target SLC-off ETM+ image:')
  envi_open_file,FileName1,r_fid=fid
  envi_file_query,fid,ns=ns,nl=nl,nb=nb,dims=dims
  map_info = envi_get_map_info(fid=fid)
  orig_ns=ns
  orig_nl=nl
  n_ns=ceil(float(ns)/patch_long)
  n_nl=ceil(float(nl)/patch_long)

  ind_patch1=intarr(4,n_ns*n_nl)           ;divide the whole scene into 1000*1000 block
  ind_patch=intarr(4,n_ns*n_nl)
  location=intarr(4,n_ns*n_nl)

  for i_ns=0,n_ns-1,1 do begin
    for i_nl=0,n_nl-1,1 do begin
        ind_patch1[0,n_ns*i_nl+i_ns]=i_ns*patch_long
        ind_patch[0,n_ns*i_nl+i_ns]=max([0,ind_patch1[0,n_ns*i_nl+i_ns]-size_wind])
        location[0,n_ns*i_nl+i_ns]=ind_patch1[0,n_ns*i_nl+i_ns]-ind_patch[0,n_ns*i_nl+i_ns]

        ind_patch1[1,n_ns*i_nl+i_ns]=min([ns-1,(i_ns+1)*patch_long-1])
        ind_patch[1,n_ns*i_nl+i_ns]=min([ns-1,ind_patch1[1,n_ns*i_nl+i_ns]+size_wind])
        location[1,n_ns*i_nl+i_ns]=ind_patch1[1,n_ns*i_nl+i_ns]-ind_patch1[0,n_ns*i_nl+i_ns]+location[0,n_ns*i_nl+i_ns]

        ind_patch1[2,n_ns*i_nl+i_ns]=i_nl*patch_long
        ind_patch[2,n_ns*i_nl+i_ns]=max([0,ind_patch1[2,n_ns*i_nl+i_ns]-size_wind])
        location[2,n_ns*i_nl+i_ns]=ind_patch1[2,n_ns*i_nl+i_ns]-ind_patch[2,n_ns*i_nl+i_ns]

        ind_patch1[3,n_ns*i_nl+i_ns]=min([nl-1,(i_nl+1)*patch_long-1])
        ind_patch[3,n_ns*i_nl+i_ns]=min([nl-1,ind_patch1[3,n_ns*i_nl+i_ns]+size_wind])
        location[3,n_ns*i_nl+i_ns]=ind_patch1[3,n_ns*i_nl+i_ns]-ind_patch1[2,n_ns*i_nl+i_ns]+location[2,n_ns*i_nl+i_ns]
    endfor
  endfor

  tempoutname=temp_file+'\temp_target'

  pos=indgen(nb)
  for isub=0,n_ns*n_nl-1,1 do begin
      dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
      envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
      out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
      envi_file_mng, id=r_fid1, /remove
  endfor

  envi_file_mng, id=fid, /remove

  ;-----------------------------------------------------------      ;input no-gap image and divide into block
  FileName2 = Dialog_PickFile(title = 'Open the input image:')
  envi_open_file,FileName2,r_fid=fid
  tempoutname=temp_file+'\temp_input'
  pos=indgen(nb)
  for isub=0,n_ns*n_nl-1,1 do begin
      dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
      envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
      out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
      envi_file_mng, id=r_fid1, /remove
  endfor
   envi_file_mng, id=fid, /remove
   
  ;-----------------------------------------------------------      ;input no-gap time-series and divide into block
  for i_timeseries=1,num_series,1 do begin
     FileName2 = Dialog_PickFile(title = 'Open the_'+strtrim(i_timeseries,1)+'_image of the time-series:')
     envi_open_file,FileName2,r_fid=fid
     tempoutname=temp_file+'\temp_series'+strtrim(i_timeseries,1)
     pos=indgen(nb)
     for isub=0,n_ns*n_nl-1,1 do begin
        dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
        envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
        out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
        envi_file_mng, id=r_fid1, /remove
     endfor
     envi_file_mng, id=fid, /remove
   endfor

;------------------------------------------------------------------
        ; begin process the gap for each block
;-------------------------------------------------------------------
tempoutname1=temp_file+'\temp_filled'
tempoutname2=temp_file+'\temp_uncertainty'

print,'there are total',n_ns*n_nl,' blocks'

for isub=0,n_ns*n_nl-1,1 do begin

       print,'Processing ',isub+1,' block'

    ;open each block image

    FileName = temp_file+'\temp_target'
    GetData,ImgData=ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,FileName = FileName+strtrim(isub+1,1),Fid = Fid1
    fine1=float(ImgData)
    fine0=fine1[location[0,isub]:location[1,isub],location[2,isub]:location[3,isub],*]    ;place the new image value
    error=fine0-fine0 ;output the prediction variance
    
    FileName = temp_file+'\temp_input'
    GetData,ImgData=ImgData,FileName = FileName+strtrim(isub+1,1),Fid = Fid2
    fine2=float(ImgData)
    ImgData=0 ;clear this variable
    
    image_series=fltarr(ns,nl,nb,num_series+1)
    image_series[*,*,*,0]=fine2
    for i_timeseries=1,num_series,1 do begin
          FileName = temp_file+'\temp_series'+strtrim(i_timeseries,1)
          GetData,ImgData=ImgData,FileName = FileName+strtrim(isub+1,1),Fid = Fid3
          image_series[*,*,*,i_timeseries]=float(ImgData)
          ImgData=0 ;clear this variable
          envi_file_mng, id=Fid3, /remove, /delete
    endfor
    
  ;classify the input image by k-means algorithm 
         fine2_1=reform(fine2,long(ns)*nl,nb)
         fine2_1=transpose(fine2_1)
         weights = CLUST_WTS(fine2_1, N_CLUSTERS =class_num)
         class=CLUSTER(fine2_1, weights, N_CLUSTERS =class_num)
         class=reform(class,ns,nl)
         fine2_1=0;clear it
         weights=0;clear it
  
     ;get the gap location
    gap=intarr(ns,nl)
    for i=0,ns-1,1 do begin
       for j=0,nl-1,1 do begin
          ind_gap=where(fine1[i,j,*] eq 0, c_gap_band)
          ind_on=where(fine2[i,j,*] eq 0, c_gap_on)   ;do not consider the background as gap
          if (c_gap_band gt 0 and c_gap_on eq 0) then begin
             gap[i,j]=1
          endif
       endfor
    endfor
    
  ; compute the residule
    non_gap=where(gap ne 1)
    class_on=class[non_gap]  
    dif=fltarr(ns,nl,nb)     
    regress_result=fltarr(nb,class_num,2)
    for iband=0,nb-1,1 do begin
       fine1_on=(fine1[*,*,iband])[non_gap]
       fine2_on=(fine2[*,*,iband])[non_gap]
       for iclass=0,class_num-1,1 do begin
          ind_class=where(class_on eq iclass,num_class)
          if (num_class ge 2) then begin
            x_variable=fine2_on[ind_class]
            y_variable=fine1_on[ind_class]
            result=REGRESS(x_variable,y_variable,SIGMA=sigma, CONST=const,ftest=fvalue)
            regress_result[iband,iclass,0]=const
            regress_result[iband,iclass,1]=result
            sig=1.0-f_pdf(fvalue,1,num_class-2)
          endif else begin
            regress_result[iband,iclass,0]=0
            regress_result[iband,iclass,1]=1
          endelse
;      PRINT, 'class',iclass,'const: ', const,'Coefficients: ', result[*]
;      PRINT, 'Standard errors: ', sigma
;      print,'p_value:', sig
        endfor
     endfor
    print,'OLS regression (matrix structure:band/class/constant and slop):'
    print,regress_result

       ;predict the value and get the residual
    for i=0,ns-1,1 do begin
     for j=0,nl-1,1 do begin
       if (fine1[i,j,0] ne 0) then begin
         classij=class[i,j]
         for iband=0,nb-1,1 do begin
           dif[i,j,iband]=fine1[i,j,iband]-(fine2[i,j,iband]*regress_result[iband,classij,1]+regress_result[iband,classij,0])
         endfor
       endif
     endfor
    endfor
    
    
    ; estimate the semi-variogram for each class
    sample_semiv=1000  ;1000 samples for estimating semi-variogram
    result_semiv=fltarr(3,nb,class_num)
    for classi=0,class_num-1,1 do begin
    
        location_s=fltarr(2,sample_semiv)    ; select the samples randomly
        variable=fltarr(nb,sample_semiv)
        ind_value=where(gap[*,*] ne 1 and class eq classi, c_value)
        if (c_value ge sample_semiv) then begin
          rand_num=randomu(seed,c_value)
          order=sort(rand_num)
          for i=0,sample_semiv-1,1 do begin
             loc_i=ind_value[order[i]]
             location_s[1,i]=fix(loc_i/ns)
             location_s[0,i]=loc_i-ns*location_s[1,i]
             for iband=0,nb-1,1 do begin
               variable[iband,i]=(dif[*,*,iband])[loc_i]
             endfor
          endfor 

         for iband=0,nb-1,1 do begin
             data1=variable[iband,*]
             type=1
             bynp=[10.0,4.0,0.0]
             vgram2, data1, location_s, type, bynp, vargram, vardist, count
             weights = count/(vargram^2)
             A1 = [-1.0*MAX(vargram),-3.0/(bynp[0]*bynp[1]),MAX(vargram)]
            ;Compute the parameters.
             yfit = CURVEFIT(vardist, vargram, weights, A1, SIGMA, FUNCTION_NAME='gfunct')
             psill=A1[2]
             nug=A1[2]+A1[0]
             range=-3.0/A1[1]
    
            if (nug ge 0) then begin
              result_semiv[0,iband,classi]=psill
              result_semiv[1,iband,classi]=nug
              result_semiv[2,iband,classi]=range
            endif else begin
               A2 = [MAX(vargram),-3.0/(bynp[0]*bynp[1])]
               yfit2 = CURVEFIT(vardist, vargram, weights, A2, SIGMA, FUNCTION_NAME='gfunct2')
               psill=A2[0]
               nug=0.0
               range=-3.0/A2[1]
               result_semiv[0,iband,classi]=psill
               result_semiv[1,iband,classi]=nug
               result_semiv[2,iband,classi]=range 
           endelse      
        endfor
       endif else begin
               result_semiv[0,*,classi]=1.0
               result_semiv[1,*,classi]=0.0
               result_semiv[2,*,classi]=40.0
       endelse
     endfor
  print,'semi-variogram parameters:sill,nugget,range(matrix structure:sill,nugget,range/band/class)'
  print,result_semiv
 
dis_window_whole=fltarr(2*size_wind+1,2*size_wind+1) ; spatial distance for a completed window
x_window_whole=intarr(2*size_wind+1,2*size_wind+1)
y_window_whole=intarr(2*size_wind+1,2*size_wind+1)
for i=0,2*size_wind,1 do begin
  for j=0,2*size_wind,1 do begin
     dis_window_whole[i,j]=((i-size_wind)^2+(j-size_wind)^2)^0.5
     x_window_whole[i,j]=i
     y_window_whole[i,j]=j
  endfor
endfor

 ;compute the threshold for selecting similar pixels
 similar_th=fltarr(num_series+1)
 for i_series=0,num_series,1 do begin
     similar_th_band=fltarr(nb)
     for iband=0,nb-1,1 do begin
        similar_th_band[iband]=stddev(image_series[*,*,iband,i_series])*2.0/float(class_num)    ;compute the threshold of similar pixel
     endfor
    similar_th[i_series]=mean(similar_th_band)
 endfor  
;print,'threshold of time-series:',similar_th
 
  
  
for i=location[0,isub],location[1,isub],1 do begin
  for j=location[2,isub],location[3,isub],1 do begin

    if (gap[i,j] eq 1) then begin  
       class_t=class[i,j]
       a1=max([0,i-size_wind])                      ;moving window
       a2=min([ns-1,i+size_wind])
       b1=max([0,j-size_wind])
       b2=min([nl-1,j+size_wind])
       sub_off=fine1[a1:a2,b1:b2,*]
       sub_dif=dif[a1:a2,b1:b2,*]
       gap_wind=gap[a1:a2,b1:b2]
       sub_class=class[a1:a2,b1:b2]
       sub_on=fine2[a1:a2,b1:b2,*]
       image_series_wind=image_series[a1:a2,b1:b2,*,*]

 
       ti=i-a1  ;location of target pixel
       tj=j-b1

       if (a2-a1 lt 2*size_wind or b2-b1 lt 2*size_wind ) then begin ; spatial distance for un-completed window
         dis_window=fltarr(a2-a1+1,b2-b1+1)
         x_window=intarr(a2-a1+1,b2-b1+1)
         y_window=intarr(a2-a1+1,b2-b1+1)
          for ii=0,a2-a1,1 do begin
            for jj=0,b2-b1,1 do begin
              dis_window[ii,jj]=((ii-ti)^2+(jj-tj)^2)^0.5
              x_window[ii,jj]=ii
              y_window[ii,jj]=jj
            endfor
          endfor
       endif else begin
          dis_window=dis_window_whole
          x_window=x_window_whole
          y_window=y_window_whole
       endelse
    
  
          
   ; find the intersection of similar pixels    
          position_intersect=intarr((a2-a1+1)*(b2-b1+1))+1  ;place the location of similar pixel in the intersection
          for i_series=0,num_series,1 do begin
              cand_band=intarr((a2-a1+1)*(b2-b1+1))      
              dis_window1=fltarr(a2-a1+1,b2-b1+1) ;spectral distance of each pixel in image1
             for ii=0,a2-a1,1 do begin
               for jj=0,b2-b1,1 do begin      
                 dis_window1[ii,jj]=(total((image_series_wind[ii,jj,*,i_series]-image_series_wind[ti,tj,*,i_series])^2)/3.0)^0.5
               endfor
             endfor
                ind_cand=where(dis_window1 lt similar_th[i_series])
                 cand_band[ind_cand]=1
                 position_intersect=position_intersect*cand_band
           endfor
          

          
     ;select the set of sample pixels     
      ind_common=where(gap_wind[*,*] ne 1 and sub_class eq class_t and position_intersect ne 0, c_gap)
      
      ; predict the gap-pixel values
       if (c_gap ge sample_size) then begin
           dis_window2=dis_window[ind_common]
           x_window2=x_window[ind_common]
           y_window2=y_window[ind_common]
           order_dis=sort(dis_window2)
           x_sample=x_window2[order_dis[0:sample_size-1]]
           y_sample=y_window2[order_dis[0:sample_size-1]]

            for iband=0,nb-1,1 do begin
                data_sample=((sub_dif[*,*,iband])[ind_common])[order_dis[0:sample_size-1]]
                sill=result_semiv[0,iband,class_t]
                range=result_semiv[2,iband,class_t]
                nugget=result_semiv[1,iband,class_t]
                result=O_kriging(x_sample,y_sample,data_sample,ti,tj,sill,range,nugget)
                fine0[i-location[0,isub],j-location[2,isub],iband]=(fine2[i,j,iband]*regress_result[iband,class_t,1]+regress_result[iband,class_t,0])+result[0]
                error[i-location[0,isub],j-location[2,isub],iband]=result[1]
            endfor
          endif else begin
              if (c_gap gt 0 and c_gap lt sample_size) then begin
                dis_window2=dis_window[ind_common]
                 x_sample=x_window[ind_common]
                 y_sample=y_window[ind_common]

                 for iband=0,nb-1,1 do begin
                    data_sample=(sub_dif[*,*,iband])[ind_common]
                    sill=result_semiv[0,iband,class_t]
                    range=result_semiv[2,iband,class_t]
                    nugget=result_semiv[1,iband,class_t]
                    result=O_kriging(x_sample,y_sample,data_sample,ti,tj,sill,range,nugget)
                    fine0[i-location[0,isub],j-location[2,isub],iband]=(fine2[i,j,iband]*regress_result[iband,class_t,1]+regress_result[iband,class_t,0])+result[0]
                    error[i-location[0,isub],j-location[2,isub],iband]=result[1]
                 endfor
               endif else begin              
                    ind_sameclass=where(gap_wind[*,*] ne 1 and sub_class eq class_t, c_sameclass)
                    if (c_sameclass gt 0) then begin                    
                       dis_window2=dis_window[ind_sameclass]
                       x_window2=x_window[ind_sameclass]
                       y_window2=y_window[ind_sameclass]
                       ind_sample=where(dis_window2 eq min(dis_window2))                   
                       x_sample=x_window2[ind_sample]
                       y_sample=y_window2[ind_sample]
                       for iband=0,nb-1,1 do begin
                         predict1=fine2[i,j,iband]*regress_result[iband,class_t,1]+regress_result[iband,class_t,0]
                         predict2=mean((sub_off[*,*,iband])[x_sample,y_sample])
                         fine0[i-location[0,isub],j-location[2,isub],iband]=(predict1+predict2)/2.0
                         error[i-location[0,isub],j-location[2,isub],iband]=0
                       endfor
                      endif else begin
                         for iband=0,nb-1,1 do begin  
                            fine0[i-location[0,isub],j-location[2,isub],iband]=fine2[i,j,iband]*regress_result[iband,class_t,1]+regress_result[iband,class_t,0] 
                            error[i-location[0,isub],j-location[2,isub],iband]=0
                          endfor
                      endelse
                  endelse
               endelse
                for iband=0,nb-1,1 do begin          ;modify the abnormal values
                    if ((fine0[i-location[0,isub],j-location[2,isub],iband]) lt DN_min or (fine0[i-location[0,isub],j-location[2,isub],iband]) gt DN_max)then begin
                        fine0[i-location[0,isub],j-location[2,isub],iband]=mean(image_series[i,j,iband,*])
                     endif
                endfor
        endif
     endfor
   endfor


uncertainty=error/fine0*100.0      ;compute the uncertainty of prediction



size_result=size(fine0)

    case Data_Type Of
        1:fine0 = Byte(fine0)    ;  BYTE  Byte
        2:fine0 = FIX(fine0)     ;  INT  Integer
        3:fine0 = LONG(fine0)    ;  LONG  Longword integer
        4:fine0 = FLOAT(fine0)   ;  FLOAT  Floating point
        5:fine0 = DOUBLE(fine0)  ;  DOUBLE  Double-precision floating
        6:fine0 = COMPLEX(fine0); complex, single-precision, floating-point
        9:fine0 = DCOMPLEX(fine0);complex, double-precision, floating-point
        12:fine0 = UINT(fine0)   ; unsigned integer vector or array
        13:fine0 = ULONG(fine0)   ;  unsigned longword integer vector or array
        14:fine0 = LONG64(fine0)   ;a 64-bit integer vector or array
        15:fine0 = ULONG64(fine0)   ;an unsigned 64-bit integer vector or array
    EndCase



         Envi_Write_Envi_File,fine0,Out_Name = tempoutname1+strtrim(isub+1,1)
         Envi_Write_Envi_File,uncertainty,Out_Name = tempoutname2+strtrim(isub+1,1)
         envi_file_mng, id=Fid1, /remove, /delete
         envi_file_mng, id=Fid2, /remove, /delete

endfor

;--------------------------------------------------------------------------------------
;mosiac all the filled patch

  mfid=intarr(n_ns*n_nl)
  mdims=intarr(5,n_ns*n_nl)
  mpos=intarr(nb,n_ns*n_nl)
  pos=indgen(nb)
  x0=intarr(n_ns*n_nl)
  y0=intarr(n_ns*n_nl)

  for isub=0,n_ns*n_nl-1,1 do begin
      envi_open_file, tempoutname1+strtrim(isub+1,1), r_fid= sub_fid
     if (sub_fid eq -1) then begin
       envi_batch_exit
       return
     endif
      envi_file_query,  sub_fid, ns=sub_ns, nl=sub_nl
      mfid[isub] = sub_fid
      mpos[*,isub] = indgen(nb)
      mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
      x0[isub]=ind_patch1[0,isub]
      y0[isub]=ind_patch1[2,isub]
  endfor

  xsize = orig_ns
  ysize = orig_nl
  pixel_size = [1.,1.]

  use_see_through = replicate(1L,n_ns*n_nl)
  see_through_val = replicate(0L,n_ns*n_nl)

    out_name=FileName1+'_filled_GNSPI'
    envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
    dims=mdims, out_name=out_name, xsize=xsize, $
    ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
    out_dt=Data_Type, pixel_size=pixel_size, $
    background=0, see_through_val=see_through_val, $
    use_see_through=use_see_through

    for i=0,n_ns*n_nl-1,1 do begin
      envi_file_mng, id=mfid[i], /remove, /delete
    endfor

;--------------------------------------------------------------------------------------
;mosiac all the uncertainty patch

  mfid=intarr(n_ns*n_nl)
  mdims=intarr(5,n_ns*n_nl)
  mpos=intarr(nb,n_ns*n_nl)
  pos=indgen(nb)
  x0=intarr(n_ns*n_nl)
  y0=intarr(n_ns*n_nl)

  for isub=0,n_ns*n_nl-1,1 do begin
      envi_open_file, tempoutname2+strtrim(isub+1,1), r_fid= sub_fid
     if (sub_fid eq -1) then begin
       envi_batch_exit
       return
     endif
      envi_file_query,  sub_fid, ns=sub_ns, nl=sub_nl
      mfid[isub] = sub_fid
      mpos[*,isub] = indgen(nb)
      mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
      x0[isub]=ind_patch1[0,isub]
      y0[isub]=ind_patch1[2,isub]
  endfor

  xsize = orig_ns
  ysize = orig_nl
  pixel_size = [1.,1.]

  use_see_through = replicate(1L,n_ns*n_nl)
  see_through_val = replicate(0L,n_ns*n_nl)

    out_name=FileName1+'_uncertainty_GNSPI'
    envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
    dims=mdims, out_name=out_name, xsize=xsize, $
    ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
    out_dt=Data_Type, pixel_size=pixel_size, $
    background=0, see_through_val=see_through_val, $
    use_see_through=use_see_through

    for i=0,n_ns*n_nl-1,1 do begin
      envi_file_mng, id=mfid[i], /remove, /delete
    endfor

print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'

END



