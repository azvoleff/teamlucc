
;---------------------------------------------------------------------------
;           GNSPI algorithm for FILLING THE SLC-OFF GAP OF ETM+ IMAGES
;                           Using TM input images
;        Can process whole ETM+ scene using block strategy
;            Developed by Xiaolin Zhu, email: zhu.381@osu.edu
;             Department of Geography,The Ohio State University
;                         update date:  3/17/2013
;                     Copyright belong to Xiaolin Zhu
; Please cite the reference:
; Zhu, X., Liu, D. and Chen, J. 2012. A new geostatistical approach for filling
;      gaps in Landsat ETM+ SLC-off images, Remote Sensing of Environment,
;      124,49-60.
;-----------------------------------------------------------------------------
;
; Modified from original by Xiaolin Zhu on Feb 6, 2014 by Alex Zvoleff for 
; inclusion in the 'teamlucc' R package. For the original version, see Xiaolin's 
; website at: http://geography.osu.edu/grads/xzhu/
;
;-----------------------------------------------------------------------------


;function for open the file

PRO getdata,ImgData = ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,$
    FileName = FileName,Map_info = map_Info, Fid = Fid

    COMPILE_OPT idl2, hidden
    e = ENVI(/HEADLESS)
    ENVI_BATCH_STATUS_WINDOW, /ON

    Filter = ['all file;*.*']
    ENVI_OPEN_FILE,FileName,R_Fid = Fid
    ENVI_FILE_QUERY,Fid,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type
    map_info = ENVI_GET_MAP_INFO(fid=Fid)
    dims = [-1,0,ns - 1 ,0,nl - 1]
    CASE Data_Type OF
        1:ImgData = BYTARR(ns,nl,nb)    ;  BYTE  Byte
        2:ImgData = INTARR(ns,nl,nb)    ;  INT  Integer
        3:ImgData = LONARR(ns,nl,nb)    ;  LONG  Longword integer
        4:ImgData = FLTARR(ns,nl,nb)    ;  FLOAT  Floating point
        5:ImgData = DBLARR(ns,nl,nb)    ;  DOUBLE  Double-precision floating
        6:ImgData = COMPLEXARR(ns,nl,nb); complex, single-precision, floating-point
        9:ImgData = DCOMPLEXARR(ns,nl,nb);complex, double-precision, floating-point
        12:ImgData = UINTARR(ns,nl,nb)   ; unsigned integer vector or array
        13:ImgData = ULONARR(ns,nl,nb)   ;  unsigned longword integer vector or array
        14:ImgData = LON64ARR(ns,nl,nb)   ;a 64-bit integer vector or array
        15:ImgData = ULON64ARR(ns,nl,nb)   ;an unsigned 64-bit integer vector or array
    ENDCASE
    FOR i = 0,nb-1 DO BEGIN
        Dt = ENVI_GET_DATA(Fid = Fid,dims = dims,pos=i)
        ImgData[*,*,i] = Dt[*,*]
    ENDFOR
END

; estimate semi-variogram function
PRO vgram2, data, locs, bintype, binparms, vargram, vardist, count
    sz = SIZE(locs)
    IF sz(0) EQ 1 THEN d = 1 ELSE d = sz(2)
    N = N_ELEMENTS(data)

    ;construct bin boundaries
    IF bintype EQ 1 THEN BEGIN
        nbins = binparms(0) & binsize = binparms(1) & minbin = binparms(2)
        bounds = binsize*FINDGEN(nbins+1) + minbin
    ENDIF ELSE IF bintype EQ 2 THEN BEGIN
        nbins = N_ELEMENTS(binparms)-1 & bounds = binparms
    ENDIF

    vargram = FLTARR(nbins)
    vardist = FLTARR(nbins)
    count = FLTARR(nbins)
    maxdist = 0.
    FOR i=1,N-1 DO BEGIN
        FOR j=0,i-1 DO BEGIN
            IF d EQ 1 THEN pdist = ABS(locs(i) - locs(j)) ELSE pdist =norm(locs(*,i) - locs(*,j))
            IF pdist GT maxdist THEN maxdist = pdist
            subs = WHERE(bounds LE pdist, cnt)
            IF cnt GT 0 AND cnt LE nbins THEN BEGIN
                sub = MAX(subs)
                vardist(sub) = vardist(sub) + pdist
                vargram(sub) = vargram(sub) + .5*(data(i)-data(j))^2
                count(sub) = count(sub) + 1.
            ENDIF
        ENDFOR
    ENDFOR

    FOR i=0,nbins-1 DO BEGIN
        IF count(i) GT 0 THEN BEGIN
            vardist(i) = vardist(i)/count(i)
            vargram(i) = vargram(i)/count(i)
        ENDIF ELSE vardist(i) = .5*(bounds(i) + bounds(i+1))
    ENDFOR

END

; exponential function of variogram

PRO gfunct, X, A, F, pder  ;fuction with nugget
    bx = EXP(A[1] * X)
    F = A[0] * bx + A[2]

    IF N_PARAMS() GE 4 THEN $
        pder = [[bx], [A[0] * X * bx], [REPLICATE(1.0, N_ELEMENTS(X))]]
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
FUNCTION o_kriging,X,Y,Z,X0,Y0,sill,range,nugget

    N_E=N_ELEMENTS(X)

    D=FLTARR(N_E,N_E)
    FOR iok=0,N_E-1,1 DO BEGIN
        FOR jok=0,N_E-1,1 DO BEGIN
            D[iok,jok]=((X[iok]-X[jok])^2+(Y[iok]-Y[jok])^2)^0.5
        ENDFOR
    ENDFOR

    Ch=sill*EXP(-3.0*D/range)
    ind_d0=WHERE(D EQ 0)
    Ch[ind_d0]=sill+nugget

    cs=FLTARR(1,N_E)
    FOR iok=0,N_E-1,1 DO BEGIN
        Ds=((X[iok]-x0)^2+(Y[iok]-y0)^2)^0.5
        cs[0,iok]=sill*EXP(-3.0*Ds/range)
    END

    Ch_plus=FLTARR(N_E+1,N_E+1)+1
    Ch_plus[N_E,N_E]=0
    Ch_plus[0:N_E-1,0:N_E-1]=Ch
    ;  print,'ch_plus',Ch_plus
    cs_plus=FLTARR(1,N_E+1)+1
    cs_plus[0,0:N_E-1]=cs
    ;   print,'cs_plus',cs_plus
    lamta=INVERT(Ch_plus)##cs_plus
    ;  print,'weight',lamta[0:N_E-1]
    result=FLTARR(2)
    result[0]=TOTAL(lamta[0:N_E-1]*Z)
    result[1]=(sill-TRANSPOSE(cs_plus)##INVERT(Ch_plus)##cs_plus)^0.5*1.96
    RETURN, result
END

;##################################################################
;                  main program
;##################################################################

PRO gnspi, slc_off_file, input_file, timeseries_files, out_base, sample_size,$
    size_wind, class_num, DN_min, DN_max, patch_long, temp_dir

    COMPILE_OPT idl2, hidden

    e = ENVI(/HEADLESS)
    ENVI_BATCH_STATUS_WINDOW, /ON

    t0=SYSTIME(1)                  ;the initial time of program running

    ;open the SLC-off ETM+ image
    ENVI_OPEN_FILE,slc_off_file,r_fid=fid
    ENVI_FILE_QUERY,fid,ns=ns,nl=nl,nb=nb,dims=dims
    map_info = ENVI_GET_MAP_INFO(fid=fid)
    orig_ns=ns
    orig_nl=nl
    n_ns=CEIL(FLOAT(ns)/patch_long)
    n_nl=CEIL(FLOAT(nl)/patch_long)

    ind_patch1=INTARR(4,n_ns*n_nl)           ;divide the whole scene into 1000*1000 block
    ind_patch=INTARR(4,n_ns*n_nl)
    location=INTARR(4,n_ns*n_nl)

    FOR i_ns=0,n_ns-1,1 DO BEGIN
        FOR i_nl=0,n_nl-1,1 DO BEGIN
            ind_patch1[0,n_ns*i_nl+i_ns]=i_ns*patch_long
            ind_patch[0,n_ns*i_nl+i_ns]=MAX([0,ind_patch1[0,n_ns*i_nl+i_ns]-size_wind])
            location[0,n_ns*i_nl+i_ns]=ind_patch1[0,n_ns*i_nl+i_ns]-ind_patch[0,n_ns*i_nl+i_ns]

            ind_patch1[1,n_ns*i_nl+i_ns]=MIN([ns-1,(i_ns+1)*patch_long-1])
            ind_patch[1,n_ns*i_nl+i_ns]=MIN([ns-1,ind_patch1[1,n_ns*i_nl+i_ns]+size_wind])
            location[1,n_ns*i_nl+i_ns]=ind_patch1[1,n_ns*i_nl+i_ns]-ind_patch1[0,n_ns*i_nl+i_ns]+location[0,n_ns*i_nl+i_ns]

            ind_patch1[2,n_ns*i_nl+i_ns]=i_nl*patch_long
            ind_patch[2,n_ns*i_nl+i_ns]=MAX([0,ind_patch1[2,n_ns*i_nl+i_ns]-size_wind])
            location[2,n_ns*i_nl+i_ns]=ind_patch1[2,n_ns*i_nl+i_ns]-ind_patch[2,n_ns*i_nl+i_ns]

            ind_patch1[3,n_ns*i_nl+i_ns]=MIN([nl-1,(i_nl+1)*patch_long-1])
            ind_patch[3,n_ns*i_nl+i_ns]=MIN([nl-1,ind_patch1[3,n_ns*i_nl+i_ns]+size_wind])
            location[3,n_ns*i_nl+i_ns]=ind_patch1[3,n_ns*i_nl+i_ns]-ind_patch1[2,n_ns*i_nl+i_ns]+location[2,n_ns*i_nl+i_ns]
        ENDFOR
    ENDFOR

    tempoutname=temp_dir+'\temp_target'

    pos=INDGEN(nb)
    FOR isub=0,n_ns*n_nl-1,1 DO BEGIN
        dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
        ENVI_DOIT, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
            out_name=tempoutname+STRTRIM(isub+1,1), r_fid=r_fid1
        ENVI_FILE_MNG, id=r_fid1, /remove
    ENDFOR

    ENVI_FILE_MNG, id=fid, /remove

    ;-----------------------------------------------------------      ;input no-gap image and divide into block
    ENVI_OPEN_FILE,input_file,r_fid=fid
    tempoutname=temp_dir+'\temp_input'
    pos=INDGEN(nb)
    FOR isub=0,n_ns*n_nl-1,1 DO BEGIN
        dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
        ENVI_DOIT, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
            out_name=tempoutname+STRTRIM(isub+1,1), r_fid=r_fid1
        ENVI_FILE_MNG, id=r_fid1, /remove
    ENDFOR
    ENVI_FILE_MNG, id=fid, /remove

    ;-----------------------------------------------------------      ;input no-gap time-series and divide into block
    num_series=SIZE(timeseries_files, /N_ELEMENTS)
    FOR i_timeseries=1,num_series,1 DO BEGIN
        timeseries_file=timeseries_files[i_timeseries-1]
        ENVI_OPEN_FILE,timeseries_file,r_fid=fid
        tempoutname=temp_dir+'\temp_series'+STRTRIM(i_timeseries,1)
        pos=INDGEN(nb)
        FOR isub=0,n_ns*n_nl-1,1 DO BEGIN
            dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
            ENVI_DOIT, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
                out_name=tempoutname+STRTRIM(isub+1,1), r_fid=r_fid1
            ENVI_FILE_MNG, id=r_fid1, /remove
        ENDFOR
        ENVI_FILE_MNG, id=fid, /remove
    ENDFOR

    ;------------------------------------------------------------------
    ; begin process the gap for each block
    ;-------------------------------------------------------------------
    tempoutname1=temp_dir+'\temp_filled'
    tempoutname2=temp_dir+'\temp_uncertainty'

    PRINT,'there are total',n_ns*n_nl,' blocks'

    FOR isub=0,n_ns*n_nl-1,1 DO BEGIN

        PRINT,'Processing ',isub+1,' block'

        ;open each block image

        FileName = temp_dir+'\temp_target'
        getdata,ImgData=ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,FileName = FileName+STRTRIM(isub+1,1),Fid = Fid1
        fine1=FLOAT(ImgData)
        fine0=fine1[location[0,isub]:location[1,isub],location[2,isub]:location[3,isub],*]    ;place the new image value
        error=fine0-fine0 ;output the prediction variance

        FileName = temp_dir+'\temp_input'
        getdata,ImgData=ImgData,FileName = FileName+STRTRIM(isub+1,1),Fid = Fid2
        fine2=FLOAT(ImgData)
        ImgData=0 ;clear this variable

        image_series=FLTARR(ns,nl,nb,num_series+1)
        image_series[*,*,*,0]=fine2
        FOR i_timeseries=1,num_series,1 DO BEGIN
            FileName = temp_dir+'\temp_series'+STRTRIM(i_timeseries,1)
            getdata,ImgData=ImgData,FileName = FileName+STRTRIM(isub+1,1),Fid = Fid3
            image_series[*,*,*,i_timeseries]=FLOAT(ImgData)
            ImgData=0 ;clear this variable
            ENVI_FILE_MNG, id=Fid3, /remove, /delete
        ENDFOR

        ;classify the input image by k-means algorithm
        fine2_1=REFORM(fine2,LONG(ns)*nl,nb)
        fine2_1=TRANSPOSE(fine2_1)
        weights = clust_wts(fine2_1, N_CLUSTERS =class_num)
        class=cluster(fine2_1, weights, N_CLUSTERS =class_num)
        class=REFORM(class,ns,nl)
        fine2_1=0;clear it
        weights=0;clear it

        ;get the gap location
        gap=INTARR(ns,nl)
        FOR i=0,ns-1,1 DO BEGIN
            FOR j=0,nl-1,1 DO BEGIN
                ind_gap=WHERE(fine1[i,j,*] EQ 0, c_gap_band)
                ind_on=WHERE(fine2[i,j,*] EQ 0, c_gap_on)   ;do not consider the background as gap
                IF (c_gap_band GT 0 AND c_gap_on EQ 0) THEN BEGIN
                    gap[i,j]=1
                ENDIF
            ENDFOR
        ENDFOR

        ; compute the residule
        non_gap=WHERE(gap NE 1)
        class_on=class[non_gap]
        dif=FLTARR(ns,nl,nb)
        regress_result=FLTARR(nb,class_num,2)
        FOR iband=0,nb-1,1 DO BEGIN
            fine1_on=(fine1[*,*,iband])[non_gap]
            fine2_on=(fine2[*,*,iband])[non_gap]
            FOR iclass=0,class_num-1,1 DO BEGIN
                ind_class=WHERE(class_on EQ iclass,num_class)
                IF (num_class GE 2) THEN BEGIN
                    x_variable=fine2_on[ind_class]
                    y_variable=fine1_on[ind_class]
                    result=regress(x_variable,y_variable,SIGMA=sigma, CONST=const,ftest=fvalue)
                    regress_result[iband,iclass,0]=const
                    regress_result[iband,iclass,1]=result
                    sig=1.0-f_pdf(fvalue,1,num_class-2) ; get underflow here
                ENDIF ELSE BEGIN
                    regress_result[iband,iclass,0]=0
                    regress_result[iband,iclass,1]=1
                ENDELSE
                ;      PRINT, 'class',iclass,'const: ', const,'Coefficients: ', result[*]
                ;      PRINT, 'Standard errors: ', sigma
                ;      print,'p_value:', sig
            ENDFOR
        ENDFOR
        PRINT,'OLS regression (matrix structure:band/class/constant and slop):'
        PRINT,regress_result

        ;predict the value and get the residual
        FOR i=0,ns-1,1 DO BEGIN
            FOR j=0,nl-1,1 DO BEGIN
                IF (fine1[i,j,0] NE 0) THEN BEGIN
                    classij=class[i,j]
                    FOR iband=0,nb-1,1 DO BEGIN
                        dif[i,j,iband]=fine1[i,j,iband]-(fine2[i,j,iband]*regress_result[iband,classij,1]+regress_result[iband,classij,0])
                    ENDFOR
                ENDIF
            ENDFOR
        ENDFOR


        ; estimate the semi-variogram for each class
        sample_semiv=1000  ;1000 samples for estimating semi-variogram
        result_semiv=FLTARR(3,nb,class_num)
        FOR classi=0,class_num-1,1 DO BEGIN

            location_s=FLTARR(2,sample_semiv)    ; select the samples randomly
            variable=FLTARR(nb,sample_semiv)
            ind_value=WHERE(gap[*,*] NE 1 AND class EQ classi, c_value)
            IF (c_value GE sample_semiv) THEN BEGIN
                rand_num=RANDOMU(seed,c_value)
                order=SORT(rand_num)
                FOR i=0,sample_semiv-1,1 DO BEGIN
                    loc_i=ind_value[order[i]]
                    location_s[1,i]=FIX(loc_i/ns)
                    location_s[0,i]=loc_i-ns*location_s[1,i]
                    FOR iband=0,nb-1,1 DO BEGIN
                        variable[iband,i]=(dif[*,*,iband])[loc_i]
                    ENDFOR
                ENDFOR

                FOR iband=0,nb-1,1 DO BEGIN
                    data1=variable[iband,*]
                    type=1
                    bynp=[10.0,4.0,0.0]
                    vgram2, data1, location_s, type, bynp, vargram, vardist, count
                    weights = count/(vargram^2)
                    A1 = [-1.0*MAX(vargram),-3.0/(bynp[0]*bynp[1]),MAX(vargram)]
                    ;Compute the parameters.
                    yfit = curvefit(vardist, vargram, weights, A1, SIGMA, FUNCTION_NAME='gfunct')
                    psill=A1[2]
                    nug=A1[2]+A1[0]
                    range=-3.0/A1[1]

                    IF (nug GE 0) THEN BEGIN
                        result_semiv[0,iband,classi]=psill
                        result_semiv[1,iband,classi]=nug
                        result_semiv[2,iband,classi]=range
                    ENDIF ELSE BEGIN
                        A2 = [MAX(vargram),-3.0/(bynp[0]*bynp[1])]
                        yfit2 = curvefit(vardist, vargram, weights, A2, SIGMA, FUNCTION_NAME='gfunct2')
                        psill=A2[0]
                        nug=0.0
                        range=-3.0/A2[1]
                        result_semiv[0,iband,classi]=psill
                        result_semiv[1,iband,classi]=nug
                        result_semiv[2,iband,classi]=range
                    ENDELSE
                ENDFOR
            ENDIF ELSE BEGIN
                result_semiv[0,*,classi]=1.0
                result_semiv[1,*,classi]=0.0
                result_semiv[2,*,classi]=40.0
            ENDELSE
        ENDFOR
        PRINT,'semi-variogram parameters:sill,nugget,range(matrix structure:sill,nugget,range/band/class)'
        PRINT,result_semiv

        dis_window_whole=FLTARR(2*size_wind+1,2*size_wind+1) ; spatial distance for a completed window
        x_window_whole=INTARR(2*size_wind+1,2*size_wind+1)
        y_window_whole=INTARR(2*size_wind+1,2*size_wind+1)
        FOR i=0,2*size_wind,1 DO BEGIN
            FOR j=0,2*size_wind,1 DO BEGIN
                dis_window_whole[i,j]=((i-size_wind)^2+(j-size_wind)^2)^0.5
                x_window_whole[i,j]=i
                y_window_whole[i,j]=j
            ENDFOR
        ENDFOR

        ;compute the threshold for selecting similar pixels
        similar_th=FLTARR(num_series+1)
        FOR i_series=0,num_series,1 DO BEGIN
            similar_th_band=FLTARR(nb)
            FOR iband=0,nb-1,1 DO BEGIN
                similar_th_band[iband]=stddev(image_series[*,*,iband,i_series])*2.0/FLOAT(class_num)    ;compute the threshold of similar pixel
            ENDFOR
            similar_th[i_series]=mean(similar_th_band)
        ENDFOR
        ;print,'threshold of time-series:',similar_th



        FOR i=location[0,isub],location[1,isub],1 DO BEGIN
            FOR j=location[2,isub],location[3,isub],1 DO BEGIN

                IF (gap[i,j] EQ 1) THEN BEGIN
                    class_t=class[i,j]
                    a1=MAX([0,i-size_wind])                      ;moving window
                    a2=MIN([ns-1,i+size_wind])
                    b1=MAX([0,j-size_wind])
                    b2=MIN([nl-1,j+size_wind])
                    sub_off=fine1[a1:a2,b1:b2,*]
                    sub_dif=dif[a1:a2,b1:b2,*]
                    gap_wind=gap[a1:a2,b1:b2]
                    sub_class=class[a1:a2,b1:b2]
                    sub_on=fine2[a1:a2,b1:b2,*]
                    image_series_wind=image_series[a1:a2,b1:b2,*,*]


                    ti=i-a1  ;location of target pixel
                    tj=j-b1

                    IF (a2-a1 LT 2*size_wind OR b2-b1 LT 2*size_wind ) THEN BEGIN ; spatial distance for un-completed window
                        dis_window=FLTARR(a2-a1+1,b2-b1+1)
                        x_window=INTARR(a2-a1+1,b2-b1+1)
                        y_window=INTARR(a2-a1+1,b2-b1+1)
                        FOR ii=0,a2-a1,1 DO BEGIN
                            FOR jj=0,b2-b1,1 DO BEGIN
                                dis_window[ii,jj]=((ii-ti)^2+(jj-tj)^2)^0.5
                                x_window[ii,jj]=ii
                                y_window[ii,jj]=jj
                            ENDFOR
                        ENDFOR
                    ENDIF ELSE BEGIN
                        dis_window=dis_window_whole
                        x_window=x_window_whole
                        y_window=y_window_whole
                    ENDELSE



                    ; find the intersection of similar pixels
                    position_intersect=INTARR((a2-a1+1)*(b2-b1+1))+1  ;place the location of similar pixel in the intersection
                    FOR i_series=0,num_series,1 DO BEGIN
                        cand_band=INTARR((a2-a1+1)*(b2-b1+1))
                        dis_window1=FLTARR(a2-a1+1,b2-b1+1) ;spectral distance of each pixel in image1
                        FOR ii=0,a2-a1,1 DO BEGIN
                            FOR jj=0,b2-b1,1 DO BEGIN
                                dis_window1[ii,jj]=(TOTAL((image_series_wind[ii,jj,*,i_series]-image_series_wind[ti,tj,*,i_series])^2)/3.0)^0.5
                            ENDFOR
                        ENDFOR
                        ind_cand=WHERE(dis_window1 LT similar_th[i_series])
                        cand_band[ind_cand]=1
                        position_intersect=position_intersect*cand_band
                    ENDFOR



                    ;select the set of sample pixels
                    ind_common=WHERE(gap_wind[*,*] NE 1 AND sub_class EQ class_t AND position_intersect NE 0, c_gap)

                    ; predict the gap-pixel values
                    IF (c_gap GE sample_size) THEN BEGIN
                        dis_window2=dis_window[ind_common]
                        x_window2=x_window[ind_common]
                        y_window2=y_window[ind_common]
                        order_dis=SORT(dis_window2)
                        x_sample=x_window2[order_dis[0:sample_size-1]]
                        y_sample=y_window2[order_dis[0:sample_size-1]]

                        FOR iband=0,nb-1,1 DO BEGIN
                            data_sample=((sub_dif[*,*,iband])[ind_common])[order_dis[0:sample_size-1]]
                            sill=result_semiv[0,iband,class_t]
                            range=result_semiv[2,iband,class_t]
                            nugget=result_semiv[1,iband,class_t]
                            result=o_kriging(x_sample,y_sample,data_sample,ti,tj,sill,range,nugget)
                            fine0[i-location[0,isub],j-location[2,isub],iband]=(fine2[i,j,iband]*regress_result[iband,class_t,1]+regress_result[iband,class_t,0])+result[0]
                            error[i-location[0,isub],j-location[2,isub],iband]=result[1]
                        ENDFOR
                    ENDIF ELSE BEGIN
                        IF (c_gap GT 0 AND c_gap LT sample_size) THEN BEGIN
                            dis_window2=dis_window[ind_common]
                            x_sample=x_window[ind_common]
                            y_sample=y_window[ind_common]

                            FOR iband=0,nb-1,1 DO BEGIN
                                data_sample=(sub_dif[*,*,iband])[ind_common]
                                sill=result_semiv[0,iband,class_t]
                                range=result_semiv[2,iband,class_t]
                                nugget=result_semiv[1,iband,class_t]
                                result=o_kriging(x_sample,y_sample,data_sample,ti,tj,sill,range,nugget)
                                fine0[i-location[0,isub],j-location[2,isub],iband]=(fine2[i,j,iband]*regress_result[iband,class_t,1]+regress_result[iband,class_t,0])+result[0]
                                error[i-location[0,isub],j-location[2,isub],iband]=result[1]
                            ENDFOR
                        ENDIF ELSE BEGIN
                            ind_sameclass=WHERE(gap_wind[*,*] NE 1 AND sub_class EQ class_t, c_sameclass)
                            IF (c_sameclass GT 0) THEN BEGIN
                                dis_window2=dis_window[ind_sameclass]
                                x_window2=x_window[ind_sameclass]
                                y_window2=y_window[ind_sameclass]
                                ind_sample=WHERE(dis_window2 EQ MIN(dis_window2))
                                x_sample=x_window2[ind_sample]
                                y_sample=y_window2[ind_sample]
                                FOR iband=0,nb-1,1 DO BEGIN
                                    predict1=fine2[i,j,iband]*regress_result[iband,class_t,1]+regress_result[iband,class_t,0]
                                    predict2=mean((sub_off[*,*,iband])[x_sample,y_sample])
                                    fine0[i-location[0,isub],j-location[2,isub],iband]=(predict1+predict2)/2.0
                                    error[i-location[0,isub],j-location[2,isub],iband]=0
                                ENDFOR
                            ENDIF ELSE BEGIN
                                FOR iband=0,nb-1,1 DO BEGIN
                                    fine0[i-location[0,isub],j-location[2,isub],iband]=fine2[i,j,iband]*regress_result[iband,class_t,1]+regress_result[iband,class_t,0]
                                    error[i-location[0,isub],j-location[2,isub],iband]=0
                                ENDFOR
                            ENDELSE
                        ENDELSE
                    ENDELSE
                    FOR iband=0,nb-1,1 DO BEGIN          ;modify the abnormal values
                        IF ((fine0[i-location[0,isub],j-location[2,isub],iband]) LT DN_min OR (fine0[i-location[0,isub],j-location[2,isub],iband]) GT DN_max)THEN BEGIN
                            fine0[i-location[0,isub],j-location[2,isub],iband]=mean(image_series[i,j,iband,*])
                        ENDIF
                    ENDFOR
                ENDIF
            ENDFOR
        ENDFOR


        uncertainty=error/fine0*100.0      ;compute the uncertainty of prediction



        size_result=SIZE(fine0)

        CASE Data_Type OF
            1:fine0 = BYTE(fine0)    ;  BYTE  Byte
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
        ENDCASE



        ENVI_WRITE_ENVI_FILE,fine0,Out_Name = tempoutname1+STRTRIM(isub+1,1)
        ENVI_WRITE_ENVI_FILE,uncertainty,Out_Name = tempoutname2+STRTRIM(isub+1,1)
        ENVI_FILE_MNG, id=Fid1, /remove, /delete
        ENVI_FILE_MNG, id=Fid2, /remove, /delete

    ENDFOR

    ;--------------------------------------------------------------------------------------
    ;mosiac all the filled patch

    mfid=INTARR(n_ns*n_nl)
    mdims=INTARR(5,n_ns*n_nl)
    mpos=INTARR(nb,n_ns*n_nl)
    pos=INDGEN(nb)
    x0=INTARR(n_ns*n_nl)
    y0=INTARR(n_ns*n_nl)

    FOR isub=0,n_ns*n_nl-1,1 DO BEGIN
        ENVI_OPEN_FILE, tempoutname1+STRTRIM(isub+1,1), r_fid= sub_fid
        IF (sub_fid EQ -1) THEN BEGIN
            ENVI_BATCH_EXIT
            RETURN
        ENDIF
        ENVI_FILE_QUERY,  sub_fid, ns=sub_ns, nl=sub_nl
        mfid[isub] = sub_fid
        mpos[*,isub] = INDGEN(nb)
        mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
        x0[isub]=ind_patch1[0,isub]
        y0[isub]=ind_patch1[2,isub]
    ENDFOR

    xsize = orig_ns
    ysize = orig_nl
    pixel_size = [1.,1.]

    use_see_through = REPLICATE(1L,n_ns*n_nl)
    see_through_val = REPLICATE(0L,n_ns*n_nl)

    out_name=out_base + '_GNSPI.envi'
    ENVI_DOIT, 'mosaic_doit', fid=mfid, pos=mpos, $
        dims=mdims, out_name=out_name, xsize=xsize, $
        ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
        out_dt=Data_Type, pixel_size=pixel_size, $
        background=0, see_through_val=see_through_val, $
        use_see_through=use_see_through

    FOR i=0,n_ns*n_nl-1,1 DO BEGIN
        ENVI_FILE_MNG, id=mfid[i], /remove, /delete
    ENDFOR

    ;--------------------------------------------------------------------------------------
    ;mosiac all the uncertainty patch

    mfid=INTARR(n_ns*n_nl)
    mdims=INTARR(5,n_ns*n_nl)
    mpos=INTARR(nb,n_ns*n_nl)
    pos=INDGEN(nb)
    x0=INTARR(n_ns*n_nl)
    y0=INTARR(n_ns*n_nl)

    FOR isub=0,n_ns*n_nl-1,1 DO BEGIN
        ENVI_OPEN_FILE, tempoutname2+STRTRIM(isub+1,1), r_fid= sub_fid
        IF (sub_fid EQ -1) THEN BEGIN
            ENVI_BATCH_EXIT
            RETURN
        ENDIF
        ENVI_FILE_QUERY,  sub_fid, ns=sub_ns, nl=sub_nl
        mfid[isub] = sub_fid
        mpos[*,isub] = INDGEN(nb)
        mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
        x0[isub]=ind_patch1[0,isub]
        y0[isub]=ind_patch1[2,isub]
    ENDFOR

    xsize = orig_ns
    ysize = orig_nl
    pixel_size = [1.,1.]

    use_see_through = REPLICATE(1L,n_ns*n_nl)
    see_through_val = REPLICATE(0L,n_ns*n_nl)

    out_name=out_base+'_GNSPI_uncertainty.envi'
    ENVI_DOIT, 'mosaic_doit', fid=mfid, pos=mpos, $
        dims=mdims, out_name=out_name, xsize=xsize, $
        ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
        out_dt=Data_Type, pixel_size=pixel_size, $
        background=0, see_through_val=see_through_val, $
        use_see_through=use_see_through

    FOR i=0,n_ns*n_nl-1,1 DO BEGIN
        ENVI_FILE_MNG, id=mfid[i], /remove, /delete
    ENDFOR

    PRINT, 'time used:', FLOOR((SYSTIME(1)-t0)/3600), 'h',FLOOR(((SYSTIME(1)-t0) MOD 3600)/60),'m',(SYSTIME(1)-t0) MOD 60,'s'

END
