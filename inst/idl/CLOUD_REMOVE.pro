
;-------------------------------------------------------------------------------------------------
;                REMOVE THICK CLOUD IN SATELLITE IMAGES
;
;               This code can be used for whole TM scene
;                  Developed and coded by Zhu Xiaolin
;           Department of Geography,the Ohio State University
;                    Email:zhuxiaolin55@gmail.com
;                           updated：2013-10-23
;     Debug history:
;         2013-10-23 correct:1) the bug of similar pixel searching
;                            2) the bug of using wrong values in cloudy image
;
;     Please cite: Zhu, X., Gao, F., Liu, D. and Chen, J. A modified
;     neighborhood similar pixel interpolator approach for removing
;     thick clouds in Landsat images. IEEE Geoscience and Remote
;     Sensing Letters, 2012, 9(3), 521-525
;     NOTE: the efficiency may be low for large image.If so, please use CLOUD_REMOVE_FAST,
;     in which some process was simplified, but the accuracy is still satisfactory
;-----------------------------------------------------------------------------
;
; Modified from original by Xiaolin Zhu on Feb 6, 2014 by Alex Zvoleff for 
; inclusion in the 'teamr' R package. For the original version, see Xiaolin's 
; website at: http://geography.osu.edu/grads/xzhu/
;
;-----------------------------------------------------------------------------

;function of open file
Pro GetData,ImgData = ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,$
    FileName = FileName,Map_info = map_Info, Fid = Fid

    COMPILE_OPT idl2, hidden
    e = ENVI(/HEADLESS)
    ENVI_BATCH_STATUS_WINDOW, /ON

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
;-------------------------------------------------------------------
;                  main body of the program
;-------------------------------------------------------------------

pro CLOUD_REMOVE, cloudy_file, clear_file, mask_file, out_name, num_class, $
    min_pixel, extent1, DN_min, DN_max, patch_long

    COMPILE_OPT idl2, hidden
    e = ENVI(/HEADLESS)
    ENVI_BATCH_STATUS_WINDOW, /ON

    t0=systime(1)

    envi_open_file,cloudy_file,r_fid=fid
    envi_file_query,fid,ns=ns,nl=nl,nb=nb,dims=dims
    map_info = envi_get_map_info(fid=fid)
    orig_ns=ns
    orig_nl=nl
    n_ns=ceil(float(ns)/patch_long)
    n_nl=ceil(float(nl)/patch_long)

    ind_patch=intarr(4,n_ns*n_nl)

    for i_ns=0,n_ns-1,1 do begin
        for i_nl=0,n_nl-1,1 do begin
            ind_patch[0,n_ns*i_nl+i_ns]=i_ns*patch_long
            ind_patch[1,n_ns*i_nl+i_ns]=min([ns-1,(i_ns+1)*patch_long-1])
            ind_patch[2,n_ns*i_nl+i_ns]=i_nl*patch_long
            ind_patch[3,n_ns*i_nl+i_ns]=min([nl-1,(i_nl+1)*patch_long-1])
        endfor
    endfor

    ; Below line is necessary to remove file extensions from the middle of temp 
    ; filenames
    cloudy_file_split=strsplit(cloudy_file, '.', /EXTRACT)
    temp_file=cloudy_file_split[0]
    tempoutname=temp_file+'_temp_cloud'

    pos=indgen(nb)
    for isub=0,n_ns*n_nl-1,1 do begin
        dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
        envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
            out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
        envi_file_mng, id=r_fid1, /remove
    endfor

    envi_file_mng, id=fid, /remove

    ;-----------------------------------------------------------
    envi_open_file,clear_file,r_fid=fid
    tempoutname=temp_file+'_temp_clear'
    pos=indgen(nb)
    for isub=0,n_ns*n_nl-1,1 do begin
        dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
        envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
            out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
        envi_file_mng, id=r_fid1, /remove
    endfor
    envi_file_mng, id=fid, /remove

    ;------------------------------------------------------------
    envi_open_file,mask_file,r_fid=fid
    tempoutname=temp_file+'_temp_cloud_mask'
    pos=indgen(1)
    for isub=0,n_ns*n_nl-1,1 do begin
        dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
        envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
            out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
        envi_file_mng, id=r_fid1, /remove
    endfor
    envi_file_mng, id=fid, /remove

    ;-------------------------------------------------------------------
    tempoutname=temp_file+'_temp_filled'

    print,'there are total',n_ns*n_nl,'patches'

    for isub=0,n_ns*n_nl-1,1 do begin

        ;open each block

        fine1_FileName = temp_file + '_temp_cloud' + strtrim(isub+1,1)
        GetData,ImgData=fine1,ns=ns1,nl=nl1,nb=nb,Data_Type=Data_Type,FileName=fine1_FileName,Fid = Fid1

        fine2_FileName = temp_file+'_temp_clear'+strtrim(isub+1,1)
        GetData,ImgData=fine2,FileName=fine2_FileName,Fid=Fid2

        cloud_FileName = temp_file+'_temp_cloud_mask'+strtrim(isub+1,1)
        GetData,ImgData=cloud,FileName=cloud_FileName,Fid=Fid3


        ;-----------------------------------------------------------------------


        cloud_posite=where(cloud ne 0,c_c)

        if (c_c gt 0) then begin

            num_cloud_path=max(cloud[cloud_posite])
            ;
            for icloud=1,num_cloud_path,1 do begin       ; remove each cloud patches
                cloud_area=where(cloud eq icloud,num_cloud1)

                if (num_cloud1 gt 0) then begin

                    cloud_pixel_locate=intarr(2,num_cloud1)
                    cloud_pixel_locate[1,*]=fix(cloud_area/ns1)
                    cloud_pixel_locate[0,*]=cloud_area-ns1*cloud_pixel_locate[1,*]
                    left_cloud=min(cloud_pixel_locate[0,*])     ;find out the cloud area
                    right_cloud=max(cloud_pixel_locate[0,*])
                    up_cloud=min(cloud_pixel_locate[1,*])
                    down_cloud=max(cloud_pixel_locate[1,*])
                    cloud_area=0 ;clear it
                    cloud_pixel_locate=0 ;clear it

                    ;neigborhood of cloud area: rectangle which can completely cover cloud
                    left_region=max([0,left_cloud-extent1])
                    right_region=min([ns1-1,right_cloud+extent1])
                    up_region=max([0,up_cloud-extent1])
                    down_region=min([nl1-1,down_cloud+extent1])


                    a_region=right_region-left_region+1
                    b_region=down_region-up_region+1
                    x_center=a_region/2.0      ; the position of cloud center
                    y_center=b_region/2.0


                    sub_fine1=fine1[left_region:right_region,up_region:down_region,*]
                    sub_fine2=fine2[left_region:right_region,up_region:down_region,*]
                    sub_cloud=cloud[left_region:right_region,up_region:down_region]

                    ;---------------------remove background-----------------
                    fine1_value=sub_cloud
                    fine1_value[*,*]=0
                    fine2_value=fine1_value
                    for ib=0,nb-1,1 do begin
                        fine1_value=mean(fine1_value+sub_fine1[*,*,ib])
                        fine2_value=mean(fine2_value+sub_fine2[*,*,ib])
                    endfor
                    ind_com=where(fine1_value eq 0 or fine2_value eq 0, c_com)
                    if (c_com gt 0) then begin
                        sub_cloud[ind_com]=-1          ;mask out background
                    endif
                    fine1_value=0
                    fine2_value=0

                    similar_th_band=fltarr(nb)
                    for iband=0,nb-1,1 do begin
                        similar_th_band[iband]=stddev(sub_fine2[*,*,iband])*2.0/float(num_class)    ;compute the threshold of similar pixel
                    endfor

                    ind_cloud=where(sub_cloud eq icloud,num_cloud)
                    ind_clear=where(sub_cloud eq 0,num_clear)
                    xy_clear=fltarr(2,num_clear)                    ;the location of each clear pixel
                    xy_clear[1,*]=fix(ind_clear/(a_region))
                    xy_clear[0,*]=ind_clear-(a_region)*xy_clear[1,*]
                    fine2_clear=fltarr(num_clear,nb)
                    fine1_clear=fltarr(num_clear,nb)
                    for ib=0,nb-1,1 do begin
                        fine2_clear[*,ib]=(sub_fine2[*,*,ib])[ind_clear]
                        fine1_clear[*,ib]=(sub_fine1[*,*,ib])[ind_clear]
                    endfor

                    for ic=0.0,num_cloud-1.0,1.0 do begin
                        ri=fix(ind_cloud[ic]/(a_region))    ; the location of target pixel
                        ci=ind_cloud[ic]-(a_region)*ri
                        center_dis=((x_center-ci)^2+(y_center-ri)^2)^0.5; the distance between target pixel and center of cloud
                        s_dis=((xy_clear[0,*]-ci)^2+(xy_clear[1,*]-ri)^2)^0.5

                        order_clear=sort(s_dis,/L64)

                        iclear=0
                        isimilar=0
                        fine1_similar=fltarr(min_pixel,nb)
                        fine2_similar=fltarr(min_pixel,nb)
                        rmse_similar=fltarr(min_pixel)
                        dis_similar=fltarr(min_pixel)
                        similar_rmse12=fltarr(min_pixel)
                        while (isimilar le min_pixel-1 and iclear le num_clear-1.0 ) do begin
                            indicate_similar=0
                            for iband=0,nb-1,1 do begin
                                if (abs(fine2_clear[order_clear[iclear],iband]-sub_fine2[ci,ri,iband])le similar_th_band[iband]) then begin
                                    indicate_similar=indicate_similar+1
                                endif
                            endfor
                            if (indicate_similar eq nb) then begin
                                for ib=0, nb-1,1 do begin
                                    fine1_similar[isimilar,ib]=(fine1_clear[*,ib])[order_clear[iclear]]
                                    fine2_similar[isimilar,ib]=(fine2_clear[*,ib])[order_clear[iclear]]
                                endfor
                                rmse_similar[isimilar]=((total((fine2_clear[order_clear[iclear],*]-sub_fine2[ci,ri,*])^2))/float(nb))^0.5
                                dis_similar[isimilar]=s_dis[order_clear[iclear]]
                                similar_rmse12[isimilar]=((total((fine2_clear[order_clear[iclear],*]-fine1_clear[order_clear[iclear],*])^2))/float(nb))^0.5
                                isimilar=isimilar+1
                            endif
                            iclear=iclear+1
                        endwhile

                        ; if similar pixel larger than 0
                        if (isimilar gt 0 ) then begin
                            ind_null=where(dis_similar gt 0)
                            rmse_similar1=(rmse_similar[ind_null]-min(rmse_similar[ind_null]))/(max(rmse_similar[ind_null])-min(rmse_similar[ind_null])+0.000001)+1.0
                            dis_similar1=(dis_similar[ind_null]-min(dis_similar[ind_null]))/(max(dis_similar[ind_null])-min(dis_similar[ind_null])+0.000001)+1.0
                            C_D=(rmse_similar1)*dis_similar1+0.0000001
                            weight=(1.0/C_D)/total(1.0/C_D)

                            ;compute the time weight
                            W_T1=center_dis/(center_dis+mean(dis_similar[ind_null]))
                            W_T2=mean(dis_similar[ind_null])/(center_dis+mean(dis_similar[ind_null]))
                            for iband=0,nb-1,1 do begin
                                predict_1=total((fine1_similar[*,iband])[ind_null]*weight)
                                predict_2=sub_fine2[ci,ri,iband]+total(((fine1_similar[*,iband])[ind_null]-(fine2_similar[*,iband])[ind_null])*weight)
                                if (predict_2 gt DN_min and predict_2 lt DN_max) then begin
                                    sub_fine1[ci,ri,iband]=W_T1*predict_1+W_T2*predict_2
                                endif else begin
                                    sub_fine1[ci,ri,iband]=predict_1
                                endelse
                            endfor
                        endif else begin
                            ; if no similar pixel, use all similar pixels
                            for iband=0,nb-1,1 do begin
                                diff=mean(fine1_clear[*,iband]-fine2_clear[*,iband])
                                sub_fine1[ci,ri,iband]=sub_fine2[ci,ri,iband]+diff
                            endfor
                        endelse

                    endfor
                    fine1[left_region:right_region,up_region:down_region,*]=sub_fine1
                endif
            endfor
        endif
        print,'finished ',isub+1,' patch'

        Envi_Write_Envi_File,fine1,Out_Name = tempoutname+strtrim(isub+1,1)
        envi_file_mng, id=Fid1, /remove, /delete
        envi_file_mng, id=Fid2, /remove, /delete
        envi_file_mng, id=Fid3, /remove, /delete

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

        envi_open_file, tempoutname+strtrim(isub+1,1), r_fid= sub_fid
        if (sub_fid eq -1) then begin
            envi_batch_exit
            return
        endif

        envi_file_query,  sub_fid, ns=sub_ns, nl=sub_nl

        mfid[isub] = sub_fid
        mpos[*,isub] = indgen(nb)
        mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
        x0[isub]=ind_patch[0,isub]
        y0[isub]=ind_patch[2,isub]
    endfor


    xsize = orig_ns
    ysize = orig_nl
    pixel_size = [1.,1.]


    use_see_through = replicate(1L,n_ns*n_nl)
    see_through_val = replicate(0L,n_ns*n_nl)

    envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
        dims=mdims, out_name=out_name, xsize=xsize, $
        ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
        out_dt=Data_Type, pixel_size=pixel_size, $
        background=0, see_through_val=see_through_val, $
        use_see_through=use_see_through


    for i=0,n_ns*n_nl-1,1 do begin
        envi_file_mng, id=mfid[i], /remove, /delete
    endfor

    print, 'time used', floor((systime(1)-t0)/3600), 'hour',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'

end
