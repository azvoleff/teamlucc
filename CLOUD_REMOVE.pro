
;-------------------------------------------------------------------------------------------------
;                REMOVE THICK CLOUD IN SATELLITE IMAGES
;
;               This code can be used for whole TM scene
;                  Developed and coded by Zhu Xiaolin
;           Department of Geography,the Ohio State University
;                    Email:zhuxiaolin55@gmail.com
;                           2013-04-24
;     Reference: Zhu, X., Gao, F., Liu, D. and Chen, J. A modified
;     neighborhood similar pixel interpolator approach for removing
;     thick clouds in Landsat images. IEEE Geoscience and Remote
;     Sensing Letters, 2012, 9(3), 521-525
;     NOTE: to improve the efficiency for large image, some process in the paper
;     was simplified a little, but the accuracy is comparable
;-------------------------------------------------------------------------------------------------

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

pro CLOUD_REMOVE, cloudy_file, clear_file, mask_file, out_name, num_class, extent1, DN_min,$
    DN_max, patch_long

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

    temp_file=cloudy_file
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

        FileName = temp_file+'_temp_cloud'
        GetData,ImgData=ImgData,ns = ns1,nl = nl1,nb = nb,Data_Type = Data_Type,FileName = FileName+strtrim(isub+1,1),Fid = Fid1
        fine1=ImgData


        FileName = temp_file+'_temp_clear'
        GetData,ImgData=ImgData,ns = ns2,nl = nl2,nb = nb,FileName = FileName+strtrim(isub+1,1),Fid = Fid2
        fine2=ImgData

        FileName = temp_file+'_temp_cloud_mask'
        GetData,ImgData=ImgData,FileName = FileName+strtrim(isub+1,1),Fid = Fid3
        cloud=ImgData


        ;-----------------------------------------------------------------------


        cloud_posite=where(cloud ne 0,c_c)

        if (c_c gt 0) then begin

            num_cloud_path=max(cloud[cloud_posite])
            ;
            for icloud=1,num_cloud_path,1 do begin       ; remove each could patches
                cloud_area=where(cloud eq icloud,num_cloud1)
                if (num_cloud1 gt 0) then begin

                    cloud_pixel_locate=intarr(2,num_cloud1)
                    for ic=0.0,num_cloud1-1.0,1.0 do begin
                        cloud_pixel_locate[1,ic]=fix(cloud_area[ic]/ns1)
                        cloud_pixel_locate[0,ic]=cloud_area[ic]-ns1*cloud_pixel_locate[1,ic]
                    endfor
                    left_cloud=min(cloud_pixel_locate[0,*])
                    right_cloud=max(cloud_pixel_locate[0,*])
                    up_cloud=min(cloud_pixel_locate[1,*])
                    down_cloud=max(cloud_pixel_locate[1,*])
                    cloud_area=0 ;clear it
                    cloud_pixel_locate=0
                    ;

                    left_region=max([0,left_cloud-extent1])            ;neigborhood of cloud area: rectangle which can completely cover cloud
                    right_region=min([ns1-1,right_cloud+extent1])
                    up_region=max([0,up_cloud-extent1])
                    down_region=min([nl1-1,down_cloud+extent1])

                    fine1_sub=fine1[left_region:right_region,up_region:down_region,*]
                    a_region=right_region-left_region+1
                    b_region=down_region-up_region+1
                    x_center=a_region/2.0      ; the position of cloud center
                    y_center=b_region/2.0
                    fine2_sub=fine2[left_region:right_region,up_region:down_region,*]
                    cloud_sub=cloud[left_region:right_region,up_region:down_region]

                    fine2_1=reform(fine2_sub,long(a_region)*b_region,nb)    ; classify the cloud-free image, used for selecting similar pixels
                    fine2_1=transpose(fine2_1)
                    weights = CLUST_WTS(fine2_1, N_CLUSTERS =num_class)
                    class_img=CLUSTER(fine2_1, weights, N_CLUSTERS =num_class)
                    class_img=reform(class_img,a_region,b_region)
                    fine2_1=0;clear it
                    weights=0;clear it


                    ;----------------------------------------------------------------------

                    indcloud_sub=where(cloud_sub eq icloud,num_cloud)

                    fine02=fine1_sub
                    ;---------------------remove background-----------------
                    fine1_value=cloud_sub
                    fine1_value[*,*]=0
                    fine2_value=fine1_value
                    for ib=0,nb-1,1 do begin
                        fine1_value=mean(fine1_value+fine1_sub[*,*,ib])
                        fine2_value=mean(fine2_value+fine2_sub[*,*,ib])
                    endfor
                    ind_com=where(fine1_value eq 0 or fine2_value eq 0, c_com)
                    if (c_com gt 0) then begin
                        cloud_sub[ind_com]=-1
                    endif
                    fine1_value=0
                    fine2_value=0


                    for ic=0,num_class-1,1 do begin
                        ind_clear=where(cloud_sub eq 0 and class_img eq ic,c_clear)
                        ind_cloud=where(cloud_sub eq icloud and class_img eq ic,c_cloud)

                        if(c_clear gt 0 and c_cloud gt 0) then begin
                            fine2_clear=fltarr(c_clear,nb)
                            fine1_clear=fltarr(c_clear,nb)
                            for ib=0,nb-1,1 do begin
                                fine2_clear[*,ib]=(fine2_sub[*,*,ib])[ind_clear]
                                fine1_clear[*,ib]=(fine2_sub[*,*,ib])[ind_clear]
                            endfor
                            xy_clear=fltarr(2,c_clear)                    ;the location of each clear pixel
                            for iclear=0.0,c_clear-1.0,1.0 do begin
                                xy_clear[1,iclear]=fix(ind_clear[iclear]/a_region)
                                xy_clear[0,iclear]=ind_clear[iclear]-a_region*xy_clear[1,iclear]
                            endfor
                            for ip_c=0,c_cloud-1.0, 1l do begin
                                ri=fix(ind_cloud[ip_c]/a_region)    ; the location of target cloudy pixel
                                ci=ind_cloud[ip_c]-a_region*ri
                                s_dis=((xy_clear[0,*]-ci)^2+(xy_clear[1,*]-ri)^2)^0.5 ; spatial distance of all clear pixels to the target pixel
                                order_dis=sort(s_dis)
                                number_similar=min([c_clear,20])
                                s_d=fltarr(number_similar)
                                s_rmsd=fltarr(number_similar)
                                for i_similar=0, number_similar-1,1 do begin
                                    ;caculate RMSD of each same-class pixel to the target cloudy pixel
                                    s_rmsd[i_similar]=((total((fine2_clear[order_dis[i_similar],*]-fine1_sub[ci,ri,*])^2))/float(nb))^0.5
                                    s_d[i_similar]=s_dis[order_dis[i_similar]]
                                endfor
                                C_D=s_d*s_rmsd+0.0000001
                                weight=(1.0/C_D)/total(1.0/C_D)
                                r1=mean(s_d)                      ; the weight of two predictions
                                r2=((ci-x_center)^2+(ri-y_center)^2)^0.5
                                W_T1=(1/r1)/(1/r1+1/r2)
                                W_T2=(1/r2)/(1/r1+1/r2)

                                for iband=0,nb-1,1 do begin
                                    fine1_similar=(fine1_clear[*,iband])[order_dis[0:(number_similar-1)]]
                                    fine2_similar=(fine2_clear[*,iband])[order_dis[0:(number_similar-1)]]
                                    predict_1=total(fine1_similar*weight)
                                    predict_2=fine2_sub[ci,ri,iband]+total((fine1_similar-fine2_similar)*weight)
                                    if (predict_2 gt DN_min and predict_2 lt DN_max) then begin
                                        fine1_sub[ci,ri,iband]=W_T1*predict_1+W_T2*predict_2
                                    endif else begin
                                        fine1_sub[ci,ri,iband]=predict_1
                                    endelse
                                endfor
                            endfor

                        endif else begin
                            if (c_cloud gt 0 and c_clear eq 0 ) then begin
                                fine1_sub[ind_cloud,*]=fine2_sub[ind_cloud,*]
                            endif
                        endelse
                    endfor

                    fine1[left_region:right_region,up_region:down_region,*]=fine1_sub
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
