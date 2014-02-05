;-------------------------------------------------------------------------------------------------
;                REMOVE THICK CLOUD IN SATELLITE IMAGES
;                          * A FSAT VERSION*
;
;               Faster VERSION: can be used for whole TM scene
;                  Developed and coded by Zhu Xiaolin
;           Department of Geography,the Ohio State University
;                    Email:zhuxiaolin55@gmail.com
;                           2013-04-24
;     Reference: Zhu, X., Gao, F., Liu, D. and Chen, J. A modified 
;     neighborhood similar pixel interpolator approach for removing 
;     thick clouds in Landsat images. IEEE Geoscience and Remote 
;     Sensing Letters, 2012, 9(3), 521-525
;     NOTE: to improve the efficiency for large image, some process in the paper 
;     was LARGELY simplified, but the accuracy is still acceptable
;---------------------------------------------------------------------------


;function of open file

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
;-------------------------------------------------------------------
;                  main body of the program
;-------------------------------------------------------------------

pro  CLOUD_REMOVE_FAST

 t0=systime(1)

;------------------------------------------------
 num_class=4  ;set the estimated number of classes in image
 extent1=1    ;set the range of cloud neighborhood
 DN_min=0     ;set the range of DN
 DN_max=255
 patch_long=1000 ;set the block size
 ;-------------------------------------------------

  FileName1 = Dialog_PickFile(title = 'Open the cloudy image:')
  envi_open_file,FileName1,r_fid=fid
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

  temp_file=FileName1
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
  FileName2 = Dialog_PickFile(title = 'Open the input clear image:')
  envi_open_file,FileName2,r_fid=fid
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
  FileName3 = Dialog_PickFile(title = 'Open the cloud mask image:')
  envi_open_file,FileName3,r_fid=fid
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
   for icloud=1,num_cloud_path,1 do begin
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

         left_region=max([0,left_cloud-extent1])            ;region of cloud area
         right_region=min([ns1-1,right_cloud+extent1])
         up_region=max([0,up_cloud-extent1])
         down_region=min([nl1-1,down_cloud+extent1])

          fine1_sub=fine1[left_region:right_region,up_region:down_region,*]
          a_region=right_region-left_region+1
          b_region=down_region-up_region+1
          fine2_sub=fine2[left_region:right_region,up_region:down_region,*]
          cloud_sub=cloud[left_region:right_region,up_region:down_region]

         fine2_1=reform(fine2_sub,long(a_region)*b_region,nb)
         fine2_1=transpose(fine2_1)
         weights = CLUST_WTS(fine2_1, N_CLUSTERS =num_class)
         class_img=CLUSTER(fine2_1, weights, N_CLUSTERS =num_class)
         class_img=reform(class_img,a_region,b_region)
          fine2_1=0;clear it
          weights=0;clear it


          ;----------------------------------------------------------------------

          indcloud_sub=where(cloud_sub eq icloud,num_cloud)

            relation=fltarr(num_class,2,nb)
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

           ;USE temporal information
           for ic=0,num_class-1,1 do begin
             ind_clear=where(cloud_sub eq 0 and class_img eq ic,c_clear)
             ind_cloud=where(cloud_sub eq icloud and class_img eq ic,c_cloud)
             if(c_clear gt 0 and c_cloud gt 0) then begin
               for ib=0,nb-1,1 do begin
                   x_fine2=float((fine2_sub[*,*,ib])[ind_clear])
                   y_fine1=float((fine1_sub[*,*,ib])[ind_clear])
                   relation[ic,*,ib]=LINFIT(x_fine2,y_fine1,YFIT=y_L_f)
                   if (abs(relation[ic,1,ib]) lt 0.33 or abs(relation[ic,1,ib]) gt 3) then begin
                      relation[ic,1,ib]=1.0
                      relation[ic,0,ib]=mean(y_fine1-x_fine2)
                   endif
                   temp=float(fine02[*,*,ib])
                   temp[ind_cloud]=relation[ic,0,ib]+relation[ic,1,ib]*(fine2_sub[*,*,ib])[ind_cloud]
                   temp1=temp[ind_cloud]
                   ind0=where(temp1 lt DN_min, c0)
                   if (c0 gt 0) then begin
                     temp1[ind0]=((fine2_sub[*,*,ib])[ind_cloud])[ind0]
                    endif
                  ind255=where(temp1 gt DN_max, c255)
                     if (c255 gt 0) then begin
                       temp1[ind255]=((fine2_sub[*,*,ib])[ind_cloud])[ind255]
                     endif
                     temp[ind_cloud]=temp1
                     fine02[*,*,ib]=temp
                 endfor
               endif else begin
                 if (c_cloud gt 0 and c_clear eq 0 ) then begin
                     fine02[ind_cloud,*]=fine2_sub[ind_cloud,*]
                  endif
               endelse
             endfor
            ind_clear=0
            x_fine2=0
            y_fine1=0
            temp=0
            temp1=0

            ;use spatial information
            for i=0.0,num_cloud-1.0,1.0 do begin

               ri=fix(indcloud_sub[i]/a_region)
               ci=indcloud_sub[i]-a_region*ri

                class_i=class_img[ci,ri]
                a1=max([0,ci-5])
                a2=min([a_region-1,ci+5])
                b1=max([0,ri-5])
                b2=min([b_region-1,ri+5])
                class_w=class_img[a1:a2,b1:b2]
                cloud_sub_w=cloud_sub[a1:a2,b1:b2]
                ind_w_sameclass=where(class_w eq class_i,c_classi)
               ind_w_clear=where(cloud_sub_w eq 0 and  class_w eq class_i,c_uncloud)
                w1=c_uncloud/float(c_classi)
                w2=1.0-w1
                ; combination of the two predictions by weight
                for ib=0,nb-1,1 do begin
                    value_fine1=(fine1_sub[*,*,ib])[a1:a2,b1:b2]
                      if (c_uncloud gt 0) then begin
                           fine01=mean(value_fine1[ind_w_clear])
                           fine1_sub[ci,ri,ib]=fine01*w1+fine02[ci,ri,ib]*w2
                       endif else begin
                            fine1_sub[ci,ri,ib]=fine02[ci,ri,ib]
                       endelse
                 endfor

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

;   out_name=Dialog_PickFile(Title = 'Enter the filename of the restored image')
    out_name=FileName1+'cloud_remove_fast'
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