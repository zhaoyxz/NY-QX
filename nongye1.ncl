begin
    filename1 = "E:/Documents/ncl/data/tmax.2020.nc"
    filename2 = "E:/Documents/ncl/data/tmin.2020.nc"
    f1 = addfile(filename1, "r")
    f2 = addfile(filename2, "r")

    lon = f1->lon
    lat = f1->lat

    tmax = f1->tmax
    tmin = f2->tmin
    tavg0 = (tmax+tmin)/2
    tavg = runave_n_Wrap(tavg0, 5, -1, 0)
    tavg1 = tavg
    tavg2 = tavg

    do i = 68,175,1
        do j = 144,273,1
            if(ismissing(tavg1(0,i,j)))then
            tavg1(:,i,j) = -5
            tavg2(:,i,j) = 10
        end if
        end do
    end do

    a = new((/360,720/),float,"No_FillValue")
    a = 0
    a!0 = "lat"
    a!1 = "lon"
    a&lat = lat
    a&lon = lon
    a@units = "First_Day"
    b = new((/360,720/),float,"No_FillValue")
    b = 0
    b!0 = "lat"
    b!1 = "lon"
    b&lat = lat
    b&lon = lon
    b@units = "total_Day"

    do i = 68,175,1
        do j = 144,273,1
            do t = 0,365,1
            if(tavg1(t,i,j).ge.5)then
                b(i,j) = b(i,j)+1
            end if
            end do
        end do
    end do

    do i = 68,175,1
        do j = 144,273,1
            do t = 0,365,1
            if(tavg2(t,i,j).lt.5)then
                a(i,j) = a(i,j)+1
            continue
            end if
            end do
        end do
    end do
 
    wks  = gsn_open_wks("pdf", "E:/Documents/ncl/save/1.2")
    map = new(2,"graphic")
    map_nanhai = new(2,"graphic")
    plotnh = new(2,"graphic")

    res                       = True
    ;res@gsnMaximize           = True
    res@gsnDraw               = False
    res@gsnFrame              = False
    
    res@mpFillOn              = True
    res@cnLinesOn             = False
    res@cnFillOn              = True 

    res@pmTickMarkDisplayMode = "Always"  ;坐标轴签上增加°符号
    res@cnLevelSelectionMode  = "ManualLevels"
    res@cnFillPalette         = "MPL_jet"
    res@cnFillDrawOrder       = "PreDraw"
    
    res@mpGridAndLimbOn       = True
    res@mpGridLineDashPattern = 15
    res@mpGridSpacingF        = 10

    
    res@cnMinLevelValF        =0
    res@cnMaxLevelValF        =360
    res@cnLevelSpacingF       =20
   
    res@mpAreaMaskingOn        = True
    res@mpOutlineBoundarySets  = "NoBoundaries"
    res@mpFillBoundarySets     = "NoBoundaries"
    res@mpDataBaseVersion      = "MediumRes"
    res@mpDataSetName          = "Earth..4"
    res@mpMaskAreaSpecifiers   = (/"China"/)
    res@mpOutlineSpecifiers    = (/"China","China:Provinces"/)
    ;res@mpFillDrawOrder        = "PostDraw" 
    res@mpFillAreaSpecifiers   = (/"China","water","land"/)
    res@mpInlandWaterFillColor = 0
    res@mpOceanFillColor       = 0
    res@mpNationalLineThicknessF = 0.2
    res@mpProvincialLineThicknessF = 0.1

    nhres                     = res

    res@mpMaxLatF             = 55
    res@mpMinLatF             = 17
    res@mpMaxLonF             = 136
    res@mpMinLonF             = 72  
    res1 = res
    res2 = res
    res1@gsnCenterString = "start_day-pass-5~F35~J~F~C"
    res2@gsnCenterString = "total_day-pass-5~F35~J~F~C"

;--- add South China Sea --- 
  
    ;nhres@gsnMaximize              = False

    nhres@vpHeightF                = 0.15    
    nhres@vpWidthF                 = 0.15
    
    nhres@mpMinLatF                =   2   
    nhres@mpMaxLatF                =  23
    nhres@mpMinLonF                = 105
    nhres@mpMaxLonF                = 123

    nhres@lbLabelBarOn             = False
    nhres@tmXBOn                   = False
    nhres@tmXTOn                   = False
    nhres@tmYLOn                   = False
    nhres@tmYROn                   = False
    nhres@gsnLeftString            = ""
    nhres@gsnRightString           = ""

    adres                          = True
    adres@amParallelPosF           = 0.495 ; -0.5 is the left edge of the plot.
    adres@amOrthogonalPosF         = 0.49  ; -0.5 is the top edge of the plot.
    adres@amJust                   = "BottomRight"

    map(0) = gsn_csm_contour_map(wks,a,res1) 
    map_nanhai(0) = gsn_csm_contour_map(wks,a,nhres)
    plotnh(0) = gsn_add_annotation(map(0),map_nanhai(0),adres)  

    draw(map(0))
    frame(wks) 

    map(1) = gsn_csm_contour_map(wks,b,res2)
    map_nanhai(1) = gsn_csm_contour_map(wks,b,nhres)
    plotnh(1) = gsn_add_annotation(map(1),map_nanhai(1),adres) 

    draw(map(1))
    frame(wks) 

end