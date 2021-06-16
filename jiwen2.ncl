begin

filename1 = "E:/Documents/ncl/data/tmax.2020.nc"
filename2 = "E:/Documents/ncl/data/tmin.2020.nc"
f1 = addfile(filename1, "r")
f2 = addfile(filename2, "r")

lon = f1->lon
lat = f1->lat
time = f1->time

tmax = f1->tmax
tmin = f2->tmin
tavg = (tmax+tmin)/2

a = new((/360,720/),float,"No_FillValue")
a = 0
a!0 = "lat"
a!1 = "lon"
a&lat = lat
a&lon = lon
a@units = "~F35~J~F~C"

do i = 68,175,1
    do j = 144,273,1
        if(ismissing(tavg(0,i,j)))then
        tavg(:,i,j) = -5
    end if
    end do
end do

do i = 68,175,1
    do j = 144,273,1
        s = 0.0
        do k = 0,365,1
        if(tavg(k,i,j).gt.10)then
            s = s + tavg(k,i,j)
        else
            continue 
        end if
        a(i,j) = s
        end do
    end do
end do
printMinMax(a, 0)
wks  = gsn_open_wks("pdf", "E:/Documents/ncl/save/1.311")
res                       = True
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnCenterString       = "10~F35~J~F~C Active Accumulated Temperature"
res@gsnCenterStringFontHeightF = 0.025
res@gsnCenterStringOrthogonalPosF = 0.06

res@mpFillOn              = True
res@cnLinesOn             = False
res@cnFillOn              = True 

res@pmTickMarkDisplayMode = "Always"  ;坐标轴签上增加°符号
res@cnLevelSelectionMode  = "ManualLevels"
res@cnFillPalette         = "WhiteBlueGreenYellowRed"
res@cnFillDrawOrder       = "PreDraw"

res@mpGridAndLimbOn       = True
res@mpGridLineDashPattern = 15
res@mpGridSpacingF        = 10

res@cnMinLevelValF        =0
res@cnMaxLevelValF        =10000
res@cnLevelSpacingF       =500

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

map = gsn_csm_contour_map(wks,a,res) 

nhres@tiMainString = ""
nhres@gsnCenterString       = ""
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

map_nanhai = gsn_csm_contour_map(wks,a,nhres)

adres                          = True
adres@amParallelPosF           = 0.495 ; -0.5 is the left edge of the plot.
adres@amOrthogonalPosF         = 0.49  ; -0.5 is the top edge of the plot.
adres@amJust                   = "BottomRight"

plotnh = gsn_add_annotation(map,map_nanhai,adres)  

draw(map)
frame(wks) 
end