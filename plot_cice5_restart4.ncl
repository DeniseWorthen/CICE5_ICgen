;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;************************************************
begin

  batch = "no"

  hemi="NH"
  ;hemi="SH"

  if(hemi .eq. "NH")then
    minlat =  50.0
    midlon = -45.0
  else
    maxlat = -50.0
    midlon =   0.0
  end if

;************************************************
; read in netCDF file
;************************************************

  ; hera
  RT1 = "/scratch2/NCEPDEV/stmp1/Denise.Worthen/"
  RT2 = "/scratch1/NCEPDEV/stmp2/Denise.Worthen/"
  RT3 = "/scratch2/NCEPDEV/stmp3/Denise.Worthen/"
  RT4 = "/scratch1/NCEPDEV/stmp4/Denise.Worthen/"

  ncat = 5

;************************************************
; read in netCDF file
;************************************************

  cyear="2012"
  cmons = (/"01","04","07","10"/)
  dates = cyear+cmons+"0100"
  
  fsrc="/scratch2/NCEPDEV/climate/Denise.Worthen/ICmx050/"+dates+"/"
  fres="cice5_model_0.50.res_"+dates+".nc"
  files=addfiles(fsrc+fres,"r")

  nems="/scratch2/NCEPDEV/climate/Denise.Worthen/NEMS_INPUT0.1/MOM6_FIX_05deg/"
  gf=addfile(nems+"ocean_static.nc","r")

  icelonname = "geolon"
  icelatname = "geolat"

  wkstitle="restarts"
;************************************************
;
;************************************************

  rhoi = 917.0
  rhos = 330.0

;************************************************
; add a poly line showing the polar seem along
; great circle path from (lon1,lat1) to (lon2,lat2)
;************************************************

  lon1=  60.00
  lon2=-120.00
  lat1 = 60.00
  lat2 = 60.00

;************************************************
; define the map
;************************************************

  if(batch .eq. "no")then
   wkstype="X11"
  else
   wkstype="png"
  end if
  ; wkstype="png"
  wkstype@wkWidth=800
  wkstype@wkHeight=800
  wks   = gsn_open_wks (wkstype,wkstitle)

  ; panels
  npanels = dimsizes(dates) 
  plot = new(npanels,graphic)
   dum = new(npanels,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@mpFillOn   = True
  res@cnFillMode           = "MeshFill"
  ;res@cnFillMode           = "RasterFill"
  res@cnFillDrawOrder      = "PreDraw"             ; draw contours first

  res@cnFillOn           = True                  ; color fill
  res@cnLinesOn          = False                 ; no contour lines
  res@lbLabelBarOn       = True                 ; turn off individual colorbars
  res@lbLabelFontHeightF = 0.015

  res@gsnStringFontHeightF = 0.05

   res@gsnPolar   =  hemi                         ; specify the hemisphere
  if(hemi .eq. "NH")then
   res@mpMinLatF    =  minlat                    ; minimum lat to plot
  else
   res@mpMaxLatF    = maxlat
  end if
   res@mpCenterLonF = midlon

  res@mpFillOn   = False
  res@mpOutlineOn          = False
  ;res@cnMissingValFillPattern = 0               ; choose fill pattern 0 (solid fill)
  ;res@cnMissingValFillColor   = "gray"          ; color fill missing areas blue
 
;************************************************
; Resources for the polyline
;************************************************

  res@mpGreatCircleLinesOn = True
  pres                  = True           ; polyline mods desired
  pres@gsLineThicknessF = 2.0            ; line thickness
  pres@gsLineColor      = "black"         ; color of lines
  pres@gsLineDashPattern =  1            ; dashes
  ;pres@tfPolyDrawOrder = "postDraw"

;************************************************
; the subplot call
;************************************************
function subplot(ii,name,pltvar,pltlevs)
begin
    res@cnLevelSelectionMode = "ManualLevels"        ; set the contour levels with the following 3 resources
   if(name .eq. "vicen")then
    res@cnFillPalette        = cryosat2      ; choose colormap
   else
    res@cnFillPalette        = "MPL_rainbow"      ; choose colormap
    ;res@cnFillPalette        = "WhBlGrYeRe"      ; choose colormap
   end if
    res@cnMinLevelValF       = pltlevs(0)            ; set the minimum contour level
    res@cnMaxLevelValF       = pltlevs(1)            ; set the maximum contour level
    res@cnLevelSpacingF      = pltlevs(2)            ; set the interval between contours

   splot = gsn_csm_contour_map_polar(wks,pltvar,res)    ; create the plot

   return(splot)
end
;************************************************
;
;************************************************

  ainame = "aicen"
  ;ailevs = (/0.85, 1.0, 0.01/)
  ailevs = (/0.00, 1.0, 0.05/)

  hiname = "vicen"
  ;hilevs = (/0.55, 5.00, 0.0275/)
  hilevs = (/0.00, 5.0, 0.0275/)

   vname = ainame
    levs = ailevs

  do ii = 0,npanels-1
   icevar = files[ii]->$vname$(:,:,:)
   icevar@long_name = " "
   icevar@units = " "

  ; to plot the sum over the categories
  plotvar = dim_sum_n(icevar,0)
  land = gf->wet
  plotvar = where(land .eq. 0.0, default_fillvalue("double"), plotvar)
  plotvar@lon2d = gf->$icelonname$
  plotvar@lat2d = gf->$icelatname$
  plot(ii) = subplot(ii,vname,plotvar,levs)
   dum(ii) = gsn_add_polyline(wks,plot(ii),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  delete(icevar)
 end do

;************************************************
; panel the plots
;************************************************

  resP                     = True                ; modify the panel plot
                                                 ; new resource added in NCL V6.4.0
  ;resP@gsnPanelLabelBar    = True                ; add common colorbar
  ;resP@lbLabelFontHeightF  = 0.007               ; make labels smaller

 ;resP@gsnPanelYWhiteSpacePercent = 5
 ;resP@gsnPanelXWhiteSpacePercent = 5
  ;resP@gsnPanelMainString = vartitle+"  "+cd_string(times(lstep),"%Y%N%D")
  resP@gsnPanelFigureStrings = cmons

  gsn_panel(wks,plot,(/2,2/),resP)               ; now draw as one plot

exit
end

