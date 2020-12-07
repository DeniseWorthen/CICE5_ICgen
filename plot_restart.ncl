;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;************************************************
begin

  ; hera
  RT1 = "/scratch2/NCEPDEV/stmp1/Denise.Worthen/"
  RT2 = "/scratch1/NCEPDEV/stmp2/Denise.Worthen/"
  RT3 = "/scratch2/NCEPDEV/stmp3/Denise.Worthen/"
  RT4 = "/scratch1/NCEPDEV/stmp4/Denise.Worthen/"

  ;hemi="SH"

  hemi="NH"
  minlat = 60.0

  ncat = 5

;************************************************
;
;************************************************

   cdate = "2016100300"
  grdsrc = "/scratch2/NCEPDEV/climate/Denise.Worthen/GRIDS/"
   icsrc = "/scratch2/NCEPDEV/climate/Denise.Worthen/NEWIC/"
   ;icsrc = "/scratch1/NCEPDEV/nems/emc.nemspara/RT/NEMSfv3gfs/input-data-20201201/CICE_IC/"

  ; in grdsrc
  latname = "latCt"
  lonname = "lonCt"
  lndname = "wet"

  ; fields to plot (these are all ncat)
  iceflds = (/"aicen",  \
              "vicen",  \
              "vsnon",  \
              "Tsfcn"/)
  nflds = dimsizes(iceflds)

       plotlevs = new((/nflds,3/), float)
  plotlevs(0,:) = (/  0.0,  1.00,  0.10/)
  plotlevs(1,:) = (/  0.0,  3.00,  0.10/)
  plotlevs(2,:) = (/  0.0,  0.50,  0.05/)
  plotlevs(3,:) = (/-20.0,  0.00,  1.00/)

;************************************************
;
;************************************************

   ; arrays of characterstrings to set locations and filenames
   ;  fsrcs => tripole grid file
   ; rtsrcs => subdirs for the rt ICs
   ;  csrcs => cice restart file
    fsrcs = (/ "mx100",  "mx050",  "mx025", "mx072"/)
   rtsrcs = (/  "100/",   "050/",   "025/",  "072/"/)
    csrcs = (/  "1.00",   "0.50",   "0.25",  "0.72"/)
   ; seam longitude for each res
     lons = (/    60.0,     60.0,     60.0,   73.75/)

  ;ires = 0 ; 1deg
  ;ires = 1 ; 1/2deg
  ;ires = 2 ; 1/4deg
  ires = 3 ; 072deg

  grdfile = "tripole."+fsrcs(ires)+".nc"
  icefile = "cice_model_"+csrcs(ires)+".res_"+cdate+".nc"
  print("using grid file "+grdfile)
  print("using IC file "+icefile)

    gf = addfile(grdsrc+grdfile,"r")
  ;icef = addfile(icsrc+rtsrcs(ires)+icefile,"r")
  icef = addfile(icsrc+icefile,"r")

  ; retrieve lat,lon,mask
  lon2d = gf->$lonname$
  lat2d = gf->$latname$
  lsmsk = gf->$lndname$
  delete(gf)

  wkstitle = cdate+"_"+fsrcs(ires)

;************************************************
; great circle path from (lon1,lat1) to (lon2,lat2)
;************************************************
  
  lon1 =  lons(ires)
  lon2 =  lon1+180.0
  lat1 =  60.00
  lat2 =  60.00

  print("Great Circle "+lon1+" "+lon2)

;************************************************
; define the map
;************************************************

   wkstype="X11"
  ; wkstype="png"
  wkstype@wkWidth=800
  wkstype@wkHeight=800
  wks   = gsn_open_wks (wkstype,wkstitle)

  ; panels
  npanels = nflds
  plot = new(npanels,graphic)
   dum = new(npanels,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@mpFillOn    = False
  res@mpOutlineOn = False
  ;res@cnFillMode           = "MeshFill"
  res@cnFillMode           = "RasterFill"
  res@cnFillDrawOrder      = "PreDraw"             ; draw contours first

  res@cnFillOn           = True                  ; color fill
  res@cnLinesOn          = False                 ; no contour lines
  res@lbLabelBarOn       = True                  ; turn on individual colorbars
  res@lbLabelFontHeightF = 0.015

  res@gsnPolar   =  hemi                         ; specify the hemisphere
  if(hemi .eq. "NH")then
   res@mpMinLatF    =  minlat                    ; minimum lat to plot
   res@mpCenterLonF = -45
  end if
  if(hemi .eq. "SH")then
   res@mpMaxLatF    = -50                        ; maximum lat to plot
   res@mpCenterLonF = 0.0
  end if

  res@cnMissingValPerimOn     = False            ; turn on the missing value perimeter
  res@cnMissingValFillColor   = "gray"           ; color fill missing areas grey

  res@cnFillPalette           = "amwg256"        ; choose colormap

;************************************************
; Resources for the polyline
;************************************************

  res@mpGreatCircleLinesOn = True
  pres                  = True           ; polyline mods desired
  pres@gsLineThicknessF = 2.0            ; line thickness
  pres@gsLineColor      = "black"        ; color of lines
  pres@gsLineDashPattern =  1            ; dashes

;************************************************
; the subplot call
;************************************************
function subplot(ii,fname,name,pltvar,pltlevs)
begin
    res@cnLevelSelectionMode = "ManualLevels"        ; set the contour levels with the following 3 resources
    res@cnMinLevelValF       = pltlevs(0)            ; set the minimum contour level
    res@cnMaxLevelValF       = pltlevs(1)            ; set the maximum contour level
    res@cnLevelSpacingF      = pltlevs(2)            ; set the interval between contours

   ;res@gsnLeftString = name
   res@gsnCenterString = name

   splot = gsn_csm_contour_map_polar(wks,pltvar,res)    ; create the plot
   return(splot)
end
;************************************************
;
;************************************************
 
  ii = 3

  do ii = 0,nflds-1

   vname = iceflds(ii)
    levs = plotlevs(ii,:)
   print("plotting variable "+vname)

   ; retrieve category vars
   if(vname .eq. "Tsfcn")then
      tmp1 = icef->aicen(:,:,:)
      tmp2 = icef->$vname$(:,:,:)
       var = dim_sum_n_Wrap(tmp1*tmp2,0)
      delete([/tmp1,tmp2/])
   else
      var = dim_sum_n_Wrap(icef->$vname$(:,:,:),0)
   end if
   icevar = mask(var,lsmsk,1)  ; return values of var where lsmsk=1
   icevar@lon2d = lon2d
   icevar@lat2d = lat2d

   plot(ii) = subplot(ii," ",vname+" "+cdate,icevar,levs)
    dum(ii) = gsn_add_polyline(wks,plot(ii),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
   delete([/icevar,var/])
  end do

;************************************************
; panel the plots
;************************************************

  resP                     = True                ; modify the panel plot
                                                 ; new resource added in NCL V6.4.0
 ;resP@gsnPanelYWhiteSpacePercent = 5
 ;resP@gsnPanelXWhiteSpacePercent = 5

  ;gsn_panel(wks,plot,(/1,1/),resP)               ; now draw as one plot
  gsn_panel(wks,plot,(/2,2/),resP)               ; now draw as one plot

;************************************************
;
;************************************************

  if(wkstype .eq. "png")then
   imcmd = "convert -trim "+wkstitle+".png "+wkstitle+".jpg"
   ;print(imcmd)
   system(imcmd)
   imcmd = "convert -bordercolor white -border 0x10 "+wkstitle+".jpg "+wkstitle+".jpg"
   ;print(imcmd)
   system(imcmd)
  end if

exit
end