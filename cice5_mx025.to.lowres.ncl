  load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;----------------------------------------------------------------------
begin

; define some constants
saltmax = 3.2d0
   nsal = 0.407d0
   msal = 0.573d0
     pi = atan(1.0d0)*4.0d0
   rhoi =  917.d0
   rhos =  330.d0
 cp_ice = 2106.d0
 cp_ocn = 4218.d0
 Lfresh = 3.34d5
   puny = 1.0d-3

; set options for poisson_grid_fill
  guess     = 1                ; use zonal means
  is_cyclic = True             ; cyclic [global]
  nscan     = 1500             ; usually much less than this
  eps       = 1.e-2            ; variable dependent
  relc      = 0.6              ; relaxation coefficient
  opt       = 0                ; not used

;----------------------------------------------------------------------
;
;----------------------------------------------------------------------

  ; specify grid and weights directory
   grdsrc = "/scratch2/NCEPDEV/climate/Denise.Worthen/grids-20210727/"
  ; specify IC source directory
  ; CPC reanalysis
  ;  icsrc = "/scratch2/NCEPDEV/climate/climpara/S2S/IC/CPC/"
    icsrc = "/scratch1/NCEPDEV/stmp4/Denise.Worthen/CPLD_V16_ICS/"
  ;  icsrc = "/scratch2/NCEPDEV/stmp1/Denise.Worthen/ForXingren/"
  ; RT IC (20161003 only)
  ;  icsrc = "/scratch1/NCEPDEV/nems/emc.nemspara/RT/NEMSfv3gfs/input-data-20201201/CICE_IC/"

  ; specify output directory
  ; dirout = "/scratch2/NCEPDEV/climate/Denise.Worthen/NEWIC/"
   dirout = "/scratch1/NCEPDEV/stmp4/Denise.Worthen/CPLD_V16_ICS/"
  ; dirout = "/scratch2/NCEPDEV/stmp1/Denise.Worthen/ForXingren/"
  ; specify a date
  ;  cdate = "2016100300"
  ;  cdate = "2011070100"
  ;  cdate = "2012091500"
    cdate = "20210322"
  ;  cdate = "00000000"

  ; specify the output resolution identifiers and grid size
  ;output = (/"0.50", "mx050" /)
  ;outsiz = (/576,720/)
  ;output = (/"1.00", "mx100" /)
  ;outsiz = (/320,360/)
  ;output = (/"0.72", "mx072" /)
  ;outsiz = (/381,500/)
  output = (/"4.00", "mx400" /)
  outsiz = (/35,72/)
  
   ores = output(0)
   cres = output(1)
     nj = outsiz(0)
     ni = outsiz(1)

  ; use filling in extra ocean regions of 3deg grid
  if(output(0) .eq. "3.00")then
   ;use_poisson = "true"
   use_poisson = "false"
  else
   use_poisson = "false"
  end if

;----------------------------------------------------------------------
;
;----------------------------------------------------------------------

   ; fout = "cice_model_"+ores+".res_"+cdate+".nc"
    fout = "cice_model_"+ores+".cpc.res_"+cdate+".nc"
 
   ; if using CPC reanalysis
   ic025 = "cice_model_0.25.cpc.res_"+cdate+".nc"
   ; f025 = addfile(icsrc+cdate+"/ice/025/"+ic025,"r")
    f025 = addfile(icsrc+ic025,"r")
   ; if RT IC
   ;f025 = addfile(icsrc+"025/"+ic025,"r")
   ;ic025 = "cice_model_0.25.res_"+cdate+".nc"
   ; f025 = addfile(icsrc+ic025,"r")

   vnames = getfilevarnames(f025)
    nvars = dimsizes(vnames)
   dNames = getfiledimnames(f025)
   dSizes = getfiledimsizes(f025)
   print(dNames+"    "+dSizes)

   ; output resolution sizes is known
   nilyr = dSizes(0)
    ncat = dSizes(1)
   nslyr = 1

;----------------------------------------------------------------------
;
;----------------------------------------------------------------------

    ; get the known land for output grid
         gr = addfile(grdsrc+"grid_cice_NEMS_"+cres+".nc","r")
      land = where(gr->kmt .eq. 0, 0.0d0, 1.0d0)
    land3d = conform_dims((/ncat,nj,ni/),land,(/1,2/))
    delete(gr)

   ; the mask contains 1's on land and 0's at valid points.
         gr = addfile(grdsrc+"grid_cice_NEMS_mx025.nc","r")
    mask025 = where(gr->kmt .eq. 0, 1.0, 0.0)
    delete(gr)
   ; regrid the mask to obtain the interpolation masks.
   ; the mask025 contain 1's on land and 0's at valid points.
   ; when remapped, any mask value > 0 identifies land values that
   ; have crept into the field. remapped model fields are then
   ; masked with this interpolation mask

   ; specify the tripole-to-tripole weights file
    wgtsfile = "tripole.mx025.Ct.to."+cres+".Ct.neareststod.nc"
      rgmask = ESMF_regrid_with_weights(mask025, grdsrc+wgtsfile,False)
    rgmask3d = conform_dims((/ncat,nj,ni/),rgmask,(/1,2/))

    ;printVarSummary(rgmask)
    ;printVarSummary(rgmask3d)

    ; use rgmask3d >0 and <1 to mask areas which are land-contaminated
    ; this makes all land areas missing in order to use the poisson_grid_fill
    ; routine
      rgmask = where(rgmask   .gt. 0.0, default_fillvalue(typeof(rgmask)),   1.0)
    rgmask3d = where(rgmask3d .gt. 0.0, default_fillvalue(typeof(rgmask3d)), 1.0)

    ;testfile = "mask_bilinear.nc"
    ;system("/bin/rm -f "+testfile)
    ; create
    ;testcdf = addfile(testfile,"c")
    ;testcdf->rgmask = rgmask
    ;testcdf->rgmask3d = rgmask3d
    ;testcdf->mask025 = mask025
    ;testcdf->land = land
    ; close
    ;delete(testcdf)

;----------------------------------------------------------------------
; define required variables
;----------------------------------------------------------------------

   ; a reasonable salinity profile
   salinz = new((/nilyr/),double)
   do l = 0,nilyr-1
           zn = (int2dble(l+1)-0.5d0)/int2dble(nilyr)
    salinz(l) = (saltmax/2.d0)*(1.d0-cos(pi*zn^(nsal/(msal+zn))))
   end do
   Tmltz = salinz / (-18.48 + (0.01848*salinz))
   print(Tmltz)
    sice = new((/nilyr,ncat,nj,ni/),double)
    qice = new((/nilyr,ncat,nj,ni/),double)
    qsno = new((/nslyr,ncat,nj,ni/),double)
     Tin = new((/nilyr,ncat,nj,ni/),double)

   ; sis2 is run with constant ice salinity, so substitue profile salinz
   do l = 0,nilyr-1
    sice(l,:,:,:) = salinz(l)
   end do
   
;----------------------------------------------------------------------
; set up the output netcdf file
;----------------------------------------------------------------------

    outfile = dirout+fout
    print(outfile)
    ;outfile = fout
    system("/bin/rm -f " + outfile)    ; remove if exists
    outcdf  = addfile (outfile, "c")  ; open output file

    ; explicitly declare file definition mode. Improve efficiency.
    setfileoption(outcdf,"DefineMode",True)

    ; create global attributes of the file
    fAtt               = True            ; assign file attributes
    fAtt@creation_date = systemfunc ("date")
    fAtt@source_file   = ic025
    fileattdef( outcdf, fAtt )           ; copy file attributes

    ; predefine the coordinate variables and their dimensionality
    dimNames = dNames
    dimSizes = (/nilyr, ncat, nj, ni/)
    dimUnlim = (/False,  False,  False, False/)
    filedimdef(outcdf,dimNames,dimSizes,dimUnlim)

    ; predefine each variable in the output file
    do nv = 0,nvars-1
      varname = vnames(nv)
       vtype := getfilevartypes(f025,varname)
      vatts  := getfilevaratts(f025,varname)       ; The := is important because these
      dnames := getfilevardimnames(f025,varname)   ; variables may be a different size
      dsizes := getfilevardimsizes(f025,varname)   ; the next time through the loop
    
    print(nv+"  "+varname+"  "+dimsizes(dsizes)+"  "+vtype)
    if(dimsizes(dsizes) .eq. 4)then
     odims = (/"nilyr", "ncat", "nj", "ni"/)
     filevardef(outcdf, varname, vtype, odims)
    end if
    if(dimsizes(dsizes) .eq. 3)then
     odims = (/"ncat", "nj", "ni"/)
     filevardef(outcdf, varname, vtype, odims)
    end if
    if(dimsizes(dsizes) .eq. 2)then
     odims = (/"nj", "ni"/)
     filevardef(outcdf, varname, vtype, odims)
    end if
     delete(odims)
    end do

    ; explicitly exit file definition mode.
    setfileoption(outcdf,"DefineMode",False)

;----------------------------------------------------------------------
; regrid each variable and place in output cdf
; this ignores any difference between angle of the two grids; the ice
; velocities could be inaccurate
;----------------------------------------------------------------------

  do nv = 0,nvars-1
      varname = vnames(nv)
        sname = str_get_cols(varname,0,3)
   ; skip the iceumask variable and layer variables
   if(varname .ne. "iceumask" .and.  \
        sname .ne.     "qice" .and.  \
        sname .ne.     "qsno" .and.  \
        sname .ne.     "Tsfc" .and.  \
        sname .ne.     "sice")then
     print("regridding variable "+varname)
          var = f025->$varname$
   rgtt = ESMF_regrid_with_weights(var,grdsrc+wgtsfile,False)
   ndims = dimsizes(dimsizes(rgtt))
   ;print(varname +"  "+ndims)

   ; mask out land contaminated areas. this skips Tin variable
   if(ndims .eq. 3)then
    if(use_poisson .eq. "true")then
     rgtt = where(ismissing(rgmask3d), rgtt@_FillValue, rgtt)
     ; fill missing areas
      poisson_grid_fill( rgtt, is_cyclic, guess, nscan, eps, relc, opt)
    end if
     ; apply known land mask
     rgtt = where(land3d .eq. 0.0, 0.0, rgtt)
   end if

   if(ndims .eq. 2)then
    if(use_poisson .eq. "true")then
      rgtt = where(ismissing(rgmask),   rgtt@_FillValue, rgtt)
     ; fill missing areas
      poisson_grid_fill( rgtt, is_cyclic, guess, nscan, eps, relc, opt)
    end if
     ; apply known land mask
     rgtt = where(land .eq. 0.0, 0.0, rgtt)
   end if

   outcdf->$varname$  = (/rgtt/)

     delete(var)
     delete(rgtt)
   end if
  end do

  ;retrieve the regridded ice concentration for masking
  ;   aicen = outcdf->aicen
  ;    aice = dim_sum_n(aicen,0)
  ;    aice = where(aice .lt. puny, 0.0d0, aice)

;----------------------------------------------------------------------
; retrieve the layer variables
;----------------------------------------------------------------------

  ; regrid the surface temperature field
  rgtt = ESMF_regrid_with_weights(f025->Tsfcn,grdsrc+wgtsfile,False)
  do k=0,nslyr-1
   qsno(k,:,:,:) = -rhos*(Lfresh - cp_ice*rgtt(:,:,:))
  end do

  do nv = 0,nvars-1
   varname = vnames(nv)
   if(str_get_cols(varname,0,3) .eq. "qice")then
    do kk = 0,nilyr-1
     qname = "qice"+sprinti("%0.3i",kk+1)
     qice(kk,:,:,:) = ESMF_regrid_with_weights(f025->$qname$,grdsrc+wgtsfile,False)
      Tin(kk,:,:,:) = (qice(kk,:,:,:) + rhoi*Lfresh)/(rhoi*cp_ice)
      Tin(kk,:,:,:) = where(Tin(kk,:,:,:) .gt. Tmltz(kk), Tmltz(kk), Tin(kk,:,:,:))
     qice(kk,:,:,:) = rhoi*cp_ice*Tin(kk,:,:,:) - rhoi*Lfresh
    end do
   end if
   ;if(str_get_cols(varname,0,3) .eq. "qsno")then
   ; do kk = 0,nslyr-1
   ;  qname = "qsno"+sprinti("%0.3i",kk+1)
   ;  qsno(kk,:,:,:) = ESMF_regrid_with_weights(f025->$qname$,grdsrc+wgtsfile,False)
   ; end do
   ;end if
  end do
 
  do k = 0,ncat-1
   do l = 0,nslyr-1
   ; qsno(l,k,:,:) = where(aice(:,:) .eq. 0.0d0, 0.0d0, qsno(l,k,:,:))
   end do
   do l = 0,nilyr-1
   ; qice(l,k,:,:) = where(aice(:,:) .eq. 0.0d0, 0.0d0, qice(l,k,:,:))
   ; sice(l,k,:,:) = where(aice(:,:) .eq. 0.0d0, 0.0d0, sice(l,k,:,:))
   end do
  end do

  ; write the qice and sice
  outcdf->qice001 = (/qice(0,:,:,:)/)
  outcdf->qice002 = (/qice(1,:,:,:)/)
  outcdf->qice003 = (/qice(2,:,:,:)/)
  outcdf->qice004 = (/qice(3,:,:,:)/)
  outcdf->qice005 = (/qice(4,:,:,:)/)
  outcdf->qice006 = (/qice(5,:,:,:)/)
  outcdf->qice007 = (/qice(6,:,:,:)/)
  outcdf->sice001 = (/sice(0,:,:,:)/)
  outcdf->sice002 = (/sice(1,:,:,:)/)
  outcdf->sice003 = (/sice(2,:,:,:)/)
  outcdf->sice004 = (/sice(3,:,:,:)/)
  outcdf->sice005 = (/sice(4,:,:,:)/)
  outcdf->sice006 = (/sice(5,:,:,:)/)
  outcdf->sice007 = (/sice(6,:,:,:)/)
  ; write the Tsfcn and qsno
  outcdf->Tsfcn   = (/rgtt/)
  outcdf->qsno001 = (/qsno(0,:,:,:)/)
  delete(rgtt)

;----------------------------------------------------------------------
; add iceumask, which must be set according the ice concentration
;----------------------------------------------------------------------

 print("setting iceumask")

     aicen = outcdf->aicen
      aice = dim_sum_n(aicen,0)
  iceumask = where(aice .ge. puny, 1., 0.)
  outcdf->iceumask = (/iceumask/)

  delete(outcdf)
exit
end
