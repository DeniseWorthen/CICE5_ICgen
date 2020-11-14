  load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;----------------------------------------------------------------------
begin

        dir = "/scratch2/NCEPDEV/climate/Denise.Worthen/TTout/"
   wgtsfile = "tripole.mx025.Ct.to.mx050.Ct.bilinear.nc"

   cdate = "2012010100"
   icdir = "/scratch1/NCEPDEV/nems/Bin.Li/S2S/FROM_HPSS/"+cdate+"/cpc/"
    g025 = "cice5_model_0.25.res_"+cdate+".nc"
     g05 = "cice5_model_0.50.res_"+cdate+".nc"

   f025 = addfile(icdir+g025,"r")

   vnames = getfilevarnames(f025)
    nvars = dimsizes(vnames)
   dNames = getfiledimnames(f025)
   dSizes = getfiledimsizes(f025)
   print(dNames+"    "+dSizes)

   ; output resolution sizes is known
   nilyr = dSizes(0)
    ncat = dSizes(1)
      nj = 576
      ni = 720

;----------------------------------------------------------------------
; set up the output netcdf file
;----------------------------------------------------------------------

    outfile = dir+g05
    ;outfile = g05
    system("/bin/rm -f " + outfile)    ; remove if exists
    outcdf  = addfile (outfile, "c")  ; open output file

    ; explicitly declare file definition mode. Improve efficiency.
    setfileoption(outcdf,"DefineMode",True)

    ; create global attributes of the file
    fAtt               = True            ; assign file attributes
    fAtt@creation_date = systemfunc ("date")
    fAtt@source_file   = g025
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

    print(varname+"  "+dimsizes(dsizes)+"  "+vtype)
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
   ; skip the iceumask variable
   if(varname .ne. "iceumask")then
     print("regridding variable "+varname)
          var = f025->$varname$
   regrid_var = ESMF_regrid_with_weights(var,dir+wgtsfile,False)
   outcdf->$varname$  = (/regrid_var/)

     delete(var)
     delete(regrid_var)
   end if
  end do

;----------------------------------------------------------------------
; add iceumask, which must be set according the ice concentration
;----------------------------------------------------------------------

 print("setting iceumask")

     aicen = outcdf->aicen
      aice = dim_sum_n(aicen,0)
  iceumask = where(aice.gt.1.0d-11,1.,0.)
  outcdf->iceumask = (/iceumask/)

  delete(outcdf)
exit
end
