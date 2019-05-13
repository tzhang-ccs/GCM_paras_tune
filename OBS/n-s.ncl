load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"

begin
	file_in="PRECT/prect_2009_Jan.nc"
    file_out="PRECT_Jan_mean.nc"
	fp_in=addfile(file_in,"r")
	var=fp_in->PRECT(::-1, :)
	lat=fp_in->lat(::-1)
	lon=fp_in->lon

	fp_out = addfile(file_out, "c")
	fp_out->PRECT = var
	fp_out->lat = lat
	fp_out->lon = lon
	 
end
