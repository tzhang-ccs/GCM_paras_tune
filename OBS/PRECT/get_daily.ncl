begin
date_start = 20090701
date_len   = 31

fid   = addfile("GPCP_1DD_v1.2_199610-201507.nc", "r")
date  = fid->date
prect = fid->PREC
lat   = fid->lat
lon   = fid->lon

index = ind(date.eq.date_start)
print(date(index:index+date_len-1))  

prect_part = prect(index:index+date_len-1,::-1,:)

prect_part_mean = dim_avg_n_Wrap(prect_part, 0)
;prect_part_mean!0 = "lat"
;prect_part_mean!1 = "lon"
;prect_part_mean&lat = lat
;prect_part_mean&lon = lon


;fid_out = addfile("PRECT_200907_5day_mean.nc", "c")
;fid_out->PRECT =  prect_part_mean

fid_out = addfile("PRECT_2009_07_daily.nc", "c")
fid_out->PRECT = prect_part

end
