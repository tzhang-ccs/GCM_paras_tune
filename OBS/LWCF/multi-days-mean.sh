nc_list=""
for i in `seq -w 1 14`
do
	nc_list=$nc_list" LWCF_200907$i.nc"
done

ncea $nc_list  LWCF_200907_14day_mean.nc
