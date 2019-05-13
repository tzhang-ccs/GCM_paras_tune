nc_list=""
for i in `seq -f "%02g" 1 5`
do
	nc_list=$nc_list" SWCF_200907$i.nc"
done

echo $nc_list
ncea $nc_list  SWCF_200907_5day_mean.nc
