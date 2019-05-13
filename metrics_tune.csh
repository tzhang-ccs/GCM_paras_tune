#!/bin/csh
set experimenter = "'zhangtao'"
set NOTE = "'E3SM_low_opt_TWP'"

set hist_prefix = "low_ne30.cam.h1.2009-07"
set hist_path = "/global/cscratch1/sd/zhangtao/acme_scratch/cori-knl/low_ne30/run/regridded/"
set hist_tune = "/global/cscratch1/sd/zhangtao/acme_scratch/cori-knl/low_ne30/run/hist/"
set hist_run = "/global/cscratch1/sd/zhangtao/acme_scratch/cori-knl/low_ne30/run/"
set uid=` date +"%y-%m-%d-%H-%M"`
setenv climo_path  climo/

#computing climate average
set file_list = ''
foreach day (`seq -w 1 31`)
    set file_list = "$file_list $hist_path/$hist_prefix-$day-00000.nc"
end

set var_list = "PRECT CLDHGH CLDLOW CLDMED LWCF SWCF"
foreach var_name ($var_list)
	ncea -O -v $var_name $file_list ${climo_path}/model/${var_name}.nc
    #ncks -O -v $var_name $hist_tune/low_ne30.cam.h1.2009-07-mean-19-02-14-20-44.nc ${climo_path}/model/${var_name}.nc
end

#save to hist tune
#cp $hist_run/atm_in $hist_tune/atm_in_$uid
#ncea -O $file_list $hist_tune/$hist_prefix-mean-$uid.nc

##computing metrics
ncl calc_metrics_tune.ncl
#ncl set_pdf.ncl
#
#
set inst_str = "$experimenter, "

set inst_str = "$inst_str `grep -w zmconv_dmpdz  ../../atm_in  | cut -d = -f 2`,"
set inst_str = "$inst_str `grep -w zmconv_c0_ocn ../../atm_in  | cut -d = -f 2`,"
set inst_str = "$inst_str `grep -w zmconv_c0_lnd ../../atm_in  | cut -d = -f 2`,"
set inst_str = "$inst_str `grep -w cldfrc_dp1    ../../atm_in  | cut -d = -f 2`,"
set inst_str = "$inst_str `grep -w clubb_c1      ../../atm_in  | cut -d = -f 2`,"
set inst_str = "$inst_str `grep -w clubb_c8      ../../atm_in  | cut -d = -f 2`,"
set inst_str = "$inst_str `grep -w clubb_c14     ../../atm_in  | cut -d = -f 2`,"
set inst_str = "$inst_str `grep -w ice_sed_ai    ../../atm_in  | cut -d = -f 2`,"

foreach i (`cat mcpi_ratio`)
    set inst_str = "$inst_str `echo $i,`"
end

set inst_str = "$inst_str `cat mcpi`",
set inst_str = "$inst_str  $NOTE"


cat >! metrics.sql << EOF
USE uq_e3sm;
INSERT INTO e3sm_tune_MOO (experimenter, dmpdz, c0_ocn, c0_lnd, dp1, c1, c8, c14, ai, 
PRECT, CLDHGH, CLDLOW, CLDMED, LWCF, SWCF, MCPI, NOTE) VALUES
($inst_str);
EOF

mysql uq_e3sm  -u  uq_e3sm_admin -h nerscdb04.nersc.gov -p3Ii3i3i2fdd_2s25j333jjdd < metrics.sql
