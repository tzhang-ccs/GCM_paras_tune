#!/bin/bash
paras="zmconv_dmpdz:zmconv_c0_ocn:zmconv_c0_lnd:cldfrc_dp1:clubb_c1:clubb_c8:clubb_c14:ice_sed_ai"
sample_dir="/global/homes/z/zhangtao/ACME/GCM_paras_tune/algorithms/downhill_simplex/model/"
run_dir="/global/cscratch1/sd/zhangtao/acme_scratch/cori-knl/low_ne30/run/"
metrics_path="$run_dir/metrics"
sampling_file="sampling_data"

#run
cd  $run_dir
paras_num=`echo $paras | awk -F ':' '{print NF}'`
for i in `seq 1 $paras_num`
do
	para=`echo $paras |cut -d : -f $i`
	var_line=`sed -n "1p" $sample_dir/$sampling_file`
    para_val="$para=`echo $var_line |cut -d ' ' -f $i`"
    echo $para_val
    sed -i "/\<$para\>/c \\  $para_val" $run_dir/atm_in
done


./run_e3sm.sh
wait

echo "e3sm has finished running!"
./remap_2d.sh low_ne30.cam.h1.2009-07*nc
./remap_2d.sh low_ne30.cam.h0*

#diag
cd $metrics_path/TWP
./metrics_tune.csh > metrics.log
cd $metrics_path/Global
./metrics_tune.csh > metrics.log
