begin
	model_path = "climo/model/"
	obs_path   = "/global/homes/z/zhangtao/ACME/GCM_paras_tune/obs/CAPT/"
	var_names = (/"PRECT","CLDHGH","CLDLOW","CLDMED","LWCF","SWCF"/)	

	deg2rad = atan(1.0)/45.0

	num_var = dimsizes(var_names)
	e_score = new((/num_var/),float)
    e_score_str = new((/num_var/), string)
	e_score@_FillValue = 1.0e33

    lat_min = -26 
    lat_max = 6 
    lon_min = 107 
    lon_max = 159 

	do k = 0, num_var - 1
		obs_file   = obs_path + var_names(k) + "/"+var_names(k)+"_Jul_mean.nc"
		model_file = model_path + var_names(k) + ".nc"
		fp_obs     = addfile(obs_file, "r")
		fp_model   = addfile(model_file, "r")
		var_obs    = fp_obs->$var_names(k)$
		var_model  = fp_model->$var_names(k)$(0,:,:)

        lat_obs    = fp_obs->lat           
        lon_obs    = fp_obs->lon          
        lat_mod    = fp_model->lat
        lon_mod    = fp_model->lon

        var_obs!0 = "lat"
        var_obs!1 = "lon"
        var_obs&lat = lat_obs
        var_obs&lon = lon_obs

        if ((var_names(k) .eq. "PRECT"))
            var_model = var_model * 1000.0*86400.0
        end if
        if (var_names(k) .eq. "Q850")then
            var_model = var_model * 1000.0
        end if

        var_mod_remap = area_conserve_remap_Wrap(lon_mod, lat_mod, var_model, lon_obs, lat_obs, False)

        var_obs_reg := var_obs({lat_min:lat_max},{lon_min:lon_max})
        var_mod_reg := var_model({lat_min:lat_max},{lon_min:lon_max})
        var_mod_remap_reg := var_mod_remap({lat_min:lat_max},{lon_min:lon_max})


		sd_mod = stddev(var_mod_reg)
		sd_obs = stddev(var_obs_reg)

		R = pattern_cor(var_mod_remap_reg, var_obs_reg, 1.0, 0)
        R0 = 1

        e_score(k) = log(((sd_mod/sd_obs + sd_obs/sd_mod)^2*16)/(4*(1+R)^4))
	    e_score_str(k) = sprintf("%10.10f", e_score(k))
		delete(var_model)
        delete(var_obs)
        delete(var_mod_remap)
        delete(lat_obs)
        delete(lon_obs)
        delete(lat_mod)
        delete(lon_mod)
	end do
    ;print(e_score_str)
	cntl_score = asciiread("cntl_score", -1, "float")
	mcpi_ratio = e_score/cntl_score
    print(mcpi_ratio)
	asciiwrite("mcpi_ratio", mcpi_ratio)
	mcpi = dim_avg_n(mcpi_ratio,0)
	asciiwrite("mcpi", mcpi)
    ;print(e_score_str)
end
