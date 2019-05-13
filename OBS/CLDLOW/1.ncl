begin
fin = addfile("CLDLOW_Jul_mean.nc", "r")
var = fin->CLDLOW(:,:)
var = var * 100.0

fout = addfile("tmp.nc", "c")
fout->CLDLOW = var
end
