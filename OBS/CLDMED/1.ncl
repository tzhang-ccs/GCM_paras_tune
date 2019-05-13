begin
fin = addfile("CLDMED_Jul_mean.nc", "r")
var = fin->CLDMED(:,:)
var = var * 100.0

fout = addfile("tmp.nc", "c")
fout->CLDMED = var
end
