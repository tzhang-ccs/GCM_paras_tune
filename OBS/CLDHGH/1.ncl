begin
fin = addfile("CLDHGH_Jul_mean.nc", "r")
var = fin->CLDHGH(::-1,:)

fout = addfile("tmp.nc", "c")
fout->CLDHGH = var
end
