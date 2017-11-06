from icetdev.latticeSite import LatticeSite

index = 1
offset = [0., 0., 0.,]
latnbr = LatticeSite(index, offset)


#Test that it is hashable

latnbr_map = {}

latnbr_map[latnbr] = 1

assert latnbr in latnbr_map

print(latnbr_map[latnbr])
