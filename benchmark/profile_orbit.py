from benchmark.benchmark_orbit import init_python_orbit

from benchmark.benchmark_orbit import time_orbit_translating
from benchmark.benchmark_orbit import time_orbit_sites_permutations




orbit_pair, orbit_triplet = init_python_orbit(100)

py_timing = time_orbit_translating(orbit_pair)
py_timing = time_orbit_sites_permutations(orbit_pair)


