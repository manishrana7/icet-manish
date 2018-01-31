Some notes about the benchmark folder:
* Files beginning with benchmark will print out timings of different functions  
  and if possible compare C++ and python versions  
* Files beginning with profile are similar to the benchmark scripts but are  
  designed to be profiled. As such they only use either the C++ or Python version.  

When making changes to core classes run the benchmark scripts before and after the  
change to make sure nothing gets broken in terms of performance.  
If something gets broken then use the profile scripts to pin point the problem.  


Read about basic profiling here:
* https://docs.python.org/3/library/profile.html#introduction-to-the-profilers


Example
-------
Profile LatticeSite by running

```
python3 -m cProfile benchmark/profile_lattice_site.py 
```

One can easily sort the output by cumulative time with sort:

```
python3 -m cProfile -s cumtime benchmark/profile_lattice_site.py
```
Which results in something like:

```
Timing for sorting: 593.750000 micro sec
Timing for hash: 656.250000 micro sec
Timing for index lookup: 109.375000 micro sec
Timing for offset lookup: 109.375000 micro sec
         2121324 function calls (2114707 primitive calls) in 3.621 seconds

   Ordered by: cumulative time

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
    734/1    0.060    0.000    3.623    3.623 {built-in method builtins.exec}
        1    0.000    0.000    3.623    3.623 profile_lattice_site.py:1(<module>)
       23    0.001    0.000    2.173    0.094 __init__.py:1(<module>)
    661/1    0.006    0.000    1.501    1.501 <frozen importlib._bootstrap>:966(_find_and_load)
    661/1    0.006    0.000    1.501    1.501 <frozen importlib._bootstrap>:939(_find_and_load_unlocked)
    815/2    0.002    0.000    1.499    0.750 <frozen importlib._bootstrap>:214(_call_with_frames_removed)
    415/1    0.001    0.000    1.499    1.499 {built-in method builtins.__import__}
    592/3    0.006    0.000    1.499    0.500 <frozen importlib._bootstrap>:659(_load_unlocked)
        2    0.000    0.000    1.498    0.749 __init__.py:2(<module>)
    489/2    0.003    0.000    1.498    0.749 <frozen importlib._bootstrap_external>:659(exec_module)
        1    0.003    0.003    1.195    1.195 profile_lattice_site.py:21(time_sorting)
     2000    0.153    0.000    1.192    0.001 {method 'sort' of 'list' objects}
   198000    0.604    0.000    1.039    0.000 lattice_site.py:63(__lt__)
2937/1955    0.007    0.000    0.866    0.000 <frozen importlib._bootstrap>:996(_handle_fromlist)
```


Benchmarks
==========

Orbit
-----
This benchmark typical orbit operations that are done
in `icet`. 
Orbit translating is for creating
a new orbit where all lattice sites have had their
unitcell offsets translated. This is useful to get all
lattice sites of a supercell if your orbit is for the
primitive cell.

orbit permutation is a property that returns
the permutated equivalent sites so that they
are in an identical order to the representative sites.
This is useful when counting element combinations on
an orbit. E.g, if we are interested in how many ABB elements
there are compared to BAB then we would need the permutated
sites.

A typical output from benchmark/benchmark_orbit.py is

```
Time for python pair orbit translating: 0.00318750s
Time for C++ pair orbit translating: 0.00001563s
Cpp speedup 204.000

Time for python pair orbit translating: 0.00484375s
Time for C++ pair orbit translating: 0.00000000s
Cpp speedup inf

Time for python pair orbit permutation: 0.00010312s
Time for C++ pair orbit permutation: 0.00014062s
Cpp speedup 0.733

Time for python triplet orbit permutation:
     0.00010938s
Time for C++ triplet orbit permutation: 0.00021875s
Cpp speedup 0.500
```
This test was run on a Intel(R) Core(TM) i7-6500U CPU @ 2.50GHz




Profile
======

Orbit
-----

Running
```
 python3 -m cProfile -s cumtime benchmark/profile_orbit.py
 ```
 will result in something like this:

 ```
          18279299 function calls (16467651 primitive calls) in 29.205 seconds

   Ordered by: cumulative time

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
    738/1    0.054    0.000   29.207   29.207 {built-in method builtins.exec}
        1    0.000    0.000   29.207   29.207 profile_orbit.py:1(<module>)
        1    0.078    0.078   24.173   24.173 benchmark_orbit.py:98(time_orbit_translating)
     1000    0.636    0.001   24.095    0.024 orbit.py:130(__add__)
1701000/1000    8.138    0.000   23.171    0.023 copy.py:137(deepcopy)
101000/1000    0.710    0.000   23.158    0.023 copy.py:214(_deepcopy_list)
   200000    2.077    0.000   17.394    0.000 copy.py:269(_reconstruct)
   200000    1.262    0.000    8.114    0.000 copy.py:239(_deepcopy_dict)
        1    0.044    0.044    3.553    3.553 benchmark_orbit.py:116(time_orbit_sites_permutations)
     5000    0.023    0.000    3.509    0.001 orbit.py:205(permutated_sites)
     5000    0.390    0.000    3.456    0.001 orbit.py:211(<listcomp>)
   500000    1.380    0.000    3.066    0.000 orbit.py:213(get_permutated_sites)
   200000    0.873    0.000    2.564    0.000 copy.py:222(_deepcopy_tuple)
   701000    1.441    0.000    2.221    0.000 copy.py:253(_keep_alive)
  3603912    2.184    0.000    2.184    0.000 {method 'get' of 'dict' objects}
       23    0.001    0.000    2.097    0.091 __init__.py:1(<module>)
  3104339    1.660    0.000    1.660    0.000 {built-in method builtins.id}
   200000    0.365    0.000    1.584    0.000 copy.py:223(<listcomp>)
   ```
showing that deepcopy and _reconstruct is what is so slow with the
constructing a translated orbit.

