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