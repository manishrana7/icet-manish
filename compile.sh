#rm -r build
#mkdir build
cd build
cmake ..
make -j 4
cd ..
#export PYTHONPATH=$PWD/build/src/:$PYTHONPATH
#python3 tests/__init__.py