rm -r build
mkdir build
cd build
cmake ..
make
cd ..
export PYTHONPATH=$PWD/build/src/:$PYTHONPATH