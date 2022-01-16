mkdir build 
cd build
cmake .. -G "Visual Studio 16 2019"
msbuild ALL_BUILD.vcxproj /p:Configuration=Debug
cd ..
