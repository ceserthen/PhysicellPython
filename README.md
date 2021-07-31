# PhysicellPython
This project has two primary goals: cmake compilation of PhysiCell and a resulting Python API for building physicell projects. 

# Build/Installation Process
## Windows
### 64-bit gcc/OpenMP
Before installing the software, I reccomend following the tutorial on how to [set up a 64-bit gcc/OpenMP environment on Windows](http://www.mathcancer.org/blog/setting-up-a-64-bit-gcc-environment-on-windows/). 

### CMake
Once you've completed that you will need to install cmake. There are several [precompliled binaries available here](https://cmake.org/download/).

### Using CMake for Install
CMake is the build system for this project. You will likely need to configure your build according to the specifics of your system. Basic CMake usage for this project is as follows:
```
>mkdir build
>cmake path_to_your_source_directory -DCMAKE_INSTALL_PREFIX=your_path_to_a_directory_for_install -DCMAKE_BUILD_TYPE=Release
>cmake --build .
>cmake --install .
```
### Running the Executable
Once you have completed installation you can go to `your_path_to_a_directory_for_install`. From here you have a folder for each sample project.

### Heterogeneity
In the `heterogeneity` folder you will have the executable (`heterogeneity.exe`), `config` folder,`output` folder, and associated `.dll` files. Settings and default initial conditions can be adjusted in the config folders. You can run the executable in the target directory by:
```
>heterogeneity.exe
```

## Linux
CMake is the build system for this project. You will likely need to configure your build according to the specifics of your system. Basic CMake usage for this project is as follows:
```
$cd YOUR_Physicell-Cmake_Directory
$mkdir build
cd build
$cmake path_to_source_code_directory -DCMAKE_INSTALL_PREFIX=your_path_to_a_directory_for_install
$make install
```
### Running the Executable
Once you have completed installation you can go to `your_path_to_a_directory_for_install`. From here you have a folder for each sample project.

### Heterogeneity
In the `heterogeneity` folder you will have the executable (`heterogeneity`), `config` folder, and `output` folder. Settings and default initial conditions can be adjusted in the config folders. You can run the executable in the target directory by:
```
$.heterogeneity
```
After running your output files will be in the `output` directory. 

## Mac

### CMake
CMake is the build system for this project. There are several [precompliled binaries available here](https://cmake.org/download/).You will likely need to configure your build according to the specifics of your system. Basic CMake usage for this project is as follows:
```
~/git/PhysicellPython$ mkdir build
~/git/PhysicellPython$ cd build
~/git/PhysicellPython/build$ mkdir installed
```
this installation provides custom compilers and paths
```
cmake  -DCMAKE_CXX_COMPILER=g++-11 -DCMAKE_C_COMPILER=gcc-11 \
-DCMAKE_INSTALL_PREFIX:PATH=installed  \
-DPYTHON_EXECUTABLE:FILEPATH=/Users/USERNAME/opt/anaconda3/bin/python \
-DPYTHON_INCLUDE_DIR:PATH=/Users/USERNAME/opt/anaconda3/include/python3.8 \
-DPYTHON_LIBRARY:FILEPATH=/Users/USERNAME/opt/anaconda3/lib/libpython3.8.dylib ..
```
After setting the cmake configuration you can make and install the project.
```
make -j2
make install   # this will put everything in ./installed (specified above)
```
