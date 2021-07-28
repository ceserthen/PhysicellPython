# Physicell-Cmake
This project is an attempt at making a Cmake-compiled version of the PhysiCell code base. 

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
# Running the Executable
Once you have completed installation you can go to `your_path_to_a_directory_for_install`. From here you have a folder for each sample project.

# Heterogeneity
In the `heterogeneity` folder you will have the executable (`heterogeneity`), `config` folder, and `output` folder. Settings and default initial conditions can be adjusted in the config folders. You can run the executable in the target directory by:
```
$.heterogeneity
```
After running your output files will be in the `output` directory. 

## Mac
Currently untested




