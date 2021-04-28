# Physicell-Cmake
This project is an attempt at making a Cmake-compiled version of the PhysiCell code base. 

# Build/Installation Process
CMake is the build system for this project. You will likely need to configure your build according to the specifics of your system. Basic CMake usage for this project is as follows:
```
$cd YOUR_Physicell-Cmake_Directory
$mkdir build
$cmake -DCMAKE_INSTALL_PREFIX=your_path_to_a_directory_for_install
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




