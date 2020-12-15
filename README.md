# wacommplusplus
WaComM++ is a steroidized version of WaComM - Water quality Community Model.

# Compiling

WaComM++ is developed using C++17. Be sure a compatible toolchain is installed and available.
* Linux CentOS:
  yum install devtoolset-9
  
* MacOS:
  CLang
    
* Windows: Visual Sutudio

## Using the command line interface

1) Clone the repository
```bash
git clone https://github.com/CCMMMA/wacommplusplus.git
```
2) Enter the project directory
```bash
cd wacomplusplus
```
3) Create the build directory
```bash
mkdir build
```
4) Enter the build directory
```bash
cd build
```
5) Invoke cmake
```bash
cmake ..
```
6) Run make and wait
```bash
make
```
## Using CLion 

1) Clone a new project from CVS
2) Enter https://github.com/CCMMMA/wacommplusplus.git as project URL
3) Press the Clone Button
4) Select *wacomplusplus* | Debug as target
5) Click on menu *Build*, option *Build Project*
6) Wait


# Testing
Download the data.

## CLion

1) Select *data | Debug*
2) Click on menu *Build*, option *Build data*
3) Wait
