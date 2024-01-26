# M-SPDZ
the implementation of the matrix SPDZ

## emp-tool
![arm](https://github.com/emp-toolkit/emp-tool/workflows/arm/badge.svg)
![x86](https://github.com/emp-toolkit/emp-tool/workflows/x86/badge.svg)
[![CodeQL](https://github.com/emp-toolkit/emp-tool/actions/workflows/codeql.yml/badge.svg)](https://github.com/emp-toolkit/emp-tool/actions/workflows/codeql.yml)

<img src="https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/art/logo-full.jpg" width=300px/>



## Installation
1. `wget https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/scripts/install.py`
2. `python install.py --deps --tool `
    1. You can use `--ot=[release]` to install a particular branch or release
    2. By default it will build for Release. `-DCMAKE_BUILD_TYPE=[Release|Debug]` option is also available.
    3. No sudo? Change [`CMAKE_INSTALL_PREFIX`](https://cmake.org/cmake/help/v2.8.8/cmake.html#variable%3aCMAKE_INSTALL_PREFIX).
   
This program is implemented based on emp-tool. The following is an introduction to each folder in the project.

## network

## offline 
writing......

## online
spdz.h cantains the basic spdz protocol

m-spdz.h is the protocol designed in the [paper](https://eprint.iacr.org/2023/1912/.)

## test 
已实现线上测试
> g++ test/mul.cpp -o mul

打开三个线程，分别执行
> ./mul 1 12345
> ./mul 2 12345
> ./mul 3 12345

对于矩阵乘法测试 类似

由于预处理数据是自己随机写的，所以check error

矩阵乘法时间快了0.9m倍。
## pre_data
测试数据，尝试用mp-spdz生成中(一直报错)