# rs@md
reactive steps @ molecular dynamics

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

rs@md is a wrapper program written in C++ 2017. It wraps around a molecular dynamics engine, i.e. it can read, manipulate and write input files, start a molecular dynamics sequence with that engine, wait for it to finish and subsequently read the output of that last MD run. Thereby, rs@md introduces reactive steps into classical molecular dynamics simulations.

Currently, the only supported molecular dynamics engine is [GROMACS](https://gitlab.com/gromacs/gromacs).

### Requirements
- Boost Program Options

### Installation guide
```bash
git clone https://github.com/myrabiedermann/rsmd.git
cd rsmd
mkdir build
cd build
cmake ..
```

### Quick start
Building rs@md results in an executable named "rsmd", which can be executed in the command line via 
```bash
./rsmd
```
Additional information are provided when using the --help flag
```bash
./rsmd --help
```



