# C++ version

## Installation

Prerequisites:

* C++14
* GoogleTest 1.11.0: ```sudo apt-get install libgtest-dev```
* clang-format: ```sudo apt install clang-format```
* cmake 3.22.1: ```sudo apt-get install cmake```
* benchmark 1.8.4: see <https://github.com/google/benchmark?tab=readme-ov-file#installation>

Compile and build from scratch:

```bash
mkdir build
cd build
cmake ..
make
```

## Tests

```bash
./phylo2vec_test
```

## Benchmarks

```bash
./phylo2vec_bench
```

## Python

```bash
pip install -e .
```

## R

```bash
source("R/build.R")
```
