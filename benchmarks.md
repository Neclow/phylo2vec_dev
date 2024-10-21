# Benchmark stats

## Python example

```ipython
from phylo2vec.base import to_newick
from phylo2vec.utils import sample
# Sample a 1024-leaf tree
v = sample(1024)
# Compile
to_newick(v)
# Run timeit
%timeit to_newick(v);
```

## C++ example

```bash
# Modify cpp/benchmarks/bench_base.cpp accordingly
cd cpp
mkdir build
cd build
cmake ..
make
./phylo2vec_bench
```

## Julia example

```julia
using BenchmarkTools

include("jl/base.jl")
include("jl/to_newick.jl")

# Sample a 1024-leaf tree
v = sample_v(1024)
# Compile
to_newick(v)

# Benchmark performance
@benchmark to_newick(v) setup=(v=$v)
```

## Summary table

Approx. execution times (ms)

|                      | size  | cpp   | python | julia | winner |
|----------------------|------ |------ |--------|-------|--------|
| sample               | 512   | 0.004 | 0.009  | 0.003 | julia  |
| sample               | 1024  | 0.007 | 0.013  | 0.006 | julia  |
| sample               | 4096  | 0.037 | 0.098  | 0.023 | julia  |
| sample               | 32768 | 0.294 | 0.728  | 0.222 | julia  |
| to_newick            | 512   | 0.19  | 0.22   | 0.24  | cpp    |
| to_newick            | 1024  | 0.34  | 0.47   | 0.64  | cpp    |
| to_newick            | 4096  | 1.89  | 2.29   | 4.6  | cpp    |
| to_newick            | 32768 | 34.0  | 27     | 117.6 | python |
| to_vector            | 512   | 0.38  | 1.7    | 0.4   | cpp    |
| to_vector            | 1024  | 1.2   | 5.7    | 1.5   | cpp    |
| to_vector            | 4096  | 19.0  |        | 12.8  | julia  |
| to_vector            | 32768 | 1005  |        | 871.0 | julia  |
| to_vector_no_parents | 512   | 1.9   | 2.2    |       | cpp    |
| to_vector_no_parents | 1024  | 6.9   | 6.8    |       | python |
| to_vector_no_parents | 4096  | 117   |        |       |        |
| to_vector_no_parents | 32768 | 8992  |        |       |        |

|                      | size     | AVL   | old    |
|----------------------|----------|------ |--------|
| getPairs             | 512      | 0.188 | 0.042  |
| getPairs             | 4096     | 1.90  | 0.446  |
| getPairs             | 32768    | 19.1  | 15.4   |
| getPairs             | 262144   | 188   | 1020   |
| getPairs             | 2097152  | 2130  | 69144  |
