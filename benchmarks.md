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
include("jl/base.jl")

# Sample a 1024-leaf tree
v = sample_v(1024)
# Compile
to_newick(v)

# Benchmark performance
using BenchmarkTools
@benchmark to_newick(v) setup=(v=$v)
```

## Summary table

Approx. execution times (ms)

|                      | size | cpp   | python | julia | winner |
|----------------------|------|------ |--------|-------|--------|
| sample               | 512  | 0.009 | 0.009  | 0.003 | julia  |
| sample               | 1024 | 0.015 | 0.013  | 0.006 | julia  |
| to_newick            | 512  | 0.18  | 0.39   | 0.24  | cpp    |
| to_newick            | 1024 | 0.38  | 0.85   | 0.64  | cpp    |
| to_vector            | 512  | 0.65  | 2.5    | 0.4   | julia  |
| to_vector            | 1024 | 2.2   | 6.5    | 1.5   | julia  |
| to_vector_no_parents | 512  | 2.66  | 3      |       | cpp    |
| to_vector_no_parents | 1024 | 9.7   | 8.9    |       | python |