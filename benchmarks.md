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

|                      | size | cpp   | python | julia | winner |
|----------------------|------|------ |--------|-------|--------|
| sample               | 512  | 0.004 | 0.009  | 0.003 | julia  |
| sample               | 1024 | 0.007 | 0.013  | 0.006 | julia  |
| to_newick            | 512  | 0.16  | 0.39   | 0.24  | cpp    |
| to_newick            | 1024 | 0.34  | 0.85   | 0.64  | cpp    |
| to_vector            | 512  | 0.38  | 1.7    | 0.4   | cpp    |
| to_vector            | 1024 | 1.2   | 5.7    | 1.5   | cpp    |
| to_vector_no_parents | 512  | 1.9   | 2.2    |       | cpp    |
| to_vector_no_parents | 1024 | 6.9   | 6.8    |       | python |
