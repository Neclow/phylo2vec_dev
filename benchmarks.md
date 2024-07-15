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
# Modify cpp/benchmarks/bench.cpp accordingly
cd cpp
mkdir build
cd build
cmake ..
make
./phylo2vec_bench
```

## Summary table

Approx. execution times (ms)

|                      | size | cpp  | python |
|----------------------|------|------|--------|
| to_newick            | 512  | 0.17 | 0.4    |
| to_newick            | 1024 | 0.4  | 0.85   |
| to_vector            | 512  | 0.7  | 2.5    |
| to_vector            | 1024 | 2.5  | 6.5    |
| to_vector_no_parents | 512  | 3.2  | 3      |
| to_vector_no_parents | 1024 | 12.5 | 8.9    |
