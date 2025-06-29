#include <benchmark/benchmark.h>

#include "../base/core.hpp"
#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"
#include "../ops/newick.hpp"
#include "../ops/vector.hpp"

// Benchmark sample
static void BM_sample(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        PhyloVec v = sample(n, false);
        benchmark::DoNotOptimize(v);
        benchmark::ClobberMemory();
    }
}

// Benchmark toNewick
static void BM_toNewick(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        PhyloVec v = sample(n, false);
        state.ResumeTiming();
        std::string newick = toNewick(v);
        benchmark::DoNotOptimize(newick);
        benchmark::ClobberMemory();
    }
}

// Benchmark toVector
static void BM_toVector(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        PhyloVec v1 = sample(n, false);
        std::string newick = toNewick(v1);
        state.ResumeTiming();
        PhyloVec v2 = toVector(newick);
        benchmark::DoNotOptimize(v2);
        benchmark::ClobberMemory();
    }
}

// Benchmark toVectorNoParents
static void BM_toVectorNoParents(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        PhyloVec v1 = sample(n, false);
        std::string newick = toNewick(v1, false);
        state.ResumeTiming();
        PhyloVec v2 = toVectorNoParents(newick);
        benchmark::DoNotOptimize(v2);
        benchmark::ClobberMemory();
    }
}

#define QUICK_RANGE Range(8 << 6, 8 << 12)->Unit(benchmark::kMillisecond)
#define BIG_RANGE DenseRange(10000, 100000, 10000)->Unit(benchmark::kMillisecond)

BENCHMARK(BM_sample)->BIG_RANGE;
BENCHMARK(BM_toNewick)->BIG_RANGE;
BENCHMARK(BM_toVector)->BIG_RANGE;
BENCHMARK(BM_toVectorNoParents)->BIG_RANGE;

// Run the benchmark
BENCHMARK_MAIN();