#include "../base/core.hpp"
#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"
#include "../utils/newick.hpp"
#include "../utils/random.hpp"

#include <benchmark/benchmark.h>

// Benchmark sample
static void BM_sample(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        PhyloVec v = sample(n);
        benchmark::DoNotOptimize(v);
        benchmark::ClobberMemory();
    }
}

// Benchmark toNewick
static void BM_toNewick(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        PhyloVec v = sample(n);
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
        PhyloVec v1 = sample(n);
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
        PhyloVec v1 = sample(n);
        std::string newick = toNewick(v1);
        removeParentLabels(newick);
        state.ResumeTiming();
        PhyloVec v2 = toVectorNoParents(newick);
        benchmark::DoNotOptimize(v2);
        benchmark::ClobberMemory();
    }
}

#define BENCHMARK_RANGE Range(8 << 6, 8 << 12)->Unit(benchmark::kMillisecond)

BENCHMARK(BM_sample)->BENCHMARK_RANGE;
BENCHMARK(BM_toNewick)->BENCHMARK_RANGE;
BENCHMARK(BM_toVector)->BENCHMARK_RANGE;
BENCHMARK(BM_toVectorNoParents)->BENCHMARK_RANGE;

// Run the benchmark
BENCHMARK_MAIN();