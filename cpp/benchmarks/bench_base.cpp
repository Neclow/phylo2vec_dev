#include "../base/core.hpp"
#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"
#include "../utils/newick.hpp"
#include "../utils/random.hpp"

#include <benchmark/benchmark.h>

constexpr int NUM_LEAVES = 1024;

// Benchmark sample
void benchSample(PhyloVec v) { std::string newick = toNewick(v); }

static void BM_sample(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        PhyloVec v = sample(n);
    }
}

// Benchmark toNewick
void benchToNewick(PhyloVec v) { std::string newick = toNewick(v); }

static void BM_toNewick(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        PhyloVec v = sample(n);
        state.ResumeTiming();
        benchToNewick(v);
    }
}

// Benchmark toVector
void benchToVector(std::string newick) { PhyloVec v = toVector(newick); }

static void BM_toVector(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        PhyloVec v = sample(n);
        std::string newick = toNewick(v);
        state.ResumeTiming();
        benchToVector(newick);
    }
}

// Benchmark toVectorNoParents
void benchToVectorNoParents(std::string newick) {
    PhyloVec v = toVectorNoParents(newick);
}

static void BM_toVectorNoParents(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        PhyloVec v = sample(n);
        std::string newick = toNewick(v);
        removeParentLabels(newick);
        state.ResumeTiming();
        benchToVectorNoParents(newick);
    }
}

BENCHMARK(BM_sample)
    ->Arg(NUM_LEAVES)
    ->Unit(benchmark::kMillisecond)
    ->Complexity();

BENCHMARK(BM_toNewick)
    ->Arg(NUM_LEAVES)
    ->Unit(benchmark::kMillisecond)
    ->Complexity();

BENCHMARK(BM_toVector)
    ->Arg(NUM_LEAVES)
    ->Unit(benchmark::kMillisecond)
    ->Complexity();

BENCHMARK(BM_toVectorNoParents)
    ->Arg(NUM_LEAVES)
    ->Unit(benchmark::kMillisecond)
    ->Complexity();

// Run the benchmark
BENCHMARK_MAIN();