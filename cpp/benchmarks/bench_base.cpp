#include "../base/core.hpp"
#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"
#include "../utils/newick.hpp"
#include "../utils/random.hpp"

#include <benchmark/benchmark.h>

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

// Register the benchmark
BENCHMARK(BM_toNewick)->Arg(1024)->Unit(benchmark::kMillisecond)->Complexity();

// BENCHMARK(BM_toNewick)->RangeMultiplier(2)->Range(2, 2 <<
// 10)->Unit(benchmark::kMillisecond);

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

// Register the benchmark
BENCHMARK(BM_toVector)->Arg(1024)->Unit(benchmark::kMillisecond)->Complexity();
// BENCHMARK(BM_toVector)->RangeMultiplier(2)->Range(2, 2 <<
// 10)->Unit(benchmark::kMillisecond);

void benchToVectorNoParents(std::string newick) {
    // PhyloVec v = toVectorNoParents(newick);

    toVectorNoParents(newick);
}

static void BM_toVectorNoParents(benchmark::State &state) {
    int n = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        PhyloVec v = sample(n);
        std::string newick = toNewick(v);
        removeParentLabels(newick);
        // Ancestry anc = reduceNoParents(newick);
        state.ResumeTiming();
        benchToVectorNoParents(newick);
    }
}

BENCHMARK(BM_toVectorNoParents)
    ->Arg(1024)
    ->Unit(benchmark::kMillisecond)
    ->Complexity();

// Run the benchmark
BENCHMARK_MAIN();