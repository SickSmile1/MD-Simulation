//
// Created by ilia on 28/06/24.
//
#include <benchmark/benchmark.h>
#include "lj_direct_summation.h"
#include "lj.h"
#include "cubicLatice.h"
#include "berendsten.h"
#include "neighbors.h"
#include "verlet.h"


// 8 27 64 125 216 343 512 729 1000

static void lj_cubes(benchmark::State& state) {
    for (auto _: state) {
        state.PauseTiming();
        state.SetComplexityN(std::pow(state.range(0), 3));
        Atoms at(std::pow(state.range(0), 3));
        const double timestep = 1e-3;
        createCubicLatice(at, state.range(0));
        state.ResumeTiming();
        for (int i = 0; i < 100; i++) {
            verlet_step1(at.positions, at.velocities, at.forces, timestep, 1);
            lj_direct_summation(at, 1., 1.);
            verlet_step2(at.velocities, at.forces, timestep);
            berendsen_thermostat(at, 0.6, timestep, 1e-2);
        }
    }
}


BENCHMARK(lj_cubes)
->ArgsProduct({
benchmark::CreateDenseRange(2, 8, /*step=*/1) // This creates a DenseRange from 1 to 8
})->Complexity(benchmark::oNSquared)->Unit(benchmark::kMillisecond);

static void lj_cubes_cutoff(benchmark::State& state) {
    for(auto _:state) {
        state.PauseTiming();
        state.SetComplexityN(std::pow(state.range(0), 3));
        Atoms at(std::pow(state.range(0), 3));
        const double timestep = 1e-3;
        createCubicLatice(at, state.range(0));
        NeighborList nl;
        nl.update(at,1.5);
        state.ResumeTiming();
        for (int i = 0; i < 100; i++) {
            verlet_step1(at.positions, at.velocities, at.forces, timestep, 1);
            lj_direct_summation(at, nl, 1., 1., 1.5);
            verlet_step2(at.velocities, at.forces, timestep);
            berendsen_thermostat(at, 0.6, timestep, 1e-2);
            nl.update(at,1.);
        }
    }
}

BENCHMARK(lj_cubes_cutoff)
    ->ArgsProduct({
    benchmark::CreateDenseRange(2, 8, /*step=*/1) // This creates a DenseRange from 1 to 8
})->Complexity(benchmark::oNSquared)->Unit(benchmark::kMillisecond);

// MY_BENCHMARK(lj_cubes);

BENCHMARK_MAIN();