//
// Created by ilia on 01/07/24.
//


#include <gtest/gtest.h>

#include "Atoms.h"
#include "ducastelle.h"
#include "neighbors.h"
#include "xyz.h"
#include "mpi.h"

std::tuple<std::vector<double>,std::vector<double>> get_max_pos(Eigen::Array3Xd &p) {
    std::vector<double> max, min;
    max.push_back(p.row(0).maxCoeff());
    max.push_back(p.row(1).maxCoeff());
    max.push_back(p.row(2).maxCoeff());
    min.push_back(p.row(0).minCoeff());
    min.push_back(p.row(1).minCoeff());
    return std::tuple(max,min);
}

int main(int argc, char** argv) {
    // Filter out Google Test arguments
    ::testing::InitGoogleTest(&argc, argv);

    // read whisker and initiate atoms and domain length
    auto [names, positions, velocities]{read_xyz_with_velocities("whisker_small.xyz")};
    Atoms at(positions);
    auto [maxPos, minPos] = get_max_pos(at.positions);
    Eigen::Array3d a(3,1);
    a << maxPos[0]*3, maxPos[1]*3, std::ceil(maxPos[2]+1);

    // initiate standard values for gold, stable timestep
    const double timestep = 8;
    const double fixed_mass = 196.96657 / 0.009649;
    NeighborList nl;
    nl.update(at,5);

    // initiate and test for set tempearture for berendsen single to multi test
    at.velocities.setRandom();
    at.velocities.colwise().normalize();
    at.velocities*= std::sqrt( (3*kB* 50)/ fixed_mass );

    double ekin = at.get_ekin(fixed_mass, at.velocities);
    double T = ( ekin / (1.5 * at.positions.cols() * kB) );
    EXPECT_NEAR(50., T, 1e-10);

    // get epot from single core for Test
    double epot_single{ducastelle(at,nl)};

    // initiate mpi-simulation
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Domain domain(MPI_COMM_WORLD,
                  {
                      a[0],
                      a[1],
                      a[2]},
                  {1,1,12},
                  {0, 0, 1});
    domain.enable(at);
    domain.update_ghosts(at,10.);
    nl.update(at,5);

    // get epot from all cores
    double epot_tot{MPI::allreduce(ducastelle(at,nl, domain),
                                   MPI_SUM, MPI_COMM_WORLD)};
    EXPECT_NEAR(epot_single, epot_tot, 1e-10);



    MPI_Finalize();
}