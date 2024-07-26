//
// Created by ilia on 01/07/24.
//


#include <gtest/gtest.h>

#include "Atoms.h"
#include "ducastelle.h"
#include "neighbors.h"
#include "xyz.h"
#include "mpi.h"
#include "berendsten.h"

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
    Atoms at2(positions);
    auto [maxPos, minPos] = get_max_pos(at.positions);
    Eigen::Array3d a(3,1);
    a << maxPos[0]*3, maxPos[1]*3, std::ceil(maxPos[2]+1);

    // initiate standard values for gold, stable timestep
    const double timestep = 2;
    const double fixed_mass = 196.96657 / 0.009649;
    NeighborList nl;
    nl.update(at,5);

    // initiate and test for set tempearture for berendsen single to multi test
    at.velocities.setRandom();
    at.velocities.colwise().normalize();
    at.velocities*= std::sqrt( (3*kB* 50)/ fixed_mass );

    // copy instance to verify results
    at2 = at;

    double ekin = at.get_ekin(fixed_mass, at.velocities);
    double T = ( ekin / (1.5 * at.positions.cols() * kB) );
    EXPECT_NEAR(50., T, 1e-10);
    for (int i = 0; i < at.velocities.cols();i++) {
        EXPECT_NEAR(at.velocities(i),at2.velocities(i),1e-5);
    }
    double lambda = berendsen_thermostat(at2,T,timestep,1000, fixed_mass);
    std::cout << at.nb_atoms() << std::endl;
    // get epot from single core for Test
    double epot_single{ducastelle(at2,nl)};

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
                  {1,1,6},
                  {0, 0, 0});
    domain.enable(at);
    domain.update_ghosts(at,10.);
    nl.update(at,5);

    // get epot from all cores
    std::cout << "numbers: "<< MPI::allreduce(domain.nb_local(),MPI_SUM, MPI_COMM_WORLD) << std::endl;
    double epot_tot{MPI::allreduce(ducastelle(at,nl, domain),
                                   MPI_SUM, MPI_COMM_WORLD)};
    double lambda2 = berendsen_thermostat(at,domain,T,timestep,1000,fixed_mass);

    MPI_Finalize();
    EXPECT_NEAR(epot_single, epot_tot, 1e-10);
    for (int i = 0; i < at.velocities.cols();i++) {
        EXPECT_NEAR(at.velocities(i),at2.velocities(i),1e-2);
    }
    EXPECT_NEAR(lambda, lambda2, 1e-1);
}