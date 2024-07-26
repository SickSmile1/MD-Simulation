//
// Created by ilia on 17/07/24.
//

#include "vector"
#include "berendsten.h"
#include "verlet.h"
#include "domain.h"
#include "ducastelle.h"
#include "neighbors.h"
#include "iostream"
#include "xyz.h"
#include "mpi.h"

int main(int argc, char** argv) {
    auto [names, positions, velocities]{read_xyz_with_velocities("whisker_small.xyz")};
    Atoms at(positions);
    const double m = 196.96657 / 0.009649;
    at.velocities.setRandom();
    at.velocities.colwise().normalize();
    at.velocities*= std::sqrt( (3*kB* 300)/ m );

    double ekin = at.get_ekin(m, at.velocities);
    double T = ( ekin / (1.5 * at.positions.cols() * kB) );
    std::cout<< "initial temp is 300?: " << T << std::endl;
    Eigen::Array3Xd one{3,1};
    Eigen::Array3Xd two{3,1};
    Eigen::Vector3d one1{1,2,3};
    Eigen::Vector3d two1{1,2,3};

    one<<1,2,3;
    two<<1,2,3;
    // multi = one*two.transpose();


    std::cout << one << std::endl;
    std::cout << two << std::endl;
    // two = two.transposeInPlace();
    std::cout << one.matrix()*two.matrix().transpose() << std::endl;

    std::cout <<"und als vektoren: \n"<< two1*one1.transpose() << std::endl;

    /*MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Domain domain(MPI_COMM_WORLD,
                  {
                          50,
                          50,
                          144.24978336},
                  {1,1,6},
                  {0, 0, 0});
    domain.enable(at);
    domain.update_ghosts(at,11.);
    NeighborList nl;
    nl.update(at,5.5);

    double ekin_total{MPI::allreduce(at.get_ekin(m, at.velocities, domain.nb_local()), MPI_SUM, MPI_COMM_WORLD)};
    double ekin1 = at.get_ekin(m, at.velocities, domain.nb_local());
    int n_atoms{MPI::allreduce(domain.nb_local(), MPI_SUM, MPI_COMM_WORLD)};
    std::cout << at.nb_atoms() << " domain atoms: " << domain.nb_local() << std::endl;
    if (rank ==0) {
        // std::cout << "globpot" << global_pot << "reduce: " << global_pot+ekin << std::endl;
        double T1 = ekin_total / (1.5 * n_atoms * kB);
        std::cout << "t1: " << T << " t2: " << T1 << std::endl;
    }
    MPI_Finalize();*/
}
