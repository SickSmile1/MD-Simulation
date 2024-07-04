//
// Created by ilia on 28/06/24.
//
#include "Atoms.h"
#include "neighbors.h"
#include "xyz.h"
#include "lj_direct_summation.h"

#include <gtest/gtest.h>

TEST(NeighborsLJTest, stillWorks) {
    constexpr int nb_atoms = 10;
    constexpr double epsilon = 0.7;  // choose different to 1 to pick up missing factors
    constexpr double sigma = 0.3;
    constexpr double delta = 0.0001;  // difference used for numerical (finite difference) computation of forces

    Atoms atoms(nb_atoms);
    atoms.positions.setRandom();  // random numbers between -1 and 1
    NeighborList nl;
    nl.update(atoms, 4.);
    // compute and store energy of the indisturbed configuration
    double e0{lj_direct_summation(atoms, nl, epsilon, sigma)};
    Forces_t forces0{atoms.forces};

    // loop over all atoms and compute forces from a finite differences approximation
    for (int i{0}; i < nb_atoms; ++i) {
        // loop over all Cartesian directions
        for (int j{0}; j < 3; ++j) {
            // move atom to the right
            atoms.positions(j, i) += delta;
            double eplus{lj_direct_summation(atoms, nl, epsilon, sigma)};
            // move atom to the left
            atoms.positions(j, i) -= 2 * delta;
            double eminus{lj_direct_summation(atoms, nl, epsilon, sigma)};
            // move atom back to original position
            atoms.positions(j, i) += delta;

            // finite-differences forces
            double fd_force{-(eplus - eminus) / (2 * delta)};

            // check whether finite-difference and analytic forces agree
            if (abs(forces0(j, i)) > 1e-10) {
                EXPECT_NEAR(abs(fd_force - forces0(j, i)) / forces0(j, i), 0, 1e-5);
            } else {
                EXPECT_NEAR(fd_force, forces0(j, i), 1e-10);
            }
        }
    }
}