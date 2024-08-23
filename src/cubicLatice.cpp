//
// Created by ilia on 25/06/24.
//

#include "Atoms.h"

void createCubicLatice(Atoms &at, const int num, const double dist = 1) {
    // simple cubic
    int counter = 0;
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < num; j++) {
            for (int k = 0; k < num; k++) {
                at.positions.col(counter) << i*dist, j*dist, k*dist;
                counter++;
            }
        }
    }
}

