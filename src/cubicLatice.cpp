//
// Created by ilia on 25/06/24.
//

#include "Atoms.h"

void createCubicLatice(Atoms &at, int num) {
    int counter = 0;
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < num; j++) {
            for (int k = 0; k < num; k++) {
                at.positions.col(counter) << i, j, k;
                counter++;
            }
        }
    }
}
