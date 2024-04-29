//
// Created by ilia on 19/04/24.
//

#include "verlet.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
    double fx, double fy, double fz, double timestep, double mass) {
// void verlet_step1(Atoms &atoms, double timestep, double mass) {
    vx += 0.5 * fx * timestep / mass;
    vy += 0.5 * fy * timestep / mass;
    vz += 0.5 * fz * timestep / mass;
    x += vx * timestep;
    y += vy * timestep;
    z += vy * timestep;
    // atoms.positions += atoms.velocities * timestep;
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass) {
//void verlet_step2(Atoms &atoms, double timestep, double mass) {
    // atoms.velocities += 0.5 * atoms.forces * timestep / mass;
    vx += 0.5 * fx * timestep/mass;
    vy += 0.5 * fy * timestep/mass;
    vz += 0.5 * fz * timestep/mass;
}