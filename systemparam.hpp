#ifndef SYSTEMPARAM_HPP_
#define SYSTEMPARAM_HPP_

const double L = 10;
const double dt = 0.01;
const double CUTOFF = 2.0;
const double CL2 = (CUTOFF*CUTOFF);
const double RC2 = 1.0 / CL2;
const double RC6 = RC2 * RC2 * RC2;
const double RC12 = RC6 * RC6;
const double C0 = - 4.0 * (RC12 - RC6);
void adjust_periodic(double &dx, double &dy, double &dz);

#endif