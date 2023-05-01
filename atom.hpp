#ifndef ATOM_HPP_
#define ATOM_HPP_

struct Atom {
    double qx, qy, qz;
    double px, py, pz;
    Atom();
    Atom(double, double, double, double, double, double);
};

#endif
