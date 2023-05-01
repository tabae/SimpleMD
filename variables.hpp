#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

#include "atom.hpp"
#include <vector>

struct Variables {
    double time;
    Variables();
    std::vector<Atom> atoms;
    void add_atoms(double x, double y, double z);
    void set_initial_velocity(const double V0);
    void export_cdview();
    int number_of_atoms() const;
};

#endif
