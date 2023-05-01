#ifndef OBSERVER_HPP_
#define OBSERVER_HPP_

#include "variables.hpp"

struct Observer {
    double kinetic_energy(const Variables&);
    double potential_energy(const Variables&);
    double temperature(const Variables&);
    double total_energy(const Variables&);
};

#endif