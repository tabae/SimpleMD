#include "observer.hpp"
#include "systemparam.hpp"

double Observer::kinetic_energy(const Variables& vars) {
    double k = 0;
    for(const auto& atom: vars.atoms) {
        k += atom.px * atom.px;
        k += atom.py * atom.py;
        k += atom.pz * atom.pz;
    }
    k /= static_cast<double>(vars.number_of_atoms());
    return k * 0.5;
}

double Observer::potential_energy(const Variables& vars) {
    double v = 0.0;
    const int pn = vars.number_of_atoms();
    const std::vector<Atom>& atoms = vars.atoms;
    for(int i = 0; i < pn; i++) {
        for(int j = i + 1; j < pn; j++) {
            double dx = atoms[j].qx - atoms[i].qx;
            double dy = atoms[j].qy - atoms[i].qy;
            double dz = atoms[j].qz - atoms[i].qz;
            adjust_periodic(dx, dy, dz);
            double r2 = dx * dx + dy * dy + dz * dz;
            if(r2 > CL2) continue;
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            v += 4.0 * (1.0 / r12 - 1.0 / r6) + C0;
        }
    }
    v /= static_cast<double>(pn);
    return v;
}

double Observer::temperature(const Variables& vars) {
    return kinetic_energy(vars) / 1.5;
}

double Observer::total_energy(const Variables& vars) {
    return kinetic_energy(vars) + potential_energy(vars);
}
