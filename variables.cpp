#include <iostream>
#include <fstream>
#include <random>
#include "systemparam.hpp"
#include "variables.hpp"

Variables::Variables() : time(0.0) { ; }

void Variables::add_atoms(double x, double y, double z) {
    atoms.emplace_back(Atom(x, y, z, 0, 0, 0));
}

void Variables::set_initial_velocity(const double V0) { 
    std::mt19937 mt(2);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    double avx = 0;
    double avy = 0;
    double avz = 0;
    for(auto& atom: atoms) {
        double z = ud(mt) * 2.0 - 1.0;
        double phi = 2.0 * ud(mt) * M_PI;
        double vx = V0 * sqrt(1 - z * z) * cos(phi);
        double vy = V0 * sqrt(1 - z * z) * sin(phi);
        double vz = V0 * z;
        atom.px = vx;
        atom.py = vy;
        atom.pz = vz;
        avx += vx;
        avy += vy;
        avz += vz;
    }
    const int pn = atoms.size();
    avx /= static_cast<double>(pn);
    avy /= static_cast<double>(pn);
    avz /= static_cast<double>(pn);
    for(auto& atom: atoms) {
        atom.px -= avx;
        atom.py -= avy;
        atom.pz -= avz;
    }
}

void Variables::export_cdview() {
    static int count = 0;
    char filename[256];
    sprintf(filename, "conf%03d.cdv", count);
    ++count;
    std::ofstream ofs(filename);
    int i = 0;
    for(auto& atom: atoms) {
        ofs << i << " " << "0" << " ";
        ofs << atom.qx << " ";
        ofs << atom.qy << " ";
        ofs << atom.qz << " ";
        ofs << std::endl;
        ++i;
    }
}

int Variables::number_of_atoms() const {
    return static_cast<int>(atoms.size());
}